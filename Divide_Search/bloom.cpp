#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>

#include "secp256k1/SECP256k1.h"
#include "secp256k1/Int.h"
#include "secp256k1/IntGroup.h"
#include "util/util.h"

using namespace std;

const int n_cores = 2;
const int POINTS_BATCH_SIZE = 1024; // Batch addition with batch inversion using IntGroup class

std::mutex mtx_lock;

auto main() -> int {
    
    auto chrono_start = std::chrono::high_resolution_clock::now();    // starting the timer
    Secp256K1* secp256k1 = new Secp256K1();
    secp256k1->Init(); // initializing secp256k1 context
    
    Int pk; pk.SetInt32(1); // generating power of two points table (2^0..2^256) 
    uint64_t mult = 2;
    vector<Point> P_table;
    Point P;
    for (int i = 0; i < 256; i++)
    {
        P = secp256k1->ScalarMultiplication(&pk);
        P_table.push_back(P);
        pk.Mult(mult);
    }
    print_time(); cout << "P_table generated" << endl;

    uint64_t range_start, range_end, block_width; // block_width = number of elements in the bloomfilter and a stride size to walk the range
    string temp, search_pub;
    ifstream inFile("settings.txt");
    getline(inFile, temp); range_start = std::stoull(temp);
    getline(inFile, temp); range_end = std::stoull(temp);
    getline(inFile, temp); block_width = std::stoull(temp);
    getline(inFile, temp); search_pub = trim(temp);
    inFile.close();
    
    print_time(); cout << "Range Start: " << range_start << " bits" << endl;
    print_time(); cout << "Range End  : " << range_end << " bits" << endl;
    print_time(); cout << "Block Width: 2^" << block_width << endl;
    print_time(); cout << "Search Pub : " << search_pub << endl;

    Point target_point = secp256k1->ParsePublicKeyHex(search_pub);

    ofstream outFile; // writing steps_sum to file (initial value = 0)
    outFile.open("steps_sum.txt");
    outFile << "0" << '\n';
    outFile.close();
    
    print_time(); cout << "Settings written to file" << endl;
    
    //uint64_t n_elements = pow(2, block_width);  // number of elements == 2^block_width
    uint64_t n_elements = 1 << block_width;  // number of elements == 2^block_width
    uint64_t keysPerThread = n_elements / n_cores; // elements per thread
    Int add_key; add_key.SetInt64(keysPerThread);
    Point Add_Point = secp256k1->ScalarMultiplication(&add_key); // helper point to calculate the starting points
    
    Point addPoints[POINTS_BATCH_SIZE]; // array for the batch addition points(1G .. 1024G)
    Point batch_Add = secp256k1->DoublePoint(secp256k1->G); // 2G
    addPoints[0] = secp256k1->G; // 1G
    addPoints[1] = batch_Add;    // 2G
    for (int i = 2; i < POINTS_BATCH_SIZE; i++) // filling in the batch addition points array with points from(3G .. 1024G)
    {
        batch_Add = secp256k1->AddPoints(batch_Add, secp256k1->G);
        addPoints[i] = batch_Add;
    }
    
    int nbBatch = keysPerThread / POINTS_BATCH_SIZE; // number of batches for the single thread

    uint64_t bloom_size = n_elements * 4;
    uint64_t bloom_mod = bloom_size * 8;
    int iterations = 4;
    
    print_time(); cout << "Bloomfilter Size : " << bytesToSize(double(bloom_size), 2) << endl;

    auto bloom_create = [&]() {
        
        char * bloomfile = "bloom.bf";
        Point P = secp256k1->SubtractPoints(target_point, secp256k1->G);
        Point starting_points[n_cores];
        for (int i = 0; i < n_cores; i++) { // calculating the starting points 
            starting_points[i] = P;
            P = secp256k1->AddPoints(P, Add_Point);
        }

        unsigned char * bloom = (unsigned char *)calloc(bloom_size, sizeof(unsigned char));
        
        auto process_chunk = [&](Point start_point) { // function for a thread
            
            IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            Int deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            modGroup.Set(deltaX); // assign array deltaX to modGroup for batch inversion
            Int pointBatchX[POINTS_BATCH_SIZE]; // X coordinates of the batch
            Int pointBatchY[POINTS_BATCH_SIZE]; // Y coordinates of the batch
            Int deltaY, slope; // values to store the results of points addition formula
                      
            Point startPoint = start_point; // start point
            
            for (int s = 0; s < nbBatch; s++) {
                
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                    deltaX[i].ModSub(&startPoint.x, &addPoints[i].x); // insert each entry into the deltaX array
                }
    
                modGroup.ModInv();    // doing batch inversion
                
                int i;
                for (i = 0; i < POINTS_BATCH_SIZE - 1; i++) { // follow points addition formula logic
                    
                    deltaY.ModSub(&startPoint.y, &addPoints[i].y);
                    slope.ModMulK1(&deltaY, &deltaX[i]); // deltaX already inverted for each entry of the batch

                    pointBatchX[i].ModSquareK1(&slope); // calculating just X coordinate
                    pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                    pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);

                }
                
                deltaY.ModSub(&startPoint.y, &addPoints[i].y); // calculating X,Y coordinates for the last of the batch entry (used as the next startPoint)
                slope.ModMulK1(&deltaY, &deltaX[i]);

                pointBatchX[i].ModSquareK1(&slope);
                pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                    
                pointBatchY[i].ModSub(&startPoint.x, &pointBatchX[i]);
                pointBatchY[i].ModMulK1(&slope, &pointBatchY[i]);
                pointBatchY[i].ModSub(&pointBatchY[i], &startPoint.y);

                startPoint.x.Set(&pointBatchX[POINTS_BATCH_SIZE - 1]); // setting the new startPoint for the next batch iteration
                startPoint.y.Set(&pointBatchY[POINTS_BATCH_SIZE - 1]);
                
                mtx_lock.lock();

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // inserting all batch points X coordinates into the bloomfilter
                    for(int a = 0; a < iterations; a++) {
                        set_bit(bloom, pointBatchX[i].bits64[a] & (bloom_mod - 1));
                    }
                }

                mtx_lock.unlock();

            }      
        };
        
        std::thread myThreads[n_cores]; // launching threads
        for (int i = 0; i < n_cores; i++) {
            myThreads[i] = std::thread(process_chunk, starting_points[i]);
        }

        for (int i = 0; i < n_cores; i++) {
            myThreads[i].join(); // waiting for threads to finish
        }

        print_time(); cout << "Writing BloomFilter to bloom.bf" << endl; 
        save_bloom_filter(bloomfile, bloom, bloom_size);
    };        

    std::thread bloom_thread(bloom_create);
    
    print_time(); cout << "Creating bloomfilter image" << '\n';
    
    bloom_thread.join();
    
    print_elapsed_time(chrono_start);

}
