#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <omp.h>
#include <format>

#include "secp256k1/SECP256k1.h"
#include "secp256k1/Int.h"
#include "secp256k1/IntGroup.h"
#include "bloom/filter.hpp"
#include "util/util.h"

using namespace std;
using filter = boost::bloom::filter<std::string, 32>; // bloomfilter settings

const double error = 0.0000000001; // errror rate for bloomfilter
const int n_cores = 2; //actual number of processing cores equal to some power of two value(2,4,8,16,32,64,...) divided by 2
//const int xC_len = 10; //X coordinate length to be inserted into the bloomfilter (should be the same for generate_bloom and point_search max=33(full length X coordinate))

static constexpr int POINTS_BATCH_SIZE = 1024; // Batch addition with batch inversion using IntGroup class

static omp_lock_t lock1;
static omp_lock_t lock2;

auto main() -> int {
    
    auto chrono_start = std::chrono::high_resolution_clock::now();    // starting the timer
    Secp256K1* secp256k1 = new Secp256K1(); secp256k1->Init(); // initializing secp256k1 context
   
    std::remove("settings1.txt"); // remove previous settings and bloom files
    std::remove("settings2.txt");
    std::remove("bloom1.bf");
    std::remove("bloom2.bf");
    
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

    Point point_05, puzzle_point, puzzle_point_05, puzzle_point_divide2; // calculation points
    Point first_point, second_point, P1, P2, Q;                     // take 2^10 .. 2^11 range for example
    Int d_05; d_05.SetBase10("57896044618658097711785492504343953926418782139537452191302581570759080747169");
    
    point_05 = secp256k1->ScalarMultiplication(&d_05); // multiplicative inverse of 2 mod N (0.5)
    puzzle_point = secp256k1->ParsePublicKeyHex(search_pub); // 1288 1289 points for example
    puzzle_point_05 = secp256k1->AddPoints(puzzle_point, point_05); // 1288 + 0.5 = 1288.5
                                                                    // 1289 + 0.5 = 1289.5
    puzzle_point_divide2 = secp256k1->PointMultiplication(puzzle_point, &d_05); // 1288 : 2 = 644
                                                                                // 1289 : 2 = 644.5
    first_point  = P_table[range_start - 1]; // 512
    second_point = P_table[range_start - 2]; // 256

    P1 = secp256k1->SubtractPoints(puzzle_point_divide2, first_point); // 644 - 512 = 132  644.5 - 512 = 132.5
    P2 = secp256k1->SubtractPoints(puzzle_point_divide2, second_point);// 644 - 256 = 388  644.5 - 256 = 388.5
    Q = secp256k1->AddPoints(P1, P2);                                  // 132 + 388 = 520  132.5 + 388.5 = 521
    Q = secp256k1->AddPoints(Q, puzzle_point_divide2);                 // 644 + 520 = 1164 (point 1164)  644.5 + 521 = 1165.5
                                                                       // if the puzzle point is in the lower range half the calculated point will behind it 1164 - > 1288(addition_search)
    string sp{secp256k1->GetPublicKeyHex(Q)};                          // if the puzzle point is in higher range half
                                                                       // the calculated point will be ahead of it// 1784 < - 1908 (subtraction_search)
    ofstream outFile1; // writing settings to files
    outFile1.open("settings1.txt", ios::app);                          // for odd number 1289 the point will be 1165.5
    outFile1 << sp <<'\n';                                             // that is why two bloomfilters are used
    outFile1 << "0" << '\n';
    outFile1.close();
    
    ofstream outFile2;
    outFile2.open("settings2.txt", ios::app);
    outFile2 << sp <<'\n';
    outFile2 << "0" << '\n';
    outFile2.close();
    
    print_time(); cout << "Settings written to file" << endl;
    
    uint64_t n_elements = pow(2, block_width);  // number of elements == 2^block_width
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
    
    omp_init_lock(&lock1);
    omp_init_lock(&lock2);
    
    auto bloom_create1 = [&]() {
        
        string bloomfile = "bloom1.bf"; //       bloomfilter for even case (1164->1288)   [1288 + block_width] will hit
        Point P = secp256k1->SubtractPoints(puzzle_point, secp256k1->G); //(1165.5->1289) [1289 + block_width] will not hit
        Point starting_points[n_cores];
        for (int i = 0; i < n_cores; i++) { // calculating the starting points 
            starting_points[i] = P;
            P = secp256k1->AddPoints(P, Add_Point);
        }

        filter bf(n_elements, error);
        
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
                
                omp_set_lock(&lock1);
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // inserting all batch points X coordinates into the bloomfilter
                    bf.insert(std::format("{:x}", pointBatchX[i].bits64[3]));
                    //bf.insert(secp256k1->GetXHex(&pointBatchX[i], xC_len)); // xC_len (length of the X coordinate hex representation)
                }
                omp_unset_lock(&lock1);
                
                startPoint.x.Set(&pointBatchX[POINTS_BATCH_SIZE - 1]); // setting the new startPoint for the next batch iteration
                startPoint.y.Set(&pointBatchY[POINTS_BATCH_SIZE - 1]);

            }
            
        };
        
        std::thread myThreads[n_cores]; // launching threads
        for (int i = 0; i < n_cores; i++) {
            myThreads[i] = std::thread(process_chunk, starting_points[i]);
        }

        for (int i = 0; i < n_cores; i++) {
            myThreads[i].join(); // waiting for threads to finish
        }
        
        omp_destroy_lock(&lock1);
        
        print_time(); cout << "Writing bloom1 image to bloom1.bf" << '\n';
        std::ofstream out(bloomfile, std::ios::binary);
        std::size_t c1= bf.capacity();
        out.write((const char*) &c1, sizeof(c1)); // save capacity (bits)
        boost::span<const unsigned char> s1 = bf.array();
        out.write((const char*) s1.data(), s1.size()); // save array
        out.close();
    };

    auto bloom_create2 = [&]() {
        
        string bloomfile = "bloom2.bf"; //           bloomfilter for odd case  (1165.5->1289.5) [1289.5 + block_width] will hit
        Point P = secp256k1->SubtractPoints(puzzle_point_05, secp256k1->G); // (1164->1288.5)   [1288.5 + block_width] will not hit
        Point starting_points[n_cores];
        for (int i = 0; i < n_cores; i++) { // calculating the starting points 
            starting_points[i] = P;
            P = secp256k1->AddPoints(P, Add_Point);
        }
        
        filter bf(n_elements, error);
        
        auto process_chunk = [&](Point start_point) {  // function for a thread
            
            IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            Int deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            modGroup.Set(deltaX); // assign array deltaX to modGroup for batch inversion
            Int pointBatchX[POINTS_BATCH_SIZE]; // X coordinates of the batch
            Int pointBatchY[POINTS_BATCH_SIZE]; // Y coordinates of the batch
            Int deltaY, slope; // values to store the results of points addition formula
                      
            Point startPoint = start_point;  // start point
            
            for (int s = 0; s < nbBatch; s++) {
                
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) {         // // we compute (x1 - x2) for each entry of the entire batch
                    deltaX[i].ModSub(&startPoint.x, &addPoints[i].x); // insert each entry into the deltaX array
                }
    
                modGroup.ModInv();    // doing batch inversion
                
                int i;
                for (i = 0; i < POINTS_BATCH_SIZE - 1; i++) { // follow points addition formula logic

                    deltaY.ModSub(&startPoint.y, &addPoints[i].y); // deltaX already inverted for each entry of the batch
                    slope.ModMulK1(&deltaY, &deltaX[i]);
                    
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
                
                omp_set_lock(&lock2);
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // inserting all batch points X coordinates into the bloomfilter
                    bf.insert(std::format("{:x}", pointBatchX[i].bits64[3]));
                    //bf.insert(secp256k1->GetXHex(&pointBatchX[i], xC_len)); // xC_len (length of the X coordinate hex representation)
                }
                omp_unset_lock(&lock2);
                
                startPoint.x.Set(&pointBatchX[POINTS_BATCH_SIZE - 1]); // setting the new startPoint for the next batch iteration
                startPoint.y.Set(&pointBatchY[POINTS_BATCH_SIZE - 1]);
                
            }
            
        };
        
        std::thread myThreads[n_cores]; // launching threads
        for (int i = 0; i < n_cores; i++) {
            myThreads[i] = std::thread(process_chunk, starting_points[i]);
        }
        
        for (int i = 0; i < n_cores; i++) {
            myThreads[i].join(); // waiting for threads to finish
        }
        
        omp_destroy_lock(&lock2);

        print_time(); cout << "Writing bloom2 image to bloom2.bf" << '\n'; 
        std::ofstream out(bloomfile, std::ios::binary);
        std::size_t c1= bf.capacity();
        out.write((const char*) &c1, sizeof(c1)); // save capacity (bits)
        boost::span<const unsigned char> s1 = bf.array();
        out.write((const char*) s1.data(), s1.size()); // save array
        out.close();
    };

    std::thread thread1(bloom_create1);
    std::thread thread2(bloom_create2);
    
    print_time(); cout << "Creating bloomfilter images" << '\n';
    
    thread1.join();
    thread2.join();
    
    print_elapsed_time(chrono_start);

}
