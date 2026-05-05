#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <omp.h>

#include "secp256k1/SECP256k1.h"
#include "secp256k1/Int.h"
#include "secp256k1/IntGroup.h"
#include "bloom/filter.hpp"
#include "util/util.h"

using namespace std;
using filter = boost::bloom::filter<std::string, 32>; // bloomfilter settings

const double error = 0.0000000001; // errror rate for bloomfilter
//const int n_cores = 2; //actual number of processing cores equal to some power of two value (2,4,8,16,32,64,...)
//const int POINTS_BATCH_SIZE = 1024; // Batch addition with batch inversion using IntGroup class

auto main() -> int {
    
    auto chrono_start = std::chrono::high_resolution_clock::now();    // starting the timer
    Secp256K1* secp256k1 = new Secp256K1(); secp256k1->Init(); // initializing secp256k1 context
    
    std::remove("points.txt"); // remove previous settings and bloom files
    std::remove("points_match.txt");
    std::remove("bloomfile_points.bf");
    std::remove("search_settings.txt");
    
    Int pk; pk.SetInt32(1); // generating power of two values (2^0..2^256) table
    uint64_t mult = 2;
    vector<Int> S_table;
    for (int i = 0; i < 256; i++)
    {
        S_table.push_back(pk);
        pk.Mult(mult);
    }
    print_time(); cout << "S_table generated" << endl;
 
    pk.SetInt32(1); // generating power of two points table (2^0..2^256)
    vector<Point> P_table;
    Point P;
    for (int i = 0; i < 256; i++)
    {
        P = secp256k1->ScalarMultiplication(&pk);
        P_table.push_back(P);
        pk.Mult(mult);
    }
    print_time(); cout << "P_table generated" << endl;

    uint64_t range_start, range_end, block_width, search_power; // block_width = number of elements in the bloomfilter and a stride size to walk the range
    string temp, search_pub;
    ifstream inFile("settings.txt");
    getline(inFile, temp); range_start = std::stoull(temp);
    getline(inFile, temp); range_end = std::stoull(temp);
    getline(inFile, temp); block_width = std::stoull(temp);
    getline(inFile, temp); search_pub = trim(temp);
    inFile.close();
    
    search_power = block_width + 1;
    
    print_time(); cout << "Range Start: " << range_start << " bits" << endl;
    print_time(); cout << "Range End  : " << range_end << " bits" << endl;
    print_time(); cout << "Block Width: 2^" << block_width << endl;
    print_time(); cout << "Search Pub : " << search_pub << endl;

    Point start_point = P_table[range_start];
    Point end_point   = P_table[range_end];
    Point inverse_find_point = secp256k1->AddPoints(start_point, end_point);
    
    Point puzzle_point = secp256k1->ParsePublicKeyHex(search_pub);
    Point puzzle_inverse_point = secp256k1->SubtractPoints(inverse_find_point, puzzle_point);
    
    vector<Point> points;
    points.push_back(puzzle_point);
    points.push_back(puzzle_inverse_point);
    
    vector<Point> subtracted_points;
    
    Point A1, A2, S1, S2;
    ofstream WriteFile1("points_match.txt", std::ios_base::app);
    for(uint64_t i = search_power; i < range_start; i++) {
        A1 = secp256k1->AddPoints(puzzle_point, P_table[i]);
        A2 = secp256k1->AddPoints(puzzle_inverse_point, P_table[i]);
        S1 = secp256k1->SubtractPoints(puzzle_point, P_table[i]);
        S2 = secp256k1->SubtractPoints(puzzle_inverse_point, P_table[i]);
        points.push_back(A1);
        points.push_back(A2);
        points.push_back(S1);
        points.push_back(S2);
        WriteFile1 << i << " A : " << secp256k1->GetPublicKeyHex(A1) << " " << secp256k1->GetPublicKeyHex(A2) << '\n';
        WriteFile1 << i << " S : " << secp256k1->GetPublicKeyHex(S1) << " " << secp256k1->GetPublicKeyHex(S2) << '\n';
    }
    WriteFile1.close();
    
    ofstream WriteFile2("points.txt", std::ios_base::app);
    for (auto& point : points) {
         WriteFile2 << secp256k1->GetPublicKeyHex(point) << '\n';
    }
    WriteFile2.close();
    
    Int pkey(uint64_t(pow(2, block_width - 1)));
    Point sub_point = secp256k1->ScalarMultiplication(&pkey);
    for (auto& point : points) {
        P = secp256k1->SubtractPoints(point, sub_point);
        subtracted_points.push_back(P);
    }
    
    ofstream WriteFile4("search_settings.txt", std::ios_base::app);
    WriteFile4 << S_table[range_end].GetBase10() << '\n';
    WriteFile4 << S_table[range_start].GetBase10() << '\n';
    WriteFile4 << 2 << '\n';
    WriteFile4.close();
    
    print_time(); cout << "All necessary files created and written" << endl;
    
    print_time(); cout << "Creating bloomfilter image" << '\n';
    uint64_t pow_exp = uint64_t(pow(2, block_width)) + 1;
    uint64_t n_elements = subtracted_points.size() * pow_exp;
    filter bf(n_elements, error);
    
    for (auto& point : subtracted_points) {
        P = point;
        for(uint64_t i = 0; i < pow_exp; i++) {
            P = secp256k1->AddPoints(P, secp256k1->G);
            bf.insert(secp256k1->GetPublicKeyHex(P));
        }
    }
    
    string bloomfile = "bloomfile_points.bf";
    print_time(); cout << "Writing image to bloomfile_points.bf" << '\n';
    std::ofstream out(bloomfile, std::ios::binary);
    std::size_t c1= bf.capacity();
    out.write((const char*) &c1, sizeof(c1)); // save capacity (bits)
    boost::span<const unsigned char> s1 = bf.array();
    out.write((const char*) s1.data(), s1.size()); // save array
    out.close();
    /*
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
    
    auto bloom_create = [&]() {

        omp_lock_t lock;
        omp_init_lock(&lock);
        
        string bloomfile = "bloomfile_points.bf"; //       bloomfilter for even case (1164->1288)   [1288 + block_width] will hit
        Point P = secp256k1->SubtractPoints(target_point, secp256k1->G); //(1165.5->1289) [1289 + block_width] will not hit
        Point starting_points[n_cores];
        for (int i = 0; i < n_cores; i++) { // calculating the starting points 
            starting_points[i] = P;
            P = secp256k1->AddPoints(P, Add_Point);
        }

        filter bf(n_elements, error);
        
        auto process_chunk = [&](Point start_point) { // function for a thread
            
            IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            Int deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            modGroup.Set(deltaX); // assign array deltaX to modGroup for batch inversion (JLP way one time for all)
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
                
                omp_set_lock(&lock);
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // inserting all batch points X coordinates into the bloomfilter
                    bf.insert(pointBatchX[i].bits64[3]);
                }
                omp_unset_lock(&lock);
                
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
        
        omp_destroy_lock(&lock);
        
        print_time(); cout << "Writing image to bloom.bf" << '\n';
        std::ofstream out(bloomfile, std::ios::binary);
        std::size_t c1= bf.capacity();
        out.write((const char*) &c1, sizeof(c1)); // save capacity (bits)
        boost::span<const unsigned char> s1 = bf.array();
        out.write((const char*) s1.data(), s1.size()); // save array
        out.close();
    };


    std::thread bloom_thread(bloom_create);
    
    print_time(); cout << "Creating bloomfilter image" << '\n';
    
    bloom_thread.join();
    */
    print_elapsed_time(chrono_start);
  
}
