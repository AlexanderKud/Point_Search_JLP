#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <thread>
#include <omp.h>

#include "secp256k1/SECP256k1.h"
#include "secp256k1/Int.h"
#include "secp256k1/IntGroup.h"
#include "bloom/filter.hpp"
#include "util/util.h"

using namespace std;
namespace fs = filesystem;

static constexpr int POINTS_BATCH_SIZE = 1024; // Batch addition with batch inversion using IntGroup class

static omp_lock_t lock1;
static omp_lock_t lock2;

auto main() -> int {
    
    auto start = std::chrono::high_resolution_clock::now();    // starting the timer
    Secp256K1* secp256k1 = new Secp256K1(); secp256k1->Init(); // initializing secp256k1 context
    int xC_len = 10; //X coordinate length to be inserted into the bloomfilter (should be the same for generate_bloom and point_search max=33(full length X coordinate))
    
    fs::path current_path = fs::current_path(); // deleting previous settings and bloom files
    auto file_list = get_files_in_directory(current_path);
    vector<string> targets = {"settings1.txt", "settings2.txt", "bloom1.bf", "bloom2.bf"};
    for (auto& i : file_list) {
        for (auto& t : targets) {
            if (i == t) { std::remove(t.c_str()); }
        }
    }

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
    getline(inFile, temp); range_start = str_to_uint64(temp);
    getline(inFile, temp); range_end = str_to_uint64(temp);
    getline(inFile, temp); block_width = str_to_uint64(temp);
    getline(inFile, temp); search_pub = trim(temp);
    inFile.close();
    
    print_time(); cout << "Range Start: " << range_start << " bits" << endl;
    print_time(); cout << "Range End  : " << range_end << " bits" << endl;
    print_time(); cout << "Block Width: 2^" << block_width << endl;
    print_time(); cout << "Search Pub : " << search_pub << endl;

    Point point_05, puzzle_point, puzzle_point_05, puzzle_point_divide2; // calculation points
    Point first_point, second_point, P1, P2, Q1, Q2;                     // take 2^10 .. 2^11 range for example
    Int stride_sum; stride_sum.SetInt32(0); // value to keep the travelled distance by stride
    Int d_05, div2;
    d_05.SetBase10("57896044618658097711785492504343953926418782139537452191302581570759080747169");
    div2.SetInt32(2);
    point_05 = secp256k1->ScalarMultiplication(&d_05); // multiplicative inverse of 2 mod N (0.5)
    puzzle_point = secp256k1->ParsePublicKeyHex(search_pub); // 1288 1289 points for example
    puzzle_point_05 = secp256k1->AddPoints(puzzle_point, point_05); // 1288 + 0.5 = 1288.5
                                                                    // 1289 + 0.5 = 1289.5
    puzzle_point_divide2 = secp256k1->PointDivision(puzzle_point, &div2); // 1288 : 2 = 644
                                                                          // 1289 : 2 = 644.5
    first_point  = P_table[range_start - 1]; // 512
    second_point = P_table[range_start - 2]; // 256

    P1 = secp256k1->SubtractPoints(puzzle_point_divide2, first_point); // 644 - 512 = 132  644.5 - 512 = 132.5
    P2 = secp256k1->SubtractPoints(puzzle_point_divide2, second_point);// 644 - 256 = 388  644.5 - 256 = 388.5
    Q1 = secp256k1->AddPoints(P1, P2);                                 // 132 + 388 = 520  132.5 + 388.5 = 521
    Q2 = secp256k1->AddPoints(puzzle_point_divide2, Q1);               // 644 + 520 = 1164 (point 1164)  644.5 + 521 = 1165.5
                                                                       // if the puzzle point is in the lower range half
    string s1, s2;                                                     // the calculated point will behind it 1164 - > 1288(addition_search)
    s1 = secp256k1->GetPublicKeyHex(Q2);                               // if the puzzle point is in the higher range half 1784
    s2 = stride_sum.GetBase10();                                       // the calculated point will be ahead of it
                                                                       // 1784 < - 1908 (subtraction_search)
    ofstream outFile1;                        // writing settings to files
    outFile1.open("settings1.txt", ios::app);                          // for odd number 1289 the point will be 1165.5
    outFile1 << s1 <<'\n';                                             // that is why two bloomfilters are used
    outFile1 << s2 << '\n';
    outFile1.close();
    
    ofstream outFile2;
    outFile2.open("settings2.txt", ios::app);
    outFile2 << s1 <<'\n';
    outFile2 << s2 << '\n';
    outFile2.close();
    
    print_time(); cout << "Settings written to file" << endl;
    
    using filter = boost::bloom::filter<std::string, 32>; // bloomfilter settings
    uint64_t n_elements = uint64_t(pow(2, block_width));  // number of elements == 2^block_width
    double error = 0.0000000001; // errror rate for bloomfilter
    int n_cores = 4; //actual number of processing cores equal to some power of two value(2,4,8,16,32,64,...) divided by 2
    uint64_t count = uint64_t(pow(2, block_width) / n_cores); // elements per thread
    Int add_key; add_key.SetInt64(count);
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
    
    int nbBatch = count / POINTS_BATCH_SIZE; // number of batches for the single thread
    
    omp_init_lock(&lock1);
    omp_init_lock(&lock2);
    
    auto bloom_create1 = [&]() {
        
        string bloomfile = "bloom1.bf"; //       bloomfilter for even case (1164->1288)   [1288 + block_width] will hit
        Point P = secp256k1->SubtractPoints(puzzle_point, secp256k1->G); //(1165.5->1289) [1289 + block_width] will not hit
        vector<Point> starting_points;
        for (int i = 0; i < n_cores; i++) { // calculating the starting points 
            starting_points.push_back(P);
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

                    pointBatchX[i].ModSquareK1(&slope);
                    pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                    pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);

                }
                
                deltaY.ModSub(&startPoint.y, &addPoints[i].y);
                slope.ModMulK1(&deltaY, &deltaX[i]);

                pointBatchX[i].ModSquareK1(&slope);
                pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                    
                pointBatchY[i].ModSub(&startPoint.x, &pointBatchX[i]);
                pointBatchY[i].ModMulK1(&slope, &pointBatchY[i]);
                pointBatchY[i].ModSub(&pointBatchY[i], &startPoint.y);
                
                omp_set_lock(&lock1);
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // inserting all batch points into the bloomfilter
                    bf.insert(secp256k1->GetXHex(&pointBatchX[i], xC_len));
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
    
        print_time(); cout << "Creating bloom1 image with " << n_cores << " threads" << '\n';

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
        vector<Point> starting_points;
        for (int i = 0; i < n_cores; i++) { // calculating the starting points 
            starting_points.push_back(P);
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

                    deltaY.ModSub(&startPoint.y, &addPoints[i].y);
                    slope.ModMulK1(&deltaY, &deltaX[i]);
                    
                    pointBatchX[i].ModSquareK1(&slope);
                    pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                    pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);

                }
                
                deltaY.ModSub(&startPoint.y, &addPoints[i].y);
                slope.ModMulK1(&deltaY, &deltaX[i]);
                    
                pointBatchX[i].ModSquareK1(&slope);
                pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                    
                pointBatchY[i].ModSub(&startPoint.x, &pointBatchX[i]);
                pointBatchY[i].ModMulK1(&slope, &pointBatchY[i]);
                pointBatchY[i].ModSub(&pointBatchY[i], &startPoint.y);
                
                omp_set_lock(&lock2);
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // inserting all batch points into the bloomfilter
                    bf.insert(secp256k1->GetXHex(&pointBatchX[i], xC_len));
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
    
        print_time(); cout << "Creating bloom2 image with " << n_cores << " threads" << '\n';

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
    
    thread1.join();
    thread2.join();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = end - start;
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    duration -= hours;
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
    duration -= minutes;
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
    print_time(); cout << "Elapsed time: (" << hours.count() << ")hours (" << minutes.count() << ")minutes (" << seconds.count() << ")seconds\n";
}
