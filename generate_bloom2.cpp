#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>

#include "secp256k1/SECP256k1.h"
#include "secp256k1/Int.h"
#include "secp256k1/IntGroup.h"
#include "util/util.h"

using namespace std;

#define CPU_GRP_SIZE 1024

const int nbCPUThread = 2;

std::mutex lock1;
std::mutex lock2;

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
    Q = secp256k1->AddPoints(puzzle_point_divide2, Q);                 // 644 + 520 = 1164 (point 1164)  644.5 + 521 = 1165.5
                                                                       // if the puzzle point is in the lower range half the calculated point will behind it 1164 - > 1288(addition_search)
    string sp{secp256k1->GetPublicKeyHex(Q)};                         // if the puzzle point is in the higher range half 1784 the calculated point will be ahead of it
                                                                       // 1784 < - 1908 (subtraction_search)
    ofstream outFile1;                        // writing settings to files
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
    
    uint64_t bsSize = pow(2, block_width);  // number of elements == 2^block_width
    uint64_t kPerThread = bsSize / nbCPUThread; // elements per thread
    uint64_t nbStep = kPerThread / CPU_GRP_SIZE;
    
    Int add_key; add_key.SetInt64(kPerThread);
    Point Add_Point = secp256k1->ScalarMultiplication(&add_key); // helper point to calculate the starting points

    Point Gn[CPU_GRP_SIZE / 2];
    Point _2Gn;
    
    Point g = secp256k1->G;
    Gn[0] = g;
    g = secp256k1->DoublePoint(g);
    Gn[1] = g;
    for (int i = 2; i < CPU_GRP_SIZE / 2; i++) {
        g = secp256k1->AddPoints(g, secp256k1->G);
        Gn[i] = g;
    }
    
    _2Gn = secp256k1->DoublePoint(Gn[CPU_GRP_SIZE / 2 - 1]);
    
    uint64_t bloom_size = bsSize * 4;
    uint64_t bloom_mod = bloom_size * 8;
    int iterations = 4;
    
    print_time(); cout << "Bloomfilter Size : " << bytesToSize(double(bloom_size), 2);
    cout << " Total(2 blooms): " << bytesToSize(double(bloom_size * 2), 2) << endl;
    
    auto bloom_create1 = [&]() {
        
        char * bloomfile = "bloom1.bf"; // bloomfilter for even case (1164->1288)   [1288 + block_width] will hit
        Point P(puzzle_point);                                   //  (1165.5->1289) [1289 + block_width] will not hit
        Point starting_points[nbCPUThread];
        for (int i = 0; i < nbCPUThread; i++) { // calculating the starting points 
            starting_points[i] = P;
            P = secp256k1->AddPoints(P, Add_Point);
        }

        unsigned char * bloom = (unsigned char *)calloc(bloom_size, sizeof(unsigned char));
        
        auto process_chunk = [&](Point start_point) { // function for a thread
            
            Point startP = secp256k1->AddPoints(start_point, Gn[CPU_GRP_SIZE / 2 - 1]);
            
            IntGroup grp(CPU_GRP_SIZE / 2 + 1);
            Int dx[CPU_GRP_SIZE / 2 + 1];
            grp.Set(dx);
            
            Int dy, dyn, _s, _p;
            Point pts[CPU_GRP_SIZE], pp, pn;
        
            for (uint64_t s = 0; s < nbStep; s++) {
        
                // Fill group
                int i;
                int hLength = (CPU_GRP_SIZE / 2 - 1);
        
                for(i = 0; i < hLength; i++) { //fill dx for batch inversion
                    dx[i].ModSub(&Gn[i].x, &startP.x);
                }
                dx[i].ModSub(&Gn[i].x, &startP.x);  // For the first point
                dx[i + 1].ModSub(&_2Gn.x, &startP.x); // For the next center point

                grp.ModInv(); // Grouped ModInv
        
                // We use the fact that P + i*G and P - i*G has the same deltax, so the same inverse
                // We compute key in the positive and negative way from the center of the group
                // center point
                pts[CPU_GRP_SIZE / 2] = startP;
        
                for(i = 0; i < hLength; i++) {
        
                    pp = startP;
                    pn = startP;
        
                    // P = startP + i*G
                    dy.ModSub(&Gn[i].y, &pp.y);
                    _s.ModMulK1(&dy, &dx[i]);        // s = (p2.y-p1.y)*inverse(p2.x-p1.x);
                    
                    _p.ModSquareK1(&_s);            // _p = pow2(s)
                    pp.x.ModSub(&_p, &pp.x);
                    pp.x.ModSub(&Gn[i].x);
                
                    // P = startP - i*G  , if (x,y) = i*G then (x,-y) = -i*G
                    dyn.Set(&Gn[i].y);
                    dyn.ModNeg();
                    dyn.ModSub(&pn.y);
                    _s.ModMulK1(&dyn, &dx[i]);      // s = (p2.y-p1.y)*inverse(p2.x-p1.x);
                    
                    _p.ModSquareK1(&_s);            // _p = pow2(s)
                    pn.x.ModSub(&_p, &pn.x);
                    pn.x.ModSub(&Gn[i].x);
                    
                    // setting to pts
                    pts[CPU_GRP_SIZE / 2 + (i + 1)] = pp;
                    pts[CPU_GRP_SIZE / 2 - (i + 1)] = pn;
                
                }
        
                // First point (startP - (CPU_GRP_SIZE/2)*G)
                pn = startP;
                dyn.Set(&Gn[i].y);
                dyn.ModNeg();
                dyn.ModSub(&pn.y);
                _s.ModMulK1(&dyn, &dx[i]);
                
                _p.ModSquareK1(&_s);
                pn.x.ModSub(&_p, &pn.x);
                pn.x.ModSub(&Gn[i].x);
            
                pts[0] = pn;
            
                // Add to bloomfilter
                lock1.lock();
                for (int i = 0; i < CPU_GRP_SIZE; i++) {
                    for(int a = 0; a < iterations; a++) {
                        set_bit(bloom, pts[i].x.bits64[a] & (bloom_mod - 1));
                    }
                }
                lock1.unlock();
            
                // Next start point (startP + CPU_GRP_SIZE*G)
                pp = startP;
                
                dy.ModSub(&_2Gn.y, &pp.y);
                _s.ModMulK1(&dy, &dx[i + 1]);
                
                _p.ModSquareK1(&_s);
                pp.x.ModSub(&_p, &pp.x);
                pp.x.ModSub(&_2Gn.x);
            
                pp.y.ModSub(&_2Gn.x, &pp.x);
                pp.y.ModMulK1(&_s);
                pp.y.ModSub(&_2Gn.y);
                
                startP = pp;
        
            }
            
        };
        
        std::thread myThreads[nbCPUThread]; // launching threads
        for (int i = 0; i < nbCPUThread; i++) {
            myThreads[i] = std::thread(process_chunk, starting_points[i]);
        }

        for (int i = 0; i < nbCPUThread; i++) {
            myThreads[i].join(); // waiting for threads to finish
        }
         
        print_time(); cout << "Writing BloomFilter to bloom1.bf" << endl; 
        save_bloom_filter(bloomfile, bloom, bloom_size);
    };
    
    auto bloom_create2 = [&]() {
        
        char * bloomfile = "bloom2.bf"; // bloomfilter for odd case  (1165.5->1289.5) [1289.5 + block_width] will hit
        Point P(puzzle_point_05);                                 // (1164->1288.5)   [1288.5 + block_width] will not hit
        Point starting_points[nbCPUThread];
        for (int i = 0; i < nbCPUThread; i++) { // calculating the starting points 
            starting_points[i] = P;
            P = secp256k1->AddPoints(P, Add_Point);
        }
        
        unsigned char * bloom = (unsigned char *)calloc(bloom_size, sizeof(unsigned char));
        
        auto process_chunk = [&](Point start_point) {  // function for a thread
            
            Point startP = secp256k1->AddPoints(start_point, Gn[CPU_GRP_SIZE / 2 - 1]);
            
            IntGroup grp(CPU_GRP_SIZE / 2 + 1);
            Int dx[CPU_GRP_SIZE / 2 + 1];
            grp.Set(dx);
            
            Int dy, dyn, _s, _p;
            Point pts[CPU_GRP_SIZE], pp, pn;
        
            for (uint64_t s = 0; s < nbStep; s++) {
        
                // Fill group
                int i;
                int hLength = (CPU_GRP_SIZE / 2 - 1);
        
                for(i = 0; i < hLength; i++) { //fill dx for batch inversion
                    dx[i].ModSub(&Gn[i].x, &startP.x);
                }
                dx[i].ModSub(&Gn[i].x, &startP.x);  // For the first point
                dx[i + 1].ModSub(&_2Gn.x, &startP.x); // For the next center point
        
                grp.ModInv(); // Grouped ModInv
        
                // We use the fact that P + i*G and P - i*G has the same deltax, so the same inverse
                // We compute key in the positive and negative way from the center of the group
                // center point
                pts[CPU_GRP_SIZE / 2] = startP;
        
                for(i = 0; i < hLength; i++) {
        
                    pp = startP;
                    pn = startP;
                    
                    // P = startP + i*G
                    dy.ModSub(&Gn[i].y, &pp.y);
                    _s.ModMulK1(&dy, &dx[i]);        // s = (p2.y-p1.y)*inverse(p2.x-p1.x);
                    
                    _p.ModSquareK1(&_s);            // _p = pow2(s)
                    pp.x.ModSub(&_p, &pp.x);
                    pp.x.ModSub(&Gn[i].x);
                
                    // P = startP - i*G  , if (x,y) = i*G then (x,-y) = -i*G
                    dyn.Set(&Gn[i].y);
                    dyn.ModNeg();
                    dyn.ModSub(&pn.y);
                    _s.ModMulK1(&dyn, &dx[i]);      // s = (p2.y-p1.y)*inverse(p2.x-p1.x);
                    
                    _p.ModSquareK1(&_s);            // _p = pow2(s)
                    pn.x.ModSub(&_p, &pn.x);
                    pn.x.ModSub(&Gn[i].x);
                    
                    // setting to pts
                    pts[CPU_GRP_SIZE / 2 + (i + 1)] = pp;
                    pts[CPU_GRP_SIZE / 2 - (i + 1)] = pn;
                
                }
        
                // First point (startP - (CPU_GRP_SIZE/2)*G)
                pn = startP;
                
                dyn.Set(&Gn[i].y);
                dyn.ModNeg();
                dyn.ModSub(&pn.y);
                _s.ModMulK1(&dyn, &dx[i]);
                
                _p.ModSquareK1(&_s);
                pn.x.ModSub(&_p, &pn.x);
                pn.x.ModSub(&Gn[i].x);
            
                pts[0] = pn;
            
                // Add to bloomfilter
                lock2.lock();
                for (int i = 0; i < CPU_GRP_SIZE; i++) {
                    for(int a = 0; a < iterations; a++) {
                        set_bit(bloom, pts[i].x.bits64[a] & (bloom_mod - 1));
                    }
                }
                lock2.unlock();
            
                // Next start point (startP + CPU_GRP_SIZE*G)
                pp = startP;
                
                dy.ModSub(&_2Gn.y, &pp.y);
                _s.ModMulK1(&dy, &dx[i + 1]);
                
                _p.ModSquareK1(&_s);
                pp.x.ModSub(&_p, &pp.x);
                pp.x.ModSub(&_2Gn.x);
            
                pp.y.ModSub(&_2Gn.x, &pp.x);
                pp.y.ModMulK1(&_s);
                pp.y.ModSub(&_2Gn.y);
                
                startP = pp;
        
            }
            
        };
        
        std::thread myThreads[nbCPUThread]; // launching threads
        for (int i = 0; i < nbCPUThread; i++) {
            myThreads[i] = std::thread(process_chunk, starting_points[i]);
        }

        for (int i = 0; i < nbCPUThread; i++) {
            myThreads[i].join(); // waiting for threads to finish
        }

        print_time(); cout << "Writing BloomFilter to bloom2.bf" << endl; 
        save_bloom_filter(bloomfile, bloom, bloom_size);
    };
    
    std::thread thread1(bloom_create1);
    std::thread thread2(bloom_create2);
    
    print_time(); cout << "Creating bloomfilter images" << '\n';
    
    thread1.join();
    thread2.join();
    
    print_elapsed_time(chrono_start);
    
}
