#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <vector>
#include <algorithm>
#include <thread>

#include "secp256k1/SECP256k1.h"
#include "secp256k1/Int.h"
#include "secp256k1/IntGroup.h"
#include "bloom/filter.hpp"
#include "util/util.h"

using namespace std;

static constexpr int POINTS_BATCH_SIZE = 1024; // Batch addition with batch inversion using IntGroup class

auto main() -> int {

    Secp256K1* secp256k1 = new Secp256K1(); secp256k1->Init(); // initialize secp256k1 context
    int cpuCores = 4; // actual number of processing cores divided by 2
    int xC_len = 15; // X coordinate length to be checked for being insererted into the bloomfilter (should be equal for generate_bloom and point_search max=33(full length X coordinate))
    
    Int pk; pk.SetInt32(1); // generating power of two values (2^0..2^256) table
    uint64_t mult = 2;
    vector<Int> S_table;
    for (int i = 0; i < 256; i++)
    {
        S_table.push_back(pk);
        pk.Mult(mult);
    }
    print_time(); cout << "S_table generated" << endl;

    uint64_t range_start, range_end, block_width; // block_width = number of elements in the bloomfilter and a stride size to walk the range
    string temp, search_pub;
    ifstream inFile("settings.txt"); // load setiings from file
    getline(inFile, temp); range_start = str_to_uint64(temp);
    getline(inFile, temp); range_end = str_to_uint64(temp);
    getline(inFile, temp); block_width = str_to_uint64(temp);
    getline(inFile, temp); search_pub = trim(temp);
    inFile.close();
    print_time(); cout << "Range Start: " << range_start << " bits" << endl;
    print_time(); cout << "Range End  : " << range_end << " bits" << endl;
    print_time(); cout << "Block Width: 2^" << block_width << endl;
    print_time(); cout << "Search Pub : " << search_pub << endl;

    Int pre_calc_sum; // precalculated sum for private key recovering 512 + 256 = 768 for range[2^10..2^11] 1288 in example
    pre_calc_sum.Add(&S_table[range_start - 1], &S_table[range_start - 2]);
    
    string bloomfile1 = "bloom1.bf"; // bloomfilter stuff
    string bloomfile2 = "bloom2.bf";
    using filter = boost::bloom::filter<std::string, 32>;
    
    print_time(); cout << "Loading Bloomfilter bloom1.bf" << endl;
    filter bf1;
    std::ifstream in1(bloomfile1, std::ios::binary);
    std::size_t c1;
    in1.read((char*) &c1, sizeof(c1));
    bf1.reset(c1); // restore capacity
    boost::span<unsigned char> s1 = bf1.array();
    in1.read((char*) s1.data(), s1.size()); // load array
    in1.close();

    
    print_time(); cout << "Loading Bloomfilter bloom2.bf" << endl;
    filter bf2;
    std::ifstream in2(bloomfile2, std::ios::binary);
    std::size_t c2;
    in2.read((char*) &c2, sizeof(c2));
    bf2.reset(c2); // restore capacity
    boost::span<unsigned char> s2 = bf2.array();
    in2.read((char*) s2.data(), s2.size()); // load array
    in2.close();
    
    auto pow10_nums = break_down_to_pow10(uint64_t(pow(2, block_width))); // decomposing the 2^block_width to the power of ten values
    vector<Point> pow10_points;                                           // to get the index of the bloomfilter element fast
    Int pow_key;
    for (auto& n : pow10_nums) { // calculating points corresponding to the decomposition components
        pow_key.SetInt64(n);     
        pow10_points.push_back(secp256k1->ScalarMultiplication(&pow_key));
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    auto addition_search = [&]() { // addition search for the case when the starting point is behind the target point after calculations
        uint64_t mult = 2;         // the closer the target point to the center of the range from either side
        int save_counter = 0;      // the faster collision will happen
        string temp, cpub;
        Point start_point, stride_point, calc_point;
        Int stride_sum, stride;
        ifstream inFile("settings1.txt");
        getline(inFile, temp);
        start_point = secp256k1->ParsePublicKeyHex(trim(temp));
        getline(inFile, temp);
        stride_sum.SetBase10(trim(temp).data());
        inFile.close();
        
        stride.SetInt64(uint64_t(pow(2, block_width)));
        stride_point = secp256k1->ScalarMultiplication(&stride);
        
        int n_cores = cpuCores; //start splitting the search initiative according to the chosen number of cpu cores
    
        Int offset_Step, int_Cores, vector_Num;
        int_Cores.SetInt32(n_cores);
        offset_Step.floor_Div(&S_table[range_start - 2], &int_Cores);
        
        vector_Num.SetInt32(0);
        vector<Int> offset_Nums;
        for (int i = 0; i < n_cores; i++) {
            offset_Nums.push_back(vector_Num);
            vector_Num.Add(&offset_Step);
        }
        
        vector<Point> offset_Points;
        offset_Points.push_back(secp256k1->G);
        for (int i = 1; i < n_cores; i++) {
            offset_Points.push_back(secp256k1->ScalarMultiplication(&offset_Nums[i]));
        }

        vector<Point> starting_Points;
        Point vector_Point(start_point);
        starting_Points.push_back(vector_Point);

        for (int i = 1; i < n_cores; i++) {
            vector_Point = secp256k1->AddPoints(start_point, offset_Points[i]);
            starting_Points.push_back(vector_Point);
        }

        Point addPoints[POINTS_BATCH_SIZE]; // array for batch addition points       
        Point batch_Add = secp256k1->DoublePoint(stride_point);
        addPoints[0] = stride_point;
        addPoints[1] = batch_Add;
        for (int i = 2; i < POINTS_BATCH_SIZE; i++) // filling in batch addition points array with points
        {
            batch_Add = secp256k1->AddPoints(batch_Add, stride_point);
            addPoints[i] = batch_Add;
        }
        // scalable lambda gets its chunk to search through
        auto scalable_addition_search = [&](Point starting_Point, int threadIdx, Int offset, Int stride_Sum) {

            Point P;
            Int stride_sum; stride_sum.Set(&stride_Sum);
            Int Int_steps, Int_temp, privkey;
            string cpub, xc, xc1, xc2;
            int index, count;
            uint64_t steps;
            vector<uint64_t> privkey_num;

            Int deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            Int pointBatchX[POINTS_BATCH_SIZE]; // X coordinates of the batch
            Int pointBatchY[POINTS_BATCH_SIZE]; // Y coordinates of the batch
                      
            Point startPoint = starting_Point; // start point

            Point BloomP; // point for insertion of the batch into the bloomfilter
            Int deltaY, slope, slopeSquared; // values to store the results of points addition formula

            Int batch_stride, batch_index;
            batch_stride.Mult(&stride, uint64_t(POINTS_BATCH_SIZE));
        
            while (true) {

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                    deltaX[i].ModSub(&startPoint.x, &addPoints[i].x); // insert each entry into the deltaX array
                }
    
                modGroup.Set(deltaX); // assign array deltaX to modGroup for batch inversion
                modGroup.ModInv();    // doing batch inversion
                
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // follow points addition formula logic
                    
                    deltaY.ModSub(&startPoint.y, &addPoints[i].y);
                    slope.ModMulK1(&deltaY, &deltaX[i]); // deltaX already inverted for each entry of the batch

                    slopeSquared.ModSquareK1(&slope);
                    pointBatchX[i].ModSub(&slopeSquared, &startPoint.x);
                    pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                    
                    pointBatchY[i].ModSub(&startPoint.x, &pointBatchX[i]);
                    pointBatchY[i].ModMulK1(&slope, &pointBatchY[i]);
                    pointBatchY[i].ModSub(&pointBatchY[i], &startPoint.y);
                }

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) {
                    
                    xc = secp256k1->GetXHex(&pointBatchX[i], xC_len);                  
                    
                    if (bf1.may_contain(xc)) {
                        
                        print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Lower Range Half]" << endl;
                        
                        BloomP.x.Set(&pointBatchX[i]);
                        BloomP.y.Set(&pointBatchY[i]);
                        BloomP.z.SetInt32(1);
                        
                        P = BloomP;
                        
                        privkey_num.clear();
                        index = 0;
                        for (auto& p : pow10_points) { // getting the index of the element in the bloomfilter
                            count = 0;
                            xc1 = secp256k1->GetXHex(&P.x, xC_len);
                            while (bf1.may_contain(xc1)) {
                                P = secp256k1->SubtractPoints(P, p);
                                xc1 = secp256k1->GetXHex(&P.x, xC_len);
                                count += 1;
                            }
                            privkey_num.push_back(pow10_nums[index] * (count - 1));
                            P = secp256k1->AddPoints(P, p);
                            index += 1;
                        }
                        steps = 0;
                        for (auto& n : privkey_num) { steps += n; } // we got here the index of the element in the bloomfilter
                        Int_steps.SetInt64(steps); // restoring the private key
                        batch_index.Mult(&stride, uint64_t(i + 1));
                        Int_temp.Add(&stride_sum, &batch_index);
                        Int_temp.Add(&offset);
                        Int_temp.Sub(&Int_steps);
                        privkey.Sub(&pre_calc_sum, &Int_temp);
                        privkey.Mult(mult); // we got here the private key
                        calc_point = secp256k1->ScalarMultiplication(&privkey);
                        if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) { // if cpubs are equal we got it
                            print_time(); cout << "Privatekey: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            auto end = std::chrono::high_resolution_clock::now();
                            auto duration = end - start;
                            auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
                            duration -= hours;
                            auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
                            duration -= minutes;
                            auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
                            print_time(); cout << "Elapsed time: (" << hours.count() << ")hours (" << minutes.count() << ")minutes (" << seconds.count() << ")seconds\n";
                            exit(0);
                        }
                        print_time(); cout << "False Positive" << endl;
                    }
                    
                    if (bf2.may_contain(xc)) {
                        
                        print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Lower Range Half]" << endl;
                        
                        BloomP.x.Set(&pointBatchX[i]);
                        BloomP.y.Set(&pointBatchY[i]);
                        BloomP.z.SetInt32(1);
                        
                        P = BloomP;
                        
                        privkey_num.clear();
                        index = 0;
                        for (auto& p : pow10_points) { // getting the index of the element in the bloomfilter
                            count = 0;
                            xc2 = secp256k1->GetXHex(&P.x, xC_len);
                            while (bf2.may_contain(xc2)) {
                                P = secp256k1->SubtractPoints(P, p);
                                xc2 = secp256k1->GetXHex(&P.x, xC_len);
                                count += 1;
                            }
                            privkey_num.push_back(pow10_nums[index] * (count - 1));
                            P = secp256k1->AddPoints(P, p);
                            index += 1;
                        }                   
                        steps = 0;
                        for (auto& n : privkey_num) { steps += n; } // we got here the index of the element in the bloomfilter
                        Int_steps.SetInt64(steps); // restoring the private key
                        batch_index.Mult(&stride, uint64_t(i + 1));
                        Int_temp.Add(&stride_sum, &batch_index);
                        Int_temp.Add(&offset);
                        Int_temp.Sub(&Int_steps);
                        privkey.Sub(&pre_calc_sum, &Int_temp);
                        privkey.Mult(mult);
                        privkey.AddOne(); // we got here the private key
                        calc_point = secp256k1->ScalarMultiplication(&privkey);
                        if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) { // if cpubs are equal we got it
                            print_time(); cout << "Privatekey: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            auto end = std::chrono::high_resolution_clock::now();
                            auto duration = end - start;
                            auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
                            duration -= hours;
                            auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
                            duration -= minutes;
                            auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
                            print_time(); cout << "Elapsed time: (" << hours.count() << ")hours (" << minutes.count() << ")minutes (" << seconds.count() << ")seconds\n";
                            exit(0);
                        }
                        print_time(); cout << "False Positive" << endl; 
                    }
                }
                
                startPoint.x.Set(&pointBatchX[POINTS_BATCH_SIZE - 1]); // setting the new startPoint for the next batch iteration
                startPoint.y.Set(&pointBatchY[POINTS_BATCH_SIZE - 1]);
                startPoint.z.SetInt32(1);
                
                stride_sum.Add(&batch_stride);
                    
                if (threadIdx == 0) {  // thread with index zero is used to save the progress
                    save_counter += 1; // all values are derived from this data after new program start
                    if (save_counter % 70000 == 0) {
                        cpub = secp256k1->GetPublicKeyHex(startPoint);
                        ofstream outFile;
                        outFile.open("settings1.txt");
                        outFile << cpub <<'\n';
                        outFile << stride_sum.GetBase10() << '\n';
                        outFile.close();
                        save_counter = 0;
                        print_time(); cout << "Save Data written to settings1.txt" << endl;
                    }
                }
            } // while (true) loop end curly brace
        };
        
        std::thread addition_Threads[n_cores];
        for (int i = 0; i < n_cores; i++) {
            addition_Threads[i] = std::thread(scalable_addition_search, starting_Points[i], i , offset_Nums[i], stride_sum);
        }

        for (int i = 0; i < n_cores; i++) {
            addition_Threads[i].join();
        }
    };
    
    auto subtraction_search = [&]() { // subtraction search for the case when the starting point is ahead of the target point after calculations
        uint64_t mult = 2;            // the closer the target point to the center of the range from either side
        int save_counter = 0;         // the faster collision will happen
        string temp;
        Point start_point, stride_point, calc_point;
        Int stride_sum, stride;
        ifstream inFile("settings2.txt");
        getline(inFile, temp);
        start_point = secp256k1->ParsePublicKeyHex(trim(temp));
        getline(inFile, temp);
        stride_sum.SetBase10(trim(temp).data());
        inFile.close();
        
        stride.SetInt64(uint64_t(pow(2, block_width)));
        stride_point = secp256k1->ScalarMultiplication(&stride);
        
        int n_cores = cpuCores; //start splitting the search according to the chosen number of cpu cores
        
        Int offset_Step, int_Cores, vector_Num;
        int_Cores.SetInt32(n_cores);
        offset_Step.floor_Div(&S_table[range_start - 2], &int_Cores);
        
        vector_Num.SetInt32(0);
        vector<Int> offset_Nums;
        for (int i = 0; i < n_cores; i++) {
            offset_Nums.push_back(vector_Num);
            vector_Num.Add(&offset_Step);
        }
        
        vector<Point> offset_Points;
        offset_Points.push_back(secp256k1->G);
        for (int i = 1; i < n_cores; i++) {
            offset_Points.push_back(secp256k1->ScalarMultiplication(&offset_Nums[i]));
        }

        vector<Point> starting_Points;
        Point vector_Point(start_point);
        starting_Points.push_back(vector_Point);

        for (int i = 1; i < n_cores; i++) {
            vector_Point = secp256k1->SubtractPoints(start_point, offset_Points[i]);
            starting_Points.push_back(vector_Point);
        }

        Point addPoints[POINTS_BATCH_SIZE]; // array for batch addition points       
        Point batch_Add = secp256k1->DoublePoint(stride_point);
        addPoints[0] = stride_point;
        addPoints[0].y.ModNeg();
        addPoints[1] = batch_Add;
        addPoints[1].y.ModNeg();
        for (int i = 2; i < POINTS_BATCH_SIZE; i++) // filling in batch addition points array with points
        {
            batch_Add = secp256k1->AddPoints(batch_Add, stride_point);
            addPoints[i] = batch_Add;
            addPoints[i].y.ModNeg();
        }
        // scalable lambda gets its chunk to search through
        auto scalable_subtraction_search = [&](Point starting_Point, int threadIdx, Int offset, Int stride_Sum) {

            Point P;
            Int stride_sum; stride_sum.Set(&stride_Sum);
            string cpub, xc, xc1, xc2;
            int index, count;
            uint64_t steps;
            vector<uint64_t> privkey_num;
            Int Int_steps, Int_temp, privkey;

            Int deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            Int pointBatchX[POINTS_BATCH_SIZE]; // X coordinates of the batch
            Int pointBatchY[POINTS_BATCH_SIZE]; // Y coordinates of the batch
                      
            Point startPoint = starting_Point; // start point

            Point BloomP; // point for insertion of the batch into the bloomfilter
            Int deltaY, slope, slopeSquared; // values to store the results of points addition formula
            
            Int batch_stride, batch_index;
            batch_stride.Mult(&stride, uint64_t(POINTS_BATCH_SIZE));
        
            while (true) {

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                    deltaX[i].ModSub(&startPoint.x, &addPoints[i].x); // insert each entry into the deltaX array
                }
    
                modGroup.Set(deltaX); // assign array deltaX to modGroup for batch inversion
                modGroup.ModInv();    // doing batch inversion
                
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // follow points addition formula logic
                    
                    deltaY.ModSub(&startPoint.y, &addPoints[i].y);
                    slope.ModMulK1(&deltaY, &deltaX[i]); // deltaX already inverted for each entry of the batch

                    slopeSquared.ModSquareK1(&slope);
                    pointBatchX[i].ModSub(&slopeSquared, &startPoint.x);
                    pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                    
                    pointBatchY[i].ModSub(&startPoint.x, &pointBatchX[i]);
                    pointBatchY[i].ModMulK1(&slope, &pointBatchY[i]);
                    pointBatchY[i].ModSub(&pointBatchY[i], &startPoint.y);
                }

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) {
                    
                    xc = secp256k1->GetXHex(&pointBatchX[i], xC_len);

                    if (bf1.may_contain(xc)) {
                        
                        print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Higher Range Half]" << endl;
                        
                        BloomP.x.Set(&pointBatchX[i]);
                        BloomP.y.Set(&pointBatchY[i]);
                        BloomP.z.SetInt32(1);
                        
                        P = BloomP;
                        
                        privkey_num.clear();
                        index = 0;
                        for (auto& p : pow10_points) { // getting the index of the element in the bloomfilter
                            count = 0;
                            xc1 = secp256k1->GetXHex(&P.x, xC_len);
                            while (bf1.may_contain(xc1)) {
                                P = secp256k1->SubtractPoints(P, p);
                                xc1 = secp256k1->GetXHex(&P.x, xC_len);
                                count += 1;
                            }
                            privkey_num.push_back(pow10_nums[index] * (count - 1));
                            P = secp256k1->AddPoints(P, p);
                            index += 1;
                        }                   
                        steps = 0;
                        for (auto& n : privkey_num) { steps += n; } // we got here the index of element in the bloomfilter
                        Int_steps.SetInt64(steps); // restoring the private key
                        batch_index.Mult(&stride, uint64_t(i + 1));
                        Int_temp.Add(&stride_sum, &batch_index);                        
                        Int_temp.Add(&offset);
                        Int_temp.Add(&Int_steps);
                        privkey.Add(&pre_calc_sum, &Int_temp);
                        privkey.Mult(mult); // we got here the private key
                        calc_point = secp256k1->ScalarMultiplication(&privkey);
                        if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) { // if cpubs are equal we got it
                            print_time(); cout << "Privatekey: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            auto end = std::chrono::high_resolution_clock::now();
                            auto duration = end - start;
                            auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
                            duration -= hours;
                            auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
                            duration -= minutes;
                            auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
                            print_time(); cout << "Elapsed time: (" << hours.count() << ")hours (" << minutes.count() << ")minutes (" << seconds.count() << ")seconds\n";
                            exit(0);
                        }
                        print_time(); cout << "False Positive" << endl;
                    }
                    
                    if (bf2.may_contain(xc)) {
                        
                        print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Higher Range Half]" << endl;
                        
                        BloomP.x.Set(&pointBatchX[i]);
                        BloomP.y.Set(&pointBatchY[i]);
                        BloomP.z.SetInt32(1);
                        
                        P = BloomP;
                        
                        privkey_num.clear();
                        index = 0;
                        for (auto& p : pow10_points) { // getting the index of the element in the bloomfilter
                            count = 0;
                            xc2 = secp256k1->GetXHex(&P.x, xC_len);
                            while (bf2.may_contain(xc2)) {
                                P = secp256k1->SubtractPoints(P, p);
                                xc2 = secp256k1->GetXHex(&P.x, xC_len);
                                count += 1;
                            }
                            privkey_num.push_back(pow10_nums[index] * (count - 1));
                            P = secp256k1->AddPoints(P, p);
                            index += 1;
                        }
                        steps = 0;
                        for (auto& n : privkey_num) { steps += n; } // we got here the index of the element in the bloomfilter
                        Int_steps.SetInt64(steps); // restoring the private key
                        batch_index.Mult(&stride, uint64_t(i + 1));
                        Int_temp.Add(&stride_sum, &batch_index);                        
                        Int_temp.Add(&offset);
                        Int_temp.Add(&Int_steps);
                        privkey.Add(&pre_calc_sum, &Int_temp);
                        privkey.Mult(mult);
                        privkey.AddOne(); // we got here the private key
                        calc_point = secp256k1->ScalarMultiplication(&privkey);
                        if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) { // if cpubs are equal we got it
                            print_time(); cout << "Privatekey: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            auto end = std::chrono::high_resolution_clock::now();
                            auto duration = end - start;
                            auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
                            duration -= hours;
                            auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
                            duration -= minutes;
                            auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
                            print_time(); cout << "Elapsed time: (" << hours.count() << ")hours (" << minutes.count() << ")minutes (" << seconds.count() << ")seconds\n";
                            exit(0);
                        }
                        print_time(); cout << "False Positive" << endl;
                    }
                }
                 
                startPoint.x.Set(&pointBatchX[POINTS_BATCH_SIZE - 1]); // setting the new startPoint for the next batch iteration
                startPoint.y.Set(&pointBatchY[POINTS_BATCH_SIZE - 1]);
                startPoint.z.SetInt32(1);
                
                stride_sum.Add(&batch_stride);
                
                if (threadIdx == 0) {  // thread with index zero is used to save the progress
                    save_counter += 1; // all values are derived from this data after new program start
                    if (save_counter % 70000 == 0) {
                        cpub = secp256k1->GetPublicKeyHex(startPoint);
                        ofstream outFile;
                        outFile.open("settings2.txt");
                        outFile << cpub <<'\n';
                        outFile << stride_sum.GetBase10() << '\n';
                        outFile.close();
                        save_counter = 0;
                        print_time(); cout << "Save Data written to settings2.txt" << endl;
                    }
                }
            }// while (true) loop end curly brace
        };
        
        std::thread subtraction_Threads[n_cores];
        for (int i = 0; i < n_cores; i++) {
            subtraction_Threads[i] = std::thread(scalable_subtraction_search, starting_Points[i], i , offset_Nums[i], stride_sum);
        }

        for (int i = 0; i < n_cores; i++) {
            subtraction_Threads[i].join();
        }
    };
    
    print_time(); cout << "Search in progress..." << endl;
    
    std::thread thread1(addition_search);
    std::thread thread2(subtraction_search);
    
    thread1.join();
    thread2.join();
}
