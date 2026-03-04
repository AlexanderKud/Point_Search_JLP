#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <cmath>

#include "secp256k1/SECP256k1.h"
#include "secp256k1/Int.h"
#include "secp256k1/IntGroup.h"
#include "util/util.h"

using namespace std;

const int cpuCores = std::thread::hardware_concurrency(); // actual number of processing cores
static constexpr int POINTS_BATCH_SIZE = 1024; // Batch addition with batch inversion using IntGroup class

auto main() -> int {

    Secp256K1* secp256k1 = new Secp256K1(); secp256k1->Init(); // initialize secp256k1 context

    Int gm; gm.SetBase16("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364140");
    Point Gm = secp256k1->ScalarMultiplication(&gm);
    
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
    getline(inFile, temp); range_start = std::stoull(temp);
    getline(inFile, temp); range_end = std::stoull(temp);
    getline(inFile, temp); block_width = std::stoull(temp);
    getline(inFile, temp); search_pub = trim(temp);
    inFile.close();

    Point TargetP1 = secp256k1->ParsePublicKeyHex(search_pub);
    Int d_05; d_05.SetBase10("57896044618658097711785492504343953926418782139537452191302581570759080747169");
    Point point_05 = secp256k1->ScalarMultiplication(&d_05);
    Point TargetP2 = secp256k1->AddPoints(TargetP1, point_05);

    print_time(); cout << "Range Start: " << range_start << " bits" << endl;
    print_time(); cout << "Range End  : " << range_end << " bits" << endl;
    print_time(); cout << "Block Width: 2^" << block_width << endl;
    print_time(); cout << "Search Pub : " << search_pub << endl;

    Int pre_calc_sum; // precalculated sum for private key recovering 512 + 256 = 768 for range[2^10..2^11] 1288 in example
    pre_calc_sum.Add(&S_table[range_start - 1], &S_table[range_start - 2]);
    uint64_t stride_bits = pow(2, block_width);

    uint64_t bloom_size = stride_bits * 4;
    uint64_t bloom_pos = bloom_size * 8;
    int iterations = 4;
    
    char * bloomfile1 = "bloom1.bf";
    char * bloomfile2 = "bloom2.bf";

    unsigned char * bloom1 = (unsigned char *)calloc(bloom_size, sizeof(unsigned char));
    unsigned char * bloom2 = (unsigned char *)calloc(bloom_size, sizeof(unsigned char));

    print_time(); cout << "Loading Bloomfilter images" << endl;
   
    auto bloom1_load = [&]() {
        load_bloom_filter(bloomfile1, bloom1, bloom_size);
    };
    
    auto bloom2_load = [&]() {
        load_bloom_filter(bloomfile2, bloom2, bloom_size);
    };
    
    std::thread bloom_load_thread1(bloom1_load);
    std::thread bloom_load_thread2(bloom2_load);
    
    bloom_load_thread1.join();
    bloom_load_thread2.join();
    
    auto pow10_nums = break_down_into_pow10(stride_bits); // decomposing the 2^block_width to the power of ten values
    size_t arr_size = pow10_nums.size();                  // to get the offset from the target point based on the bloomfilter hits fast
    Point pow10_points_Pos[arr_size];
    Point pow10_points_Neg[arr_size];                           
    Int pow_key;
    Point Pm;
    int arr_index = 0;
    for (auto& n : pow10_nums) { // calculating points corresponding to the decomposition components
        pow_key.SetInt64(n);
        Pm = secp256k1->ScalarMultiplication(&pow_key);
        pow10_points_Pos[arr_index] = Pm;
        Pm.y.ModNeg();   
        pow10_points_Neg[arr_index] = Pm;
        arr_index += 1;
    }
    
    auto chrono_start = std::chrono::high_resolution_clock::now();
    
    auto addition_search = [&]() { // addition search for the case when the starting point is behind the target point after calculations
                                   // the closer the target point to the center of the range from either side
        int save_counter = 0;      // the faster collision will happen
        string temp;
        Point start_point, stride_point, calc_point;
        Int stride_sum, stride;
        ifstream inFile("settings1.txt");
        getline(inFile, temp);
        start_point = secp256k1->ParsePublicKeyHex(trim(temp));
        getline(inFile, temp);
        stride_sum.SetBase10(trim(temp).data());
        inFile.close();
        
        stride.SetInt64(stride_bits);
        stride_point = secp256k1->ScalarMultiplication(&stride);
        
        //start splitting the search according to the chosen number of cpu cores
        Int offset_Step, int_Cores, vector_Num;
        int_Cores.SetInt32(cpuCores);
        offset_Step.floor_Div(&S_table[range_start - 2], &int_Cores);
        
        vector_Num.SetInt32(0);
        vector<Int> offset_Nums;
        for (int i = 0; i < cpuCores; i++) {
            offset_Nums.push_back(vector_Num);
            vector_Num.Add(&offset_Step);
        }
        
        vector<Point> offset_Points;
        offset_Points.push_back(secp256k1->G);
        for (int i = 1; i < cpuCores; i++) {
            offset_Points.push_back(secp256k1->ScalarMultiplication(&offset_Nums[i]));
        }

        vector<Point> starting_Points;
        Point vector_Point(start_point);
        starting_Points.push_back(vector_Point);

        for (int i = 1; i < cpuCores; i++) {
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
        // scalable_addition_search_save
        auto scalable_addition_search_save = [&](Point starting_Point, Int offset, Int stride_Sum) {

            Int stride_sum; stride_sum.Set(&stride_Sum);
            Int Int_steps, Int_temp, privkey;
            string cpub;
            int index, count;
            uint64_t steps;
            vector<uint64_t> privkey_num;
            
            IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            Int deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            modGroup.Set(deltaX); // assign array deltaX to modGroup for batch inversion (JLP way set it once)
            Int pointBatchX[POINTS_BATCH_SIZE]; // X coordinates of the batch
            Int pointBatchY[POINTS_BATCH_SIZE]; // Y coordinates of the batch
            Int deltaY; // values to store the results of points addition formula
            Int slope[POINTS_BATCH_SIZE];
            
            Point startPoint = starting_Point; // start point
            Point BloomP; // point for insertion of the batch into the bloomfilter
            Point CheckP;
        
            Int batch_stride, batch_index;
            batch_stride.Mult(&stride, uint64_t(POINTS_BATCH_SIZE));
            bool in_bloom;
        
            while (true) {

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                    deltaX[i].ModSub(&startPoint.x, &addPoints[i].x); // insert each entry into the deltaX array
                }
    
                modGroup.ModInv();    // doing batch inversion
                
                int i;
                for (i = 0; i < POINTS_BATCH_SIZE - 1; i++) { // follow points addition formula logic
                    
                    deltaY.ModSub(&startPoint.y, &addPoints[i].y);
                    slope[i].ModMulK1(&deltaY, &deltaX[i]); // deltaX already inverted for each entry of the batch

                    pointBatchX[i].ModSquareK1(&slope[i]); // computing just x coordinate for the (batch_size - 1)
                    pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                    pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                    
                }
                
                deltaY.ModSub(&startPoint.y, &addPoints[i].y); // computing the last entry of the batch full (x,y) coordinates
                slope[i].ModMulK1(&deltaY, &deltaX[i]); // deltaX already inverted for each entry of the batch

                pointBatchX[i].ModSquareK1(&slope[i]);
                pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                    
                pointBatchY[i].ModSub(&startPoint.x, &pointBatchX[i]);
                pointBatchY[i].ModMulK1(&slope[i], &pointBatchY[i]);
                pointBatchY[i].ModSub(&pointBatchY[i], &startPoint.y);

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) {

                    in_bloom = true;
                    for (int a = 0; a < iterations; a++) {
                        if (!check_bit(bloom1, pointBatchX[i].bits64[a] % bloom_pos)) {
                            in_bloom = false;
                            break;
                        }
                    }

                    if (in_bloom) {
                        
                        BloomP.x.Set(&pointBatchX[i]);
                        BloomP.y.ModSub(&startPoint.x, &pointBatchX[i]);
                        BloomP.y.ModMulK1(&slope[i], &BloomP.y);
                        BloomP.y.ModSub(&BloomP.y, &startPoint.y);

                        if (BloomP.x_equals(TargetP1)) {
                            Int_steps.SetInt64(0); // restoring the private key
                            batch_index.Mult(&stride, uint64_t(i + 1));
                            Int_temp.Add(&stride_sum, &batch_index);
                            Int_temp.Add(&offset);
                            Int_temp.Sub(&Int_steps);
                            privkey.Sub(&pre_calc_sum, &Int_temp);
                            privkey.Mult(mult); // we got here the private key
                            
                            print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Lower Range Half]" << endl;
                            print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }

                        CheckP = secp256k1->AddPoints(BloomP, Gm);
                        in_bloom = true;
                        for (int a = 0; a < iterations; a++) {
                            if (!check_bit(bloom1, CheckP.x.bits64[a] % bloom_pos)) {
                                in_bloom = false;
                                break;
                            }
                        }
                        
                        if (in_bloom) {

                            privkey_num.clear();
                            index = 0;
                            for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                                count = 0;
                                do {
                                    BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                                    count += 1;
                                    in_bloom = true;
                                    for (int c = 0; c < iterations; c++) {
                                        if (!check_bit(bloom1, BloomP.x.bits64[c] % bloom_pos)) {
                                            in_bloom = false;
                                            break;
                                        } 
                                    }
                                } while(in_bloom);
                                privkey_num.push_back(pow10_nums[index] * (count - 1));
                                BloomP = secp256k1->AddPoints(BloomP, pow10_points_Pos[i]);
                                index += 1;
                            }
                            
                            steps = 0;
                            for (auto& n : privkey_num) { steps += n; } // we got here the offset
                            Int_steps.SetInt64(steps); // restoring the private key
                            batch_index.Mult(&stride, uint64_t(i + 1));
                            Int_temp.Add(&stride_sum, &batch_index);
                            Int_temp.Add(&offset);
                            Int_temp.Sub(&Int_steps);
                            privkey.Sub(&pre_calc_sum, &Int_temp);
                            privkey.Mult(mult); // we got here the private key
                            calc_point = secp256k1->ScalarMultiplication(&privkey);
                            
                            if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) { // if cpubs are equal we got it
                                print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Lower Range Half]" << endl;
                                print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                                ofstream outFile;
                                outFile.open("found.txt", ios::app);
                                outFile << privkey.GetBase10() << '\n';
                                outFile.close();
                                print_elapsed_time(chrono_start);
                                exit(0);
                            }
                        }
                    }

                    in_bloom = true;
                    for (int b = 0; b < iterations; b++) {
                        if (!check_bit(bloom2, pointBatchX[i].bits64[b] % bloom_pos)) {
                            in_bloom = false;
                            break;
                        }
                    }
                    
                    if (in_bloom) {
                        
                        BloomP.x.Set(&pointBatchX[i]);
                        BloomP.y.ModSub(&startPoint.x, &pointBatchX[i]);
                        BloomP.y.ModMulK1(&slope[i], &BloomP.y);
                        BloomP.y.ModSub(&BloomP.y, &startPoint.y);

                        if (BloomP.x_equals(TargetP2)) {
                            Int_steps.SetInt64(0); // restoring the private key
                            batch_index.Mult(&stride, uint64_t(i + 1));
                            Int_temp.Add(&stride_sum, &batch_index);
                            Int_temp.Add(&offset);
                            Int_temp.Sub(&Int_steps);
                            privkey.Sub(&pre_calc_sum, &Int_temp);
                            privkey.Mult(mult);
                            privkey.AddOne(); // we got here the private key
                            
                            print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Lower Range Half]" << endl;
                            print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }

                        CheckP = secp256k1->AddPoints(BloomP, Gm);
                        in_bloom = true;
                        for (int a = 0; a < iterations; a++) {
                            if (!check_bit(bloom2, CheckP.x.bits64[a] % bloom_pos)) {
                                in_bloom = false;
                                break;
                            }
                        }

                        if (in_bloom) {
                        
                            privkey_num.clear();
                            index = 0;
                            for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                                count = 0;
                                do {
                                    BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                                    count += 1;
                                    in_bloom = true;
                                    for (int c = 0; c < iterations; c++) {
                                        if (!check_bit(bloom2, BloomP.x.bits64[c] % bloom_pos)) {
                                            in_bloom = false;
                                            break;
                                        } 
                                    }
                                } while(in_bloom);
                                privkey_num.push_back(pow10_nums[index] * (count - 1));
                                BloomP = secp256k1->AddPoints(BloomP, pow10_points_Pos[i]);
                                index += 1;
                            }
                                              
                            steps = 0;
                            for (auto& n : privkey_num) { steps += n; } // we got here the offset
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
                                print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Lower Range Half]" << endl;
                                print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                                ofstream outFile;
                                outFile.open("found.txt", ios::app);
                                outFile << privkey.GetBase10() << '\n';
                                outFile.close();
                                print_elapsed_time(chrono_start);
                                exit(0);
                            }
                        }
                    }
                }
                
                startPoint.x.Set(&pointBatchX[POINTS_BATCH_SIZE - 1]); // last batch entry as the new startPoint for the next batch iteration
                startPoint.y.Set(&pointBatchY[POINTS_BATCH_SIZE - 1]);
                
                stride_sum.Add(&batch_stride);
                save_counter += 1;
                
                if (save_counter % 300000 == 0) {
                    cpub = secp256k1->GetPublicKeyHex(startPoint);
                    ofstream outFile;
                    outFile.open("settings1.txt");
                    outFile << cpub <<'\n';
                    outFile << stride_sum.GetBase10() << '\n';
                    outFile.close();
                    save_counter = 0;
                    print_time(); cout << "Save Data written to settings1.txt" << endl;
                }
            } // while (true) loop end curly brace
        }; //scalable_addition_search_save
        // scalable_addition_search
        auto scalable_addition_search = [&](Point starting_Point, Int offset, Int stride_Sum) {

            Int stride_sum; stride_sum.Set(&stride_Sum);
            Int Int_steps, Int_temp, privkey;
            string cpub;
            int index, count;
            uint64_t steps;
            vector<uint64_t> privkey_num;
            
            IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            Int deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            modGroup.Set(deltaX); // assign array deltaX to modGroup for batch inversion (JLP way set it once)
            Int pointBatchX[POINTS_BATCH_SIZE]; // X coordinates of the batch
            Int pointBatchY[POINTS_BATCH_SIZE]; // Y coordinates of the batch
            Int deltaY; // values to store the results of points addition formula
            Int slope[POINTS_BATCH_SIZE];
            
            Point startPoint = starting_Point; // start point
            Point BloomP; // point for insertion of the batch into the bloomfilter
            Point CheckP;
        
            Int batch_stride, batch_index;
            batch_stride.Mult(&stride, uint64_t(POINTS_BATCH_SIZE));
            bool in_bloom;
        
            while (true) {

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                    deltaX[i].ModSub(&startPoint.x, &addPoints[i].x); // insert each entry into the deltaX array
                }
    
                modGroup.ModInv();    // doing batch inversion
                
                int i;
                for (i = 0; i < POINTS_BATCH_SIZE - 1; i++) { // follow points addition formula logic
                    
                    deltaY.ModSub(&startPoint.y, &addPoints[i].y);
                    slope[i].ModMulK1(&deltaY, &deltaX[i]); // deltaX already inverted for each entry of the batch

                    pointBatchX[i].ModSquareK1(&slope[i]); // computing just x coordinate for the (batch_size - 1)
                    pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                    pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                    
                }
                
                deltaY.ModSub(&startPoint.y, &addPoints[i].y); // computing the last entry of the batch full (x,y) coordinates
                slope[i].ModMulK1(&deltaY, &deltaX[i]); // deltaX already inverted for each entry of the batch

                pointBatchX[i].ModSquareK1(&slope[i]);
                pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                    
                pointBatchY[i].ModSub(&startPoint.x, &pointBatchX[i]);
                pointBatchY[i].ModMulK1(&slope[i], &pointBatchY[i]);
                pointBatchY[i].ModSub(&pointBatchY[i], &startPoint.y);

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) {
                    
                    in_bloom = true;
                    for (int a = 0; a < iterations; a++) {
                        if (!check_bit(bloom1, pointBatchX[i].bits64[a] % bloom_pos)) {
                            in_bloom = false;
                            break;
                        }
                    }

                    if (in_bloom) {
                        
                        BloomP.x.Set(&pointBatchX[i]);
                        BloomP.y.ModSub(&startPoint.x, &pointBatchX[i]);
                        BloomP.y.ModMulK1(&slope[i], &BloomP.y);
                        BloomP.y.ModSub(&BloomP.y, &startPoint.y);

                        if (BloomP.x_equals(TargetP1)) {
                            Int_steps.SetInt64(0); // restoring the private key
                            batch_index.Mult(&stride, uint64_t(i + 1));
                            Int_temp.Add(&stride_sum, &batch_index);
                            Int_temp.Add(&offset);
                            Int_temp.Sub(&Int_steps);
                            privkey.Sub(&pre_calc_sum, &Int_temp);
                            privkey.Mult(mult); // we got here the private key
                            
                            print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Lower Range Half]" << endl;
                            print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }

                        CheckP = secp256k1->AddPoints(BloomP, Gm);
                        in_bloom = true;
                        for (int a = 0; a < iterations; a++) {
                            if (!check_bit(bloom1, CheckP.x.bits64[a] % bloom_pos)) {
                                in_bloom = false;
                                break;
                            }
                        }

                        if (in_bloom) {
                        
                            privkey_num.clear();
                            index = 0;
                            for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                                count = 0;
                                do {
                                    BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                                    count += 1;
                                    in_bloom = true;
                                    for (int c = 0; c < iterations; c++) {
                                        if (!check_bit(bloom1, BloomP.x.bits64[c] % bloom_pos)) {
                                            in_bloom = false;
                                            break;
                                        } 
                                    }
                                } while(in_bloom);
                                privkey_num.push_back(pow10_nums[index] * (count - 1));
                                BloomP = secp256k1->AddPoints(BloomP, pow10_points_Pos[i]);
                                index += 1;
                            }
                            
                            steps = 0;
                            for (auto& n : privkey_num) { steps += n; } // we got here the offset
                            Int_steps.SetInt64(steps); // restoring the private key
                            batch_index.Mult(&stride, uint64_t(i + 1));
                            Int_temp.Add(&stride_sum, &batch_index);
                            Int_temp.Add(&offset);
                            Int_temp.Sub(&Int_steps);
                            privkey.Sub(&pre_calc_sum, &Int_temp);
                            privkey.Mult(mult); // we got here the private key
                            calc_point = secp256k1->ScalarMultiplication(&privkey);
                            
                            if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) { // if cpubs are equal we got it
                                print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Lower Range Half]" << endl;
                                print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                                ofstream outFile;
                                outFile.open("found.txt", ios::app);
                                outFile << privkey.GetBase10() << '\n';
                                outFile.close();
                                print_elapsed_time(chrono_start);
                                exit(0);
                            }
                        }
                    }

                    in_bloom = true;
                    for (int b = 0; b < iterations; b++) {
                        if (!check_bit(bloom2, pointBatchX[i].bits64[b] % bloom_pos)) {
                            in_bloom = false;
                            break;
                        }
                    }
                    
                    if (in_bloom) {
                        
                        BloomP.x.Set(&pointBatchX[i]);
                        BloomP.y.ModSub(&startPoint.x, &pointBatchX[i]);
                        BloomP.y.ModMulK1(&slope[i], &BloomP.y);
                        BloomP.y.ModSub(&BloomP.y, &startPoint.y);

                        if (BloomP.x_equals(TargetP2)) {
                            Int_steps.SetInt64(0); // restoring the private key
                            batch_index.Mult(&stride, uint64_t(i + 1));
                            Int_temp.Add(&stride_sum, &batch_index);
                            Int_temp.Add(&offset);
                            Int_temp.Sub(&Int_steps);
                            privkey.Sub(&pre_calc_sum, &Int_temp);
                            privkey.Mult(mult);
                            privkey.AddOne(); // we got here the private key
                            
                            print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Lower Range Half]" << endl;
                            print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }

                        CheckP = secp256k1->AddPoints(BloomP, Gm);
                        in_bloom = true;
                        for (int a = 0; a < iterations; a++) {
                            if (!check_bit(bloom2, CheckP.x.bits64[a] % bloom_pos)) {
                                in_bloom = false;
                                break;
                            }
                        }

                        if (in_bloom) {
                        
                            privkey_num.clear();
                            index = 0;
                            for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                                count = 0;
                                do {
                                    BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                                    count += 1;
                                    in_bloom = true;
                                    for (int c = 0; c < iterations; c++) {
                                        if (!check_bit(bloom2, BloomP.x.bits64[c] % bloom_pos)) {
                                            in_bloom = false;
                                            break;
                                        } 
                                    }
                                } while(in_bloom);
                                privkey_num.push_back(pow10_nums[index] * (count - 1));
                                BloomP = secp256k1->AddPoints(BloomP, pow10_points_Pos[i]);
                                index += 1;
                            }
                                              
                            steps = 0;
                            for (auto& n : privkey_num) { steps += n; } // we got here the offset
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
                                print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Lower Range Half]" << endl;
                                print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                                ofstream outFile;
                                outFile.open("found.txt", ios::app);
                                outFile << privkey.GetBase10() << '\n';
                                outFile.close();
                                print_elapsed_time(chrono_start);
                                exit(0);
                            } 
                        }
                    }
                }
                
                startPoint.x.Set(&pointBatchX[POINTS_BATCH_SIZE - 1]); // last batch entry as the new startPoint for the next batch iteration
                startPoint.y.Set(&pointBatchY[POINTS_BATCH_SIZE - 1]);
                
                stride_sum.Add(&batch_stride);

            } // while (true) loop end curly brace
        }; // scalable_addition_search
        
        std::thread addition_Threads[cpuCores];
        
        addition_Threads[0] = std::thread(scalable_addition_search_save, starting_Points[0], offset_Nums[0], stride_sum);
        for (int i = 1; i < cpuCores; i++) {
            addition_Threads[i] = std::thread(scalable_addition_search, starting_Points[i], offset_Nums[i], stride_sum);
        }

        for (int i = 0; i < cpuCores; i++) {
            addition_Threads[i].join();
        }
    };
    
    auto subtraction_search = [&]() { // subtraction search for the case when the starting point is ahead of the target point after calculations
                                      // the closer the target point to the center of the range from either side
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
        
        stride.SetInt64(stride_bits);
        stride_point = secp256k1->ScalarMultiplication(&stride);
        
        //start splitting the search according to the chosen number of cpu cores
        Int offset_Step, int_Cores, vector_Num;
        int_Cores.SetInt32(cpuCores);
        offset_Step.floor_Div(&S_table[range_start - 2], &int_Cores);
        
        vector_Num.SetInt32(0);
        vector<Int> offset_Nums;
        for (int i = 0; i < cpuCores; i++) {
            offset_Nums.push_back(vector_Num);
            vector_Num.Add(&offset_Step);
        }
        
        vector<Point> offset_Points;
        offset_Points.push_back(secp256k1->G);
        for (int i = 1; i < cpuCores; i++) {
            offset_Points.push_back(secp256k1->ScalarMultiplication(&offset_Nums[i]));
        }

        vector<Point> starting_Points;
        Point vector_Point(start_point);
        starting_Points.push_back(vector_Point);

        for (int i = 1; i < cpuCores; i++) {
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
        // scalable_subtraction_search_save
        auto scalable_subtraction_search_save = [&](Point starting_Point, Int offset, Int stride_Sum) {

            Int stride_sum; stride_sum.Set(&stride_Sum);
            string cpub;
            int index, count;
            uint64_t steps;
            vector<uint64_t> privkey_num;
            Int Int_steps, Int_temp, privkey;
            
            IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            Int deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            modGroup.Set(deltaX); // assign array deltaX to modGroup for batch inversion (JLP way set it once)
            Int pointBatchX[POINTS_BATCH_SIZE]; // X coordinates of the batch
            Int pointBatchY[POINTS_BATCH_SIZE]; // Y coordinates of the batch
            Int deltaY; // values to store the results of points addition formula
            Int slope[POINTS_BATCH_SIZE];
            
            Point startPoint = starting_Point; // start point
            Point BloomP; // point for insertion of the batch into the bloomfilter
            Point CheckP;
            
            Int batch_stride, batch_index;
            batch_stride.Mult(&stride, uint64_t(POINTS_BATCH_SIZE));
            bool in_bloom;
        
            while (true) {

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                    deltaX[i].ModSub(&startPoint.x, &addPoints[i].x); // insert each entry into the deltaX array
                }
    
                modGroup.ModInv();    // doing batch inversion
                
                int i;
                for (i = 0; i < POINTS_BATCH_SIZE - 1; i++) { // follow points addition formula logic
                    
                    deltaY.ModSub(&startPoint.y, &addPoints[i].y);
                    slope[i].ModMulK1(&deltaY, &deltaX[i]); // deltaX already inverted for each entry of the batch

                    pointBatchX[i].ModSquareK1(&slope[i]); // computing just x coordinate for the (batch_size - 1)
                    pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                    pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                    
                }
                
                deltaY.ModSub(&startPoint.y, &addPoints[i].y); // computing the last entry of the batch full (x,y) coordinates
                slope[i].ModMulK1(&deltaY, &deltaX[i]); // deltaX already inverted for each entry of the batch

                pointBatchX[i].ModSquareK1(&slope[i]);
                pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                    
                pointBatchY[i].ModSub(&startPoint.x, &pointBatchX[i]);
                pointBatchY[i].ModMulK1(&slope[i], &pointBatchY[i]);
                pointBatchY[i].ModSub(&pointBatchY[i], &startPoint.y);

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) {
                    
                    in_bloom = true;
                    for (int a = 0; a < iterations; a++) {
                        if (!check_bit(bloom1, pointBatchX[i].bits64[a] % bloom_pos)) {
                            in_bloom = false;
                            break;
                        }
                    }

                    if (in_bloom) {
                        
                        BloomP.x.Set(&pointBatchX[i]);
                        BloomP.y.ModSub(&startPoint.x, &pointBatchX[i]);
                        BloomP.y.ModMulK1(&slope[i], &BloomP.y);
                        BloomP.y.ModSub(&BloomP.y, &startPoint.y);

                        if (BloomP.x_equals(TargetP1)) {
                            Int_steps.SetInt64(0); // restoring the private key
                            batch_index.Mult(&stride, uint64_t(i + 1));
                            Int_temp.Add(&stride_sum, &batch_index);                        
                            Int_temp.Add(&offset);
                            Int_temp.Add(&Int_steps);
                            privkey.Add(&pre_calc_sum, &Int_temp);
                            privkey.Mult(mult); // we got here the private key

                            print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Higher Range Half]" << endl;
                            print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }

                        CheckP = secp256k1->AddPoints(BloomP, Gm);
                        in_bloom = true;
                        for (int a = 0; a < iterations; a++) {
                            if (!check_bit(bloom1, CheckP.x.bits64[a] % bloom_pos)) {
                                in_bloom = false;
                                break;
                            }
                        }

                        if (in_bloom) {
                        
                            privkey_num.clear();
                            index = 0;
                            for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                                count = 0;
                                do {
                                    BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                                    count += 1;
                                    in_bloom = true;
                                    for (int c = 0; c < iterations; c++) {
                                        if (!check_bit(bloom1, BloomP.x.bits64[c] % bloom_pos)) {
                                            in_bloom = false;
                                            break;
                                        } 
                                    }
                                } while(in_bloom);
                                privkey_num.push_back(pow10_nums[index] * (count - 1));
                                BloomP = secp256k1->AddPoints(BloomP, pow10_points_Pos[i]);
                                index += 1;
                            }
                                               
                            steps = 0;
                            for (auto& n : privkey_num) { steps += n; } // we got here the offset
                            Int_steps.SetInt64(steps); // restoring the private key
                            batch_index.Mult(&stride, uint64_t(i + 1));
                            Int_temp.Add(&stride_sum, &batch_index);                        
                            Int_temp.Add(&offset);
                            Int_temp.Add(&Int_steps);
                            privkey.Add(&pre_calc_sum, &Int_temp);
                            privkey.Mult(mult); // we got here the private key
                            calc_point = secp256k1->ScalarMultiplication(&privkey);
                            
                            if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) { // if cpubs are equal we got it
                                print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Higher Range Half]" << endl;
                                print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                                ofstream outFile;
                                outFile.open("found.txt", ios::app);
                                outFile << privkey.GetBase10() << '\n';
                                outFile.close();
                                print_elapsed_time(chrono_start);
                                exit(0);
                            }
                        }
                    }

                    in_bloom = true;
                    for (int b = 0; b < iterations; b++) {
                        if (!check_bit(bloom2, pointBatchX[i].bits64[b] % bloom_pos)) {
                            in_bloom = false;
                            break;
                        }
                    }
                    
                    if (in_bloom) {
                        
                        BloomP.x.Set(&pointBatchX[i]);
                        BloomP.y.ModSub(&startPoint.x, &pointBatchX[i]);
                        BloomP.y.ModMulK1(&slope[i], &BloomP.y);
                        BloomP.y.ModSub(&BloomP.y, &startPoint.y);

                        if (BloomP.x_equals(TargetP2)) {
                            Int_steps.SetInt64(0); // restoring the private key
                            batch_index.Mult(&stride, uint64_t(i + 1));
                            Int_temp.Add(&stride_sum, &batch_index);                        
                            Int_temp.Add(&offset);
                            Int_temp.Add(&Int_steps);
                            privkey.Add(&pre_calc_sum, &Int_temp);
                            privkey.Mult(mult);
                            privkey.AddOne(); // we got here the private key

                            print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Higher Range Half]" << endl;
                            print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }

                        CheckP = secp256k1->AddPoints(BloomP, Gm);
                        in_bloom = true;
                        for (int a = 0; a < iterations; a++) {
                            if (!check_bit(bloom2, CheckP.x.bits64[a] % bloom_pos)) {
                                in_bloom = false;
                                break;
                            }
                        }

                        if (in_bloom) {
                        
                            privkey_num.clear();
                            index = 0;
                            for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                                count = 0;
                                do {
                                    BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                                    count += 1;
                                    in_bloom = true;
                                    for (int c = 0; c < iterations; c++) {
                                        if (!check_bit(bloom2, BloomP.x.bits64[c] % bloom_pos)) {
                                            in_bloom = false;
                                            break;
                                        } 
                                    }
                                } while(in_bloom);
                                privkey_num.push_back(pow10_nums[index] * (count - 1));
                                BloomP = secp256k1->AddPoints(BloomP, pow10_points_Pos[i]);
                                index += 1;
                            }
                            
                            steps = 0;
                            for (auto& n : privkey_num) { steps += n; } // we got here the offset
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
                                print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Higher Range Half]" << endl;
                                print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                                ofstream outFile;
                                outFile.open("found.txt", ios::app);
                                outFile << privkey.GetBase10() << '\n';
                                outFile.close();
                                print_elapsed_time(chrono_start);
                                exit(0);
                            }
                        }
                    }
                }
                 
                startPoint.x.Set(&pointBatchX[POINTS_BATCH_SIZE - 1]); // last batch entry as the new startPoint for the next batch iteration
                startPoint.y.Set(&pointBatchY[POINTS_BATCH_SIZE - 1]); 
                
                stride_sum.Add(&batch_stride);
                save_counter += 1; // all values are derived from this data after new program start
                
                if (save_counter % 300000 == 0) {
                    cpub = secp256k1->GetPublicKeyHex(startPoint);
                    ofstream outFile;
                    outFile.open("settings2.txt");
                    outFile << cpub <<'\n';
                    outFile << stride_sum.GetBase10() << '\n';
                    outFile.close();
                    save_counter = 0;
                    print_time(); cout << "Save Data written to settings2.txt" << endl;
                }
            }// while (true) loop end curly brace
        };// scalable_subtraction_search_save
        // scalable_subtraction_search
        auto scalable_subtraction_search = [&](Point starting_Point, Int offset, Int stride_Sum) {

            Int stride_sum; stride_sum.Set(&stride_Sum);
            string cpub;
            int index, count;
            uint64_t steps;
            vector<uint64_t> privkey_num;
            Int Int_steps, Int_temp, privkey;
            
            IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            Int deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            modGroup.Set(deltaX); // assign array deltaX to modGroup for batch inversion (JLP way set it once)
            Int pointBatchX[POINTS_BATCH_SIZE]; // X coordinates of the batch
            Int pointBatchY[POINTS_BATCH_SIZE]; // Y coordinates of the batch
            Int deltaY; // values to store the results of points addition formula
            Int slope[POINTS_BATCH_SIZE];
            
            Point startPoint = starting_Point; // start point
            Point BloomP; // point for insertion of the batch into the bloomfilter
            Point CheckP;
            
            Int batch_stride, batch_index;
            batch_stride.Mult(&stride, uint64_t(POINTS_BATCH_SIZE));
            bool in_bloom;
        
            while (true) {

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                    deltaX[i].ModSub(&startPoint.x, &addPoints[i].x); // insert each entry into the deltaX array
                }
    
                modGroup.ModInv();    // doing batch inversion
                
                int i;
                for (i = 0; i < POINTS_BATCH_SIZE - 1; i++) { // follow points addition formula logic
                    
                    deltaY.ModSub(&startPoint.y, &addPoints[i].y);
                    slope[i].ModMulK1(&deltaY, &deltaX[i]); // deltaX already inverted for each entry of the batch

                    pointBatchX[i].ModSquareK1(&slope[i]); // computing just x coordinate for the (batch_size - 1)
                    pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                    pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                    
                }
                
                deltaY.ModSub(&startPoint.y, &addPoints[i].y); // computing the last entry of the batch full (x,y) coordinates
                slope[i].ModMulK1(&deltaY, &deltaX[i]); // deltaX already inverted for each entry of the batch

                pointBatchX[i].ModSquareK1(&slope[i]);
                pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                    
                pointBatchY[i].ModSub(&startPoint.x, &pointBatchX[i]);
                pointBatchY[i].ModMulK1(&slope[i], &pointBatchY[i]);
                pointBatchY[i].ModSub(&pointBatchY[i], &startPoint.y);

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) {

                    in_bloom = true;
                    for (int a = 0; a < iterations; a++) {
                        if (!check_bit(bloom1, pointBatchX[i].bits64[a] % bloom_pos)) {
                            in_bloom = false;
                            break;
                        }
                    }

                    if (in_bloom) {
                        
                        BloomP.x.Set(&pointBatchX[i]);
                        BloomP.y.ModSub(&startPoint.x, &pointBatchX[i]);
                        BloomP.y.ModMulK1(&slope[i], &BloomP.y);
                        BloomP.y.ModSub(&BloomP.y, &startPoint.y);

                        if (BloomP.x_equals(TargetP1)) {
                            Int_steps.SetInt64(0); // restoring the private key
                            batch_index.Mult(&stride, uint64_t(i + 1));
                            Int_temp.Add(&stride_sum, &batch_index);                        
                            Int_temp.Add(&offset);
                            Int_temp.Add(&Int_steps);
                            privkey.Add(&pre_calc_sum, &Int_temp);
                            privkey.Mult(mult); // we got here the private key

                            print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Higher Range Half]" << endl;
                            print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }

                        CheckP = secp256k1->AddPoints(BloomP, Gm);
                        in_bloom = true;
                        for (int a = 0; a < iterations; a++) {
                            if (!check_bit(bloom1, CheckP.x.bits64[a] % bloom_pos)) {
                                in_bloom = false;
                                break;
                            }
                        }

                        if (in_bloom) {
                        
                            privkey_num.clear();
                            index = 0;
                            for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                                count = 0;
                                do {
                                    BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                                    count += 1;
                                    in_bloom = true;
                                    for (int c = 0; c < iterations; c++) {
                                        if (!check_bit(bloom1, BloomP.x.bits64[c] % bloom_pos)) {
                                            in_bloom = false;
                                            break;
                                        } 
                                    }
                                } while(in_bloom);
                                privkey_num.push_back(pow10_nums[index] * (count - 1));
                                BloomP = secp256k1->AddPoints(BloomP, pow10_points_Pos[i]);
                                index += 1;
                            }
                                               
                            steps = 0;
                            for (auto& n : privkey_num) { steps += n; } // we got here the offset
                            Int_steps.SetInt64(steps); // restoring the private key
                            batch_index.Mult(&stride, uint64_t(i + 1));
                            Int_temp.Add(&stride_sum, &batch_index);                        
                            Int_temp.Add(&offset);
                            Int_temp.Add(&Int_steps);
                            privkey.Add(&pre_calc_sum, &Int_temp);
                            privkey.Mult(mult); // we got here the private key
                            calc_point = secp256k1->ScalarMultiplication(&privkey);
                            
                            if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) { // if cpubs are equal we got it
                                print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Higher Range Half]" << endl;
                                print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                                ofstream outFile;
                                outFile.open("found.txt", ios::app);
                                outFile << privkey.GetBase10() << '\n';
                                outFile.close();
                                print_elapsed_time(chrono_start);
                                exit(0);
                            }
                        }
                    }

                    in_bloom = true;
                    for (int b = 0; b < iterations; b++) {
                        if (!check_bit(bloom2, pointBatchX[i].bits64[b] % bloom_pos)) {
                            in_bloom = false;
                            break;
                        }
                    }
                    
                    if (in_bloom) {
                        
                        BloomP.x.Set(&pointBatchX[i]);
                        BloomP.y.ModSub(&startPoint.x, &pointBatchX[i]);
                        BloomP.y.ModMulK1(&slope[i], &BloomP.y);
                        BloomP.y.ModSub(&BloomP.y, &startPoint.y);

                        if (BloomP.x_equals(TargetP2)) {
                            Int_steps.SetInt64(0); // restoring the private key
                            batch_index.Mult(&stride, uint64_t(i + 1));
                            Int_temp.Add(&stride_sum, &batch_index);                        
                            Int_temp.Add(&offset);
                            Int_temp.Add(&Int_steps);
                            privkey.Add(&pre_calc_sum, &Int_temp);
                            privkey.Mult(mult);
                            privkey.AddOne(); // we got here the private key

                            print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Higher Range Half]" << endl;
                            print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }

                        CheckP = secp256k1->AddPoints(BloomP, Gm);
                        in_bloom = true;
                        for (int a = 0; a < iterations; a++) {
                            if (!check_bit(bloom2, CheckP.x.bits64[a] % bloom_pos)) {
                                in_bloom = false;
                                break;
                            }
                        }

                        if (in_bloom) {
                        
                            privkey_num.clear();
                            index = 0;
                            for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                                count = 0;
                                do {
                                    BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                                    count += 1;
                                    in_bloom = true;
                                    for (int c = 0; c < iterations; c++) {
                                        if (!check_bit(bloom2, BloomP.x.bits64[c] % bloom_pos)) {
                                            in_bloom = false;
                                            break;
                                        } 
                                    }
                                } while(in_bloom);
                                privkey_num.push_back(pow10_nums[index] * (count - 1));
                                BloomP = secp256k1->AddPoints(BloomP, pow10_points_Pos[i]);
                                index += 1;
                            }
                            
                            steps = 0;
                            for (auto& n : privkey_num) { steps += n; } // we got here the offset
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
                                print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Higher Range Half]" << endl;
                                print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                                ofstream outFile;
                                outFile.open("found.txt", ios::app);
                                outFile << privkey.GetBase10() << '\n';
                                outFile.close();
                                print_elapsed_time(chrono_start);
                                exit(0);
                            }
                        }
                    }
                }
                 
                startPoint.x.Set(&pointBatchX[POINTS_BATCH_SIZE - 1]); // last batch entry as the new startPoint for the next batch iteration
                startPoint.y.Set(&pointBatchY[POINTS_BATCH_SIZE - 1]); 
                
                stride_sum.Add(&batch_stride);
                
            } // while (true) loop end curly brace
        }; // scalable_subtraction_search
        
        std::thread subtraction_Threads[cpuCores];
        
        subtraction_Threads[0] = std::thread(scalable_subtraction_search_save, starting_Points[0], offset_Nums[0], stride_sum);
        for (int i = 1; i < cpuCores; i++) {
            subtraction_Threads[i] = std::thread(scalable_subtraction_search, starting_Points[i], offset_Nums[i], stride_sum);
        }

        for (int i = 0; i < cpuCores; i++) {
            subtraction_Threads[i].join();
        }
    };
    
    print_time(); cout << "Search in progress..." << endl;
    
    std::thread thread1(addition_search);
    std::thread thread2(subtraction_search);
    
    thread1.join();
    thread2.join();
}
