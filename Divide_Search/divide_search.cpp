#include <iostream>
#include <fstream>
#include <thread>
#include <utility>

#include "secp256k1/SECP256k1.h"
#include "secp256k1/Int.h"
#include "secp256k1/IntGroup.h"
#include "util/util.h"

using namespace std;

const int cpuCores = std::thread::hardware_concurrency();

auto main() -> int {

    Secp256K1* secp256k1 = new Secp256K1(); 
    secp256k1->Init();

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

    Point TargetP = secp256k1->ParsePublicKeyHex(search_pub);

    print_time(); cout << "Range Start: " << range_start << " bits" << endl;
    print_time(); cout << "Range End  : " << range_end << " bits" << endl;
    print_time(); cout << "Block Width: 2^" << block_width << endl;
    print_time(); cout << "Search Pub : " << search_pub << endl;

    uint64_t stride_bits = 1 << block_width;

    uint64_t bloom_size = stride_bits * 4;
    uint64_t bloom_mod = bloom_size * 8;
    int iterations = 4;
    
    char * bloomfile = "bloom.bf";

    unsigned char * bloom = (unsigned char *)calloc(bloom_size, sizeof(unsigned char));

    print_time(); cout << "Loading Bloomfilter image" << endl;
   
    load_bloom_filter(bloomfile, bloom, bloom_size);
    
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

    auto divide_search = [&]() {

        int exp = block_width / 2;

        string temp;

        Int steps_sum; // reading steps_sum from file
        ifstream inFile("steps_sum.txt");
        getline(inFile, temp);
        steps_sum.SetBase10(trim(temp).data());
        inFile.close();

        Int range_start_Int, divide, step, rem;

        range_start_Int.Set(&S_table[range_start]); // range_start
        divide.SetInt32(cpuCores); // divide = cpuCores                 
        step.Set(&range_start_Int);
        step.Div(&divide, &rem); // step = range_start / cpuCores


        std::pair<Int , Int> range_nums[cpuCores]; // whole range divided into pairs (range start + step)

        for (int i = 0; i < cpuCores; i++) {
            range_nums[i].first.Set(&range_start_Int);
            range_start_Int.Add(&step);
            range_nums[i].second.Set(&range_start_Int);
        }
        
        Int stride;
        stride.SetInt64(stride_bits);
        Point stride_point = secp256k1->ScalarMultiplication(&stride);

        Int offset;
        offset.Mult(&stride, &steps_sum);

        for (int i = 0; i < cpuCores; i++) {
            range_nums[i].first.Add(&offset);
            range_nums[i].second.Sub(&offset);
        }

        std::pair<Point , Point> range_pts[cpuCores];
        Point P, Q;
        for (int i = 0; i < cpuCores; i++) {
            P = secp256k1->ScalarMultiplication(&range_nums[i].first);
            Q = secp256k1->ScalarMultiplication(&range_nums[i].second);
            range_pts[i].first = P;
            range_pts[i].second = Q;
        }

        auto scalable_divide_search_save = [&](const std::pair<Int, Int>& keys, const std::pair<Point, Point>& pts) {

            int save_counter = 0;
            
            Int min_value(keys.first);
            Int max_value(keys.second);

            Int between, walker, walker_stride, rem, loop;
            Point walker_stride_point, walker_point;
            
            Point min_point = pts.first;
            Point max_point = pts.second;

            Int divide((uint64_t) 1 << exp);

            bool in_bloom;
            uint64_t steps;
            vector<uint64_t> privkey_num;
            Point BloomP, calc_point;
            int index, count;
            Int privkey, Int_steps;
            
            while (true) {

                in_bloom = true;
                for (int a = 0; a < iterations; a++) {
                    if (!check_bit(bloom, min_point.x.bits64[a] & (bloom_mod - 1))) {
                        in_bloom = false;
                        break;
                    }
                }

                if (in_bloom) {

                    BloomP = min_point;

                    privkey_num.clear();
                    index = 0;
                    for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                        count = 0;
                        do {
                            BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                            count += 1;
                            in_bloom = true;
                            for (int c = 0; c < iterations; c++) {
                                if (!check_bit(bloom, BloomP.x.bits64[c] & (bloom_mod - 1))) {
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
                    Int_steps.SetInt64(steps);
                    privkey.Sub(&min_value, &Int_steps);

                    //print_time(); cout << "Private key(save min): " << privkey.GetBase10() << endl;

                    calc_point = secp256k1->ScalarMultiplication(&privkey);
                            
                    if (calc_point.x_equals(TargetP)) {
                        print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                        ofstream outFile;
                        outFile.open("found.txt", ios::app);
                        outFile << privkey.GetBase10() << '\n';
                        outFile.close();
                        print_elapsed_time(chrono_start);
                        exit(0);
                    }

                }

                in_bloom = true;
                for (int a = 0; a < iterations; a++) {
                    if (!check_bit(bloom, max_point.x.bits64[a] & (bloom_mod - 1))) {
                        in_bloom = false;
                        break;
                    }
                }

                if (in_bloom) {

                    BloomP = max_point;

                    privkey_num.clear();
                    index = 0;
                    for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                        count = 0;
                        do {
                            BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                            count += 1;
                            in_bloom = true;
                            for (int c = 0; c < iterations; c++) {
                                if (!check_bit(bloom, BloomP.x.bits64[c] & (bloom_mod - 1))) {
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
                    Int_steps.SetInt64(steps);
                    privkey.Sub(&max_value, &Int_steps);

                    //print_time(); cout << "Private key(save max): " << privkey.GetBase10() << endl;

                    calc_point = secp256k1->ScalarMultiplication(&privkey);
                            
                    if (calc_point.x_equals(TargetP)) {
                        print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                        ofstream outFile;
                        outFile.open("found.txt", ios::app);
                        outFile << privkey.GetBase10() << '\n';
                        outFile.close();
                        print_elapsed_time(chrono_start);
                        exit(0);
                    }
                    
                }



                between.Sub(&max_value, &min_value);

                walker.Set(&min_value);
                walker_stride.Set(&between);
                walker_stride.Div(&divide, &rem);
                walker_point = secp256k1->ScalarMultiplication(&walker);
                walker_stride_point = secp256k1->ScalarMultiplication(&walker_stride);

                for(loop.SetBase10("0"); loop.IsLower(&divide); loop.AddOne()) {

                    walker_point = secp256k1->AddPoints(walker_point, walker_stride_point);
                    walker.Add(&walker_stride);

                    in_bloom = true; // max_point check
                    for (int a = 0; a < iterations; a++) {
                        if (!check_bit(bloom, walker_point.x.bits64[a] & (bloom_mod - 1))) {
                            in_bloom = false;
                            break;
                        }
                    }

                    if (in_bloom) {

                        BloomP = walker_point;

                        privkey_num.clear();
                        index = 0;
                        for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                            count = 0;
                            do {
                                BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                                count += 1;
                                in_bloom = true;
                                for (int c = 0; c < iterations; c++) {
                                    if (!check_bit(bloom, BloomP.x.bits64[c] & (bloom_mod - 1))) {
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
                        Int_steps.SetInt64(steps);
                        privkey.Sub(&walker, &Int_steps);

                        //print_time(); cout << "Private key(save walker): " << privkey.GetBase10() << endl;

                        calc_point = secp256k1->ScalarMultiplication(&privkey);
                            
                        if (calc_point.x_equals(TargetP)) {
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

                min_value.Add(&stride);
                max_value.Sub(&stride);
                min_point = secp256k1->AddPoints(min_point, stride_point);
                max_point = secp256k1->SubtractPoints(max_point, stride_point);

                steps_sum.AddOne();
                save_counter += 1;

                if (save_counter % 3000 == 0) {
                        ofstream outFile;
                        outFile.open("steps_sum.txt");
                        outFile << steps_sum.GetBase10() << '\n';
                        outFile.close();
                        save_counter = 0;
                        print_time(); cout << "Save Data written to steps_sum.txt" << endl;
                }
            }

        };

        auto scalable_divide_search = [&](const std::pair<Int, Int>& keys, const std::pair<Point, Point>& pts) {

            Int min_value(keys.first);
            Int max_value(keys.second);

            Int between, walker, walker_stride, rem, loop;
            Point walker_stride_point, walker_point;
            
            Point min_point = pts.first;
            Point max_point = pts.second;

            Int divide((uint64_t) 1 << exp);

            bool in_bloom;
            uint64_t steps;
            vector<uint64_t> privkey_num;
            Point BloomP, calc_point;
            int index, count;
            Int privkey, Int_steps;

            while (true) {

                in_bloom = true; // min_point check
                for (int a = 0; a < iterations; a++) {
                    if (!check_bit(bloom, min_point.x.bits64[a] & (bloom_mod - 1))) {
                        in_bloom = false;
                        break;
                    }
                }

                if (in_bloom) {

                    BloomP = min_point;

                    privkey_num.clear();
                    index = 0;
                    for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                        count = 0;
                        do {
                            BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                            count += 1;
                            in_bloom = true;
                            for (int c = 0; c < iterations; c++) {
                                if (!check_bit(bloom, BloomP.x.bits64[c] & (bloom_mod - 1))) {
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
                    Int_steps.SetInt64(steps);
                    privkey.Sub(&min_value, &Int_steps);

                    //print_time(); cout << "Private key(min): " << privkey.GetBase10() << endl;

                    calc_point = secp256k1->ScalarMultiplication(&privkey);
                            
                    if (calc_point.x_equals(TargetP)) {
                        print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                        ofstream outFile;
                        outFile.open("found.txt", ios::app);
                        outFile << privkey.GetBase10() << '\n';
                        outFile.close();
                        print_elapsed_time(chrono_start);
                        exit(0);
                    }

                }

                in_bloom = true; // max_point check
                for (int a = 0; a < iterations; a++) {
                    if (!check_bit(bloom, max_point.x.bits64[a] & (bloom_mod - 1))) {
                        in_bloom = false;
                        break;
                    }
                }

                if (in_bloom) {

                    BloomP = max_point;

                    privkey_num.clear();
                    index = 0;
                    for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                        count = 0;
                        do {
                            BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                            count += 1;
                            in_bloom = true;
                            for (int c = 0; c < iterations; c++) {
                                if (!check_bit(bloom, BloomP.x.bits64[c] & (bloom_mod - 1))) {
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
                    Int_steps.SetInt64(steps);
                    privkey.Sub(&max_value, &Int_steps);

                    //print_time(); cout << "Private key(max): " << privkey.GetBase10() << endl;

                    calc_point = secp256k1->ScalarMultiplication(&privkey);
                            
                    if (calc_point.x_equals(TargetP)) {
                        print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                        ofstream outFile;
                        outFile.open("found.txt", ios::app);
                        outFile << privkey.GetBase10() << '\n';
                        outFile.close();
                        print_elapsed_time(chrono_start);
                        exit(0);
                    }
                    
                }



                between.Sub(&max_value, &min_value);

                walker.Set(&min_value);
                walker_stride.Set(&between);
                walker_stride.Div(&divide, &rem);
                walker_point = secp256k1->ScalarMultiplication(&walker);
                walker_stride_point = secp256k1->ScalarMultiplication(&walker_stride);

                for(loop.SetBase10("0"); loop.IsLower(&divide); loop.AddOne()) {

                    walker_point = secp256k1->AddPoints(walker_point, walker_stride_point);
                    walker.Add(&walker_stride);

                    in_bloom = true; // max_point check
                    for (int a = 0; a < iterations; a++) {
                        if (!check_bit(bloom, walker_point.x.bits64[a] & (bloom_mod - 1))) {
                            in_bloom = false;
                            break;
                        }
                    }

                    if (in_bloom) {

                        BloomP = walker_point;

                        privkey_num.clear();
                        index = 0;
                        for (size_t i = 0; i < arr_size; i++) { // getting the offset from the target point based on bloomfilter hits
                            count = 0;
                            do {
                                BloomP = secp256k1->AddPoints(BloomP, pow10_points_Neg[i]);
                                count += 1;
                                in_bloom = true;
                                for (int c = 0; c < iterations; c++) {
                                    if (!check_bit(bloom, BloomP.x.bits64[c] & (bloom_mod - 1))) {
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
                        Int_steps.SetInt64(steps);
                        privkey.Sub(&walker, &Int_steps);

                        //print_time(); cout << "Private key(walker): " << privkey.GetBase10() << endl;

                        calc_point = secp256k1->ScalarMultiplication(&privkey);
                            
                        if (calc_point.x_equals(TargetP)) {
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

                min_value.Add(&stride);
                max_value.Sub(&stride);
                min_point = secp256k1->AddPoints(min_point, stride_point);
                max_point = secp256k1->SubtractPoints(max_point, stride_point);


            }
        };

        std::thread divide_search_Threads[cpuCores];
        
        divide_search_Threads[0] = std::thread(scalable_divide_search_save, range_nums[0], range_pts[0]);
        for (int i = 1; i < cpuCores; i++) {
            divide_search_Threads[i] = std::thread(scalable_divide_search, range_nums[i], range_pts[i] );
        }

        for (int i = 0; i < cpuCores; i++) {
            divide_search_Threads[i].join();
        }

    };

    print_time(); cout << "Search in progress..." << endl;
    
    std::thread search_thread(divide_search);
    
    search_thread.join();
}

/*
    Int max_value, min_value, divide, max_min_stride, between, walker, walker_stride, loop;
    Point max_point, min_point, max_min_stride_point, walker_point, walker_stride_point;
    string temp;
    int counter{0};
    
    ifstream inFile("settings.txt");
    getline(inFile, temp); max_value.SetBase10(trim(const_cast<char*>(temp.c_str()), NULL));
    getline(inFile, temp); min_value.SetBase10(trim(const_cast<char*>(temp.c_str()), NULL));
    getline(inFile, temp); divide.SetBase10(trim(const_cast<char*>(temp.c_str()), NULL));
    inFile.close();
    
    max_point = secp256k1->ComputePublicKey(&max_value);
    min_point = secp256k1->ComputePublicKey(&min_value);
    max_min_stride.SetBase10("1");
    max_min_stride_point = secp256k1->ComputePublicKey(&max_min_stride);
    between.Sub(&max_value, &min_value);

    while(divide.IsLower(&between)) {
        seconds = time(NULL);
        timeStruct = localtime(&seconds);
        printf("[%02d:%02d:%02d] [%d] %s\n", timeStruct->tm_hour, timeStruct->tm_min, timeStruct->tm_sec, counter, divide.GetBase10().c_str());
        // max_point
        bitAddr_U = secp256k1->GetAddress(0, false, max_point);
        if (bloom_filter_check_string(&bfd, bitAddr_U.data()) == BLOOM_SUCCESS) {
            wif_U = secp256k1->GetPrivAddress(false, max_value);
            printf("[--------] Found address: %s %s \n", bitAddr_U.c_str(), wif_U.c_str());
            ofstream outFile;
            outFile.open("found.txt", ios::app);
            outFile << bitAddr_U << " " << wif_U << '\n';
            outFile.close();
            counter++;
        }
        bitAddr_C = secp256k1->GetAddress(0, true, max_point);
        if (bloom_filter_check_string(&bfd, bitAddr_C.data()) == BLOOM_SUCCESS) {
            wif_C = secp256k1->GetPrivAddress(true, max_value);
            printf("[--------] Found address: %s %s \n", bitAddr_C.c_str(), wif_C.c_str());
            ofstream outFile;
            outFile.open("found.txt", ios::app);
            outFile << bitAddr_C << " " << wif_C << '\n';
            outFile.close();
            counter++;
        }
        addrBech32 = secp256k1->GetAddress(2, true, max_point);
        if (bloom_filter_check_string(&bfd, addrBech32.data()) == BLOOM_SUCCESS) {
            wif_C = secp256k1->GetPrivAddress(true, max_value);
            printf("[--------] Found address: %s %s \n", addrBech32.c_str(), wif_C.c_str());
            ofstream outFile;
            outFile.open("found.txt", ios::app);
            outFile << addrBech32 << " " << wif_C << '\n';
            outFile.close();
            counter++;
        }
        // min_point
        bitAddr_U = secp256k1->GetAddress(0, false, min_point);
        if (bloom_filter_check_string(&bfd, bitAddr_U.data()) == BLOOM_SUCCESS) {
            wif_U = secp256k1->GetPrivAddress(false, min_value);
            printf("[--------] Found address: %s %s \n", bitAddr_U.c_str(), wif_U.c_str());
            ofstream outFile;
            outFile.open("found.txt", ios::app);
            outFile << bitAddr_U << " " << wif_U << '\n';
            outFile.close();
            counter++;
        }
        bitAddr_C = secp256k1->GetAddress(0, true, min_point);
        if (bloom_filter_check_string(&bfd, bitAddr_C.data()) == BLOOM_SUCCESS) {
            wif_C = secp256k1->GetPrivAddress(true, min_value);
            printf("[--------] Found address: %s %s \n", bitAddr_C.c_str(), wif_C.c_str());
            ofstream outFile;
            outFile.open("found.txt", ios::app);
            outFile << bitAddr_C << " " << wif_C << '\n';
            outFile.close();
            counter++;
        }
        addrBech32 = secp256k1->GetAddress(2, true, min_point);
        if (bloom_filter_check_string(&bfd, addrBech32.data()) == BLOOM_SUCCESS) {
            wif_C = secp256k1->GetPrivAddress(true, min_value);
            printf("[--------] Found address: %s %s \n", addrBech32.c_str(), wif_C.c_str());
            ofstream outFile;
            outFile.open("found.txt", ios::app);
            outFile << addrBech32 << " " << wif_C << '\n';
            outFile.close();
            counter++;
        }
        // walker_point
        between.Sub(&max_value, &min_value);
        walker.Set(&min_value);
        walker_stride.floor_Div(&between, &divide);
        walker_point = secp256k1->ComputePublicKey(&walker);
        walker_stride_point = secp256k1->ComputePublicKey(&walker_stride);
        for(loop.SetBase10("0"); loop.IsLower(&divide); loop.AddOne()) {
            walker_point = secp256k1->Add2(walker_point, walker_stride_point);
            walker_point.Reduce();
            walker.Add(&walker_stride);
            bitAddr_U = secp256k1->GetAddress(0, false, walker_point);
            if (bloom_filter_check_string(&bfd, bitAddr_U.data()) == BLOOM_SUCCESS) {
                wif_U = secp256k1->GetPrivAddress(false, walker);
                printf("[--------] Found address: %s %s \n", bitAddr_U.c_str(), wif_U.c_str());
                ofstream outFile;
                outFile.open("found.txt", ios::app);
                outFile << bitAddr_U << " " << wif_U << '\n';
                outFile.close();
                counter++;
            }
            bitAddr_C = secp256k1->GetAddress(0, true, walker_point);
            if (bloom_filter_check_string(&bfd, bitAddr_C.data()) == BLOOM_SUCCESS) {
                wif_C = secp256k1->GetPrivAddress(true, walker);
                printf("[--------] Found address: %s %s \n", bitAddr_C.c_str(), wif_C.c_str());
                ofstream outFile;
                outFile.open("found.txt", ios::app);
                outFile << bitAddr_C << " " << wif_C << '\n';
                outFile.close();
                counter++;
            }
            addrBech32 = secp256k1->GetAddress(2, true, walker_point);
            if (bloom_filter_check_string(&bfd, addrBech32.data()) == BLOOM_SUCCESS) {
                wif_C = secp256k1->GetPrivAddress(true, walker);
                printf("[--------] Found address: %s %s \n", addrBech32.c_str(), wif_C.c_str());
                ofstream outFile;
                outFile.open("found.txt", ios::app);
                outFile << addrBech32 << " " << wif_C << '\n';
                outFile.close();
                counter++;
            }
        }
        
        divide.AddOne();
        min_value.Add(&max_min_stride);
        max_value.Sub(&max_min_stride);
        min_point = secp256k1->Add2(min_point, max_min_stride_point);
        min_point.Reduce();
        max_point = secp256k1->Subtract(max_point, max_min_stride_point);
        max_point.Reduce();
        
        ofstream outFile("settings.txt");
        outFile << max_value.GetBase10() << '\n';
        outFile << min_value.GetBase10() << '\n';
        outFile << divide.GetBase10() << '\n';
        outFile.close();
    }
   
    seconds = time(NULL);
    timeStruct = localtime(&seconds);
    printf("[%02d:%02d:%02d] Finished...Found %d targets\n", timeStruct->tm_hour, timeStruct->tm_min, timeStruct->tm_sec, counter);
    cin.clear();
    cin.get();
    return 0;
*/
