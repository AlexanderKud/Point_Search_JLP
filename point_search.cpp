#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <vector>
#include <algorithm>
#include <thread>
#include <cmath>

#include "secp256k1/SECP256k1.h"
#include "secp256k1/Int.h"
#include "bloom/bloom.h"

using namespace std;

inline uint64_t str_to_uint64(std::string const& value) {
  uint64_t result = 0;
  size_t const length = value.size();
  switch (length) {
    case 20:    result += (value[length - 20] - '0') * 10000000000000000000ULL;
    case 19:    result += (value[length - 19] - '0') * 1000000000000000000ULL;
    case 18:    result += (value[length - 18] - '0') * 100000000000000000ULL;
    case 17:    result += (value[length - 17] - '0') * 10000000000000000ULL;
    case 16:    result += (value[length - 16] - '0') * 1000000000000000ULL;
    case 15:    result += (value[length - 15] - '0') * 100000000000000ULL;
    case 14:    result += (value[length - 14] - '0') * 10000000000000ULL;
    case 13:    result += (value[length - 13] - '0') * 1000000000000ULL;
    case 12:    result += (value[length - 12] - '0') * 100000000000ULL;
    case 11:    result += (value[length - 11] - '0') * 10000000000ULL;
    case 10:    result += (value[length - 10] - '0') * 1000000000ULL;
    case  9:    result += (value[length -  9] - '0') * 100000000ULL;
    case  8:    result += (value[length -  8] - '0') * 10000000ULL;
    case  7:    result += (value[length -  7] - '0') * 1000000ULL;
    case  6:    result += (value[length -  6] - '0') * 100000ULL;
    case  5:    result += (value[length -  5] - '0') * 10000ULL;
    case  4:    result += (value[length -  4] - '0') * 1000ULL;
    case  3:    result += (value[length -  3] - '0') * 100ULL;
    case  2:    result += (value[length -  2] - '0') * 10ULL;
    case  1:    result += (value[length -  1] - '0');
  }
  return result;
}

std::string trim(const std::string& str) {
    auto start = str.begin();
    while (start != str.end() && std::isspace(*start)) ++start;
    auto end = str.end();
    do { --end; } while (end != start && std::isspace(*end));
    return std::string(start, end + 1);
}

void print_time() {
    time_t timestamp = time(NULL);
    struct tm datetime = *localtime(&timestamp);
    char output[50];
    strftime(output, 50, "%H:%M:%S", &datetime);
    cout << "[" << output << "] ";
}

vector<uint64_t> break_down_to_pow10(uint64_t num) {
    vector<uint64_t> nums;
    string stri = to_string(num);
    int num_len = stri.length() - 2;
    for (int pw = num_len; pw >= 0; pw--) { nums.push_back(uint64_t(pow(10, pw))); }
    return nums;
}

int main() {

    Secp256K1* secp256k1 = new Secp256K1(); secp256k1->Init();
    
    Int pk; pk.SetBase10("1");
    uint64_t mult = 2;
    vector<Int> S_table;
    for (int i = 0; i < 256; i++)
    {
        S_table.push_back(pk);
        pk.Mult(mult);
    }
    print_time(); cout << "S_table generated" << endl;

    uint64_t range_start, range_end, block_width;
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

    Int first_scalar, second_scalar, pre_calc_sum;
    first_scalar.Set(&S_table[range_start - 1]);
    second_scalar.Set(&S_table[range_start - 2]);
    pre_calc_sum.Add(&first_scalar, &second_scalar);
    
    const char *bloomfile1 = "bloom1.bf";
    const char *bloomfile2 = "bloom2.bf";
    
    print_time(); cout << "Loading Bloomfilter bloom1.bf" << endl;
    BloomFilter bf1;
    bloom_filter_import(&bf1, bloomfile1);
    
    print_time(); cout << "Loading Bloomfilter bloom2.bf" << endl;
    BloomFilter bf2;
    bloom_filter_import(&bf2, bloomfile2);
    
    auto pow10_nums = break_down_to_pow10(uint64_t(pow(2, block_width)));
    vector<Point> pow10_points;
    Int pow_key;
    for (auto n : pow10_nums) {
        pow_key.SetInt64(n);
        pow10_points.push_back(secp256k1->ComputePublicKey(&pow_key));
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    auto addition_search = [&]() {
        uint64_t mult = 2;
        int save_counter = 0;
        string temp, cpub;
        Point starting_point, stride_point, calc_point;
        Int stride_sum, stride;
        ifstream inFile("settings1.txt");
        getline(inFile, temp);
        starting_point = secp256k1->ParsePublicKeyHex2(trim(temp));
        getline(inFile, temp);
        stride_sum.SetBase10(trim(temp).data());
        inFile.close();
        
        stride.SetInt64(uint64_t(pow(2, block_width)));
        stride_point = secp256k1->ComputePublicKey(&stride);
        
        while (true) {
            cpub = secp256k1->GetPublicKeyHex(true, starting_point);
            if (bloom_filter_check_string(&bf1, cpub.data()) == BLOOM_SUCCESS) {
                print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Lower Range Half]" << endl;
                Point P(starting_point);
                vector<uint64_t> privkey_num;
                int index = 0;
                string cpub1;
                for (auto p : pow10_points) {
                    int count = 0;
                    cpub1 = secp256k1->GetPublicKeyHex(true, P);
                    while (bloom_filter_check_string(&bf1, cpub1.data()) == BLOOM_SUCCESS) {
                        P = secp256k1->Subtract(P, p);
                        cpub1 = secp256k1->GetPublicKeyHex(true, P);
                        count += 1;
                    }
                    privkey_num.push_back(pow10_nums[index] * (count - 1));
                    P = secp256k1->Add2(P, p);
                    P.Reduce();
                    index += 1;
                }
                Int Int_steps, Int_temp, privkey;
                uint64_t steps = 0;
                for (auto i : privkey_num) { steps += i; }
                Int_steps.SetInt64(steps);
                Int_temp.Sub(&stride_sum, &Int_steps);
                privkey.Sub(&pre_calc_sum, &Int_temp);
                privkey.Mult(mult);
                calc_point = secp256k1->ComputePublicKey(&privkey);
                if (secp256k1->GetPublicKeyHex(true, calc_point) == search_pub) {
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
            
            if (bloom_filter_check_string(&bf2, cpub.data()) == BLOOM_SUCCESS) {
                print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Lower Range Half]" << endl;
                Point P(starting_point);
                vector<uint64_t> privkey_num;
                int index = 0;
                string cpub1;
                for (auto p : pow10_points) {
                    int count = 0;
                    cpub1 = secp256k1->GetPublicKeyHex(true, P);
                    while (bloom_filter_check_string(&bf2, cpub1.data()) == BLOOM_SUCCESS) {
                        P = secp256k1->Subtract(P, p);
                        cpub1 = secp256k1->GetPublicKeyHex(true, P);
                        count += 1;
                    }
                    privkey_num.push_back(pow10_nums[index] * (count - 1));
                    P = secp256k1->Add2(P, p);
                    P.Reduce();
                    index += 1;
                }
                Int Int_steps, Int_temp, privkey;
                uint64_t steps = 0;
                for (auto i : privkey_num) { steps += i; }
                Int_steps.SetInt64(steps);
                Int_temp.Sub(&stride_sum, &Int_steps);
                privkey.Sub(&pre_calc_sum, &Int_temp);
                privkey.Mult(mult);
                privkey.AddOne();
                calc_point = secp256k1->ComputePublicKey(&privkey);
                if (secp256k1->GetPublicKeyHex(true, calc_point) == search_pub) {
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
            starting_point = secp256k1->Add2(starting_point, stride_point);
            starting_point.Reduce();
            stride_sum.Add(&stride);
            save_counter += 1;
            if (save_counter % 20000000 == 0) {
                cpub = secp256k1->GetPublicKeyHex(true, starting_point);
                ofstream outFile;
                outFile.open("settings1.txt");
                outFile << cpub <<'\n';
                outFile << stride_sum.GetBase10() << '\n';
                outFile.close();
                save_counter = 0;
                print_time(); cout << "Save Data written to settings1.txt" << endl;
            }
        }
    };
    
    auto subtraction_search = [&]() {
        uint64_t mult = 2;
        int save_counter = 0;
        string temp, cpub;
        Point starting_point,stride_point, calc_point;
        Int stride_sum, stride;
        ifstream inFile("settings2.txt");
        getline(inFile, temp);
        starting_point = secp256k1->ParsePublicKeyHex2(trim(temp));
        getline(inFile, temp);
        stride_sum.SetBase10(trim(temp).data());
        inFile.close();
        
        stride.SetInt64(uint64_t(pow(2, block_width)));
        stride_point = secp256k1->ComputePublicKey(&stride);
        
        while (true) {
            cpub = secp256k1->GetPublicKeyHex(true, starting_point);
            if (bloom_filter_check_string(&bf1, cpub.data()) == BLOOM_SUCCESS) {
                print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Higher Range Half]" << endl;
                Point P(starting_point);
                vector<uint64_t> privkey_num;
                int index = 0;
                string cpub1;
                for (auto p : pow10_points) {
                    int count = 0;
                    cpub1 = secp256k1->GetPublicKeyHex(true, P);
                    while (bloom_filter_check_string(&bf1, cpub1.data()) == BLOOM_SUCCESS) {
                        P = secp256k1->Subtract(P, p);
                        cpub1 = secp256k1->GetPublicKeyHex(true, P);
                        count += 1;
                    }
                    privkey_num.push_back(pow10_nums[index] * (count - 1));
                    P = secp256k1->Add2(P, p);
                    P.Reduce();
                    index += 1;
                }
                Int Int_steps, Int_temp, privkey;
                uint64_t steps = 0;
                for (auto i : privkey_num) { steps += i; }
                Int_steps.SetInt64(steps);
                Int_temp.Add(&stride_sum, &Int_steps);
                privkey.Add(&pre_calc_sum, &Int_temp);
                privkey.Mult(mult);
                calc_point = secp256k1->ComputePublicKey(&privkey);
                if (secp256k1->GetPublicKeyHex(true, calc_point) == search_pub) {
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
            
            if (bloom_filter_check_string(&bf2, cpub.data()) == BLOOM_SUCCESS) {
                print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Higher Range Half]" << endl;
                Point P(starting_point);
                vector<uint64_t> privkey_num;
                int index = 0;
                string cpub1;
                for (auto p : pow10_points) {
                    int count = 0;
                    cpub1 = secp256k1->GetPublicKeyHex(true, P);
                    while (bloom_filter_check_string(&bf2, cpub1.data()) == BLOOM_SUCCESS) {
                        P = secp256k1->Subtract(P, p);
                        cpub1 = secp256k1->GetPublicKeyHex(true, P);
                        count += 1;
                    }
                    privkey_num.push_back(pow10_nums[index] * (count - 1));
                    P = secp256k1->Add2(P, p);
                    P.Reduce();
                    index += 1;
                }
                Int Int_steps, Int_temp, privkey;
                uint64_t steps = 0;
                for (auto i : privkey_num) { steps += i; }
                Int_steps.SetInt64(steps);
                Int_temp.Add(&stride_sum, &Int_steps);
                privkey.Add(&pre_calc_sum, &Int_temp);
                privkey.Mult(mult);
                privkey.AddOne();
                calc_point = secp256k1->ComputePublicKey(&privkey);
                if (secp256k1->GetPublicKeyHex(true, calc_point) == search_pub) {
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
            starting_point = secp256k1->Subtract(starting_point, stride_point);
            stride_sum.Add(&stride);
            save_counter += 1;
            if (save_counter % 20000000 == 0) {
                cpub = secp256k1->GetPublicKeyHex(true, starting_point);
                ofstream outFile;
                outFile.open("settings2.txt");
                outFile << cpub <<'\n';
                outFile << stride_sum.GetBase10() << '\n';
                outFile.close();
                save_counter = 0;
                print_time(); cout << "Save Data written to settings2.txt" << endl;
            }
        }
    };
    
    print_time(); cout << "Search in progress..." << endl;
    
    std::thread thread1(addition_search);
    std::thread thread2(subtraction_search);
    
    thread1.join();
    thread2.join();
}
