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

auto main() -> int {
    
    auto chrono_start = std::chrono::high_resolution_clock::now();    // starting the timer
    Secp256K1* secp256k1 = new Secp256K1(); secp256k1->Init(); // initializing secp256k1 context
    
    Int pk; pk.SetInt32(1); // generating power of two values (2^0..2^256) table
    uint64_t mult = 2;
    vector<Int> S_table;
    for (int i = 0; i < 256; i++)
    {
        S_table.push_back(pk);
        pk.Mult(mult);
    }
    print_time(); cout << "S_table generated" << endl;
    
    string bloomfile = "bloomfile_points.bf";
    string settingsFile = "search_settings.txt";
    
    uint64_t range_start, range_end, block_width, search_power; // block_width = number of elements in the bloomfilter and a stride size to walk the range
    string temp, search_pub;
    ifstream inFile1("settings.txt");
    getline(inFile1, temp); range_start = std::stoull(temp);
    getline(inFile1, temp); range_end = std::stoull(temp);
    getline(inFile1, temp); block_width = std::stoull(temp);
    getline(inFile1, temp); search_pub = trim(temp);
    inFile1.close();
    
    print_time(); cout << "Range Start: " << range_start << " bits" << endl;
    print_time(); cout << "Range End  : " << range_end << " bits" << endl;
    print_time(); cout << "Block Width: 2^" << block_width << endl;
    print_time(); cout << "Search Pub : " << search_pub << endl;
    
    search_power = block_width + 1;
    Int inverse_find; 
    inverse_find.Add(&S_table[range_start], &S_table[range_end]);
    
    Point max_point, min_point, max_min_stride_point;
    Int max_value, min_value, divide;
    
    ifstream inFile2(settingsFile);
    getline(inFile2, temp); max_value.SetBase10(trim(temp).data());
    getline(inFile2, temp); min_value.SetBase10(trim(temp).data());
    getline(inFile2, temp); divide.SetBase10(trim(temp).data());
    inFile2.close();
    
    max_point = secp256k1->ScalarMultiplication(&max_value);
    min_point = secp256k1->ScalarMultiplication(&min_value);

    Int max_min_stride(uint64_t(pow(2, search_power - 1)));
    max_min_stride_point = secp256k1->ScalarMultiplication(&max_min_stride);
    
    print_time(); cout << "Loading Bloomfilter image" << endl;
    filter bf;
    std::ifstream in1(bloomfile, std::ios::binary);
    std::size_t c1;
    in1.read((char*) &c1, sizeof(c1));
    bf.reset(c1); // restore capacity
    boost::span<unsigned char> s1 = bf.array();
    in1.read((char*) s1.data(), s1.size()); // load array
    in1.close();
    
    Int between, walker, walker_stride; between.Sub(&max_value, &min_value);
    Point walker_point, walker_stride_point, P;
    string cpub;
    
    vector<std::string> points, points_match;
    
    ifstream inFile3("points.txt");
    while (getline(inFile3, temp)) {
        points.push_back(trim(temp).data());
    }
    inFile3.close();
    
    ifstream inFile4("points_match.txt");
    while (getline(inFile4, temp)) {
        points_match.push_back(trim(temp).data());
    }
    inFile4.close();
    
    Int loop, scalar;
    uint64_t pow_exp = uint64_t(pow(2, search_power));
    
    while(divide.IsLower(&between)) {
        
        print_time(); cout << divide.GetBase10() << endl;
        
        cpub = secp256k1->GetPublicKeyHex(max_point);
        if (bf.may_contain(cpub)) {
            print_time(); cout << "Bloomfilter hit (Max point). Start searching for key..." << endl;
            scalar = max_value;
            scalar.Sub(uint64_t(pow(2, search_power - 2)));
            for(uint64_t i = 0; i < pow_exp; i++) {
                P = secp256k1->ScalarMultiplication(&scalar);
                for (auto& str : points_match) {
                    if (str.find(secp256k1->GetPublicKeyHex(P)) != std::string::npos) {
                        int ind = str.find(' ');
                        uint64_t power = std::stoull(str.substr(0, ind));
                        string op = str.substr(ind + 1, 1);
                        if (op == "S") { scalar.Add(uint64_t(pow(2, power))); }
                        else { scalar.Sub(uint64_t(pow(2, power))); }
                        Point fP = secp256k1->ScalarMultiplication(&scalar);
                        if (secp256k1->GetPublicKeyHex(fP) == search_pub) {
                            print_time(); cout << "Private key found: " << scalar.GetBase10() << endl;
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }
                        else {
                            scalar.Sub(&inverse_find, &scalar);
                            print_time(); cout << "Private key found: " << scalar.GetBase10() << endl;
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }
                    }
                }
                scalar.AddOne();
            }
            print_time(); cout << "False Positive" << endl;
        }
        
        cpub = secp256k1->GetPublicKeyHex(min_point);
        if (bf.may_contain(cpub)) {
            print_time(); cout << "Bloomfilter hit (Min point). Start searching for key..." << endl;
            scalar = min_value;
            scalar.Sub(uint64_t(pow(2, search_power - 2)));
            for(uint64_t i = 0; i < pow_exp; i++) {
                P = secp256k1->ScalarMultiplication(&scalar);
                for (auto& str : points_match) {
                    if (str.find(secp256k1->GetPublicKeyHex(P)) != std::string::npos) {
                        int ind = str.find(' ');
                        uint64_t power = std::stoull(str.substr(0, ind));
                        string op = str.substr(ind + 1, 1);
                        if (op == "S") { scalar.Add(uint64_t(pow(2, power))); }
                        else { scalar.Sub(uint64_t(pow(2, power))); }
                        Point fP = secp256k1->ScalarMultiplication(&scalar);
                        if (secp256k1->GetPublicKeyHex(fP) == search_pub) {
                            print_time(); cout << "Private key found: " << scalar.GetBase10() << endl;
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }
                        else {
                            scalar.Sub(&inverse_find, &scalar);
                            print_time(); cout << "Private key found: " << scalar.GetBase10() << endl;
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }
                    }
                }
                scalar.AddOne();
            }
            print_time(); cout << "False Positive" << endl;
        }
        
        between.Sub(&max_value, &min_value);
        walker.Set(&min_value);
        walker_stride.floor_Div(&between, &divide);
        walker_point = secp256k1->ScalarMultiplication(&walker);
        walker_stride_point = secp256k1->ScalarMultiplication(&walker_stride);
        
        for(loop.SetBase10("0"); loop.IsLower(&divide); loop.AddOne()) {
            
            walker_point = secp256k1->AddPoints(walker_point, walker_stride_point);
            walker.Add(&walker_stride);
            cpub = secp256k1->GetPublicKeyHex(walker_point);
            if (bf.may_contain(cpub)) {
                print_time(); cout << "Bloomfilter hit (Walker point). Start searching for key..." << endl;
                scalar = walker;
                scalar.Sub(uint64_t(pow(2, search_power - 2)));
                for(uint64_t i = 0; i < pow_exp; i++) {
                    P = secp256k1->ScalarMultiplication(&scalar);
                    for (auto& str : points_match) {
                        if (str.find(secp256k1->GetPublicKeyHex(P)) != std::string::npos) {
                            int ind = str.find(' ');
                            uint64_t power = std::stoull(str.substr(0, ind));
                            string op = str.substr(ind + 1, 1);
                            if (op == "S") { scalar.Add(uint64_t(pow(2, power))); }
                            else { scalar.Sub(uint64_t(pow(2, power))); }
                            Point fP = secp256k1->ScalarMultiplication(&scalar);
                            if (secp256k1->GetPublicKeyHex(fP) == search_pub) {
                                print_time(); cout << "Private key found: " << scalar.GetBase10() << endl;
                                print_elapsed_time(chrono_start);
                                exit(0);
                            }
                            else {
                                scalar.Sub(&inverse_find, &scalar);
                                print_time(); cout << "Private key found: " << scalar.GetBase10() << endl;
                                print_elapsed_time(chrono_start);
                                exit(0);
                            }
                        }
                    }
                    scalar.AddOne();
                }
                print_time(); cout << "False Positive" << endl;
            }
        }
        
        divide.AddOne();
        min_value.Add(&max_min_stride);
        max_value.Sub(&max_min_stride);
        min_point = secp256k1->AddPoints(min_point, max_min_stride_point);
        max_point = secp256k1->SubtractPoints(max_point, max_min_stride_point);
        
        ofstream WriteFile1(settingsFile);
        WriteFile1 << max_value.GetBase10() << '\n';
        WriteFile1 << min_value.GetBase10() << '\n';
        WriteFile1 << divide.GetBase10() << '\n';
        WriteFile1.close();
        
    }
  
}
