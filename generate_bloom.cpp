#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <thread>
#include <cmath>

#include "secp256k1/SECP256k1.h"
#include "secp256k1/Int.h"
#include "bloom/filter.hpp"


using namespace std;
namespace fs = filesystem;

vector<string> get_files_in_directory(const string& directory_path) {
    vector<string> files;
    
    try {
        for (const auto& entry : fs::directory_iterator(directory_path)) {
            if (entry.is_regular_file()) {
                files.push_back(entry.path().filename().string());
            }
        }
    } catch (const fs::filesystem_error& e) {
        cerr << "Error accessing directory: " << e.what() << std::endl;
    }
    
    return files;
}

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

auto main() -> int {
    
    auto start = std::chrono::high_resolution_clock::now();
    Secp256K1* secp256k1 = new Secp256K1(); secp256k1->Init();
    
    fs::path current_path = fs::current_path();
    auto file_list = get_files_in_directory(current_path);
    vector<string> targets = {"settings1.txt", "settings2.txt", "bloom1.bf", "bloom2.bf"};
    for (auto i : file_list) {
        for (auto t : targets) {
            if (i == t) { std::remove(t.c_str()); }
        }
    }

    Int pk; pk.SetBase10("1");
    uint64_t mult = 2;
    vector<Point> P_table;
    Point P;
    for (int i = 0; i < 256; i++)
    {
        P = secp256k1->ComputePublicKey(&pk);
        P_table.push_back(P);
        pk.Mult(mult);
    }
    print_time(); cout << "P_table generated" << endl;

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

    Point start_point, end_point, point_05, puzzle_point, puzzle_point_05, puzzle_point_divide2;
    Point first_point, second_point, P1, P2, Q1, Q2;
    Int stride_sum; stride_sum.SetInt32(0);
    start_point = P_table[range_start];
    end_point   = P_table[range_end];
    Int d_05, div2;
    d_05.SetBase10("57896044618658097711785492504343953926418782139537452191302581570759080747169");
    div2.SetBase10("2");
    point_05 = secp256k1->ComputePublicKey(&d_05);
    puzzle_point = secp256k1->ParsePublicKeyHex(search_pub);
    puzzle_point_05 = secp256k1->AddDirect(puzzle_point, point_05);

    puzzle_point_divide2 = secp256k1->PointDivision2(puzzle_point, &div2);

    first_point  = P_table[range_start - 1];
    second_point = P_table[range_start - 2];

    P1 = secp256k1->Subtract(puzzle_point_divide2, first_point);
    P2 = secp256k1->Subtract(puzzle_point_divide2, second_point);
    Q1 = secp256k1->AddDirect(P1, P2);
    Q2 = secp256k1->AddDirect(puzzle_point_divide2, Q1);
    
    ofstream outFile1;
    outFile1.open("settings1.txt", ios::app);
    outFile1 << secp256k1->GetPublicKeyHex(true, Q2) <<'\n';
    outFile1 << stride_sum.GetBase10() << '\n';
    outFile1.close();
    
    ofstream outFile2;
    outFile2.open("settings2.txt", ios::app);
    outFile2 << secp256k1->GetPublicKeyHex(true, Q2) <<'\n';
    outFile2 << stride_sum.GetBase10() << '\n';
    outFile2.close();
    
    print_time(); cout << "Settings written to file" << endl;
    
    using filter = boost::bloom::filter<std::string, 66>;
    
    auto bloom_create1 = [&]() {
        string bloomfile = "bloom1.bf";
        uint64_t n_elements = uint64_t(pow(2, block_width) * 1.5);
        double error = 0.00000001;
        Point P(puzzle_point);
        filter bf(n_elements, error);
        print_time(); cout << "Creating BloomFile1" << '\n';
        string cpub;
        for (int i = 0; i < int(n_elements); i++) {
             bf.insert(secp256k1->GetPublicKeyHex(true, P));
            P = secp256k1->AddDirect(P, secp256k1->G);
        }
        print_time(); cout << "Writing BloomFile1 to bloom1.bf" << '\n';
        std::ofstream out(bloomfile, std::ios::binary);
        std::size_t c1= bf.capacity();
        out.write((const char*) &c1, sizeof(c1)); // save capacity (bits)
        boost::span<const unsigned char> s1 = bf.array();
        out.write((const char*) s1.data(), s1.size()); // save array
        out.close();
    };
    
    auto bloom_create2 = [&]() {
        string bloomfile = "bloom2.bf";
        uint64_t n_elements = uint64_t(pow(2, block_width) * 1.5);
        double error = 0.00000001;
        Point P(puzzle_point_05);
        filter bf(n_elements, error);
        print_time(); cout << "Creating BloomFile2" << '\n';
        for (int i = 0; i < int(n_elements); i++) {
            bf.insert(secp256k1->GetPublicKeyHex(true, P));
            P = secp256k1->AddDirect(P, secp256k1->G);
        }
        print_time(); cout << "Writing BloomFile2 to bloom2.bf" << '\n'; 
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
