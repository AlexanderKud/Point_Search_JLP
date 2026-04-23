#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <iostream>
#include <ctime>
#include <chrono>
#include <string>

#include "util.h"

using namespace std;

void substr(char *dst, char *src, int position, int length) {
    int c = 0;
    while (c < length) {
      dst[c] = src[position+c-1];
      c++;
    }
    dst[c] = '\0';
}

bool startsWith(const char *pre, const char *str)
{
    size_t lenpre = strlen(pre),
           lenstr = strlen(str);
    return lenstr < lenpre ? false : memcmp(pre, str, lenpre) == 0;
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
    std::cout << "[" << output << "] ";
}

vector<uint64_t> break_down_into_pow10(uint64_t num) {
    vector<uint64_t> nums;
    string stri = to_string(num);
    int num_len = stri.length() - 2;
    for (int pw = num_len; pw >= 0; pw--) { nums.push_back(uint64_t(pow(10, pw))); }
    return nums;
}

void print_elapsed_time(std::chrono::time_point<std::chrono::system_clock> start) {
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = end - start;
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    duration -= hours;
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
    duration -= minutes;
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
    print_time(); cout << "Elapsed time: (" << hours.count() << ")hours (" << minutes.count() << ")minutes (" << seconds.count() << ")seconds\n";
}

void set_bit(unsigned char *bloom, size_t pos) {
    bloom[pos >> 3] |= (1 << (pos & 7));
}

int check_bit(unsigned char *bloom, size_t pos) {
    return (bloom[pos >> 3] >> (pos & 7)) & 1;
}

void save_bloom_filter(const char *filename, unsigned char *bloom, size_t size) {
    FILE *fp = fopen(filename, "wb");
    fwrite(bloom, 1, size, fp);
    fclose(fp);
}

void load_bloom_filter(const char *filename, unsigned char *bloom, size_t size) {
    FILE *fp = fopen(filename, "rb");
    fread(bloom, 1, size, fp);
    fclose(fp);
}

std::string bytesToSize(double bytes, int precision) {
    std::string sizes[] = {"Bytes", "KB", "MB", "GB", "TB"};
    int posttxt = 0;
    while (bytes >= 1024 && posttxt < 4) {
        posttxt++;
        bytes = bytes / 1024;
    }
    return std::to_string(bytes).substr(0, std::to_string(bytes).find(".") + precision + 1) + " " + sizes[posttxt];
}
