#include <string.h>
#include <ctime>
#include <string>
#include <iostream>
#include <cmath>
#include <cstdint>

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

vector<uint64_t> break_down_to_pow10(uint64_t num) {
    vector<uint64_t> nums;
    string stri = to_string(num);
    int num_len = stri.length() - 2;
    for (int pw = num_len; pw >= 0; pw--) { nums.push_back(uint64_t(pow(10, pw))); }
    return nums;
}
