#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <cstdint>
#include <cmath>

void substr(char *dst, char *src, int position, int length);
bool startsWith(const char *pre, const char *str);
std::string trim(const std::string& str);
void print_time();
std::vector<uint64_t> break_down_into_pow10(uint64_t num);
void print_elapsed_time(std::chrono::time_point<std::chrono::system_clock> start);
void set_bit(unsigned char *bloom, size_t pos);
int check_bit(unsigned char *bloom, size_t pos);
void save_bloom_filter(const char *filename, unsigned char *bloom, size_t size);
void load_bloom_filter(const char *filename, unsigned char *bloom, size_t size);
std::string bytesToSize(double bytes, int precision);

#endif
