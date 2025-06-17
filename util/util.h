#ifndef UTIL_H
#define UTIL_H

#include <vector>

void substr(char *dst, char *src, int position, int length);
bool startsWith(const char *pre, const char *str);
std::string trim(const std::string& str);
void print_time();
std::vector<uint64_t> break_down_to_pow10(uint64_t num);

#endif
