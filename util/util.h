#ifndef UTIL_H
#define UTIL_H

#include <vector>

std::vector<std::string> get_files_in_directory(const std::string& directory_path);
void substr(char *dst, char *src, int position, int length);
bool startsWith(const char *pre, const char *str);
std::string trim(const std::string& str);
void print_time();
std::vector<uint64_t> break_down_to_pow10(uint64_t num);
uint64_t str_to_uint64(std::string const& value);

#endif
