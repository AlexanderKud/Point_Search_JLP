#ifndef BLOOM_FILTER_H
#define BLOOM_FILTER_H

#include <vector>
#include <string>
#include <cstdint>
#include <functional>
#include <fstream>

class BloomFilter {
private:
    std::vector<bool> bits;
    uint64_t size;
    int hash_functions;
    
    uint64_t fnv1a_hash(const std::string& str) const;
    uint64_t murmur3_hash(const std::string& str, uint64_t seed) const;
    uint64_t get_hash_position(const std::string& item, int hash_function_index) const;

public:
    BloomFilter(size_t capacity, double error_rate);
    BloomFilter(const std::string& filename);
    ~BloomFilter() = default;
    
    void add(const std::string& item);
    bool contains(const std::string& item) const;
    bool contains_batch(const std::vector<std::string>& items, std::vector<bool>& results) const;
    
    bool save_to_file(const std::string& filename) const;
    bool load_from_file(const std::string& filename);
    
    uint64_t get_size() const { return size; }
    int get_hash_count() const { return hash_functions; }
};
#endif
