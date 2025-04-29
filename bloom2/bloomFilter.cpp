#include "bloomFilter.h"
#include <cmath>
#include <stdexcept>
#include <fstream>

BloomFilter::BloomFilter(uint64_t capacity, double error_rate) {
    size = static_cast<uint64_t>(std::ceil(-(capacity * std::log(error_rate)) / (std::log(2) * std::log(2))));
    hash_functions = static_cast<int>(std::ceil((size / static_cast<double>(capacity)) * std::log(2)));
    bits.resize(size, false);
}

BloomFilter::BloomFilter(const std::string& filename) {
    if (!load_from_file(filename)) {
        throw std::runtime_error("Failed to load Bloom filter from file");
    }
}

uint64_t BloomFilter::fnv1a_hash(const std::string& str) const {
    const uint64_t prime = 1099511628211ULL;
    uint64_t hash = 14695981039346656037ULL;
    
    for (char c : str) {
        hash ^= static_cast<uint8_t>(c);
        hash *= prime;
    }
    
    return hash % size;
}

uint64_t BloomFilter::murmur3_hash(const std::string& str, uint64_t seed) const {
    const uint64_t c1 = 0x87c37b91114253d5ULL;
    const uint64_t c2 = 0x4cf5ad432745937fULL;
    const uint32_t r1 = 31;
    const uint32_t r2 = 27;
    const uint32_t m = 5;
    const uint64_t n = 0x52dce729;

    uint64_t hash = seed;
    const uint8_t* data = reinterpret_cast<const uint8_t*>(str.c_str());
    const int nblocks = str.length() / 8;

    const uint64_t* blocks = reinterpret_cast<const uint64_t*>(data);
    for (int i = 0; i < nblocks; i++) {
        uint64_t k = blocks[i];
        k *= c1;
        k = (k << r1) | (k >> (64 - r1));
        k *= c2;

        hash ^= k;
        hash = ((hash << r2) | (hash >> (64 - r2))) * m + n;
    }

    const uint8_t* tail = data + nblocks * 8;
    uint64_t k1 = 0;

    switch (str.length() & 7) {
        case 7: k1 ^= static_cast<uint64_t>(tail[6]) << 48;
        case 6: k1 ^= static_cast<uint64_t>(tail[5]) << 40;
        case 5: k1 ^= static_cast<uint64_t>(tail[4]) << 32;
        case 4: k1 ^= static_cast<uint64_t>(tail[3]) << 24;
        case 3: k1 ^= static_cast<uint64_t>(tail[2]) << 16;
        case 2: k1 ^= static_cast<uint64_t>(tail[1]) << 8;
        case 1:
            k1 ^= static_cast<uint64_t>(tail[0]);
            k1 *= c1;
            k1 = (k1 << r1) | (k1 >> (64 - r1));
            k1 *= c2;
            hash ^= k1;
    }

    hash ^= str.length();
    hash ^= (hash >> 33);
    hash *= 0xff51afd7ed558ccdULL;
    hash ^= (hash >> 33);
    hash *= 0xc4ceb9fe1a85ec53ULL;
    hash ^= (hash >> 33);

    return hash % size;
}

uint64_t BloomFilter::get_hash_position(const std::string& item, int hash_function_index) const {
    switch (hash_function_index % 2) {
        case 0: return fnv1a_hash(item);
        case 1: return murmur3_hash(item, hash_function_index * 0x9e3779b97f4a7c15ULL);
        default: return 0;
    }
}

void BloomFilter::add(const std::string& item) {
    for (int i = 0; i < hash_functions; ++i) {
        uint64_t position = get_hash_position(item, i);
        bits[position] = true;
    }
}

bool BloomFilter::contains(const std::string& item) const {
    for (int i = 0; i < hash_functions; ++i) {
        uint64_t position = get_hash_position(item, i);
        if (!bits[position]) {
            return false;
        }
    }
    return true;
}

bool BloomFilter::contains_batch(const std::vector<std::string>& items, std::vector<bool>& results) const {
    results.resize(items.size());
    
    for (uint64_t i = 0; i < items.size(); i++) {
        results[i] = true;
        for (int j = 0; j < hash_functions; j++) {
            uint64_t position = get_hash_position(items[i], j);
            if (!bits[position]) {
                results[i] = false;
                break;
            }
        }
    }
    
    return true;
} 

bool BloomFilter::save_to_file(const std::string& filename) const {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) return false;

    // Write metadata
    out.write(reinterpret_cast<const char*>(&size), sizeof(size));
    out.write(reinterpret_cast<const char*>(&hash_functions), sizeof(hash_functions));

    // Write bits - we need to pack them since vector<bool> is special
    for (uint64_t i = 0; i < bits.size(); i += 8) {
        uint8_t byte = 0;
        for (uint64_t j = 0; j < 8 && (i + j) < bits.size(); ++j) {
            if (bits[i + j]) {
                byte |= (1 << j);
            }
        }
        out.write(reinterpret_cast<const char*>(&byte), sizeof(byte));
    }

    return out.good();
}

bool BloomFilter::load_from_file(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) return false;

    // Read metadata
    uint64_t loaded_size;
    int loaded_hash_functions;

    in.read(reinterpret_cast<char*>(&loaded_size), sizeof(loaded_size));
    in.read(reinterpret_cast<char*>(&loaded_hash_functions), sizeof(loaded_hash_functions));

    if (!in.good()) return false;

    // Initialize the filter
    size = loaded_size;
    hash_functions = loaded_hash_functions;
    bits.resize(size);

    // Read bits
    for (uint64_t i = 0; i < bits.size(); i += 8) {
        uint8_t byte;
        in.read(reinterpret_cast<char*>(&byte), sizeof(byte));
        if (!in.good() && !in.eof()) return false;

        for (uint64_t j = 0; j < 8 && (i + j) < bits.size(); ++j) {
            bits[i + j] = (byte & (1 << j)) != 0;
        }
    }

    return true;
}
