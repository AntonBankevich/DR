//
// Created by anton on 7/24/20.
//

#pragma once
typedef unsigned __int128 htype128;

void writeHashs(std::ostream &os, const std::vector<htype128> &hash_list) {
    os << hash_list.size() << std::endl;
    for(htype128 h : hash_list) {
        size_t * tmp = reinterpret_cast<size_t*>(&h);
        os << *tmp << " " << *(tmp + 1) << std::endl;
    }
}

void readHashs(std::istream &is, std::vector<htype128> &hash_list) {
    size_t len;
    is >> len;
    size_t a[2];
    for(size_t i = 0; i < len; i++) {
        is >> a[0] >> a[1];
        htype128 * tmp = reinterpret_cast<htype128 *>(a);
        hash_list.push_back(*tmp);
    }
}

std::string decimal_string(htype128 n) {
    std::vector<char> res;
    while(n) {
        res.push_back('0' + n%10);
        n = n / 10;
    }
    return std::string(res.rbegin(), res.rend());
}