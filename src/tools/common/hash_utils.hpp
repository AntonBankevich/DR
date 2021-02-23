#pragma once
typedef unsigned __int128 htype;
#include <iostream>
#include <vector>

template<class Key>
struct alt_hasher {
    size_t operator()( const Key& k ) const;
};

template <>
struct alt_hasher<htype>
{
    size_t operator()( const htype& x ) const
    {
        return (size_t(x) * 31) ^ size_t(x >> 64u);
    }
};

inline std::ostream &operator<<(std::ostream &os, htype val) {
    std::vector<size_t> res;
    while(val != 0) {
        res.push_back(val%10);
        val /= 10;
    }
    for(auto it = res.rbegin(); it != res.rend(); ++it) {
        os << *it;
    }
    return os;
}

inline std::istream &operator>>(std::istream &is, htype &val) {
    val = 0;
    std::string tmp;
    is >> tmp;
    for(char c : tmp) {
        val = val * 10 + c - '0';
    }
    return is;
}

inline void writeHashs(std::ostream &os, const std::vector<htype> &hash_list) {
    os << hash_list.size() << std::endl;
    for(htype h : hash_list) {
        size_t * tmp = reinterpret_cast<size_t*>(&h);
        os << *tmp << " " << *(tmp + 1) << std::endl;
    }
}

inline std::vector<htype> readHashs(std::istream &is) {
    std::vector<htype> result;
    size_t len;
    is >> len;
    size_t a[2];
    for(size_t i = 0; i < len; i++) {
        is >> a[0] >> a[1];
        htype * tmp = reinterpret_cast<htype *>(a);
        result.push_back(*tmp);
    }
    return std::move(result);
}

inline std::string decimal_string(htype n) {
    std::vector<char> res;
    while(n) {
        res.push_back('0' + n%10);
        n = n / 10;
    }
    return std::string(res.rbegin(), res.rend());
}

inline htype string128(const std::string &s) {
    htype res = 0;
    for(char c : s)
        res = res * 10 + c;
    return res;
}