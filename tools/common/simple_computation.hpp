//
// Created by anton on 23.01.2020.
//

#pragma once

template<class T>
T CountSum(T const * from, const T * const to) {
    T res(0);
    while(from != to) {
        res += *from;
        from += 1;
    }
    return res;
}

template<class T>
size_t total_size(const std::vector<T> &data) {
    size_t res = 0;
    for(const T & val : data) {
        res += val.size();
    }
    return res;
}