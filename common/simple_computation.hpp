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