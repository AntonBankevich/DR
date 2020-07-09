//
// Created by anton on 08.07.2020.
//

#pragma once
#include <vector>
#include <functional>
#include <algorithm>

namespace oneline {
    template<class U, class V, class I>
    std::vector<V> map(I begin, I end, std::function<V(U)> f) {
        std::vector<V> result;
        std::for_each(begin, end, [&](const U& param){ result.push_back(f(param));});
        return std::move(result);
    }

    template<class U, class V, class I>
    std::vector<V> initialize(I begin, const I &end) {
        std::vector<V> result;
        std::for_each(begin, end, [&](const U& param){ result.emplace_back(param);});
        return std::move(result);
    }

    template<class U, class I>
    std::vector<U> filter(I begin, I end, std::function<bool(U)> f) {
        std::vector<U> result;
        std::for_each(begin, end, [&](const U& param){ if (f(param)) {result.push_back(param);}});
        return std::move(result);
    }

}