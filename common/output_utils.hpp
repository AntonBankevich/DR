//
// Created by anton on 23.01.2020.
//

#pragma once
using std::cout;
using std::endl;

template<class U, class V>
std::ostream& operator<<(std::ostream& out, const std::pair<U, V>& item) {
    return out << "(" << item.first << ", " << item.second << ")";
}

template<class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& tree) {
    if(tree.size() == 0) {
        return out << "[]" << std::endl;
    }
    out << "[";
    for(size_t i = 0; i + 1 < tree.size(); i += 1) {
        out << tree[i] << ", ";
    }
    return out << tree[tree.size() - 1] << "]";
}