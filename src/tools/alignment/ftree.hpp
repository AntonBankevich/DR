//
// Created by anton on 10.01.2020.
//

#pragma once

#include <vector>
#include "common/verify.hpp"
#include "common/output_utils.hpp"

const size_t LOG_MAX_SEQUENCE_LENGTH = 30;

template<class Item>
struct OneValue {
    typedef Item type;
    static Item value(const Item & item) {
        return item;
    }
};

template <class U, class V>
struct PairedValue {
    typedef std::pair<U, V> type;
    struct FirstGetter {
        typedef U type;
        static U value(const PairedValue::type & item) {
            return item.first;
        }
    };
    struct SecondGetter {
        typedef V type;
        static V value(const PairedValue::type & item) {
            return item.second;
        }
    };
};

template <class Item>
class FTree {
private:
    std::vector<Item> content;
    size_t deg2 = 1;
public:
    template<class T>
    friend std::ostream& operator<<(std::ostream&, const FTree<T>&);

    class Iterator {
    private:
        FTree<Item> const * data;
        size_t pos;
        Item val; // val = value[pos - 1]
    public:
        explicit Iterator(FTree<Item> const * _data, size_t i = 0): data(_data), pos(-1) {
            pos = i;
            if (pos == 0) {
                val = Item(0);
            } else {
                val = data->get(pos - 1);
            }
        }

        bool hasNext() const {
            return pos < data->content.size();
        }

        Item next() {
            if ((pos & 1u) == 0) {
                val += data->content[pos];
            } else {
                val += data->content[pos];
                size_t back = (pos & (pos + 1u)) - 1;
                size_t tmp = pos - 1;
                while(tmp != back and tmp != -1) {
                    val -= data->content[tmp];
                    tmp = (tmp & (tmp + 1u)) - 1;
                }
                VERIFY(tmp == back);
            }
            pos += 1;
            return val;
        }
    };

    explicit FTree(std::vector<Item> && values): content(values) {
        VERIFY(values.size() > 0)
        while (deg2 * 2 < values.size())
            deg2 = deg2 * 2;
        for(size_t i = content.size() - 1; i != 0; i -= 1) {
            if ((i & (i + 1)) != 0)
                content[i] = content[i] - content[(i & (i + 1)) - 1];
        }
    }

    explicit FTree(const std::vector<Item> & values): content(values) {
        VERIFY(values.size() > 0)
        while (deg2 * 2 < values.size())
            deg2 = deg2 * 2;
        for(size_t i = content.size() - 1; i != 0; i -= 1) {
            if ((i & (i + 1)) != 0)
                content[i] = content[i] - content[(i & (i + 1)) - 1];
        }
    }

    void push_back(const Item &item) {
        size_t n = content.size() - 1;
        size_t left = (n & (n + 1)) - 1;
        n = ((n + 1) & (n + 2)) - 1;
        Item res = item;
        while (n > left) {
            res += content[n];
            n = (n & (n + 1)) - 1;
        }
        if (deg2 == content.size())
            deg2 = deg2 * 2;
        content.push_back(res + item);
    }

    void extend(const std::vector<Item> & items) {
        for(const Item & item: items) {
            push_back(item);
        }
    }

    Iterator iterator(size_t pos = 0) const {
        return Iterator(this, pos);
    }

//    Find value at index n
//TODO: check if version without ifs works faster
    Item get(size_t n) const {
        VERIFY(n < content.size())
        Item res = 0;
        while (n + 1 != 0) {
            res += content[n];
            n = (n & (n + 1)) - 1;
        }
        return res;
    }

//    val = get(cur - 1); [left,right) overlaps [cur, cur + 2d2)
    void values(std::vector<Item> & result, size_t left, size_t right, size_t cur, size_t d2, Item val) const {
        VERIFY(left <= right and right <= content.size())
        if (d2 == 1) {
            if (left + 1 <= cur + 1)
                result.push_back(val);
            if (right > cur + 1)
                result.push_back(val + content[cur + 1]);
        } else {
            if (left < cur + d2) {// go into [cur , cur + d2)
                values(result, left, right, cur, d2 >> 1u, val);
            }
            if (right > cur + d2) {// go into [cur + d2, cur + 2d2)
                values(result, left, right, cur + d2, d2 >> 1u, val + content[cur + d2]);
            }
        }
    }

//    Get values at indices in [left, right)
    std::vector<Item> get(size_t left, size_t right) const {
        VERIFY(left <= right and right <= content.size())
        std::vector<Item> result;
        result.reserve(right - left + 1);
        values(result, left, right, size_t(-1), deg2, Item(0));
        return result;
    }

    void addToSuffix(size_t pos, Item val) {
        while(pos < size()) {
            content[pos] += val;
            pos = pos | (pos + 1);
        }
    }

    void addToPrefix(size_t pos, Item val) {
        addToSuffix(0, val);
        addToSuffix(pos, -val);
    }

    void addToSegment(size_t from, size_t to, Item val) {
        addToSuffix(from, val);
        addToSuffix(to, -val);
    }

    size_t size() const {
        return content.size();
    }

    //    Use binary search to find rightmost index with prefix sum lesser or equal to item. Result value between -1 and content.size() - 1.
    template<class F = OneValue<Item>>
    std::pair<size_t, typename F::type > binFindLE(const Item &item) const {
        size_t res = 0;
        typename F::type val;
        VERIFY(2 * deg2 > content.size())
//        res - 1: Yes; res + 2d2 - 1: No
        for(size_t d2 = deg2; d2 > 1; d2 = d2 >> 1u) {
            if (res + d2 - 1 >= content.size())
                continue;
            typename F::type new_value = val + F::value(content[res + d2 - 1]);
            if (new_value <= item) {
                res += d2;
                val = new_value;
            }
        }
        return std::make_pair(res - 1, val);
    }
};

template<class Item>
std::ostream& operator<<(std::ostream& out, const FTree<Item>& tree) {
    return out << tree.content << std::endl;
}

