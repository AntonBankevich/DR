//
// Created by anton on 7/27/20.
//
#pragma once
#include <omp.h>

template<class T>
class ParallelRecordCollector {
    std::vector<std::vector<T>> recs;
public:
    friend class Iterator;
    class Iterator : public std::iterator<std::forward_iterator_tag, T, size_t,  T*, T&>{
    private:
        ParallelRecordCollector<T> &data;
        size_t row;
        size_t col;
    public:
        Iterator(ParallelRecordCollector<T> &_data, size_t _row = 0, size_t _col = 0) : data(_data), row(_row), col(_col) {
            while(row < data.recs.size() && col == data.recs[row].size()) {
                row += 1;
                col = 0;
            }
        }

        void operator++() {
            col += 1;
            while(row < data.recs.size() && col == data.recs[row].size()) {
                row += 1;
                col = 0;
            }
        }

        T &operator *() {
            return data.recs[row][col];
        }

        bool operator==(const Iterator &other) {
            return row == other.row && col == other.col;
        }
        bool operator!=(const Iterator &other) {
            return row != other.row || col != other.col;
        }

    };
    explicit ParallelRecordCollector(size_t thread_num) : recs(thread_num){
    }

    void add(const T &rec) {
        recs[omp_get_thread_num()].emplace_back(rec);
    }

    template< class... Args >
    void emplace_back( Args&&... args ) {
        recs[omp_get_thread_num()].emplace_back(args...);
    }

    Iterator begin() {
        return Iterator(*this, 0, 0);
    }

    Iterator end() {
        return Iterator(*this, recs.size(), 0);
    }

    size_t size() const {
        size_t res = 0;
        for (const std::vector<T> & row : recs) {
            res += row.size();
        }
        return res;
    }

    std::vector<T> collect() {
        std::vector<T> res;
        for(std::vector<T> &row : recs) {
            res.insert(res.end(), row.begin(), row.end());
            row.clear();
        }
        return std::move(res);
    }
};

