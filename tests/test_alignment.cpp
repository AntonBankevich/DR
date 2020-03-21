//
// Created by anton on 23.01.2020.
//

#include "gtest/gtest.h"
#include "../alignment/ftree.hpp"
#include "../common/output_utils.hpp"
#include "../common/simple_computation.hpp"


template<class T>
void assert_str(const T &item, const T &expected_value) {
    std::stringstream ss1, ss2;
    ss1 << item;
    ss2 << expected_value;
    ASSERT_EQ(ss1.str(), ss2.str());
}

template<class T>
void assert_str(const T &item, const std::string &expected_value) {
    std::stringstream ss;
    ss << item;
    ASSERT_EQ(ss.str(), expected_value);
}

TEST(FTreeTests, Construction) {
    std::vector<size_t> val = {1, 2, 5, 6, 50, 100, 200, 201, 202, 204};
    FTree<size_t> tree(val);
    ASSERT_EQ(tree.get(0, tree.size()), val);
}

TEST(FTreeTests, Get) {
    std::vector<size_t> val = {1, 2, 5, 6, 50, 100, 200, 201, 202, 204};
    FTree<size_t> tree(val);
    for(size_t i = 0; i < val.size(); i++) {
        ASSERT_EQ(tree.get(i), val[i]);
    }
    for(size_t i = 0; i < val.size(); i++) {
        for(size_t j = i + 1; j <= val.size(); j++) {
            ASSERT_EQ(tree.get(i, j), std::vector<size_t>(val.data() + i, val.data() + j));
        }
    }
}

TEST(FTreeTests, Iterator) {
    std::vector<size_t> val = {1, 2, 5, 6, 50, 100, 200, 201, 202, 204};
    FTree<size_t> tree(val);
    auto it = tree.iterator();
    for(size_t v : val) {
        ASSERT_EQ(it.next(), v);
    }
    ASSERT_FALSE(it.hasNext());
    for(size_t start = 0; start < val.size(); start++) {
        it = tree.iterator(start);
        for (size_t i = start; i < val.size(); i++) {
            ASSERT_EQ(it.next(), val[i]);
        }
        ASSERT_FALSE(it.hasNext());
    }
}

TEST(FTreeTests, Add) {
    std::vector<size_t> val = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    FTree<size_t> tree(val);
    tree.addToSuffix(3, 1);
    ASSERT_EQ(tree.get(0, tree.size()), std::vector<size_t>({0, 0, 0, 1, 1, 1, 1, 1, 1, 1}));
    tree.addToSuffix(3, -1);
    for(size_t i = 0; i < tree.size(); i++) {
        tree.addToSuffix(i, 1);
    }
    ASSERT_EQ(tree.get(0, tree.size()), std::vector<size_t>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
}

TEST(FTreeTests, Time) {
    std::vector<size_t> val({0});
    for(size_t i = 0; i < 100000; i++)
        for(size_t j = 0; j < 10; j++) {
            val.push_back(val.back() + j);
        }
    clock_t t = clock();
    FTree<size_t> tree(val);
    double work_time = double(clock() - t) / double(CLOCKS_PER_SEC);
    cout << "Testing size " << 1000000 << endl;
    cout << "Construction: " << work_time << endl;
    ASSERT_TRUE(work_time <= 0.006);
    t = clock();
    size_t sum = 0;
    for(size_t i = 0; i < val.size(); i++) {
        sum += tree.get(i);
    }
    work_time = double(clock() - t) / double(CLOCKS_PER_SEC);
    cout << "Random access: " << work_time << endl;
    ASSERT_TRUE(work_time <= 0.01);
    t = clock();
    sum = 0;
    for(auto it = tree.iterator(); it.hasNext();) {
        sum += it.next();
    }
    work_time = double(clock() - t) / double(CLOCKS_PER_SEC);
    cout << "Sequential access: " << work_time << endl;
    ASSERT_TRUE(work_time <= 0.003);
    t = clock();
    for(size_t i = 0; i < val.size(); i++) {
        tree.addToSuffix(i, 1);
    }
    work_time = double(clock() - t) / double(CLOCKS_PER_SEC);
    cout << "Random edit: " << work_time << endl;
    ASSERT_TRUE(work_time <= 0.015);
}