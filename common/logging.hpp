//
// Created by anton on 7/27/20.
//
#pragma once
#include "sys/types.h"
#include "sys/sysinfo.h"
#include <sys/resource.h>
#include "common/output_utils.hpp"
#include <experimental/filesystem>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <ostream>
#include <fstream>

namespace logging {

const std::string endl = "\n";

std::string itos(size_t val) {
    std::stringstream ss;
    if(val < 10)
        ss << "0";
    ss << val;
    return ss.str();
}

class TimeSpace {
private:
    timespec start{};
public:
    TimeSpace() {
        clock_gettime(CLOCK_MONOTONIC, &start);
    }

    std::string get() const {
        timespec finish{};
        clock_gettime(CLOCK_MONOTONIC, &finish);
        auto worktime = size_t(double(finish.tv_sec - start.tv_sec) + double(finish.tv_nsec - start.tv_nsec) / 1000000000.0);
        struct sysinfo memInfo;
        sysinfo (&memInfo);
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        std::stringstream ss;
        double mem = size_t(usage.ru_maxrss * 0.001);
        std::string t = "Mb";
        if (mem > 500) {
            mem = (size_t(mem) / 100) * 0.1;
            t = "Gb";
        }
        ss << itos(worktime / 60 / 60) << ":" << itos(worktime / 60 % 60) << ":" << itos(worktime % 60) << " " << mem << t << " ";
        return ss.str();
    }
};

class Logger {
private:
    std::vector<std::ofstream *> oss;
    TimeSpace time;
    class DummyLogger {
    private:
        Logger & logger_;
    public:
        explicit DummyLogger(Logger & logger): logger_(logger) {}

        template<class T>
        DummyLogger &operator<<(const T &val) {
            std::cout << val;
            for(std::ofstream *os : logger_.oss) {
                *os << val;
            }
            return *this;
        }

        DummyLogger&
        operator<<(std::ostream& (*__pf)(std::ostream &))
        {
            std::cout << std::endl;
            for(std::ofstream *os : logger_.oss) {
                *os << "\n";
            }
            return *this;
        }

    };
    friend class DummyLogger;
    DummyLogger dummyLogger;
public:
    Logger() : dummyLogger(*this) {
    }

    void addLogFile(const std::experimental::filesystem::path &fn) {
        oss.push_back(new std::ofstream());
        oss.back()->open(fn.c_str());
    }

    template<class T>
    DummyLogger &operator<<(const T &val) {
        std::cout << time.get() << val;
        for(std::ofstream *os : oss) {
            *os << time.get() << val;
        }
        return dummyLogger;
    }

    template<class T>
    void log(const T &val) {
        std::cout << time.get() << val << std::endl;
        for(std::ofstream *os : oss) {
            *os << time.get() << val << std::endl;
        }
    }

    void closeAll() {
        for(std::ofstream *os : oss) {
            os->close();
            delete os;
        }
        oss.clear();
    }
    ~Logger() {
        closeAll();
    }
};
//
//    template<class U, class V>
//    std::ostream& operator<<(std::ostream& out, const std::pair<U, V>& item) {
//        return out << "(" << item.first << ", " << item.second << ")";
//    }
//
//    std::ostream& operator<<(std::ostream& out, const unsigned __int128& item) {
//        std::vector<char> res;
//        unsigned __int128 tmp = item;
//        while(tmp != 0) {
//            res.push_back(char((tmp % 10) + '0'));
//            tmp /= 10;
//        }
//        return out << std::string(res.rbegin(), res.rend());
//    }
//
//
//    template<class T>
//    std::ostream& operator<<(std::ostream& out, const std::vector<T>& tree) {
//        if(tree.size() == 0) {
//            return out << "[]" << std::endl;
//        }
//        out << "[";
//        for(size_t i = 0; i + 1 < tree.size(); i += 1) {
//            out << tree[i] << ", ";
//        }
//        return out << tree[tree.size() - 1] << "]";
//    }
}
