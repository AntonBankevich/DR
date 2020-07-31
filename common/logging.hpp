//
// Created by anton on 7/27/20.
//
#pragma once
#include "sys/types.h"
#include "sys/sysinfo.h"
#include <sys/time.h>
#include <sys/resource.h>


class Time {
private:
    timespec start{};
public:
    Time() {
        clock_gettime(CLOCK_MONOTONIC, &start);
    }

    string get() const {
        timespec finish{};
        clock_gettime(CLOCK_MONOTONIC, &finish);
        auto worktime = size_t(double(finish.tv_sec - start.tv_sec) + double(finish.tv_nsec - start.tv_nsec) / 1000000000.0);
        struct sysinfo memInfo;
        sysinfo (&memInfo);
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        std::stringstream ss;
        double mem = size_t(usage.ru_maxrss * 0.001);
        string t = "Mb";
        if (mem > 500) {
            mem = (size_t(mem) / 100) * 0.1;
            t = "Gb";
        }
        ss << worktime / 60 / 60 << ":" << worktime / 60 % 60 << ":" << worktime % 60 << " " << mem << t << " ";
        return ss.str();
    }
};

