//
// Created by anton on 7/27/20.
//
#pragma once
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
        std::stringstream ss;
        ss << worktime / 60 / 60 << ":" << worktime / 60 % 60 << ":" << worktime % 60 << " ";
        return ss.str();
    }
};

