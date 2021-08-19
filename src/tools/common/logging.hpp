//
// Created by anton on 7/27/20.
//
#pragma once
#include "output_utils.hpp"
#include "dir_utils.hpp"
#include "string_utils.hpp"
#include "sys/types.h"
#include "sys/sysinfo.h"
#include <sys/resource.h>
#include <experimental/filesystem>
#include <ctime>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <iostream>
#include <ostream>
#include <fstream>

namespace logging {

const std::string endl = "\n";

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
        ss << itos(worktime / 60 / 60, 2) << ":" << itos(worktime / 60 % 60, 2) << ":"
                << itos(worktime % 60, 2) << " " << mem << t << " ";
        return ss.str();
    }
};

class LoggerStorage {
private:
    const std::experimental::filesystem::path dir;
    const std::experimental::filesystem::path logFile;
    const std::experimental::filesystem::path backupDir;
public:
    explicit LoggerStorage(std::experimental::filesystem::path _dir, const std::string &_programName):
            dir(std::move(_dir)), logFile(dir / (_programName + ".log")), backupDir(dir / "old_logs") {
    }

    std::experimental::filesystem::path backup() const {
        if(!std::experimental::filesystem::is_regular_file(logFile)) {
            return {};
        }
        ensure_dir_existance(backupDir);
        size_t max = 0;
        for (const std::experimental::filesystem::path & file : std::experimental::filesystem::directory_iterator(backupDir)) {
            std::string fname = file.filename().string();
            if (fname.size() < 5 || fname.substr(fname.size() - 4) != ".log")
                continue;
            try {
                max = std::max<size_t>(max, std::stoi(fname.substr(0, fname.size() - 4)));
            } catch (const std::invalid_argument& ia) {
            }
        }
        std::experimental::filesystem::path backup = backupDir / (itos(max + 1, 2) + ".log");
        std::experimental::filesystem::copy_file(logFile, backup);
        std::experimental::filesystem::remove(logFile);
        return std::move(backup);
    }

    std::experimental::filesystem::path newLoggerFile() const {
        backup();
        return logFile;
    }
};

//info goes to console, trace goes to log file. Debug goes to log file if debug is enabled
//stage like info but for large stage declaration
enum LogLevel {stage, info, trace, debug};

class Logger : public std::streambuf , public std::ostream {
private:
    struct LogStream {
        std::ofstream * os;
        LogLevel level;
        LogStream(const std::experimental::filesystem::path &fn, LogLevel level) :
                os(new std::ofstream()), level(level) {
            os->open(fn);
        }

        ~LogStream() {
            os->close();
            delete os;
            os = nullptr;
        }
    };

    std::vector<LogStream> oss;
    TimeSpace time;
    Logger *empty_logger = nullptr;
    LogLevel curlevel;
    bool add_cout;
public:
    explicit Logger(bool _add_cout = true) :
                std::ostream(this), curlevel(LogLevel::info), add_cout(_add_cout) {
    }

    Logger(const Logger &) = delete;

    void addLogFile(const std::experimental::filesystem::path &fn, LogLevel level = LogLevel::trace) {
        oss.emplace_back(fn, level);
    }

//    template<class T>
//    DummyLogger &operator<<(const T &val) {
//        std::cout << time.get() << val;
//        for(std::ofstream *os : oss) {
//            *os << time.get() << val;
//        }
//        return dummyLogger;
//    }

    int overflow(int c) override {
        if(curlevel <= LogLevel::info)
            std::cout.put(c);
        for(LogStream &os : oss) {
            if(curlevel <= os.level)
                os.os->put(c);
        }
        if(c == '\n') {
            forceFlush();
        }
        return 0;
    }

    void forceFlush() {
        if(curlevel <= LogLevel::info)
            std::cout.flush();
        for(LogStream &os : oss) {
            if(curlevel <= os.level)
                os.os->flush();
        }
    }

    Logger & stage() {
        curlevel = LogLevel::stage;
        *this << "==============================================================\n";
        *this << time.get() << " NEW STAGE: ";
        return *this;
    }

    Logger & info() {
        curlevel = LogLevel::info;
        *this << time.get() << " INFO: ";
        return *this;
    }

    Logger & trace() {
        curlevel = LogLevel::trace;
        *this << time.get() << " TRACE: ";
        return *this;
    }

    Logger & dubug() {
        curlevel = LogLevel::debug;
        *this << time.get() << " TRACE: ";
        return *this;
    }

    ~Logger() override {
        delete empty_logger;
    }
};
}
