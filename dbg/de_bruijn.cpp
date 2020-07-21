//
// Created by anton on 17.07.2020.
//
#define _GLIBCXX_PARALLEL
#include <iostream>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>
#include <common/dir_utils.hpp>
#include <queue>
#include <omp.h>
#include "common/cl_parser.hpp"
#include "rolling_hash.hpp"
#include<parallel/algorithm>
#include <sstream>


typedef unsigned __int128 htype128;

class Time {
private:
    timespec start{};
public:
    Time() {
        clock_gettime(CLOCK_MONOTONIC, &start);
    }

    string get() {
        timespec finish{};
        clock_gettime(CLOCK_MONOTONIC, &finish);
        size_t worktime = size_t(double(finish.tv_sec - start.tv_sec) + double(finish.tv_nsec - start.tv_nsec) / 1000000000.0);
        std::stringstream ss;
        ss << worktime / 60 / 60 << ":" << worktime / 60 % 60 << ":" << worktime % 60 << " ";
        return ss.str();
    }
};

int main(int argc, char **argv) {
    Time time;
    CLParser parser({"reads=", "output-dir=", "threads=8"}, {"o=output-dir", "t=threads"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    std::cout << "Reading reads" <<std::endl;
    io::SeqReader reader(parser.getValue("reads"));
    std::vector<std::vector<htype128>> hashs;
    size_t threads = std::stoi(parser.getValue("threads"));
    hashs.resize(threads);
    const size_t buffer_size = 1000000000;
    const RollingHash<htype128> hasher(7000, 239);
    omp_set_num_threads(threads);
    std::cout << time.get() << "Extracting minimizers" << std::endl;
    while(not reader.eof()) {
        size_t tlen = 0;
        std::cout << time.get() << "Starting new round" << std::endl;
#pragma omp parallel default(none) shared(hasher, hashs, reader, tlen, buffer_size, std::cout, time)
        {
#pragma omp single
            {
                while (not reader.eof() and tlen < buffer_size) {
                    const Contig read = reader.read();
                    tlen += read.size();
#pragma omp task default(none) shared(read, hasher, hashs, std::cout)
                    {
                        MinimizerCalculator<htype128> calc(read.seq, hasher, 4000);
                        std::vector<htype128> k_mers(calc.minimizers());
                        std::vector<htype128> & thread_hashs = hashs[omp_get_thread_num()];
                        thread_hashs.insert(thread_hashs.end(), k_mers.begin(), k_mers.end());
                    }
                }
                std::cout << time.get() << tlen  << " nucleotides read from disk. Processing in progress" << std::endl;
            }
        }
        std::cout << time.get() << "Processing finished" << std::endl;
    }
    std::cout << time.get() << "Finished read processing" << std::endl;
    std::vector<htype128> hash_list;
    for (auto & hash : hashs) {
        hash_list.insert(hash_list.end(), hash.begin(), hash.end());
    }
    std::cout << time.get() << "Finished merginf results: " << hash_list.size() << " hashs collected. Starting sorting." << std::endl;
              //    TODO replace with parallel std::sort
    __gnu_parallel::sort(hash_list.begin(), hash_list.end());
    hash_list.erase(std::unique(hash_list.begin(), hash_list.end()), hash_list.end());
    std::cout << time.get() << "Finished sorting. Total distinct minimizers: " << hash_list.size() << std::endl;
    return 0;
}