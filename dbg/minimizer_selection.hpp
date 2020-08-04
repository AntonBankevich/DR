//
// Created by anton on 8/3/20.
//

#pragma once
#include "rolling_hash.hpp"
#include "sequences/seqio.hpp"
#include "common/logging.hpp"

template<typename htype>
std::vector<htype> constructMinimizers(logging::Logger &logger, const std::string &reads_file, size_t threads, const RollingHash<htype> &hasher, const size_t w) {
    logger << "Reading reads" << std::endl;
    std::vector<std::vector<htype>> prev;
    prev.resize(threads);
    const size_t buffer_size = 1000000000;
    logger.log("Extracting minimizers");
    size_t total_len = 0;
    size_t read_num = 0;
    std::vector<htype> hash_list;
    io::SeqReader reader(reads_file);
    while(not reader.eof()) {
        size_t tlen = 0;
        logger.log("Starting new round");
        std::vector<Contig> reads;
        reads.reserve(1000000);
        std::vector<std::vector<htype>> hashs;
        hashs.resize(threads);
#pragma omp parallel default(none) shared(hasher, prev, hashs, hash_list, reader, tlen, buffer_size, logger, reads, w, cout)
        {
#pragma omp single
            {
                while (not reader.eof() and tlen < buffer_size) {
                    reads.push_back(reader.read());
                    tlen += reads.back().size();
                    if(reads.back().size() < hasher.k + w - 1) {
                        continue;
                    }
                    size_t index = reads.size() - 1;
                    const Contig &read = reads[index];
#pragma omp task default(none) shared(read, hasher, hashs, w)
                    {
                        MinimizerCalculator<htype> calc(read.seq, hasher, w);
                        std::vector<htype> minimizers(calc.minimizerHashs());
                        if (minimizers.size() > 10) {
                            std::sort(minimizers.begin(), minimizers.end());
                            minimizers.erase(std::unique(minimizers.begin(), minimizers.end()), minimizers.end());
                        }
                        std::vector<htype> & thread_hashs = hashs[omp_get_thread_num()];
                        thread_hashs.insert(thread_hashs.end(), minimizers.begin(), minimizers.end());
                    }
                }
                for(std::vector<htype> & tmp: prev) {
                    hash_list.insert(hash_list.end(), tmp.begin(), tmp.end());
                }
                logger << tlen  << " nucleotides in " << reads.size() <<
                          " sequences were read from disk. Processing in progress  " << std::endl;
            }
        }
        std::swap(hashs, prev);
        total_len += tlen;
        read_num += reads.size();
        reads.clear();
        logger << "Processing finished. Total nucleotides processed: " << total_len << std::endl;
        logger << "Total reads processed: " << read_num << std::endl;
    }
    for(std::vector<htype> & tmp: prev) {
        hash_list.insert(hash_list.end(), tmp.begin(), tmp.end());
    }
    logger << "Finished read processing" << std::endl;
    logger << hash_list.size() << " hashs collected. Starting sorting." << std::endl;
    //    TODO replace with parallel std::sort
    __gnu_parallel::sort(hash_list.begin(), hash_list.end());
    hash_list.erase(std::unique(hash_list.begin(), hash_list.end()), hash_list.end());
    logger << "Finished sorting. Total distinct minimizers: " << hash_list.size() << std::endl;
    return hash_list;
}

