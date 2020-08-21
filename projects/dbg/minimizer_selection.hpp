//
// Created by anton on 8/3/20.
//

#pragma once
#include "rolling_hash.hpp"
#include "sequences/seqio.hpp"
#include "common/logging.hpp"
#include "common/omp_utils.hpp"

template<typename htype>
std::vector<htype> constructMinimizers(logging::Logger &logger, const io::Library &reads_file, size_t threads, const RollingHash<htype> &hasher, const size_t w) {
    logger << "Reading reads" << std::endl;
    std::vector<std::vector<htype>> prev;
    prev.resize(threads);
    const size_t buffer_size = 1000000000;
    logger.log("Extracting minimizers");
    size_t min_read_size = hasher.k + w - 1;
    ParallelRecordCollector<htype> hashs(threads);
    std::function<void(StringContig &)> task = [min_read_size, w, &hasher, &hashs](StringContig & contig) {
        Sequence seq = contig.makeCompressedSequence();
        if(seq.size() >= min_read_size) {
            MinimizerCalculator<htype> calc(seq, hasher, w);
            std::vector<htype> minimizers(calc.minimizerHashs());
            if (minimizers.size() > 10) {
                std::sort(minimizers.begin(), minimizers.end());
                minimizers.erase(std::unique(minimizers.begin(), minimizers.end()), minimizers.end());
            }
            hashs.addAll(minimizers.begin(), minimizers.end());
        }
    };
    io::SeqReader reader(reads_file);
    processRecords(reader.begin(), reader.end(), logger, threads, task);

    logger << "Finished read processing" << std::endl;
    logger << hashs.size() << " hashs collected. Starting sorting." << std::endl;
    std::vector<htype> hash_list = hashs.collectUnique();
    //    TODO replace with parallel std::sort
    __gnu_parallel::sort(hash_list.begin(), hash_list.end());
    hash_list.erase(std::unique(hash_list.begin(), hash_list.end()), hash_list.end());
    logger << "Finished sorting. Total distinct minimizers: " << hash_list.size() << std::endl;
    return hash_list;
}

