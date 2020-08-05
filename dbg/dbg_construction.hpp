//
// Created by anton on 8/3/20.
//
#pragma once

#include "sparse_dbg.hpp"
#include "rolling_hash.hpp"
#include "sequences/sequence.hpp"
#include "common/bloom_filter.hpp"
#include "common/output_utils.hpp"
#include "common/logging.hpp"
#include "common/simple_computation.hpp"
#include "common/omp_utils.hpp"

template<typename htype>
std::vector<htype> findJunctions(logging::Logger & logger, const std::vector<Sequence>& disjointigs,
                                 const RollingHash<htype> &hasher, size_t threads) {
    bloom_parameters parameters;
    parameters.projected_element_count = total_size(disjointigs);
    parameters.false_positive_probability = 0.0001;
    VERIFY(!!parameters);
    parameters.compute_optimal_parameters();
    BloomFilter filter(parameters);
    const RollingHash<htype> ehasher = hasher.extensionHash();
    logger << "Filling bloom filter with k+1-mers." << std::endl;
#pragma omp parallel for default(none) shared(filter, disjointigs, ehasher, hasher, logger)
    for(size_t i = 0; i < disjointigs.size(); i++) {
        const Sequence &seq = disjointigs[i];
        KWH<htype> kmer(ehasher, seq, 0);
//        KWH<htype> kmer1(hasher, seq, 0);
        while(true) {
            filter.insert(kmer.hash());
            if(!kmer.hasNext())
                break;
            kmer = kmer.next();
//            kmer1 = kmer1.next();
//            logger << kmer1.extendRight(seq[kmer1.pos + hasher.k]) << " "
//                      << kmer1.extendLeft(seq[kmer1.pos - 1]) << " " << kmer.hash << std::endl;
        }
    }
    logger << filter.count_bits() << " " << total_size(disjointigs) << std::endl;
    logger << "Finished filling bloom filter. Selecting junctions." << std::endl;
    ParallelRecordCollector<htype> junctions(threads);
    std::vector<size_t> cnt(25);
#pragma omp parallel for default(none) shared(filter, disjointigs, hasher, junctions, logger, cnt)
    for(size_t i = 0; i < disjointigs.size(); i++) {
        const Sequence &seq = disjointigs[i];
        KWH<htype> kmer(hasher, seq, 0);
        while(true) {
            size_t cnt1 = 0;
            size_t cnt2 = 0;
            for(unsigned char c = 0; c < 4u; c++) {
                cnt1 += filter.contains(kmer.extendRight(c));
                cnt2 += filter.contains(kmer.extendLeft(c));
            }
            if(cnt1 != 1 || cnt2 != 1) {
                junctions.emplace_back(kmer.hash());
//                logger << cnt1 << " " << cnt2 << std::endl;
//                logger << kmer.seq.Subseq(kmer.pos, kmer.pos + hasher.k).str() << std::endl;
            }
            VERIFY(cnt1 <= 4 && cnt2 <= 4);
#pragma omp atomic update
            cnt[cnt1 * 5 + cnt2]++;
            if(!kmer.hasNext())
                break;
            kmer = kmer.next();
        }
    }
    logger << cnt << std::endl;

    std::vector<htype> res = junctions.collect();
    __gnu_parallel::sort(res.begin(), res.end());
    res.erase(std::unique(res.begin(), res.end()), res.end());
    logger << "Collected " << res.size() << " junctions." << std::endl;
    return res;
}

template<typename htype>
SparseDBG<htype> constructDBG(logging::Logger & logger, const std::vector<htype> &vertices, const std::vector<Sequence> &disjointigs,
                                 const RollingHash<htype> &hasher) {
    logger << "Starting DBG construction." << std::endl;
    SparseDBG<htype> dbg(vertices, hasher);
#pragma omp parallel for default(none) shared(dbg, disjointigs)
    for(size_t i = 0; i < disjointigs.size(); i++) {
        dbg.processDisjointig(disjointigs[i]);
    }
    logger << "Constructed dbg of size " << dbg.size() << std::endl;
    logger << "Merging edges " << std::endl;
    mergeAll(logger, dbg);
    logger << "Ended merging edges. Resulting size " << dbg.size() << std::endl;
    return std::move(dbg);
}


