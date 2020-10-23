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
    parameters.projected_element_count = total_size(disjointigs) - hasher.k * disjointigs.size();
    parameters.false_positive_probability = 0.0001;
    VERIFY(!!parameters);
    parameters.compute_optimal_parameters();
    BloomFilter filter(parameters);
    const RollingHash<htype> ehasher = hasher.extensionHash();
    std::function<void(const Sequence &)> task = [&filter, &ehasher](const Sequence & seq) {
        if (seq.size() < ehasher.k) {
            return;
        }
        KWH<htype> kmer(ehasher, seq, 0);
        while (true) {
            filter.insert(kmer.hash());
            if (!kmer.hasNext())
                break;
            kmer = kmer.next();
        }
    };
    logger << "Filling bloom filter with k+1-mers." << std::endl;
    processRecords(disjointigs.begin(), disjointigs.end(), logger, threads, task);
    std::pair<size_t, size_t> bits = filter.count_bits();
    logger << "Filled " << bits.first << " bits out of " << bits.second << std::endl;
    logger << "Finished filling bloom filter. Selecting junctions." << std::endl;
    ParallelRecordCollector<htype> junctions(threads);
    std::function<void(const Sequence &)> junk_task = [&filter, &hasher, &junctions](const Sequence & seq) {
        KWH<htype> kmer(hasher, seq, 0);
        size_t cnt = 0;
        while (true) {
            size_t cnt1 = 0;
            size_t cnt2 = 0;
            for (unsigned char c = 0; c < 4u; c++) {
                cnt1 += filter.contains(kmer.extendRight(c));
                cnt2 += filter.contains(kmer.extendLeft(c));
            }
            if (cnt1 != 1 || cnt2 != 1) {
                cnt += 1;
                junctions.emplace_back(kmer.hash());
            }
            VERIFY(cnt1 <= 4 && cnt2 <= 4);
            if (!kmer.hasNext())
                break;
            kmer = kmer.next();
        }
        if (cnt == 0) {
            junctions.emplace_back(KWH<htype>(hasher, seq, 0).hash());
        }
    };

    processRecords(disjointigs.begin(), disjointigs.end(), logger, threads, junk_task);
    std::vector<htype> res = junctions.collect();
    __gnu_parallel::sort(res.begin(), res.end());
    res.erase(std::unique(res.begin(), res.end()), res.end());
    logger << "Collected " << res.size() << " junctions." << std::endl;
    return res;
}

template<typename htype>
SparseDBG<htype> constructDBG(logging::Logger & logger, const std::vector<htype> &vertices, const std::vector<Sequence> &disjointigs,
                                 const RollingHash<htype> &hasher, size_t threads) {
    logger << "Starting DBG construction." << std::endl;
    SparseDBG<htype> dbg(vertices.begin(), vertices.end(), hasher);
    logger << "Vertices created." << std::endl;
    std::function<void(Sequence &)> edge_filling_task = [&dbg](Sequence & seq) {
        dbg.processRead(seq);
    };
    processRecords(disjointigs.begin(), disjointigs.end(), logger, threads, edge_filling_task);

    logger << "Filled dbg edges. Adding hanging vertices " << std::endl;
    ParallelRecordCollector<std::pair<Vertex<htype>*, Edge<htype> *>> tips(threads);

    std::function<void(std::pair<const htype, Vertex<htype>> &)> task =
            [&tips](std::pair<const htype, Vertex<htype>> & pair) {
                Vertex<htype> &rec = pair.second;
                for (Edge<htype> &edge : rec.getOutgoing()) {
                    if(edge.end() == nullptr) {
                        tips.emplace_back(&rec, &edge);
                    }
                }
                for (Edge<htype> &edge : rec.rc().getOutgoing()) {
                    if(edge.end() == nullptr) {
                        tips.emplace_back(&rec.rc(), &edge);
                    }
                }
            };
    processObjects(dbg.begin(), dbg.end(), logger, threads, task);
    for(std::pair<Vertex<htype>*, Edge<htype> *> edge : tips) {
        Vertex<htype> & vertex = dbg.bindTip(*edge.first, *edge.second);
    }
    logger << "Added " << tips.size() << " hanging vertices" << std::endl;

    logger << "Constructed dbg of size " << dbg.size() << std::endl;
//    dbg.checkConsistency(threads, logger);
//    dbg.printStats(logger);
    logger << "Merging edges " << std::endl;
    mergeAll(logger, dbg, threads);
//    dbg.checkConsistency(threads, logger);
    logger << "Ended merging edges. Resulting size " << dbg.size() << std::endl;
    logger << "Statistics for de Bruijn graph:" << std::endl;
    dbg.printStats(logger);
    return std::move(dbg);
}


