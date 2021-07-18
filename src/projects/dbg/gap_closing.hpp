#pragma once

//#include "sparse_dbg.hpp"
//#include <sequences/sequence.hpp>

//std::vector<Sequence> CloseGaps(dbg::SparseDBG &dbg, logging::Logger &logger, size_t threads, size_t min_len, size_t smallK, double allowed_divergence) {
//    logger.info() << "Started gap closing procedure" << std::endl;
//    std::vector<dbg::Edge *> tips;
//    for(dbg::Edge &edge : dbg.edges()) {
//        if(edge.size() > min_len + dbg.hasher().getK() && edge.start()->inDeg() == 0 && edge.start()->outDeg() == 1)
//            tips.emplace_back(&edge);
//    }
//    ParallelRecordCollector<std::pair<htype, dbg::EdgePosition>> candidates(threads);
//    RollingHash smallHasher(smallK, 239);
////#pragma omp parallel for default(none) shared(tips, candidates, dbg, smallHasher)
//    for(size_t i = 0; i < tips.size(); i++) {
//        KWH kwh(smallHasher, tips[i]->seq, )
//        for()
//    }
//    return {};
//}