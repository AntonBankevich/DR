//
// Created by anton on 8/3/20.
//
#pragma once

#include <common/bloom_filter.hpp>
#include <common/cl_parser.hpp>
#include "sparse_dbg.hpp"
#include "common/logging.hpp"

template<typename htype>
Sequence buildDisjointig(const Vertex<htype> &rec,
                         const std::vector<Edge<htype>> &path) {
    SequenceBuilder sb;
    sb.append(rec.seq);
    for(const Edge<htype> &e : path) {
        sb.append(e.seq());
    }
    Sequence disjointig = sb.BuildSequence();
    if(disjointig < !disjointig)
        return Sequence{};
    if(rec.inDeg() > 0) {
        disjointig = !(rec.rc().getOutgoing()[0].seq()) + disjointig;
    }
    if(path.back().end()->outDeg() > 0) {
        disjointig = disjointig + path.back().end()->getOutgoing()[0].seq();
    }
    return disjointig;
}

template<typename htype>
std::vector<Sequence> extractDisjointigs(logging::Logger & logger, SparseDBG<htype> &sdbg, size_t threads) {
    logger << "Starting to extract disjointigs." << std::endl;
    std::vector<std::vector<Sequence>> res;
    res.resize(threads);
#pragma omp parallel default(none) shared(logger, sdbg, res)
    {
#pragma omp single
        {
            for(auto it = sdbg.begin(); it != sdbg.end(); ++it) {
                const Vertex<htype> &rec = it->second;
                htype hash = it->first;
                if(!rec.isJunction())
                    continue;
#pragma omp task default(none) shared(sdbg, rec, res, hash, logger)
                {
                    for(const Edge<htype> & edge : rec.getOutgoing()) {
                        VERIFY(edge.end() != nullptr);
                        std::vector<Edge<htype>> path = sdbg.walkForward(rec, edge);
                        Sequence disjointig = buildDisjointig(rec, path);
                        if(!disjointig.empty())
                            res[omp_get_thread_num()].push_back(disjointig);
                    }
                    for(const Edge<htype> & edge : rec.rc().getOutgoing()) {
                        VERIFY(edge.end() != nullptr);
                        std::vector<Edge<htype>> path = sdbg.walkForward(rec.rc(), edge);
                        Sequence disjointig = buildDisjointig(rec.rc(), path);
                        if(!disjointig.empty())
                            res[omp_get_thread_num()].push_back(disjointig);
                    }
                }
            }
        }
    }
    std::vector<Sequence> rres;
    for(const std::vector<Sequence> &r : res) {
        rres.insert(rres.end(), r.begin(), r.end());
    }
    size_t sz = 0;
    for (Sequence & seq: rres) {
        sz += seq.size();
//        logger << seq << std::endl;
    }
    logger << "Finished extracting " << rres.size() << " disjointigs of total size " << sz << std::endl;
    return rres;
}

