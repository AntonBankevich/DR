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
    Sequence disjointig = rec.pathSeq(path);
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
void extractLinearDisjointigs(SparseDBG<htype> &sdbg, std::vector<std::vector<Sequence>> &res) {
#pragma omp parallel default(none) shared(sdbg, res)
    {
#pragma omp single
        {
            for(auto it = sdbg.begin(); it != sdbg.end(); ++it) {
                const Vertex<htype> &rec = it->second;
                htype hash = it->first;
                if(!rec.isJunction())
                    continue;
#pragma omp task default(none) shared(sdbg, rec, res, hash)
                {
                    for(const Edge<htype> & edge : rec.getOutgoing()) {
                        VERIFY(edge.end() != nullptr);
                        std::vector<Edge<htype>> path = sdbg.walkForward(rec, edge);
                        Sequence disjointig = buildDisjointig(rec, path);
                        if(!disjointig.empty()) {
                            for(size_t i = 0; i + 1 < path.size(); i++) {
                                path[i].end()->clearSequence();
                            }
                            res[omp_get_thread_num()].push_back(disjointig);
                        }
                    }
                    for(const Edge<htype> & edge : rec.rc().getOutgoing()) {
                        VERIFY(edge.end() != nullptr);
                        std::vector<Edge<htype>> path = sdbg.walkForward(rec.rc(), edge);
                        Sequence disjointig = buildDisjointig(rec.rc(), path);
                        if(!disjointig.empty()) {
                            for(size_t i = 0; i + 1 < path.size(); i++) {
                                path[i].end()->clearSequence();
                            }
                            res[omp_get_thread_num()].push_back(disjointig);
                        }
                    }
                }
            }
        }
    }
}

template<typename htype>
void extractCircularDisjointigs(SparseDBG<htype> &sdbg, std::vector<std::vector<Sequence>> &res) {
#pragma omp parallel default(none) shared(sdbg, res)
    {
#pragma omp single
        {
            for(auto it = sdbg.begin(); it != sdbg.end(); ++it) {
                Vertex<htype> &rec = it->second;
                htype hash = it->first;
                if(rec.isJunction() || rec.seq.empty())
                    continue;
#pragma omp task default(none) shared(sdbg, rec, res, hash)
                {
                    const Edge<htype> &edge = rec.getOutgoing()[0];
                    VERIFY(edge.end() != nullptr);
                    std::vector<Edge<htype>> path = sdbg.walkForward(rec, edge);
                    VERIFY(path.back().end() == &rec);
                    bool isMinimal = true;
                    for(size_t i = 0; i + 1 < path.size(); i++) {
                        if(path[i].end()->hash() < rec.hash()) {
                            isMinimal = false;
                        }
                    }
                    if(isMinimal) {
                        Sequence tmp = rec.seq;
                        rec.clearSequence();
                        Sequence disjointig = rec.pathSeq(path);
                        res[omp_get_thread_num()].push_back(tmp + disjointig + disjointig);
                    }
                }
            }
        }
    }
}

template<typename htype>
std::vector<Sequence> extractDisjointigs(logging::Logger & logger, SparseDBG<htype> &sdbg, size_t threads) {
    logger << "Starting to extract disjointigs." << std::endl;
    std::vector<std::vector<Sequence>> res;
    res.resize(threads);
    logger << "Extracting linear disjointigs." << std::endl;
    extractLinearDisjointigs(sdbg, res);
    logger << "Finished extracting linear disjointigs." << std::endl;
    logger << "Extracting circular disjointigs." << std::endl;
    extractCircularDisjointigs(sdbg, res);
    logger << "Finished extracting circular disjointigs." << std::endl;
    std::vector<Sequence> rres;
    for(const std::vector<Sequence> &r : res) {
        rres.insert(rres.end(), r.begin(), r.end());
    }
    std::sort(rres.begin(), rres.end(), [] (const Sequence& lhs, const Sequence& rhs) {
        return lhs.size() > rhs.size();
    });
    logger << "Finished extracting " << rres.size() << " disjointigs of total size " << total_size(rres) << std::endl;
    return rres;
}

