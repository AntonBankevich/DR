#pragma once
#include "sparse_dbg.hpp"
using namespace dbg;

template<class Iterator>
void fillCoverage(SparseDBG &sdbg, logging::Logger &logger, Iterator begin, Iterator end, size_t threads,
                  const RollingHash &hasher, size_t min_read_size);

SparseDBG constructSparseDBGFromReads(logging::Logger & logger, const io::Library &reads_file, size_t threads, const RollingHash &hasher,
                                      const std::vector<htype> &hash_list, size_t w);

void tieTips(logging::Logger &logger, SparseDBG &sdbg, size_t w, size_t threads);

void UpdateVertexTips(Vertex &rec, ParallelRecordCollector<Vertex *> &queue);

void findTips(logging::Logger &logger, SparseDBG &sdbg, size_t threads);

void mergeLoop(Path path);

void MergeEdge(SparseDBG &sdbg, Vertex &start, Edge &edge);

void mergeLinearPaths(logging::Logger & logger, SparseDBG &sdbg, size_t threads);

void mergeCyclicPaths(logging::Logger & logger, SparseDBG &sdbg, size_t threads);

void mergeAll(logging::Logger & logger, SparseDBG &sdbg, size_t threads);

void CalculateCoverage(const std::experimental::filesystem::path &dir, const RollingHash &hasher, const size_t w,
                       const io::Library &lib, size_t threads, logging::Logger &logger, SparseDBG &dbg);

std::experimental::filesystem::path alignLib(logging::Logger &logger, SparseDBG &dbg, const io::Library &align_lib, const RollingHash &hasher,
                                             const size_t w, const std::experimental::filesystem::path &dir, size_t threads);