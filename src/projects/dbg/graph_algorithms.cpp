#include "graph_algorithms.hpp"

template<class Iterator>
void fillCoverage(SparseDBG &sdbg, logging::Logger &logger, Iterator begin, Iterator end, size_t threads,
                  const RollingHash &hasher, const size_t min_read_size) {
    typedef typename Iterator::value_type ContigType;
    logger.info() << "Starting to fill edge coverages" << std::endl;
    ParallelRecordCollector<size_t> lens(threads);
    std::function<void(ContigType &)> task = [&sdbg, &lens, min_read_size](ContigType & contig) {
        Sequence seq = std::move(contig.makeSequence());
        if(seq.size() >= min_read_size) {
            GraphAlignment path = sdbg.align(seq);
            lens.add(path.size());
            for(Segment<Edge> &seg : path) {
                seg.contig().incCov(seg.size());
            }
            path = sdbg.align(!seq);
            for(Segment<Edge> &seg : path) {
                seg.contig().incCov(seg.size());
            }
        }
    };
    processRecords(begin, end, logger, threads, task);
    logger.info() << "Edge coverage calculated." << std::endl;
    std::vector<size_t> lens_distr(1000);
    for(size_t l : lens) {
        lens_distr[std::min(l, lens_distr.size() - 1)] += 1;
    }
}

SparseDBG constructSparseDBGFromReads(logging::Logger &logger, const io::Library &reads_file, size_t threads,
                                      const RollingHash &hasher, const std::vector<htype> &hash_list, const size_t w) {
    logger.info() << "Starting construction of sparse de Bruijn graph" << std::endl;
    SparseDBG sdbg(hash_list.begin(), hash_list.end(), hasher);
    logger.info() << "Vertex map constructed." << std::endl;
    io::SeqReader reader(reads_file, (hasher.k + w) * 20, (hasher.k + w) * 4);
    sdbg.fillSparseDBGEdges(reader.begin(), reader.end(), logger, threads, w + hasher.k - 1);
    return std::move(sdbg);
}

void tieTips(logging::Logger &logger, SparseDBG &sdbg, size_t w, size_t threads) {
    logger.info() << " Collecting tips " << std::endl;
//    TODO reduce memory consumption!! A lot of duplicated k-mer storing
    ParallelRecordCollector<std::pair<Vertex *, Sequence>> old_edges(threads);
    ParallelRecordCollector<Sequence> new_edges(threads);
    ParallelRecordCollector<htype> new_minimizers(threads);
    std::function<void(std::pair<const htype, Vertex> &)> task =
            [&sdbg, &old_edges, &new_minimizers, &new_edges](std::pair<const htype, Vertex> & pair) {
                Vertex &cvertex = pair.second;
                for(auto *vit : {&cvertex, &cvertex.rc()}) {
                    Vertex &vertex = *vit;
                    VERIFY(!vertex.seq.empty());
                    for (const Edge & ext : vertex) {
                        if (ext.end() == nullptr) {
                            Sequence seq = vertex.seq + ext.seq;
                            KWH kwh(sdbg.hasher(), seq, ext.size());
                            new_edges.add(seq);
                            new_minimizers.emplace_back(kwh.hash());
                        } else {
                            old_edges.add({&vertex, ext.seq});
                        }
                    }
                }
                cvertex.clear();
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);
    logger.info() << "Added " << new_minimizers.size() << " artificial minimizers from tips." << std::endl;
    logger.info() << "Collected " << old_edges.size() << " old edges." << std::endl;
    for(auto it = new_minimizers.begin(); it != new_minimizers.end(); ++it) {
        sdbg.addVertex(*it);
    }
    logger.info() << "New minimizers added to sparse graph." << std::endl;
    logger.info() << "Refilling graph with old edges." << std::endl;
    sdbg.refillSparseDBGEdges(old_edges.begin(), old_edges.end(), logger, threads);
    logger.info() << "Filling graph with new edges." << std::endl;
    sdbg.fillSparseDBGEdges(new_edges.begin(), new_edges.end(), logger, threads, sdbg.hasher().k + 1);
    logger.info() << "Finished fixing sparse de Bruijn graph." << std::endl;
}

void UpdateVertexTips(Vertex &rec, ParallelRecordCollector<Vertex *> &queue) {
    bool ok = true;
    for (const Edge &edge : rec) {
        if (edge.getTipSize() == size_t(-1)) {
            edge.updateTipSize();
        }
        if (edge.getTipSize() == size_t(-1)) {
            ok = false;
        }
    }
    if(ok && rec.inDeg() == 1) {
        queue.add(&(rec.rc()[0].end()->rc()));
    }
}

void findTips(logging::Logger &logger, SparseDBG &sdbg, size_t threads) {
    logger.info() << " Finding tips " << std::endl;
//    TODO reduce memory consumption!! A lot of duplicated k-mer storing
    ParallelRecordCollector<Vertex *> queue(threads);
#pragma omp parallel default(none) shared(sdbg, logger, queue)
    {
#pragma omp single
        {
            for (auto &it: sdbg) {
                Vertex &rec = it.second;
                VERIFY_OMP(!rec.seq.empty());
#pragma omp task default(none) shared(sdbg, rec, logger, queue)
                {
                    UpdateVertexTips(rec, queue);
                    UpdateVertexTips(rec.rc(), queue);
                }
            }
        }
    }
    logger.info() << "Found initial tips. Looking for iterative tips" << std::endl;
    size_t cnt = 0;
    while(!queue.empty()) {
        logger.info() << "Iteration " << cnt << ". Queue size " << queue.size() << std::endl;
        std::vector<Vertex *> prev_queue = queue.collectUnique();
        queue.clear();
#pragma omp parallel default(none) shared(sdbg, logger, prev_queue, queue)
        {
#pragma omp single
            {
                for (auto &it: prev_queue) {
                    Vertex &rec = *it;
                    VERIFY_OMP(!rec.seq.empty());
#pragma omp task default(none) shared(sdbg, rec, logger, queue)
                    {
                        UpdateVertexTips(rec, queue);
                    }
                }
            }
        }
    }
    logger.info() << "Tip finding finished" << std::endl;
}

void mergeLoop(Path path) {
    VERIFY(path.start() == path.finish())
    if(path.size() % 2 == 0 && path.getVertex(path.size() / 2) == path.start().rc()) {
        path = path.subPath(0, path.size() / 2);
    }
    Sequence newSeq = path.Seq();
    size_t cov = 0;
    for(const Edge *e : path) {
        if (e->end()->hash() != path.start().hash()) {
            e->end()->mark();
            e->end()->rc().mark();
        }
        cov += e->intCov();
    }
    size_t k = path.start().seq.size();
    Edge &new_edge = path.start().addEdgeLockFree(Edge(&path.start(), &path.finish(), newSeq.Subseq(k)));
    new_edge.incCov(cov - new_edge.intCov());
    Edge &rc_new_edge = path.finish().rc().addEdgeLockFree(Edge(&path.finish().rc(), &path.start().rc(), (!newSeq).Subseq(k)));
    rc_new_edge.incCov(cov - rc_new_edge.intCov());
}

void MergeEdge(SparseDBG &sdbg, Vertex &start, Edge &edge) {
    Path path = edge.walkForward();
    Vertex &end = path.finish().rc();
    if (path.size() > 1 && end.hash() >= start.hash()) {
        VERIFY(start.seq.size() > 0)
        VERIFY(end.seq.size() > 0);
        Sequence newSeq(path.Seq());
        if (start != end)
            end.lock();
        size_t cov = 0;
        for(size_t i = 0; i + 1 < path.size(); i++) {
//            path[i].end()->clear();
            path[i].end()->mark();
            path[i].end()->rc().mark();
            cov += path[i].intCov();
        }
        cov += path.back().intCov();
        Edge &new_edge = start.addEdgeLockFree(Edge(&start, &end.rc(), newSeq.Subseq(start.seq.size())));
        Edge &rc_new_edge = end.addEdgeLockFree(Edge(&end, &start.rc(), (!newSeq).Subseq(start.seq.size())));
        new_edge.incCov(cov - new_edge.intCov());
        rc_new_edge.incCov(cov - rc_new_edge.intCov());
        if (start != end)
            end.unlock();
    }
}

void mergeLinearPaths(logging::Logger &logger, SparseDBG &sdbg, size_t threads) {
    logger.info() << "Merging linear unbranching paths" << std::endl;
    std::function<void(std::pair<const htype, Vertex> &)> task =
            [&sdbg](std::pair<const htype, Vertex> & pair) {
                Vertex &start = pair.second;
                if (!start.isJunction())
                    return;
                start.lock();
                for (Edge &edge: start) {
                    MergeEdge(sdbg, start, edge);
                }
                start.unlock();
                start.rc().lock();
                for (Edge &edge: start.rc()) {
                    MergeEdge(sdbg, start.rc(), edge);
                }
                start.rc().unlock();
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);
    logger.info() << "Finished merging linear unbranching paths" << std::endl;
}

void mergeCyclicPaths(logging::Logger &logger, SparseDBG &sdbg, size_t threads) {
    logger.info() << "Merging cyclic paths" << std::endl;
    ParallelRecordCollector<htype> loops(threads);
    std::function<void(std::pair<const htype, Vertex> &)> task =
            [&sdbg, &loops](std::pair<const htype, Vertex> & pair) {
                Vertex &start = pair.second;
                if(start.isJunction() || start.marked()) {
                    return;
                }
                Path path = start[0].walkForward();
                VERIFY(path.finish() == start);
                bool ismin = true;
                for (const Edge *e : path) {
                    if (e->end()->hash() < start.hash()) {
                        ismin = false;
                        break;
                    }
                }
                if(ismin) {
                    loops.emplace_back(start.hash());
                }
                start.unlock();
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);
    logger.info() << "Found " << loops.size() << " perfect loops" << std::endl;
    for(htype loop : loops) {
        Vertex &start = sdbg.getVertex(loop);
        Path path = start[0].walkForward();
        mergeLoop(path);
    }
    logger.info() << "Finished merging cyclic paths" << std::endl;
}

void mergeAll(logging::Logger &logger, SparseDBG &sdbg, size_t threads) {
    logger.info() << "Merging unbranching paths" << std::endl;
    mergeLinearPaths(logger, sdbg, threads);
//    sdbg.checkConsistency(threads, logger);
    mergeCyclicPaths(logger, sdbg, threads);
//    sdbg.checkConsistency(threads, logger);
    logger.info() << "Removing isolated vertices" << std::endl;
    sdbg.removeMarked();
    logger.info() << "Finished removing isolated vertices" << std::endl;
//    sdbg.checkConsistency(threads, logger);
}

void CalculateCoverage(const std::experimental::filesystem::path &dir, const RollingHash &hasher, const size_t w,
                       const io::Library &lib, size_t threads, logging::Logger &logger, SparseDBG &dbg) {
    logger.info() << "Calculating edge coverage." << std::endl;
    io::SeqReader reader(lib);
    fillCoverage(dbg, logger, reader.begin(), reader.end(), threads, hasher, w + hasher.k - 1);
    std::ofstream os;
    os.open(dir / "coverages.save");
    os << dbg.size() << std::endl;
    for (std::pair<const htype, Vertex> &pair : dbg) {
        Vertex &v = pair.second;
        os << v.hash() << " " << v.outDeg() << " " << v.inDeg() << std::endl;
        for (const Edge &edge : v) {
            os << size_t(edge.seq[0]) << " " << edge.intCov() << std::endl;
        }
        for (const Edge &edge : v.rc()) {
            os << size_t(edge.seq[0]) << " " << edge.intCov() << std::endl;
        }
    }
    dbg.printCoverageStats(logger);
    os.close();
}

std::experimental::filesystem::path
alignLib(logging::Logger &logger, SparseDBG &dbg, const io::Library &align_lib, const RollingHash &hasher,
         const size_t w, const std::experimental::filesystem::path &dir, size_t threads) {
    logger.info() << "Aligning reads" << std::endl;
    ParallelRecordCollector<std::string> alignment_results(threads);
    std::string acgt = "ACGT";

    std::function<void(StringContig &)> task = [&dbg, &alignment_results, &hasher, w, acgt](StringContig & contig) {
        Contig read = contig.makeContig();
        if(read.size() < w + hasher.k - 1)
            return;
        Path path = dbg.align(read.seq).path();
        std::stringstream ss;
        ss << read.id << " " << path.start().hash() << int(path.start().isCanonical()) << " ";
        for (size_t i = 0; i < path.size(); i++) {
            ss << acgt[path[i].seq[0]];
        }
        alignment_results.emplace_back(ss.str());
        Contig rc_read = read.RC();
        Path rc_path = dbg.align(rc_read.seq).path();
        std::stringstream rc_ss;
        rc_ss << rc_read.id << " " << rc_path.start().hash() << int(rc_path.start().isCanonical()) << " ";
        for (size_t i = 0; i < rc_path.size(); i++) {
            rc_ss << acgt[rc_path[i].seq[0]];
        }
        alignment_results.emplace_back(rc_ss.str());
    };
    std::experimental::filesystem::path alignments_file = dir / "alignments.txt";
    std::ofstream os(alignments_file);
    io::SeqReader reader(align_lib);
    processRecords(reader.begin(), reader.end(), logger, threads, task);
    for(std::string & rec : alignment_results) {
        os << rec << "\n";
    }
    os.close();
    logger.info() << "Finished read alignment. Results are in " << (dir / "alignments.txt") << std::endl;
    return alignments_file;
}

