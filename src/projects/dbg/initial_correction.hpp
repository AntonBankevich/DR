#pragma once
#include "compact_path.hpp"

size_t edit_distance(Sequence s1, Sequence s2) {
    size_t left_skip = 0;
    while(left_skip < s1.size() && left_skip < s2.size() && s1[left_skip] == s2[left_skip]) {
        left_skip++;
    }
    s1 = s1.Subseq(left_skip, s1.size());
    s2 = s2.Subseq(left_skip, s2.size());
    size_t right_skip = 0;
    while(right_skip < s1.size() && right_skip < s2.size() && s1[s1.size() - 1 - right_skip] == s2[s2.size() - 1 - right_skip]) {
        right_skip++;
    }
    s1 = s1.Subseq(0, s1.size() - right_skip);
    s2 = s2.Subseq(0, s2.size() - right_skip);
    std::vector<std::vector<size_t>> d(s1.size() + 1, std::vector<size_t>(s2.size() + 1));
    d[0][0] = 0;
    for(unsigned int i = 1; i <= s1.size(); ++i) d[i][0] = i;
    for(unsigned int i = 1; i <= s2.size(); ++i) d[0][i] = i;

    for(unsigned int i = 1; i <= s1.size(); ++i)
        for(unsigned int j = 1; j <= s2.size(); ++j)
            d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
    return d[s1.size()][s2.size()];
}

size_t bestPrefix(const Sequence &s1, const Sequence &s2) {
    std::vector<std::vector<size_t>> d(s1.size() + 1, std::vector<size_t>(s2.size() + 1));
    d[0][0] = 0;
    for(unsigned int i = 1; i <= s1.size(); ++i) d[i][0] = i;
    for(unsigned int i = 1; i <= s2.size(); ++i) d[0][i] = i;

    for(unsigned int i = 1; i <= s1.size(); ++i)
        for(unsigned int j = 1; j <= s2.size(); ++j)
            d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
    size_t res = s2.size();
    for(size_t j = 0; j < s2.size(); j++)
        if(d[s1.size()][j] < d[s1.size()][res])
            res = j;
    return res;
}

size_t tournament(const Sequence &bulge, const std::vector<Sequence> &candidates) {
    size_t winner = 0;
    std::vector<size_t> dists;
    for(size_t i = 0; i < candidates.size(); i++) {
        dists.push_back(edit_distance(bulge, candidates[i]));
        if (dists.back() < dists[winner])
            winner = i;
    }
    size_t max_dist = std::max<size_t>(20, bulge.size() / 100);
    if(dists[winner] > max_dist)
        return -1;
    for(size_t i = 0; i < candidates.size(); i++) {
        if(i != winner) {
            size_t diff = edit_distance(candidates[winner], candidates[i]);
            VERIFY(dists[winner] <= dists[i] + diff);
            VERIFY(dists[i] <= dists[winner] + diff);
            if(dists[i] < max_dist && dists[i] != dists[winner] + diff)
                return -1;
        }
    }
    return winner;
}

std::vector<Path> FindBulgeAlternatives(const Path &path, size_t max_diff) {
    size_t k = path.start().seq.size();
    std::vector<GraphAlignment> als = GraphAlignment(path.start()).allExtensions(max_diff);
    max_diff = std::min(max_diff, path.len());
    std::vector<Path> res;
    Sequence path_seq = path.truncSeq();
    for(GraphAlignment &diff_al : als) {
        size_t path_pos = 0;
        size_t edge_pos = size_t (-1);
        for(size_t i = 0; i < max_diff; i++) {
            GraphAlignment al = diff_al;
            if(i > 0 && al.size() > 0 && al.lastNucl() == path[path_pos].seq[edge_pos])
                continue;
            Sequence seq = path_seq.Subseq(i, path_seq.size());
            al.extend(seq);
            if(al.valid() && al.endClosed() && al.back().contig().end() == &path.finish()){
                res.emplace_back(al.path());
            }
            edge_pos += 1;
            if(edge_pos == path[path_pos].size()) {
                path_pos += 1;
                edge_pos = 0;
            }
        }
    }
    return res;
}

std::unordered_map<Vertex *, size_t> findReachable(Vertex &start, size_t max_dist) {
    std::priority_queue<std::pair<size_t, Vertex *>> queue;
    std::unordered_map<Vertex *, size_t> res;
    queue.emplace(0, &start);
    while(!queue.empty()) {
        std::pair<size_t, Vertex *> next = queue.top();
        queue.pop();
        if(res.find(next.second) == res.end()) {
            res[next.second] = next.first;
            for(Edge &edge : next.second->getOutgoing()) {
                size_t new_len = next.first + edge.size();
                if(new_len <= max_dist) {
                    queue.emplace(new_len, edge.end());
                }
            }
        }
    }
    return std::move(res);
}

std::vector<GraphAlignment> FindPlausibleBulgeAlternatives(const GraphAlignment &path, size_t max_diff, double min_cov) {
    size_t k = path.start().seq.size();
    size_t max_len = path.len() + max_diff;
    std::unordered_map<Vertex *, size_t> reachable = findReachable(path.finish().rc(), max_len);
    std::vector<GraphAlignment> res;
    GraphAlignment alternative(path.start());
    size_t len = 0;
    while(len <= max_len) {
        if(alternative.finish() == path.finish() && len + max_diff > path.len()) {
            res.emplace_back(alternative);
        }
        Edge * next = nullptr;
        for(Edge &edge : path.finish().getOutgoing()) {
            if(edge.getCoverage() >= min_cov && reachable.find(&edge.end()->rc()) != reachable.end() &&
                reachable[&edge.end()->rc()] + edge.size() + len <= max_len) {
                if(next == nullptr)
                    next = &edge;
                else {
                    return {};
                }
            }
        }
        if(next == nullptr)
            break;
        else {
            alternative.push_back(Segment<Edge>(*next, 0, next->size()));
            len += next->size();
        }
    }
    return std::move(res);
}


std::vector<GraphAlignment> FilterAlternatives(logging::Logger &logger1, const GraphAlignment &initial, std::vector<GraphAlignment> &als,
                                            size_t max_diff, double threshold) {
    size_t len = initial.len();
    std::vector<GraphAlignment> res;
    size_t k = initial.getVertex(0).seq.size();
    for(GraphAlignment &al : als) {
        CompactPath cpath(al);
        bool ok = true;
        for(size_t i = 0; i < al.size(); i++) {
            if(al[i].contig().getCoverage() < threshold) {
                ok = false;
                break;
            }
        }
        if(!ok) {
            continue;
        }
        size_t al_len = al.len();
        if(len > al_len + max_diff || al_len > len + max_diff) {
            continue;
        }
        res.emplace_back(al);
    }
    return res;
}

GraphAlignment chooseBulgeCandidate(logging::Logger &logger, std::ostream &out, const GraphAlignment &bulge,
                                    const RecordStorage &reads_storage, const RecordStorage &ref_storage, double threshold,
                                    std::vector<GraphAlignment> &read_alternatives, bool dump) {
    size_t size = bulge.len();
    out << size << " bulge " << bulge.size() << " " << bulge.minCoverage();
    if(dump) {
        logger  << "New bulge from "  << bulge.start().hash() << bulge.start().isCanonical() << " "
                << bulge.start().outDeg() << " " << bulge.start().inDeg() << std::endl
                << "To " << bulge.finish().hash() << bulge.finish().isCanonical() << " "
                << bulge.finish().outDeg() << " " << bulge.finish().inDeg() << std::endl;
        logger << "Record " << bulge.start().hash() << bulge.start().isCanonical() << std::endl;
        logger << reads_storage.getRecord(bulge.start()).str() << std::endl;
    }
    if(dump) {
        logger << "Alternatives" << std::endl;
        for(const auto& it : read_alternatives) {
            logger << CompactPath(it).cpath() << std::endl;
        }
    }
    std::vector<GraphAlignment> read_alternatives_filtered = FilterAlternatives(logger, bulge, read_alternatives,
                                                                                       std::max<size_t>(100, bulge.len() / 100), threshold);
    if(dump) {
        logger << "Filtered alternatives" << std::endl;
        for(const auto& it : read_alternatives) {
            logger << CompactPath(it).cpath() << std::endl;
        }
    }
    const VertexRecord &rec = ref_storage.getRecord(bulge.start());
    size_t genome_support = 0;
    for(const GraphAlignment & al : read_alternatives) {
        if(rec.countStartsWith(CompactPath(al).seq()) > 0) {
            genome_support += 1;
        }
    }
    size_t filtered_genome_support = 0;
    for(const GraphAlignment & al : read_alternatives_filtered) {
        if(rec.countStartsWith(CompactPath(al).seq()) > 0) {
            filtered_genome_support += 1;
        }
    }
    if(dump) {
        logger << read_alternatives.size() << " " << read_alternatives_filtered.size() << std::endl;
        logger << reads_storage.getRecord(bulge.start());
        logger << "Read alternatives" << std::endl;
        for (GraphAlignment &candidate : read_alternatives) {
            logger << CompactPath(candidate) << std::endl;
        }
        logger << "Result " << read_alternatives_filtered.size() << std::endl;
    }
    out << " " << read_alternatives.size() << " " << genome_support << " "
        << read_alternatives_filtered.size() << " " << filtered_genome_support;
    if(read_alternatives_filtered.size() > 1) {
        if(dump)
            logger << "Multiple choice in bulge " << read_alternatives_filtered.size() << std::endl << bulge.truncSeq() << std::endl;
        Sequence old = bulge.truncSeq();
        std::vector<Sequence> candidates;
        for(GraphAlignment &cand : read_alternatives_filtered) {
            if(dump)
                logger << cand.truncSeq() << std::endl;
            candidates.push_back(cand.truncSeq());
        }
        size_t winner = tournament(old, candidates);
        if(dump)
            logger << "Winner found " << winner << std::endl;
        if(winner != size_t(-1)) {
            read_alternatives_filtered = {read_alternatives_filtered[winner]};
        }
    }
    filtered_genome_support = 0;
    for(const GraphAlignment & al : read_alternatives_filtered) {
        if(rec.countStartsWith(CompactPath(al).seq()) > 0) {
            filtered_genome_support += 1;
        }
    }
    out << " " << read_alternatives_filtered.size() << " " << filtered_genome_support;
    out << " " << (read_alternatives_filtered.size() == 1 ? "+" : "-");
    out << std::endl;
    if(read_alternatives_filtered.size() == 1) {
        return std::move(read_alternatives_filtered[0]);
    } else {
        return bulge;
    }
}

class AbstractAlternativeGenerator {
public:
    virtual std::vector<GraphAlignment> generate(const GraphAlignment &path) = 0;
};

class BulgeAlternativeGenerator : public AbstractAlternativeGenerator {
private:
    double threshold;
    const RecordStorage &storage;
public:
    BulgeAlternativeGenerator(const RecordStorage &_storage, double cov_threshold) :
            storage(_storage), threshold(cov_threshold) {

    }

    std::vector<GraphAlignment> generate(const GraphAlignment &path) {
        return storage.getRecord(path.start()).getBulgeAlternatives(path.finish(), threshold);
    }
};

class TipAlternativeGenerator : public AbstractAlternativeGenerator {
private:
    double threshold;
    const RecordStorage &storage;
public:
    TipAlternativeGenerator(const RecordStorage &_storage, double cov_threshold) :
            storage(_storage), threshold(cov_threshold) {
    }

    std::vector<GraphAlignment> generate(const GraphAlignment &path) override {
        return storage.getRecord(path.start()).getTipAlternatives(path.len(), threshold);
    }
};

class AbstractAlternativeFilter {
public:
    virtual void filter(const GraphAlignment &path, std::vector<GraphAlignment> &alternatives) = 0;
};

class DiffAlternativeFilter : AbstractAlternativeFilter {
private:
    double threshold;
    size_t max_diff;
public:
    DiffAlternativeFilter(double cov_threshold, size_t _max_diff) :  threshold(cov_threshold), max_diff(_max_diff) {
    }

    void filter(const GraphAlignment &path, std::vector<GraphAlignment> &alternatives) {
        size_t len = path.len();
        std::vector<GraphAlignment> res;
        size_t k = path.getVertex(0).seq.size();
        for(GraphAlignment &al : alternatives) {
            CompactPath cpath(al);
            bool ok = true;
            for(size_t i = 0; i < al.size(); i++) {
                if(al[i].contig().getCoverage() < threshold) {
                    ok = false;
                    break;
                }
            }
            if(!ok) {
                continue;
            }
            size_t al_len = al.len();
            if(len > al_len + max_diff || al_len > len + max_diff) {
                continue;
            }
            res.emplace_back(std::move(al));
        }
        return std::swap(alternatives, res);
    }
};


GraphAlignment processTip(logging::Logger &logger, std::ostream &out, const GraphAlignment &tip,
                         const RecordStorage &reads_storage, const RecordStorage &ref_storage,
                         double threshold, bool dump = false) {
    size_t size = tip.len();
    out << size << " tip " << tip.size() << " " << tip.minCoverage();
    if(dump) {
        logger << "New tip from " << tip.start().hash() << tip.start().isCanonical() << " "
               << tip.start().outDeg() << " " << tip.start().inDeg() << std::endl
               << "To " << tip.finish().hash() << tip.finish().isCanonical() << " "
               << tip.finish().outDeg() << " " << tip.finish().inDeg() << std::endl;
    }
    std::vector<GraphAlignment> read_alternatives = reads_storage.getRecord(tip.start()).getTipAlternatives(tip.len(), threshold);
    std::vector<GraphAlignment> read_alternatives_filtered =
            FilterAlternatives(logger, tip, read_alternatives, size_t(-1) / 2, threshold);
    const VertexRecord &rec = ref_storage.getRecord(tip.start());
    size_t genome_support = 0;
    for(const GraphAlignment & al : read_alternatives) {
        if(rec.countStartsWith(CompactPath(al).seq()) > 0) {
            genome_support += 1;
        }
    }
    size_t filtered_genome_support = 0;
    for(const GraphAlignment & al : read_alternatives_filtered) {
        if(rec.countStartsWith(CompactPath(al).seq()) > 0) {
            filtered_genome_support += 1;
        }
    }
    if(dump) {
        logger << read_alternatives.size() << " " << read_alternatives_filtered.size() << std::endl;
        logger << reads_storage.getRecord(tip.start());
        logger << "Read alternatives" << std::endl;
        for (GraphAlignment &candidate : read_alternatives) {
            logger << CompactPath(candidate) << std::endl;
        }
    }
    out << " " << read_alternatives.size() << " " << genome_support << " "
        << read_alternatives_filtered.size() << " " << filtered_genome_support;
    if(read_alternatives_filtered.size() > 1) {
        if(dump)
            logger << "Multiple choice for tip " << read_alternatives_filtered.size() << std::endl << tip.truncSeq() << std::endl;
        Sequence old = tip.truncSeq();
        std::vector<Sequence> candidates;
        for(GraphAlignment &cand : read_alternatives_filtered) {
            if(dump)
                logger << cand.truncSeq() << std::endl;
            Sequence candSeq = cand.truncSeq();
            Sequence prefix = candSeq.Subseq(0, bestPrefix(old, candSeq));
            if(dump)
                logger << prefix << std::endl;
            candidates.push_back(prefix);
            cand.cutBack(cand.len() - prefix.size());
        }
        size_t winner = tournament(old, candidates);
        if(dump)
            logger << "Winner found " << winner << std::endl;
        if(winner != size_t(-1)) {
            read_alternatives_filtered = {read_alternatives_filtered[winner]};
        }
    }
    if(dump)
        logger << "Result " << read_alternatives_filtered.size() << std::endl;
    filtered_genome_support = 0;
    for(const GraphAlignment & al : read_alternatives_filtered) {
        if(rec.countStartsWith(CompactPath(al).seq()) > 0) {
            filtered_genome_support += 1;
        }
    }
    out << " " << read_alternatives_filtered.size() << " " << filtered_genome_support;
    out << " " << (read_alternatives_filtered.size() == 1 ? "+" : "-");
    out << std::endl;
    if(read_alternatives_filtered.size() == 1) {
        return std::move(read_alternatives_filtered[0]);
    } else {
        return std::move(tip);
    }
}

size_t correctLowCoveredRegions(logging::Logger &logger, RecordStorage &reads_storage,
                                RecordStorage &ref_storage,
                                const std::experimental::filesystem::path &out_file,
                                double threshold, size_t k, size_t threads, bool dump) {
    ParallelRecordCollector<std::string> results(threads);
    ParallelCounter simple_bulge_cnt(threads);
    ParallelCounter bulge_cnt(threads);
    logger.info() << "Correcting low covered regions in reads" << std::endl;
    if(dump)
        omp_set_num_threads(1);
    size_t max_size = 5 * k;
#pragma omp parallel for default(none) shared(reads_storage, ref_storage, results, threshold, k, max_size, logger, simple_bulge_cnt, bulge_cnt, dump)
    for(size_t read_ind = 0; read_ind < reads_storage.size(); read_ind++) {
        std::stringstream ss;
        AlignedRead &alignedRead = reads_storage[read_ind];
        if(dump)
            logger << "Processing read " << alignedRead.id << std::endl;
        CompactPath &initial_cpath = alignedRead.path;
        GraphAlignment path = initial_cpath.getAlignment();
        GraphAlignment corrected_path(path.start());
        bool corrected = false;
        for(size_t path_pos = 0; path_pos < path.size(); path_pos++) {
            VERIFY_OMP(corrected_path.finish() == path.getVertex(path_pos));
            Edge &edge = path[path_pos].contig();
            if (edge.getCoverage() >= threshold) {
                corrected_path.push_back(path[path_pos]);
                continue;
            }
            size_t step_back = 0;
            size_t step_front = 0;
            size_t size = edge.size();
            while(step_back < corrected_path.size() && corrected_path[corrected_path.size() - step_back - 1].contig().getCoverage() < 10) {
                size += corrected_path[corrected_path.size() - step_back - 1].size();
                step_back += 1;
            }
            while(step_front + path_pos + 1 < path.size() && path[step_front + path_pos + 1].contig().getCoverage() < 10) {
                size += path[step_front + path_pos + 1].size();
                step_front += 1;
            }
            Vertex &start = corrected_path.getVertex(corrected_path.size() - step_back);
            Vertex &end = path.getVertex(path_pos + 1 + step_front);
            GraphAlignment badPath = corrected_path.subPath(corrected_path.size() - step_back, corrected_path.size())
                                            + path.subPath(path_pos, path_pos + 1 + step_front);
            corrected_path.pop_back(step_back);
            if(dump) {
                logger << "Bad read segment " << alignedRead.id << " " << path_pos << " " << step_back << " "
                       << step_front << " " << path.size()
                       << " " << size << " " << edge.getCoverage() << " size " << step_back + step_front + 1
                       << std::endl;
                if (step_back < corrected_path.size()) {
                    logger << "Start stop " << step_back << " "
                           << corrected_path.getVertex(corrected_path.size() - step_back).hash() << " "
                           << corrected_path.getVertex(corrected_path.size() - step_back).isCanonical()
                           << " " << corrected_path[corrected_path.size() - step_back - 1].contig().getCoverage()
                           << corrected_path[corrected_path.size() - step_back - 1].contig().rc().getCoverage()
                           << std::endl;
                    Vertex &tmpv = corrected_path.getVertex(corrected_path.size() - step_back);
                    logger << tmpv.outDeg() << " " << tmpv.inDeg() << std::endl;
                    for (Edge &e : tmpv.getOutgoing())
                        logger << "Edge out " << e.size() << " " << e.getCoverage() << std::endl;
                    for (Edge &e : tmpv.rc().getOutgoing())
                        logger << "Edge in " << e.size() << " " << e.getCoverage() << std::endl;
                }
                if (path_pos + 1 + step_front < path.size()) {
                    logger << "End stop " << step_front << " " << path.getVertex(path_pos + 1 + step_front).hash()
                           << " " << path.getVertex(path_pos + 1 + step_front).isCanonical()
                           << " " << path[path_pos + step_front].contig().getCoverage() << std::endl;
                    Vertex &tmpv = path.getVertex(path_pos + step_front + 1);
                    logger << tmpv.outDeg() << " " << tmpv.inDeg() << std::endl;
                    for (Edge &e : tmpv.getOutgoing())
                        logger << "Edge out " << e.size() << " " << e.getCoverage() << std::endl;
                    for (Edge &e : tmpv.rc().getOutgoing())
                        logger << "Edge in " << e.size() << " " << e.getCoverage() << std::endl;
                }
            }
            if(step_back == corrected_path.size() && step_front == path.size() - path_pos - 1) {
                if(dump)
                    logger << "Whole read has low coverage. Skipping." << std::endl;
                for(const Segment<Edge> &seg : badPath) {
                    corrected_path.push_back(seg);
                }
            } else if(size > 5 * k) {
                if(dump)
                    logger << "Very long read path segment with low coverage. Skipping." << std::endl;
                for(const Segment<Edge> &seg : badPath) {
                    corrected_path.push_back(seg);
                }
            } else if(step_back == corrected_path.size()) {
                if (size < max_size) {
                    if (dump)
                        logger << "Processing incoming tip" << std::endl;
                    GraphAlignment rcBadPath = badPath.RC();
                    GraphAlignment substitution = processTip(logger, ss, rcBadPath, reads_storage, ref_storage,
                                                             threshold, dump);
                    GraphAlignment rcSubstitution = substitution.RC();
                    corrected_path = std::move(rcSubstitution);
                } else {
                    if (dump)
                        logger << "Very long incoming tip" << std::endl;
                }
            } else if(step_front == path.size() - path_pos - 1) {
                if (size < max_size) {
                    if (dump)
                        logger << "Processing outgoing tip" << std::endl;
                    GraphAlignment substitution = processTip(logger, ss, badPath, reads_storage, ref_storage, threshold,
                                                             dump);
                    for (const Segment<Edge> &seg : substitution) {
                        corrected_path.push_back(seg);
                    }
                } else {
                    if (dump)
                        logger << "Processing outgoing tip" << std::endl;
                }
            } else {
                std::vector<GraphAlignment> read_alternatives;
                if(size < max_size)
                    read_alternatives = reads_storage.getRecord(badPath.start()).getBulgeAlternatives(badPath.finish(), threshold);
                else
                    read_alternatives = FindPlausibleBulgeAlternatives(badPath, std::max<size_t>(size / 10, 100), threshold);
                GraphAlignment substitution = chooseBulgeCandidate(logger, ss, badPath, reads_storage, ref_storage, threshold, read_alternatives, dump);
                for (const Segment<Edge> &seg : substitution) {
                    corrected_path.push_back(seg);
                }
                if(badPath.size() == 1 && corrected_path.size() == 1 && badPath[0] != corrected_path[0]) {
                    simple_bulge_cnt += 1;
                }
                bulge_cnt += 1;
            }
            path_pos = path_pos + step_front;
        }
        if(path != corrected_path) {
            if(alignedRead.id == "m64062_190806_063919/124783250/ccs")
                logger << "Corrected read " << alignedRead.id << std::endl;
            reads_storage.reroute(alignedRead, path, corrected_path);
        }
        results.emplace_back(ss.str());
    }
    logger << "Corrected " << simple_bulge_cnt.get() << " simple bulges" << std::endl;
    logger << "Total " << bulge_cnt.get() << " bulges" << std::endl;
    std::ofstream out;
    out.open(out_file);
    size_t res = 0;
    for(std::string &s : results) {
        if(!s.empty() && s.find("+") != size_t(-1))
            res += 1;
        out << s;
    }
    out.close();
    return res;
}

GraphAlignment findAlternative(logging::Logger &logger, std::ostream &out, const GraphAlignment &bulge,
                                   const RecordStorage &reads_storage) {
    std::vector<GraphAlignment> read_alternatives = reads_storage.getRecord(bulge.start()).getBulgeAlternatives(bulge.finish(), 1);
    std::vector<GraphAlignment> read_alternatives_filtered = FilterAlternatives(logger, bulge, read_alternatives,
                                                                                       std::max<size_t>(100, bulge.len() / 100), 1);
    if(read_alternatives_filtered.size() != 2)
        return bulge;
    for(GraphAlignment &al : read_alternatives_filtered) {
        if(al == bulge)
            continue;
        return std::move(al);
    }
    return bulge;
}


size_t collapseBulges(logging::Logger &logger, RecordStorage &reads_storage,
                                RecordStorage &ref_storage,
                                const std::experimental::filesystem::path &out_file,
                                double threshold, size_t k, size_t threads) {
    ParallelRecordCollector<std::string> results(threads);
    ParallelRecordCollector<Edge*> bulge_cnt(threads);
    ParallelRecordCollector<Edge*> collapsable_cnt(threads);
    ParallelRecordCollector<Edge*> genome_cnt(threads);
    ParallelRecordCollector<Edge*> corruption_cnt(threads);
    ParallelRecordCollector<Edge*> heavy_cnt(threads);
    logger.info() << "Collapsing bulges" << std::endl;
#pragma omp parallel for default(none) shared(reads_storage, ref_storage, results, threshold, k, logger, bulge_cnt, genome_cnt, corruption_cnt, collapsable_cnt)
    for(size_t read_ind = 0; read_ind < reads_storage.size(); read_ind++) {
        std::stringstream ss;
        AlignedRead &alignedRead = reads_storage[read_ind];
        CompactPath &initial_cpath = alignedRead.path;
        GraphAlignment path = initial_cpath.getAlignment();
        bool corrected = false;
        for(size_t path_pos = 0; path_pos < path.size(); path_pos++) {
            Edge &edge = path[path_pos].contig();
            if (path[path_pos].left > 0 || path[path_pos].right < path[path_pos].size()) {
                continue;
            }
            Vertex &start = path.getVertex(path_pos);
            Vertex &end = path.getVertex(path_pos + 1);
            if(start.outDeg() != 2 || start.getOutgoing()[0].end() != start.getOutgoing()[1].end()) {
                continue;
            }
            Edge & alt = edge == start.getOutgoing()[0] ? start.getOutgoing()[1] : start.getOutgoing()[0];

            const VertexRecord &rec = ref_storage.getRecord(start);
            if(edge.getCoverage() < 1 || alt.getCoverage() < 1) {
                continue;
            }
            if(edge.getCoverage() > alt.getCoverage()) {
                continue;
            }
            Edge &rcEdge = edge.rc();
            bulge_cnt.emplace_back(&edge);
            bulge_cnt.emplace_back(&rcEdge);
            if(edge.getCoverage() + alt.getCoverage() > 35 || edge.getCoverage() > alt.getCoverage()) {
                continue;
            }
            collapsable_cnt.emplace_back(&edge);
            collapsable_cnt.emplace_back(&rcEdge);
            bool edge_supp = rec.countStartsWith(Sequence(std::vector<char>({char(edge.seq[0])}))) > 0;
            bool alt_supp = rec.countStartsWith(Sequence(std::vector<char>({char(alt.seq[0])}))) > 0;
#pragma omp critical
            {
                logger << edge.size() << " " << edge.getCoverage() << " " << edge_supp << std::endl;
                logger << alt.size() << " " << alt.getCoverage() << " " << alt_supp << std::endl;
            };
            if(edge_supp != alt_supp) {
                genome_cnt.emplace_back(&edge);
                genome_cnt.emplace_back(&rcEdge);
                if(edge_supp) {
                    corruption_cnt.emplace_back(&edge);
                    corruption_cnt.emplace_back(&rcEdge);
                }
            }
            corrected = true;
            path[path_pos] = {alt, 0, alt.size()};
        }
        if(corrected) {
            GraphAlignment path0 = initial_cpath.getAlignment();
            reads_storage.reroute(alignedRead, path0, path);
        }
        results.emplace_back(ss.str());
    }
    size_t bulges = std::unordered_set<Edge*>(bulge_cnt.begin(), bulge_cnt.end()).size();
    size_t collapsable = std::unordered_set<Edge*>(collapsable_cnt.begin(), bulge_cnt.end()).size();
    size_t genome = std::unordered_set<Edge*>(genome_cnt.begin(), bulge_cnt.end()).size();
    size_t corruption = std::unordered_set<Edge*>(corruption_cnt.begin(), bulge_cnt.end()).size();
    logger << "Bulge collapsing results " << bulges << " " << collapsable << " " << genome << " " << corruption << std::endl;
    return bulges;
}

size_t correctAT(logging::Logger &logger, RecordStorage &reads_storage, size_t k, size_t threads) {
    ParallelRecordCollector<std::string> results(threads);
    logger.info() << "Correcting dinucleotide errors in reads" << std::endl;
    ParallelCounter cnt(threads);
#pragma omp parallel for default(none) shared (reads_storage, results, k, logger, cnt, std::cout)
    for(size_t read_ind = 0; read_ind < reads_storage.size(); read_ind++) {
        AlignedRead &alignedRead = reads_storage[read_ind];
        CompactPath &initial_cpath = alignedRead.path;
        GraphAlignment path = initial_cpath.getAlignment();
        size_t corrected = 0;
        for (size_t path_pos = 0; path_pos < path.size(); path_pos++) {
            if(path[path_pos].left > 0 || path[path_pos].right < path[path_pos].contig().size())
                continue;
//            std::cout << alignedRead.id << " " << path_pos << " " << path.size() << std::endl;
            Sequence seq = path.getVertex(path_pos).seq;
            size_t at_cnt1 = 2;
            while (at_cnt1 < seq.size() && seq[seq.size() - at_cnt1 - 1] == seq[seq.size() - at_cnt1 + 1])
                at_cnt1 += 1;
            if (at_cnt1 == seq.size())
                continue;
            if(at_cnt1 < 4)
                continue;
            Sequence extension = path.truncSeq(path_pos, k * 2);
            size_t at_cnt2 = 0;
            while (at_cnt2 < extension.size() && extension[at_cnt2] == seq[seq.size() - 2 + at_cnt2 % 2])
                at_cnt2 += 1;
            if(at_cnt2 >= k)
                continue;
            size_t max_variation = std::max<size_t>(3, (at_cnt1 + at_cnt2) / 50);
            if (extension.size() < k + max_variation * 2 || at_cnt2 > max_variation * 2) {
                continue;
            }
            max_variation = std::min(max_variation, at_cnt1 / 2);
            extension = extension.Subseq(0, k + max_variation * 2);
//            Sequence extension = seq.Subseq(seq.size() - 2 * max_variation, seq.size()) + path.truncSeq(path_pos, k + max_variation * 2);
            SequenceBuilder sb;
            for (size_t i = 0; i < max_variation; i++) {
                sb.append(seq.Subseq(seq.size() - 2, seq.size()));
            }
            sb.append(extension);
            Sequence longest_candidate = sb.BuildSequence();
            Sequence best_seq;
            size_t best_support = 0;
            const VertexRecord &rec = reads_storage.getRecord(path.getVertex(path_pos));
            size_t initial_support = 0;
//            for (size_t i = 0; i <= std::min(2 * max_variation, max_variation + at_cnt2 / 2); i++) {
            size_t step = 2;
            if(seq[seq.size() - 2] == seq[seq.size() - 1])
                step = 1;
            for (size_t skip = 0; skip <= std::min(4 * max_variation, 2 * max_variation + at_cnt2); skip+=step) {
//                size_t skip = i * 2;
                Sequence candidate_seq = longest_candidate.Subseq(skip, skip + k);
                GraphAlignment candidate(path.getVertex(path_pos));
                candidate.extend(candidate_seq);
                if (!candidate.valid())
                    continue;
                CompactPath ccandidate(candidate);
                size_t support = rec.countStartsWith(ccandidate.seq());
                if(skip == 2 * max_variation) {
                    initial_support = support;
                    VERIFY_OMP(support > 0);
                }
                if (support > best_support) {
                    best_seq = longest_candidate.Subseq(skip, longest_candidate.size());
                    best_support = support;
                }
                if (candidate_seq[0] != seq[seq.size() - 2] || candidate_seq[step - 1] != seq[seq.size() - 3 + step])
                    break;
            }
            if (extension.startsWith(best_seq))
                continue;
            logger  << "Correcting ATAT " << best_support << " " << initial_support  << " "
                    << at_cnt1 << " " << at_cnt2 << " " << max_variation << " "
                    << "ACGT"[seq[seq.size() - 2]] << "ACGT"[seq[seq.size() - 1]] << " "
                    << extension.size() << " " << best_seq.size() << std::endl;
            VERIFY_OMP(best_support > 0);
            GraphAlignment rerouting(path.getVertex(path_pos));
            rerouting.extend(best_seq);
            GraphAlignment old_path(path.getVertex(path_pos));
            old_path.extend(extension);
            VERIFY_OMP(old_path.valid());
            VERIFY_OMP(rerouting.valid());
            VERIFY_OMP(rerouting.back() == old_path.back());
            if (rerouting.back().right != rerouting.back().contig().size()) {
                rerouting.pop_back();
                old_path.pop_back();
            }
            GraphAlignment prev_path = path;
            path = path.reroute(path_pos, path_pos + old_path.size(), rerouting.path());
            reads_storage.reroute(alignedRead, prev_path, path);
//            std::cout << "Rerouted " << alignedRead.id << " " << initial_support << " " << best_support << std::endl;
            if(path.size() > 10000) {
                std::cout << extension << "\n" <<best_seq << std::endl;
            }
            corrected += std::max(prev_path.len(), path.len()) - std::min(prev_path.len(), path.len());
        }
        if(corrected > 0) {
#pragma omp critical
            {
                logger << "ATAT " << alignedRead.id << " " << corrected << std::endl;
            }
            ++cnt;
        }
    }
    return cnt.get();
}

void initialCorrect(SparseDBG &sdbg, logging::Logger &logger,
                    const std::experimental::filesystem::path &out_file,
                    const std::experimental::filesystem::path &out_reads,
                    const std::experimental::filesystem::path &bad_reads,
                    const io::Library &reads_lib,
                    const std::experimental::filesystem::path &ref,
                    double threshold, size_t threads, const size_t min_read_size, size_t extension_size, bool dump) {
    size_t k = sdbg.hasher().k;
    logger.info() << "Collecting info from reads" << std::endl;
//    size_t extension_size = std::max(std::min(min_read_size * 3 / 4, sdbg.hasher().k * 11 / 2), sdbg.hasher().k * 3 / 2);
    size_t min_extension = sdbg.hasher().k * 2 / 3;
    RecordStorage reads_storage(sdbg, min_extension, extension_size, true);
    io::SeqReader readReader(reads_lib);
    reads_storage.fill(readReader.begin(), readReader.end(), min_read_size, logger, threads);
    logger.info() << "Collecting info from reference" << std::endl;
    RecordStorage ref_storage(sdbg, min_extension, extension_size, false);
    io::SeqReader refReader(ref);
    ref_storage.fill(refReader.begin(), refReader.end(), min_read_size, logger, threads);
    size_t clow = 0;
    size_t cerr = 0;
    size_t both = 0;
    for(auto &it :sdbg) {
        for(Vertex * vit : {&it.second, &it.second.rc()}) {
            Vertex &v = *vit;
            for(Edge &e : v.getOutgoing()) {
                bool low = e.getCoverage() < threshold;
                bool err = ref_storage.getRecord(v).countStartsWith(Sequence(std::vector<char>({char(e.seq[0])}))) == 0;
                if(low)
                    clow += 1;
                if(err)
                    cerr += 1;
                if(low && err)
                    both += 1;
            }
        }
    }
    logger << "Edge stats " << clow << " " << cerr << " " << both << std::endl;
    {
        size_t correctedAT = correctAT(logger, reads_storage, k, threads);
        logger.info() << "Corrected " << correctedAT << " dinucleotide sequences" << std::endl;
    }
    {
        size_t corrected_low = correctLowCoveredRegions(logger, reads_storage, ref_storage, out_file, threshold, k, threads, dump);
        logger.info() << "Corrected low covered regions in " << corrected_low << " reads" << std::endl;
    }
    size_t corrected_bulges = collapseBulges(logger, reads_storage, ref_storage, out_file, threshold, k, threads);
    logger.info() << "Collapsed bulges in " << corrected_bulges << " reads" << std::endl;
    {
        size_t correctedAT = correctAT(logger, reads_storage, k, threads);
        logger.info() << "Corrected " << correctedAT << " dinucleotide sequences" << std::endl;
    }
    logger.info() << "Printing reads to disk" << std::endl;
    std::ofstream ors;
    std::ofstream brs;
    ors.open(out_reads);
    brs.open(bad_reads);
    for(auto it = reads_storage.begin(); it != reads_storage.end(); ++it) {
        AlignedRead &alignedRead = *it;
        ors << ">" << alignedRead.id << "\n" << alignedRead.path.getAlignment().Seq() << "\n";
        for(auto edge_it : alignedRead.path.getAlignment().path()) {
            if(edge_it->getCoverage() < threshold) {
                brs << ">" << alignedRead.id << "\n" << alignedRead.path.getAlignment().Seq() << "\n";
                break;
            }
        }
    }
    ors.close();
    brs.close();
}