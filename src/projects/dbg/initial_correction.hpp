#pragma once
#include "compact_path.hpp"
#include "multiplicity_estimation.hpp"

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
    std::vector<size_t> prev(s2.size() + 1);
    std::vector<size_t> cur(s2.size() + 1);
    for(unsigned int i = 0; i <= s1.size(); ++i) cur[i] = i;
    for(unsigned int i = 1; i <= s1.size(); ++i) {
        std::swap(prev, cur);
        cur[0] = i;
        for(unsigned int j = 1; j <= s2.size(); ++j)
            cur[j] = std::min({ prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
    }
    size_t res = s2.size();
    for(size_t j = 0; j < s2.size(); j++)
        if(cur[j] < cur[res])
            res = j;
    return res;
}

size_t tournament(const Sequence &bulge, const std::vector<Sequence> &candidates, bool dump = false) {
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

std::unordered_map<Vertex *, size_t> findReachable(Vertex &start, double min_cov, size_t max_dist) {
    typedef std::pair<size_t, Vertex*> StoredValue;
    std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
    std::unordered_map<Vertex *, size_t> res;
    queue.emplace(0, &start);
    while(!queue.empty()) {
        StoredValue next = queue.top();
        queue.pop();
        if(res.find(next.second) == res.end()) {
            res[next.second] = next.first;
            for(Edge &edge : *next.second) {
                size_t new_len = next.first + edge.size();
                if((edge.getCoverage() >= min_cov || edge.is_reliable) && new_len <= max_dist) {
                    queue.emplace(new_len, edge.end());
                }
            }
        }
    }
    return std::move(res);
}

std::vector<GraphAlignment> FindPlausibleBulgeAlternatives(logging::Logger &logger, const GraphAlignment &path,
                                                           size_t max_diff, double min_cov) {
    logger << "Looking for plausible alternative bulge paths" << std::endl;
    size_t k = path.start().seq.size();
    size_t max_len = path.len() + max_diff;
    std::unordered_map<Vertex *, size_t> reachable = findReachable(path.finish().rc(), min_cov, max_len);
    logger << "Calculated reachable vertices " << path.len() << std::endl;
    if(reachable.find(&path.start().rc()) != reachable.end()) {
        logger << "Calculated reachable vertices. Start is reachable using " << reachable[&path.start().rc()] << std::endl;
    } else {
        logger << "Start is not reachable from end" << std::endl;
    }
    std::vector<GraphAlignment> res;
    GraphAlignment alternative(path.start());
    size_t iter_cnt = 0;
    size_t len = 0;
    bool forward = true;
    while(true) {
        iter_cnt += 1;
        if(iter_cnt > 10000)
            return {path};
        if(forward) {
            if(alternative.finish() == path.finish() && len + max_diff >= path.len()) {
                res.emplace_back(alternative);
                logger << "Found plausible alternative path " << alternative.len() << std::endl;
                if(res.size() > 30) {
                    logger << "Too many plausible alternatives. Aborting." << std::endl;
                    return {path};
                }
            }
            forward = false;
            for(Edge &edge : alternative.finish()) {
                if((edge.getCoverage() >= min_cov || edge.is_reliable) && reachable.find(&edge.end()->rc()) != reachable.end() &&
                   reachable[&edge.end()->rc()] + edge.size() + len <= max_len) {
                    len += edge.size();
                    alternative.push_back(Segment<Edge>(edge, 0, edge.size()));
                    forward = true;
                    break;
                }
            }
        } else {
            if (alternative.size() == 0)
                break;
            Edge &old_edge = alternative.back().contig();
            alternative.pop_back();
            len -= old_edge.size();
            bool found = false;
            for(Edge &edge : alternative.finish()) {
                if((edge.getCoverage() >= min_cov || edge.is_reliable) &&
                        reachable.find(&edge.end()->rc()) != reachable.end() &&
                        reachable[&edge.end()->rc()] + edge.size() + len <= max_len) {
                    if(found) {
                        len += edge.size();
                        alternative.push_back(Segment<Edge>(edge, 0, edge.size()));
                        forward = true;
                        break;
                    } else if (&edge == &old_edge) {
                        found - true;
                    }
                }
            }
        }
    }
    logger << "Found " << res.size() << " plausible alternative paths" << std::endl;
    return std::move(res);
}

std::vector<GraphAlignment> FindPlausibleTipAlternatives(logging::Logger &logger, const GraphAlignment &path,
                                                           size_t max_diff, double min_cov) {
    logger << "Looking for plausible alternative tip paths" << std::endl;
    size_t k = path.start().seq.size();
    size_t max_len = path.len() + max_diff;
    std::vector<GraphAlignment> res;
    GraphAlignment alternative(path.start());
    size_t iter_cnt = 0;
    size_t len = 0;
    size_t tip_len = path.len();
    bool forward = true;
    while(true) {
        iter_cnt += 1;
        if(iter_cnt > 10000)
            return {path};
        if(forward) {
            forward = false;
            if(len >= tip_len + max_diff) {
                res.emplace_back(alternative);
                logger << "Found plausible alternative path " << alternative.len() << std::endl;
                if(res.size() > 10) {
                    logger << "Too many plausible alternatives. Aborting." << std::endl;
                    return {path};
                }
            } else {
                for (Edge &edge : alternative.finish()) {
                    if (edge.getCoverage() >= min_cov || edge.is_reliable) {
                        len += edge.size();
                        alternative.push_back(Segment<Edge>(edge, 0, edge.size()));
                        forward = true;
                        break;
                    }
                }
            }
        } else {
            if (alternative.size() == 0)
                break;
            Edge &old_edge = alternative.back().contig();
            alternative.pop_back();
            len -= old_edge.size();
            bool found = false;
            for(Edge &edge : alternative.finish()) {
                if(edge.getCoverage() >= min_cov || edge.is_reliable) {
                    if(found) {
                        len += edge.size();
                        alternative.push_back(Segment<Edge>(edge, 0, edge.size()));
                        forward = true;
                        break;
                    } else if (&edge == &old_edge) {
                        found - true;
                    }
                }
            }
        }
    }
    logger << "Found " << res.size() << " plausible alternative paths" << std::endl;
    return std::move(res);
}

std::vector<GraphAlignment> FilterAlternatives(logging::Logger &logger1, const GraphAlignment &initial, const std::vector<GraphAlignment> &als,
                                            size_t max_diff, double threshold) {
    size_t len = initial.len();
    std::vector<GraphAlignment> res;
    size_t k = initial.getVertex(0).seq.size();
    for(const GraphAlignment &al : als) {
        CompactPath cpath(al);
        bool ok = true;
        for(size_t i = 0; i < al.size(); i++) {
            if(al[i].contig().getCoverage() < threshold && !al[i].contig().is_reliable) {
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
                                                                std::max<size_t>(100, bulge.len() * 3 / 100), threshold);
    if(dump) {
        logger << "Filtered alternatives" << std::endl;
        for(const auto& it : read_alternatives_filtered) {
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
        if(winner != size_t(-1)) {
            read_alternatives_filtered = {read_alternatives_filtered[winner]};
            if(dump)
                logger << "Winner found " << winner << std::endl;
        } else if(dump)
            logger << "Winner not found" << std::endl;
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

void InitialCorrectionPipeline(logging::Logger &logger, SparseDBG &sdbg, RecordStorage &reads_storage,
                               RecordStorage &ref_storage, size_t threads,
                               const std::experimental::filesystem::path &out_file,
                               double reliable_coverage,
                               double threshold, double bulge_threshold, bool dump);

GraphAlignment BestAlignmentPrefix(const GraphAlignment &al, const Sequence & seq) {
    Sequence candSeq = al.truncSeq();
    Sequence prefix = candSeq.Subseq(0, bestPrefix(seq, candSeq));
    GraphAlignment res(al.start());
    res.extend(prefix);
    return res;
}

GraphAlignment
processTip(logging::Logger &logger, std::ostream &out, const GraphAlignment &tip,
                         const std::vector<GraphAlignment> & alternatives,
                         const RecordStorage &ref_storage,
                         double threshold, bool dump = false) {
    size_t size = tip.len();
    out << size << " tip " << tip.size() << " " << tip.minCoverage();
    if(dump) {
        logger << "New tip from " << tip.start().hash() << tip.start().isCanonical() << " "
               << tip.start().outDeg() << " " << tip.start().inDeg() << std::endl
               << "To " << tip.finish().hash() << tip.finish().isCanonical() << " "
               << tip.finish().outDeg() << " " << tip.finish().inDeg() << std::endl;
    }
    std::vector<GraphAlignment> read_alternatives_filtered =
            FilterAlternatives(logger, tip, alternatives, size_t(-1) / 2, threshold);
    const VertexRecord &rec = ref_storage.getRecord(tip.start());
    size_t genome_support = 0;
    for(const GraphAlignment & al : alternatives) {
        CompactPath compactPath(al);
        if(rec.countStartsWith(compactPath.seq()) > 0) {
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
        logger << alternatives.size() << " " << read_alternatives_filtered.size() << std::endl;
        logger << "Read alternatives" << std::endl;
        for (const GraphAlignment &candidate : alternatives) {
            logger << CompactPath(candidate) << std::endl;
        }
    }
    out << " " << alternatives.size() << " " << genome_support << " "
        << read_alternatives_filtered.size() << " " << filtered_genome_support;
    std::vector<GraphAlignment> trunc_alignments;
    Sequence old = tip.truncSeq();
    for(const GraphAlignment &al : read_alternatives_filtered) {
        trunc_alignments.emplace_back(BestAlignmentPrefix(al, old));
    }
    if(trunc_alignments.size() > 1) {
        if(dump)
            logger << "Multiple choice for tip " << trunc_alignments.size() << std::endl << tip.truncSeq() << std::endl;
        std::vector<Sequence> candidates;
        for(GraphAlignment &cand : trunc_alignments) {
            if(dump)
                logger << cand.truncSeq() << std::endl;
            Sequence candSeq = cand.truncSeq();
            candidates.push_back(candSeq);
        }
        size_t winner = tournament(old, candidates);
        if(dump)
            logger << "Winner found " << winner << std::endl;
        if(winner != size_t(-1)) {
            trunc_alignments = {trunc_alignments[winner]};
        }
    }
    if(dump)
        logger << "Result " << trunc_alignments.size() << std::endl;
    filtered_genome_support = 0;
    for(const GraphAlignment & al : trunc_alignments) {
        if(rec.countStartsWith(CompactPath(al).seq()) > 0) {
            filtered_genome_support += 1;
        }
    }
    out << " " << trunc_alignments.size() << " " << filtered_genome_support;
    out << " " << (trunc_alignments.size() == 1 ? "+" : "-");
    out << std::endl;
    if(trunc_alignments.size() == 1) {
        return std::move(trunc_alignments[0]);
    } else {
        return tip;
    }
}

Edge * checkBorder(Vertex &v) {
    Edge * res = nullptr;
    size_t out_rel = 0;
    for(Edge &edge : v.rc()) {
        if(edge.is_reliable) {
            if (res == nullptr)
                res = &edge.rc();
            else
                return nullptr;
        }
    }
    for(Edge &edge : v) {
        if(edge.is_reliable)
            out_rel += 1;
    }
    if(out_rel == 0)
        return res;
    else
        return nullptr;
}

bool checkInner(Vertex &v) {
    size_t inc_rel = 0;
    size_t out_rel = 0;
    for(Edge &edge : v.rc()) {
        if(edge.is_reliable)
            inc_rel += 1;
    }
    for(Edge &edge : v) {
        if(edge.is_reliable)
            out_rel += 1;
    }
    return out_rel >= 1 && inc_rel >= 1;
}

void FillReliableWithConnections(logging::Logger &logger, SparseDBG &sdbg, double threshold) {
    logger << "Remarking reliable edges" << std::endl;
    for(auto &vit : sdbg) {
        for(Vertex * vp : {&vit.second, &vit.second.rc()}) {
            Vertex &v = *vp;
            for(Edge &edge : v) {
                edge.is_reliable = edge.getCoverage() >= threshold;
            }
        }
    }
    size_t cnt_paths = 0;
    std::vector<Edge *> new_reliable;
    for(auto &vit : sdbg) {
        for(Vertex * vp : {&vit.second, &vit.second.rc()}) {
            Vertex &v = *vp;
            Edge *last = checkBorder(v);
            if(last == nullptr)
                continue;
            typedef std::pair<double, Edge *> StoredValue;
            std::unordered_map<Vertex *, std::pair<double, Edge *>> res;
            std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
            queue.emplace(0, last);
            size_t cnt = 0;
            while(!queue.empty()) {
                cnt += 1;
                if(cnt > 10000) {
                    logger << "Dijkstra too long" << std::endl;
                    break;
                }
                Edge *next = queue.top().second;
                double dist = queue.top().first;
                queue.pop();
                if(res.find(next->end()) != res.end() || checkInner(*next->end()))
                    continue;
                res.emplace(next->end(), std::make_pair(dist, next));
                if (checkBorder(next->end()->rc()) != nullptr) {
                    GraphAlignment al(next->end()->rc());
                    while(next != last) {
                        al.push_back(Segment<Edge>(next->rc(), 0, next->size()));
                        new_reliable.emplace_back(next);
                        logger << "New edge marked as reliable " << next->size() << "(" << next->getCoverage() << ")" << std::endl;
                        next = res[next->start()].second;
                    }
                    logger << "Path of size " << al.size() << " and length " << al.len() << " marked as reliable." << std::endl;
                    cnt_paths += 1;
                    break;
                } else {
                    for(Edge &edge : *next->end()) {
                        double score = edge.size() / std::max<double>(std::min(edge.getCoverage(), threshold), 1);
                        queue.emplace(dist + score, &edge);
                    }
                }
            }
        }
    }
    logger << "Marked " << new_reliable.size() << " edges in " << cnt_paths << " paths as reliable" << std::endl;
    for(Edge *edge : new_reliable) {
        edge->is_reliable = true;
        edge->rc().is_reliable = true;
    }
}

size_t correctLowCoveredRegions(logging::Logger &logger, SparseDBG &sdbg,RecordStorage &reads_storage,
                                RecordStorage &ref_storage,
                                const std::experimental::filesystem::path &out_file,
                                double threshold, double reliable_threshold, size_t k, size_t threads, bool dump) {
    if(dump)
        threads = 1;
    FillReliableWithConnections(logger, sdbg, reliable_threshold);
    ParallelRecordCollector<std::string> results(threads);
    ParallelCounter simple_bulge_cnt(threads);
    ParallelCounter bulge_cnt(threads);
    logger.info() << "Correcting low covered regions in reads" << std::endl;
    omp_set_num_threads(threads);
    size_t max_size = std::min(reads_storage.max_len * 9 / 10, std::max<size_t>(k * 2, 1000));
#pragma omp parallel for default(none) shared(std::cout, reads_storage, ref_storage, results, threshold, k, max_size, logger, simple_bulge_cnt, bulge_cnt, dump, reliable_threshold)
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
            if (edge.getCoverage() >= reliable_threshold || edge.is_reliable ||
                    (edge.start()->inDeg() > 0 && edge.end()->outDeg() > 0 && edge.getCoverage() > threshold) ) {
//              Tips need to pass reliable threshold to avoid being corrected.
                corrected_path.push_back(path[path_pos]);
                continue;
            }
            size_t step_back = 0;
            size_t step_front = 0;
            size_t size = edge.size();
            while(step_back < corrected_path.size() &&
                    (corrected_path[corrected_path.size() - step_back - 1].contig().getCoverage() < reliable_threshold &&
                    !corrected_path[corrected_path.size() - step_back - 1].contig().is_reliable)) {
                size += corrected_path[corrected_path.size() - step_back - 1].size();
                step_back += 1;
            }
            while(step_front + path_pos + 1 < path.size() &&
                    (path[step_front + path_pos + 1].contig().getCoverage() < reliable_threshold &&
                    !path[step_front + path_pos + 1].contig().is_reliable)) {
                size += path[step_front + path_pos + 1].size();
                step_front += 1;
            }
            Vertex &start = corrected_path.getVertex(corrected_path.size() - step_back);
            Vertex &end = path.getVertex(path_pos + 1 + step_front);
            GraphAlignment badPath =
                    corrected_path.subalignment(corrected_path.size() - step_back, corrected_path.size())
                                            + path.subalignment(path_pos, path_pos + 1 + step_front);
            corrected_path.pop_back(step_back);
            if(dump) {
                logger << "Bad read segment " <<    alignedRead.id << " " << path_pos << " " << step_back << " "
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
                    for (Edge &e : tmpv)
                        logger << "Edge out " << e.size() << " " << e.getCoverage() << std::endl;
                    for (Edge &e : tmpv.rc())
                        logger << "Edge in " << e.size() << " " << e.getCoverage() << std::endl;
                }
                if (path_pos + 1 + step_front < path.size()) {
                    logger << "End stop " << step_front << " " << path.getVertex(path_pos + 1 + step_front).hash()
                           << " " << path.getVertex(path_pos + 1 + step_front).isCanonical()
                           << " " << path[path_pos + step_front].contig().getCoverage() << std::endl;
                    Vertex &tmpv = path.getVertex(path_pos + step_front + 1);
                    logger << tmpv.outDeg() << " " << tmpv.inDeg() << std::endl;
                    for (Edge &e : tmpv)
                        logger << "Edge out " << e.size() << " " << e.getCoverage() << std::endl;
                    for (Edge &e : tmpv.rc())
                        logger << "Edge in " << e.size() << " " << e.getCoverage() << std::endl;
                }
            }
            if(step_back == corrected_path.size() && step_front == path.size() - path_pos - 1) {
                if(dump)
                    logger << "Whole read has low coverage. Skipping." << std::endl;
                for(const Segment<Edge> &seg : badPath) {
                    corrected_path.push_back(seg);
                }
            } else if(step_back == corrected_path.size()) {
                if (dump)
                    logger << "Processing incoming tip" << std::endl;
                GraphAlignment tip = badPath.RC();
                std::vector<GraphAlignment> alternatives;
                if(tip.len() < max_size)
                    alternatives = reads_storage.getRecord(tip.start()).getTipAlternatives(tip.len(), threshold);
                if (alternatives.empty())
                    alternatives = FindPlausibleTipAlternatives(logger, tip, std::max<size_t>(size * 3 / 100, 100), 3);
                GraphAlignment substitution = processTip(logger, ss, tip, alternatives, ref_storage,
                                                         threshold, dump);
                GraphAlignment rcSubstitution = substitution.RC();
                corrected_path = std::move(rcSubstitution);
            } else if(step_front == path.size() - path_pos - 1) {
                if (dump)
                    logger << "Processing outgoing tip" << std::endl;
                GraphAlignment tip = badPath;
                std::vector<GraphAlignment> alternatives;
                if(tip.len() < max_size)
                    alternatives = reads_storage.getRecord(tip.start()).getTipAlternatives(tip.len(), threshold);
                if (alternatives.empty())
                    alternatives = FindPlausibleTipAlternatives(logger, tip, std::max<size_t>(size * 3 / 100, 100), 3);
                GraphAlignment substitution = processTip(logger, ss, tip, alternatives, ref_storage,
                                                         threshold, dump);
                for (const Segment<Edge> &seg : substitution) {
                    corrected_path.push_back(seg);
                }
            } else {
                std::vector<GraphAlignment> read_alternatives;
                if(size < max_size)
                    read_alternatives = reads_storage.getRecord(badPath.start()).getBulgeAlternatives(badPath.finish(), threshold);
                if(read_alternatives.empty())
                    read_alternatives = FindPlausibleBulgeAlternatives(logger, badPath, std::max<size_t>(size * 3 / 100, 100), 3);
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
            reads_storage.reroute(alignedRead, path, corrected_path, "low coverage correction");
        }
        results.emplace_back(ss.str());
    }
    logger << "Corrected " << simple_bulge_cnt.get() << " simple bulges" << std::endl;
    logger << "Total " << bulge_cnt.get() << " bulges" << std::endl;
    std::ofstream out;
    out.open(out_file);
    size_t res = 0;
    for(std::string &s : results) {
        if(!s.empty() && s.find('+') != size_t(-1))
            res += 1;
        out << s;
    }
    out.close();
    logger.info() << "Corrected low covered regions in " << res << " reads" << std::endl;
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
            if(start.outDeg() != 2 || start[0].end() != start[1].end()) {
                continue;
            }
            Edge & alt = edge == start[0] ? start[1] : start[0];

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
            if(edge.getCoverage() + alt.getCoverage() > threshold || edge.getCoverage() > alt.getCoverage()) {
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
            reads_storage.reroute(alignedRead, path0, path, "simple bulge corrected");
        }
        results.emplace_back(ss.str());
    }
    size_t bulges = std::unordered_set<Edge*>(bulge_cnt.begin(), bulge_cnt.end()).size();
    size_t collapsable = std::unordered_set<Edge*>(collapsable_cnt.begin(), bulge_cnt.end()).size();
    size_t genome = std::unordered_set<Edge*>(genome_cnt.begin(), bulge_cnt.end()).size();
    size_t corruption = std::unordered_set<Edge*>(corruption_cnt.begin(), bulge_cnt.end()).size();
//    logger << "Bulge collapsing results " << bulges << " " << collapsable << " " << genome << " " << corruption << std::endl;
    logger.info() << "Collapsed bulges in " << bulges << " reads" << std::endl;
    return bulges;
}

size_t correctAT(logging::Logger &logger, RecordStorage &reads_storage, size_t k, size_t threads) {
    ParallelRecordCollector<std::string> results(threads);
    logger.info() << "Correcting dinucleotide errors in reads" << std::endl;
    ParallelCounter cnt(threads);
    omp_set_num_threads(threads);
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
            size_t max_variation = std::max<size_t>(3, (at_cnt1 + at_cnt2) / 6);
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
//            logger  << "Correcting ATAT " << best_support << " " << initial_support  << " "
//                    << at_cnt1 << " " << at_cnt2 << " " << max_variation << " "
//                    << "ACGT"[seq[seq.size() - 2]] << "ACGT"[seq[seq.size() - 1]] << " "
//                    << extension.size() << " " << best_seq.size() << std::endl;
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
            reads_storage.reroute(alignedRead, prev_path, path, "AT corrected");
//            std::cout << "Rerouted " << alignedRead.id << " " << initial_support << " " << best_support << std::endl;
            if(path.size() > 10000) {
                std::cout << extension << "\n" <<best_seq << std::endl;
            }
            corrected += std::max(prev_path.len(), path.len()) - std::min(prev_path.len(), path.len());
        }
        if(corrected > 0) {
//#pragma omp critical
//            {
//                logger << "ATAT " << alignedRead.id << " " << corrected << std::endl;
//            }
            ++cnt;
        }
    }
    logger.info() << "Corrected " << cnt.get() << " dinucleotide sequences" << std::endl;
    return cnt.get();
}

void initialCorrect(SparseDBG &sdbg, logging::Logger &logger,
                    const std::experimental::filesystem::path &out_file,
                    RecordStorage &reads_storage,
                    RecordStorage &ref_storage,
                    double threshold, double bulge_threshold, double reliable_coverage, bool remove_bad, size_t threads, bool dump) {
    size_t k = sdbg.hasher().getK();
    correctAT(logger, reads_storage, k, threads);
    correctLowCoveredRegions(logger,sdbg, reads_storage, ref_storage, out_file, threshold, reliable_coverage, k, threads, dump);
    collapseBulges(logger, reads_storage, ref_storage, out_file, bulge_threshold, k, threads);
    logger.info() << "Applying changes to the graph" << std::endl;
    RemoveUncovered(logger, threads, sdbg, {&reads_storage, &ref_storage});
    sdbg.checkConsistency(threads, logger);
    logger << "Running second round of error correction" << std::endl;
    correctAT(logger, reads_storage, k, threads);
    correctLowCoveredRegions(logger,sdbg, reads_storage, ref_storage, out_file, threshold, reliable_coverage, k, threads, dump);
    collapseBulges(logger, reads_storage, ref_storage, out_file, bulge_threshold, k, threads);
    reads_storage.invalidateBad(logger, threshold);
    logger.info() << "Applying changes to the graph" << std::endl;
    RemoveUncovered(logger, threads, sdbg, {&reads_storage, &ref_storage});
    logger.info() << "Printing reads to disk" << std::endl;
}