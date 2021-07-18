#pragma once

#include "sparse_dbg.hpp"
#include <sequences/sequence.hpp>

class Connection {
private:
    dbg::EdgePosition pos1;
    dbg::EdgePosition pos2;
    Sequence connection;
public:
    Connection(dbg::EdgePosition pos1, dbg::EdgePosition pos2, Sequence connection) :
            pos1(pos1), pos2(pos2), connection(connection) {
        VERIFY(connection.startsWith(pos1.kmerSeq()));
        VERIFY(!connection.startsWith(pos2.RC().kmerSeq()));
    }

    Connection shrink() {
        size_t k = pos1.edge->start()->seq.size();
        size_t left = 0;
        size_t right = 0;
        while(left + k < connection.size() && pos1.pos + left < pos1.edge->seq.size() &&
            pos1.edge->seq[pos1.pos + left] == connection[left + k]) {
            left++;
        }
        dbg::EdgePosition rc = pos2.RC();
        Sequence rcSeq = !connection;
        while(right + k < connection.size() && pos2.pos + right < pos2.edge->seq.size() &&
              pos2.edge->seq[pos2.pos + right] == rcSeq[right + k]) {
            right++;
        }
        VERIFY(left + right + k < connection.size());
        return {dbg::EdgePosition(*pos1.edge, pos1.pos + left),
                dbg::EdgePosition(*pos2.edge, pos2.pos - right),
                connection.Subseq(left, connection.size() - right)};
    }
};

class GapCloser {
private:
    size_t min_overlap;
    size_t max_overlap;
    size_t smallK;
    double allowed_divergence;

    struct OverlapRecord {
        OverlapRecord(size_t from, size_t to, size_t matchSizeFrom, size_t matchSizeTo) : from(from), to(to),
                                                                                          match_size_from(
                                                                                                  matchSizeFrom),
                                                                                          match_size_to(matchSizeTo) {}

        size_t from;
        size_t to;
        size_t match_size_from;
        size_t match_size_to;
    };
public:
    GapCloser(size_t min_overlap, size_t max_overlap, size_t smallK, double allowed_divergence) :
            min_overlap(min_overlap), max_overlap(max_overlap), smallK(smallK), allowed_divergence(allowed_divergence){}


    std::pair<size_t, size_t> CheckOverlap(const Sequence &s1, const Sequence &s2) {
        Sequence a = s1.Subseq(s1.size() - std::min(s1.size(), max_overlap));
        Sequence b = s2.Subseq(s2.size() - std::min(s2.size(), max_overlap));
        int64_t mult = a.size() + 1;
        int64_t match = 2 * mult;
        int64_t mismatch = 1 * mult;
        int64_t indel = 1 * mult;
        std::vector<int64_t> res(a.size() + 1);
        for(size_t i = 0; i <= a.size(); i++) {
            res[i] = i;
        }
        std::vector<int64_t> prev(a.size() + 1);
        size_t best = 0;
        size_t best_val = res[a.size()];
        for(size_t j = 1; j <= b.size(); j++) {
            prev = std::move(res);
            res[0] = prev[0] - indel; // NOLINT(bugprone-use-after-move)
            for(size_t i = 1; i < a.size(); i++) {
                if(a[i - 1] == b[j - 1]) {
                    res[i] = prev[j - 1] + match;
                } else {
                    res[i] = std::max(res[i - 1] - indel, std::max(prev[i] - indel, prev[i - 1] - mismatch));
                }
            }
            if(best_val < res[a.size()]) {
                best = j;
                best_val = res[a.size()];
            }
        }
        size_t l1 = a.size() - (best_val % mult);
        size_t l2 = best;
        best_val = best_val / mult * mult;
        double min_val = match * (1 - allowed_divergence) - allowed_divergence * std::max(indel, mismatch);
        if(l1 < min_overlap || l2 < min_overlap || best_val < std::max(l1, l2) * min_val)
            return {0, 0};
        return {l1, l2};
    }

    std::vector<Connection> GapPatches(dbg::SparseDBG &dbg, logging::Logger &logger, size_t threads) {
        logger.info() << "Started gap closing procedure" << std::endl;
        size_t k = dbg.hasher().getK();
        std::vector<dbg::Edge *> tips;
        for (dbg::Edge &edge : dbg.edges()) {
            if (edge.size() > min_overlap && edge.end()->outDeg() == 0 && edge.end()->inDeg() == 1)
                tips.emplace_back(&edge);
        }
        ParallelRecordCollector<std::pair<hashing::htype, size_t>> candidates(threads);
        hashing::RollingHash smallHasher(smallK, 239);
        omp_set_num_threads(threads);
        logger.info() << "Collecting k-mers from tips" << std::endl;
#pragma omp parallel for default(none) shared(tips, candidates, dbg, smallHasher)
        for (size_t i = 0; i < tips.size(); i++) {
            size_t max_len = std::min(tips[i]->size(), max_overlap);
            hashing::KWH kwh(smallHasher, tips[i]->seq, tips[i]->size() - max_len);
            while (true) {
                candidates.emplace_back(kwh.hash(), i);
                if (!kwh.hasNext())
                    break;
                kwh = kwh.next();
            }
        }
        logger.info() << "Sorting k-mers from tips" << std::endl;
        std::vector<std::pair<hashing::htype, size_t>> candidates_list = candidates.collect();
        __gnu_parallel::sort(candidates_list.begin(), candidates_list.end());
        std::vector<std::pair<size_t, size_t>> pairs;
        std::vector<size_t> tmp;
        for (size_t i = 0; i < candidates_list.size(); i++) {
            tmp.emplace_back(candidates_list[i].second);
            if (i + 1 == candidates_list.size() || candidates_list[i + 1].first != candidates_list[i].first) {
                for (size_t t1 : tmp)
                    for (size_t t2 : tmp)
                        if (t1 < t2)
                            pairs.emplace_back(t1, t2);
                tmp = {};
            }
        }
        __gnu_parallel::sort(pairs.begin(), pairs.end());
        pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
        logger.info() << "Found " << pairs.size() << " potential overlaps. Aligning." << std::endl;
        ParallelRecordCollector<OverlapRecord> filtered_pairs(threads);
#pragma omp parallel for default(none) shared(pairs, tips, filtered_pairs)
        for (size_t i = 0; i < pairs.size(); i++) {
            Sequence s1 = tips[pairs[i].first]->seq;
            Sequence s2 = tips[pairs[i].second]->seq;
            std::pair<size_t, size_t> overlap = CheckOverlap(s1, !s2);
            if(overlap.first > 0) {
                filtered_pairs.emplace_back(pairs[i].first, pairs[i].second, overlap.first, overlap.second);
            }
        }
        std::vector<size_t> degs(tips.size(), 0);
        for(OverlapRecord &rec : filtered_pairs) {
            degs[rec.from]++;
            degs[rec.to]++;
        }
        std::vector<Connection> res;
        for(OverlapRecord &rec : filtered_pairs) {
            if(degs[rec.from] == 1 && degs[rec.to] == 1) {
                Sequence seq1 = tips[rec.from]->suffix(tips[rec.from]->size() - rec.match_size_from);
                Sequence seq2 = tips[rec.to]->kmerSeq(tips[rec.to]->size() - rec.match_size_to);
                Connection gap(dbg::EdgePosition(*tips[rec.from], tips[rec.from]->size() - rec.match_size_from - k),
                               dbg::EdgePosition(tips[rec.to]->rc(), rec.match_size_to + k),
                               seq1 +!seq2);
                res.emplace_back(gap.shrink());
            }
        }
        return std::move(res);
    }
};