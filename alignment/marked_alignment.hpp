//
// Created by anton on 08.07.2020.
//

#pragma once

#include "hmm.hpp"
#include "cigar_alignment.hpp"

template<class U, class V>
class MarkedAlignment {
private:
    const PositionalAlignment<U, V> al;
    std::vector<size_t> region_types;
    std::vector<size_t> bound_index;
public:
    explicit MarkedAlignment(PositionalAlignment<U, V> _al, std::vector<size_t> &&states, std::vector<size_t> &&matches): al(std::move(_al)) {
        bound_index.push_back(0);
        size_t cur = 1;
        size_t cur_state = states[0];
        auto prev_state = size_t(-1);
        for(size_t i = 0; i < states.size(); i++) {
            cur_state = std::min(cur_state, states[i]);
            if (i == matches[cur]) {
                if(cur_state != prev_state) {
                    region_types.push_back(cur_state);
                    bound_index.push_back(cur);
                }
                prev_state = cur_state;
                cur_state = states[i];
                cur += 1;
            }
        }
    }

    size_t regionSize(size_t index) const {
        return std::max(al.positions_from[bound_index[index + 1]] - al.positions_from[bound_index[index]],
                al.positions_to[bound_index[index + 1]] - al.positions_to[bound_index[index]]);
    }

    size_t regionSize(size_t from, size_t to) const {
        return std::max(al.positions_from[bound_index[to]] - al.positions_from[bound_index[from]],
                        al.positions_to[bound_index[to]] - al.positions_to[bound_index[from]]);
    }

    std::vector<PositionalAlignment<U, V>> ExtractRegions(alignment_state::State rtype, size_t min_rsize = 0, size_t min_len = 0) const {
        size_t last = 0;
        std::vector<PositionalAlignment<U, V>> res;
        for(size_t i = 0; i < region_types.size(); i++) {
            if (region_types[i] == rtype and regionSize(i) >= min_rsize) {
                if(regionSize(last, i) >= min_len) {
                    res.push_back(al.subalignment(last, min_len));
                }
            }
            last = i + 1;
        }
        return res;
    }

    std::vector<string> fourstringForm() const {
        std::vector<char> stateCodes;
        for(size_t i = 0; i < region_types.size(); i++) {
            for(size_t j = bound_index[i]; j < bound_index[i + 1]; j++) {
                size_t sz = std::max(al.positions_from[j + 1] - al.positions_from[j], al.positions_to[j + 1] - al.positions_to[j]);
                for(size_t t = 0; t < sz; t++) {
                    stateCodes.push_back(alignment_state::codes[region_types[i]]);
                }
            }
        }
        stateCodes.push_back(alignment_state::codes[region_types.back()]);
        std::vector<string> res = al.treestringForm();
        res.insert(res.begin(), std::string(stateCodes.begin(), stateCodes.end()));
        return res;
    }

    friend std::ostream& operator<<(std::ostream&, const MarkedAlignment<U, V> &);
};

template<class U, class V>
class MarkingTranslator : public AlignmentTranslator<U, V, MarkedAlignment<U, V>> {
    void mark(const PositionalAlignment<U, V> &al, std::vector<size_t> & observed, std::vector<size_t> & matches) const {
        for(size_t t = 0; t < al.positions_from.size(); t++) {
            matches.push_back(observed.size());
            observed.push_back(alignment_event::match);
            Sequence from = al.seg_from.contig.seq.Subseq(al.positions_from[t], al.positions_from[t + 1] + 1);
            Sequence to = al.seg_to.contig.seq.Subseq(al.positions_to[t], al.positions_to[t + 1] + 1);
            size_t match = std::min(from.size(), to.size()) - 1;
            size_t diff = std::max(from.size(), to.size()) - match - 1;
            for (size_t i = 1; i < match; i++) {
                observed.push_back(alignment_event::sub);
            }
            if (from.size() > to.size()) {
                size_t homo = 0;
                for (size_t j = 0; j + 1 < from.size(); j++) {
                    homo += (from[j] == from[j + 1]);
                }
                homo = std::min(homo, diff);
                for (size_t i = 1; i < homo; i++) {
                    observed.push_back(alignment_event::homo_ins);
                }
                for (size_t i = 1; i < diff - homo; i++) {
                    observed.push_back(alignment_event::ins);
                }
            } else if (from.size() < to.size()) {
                size_t homo = 0;
                for (size_t j = 0; j + 1 < to.size(); j++) {
                    homo += (to[j] == to[j + 1]);
                }
                homo = std::min(homo, diff);
                for (size_t i = 1; i < homo; i++) {
                    observed.push_back(alignment_event::homo_del);
                }
                for (size_t i = 1; i < diff - homo; i++) {
                    observed.push_back(alignment_event::del);
                }
            }
        }
        matches.push_back(observed.size());
        observed.push_back(alignment_event::match);
    }

    const HMM &hmm;
public:

    explicit MarkingTranslator(const HMM &_hmm) : hmm(_hmm){
    }

    std::vector<size_t> predict(const std::vector<size_t> &events) const {
        std::vector<size_t> start_penalties = {0, 1000000000, 1000000000, 1000000000, 1000000000, 0, 1000000000};
        std::vector<size_t> end_penalties = {0, 1000000000, 1000000000, 1000000000, 1000000000, 1000000000, 0};
        return hmm.states(events, start_penalties, end_penalties);
    }

    MarkedAlignment<U, V> translateOne(PositionalAlignment<U, V> &&al) const {
        std::vector<size_t> observed;
        std::vector<size_t> matches;
        mark(al, observed, matches);
        std::vector<size_t> states(this->predict(observed));
        return MarkedAlignment<U, V>(al, std::move(states), std::move(matches));
    }
};

template<class U, class V>
std::ostream &operator<< (std::ostream & os, const MarkedAlignment<U, V> &marking) {
    os << "(";
    for(size_t i = 0; i < marking.region_types.size(); i++) {
        os << marking.region_types[i] << "[" << marking.bound_index[i] << ", " << marking.bound_index[i + 1] << "]";
        if (i + 1 < marking.region_types.size()) {
            os << ", ";
        }
    }
    return os << ")";
}