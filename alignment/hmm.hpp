#pragma once

#include <cmath>
#include <vector>
#include <istream>
#include "alignment_piece.hpp"

size_t logp(double val) {
    if (val == 0)
        return 1000000000;
    else
        return -100 * log(val);
}

namespace alignment_event {
    enum Event : size_t {
        match, sub, ins, del, homo_ins, homo_del
    };
}

namespace alignment_state {
    static const std::string codes = "+BHIDSE";
    enum State : size_t {
        dflt, bad, homo_ins, ins, del, bad_start, bad_end
    };
}

//Class stores parameters of HMM: number of states, size of alphabet, transition penalty score matrix, observation penalty matrix
class HMM {
private:
    struct Record {
        size_t prev;
        size_t cost;
        Record(size_t prev, size_t cost) : prev(prev), cost(cost) {}
    };

    std::vector<std::vector<size_t>> transitions;
    std::vector<std::vector<size_t>> observations;
    size_t event_num;
    size_t state_num;

public:
    HMM(size_t _event_num, size_t _state_num, std::vector<std::vector<size_t>> &&_transitions,
            std::vector<std::vector<size_t>> &&_observations):
                event_num(_event_num), state_num(_state_num), transitions(_transitions), observations(_observations) {
    }

    HMM(const HMM &) = default;

//    Given observed text and penalties for start and end states this method returns the sequence of hidden states with the lowest penalty.
    std::vector<size_t> states(const std::vector<size_t> &events, const std::vector<size_t> &start_penalties, const std::vector<size_t> &end_penalties) const;

    static HMM load(std::istream stream);
};