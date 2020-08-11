//
// Created by anton on 07.07.2020.
//
#include "hmm.hpp"

std::vector<size_t> HMM::states(const std::vector<size_t> &events, const std::vector<size_t> &start_penalties,
                                const std::vector<size_t> &end_penalties) const {
    std::vector<Record> res;
    res.reserve(state_num * events.size());
    for(size_t st = 0; st < state_num; st++) {
        res.emplace_back(st, start_penalties[st] + observations[st][events[0]]);
    }
    size_t rec_pos = 0;
    for(size_t state_pos = 1; state_pos < events.size(); state_pos++) {
        for(size_t j = 0; j < state_num; j++) {
            res.emplace_back(0, res[rec_pos].cost + transitions[0][j]);
            for(size_t i = 1; i < state_num; i++) {
                size_t tmp = res[rec_pos + i].cost + transitions[i][j];
                if (tmp < res.back().cost) {
                    res.back().cost = tmp;
                    res.back().prev = i;
                }
            }
            res.back().cost += observations[j][events[state_pos]];
        }
        rec_pos += state_num;
    }
    std::vector<size_t> states(events.size());
    size_t cost = res.back().cost + end_penalties[state_num - 1];
    states.back() = state_num - 1;
    for(size_t i = 0; i + 1 < state_num; i++) {
        size_t endcost = res[states.size() - state_num + i].cost + end_penalties[i];
        if (endcost < cost) {
            states.back() = i;
            cost = endcost;
        }
    }
    rec_pos = res.size() - state_num;
    for(size_t i = events.size() - 1; i > 0; i--) {
        states[i - 1] = res[rec_pos + states[i]].prev;
        rec_pos -= state_num;
    }
    return states;
}

HMM HMM::load(std::istream &stream) {
    size_t n = 0, m = 0;
    double val = 0;
    stream >> n;
    stream >> m;
    std::vector<std::vector<size_t>> transitions(m);
    std::vector<std::vector<size_t>> observations(m);
    for(size_t i = 0; i < m; i++)
        for (size_t j = 0; j < m; j++) {
            stream >> val;
            transitions[i].push_back(logp(val));
        }
    for(size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++) {
            stream >> val;
            observations[i].push_back(logp(val));
        }
    return HMM(std::move(transitions), std::move(observations));
}
