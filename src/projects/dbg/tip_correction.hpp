#pragma once

void FillReliableTips(logging::Logger &logger, dbg::SparseDBG &sdbg, double reliable_threshold) {
    logger << "Remarking reliable edges" << std::endl;
    for(auto &vit : sdbg) {
        for(Vertex * vp : {&vit.second, &vit.second.rc()}) {
            Vertex &v = *vp;
            for(Edge &edge : v) {
                edge.is_reliable = true;
            }
        }
    }
    size_t infty = 1000000000;
    std::unordered_map<Vertex *, size_t> max_tip;
    std::vector<Edge*> queue;
    for(Edge &edge : sdbg.edges()) {
        if(edge.end()->outDeg() == 0 && edge.end()->inDeg() == 1 && edge.size() < 10000 && edge.getCoverage() < reliable_threshold) {
            max_tip[edge.end()] = 0;
            queue.emplace_back(&edge);
        }
    }
    while(!queue.empty()) {
        Edge &new_edge = *queue.back();
        queue.pop_back();
        Vertex &v = *new_edge.start();
        bool good = true;
        size_t val = 0;
        for(Edge &out : v) {
            if(out.size() >= 10000 || out.getCoverage() > reliable_threshold || max_tip.find(out.end()) == max_tip.end()) {
                good = false;
                break;
            } else {
                val = std::max(val, max_tip[out.end()] + out.size());
            }
        }
    }
}
