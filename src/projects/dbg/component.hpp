#pragma once

#include <utility>


#include "sparse_dbg.hpp"

namespace dbg {
    class Component {
    public:
        SparseDBG &graph;
        std::unordered_set<htype, alt_hasher<htype>> v;
        struct EdgeRec {
            Vertex *start;
            Vertex *end;
            size_t size;
            size_t cov;
        };

        size_t outDeg(const Vertex &vert, size_t min_cov) const {
            size_t res = 0;
            for (const Edge &edge : vert) {
                if (edge.getCoverage() >= min_cov) {
                    res += 1;
                }
            }
            return res;
        }

        bool isUnbranching(const Vertex &vert, size_t min_cov) const {
            return v.find(vert.hash()) != v.end() && outDeg(vert, min_cov) == 1 && outDeg(vert.rc(), min_cov) == 1;
        }

        Edge &getOut(Vertex &vert, size_t min_cov) const {
            for (Edge &edge : vert) {
                if (edge.getCoverage() >= min_cov) {
                    return edge;
                }
            }
            VERIFY(false);
        }

        bool contains(const Vertex &vert) const {
            return v.find(vert.hash()) != v.end();
        }

        Path unbranching(Vertex &vert, Edge &edge, size_t minCov) const {
            std::vector<Edge *> res;
            res.push_back(&edge);
            Vertex *cur = edge.end();
            while (cur != &vert && isUnbranching(*cur, minCov)) {
                res.push_back(&getOut(*cur, minCov));
                cur = res.back()->end();
            }
            return Path(vert, res);
        }

        template<class I>
        Component(SparseDBG &_graph, I begin, I end) : graph(_graph), v(begin, end) {
        }

        Component(SparseDBG &_graph) : graph(_graph) {
            for (auto &vert : graph)
                v.emplace(vert.second.hash());
        }

        template<class I>
        static Component neighbourhood(SparseDBG &graph, I begin, I end, size_t radius, size_t min_coverage = 0) {
            std::unordered_set<htype, alt_hasher<htype>> v;
            typedef std::pair<size_t, htype> StoredValue;
            std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
            while (begin != end) {
                queue.emplace(0, *begin);
                ++begin;
            }
            while (!queue.empty()) {
                StoredValue val = queue.top();
                queue.pop();
                if (v.find(val.second) != v.end())
                    continue;
                v.insert(val.second);
                if (val.first > radius)
                    continue;
                Vertex &vert = graph.getVertex(val.second);
                for (Edge &edge : vert) {
                    if (edge.getCoverage() >= min_coverage)
                        queue.emplace(val.first + edge.size(), edge.end()->hash());
                }
                for (Edge &edge : vert.rc()) {
                    if (edge.getCoverage() >= min_coverage)
                        queue.emplace(val.first + edge.size(), edge.end()->hash());
                }
            }
            return Component(graph, v.begin(), v.end());
        }

        static Path findPath(Vertex &from, Vertex &to, size_t max_len = size_t(-1)) {
            typedef std::tuple<size_t, Edge *, Vertex *> StoredValue;
            std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
            queue.emplace(0, nullptr, &from);
            std::unordered_map<Vertex *, Edge *> prev;
            while (!queue.empty()) {
                size_t dist = std::get<0>(queue.top());
                Edge *last_edge = std::get<1>(queue.top());
                Vertex &v2 = *std::get<2>(queue.top());
                queue.pop();
                if (prev.find(&v2) != prev.end())
                    continue;
                prev[&v2] = last_edge;
                for (auto &edge : v2) {
                    if (dist + edge.size() <= max_len)
                        queue.emplace(dist + edge.size(), &edge.rc(), edge.end());
                }
            }
            if (prev.find(&to) == prev.end()) {
                Path res(to.rc());
                Vertex *cur = &to;
                while (cur != &from) {
                    Edge &edge = *prev[cur];
                    res += edge;
                    cur = &edge.end()->rc();
                }
                return res.RC();
            } else
                return Path(from);
        }

        static Component neighbourhood(SparseDBG &graph, Contig &contig, size_t radius, size_t min_coverage = 0) {
            std::unordered_set<htype, alt_hasher<htype>> v;
            typedef std::pair<size_t, htype> StoredValue;
            std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<StoredValue>> queue;
            std::vector<PerfectAlignment<Contig, Edge>> als1 = graph.carefulAlign(contig);
            Contig rc_contig = contig.RC();
            std::vector<PerfectAlignment<Contig, Edge>> als2 = graph.carefulAlign(rc_contig);
//        TODO Every edge must have a full sequence stored as a composite sequence
            for (PerfectAlignment<Contig, Edge> &al : als1) {
                if(al.seg_to.left < radius)
                    queue.emplace(0, al.seg_to.contig().start()->hash());
                VERIFY(al.seg_to.right <= al.seg_to.contig().size());
                if(al.seg_to.contig().size() - al.seg_to.right < radius)
                    queue.emplace(0, al.seg_to.contig().end()->hash());
            }
            if(queue.empty()) {
                for (PerfectAlignment<Contig, Edge> &al : als1) {
                    queue.emplace(0, al.seg_to.contig().start()->hash());
                    queue.emplace(0, al.seg_to.contig().end()->hash());
                }
            }
            while (!queue.empty()) {
                std::pair<size_t, htype> val = queue.top();
                queue.pop();
                if (v.find(val.second) != v.end() || val.first > radius)
                    continue;
                v.insert(val.second);
                for(Vertex *vit : graph.getVertices(val.second))
                    for (Edge &edge : *vit) {
                        if (edge.getCoverage() >= min_coverage)
                            queue.emplace(val.first + edge.size(), edge.end()->hash());
                    }
            }
            return Component(graph, v.begin(), v.end());
        }

        size_t size() const {
            return v.size();
        }
    };

    class AbstractSplitter {
    public:
        virtual std::vector<Component> split(Component component) const = 0;
        std::vector<Component> split(SparseDBG &dbg) const {
            return split(Component(dbg));
        }
    };

    class ConditionSplitter : public AbstractSplitter {
    private:
        std::function<bool(const Edge &)> splitEdge;

    public:
        explicit ConditionSplitter(std::function<bool(const Edge &)> splitEdge) : splitEdge(std::move(splitEdge)) {
        }

        std::vector<Component> split(Component comp) const override {
            SparseDBG &dbg = comp.graph;
            std::vector<Component> res;
            std::unordered_set<htype, alt_hasher<htype>> visited;
            for (const htype &vid : comp.v) {
                std::vector<htype> queue;
                if (visited.find(vid) != visited.end())
                    continue;
                queue.push_back(vid);
                std::vector<htype> component;
                while (!queue.empty()) {
                    htype val = queue.back();
                    queue.pop_back();
                    if (visited.find(val) != visited.end())
                        continue;
                    visited.insert(val);
                    component.emplace_back(val);
                    for(Vertex *vert : dbg.getVertices(val)) {
                        for (Edge &edge : *vert) {
                            if (!splitEdge(edge) && comp.contains(*edge.end())) {
                                queue.emplace_back(edge.end()->hash());
                                component.emplace_back(edge.end()->hash());
                            }
                        }
                    }
                }
                res.emplace_back(dbg, component.begin(), component.end());
            }
            return std::move(res);
        }
    };

    class LengthSplitter : public ConditionSplitter {
    public:
        explicit LengthSplitter(size_t min_len) :
                    ConditionSplitter([min_len](const Edge& edge){return edge.size() > min_len;}){
        }
    };
}