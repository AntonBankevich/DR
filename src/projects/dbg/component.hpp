#pragma once

#include "paths.hpp"
#include "sparse_dbg.hpp"
#include "common/iterator_utils.hpp"
#include <utility>

namespace dbg {
    class Component {
    public:
        SparseDBG &graph;
        std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> v;
        typedef std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>>::const_iterator iterator;
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

        size_t borderEdges() const {
            size_t res = 0;
            for (hashing::htype hash : v) {
                const Vertex &vert = graph.getVertex(hash);
                for(Edge &edge : vert) {
                    if(!contains(*edge.end()))
                        res++;
                }
                for(Edge &edge : vert.rc()) {
                    if(!contains(*edge.end()))
                        res++;
                }
            }
            return res;
        }

        size_t isAcyclic() const {
            std::unordered_set<Vertex *> visited;
            for(hashing::htype hash : v) {
                Vertex & compStart = graph.getVertex(hash);
                for(Vertex *vit : {&compStart, &compStart.rc()}) {
                    if(visited.find(vit) != visited.end())
                        continue;
                    std::vector<Vertex *> queue;
                    queue.emplace_back(vit);
                    while(!queue.empty()) {
                        Vertex &vert = *queue.back();
                        queue.pop_back();
                        if(visited.find(&vert) != visited.end())
                            continue;
                        bool ok = true;
                        for(Edge &edge : vert.rc()) {
                            Vertex &prev = edge.end()->rc();
                            if(contains(prev) && visited.find(&prev) == visited.end())
                                ok = false;
                        }
                        if(!ok)
                            continue;
                        visited.emplace(&vert);
                        for(Edge &edge : vert) {
                            if(contains(*edge.end()))
                                queue.emplace_back(edge.end());
                        }
                    }
                }
            }
            return visited.size() == v.size() * 2;
        }

        size_t realCC() const {
            std::unordered_set<Vertex *> visited;
            size_t cnt = 0;
            for(hashing::htype hash : v) {
                Vertex & compStart = graph.getVertex(hash);
                for(Vertex *vit : {&compStart, &compStart.rc()}) {
                    if(visited.find(vit) != visited.end())
                        continue;
                    cnt += 1;
                    std::vector<Vertex *> queue;
                    queue.emplace_back(vit);
                    while(!queue.empty()) {
                        Vertex &vert = *queue.back();
                        queue.pop_back();
                        if(visited.find(&vert) != visited.end())
                            continue;
                        visited.emplace(&vert);
                        for(Edge &edge : vert) {
                            if(contains(*edge.end()))
                                queue.emplace_back(edge.end());
                        }
                        for(Edge &edge : vert.rc()) {
                            if(contains(*edge.end()))
                                queue.emplace_back(&edge.end()->rc());
                        }
                    }
                }
            }
            return cnt;
        }

        Edge &getOut(Vertex &vert, size_t min_cov) const {
            for (Edge &edge : vert) {
                if (edge.getCoverage() >= min_cov) {
                    return edge;
                }
            }
            VERIFY(false);
        }

        IterableStorage<ApplyingIterator<iterator, Vertex, 2>> vertices(bool unique = false) const {
            std::function<std::array<Vertex*, 2>(const hashing::htype &)> apply = [this, unique](const hashing::htype &hash) -> std::array<Vertex*, 2> {
                size_t cur = 0;
                Vertex &vertex = graph.getVertex(hash);
                if(unique)
                    return {&graph.getVertex(hash)};
                else
                    return graph.getVertices(hash);
            };
            ApplyingIterator<iterator, Vertex, 2> begin(v.begin(), v.end(), apply);
            ApplyingIterator<iterator, Vertex, 2> end(v.end(), v.end(), apply);
            return {begin, end};
        }

        IterableStorage<ApplyingIterator<iterator, Vertex, 2>> verticesUnique(bool unique = false) const {
            return vertices(true);
        }

        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edges(bool inner = false, bool unique = false) const {
            std::function<std::array<Edge*, 16>(const hashing::htype &)> apply = [this, inner, unique](const hashing::htype &hash) {
                std::array<Edge*, 16> res = {};
                size_t cur = 0;
                for(Vertex * vertex : graph.getVertices(hash)) {
                    for (Edge &edge : *vertex) {
                        bool isInner = contains(*edge.end());
                        if(inner && !isInner)
                            continue;
                        if(!isInner || edge <= edge.rc()) {
                            res[cur] = &edge;
                            cur++;
                            if(!unique) {
                                res[cur] = &edge.rc();
                                cur++;
                            }
                        }
                    }
                }
                return res;
            };
            ApplyingIterator<iterator, Edge, 16> begin(v.begin(), v.end(), apply);
            ApplyingIterator<iterator, Edge, 16> end(v.end(), v.end(), apply);
            return {begin, end};
        }

        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edgesInner() const {
            return edges(true, false);
        }

        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edgesUnique() const {
            return edges(false, true);
        }

        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edgesInnerUnique() const {
            return edges(true, true);
        }

        bool contains(const Vertex &vert) const {
            return v.find(vert.hash()) != v.end();
        }

        bool covers(const Vertex &vert) const {
            if(contains(vert))
                return true;
            for(Edge &edge :vert)
                if(!contains(*edge.end()))
                    return false;
            for(Edge &edge :vert.rc())
                if(!contains(*edge.end()))
                    return false;
            return true;
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
            std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> v;
            typedef std::pair<size_t, hashing::htype> StoredValue;
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

        static Component neighbourhood(SparseDBG &graph, Contig &contig, size_t radius) {
            std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> v;
            typedef std::pair<size_t, hashing::htype> StoredValue;
            std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<StoredValue>> queue;
            std::vector<PerfectAlignment<Contig, Edge>> als1 = GraphAligner(graph).carefulAlign(contig);
            Contig rc_contig = contig.RC();
//        TODO Every edge must have a full sequence stored as a composite sequence
            for (PerfectAlignment<Contig, Edge> &al : als1) {
                if(al.seg_to.left < radius)
                    queue.emplace(0, al.seg_to.contig().start()->hash());
                VERIFY(al.seg_to.right <= al.seg_to.contig().size());
                if(al.seg_to.contig().size() < radius + al.seg_to.right)
                    queue.emplace(0, al.seg_to.contig().end()->hash());
            }
            if(queue.empty()) {
                for (PerfectAlignment<Contig, Edge> &al : als1) {
                    queue.emplace(0, al.seg_to.contig().start()->hash());
                    queue.emplace(0, al.seg_to.contig().end()->hash());
                }
            }
            while (!queue.empty()) {
                std::pair<size_t, hashing::htype> val = queue.top();
                queue.pop();
                if (v.find(val.second) != v.end() || val.first > radius)
                    continue;
                v.insert(val.second);
                for(Vertex *vit : graph.getVertices(val.second))
                    for (Edge &edge : *vit) {
                        queue.emplace(val.first + edge.size(), edge.end()->hash());
                    }
            }
            return Component(graph, v.begin(), v.end());
        }

        static Component longEdgeNeighbourhood(SparseDBG &graph, Contig &contig, size_t long_edge_threshold) {
            std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> v;
            typedef std::pair<size_t, hashing::htype> StoredValue;
            std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
            std::vector<PerfectAlignment<Contig, Edge>> als1 = GraphAligner(graph).carefulAlign(contig);
            Contig rc_contig = contig.RC();
//        TODO Every edge must have a full sequence stored as a composite sequence
            for (PerfectAlignment<Contig, Edge> &al : als1) {
                queue.emplace(0, al.seg_to.contig().start()->hash());
                queue.emplace(0, al.seg_to.contig().end()->hash());
            }
            while (!queue.empty()) {
                std::pair<size_t, hashing::htype> val = queue.top();
                queue.pop();
                if (v.find(val.second) != v.end())
                    continue;
                v.insert(val.second);
                for(Vertex *vit : graph.getVertices(val.second))
                    for (Edge &edge : *vit) {
                        if (edge.size() < long_edge_threshold)
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
            std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> visited;
            for (const hashing::htype &vid : comp.v) {
                std::vector<hashing::htype> queue;
                if (visited.find(vid) != visited.end())
                    continue;
                queue.push_back(vid);
                std::vector<hashing::htype> component;
                while (!queue.empty()) {
                    hashing::htype val = queue.back();
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