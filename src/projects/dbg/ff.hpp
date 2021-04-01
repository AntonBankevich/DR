#pragma once

#include <queue>
#include <unordered_map>
#include <sequences/verify.hpp>
#include <unordered_set>
#include "vector"
class Network {
public:
    struct Vertex;
    struct Edge {
        Edge(int id, size_t start, size_t end, size_t capacity) : id(id), start(start), end(end),
                                                                                      capacity(capacity) {}

        int id;
        size_t start;
        size_t end;
        size_t capacity;
        size_t flow = 0;
    };

    struct Vertex {
        explicit Vertex(size_t _id) : id(_id) {
        }

        size_t id;
        std::vector<int> out;
        std::vector<int> inc;
    };

    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<Edge> back_edges;
    size_t source;
    size_t sink;

    Edge &getEdge(int id) {
        if (id > 0)
            return edges[id - 1];
        else
            return back_edges[-id - 1];
    }

    Network() {
        source = addVertex();
        sink = addVertex();
    }

    size_t addVertex() {
        vertices.emplace_back(vertices.size());
        return vertices.size() - 1;
    }

    int addEdge(size_t from, size_t to, size_t capacity) {
        int eid = edges.size() + 1;
        edges.emplace_back(eid, from, to, capacity);
        back_edges.emplace_back(-eid, to, from, 0);
        vertices[from].out.push_back(eid);
        vertices[to].inc.push_back(eid);
        vertices[to].out.push_back(-eid);
        vertices[from].inc.push_back(-eid);
        return eid;
    }

    int addEdge(size_t from, size_t to, size_t min_capacity, size_t max_capacity) {
        if(min_capacity > 0) {
            addSource(to, min_capacity);
            addSink(from, min_capacity);
            max_capacity -= min_capacity;
        }
        return addEdge(from, to, max_capacity);
    }

    void addSource(size_t id, size_t capacity) {
        addEdge(source, id, capacity);
    }

    void addSink(size_t id, size_t capacity) {
        addEdge(id, sink, capacity);
    }

    std::vector<int> bfs(int startId, int endId) {
        std::queue<size_t> queue;
        std::unordered_map<size_t, int> prev;
        prev[startId] = 0;
        queue.push(startId);
        while(!queue.empty()) {
            size_t next = queue.front();
            queue.pop();
            if (next == endId) {
                std::vector<int> res;
                while(next != startId) {
                    res.emplace_back(prev[next]);
                    next = getEdge(prev[next]).start;
                }
                return {res.rbegin(), res.rend()};
            }
            for(int eid : vertices[next].out) {
                Edge &edge = getEdge(eid);
                size_t end = edge.end;
                if(edge.capacity > 0 && prev.find(end) == prev.end()) {
                    prev[end] = eid;
                    queue.emplace(end);
                }
            }
        }
        return {};
    }

    void pushFlow(int edgeId, size_t val = 1) {
        VERIFY(val <= getEdge(edgeId).capacity);
        getEdge(edgeId).capacity -= val;
        getEdge(-edgeId).capacity += val;
    }

    void pushFlow(const std::vector<int> &path, size_t val = 1) {
        for(int eid : path)
            pushFlow(eid, val);
    }

    size_t outCapasity(size_t vId) {
        size_t outdeg = 0;
        for(int eid : vertices[vId].out) {
            outdeg += getEdge(eid).capacity;
        }
        return outdeg;
    }

    size_t inCapasity(size_t vId) {
        size_t indeg = 0;
        for(int eid : vertices[vId].inc) {
            indeg += getEdge(eid).capacity;
        }
        return indeg;
    }

    bool fillNetwork() {
        size_t outdeg = outCapasity(source);
        for(size_t i = 0; i < outdeg; i++) {
            std::vector<int> path = bfs(source, sink);
            if(path.empty())
                return false;
            std::cout << "New path" << std::endl;
            for(int x : path) {
                std::cout << x << " ";
            }
            std::cout << std::endl;
            pushFlow(path, 1);
        }
        return true;
    }

    bool isInLoop(int edgeId) {
        Edge edge = getEdge(edgeId);
        if(edge.capacity == 0)
            return false;
        Edge &redge = getEdge(-edgeId);
        size_t rcap = redge.capacity;
        redge.capacity = 0;
        bool result = !bfs(redge.start, redge.end).empty();
        redge.capacity = rcap;
        return result;
    }

    std::unordered_map<int, size_t> findFixedMultiplicities() {
        std::unordered_map<int, size_t> res;
        for(Edge &e : edges) {
            if(e.start != source && e.end != sink && !isInLoop(e.id) && !isInLoop(-e.id)) {
                res[e.id] = getEdge(-e.id).capacity;
            }
        }
        return std::move(res);
    }

    std::vector<size_t> topSort() {
        std::unordered_set<size_t> visited;
        std::vector<size_t> stack;
        std::vector<size_t> result;
        for(Vertex &v : vertices)
            stack.push_back(v.id);
        while(!stack.empty()) {
            size_t vid = stack.back();
            stack.pop_back();
            if(vid >= vertices.size()) {
                result.push_back(vid - vertices.size());
                continue;
            }
            if(visited.count(vid) > 0)
                continue;
            visited.emplace(vid);
            stack.push_back(vid + vertices.size());
            for(int eid : vertices[vid].out) {
                Edge & edge = getEdge(eid);
                if(edge.capacity > 0 && visited.count(edge.end) == 0) {
                    stack.push_back(edge.end);
                }
            }
        }
        return {result.rbegin(), result.rend()};
    }

    std::unordered_map<size_t, size_t> strongComponents() {
        std::vector<size_t> order = topSort();
        std::unordered_map<size_t, size_t> res;
        std::vector<std::pair<size_t, size_t>> stack;
        for(size_t vid : order)
            stack.emplace_back(vid, vid);
        while(!stack.empty()) {
            size_t vid = stack.back().first;
            size_t color = stack.back().second;
            stack.pop_back();
            if(res.find(vid) != res.end())
                continue;
            res[vid] = color;
            for(int eid : vertices[vid].inc) {
                Edge & edge = getEdge(eid);
                if(edge.capacity > 0 && res.find(edge.start) == res.end()) {
                    stack.emplace_back(edge.start, color);
                }
            }
        }
        return std::move(res);
    }
};