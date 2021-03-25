#pragma once
#include "sparse_dbg.hpp"
#include "sequences/contigs.hpp"
#include <unordered_map>
#include <utility>


template<class U, class V>
class PerfectAlignment {
public:
    Segment<U> seg_from;
    Segment<V> seg_to;
    PerfectAlignment(const Segment<U> & seg_from_, const Segment<V> &seg_to_) : seg_from(seg_from_), seg_to(seg_to_) {
        VERIFY(seg_from_.size() == seg_to_.size());
    }

    size_t size() {
        return seg_from.size();
    }
};

class CompositeSequence {
private:
    std::vector<Sequence> sequences_;
    size_t left_;
    size_t right_;
    size_t size_;
public:
    CompositeSequence(std::vector<Sequence> sequences, size_t from, size_t to) :
        sequences_(std::move(sequences)), left_(from), right_(0) {
        size_ = 0;
        for(Sequence & sequence : sequences_) {
            size_ += sequence.size();
        }
        VERIFY(left_ + right_ <= size_);
        size_ -= left_ + right_;
        right_ = size_ - to;
    }

    explicit CompositeSequence(std::vector<Sequence> sequences) :
            sequences_(std::move(sequences)), left_(0), right_(0) {
        size_ = 0;
        for(Sequence & sequence : sequences_) {
            size_ += sequence.size();
        }
        VERIFY(left_ + right_ <= size_);
        size_ -= left_ + right_;
    }

    size_t size() const {
        return size_;
    }

    unsigned char operator[](size_t index) const {
        VERIFY(index < size_);
        index += left_;
        size_t cur = 0;
        while(cur < sequences_.size() && index >= sequences_[cur].size()) {
            index -= sequences_[cur].size();
            cur += 1;
        }
        return sequences_[cur][index];
    }

    CompositeSequence operator!() const {
        std::function<Sequence(const Sequence &)> rc = [this](const Sequence &seq) {
            return !seq;
        };
        return {oneline::map(sequences_.rbegin(), sequences_.rend(), rc), right_, left_};
    }

//    TODO reduce sequence vector size
    CompositeSequence Subseq(size_t from, size_t to) const {
        size_t left = from + left_;
        size_t right = size_ - to + right_;
        return {sequences_, left, right};
    }
};

std::vector<PerfectAlignment<Contig, Edge>> align(SparseDBG &dbg, Contig &contig) {
    Sequence seq = contig.seq;
    std::vector<PerfectAlignment<Contig, Edge>> res;
    KWH kwh(dbg.hasher(), seq, 0);
    size_t k = dbg.hasher().k;
    while(true) {
        if (res.empty() || kwh.pos >= res.back().seg_from.right) {
            if(dbg.containsVertex(kwh.hash())) {
                Vertex &vertex = dbg.getVertex(kwh);
                Vertex &rcVertex = vertex.rc();
                if((res.empty() || kwh.pos > res.back().seg_from.right)
                   && kwh.pos > 0 && rcVertex.hasOutgoing(seq[kwh.pos - 1] ^ 3)) {
                    Edge &edge = rcVertex.getOutgoing(seq[kwh.pos - 1] ^ 3);
                    size_t len = 1;
                    while(len < edge.size() && len < kwh.pos && edge.seq[len] == (seq[kwh.pos - len - 1] ^ 3))
                        len += 1;
                    res.emplace_back(Segment<Contig>(contig, kwh.pos - len, kwh.pos),
                                     Segment<Edge>(edge.rc(), edge.size() - len, edge.size()));
                }
                if(kwh.pos + k < seq.size() && vertex.hasOutgoing(seq[kwh.pos + k])) {
                    Edge &edge = vertex.getOutgoing(seq[kwh.pos + k]);
                    size_t len = 1;
                    while(len < edge.size() && kwh.pos + k + len < seq.size() && edge.seq[len] == seq[kwh.pos + k + len])
                        len += 1;
                    res.emplace_back(Segment<Contig>(contig, kwh.pos, kwh.pos + len), Segment<Edge>(edge, 0, len));
                }
            } else if((res.empty() || kwh.pos > res.back().seg_from.right) && dbg.isAnchor(kwh.hash())) {
                typename  SparseDBG::EdgePosition pos = dbg.getAnchor(kwh);
//                TODO replace this code with a call to expand method of PerfectAlignment class after each edge is marked by its full sequence
                Edge &edge = *pos.edge;
                Vertex &start = *pos.start;
                CompositeSequence edge_seq({start.seq, edge.seq});
                size_t left_from = kwh.pos;
                size_t right_from = kwh.pos + k;
                size_t left_to = pos.pos;
                size_t right_to = pos.pos + k;
                while(left_from > 0 && left_to > 0 && edge_seq[left_to - 1] == seq[left_from - 1]) {
                    left_from -= 1;
                    left_to -= 1;
                }
                while(right_from < seq.size() && right_to < edge_seq.size() && seq[right_from] == edge_seq[right_to]) {
                    right_from += 1;
                    right_to += 1;
                }
                if(left_to - left_from > k) {
                    res.emplace_back(Segment<Contig>(contig, left_from, right_from - k),
                                     Segment<Edge>(edge, left_to, right_to - k));
                }
            }
        }
        if (!kwh.hasNext())
            break;
        kwh = kwh.next();
    }
    return std::move(res);
}

class GraphAlignmentStorage {
private:
    std::unordered_map<const Edge *, std::vector<PerfectAlignment<Contig, Edge>>> alignments;
    std::vector<Contig*> stored_contigs;
    SparseDBG & dbg;


    void innerFill(const Contig &old_contig) {
        stored_contigs.emplace_back(new Contig(old_contig));
        Contig &contig = *stored_contigs.back();
        std::vector<PerfectAlignment<Contig, Edge>> path = align(dbg, contig);
        for(PerfectAlignment<Contig, Edge> &al : path) {
            alignments[&al.seg_to.contig()].emplace_back(al);
        }
    }

    void printEdgeDot(std::ostream &os, Vertex & start, Edge &edge) {
        Vertex &end = *edge.end();
        os << "\"";
        if (!start.isCanonical())
            os << "-";
        os << start.hash() % 100000 << "\" -> \"";
        if (!end.isCanonical())
            os << "-";
        os << end.hash() % 100000  << "\" [label=\"" << edge.size() << "(" << edge.getCoverage() << ")";
        std::vector<PerfectAlignment<Contig, Edge>> &als = alignments[&edge];
        os << "\\n" << this->operator()(edge);
        os << "\"]\n";
    }

    void printEdge(std::ostream &os, Vertex & start, Edge &edge) {
        Vertex &end = *edge.end();
        if (!start.isCanonical())
            os << "-";
        os << start.hash() % 100000 << " -> ";
        if (!end.isCanonical())
            os << "-";
        os << end.hash() % 100000 << "\n";
        std::vector<PerfectAlignment<Contig, Edge>> &als = alignments[&edge];
        for(size_t i = 0; i < als.size(); i++) {
            PerfectAlignment<Contig, Edge> &al = als[i];
            os << "\n" << al.seg_from << "->" << al.seg_to;
        }
        os << "\n";
    }

public:
    explicit GraphAlignmentStorage(SparseDBG & dbg_) : dbg(dbg_) {
    }

    GraphAlignmentStorage(const GraphAlignmentStorage &) = delete;

    GraphAlignmentStorage(GraphAlignmentStorage &&other) : dbg(other.dbg) {
        std::swap(other.alignments, alignments);
        std::swap(other.stored_contigs, stored_contigs);
    }

    ~GraphAlignmentStorage() {
        for(Contig * contig : stored_contigs) {
            delete contig;
        }
    }

    void fill(const Contig &contig) {
        innerFill(contig);
        innerFill(contig.RC());
    }


    void printDot(std::ostream &os) {
        os << "digraph {\nnodesep = 0.5;\n";
        for(std::pair<const htype, Vertex> & it : dbg) {
            Vertex &start = it.second;
            for(Edge &edge : start.getOutgoing()) {
                Vertex &end = *edge.end();
                printEdgeDot(os, start, edge);
            }
            for(Edge &edge : start.rc().getOutgoing()) {
                Vertex &end = *edge.end();
                printEdgeDot(os, start.rc(), edge);
            }
        }
        os << "}\n";
    }

    std::string operator()(Edge &edge) const {
        if(alignments.find(&edge) == alignments.end())
            return "";
        const std::vector<PerfectAlignment<Contig, Edge>> &als = alignments.find(&edge)->second;
        if(als.empty()) {
            return "";
        }
        std::stringstream ss;
        size_t num = std::min<size_t>(10, als.size());
        ss << als[0].seg_from << "->" << als[0].seg_to;
        for(size_t i = 1; i < num; i++) {
            const PerfectAlignment<Contig, Edge> &al = als[i];
            ss << "\\n" << al.seg_from << "->" << al.seg_to;
        }
        return ss.str();
    }

    void print(std::ostream &os) {
        os << "digraph {\nnodesep = 0.5;\n";
        for(std::pair<const htype, Vertex> & it : dbg) {
            Vertex &start = it.second;
            for(Edge &edge : start.getOutgoing()) {
                Vertex &end = *edge.end();
                printEdge(os, start, edge);
            }
            for(Edge &edge : start.rc().getOutgoing()) {
                Vertex &end = *edge.end();
                printEdge(os, start.rc(), edge);
            }
        }
        os << "}\n";
    }
};

class Component {
private:
    SparseDBG & graph;
    std::unordered_set<htype, alt_hasher<htype>> v;
    struct EdgeRec {
        Vertex * start;
        Vertex * end;
        size_t size;
        size_t cov;
    };

    size_t outDeg(const Vertex &vert, size_t min_cov) const {
        size_t res = 0;
        for(const Edge &edge : vert.getOutgoing()) {
            if(edge.getCoverage() >= min_cov) {
                res += 1;
            }
        }
        return res;
    }

    bool isUnbranching(const Vertex &vert, size_t min_cov) const {
        return v.find(vert.hash()) != v.end() && outDeg(vert, min_cov) == 1 && outDeg(vert.rc(), min_cov) == 1;
    }

    Edge &getOut(Vertex &vert, size_t min_cov) const {
        for(Edge &edge : vert.getOutgoing()) {
            if(edge.getCoverage() >= min_cov) {
                return edge;
            }
        }
        VERIFY(false);
    }

    Path unbranching(Vertex &vert, Edge &edge, size_t minCov) const {
        std::vector<Edge*> res;
        res.push_back(&edge);
        Vertex *cur = edge.end();
        while (cur != &vert && isUnbranching(*cur, minCov)) {
            res.push_back(&getOut(*cur, minCov));
            cur = res.back()->end();
        }
        return Path(vert, res);
    }

public:
    template<class I>
    Component(SparseDBG &_graph, I begin, I end) : graph(_graph), v(begin, end) {
    }

    Component(SparseDBG &_graph) : graph(_graph) {
        for(auto & vert : graph)
            v.emplace(vert.second.hash());
    }

    template<class I>
    static Component neighbourhood(SparseDBG &graph, I begin, I end, size_t radius, size_t min_coverage = 0) {
        std::unordered_set<htype, alt_hasher<htype>> v;
        typedef std::pair<size_t, htype> StoredValue;
        std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<StoredValue>> queue;
        while(begin != end) {
            queue.emplace(0, *begin);
            ++begin;
        }
        while(!queue.empty()) {
            StoredValue val = queue.top();
            queue.pop();
            if(v.find(val.second) != v.end())
                continue;
            v.insert(val.second);
            if(val.first > radius)
                continue;
            Vertex &vert = graph.getVertex(val.second);
            for(Edge & edge : vert.getOutgoing()) {
                if(edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
            for(Edge & edge : vert.rc().getOutgoing()) {
                if(edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
        }
        return Component(graph, v.begin(), v.end());
    }

    static Path findPath(Vertex &from, Vertex &to, size_t max_len = size_t(-1)) {
        typedef std::tuple<size_t, Edge *, Vertex *> StoredValue;
        std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<StoredValue>> queue;
        queue.emplace(0, nullptr, &from);
        std::unordered_map<Vertex *, Edge *> prev;
        while(!queue.empty()) {
            size_t dist = std::get<0>(queue.top());
            Edge *last_edge = std::get<1>(queue.top());
            Vertex &v2 = *std::get<2>(queue.top());
            queue.pop();
            if(prev.find(&v2) != prev.end())
                continue;
            prev[&v2] = last_edge;
            for(auto & edge : v2.getOutgoing()) {
                if(dist + edge.size() <= max_len)
                    queue.emplace(dist + edge.size(), &edge.rc(), edge.end());
            }
        }
        if(prev.find(&to) == prev.end()) {
            Path res(to.rc());
            Vertex *cur = &to;
            while(cur != &from) {
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
        std::vector<PerfectAlignment<Contig, Edge>> als1 = align(graph, contig);
        Contig rc_contig = contig.RC();
        std::vector<PerfectAlignment<Contig, Edge>> als2 = align(graph, rc_contig);
//        TODO Every edge must have a full sequence stored as a composite sequence
        for(PerfectAlignment<Contig, Edge> &al : als1) {
            queue.emplace(0, al.seg_to.contig().end()->hash());
        }
        for(PerfectAlignment<Contig, Edge> &al : als2) {
            queue.emplace(0, al.seg_to.contig().end()->hash());
        }
        while(!queue.empty()) {
            std::pair<size_t, htype> val = queue.top();
            queue.pop();
            if(v.find(val.second) != v.end())
                continue;
            v.insert(val.second);
            if(val.first > radius)
                continue;
            Vertex &vert = graph.getVertex(val.second);
            for(Edge & edge : vert.getOutgoing()) {
                if(edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
            for(Edge & edge : vert.rc().getOutgoing()) {
                if(edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
        }
        return Component(graph, v.begin(), v.end());
    }

    std::vector<Component> split(size_t len = 200000) const {
        std::vector<Component> res;
        std::unordered_set<htype, alt_hasher<htype>> visited;
        std::vector<htype> component;
        for(const htype &vid : v) {
            std::vector<htype> queue;
            if(visited.find(vid) != visited.end())
                continue;
            queue.push_back(vid);
            while (!queue.empty()) {
                const htype &val = queue.back();
                queue.pop_back();
                if (visited.find(val) != visited.end())
                    continue;
                visited.insert(val);
                component.emplace_back(val);
                Vertex &vert = graph.getVertex(val);
                for (Edge &edge : vert.getOutgoing()) {
                    if (edge.size() < len) {
                        queue.emplace_back(edge.end()->hash());
                        component.emplace_back(edge.end()->hash());
                    }
                }
                for (Edge &edge : vert.rc().getOutgoing()) {
                    if (edge.size() < len) {
                        queue.emplace_back(edge.end()->hash());
                        component.emplace_back(edge.end()->hash());
                    }
                }
            }
            res.emplace_back(graph, component.begin(), component.end());
        }
        return std::move(res);
    }

    std::string vertexLabel(const Vertex &vert) const {
        std::stringstream res;
        if (!vert.isCanonical())
            res << "-";
        res << vert.hash() % 1000000;
        return res.str();
    }

    void printEdge(std::ostream &os, Vertex & start, Edge &edge, const std::string &extra_label = "",
                   const std::string &color = "black") const {
        Vertex &end = *edge.end();
        os << "\"" << vertexLabel(start) << "\" -> \"" << vertexLabel(end) <<
                "\" [label=\"" << edge.size() << "(" << edge.getCoverage() << ")";
        if(!extra_label.empty()) {
            os << "\\n"<<extra_label;
        }
        os << "\", color=\"" + color + "\"]\n";
    }

    void printEdge(std::ostream &os, Path & path) const {
        size_t len = 0;
        size_t cov = 0;
        for(size_t i = 0; i < path.size(); i++) {
            len += path[i].size();
            cov += path[i].intCov();
        }
        Vertex &start = path.start();
        Vertex &end = *path.back().end();
        os << "\"" << vertexLabel(start) << "\" -> \"" << vertexLabel(end)
           << "\" [label=\"" << len << "(" << double(cov) / len << ")\"]\n";
    }

    size_t size() const {
        return v.size();
    }

    void printDot(std::ostream &os, const std::function<std::string(Edge &)> &labeler, size_t min_cov = 0) const {
        os << "digraph {\nnodesep = 0.5;\n";
        std::unordered_set<htype, alt_hasher<htype>> extended;
        for(htype vid : v)
            for(Vertex * vit : graph.getVertices(vid)) {
                Vertex &start = *vit;
                for (Edge &edge : start.getOutgoing()) {
                    extended.emplace(edge.end()->hash());
                }
            }
        for(htype vid : extended) {
            std::string color = v.find(vid) == v.end() ? "yellow" : "white";
            for(Vertex * vit : graph.getVertices(vid)) {
                Vertex &vert = *vit;
                os << vertexLabel(vert) << " [fillcolor=\"" + color + "\"]\n";
            }
        }
        for(htype vid : v) {
            for(Vertex * vit : graph.getVertices(vid)) {
                Vertex &start = *vit;
                for (Edge &edge : start.getOutgoing()) {
                    if (edge.getCoverage() < min_cov)
                        continue;
                    Vertex &end = *edge.end();
                    if (v.find(end.hash()) == v.end()) {
                        printEdge(os, start, edge, labeler(edge), "red");
                        printEdge(os, end.rc(), edge.rc(), labeler(edge.rc()), "red");
                    } else {
                        printEdge(os, start, edge, labeler(edge));
                    }
                }
            }
        }
        os << "}\n";
    }

    void printDot(std::ostream &os, size_t min_cov = 0) const {
        const std::function<std::string(Edge &)> labeler = [](Edge &) {return "";};
        printDot(os, labeler, min_cov);
    }

    void printCompressedDot(std::ostream &os, size_t min_cov = 0) const {
        os << "digraph {\nnodesep = 0.5;\n";
        for(htype vid : v) {
            Vertex &start = graph.getVertex(vid);
            if(isUnbranching(start, min_cov))
                continue;
            for(Edge &edge : start.getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Path path = unbranching(start, edge, min_cov);
                printEdge(os, path);
                Vertex &end = *path.back().end();
                if(v.find(end.hash()) == v.end()) {
                    Path rcpath = path.RC();
                    printEdge(os, rcpath);
                }
            }
            for(Edge &edge : start.rc().getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Path path = unbranching(start.rc(), edge, min_cov);
                printEdge(os, path);
                Vertex &end = *path.back().end();
                if(v.find(end.hash()) == v.end()) {
                    Path rcpath = path.RC();
                    printEdge(os, rcpath);
                }
            }
        }
        os << "}\n";
    }
};

void DrawSplit(const Component &component, const std::experimental::filesystem::path &dir,
               const std::function<std::string(Edge &)> &labeler, size_t len = 200000) {
    ensure_dir_existance(dir);
    std::vector<Component> split = component.split(len);
    for(size_t i = 0; i < split.size(); i++) {
        std::experimental::filesystem::path f = dir / (std::to_string(i) + ".dot");
        std::ofstream os;
        os.open(f);
        component.printDot(os, labeler);
        os.close();
    }
}

void DrawSplit(const Component &component, const std::experimental::filesystem::path &dir, size_t len = 200000) {
    std::function<std::string(Edge &)> labeler = [](Edge &){return "black";};
    DrawSplit(component, dir, labeler, len);
}