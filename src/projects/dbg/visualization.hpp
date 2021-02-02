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

template<typename htype>
std::vector<PerfectAlignment<Contig, Edge<htype>>> align(SparseDBG<htype> &dbg, Contig &contig) {
    Sequence seq = contig.seq;
    std::vector<PerfectAlignment<Contig, Edge<htype>>> res;
    KWH<htype> kwh(dbg.hasher(), seq, 0);
    size_t k = dbg.hasher().k;
    while(true) {
        if (res.empty() || kwh.pos >= res.back().seg_from.right) {
            if(dbg.containsVertex(kwh.hash())) {
                Vertex<htype> &vertex = dbg.getVertex(kwh);
                Vertex<htype> &rcVertex = vertex.rc();
                if((res.empty() || kwh.pos > res.back().seg_from.right)
                   && kwh.pos > 0 && rcVertex.hasOutgoing(seq[kwh.pos - 1] ^ 3)) {
                    Edge<htype> &edge = rcVertex.getOutgoing(seq[kwh.pos - 1] ^ 3);
                    size_t len = 1;
                    while(len < edge.size() && len < kwh.pos && edge.seq[len] == (seq[kwh.pos - len - 1] ^ 3))
                        len += 1;
                    res.emplace_back(Segment<Contig>(contig, kwh.pos - len, kwh.pos),
                                     Segment<Edge<htype>>(rcVertex.rcEdge(edge), edge.size() - len, edge.size()));
                }
                if(kwh.pos + k < seq.size() && vertex.hasOutgoing(seq[kwh.pos + k])) {
                    Edge<htype> &edge = vertex.getOutgoing(seq[kwh.pos + k]);
                    size_t len = 1;
                    while(len < edge.size() && kwh.pos + k + len < seq.size() && edge.seq[len] == seq[kwh.pos + k + len])
                        len += 1;
                    res.emplace_back(Segment<Contig>(contig, kwh.pos, kwh.pos + len), Segment<Edge<htype>>(edge, 0, len));
                }
            } else if((res.empty() || kwh.pos > res.back().seg_from.right) && dbg.isAnchor(kwh.hash())) {
                typename  SparseDBG<htype>::EdgePosition pos = dbg.getAnchor(kwh);
//                TODO replace this code with a call to expand method of PerfectAlignment class after each edge is marked by its full sequence
                Edge<htype> &edge = *pos.edge;
                Vertex<htype> &start = *pos.start;
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
                                     Segment<Edge<htype>>(edge, left_to, right_to - k));
                }
            }
        }
        if (!kwh.hasNext())
            break;
        kwh = kwh.next();
    }
    return std::move(res);
}

template<typename htype>
class GraphAlignmentStorage {
private:
    std::unordered_map<const Edge<htype> *, std::vector<PerfectAlignment<Contig, Edge<htype>>>> alignments;
    std::vector<Contig*> stored_contigs;
    SparseDBG<htype> & dbg;


    void innerFill(const Contig &old_contig) {
        stored_contigs.emplace_back(new Contig(old_contig));
        Contig &contig = *stored_contigs.back();
        std::vector<PerfectAlignment<Contig, Edge<htype>>> path = align(dbg, contig);
        for(PerfectAlignment<Contig, Edge<htype>> &al : path) {
            alignments[&al.seg_to.contig()].emplace_back(al);
        }
    }

    void printEdgeDot(std::ostream &os, Vertex<htype> & start, Edge<htype> &edge) {
        Vertex<htype> &end = *edge.end();
        os << "\"";
        if (!start.isCanonical())
            os << "-";
        os << start.hash() % 100000 << "\" -> \"";
        if (!end.isCanonical())
            os << "-";
        os << end.hash() % 100000  << "\" [label=\"" << edge.size() << "(" << edge.getCoverage() << ")";
        std::vector<PerfectAlignment<Contig, Edge<htype>>> &als = alignments[&edge];
        size_t num = std::min<size_t>(10, als.size());
        for(size_t i = 0; i < num; i++) {
            PerfectAlignment<Contig, Edge<htype>> &al = als[i];
            os << "\\n" << al.seg_from << "->" << al.seg_to;
        }
        if(num < als.size())
            os << "\\n" << "... (" << als.size() - num << ")";

        os << "\"]\n";
    }

    void printEdge(std::ostream &os, Vertex<htype> & start, Edge<htype> &edge) {
        Vertex<htype> &end = *edge.end();
        if (!start.isCanonical())
            os << "-";
        os << start.hash() % 100000 << " -> ";
        if (!end.isCanonical())
            os << "-";
        os << end.hash() % 100000 << "\n";
        std::vector<PerfectAlignment<Contig, Edge<htype>>> &als = alignments[&edge];
        for(size_t i = 0; i < als.size(); i++) {
            PerfectAlignment<Contig, Edge<htype>> &al = als[i];
            os << "\n" << al.seg_from << "->" << al.seg_to;
        }
        os << "\n";
    }

public:
    explicit GraphAlignmentStorage(SparseDBG<htype> & dbg_) : dbg(dbg_) {
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
        for(std::pair<const htype, Vertex<htype>> & it : dbg) {
            Vertex<htype> &start = it.second;
            for(Edge<htype> &edge : start.getOutgoing()) {
                Vertex<htype> &end = *edge.end();
                printEdgeDot(os, start, edge);
            }
            for(Edge<htype> &edge : start.rc().getOutgoing()) {
                Vertex<htype> &end = *edge.end();
                printEdgeDot(os, start.rc(), edge);
            }
        }
        os << "}\n";
    }

    std::string operator()(Edge<htype> &edge) const {
        if(alignments.find(&edge) == alignments.end())
            return "";
        const std::vector<PerfectAlignment<Contig, Edge<htype>>> &als = alignments.find(&edge)->second;
        if(als.empty()) {
            return "";
        }
        std::stringstream ss;
        ss << als[0].seg_from << "->" << als[1].seg_to;
        for(size_t i = 1; i < als.size(); i++) {
            const PerfectAlignment<Contig, Edge<htype>> &al = als[i];
            ss << "\n" << al.seg_from << "->" << al.seg_to;
        }
        return ss.str();
    }

    void print(std::ostream &os) {
        os << "digraph {\nnodesep = 0.5;\n";
        for(std::pair<const htype, Vertex<htype>> & it : dbg) {
            Vertex<htype> &start = it.second;
            for(Edge<htype> &edge : start.getOutgoing()) {
                Vertex<htype> &end = *edge.end();
                printEdge(os, start, edge);
            }
            for(Edge<htype> &edge : start.rc().getOutgoing()) {
                Vertex<htype> &end = *edge.end();
                printEdge(os, start.rc(), edge);
            }
        }
        os << "}\n";
    }
};

template<class htype>
class Component {
private:
    SparseDBG<htype> & graph;
    std::unordered_set<htype, alt_hasher<htype>> v;
    struct EdgeRec {
        Vertex<htype> * start;
        Vertex<htype> * end;
        size_t size;
        size_t cov;
    };

    size_t outDeg(const Vertex<htype> &vert, size_t min_cov) const {
        size_t res = 0;
        for(const Edge<htype> &edge : vert.getOutgoing()) {
            if(edge.getCoverage() >= min_cov) {
                res += 1;
            }
        }
        return res;
    }

    bool isUnbranching(const Vertex<htype> &vert, size_t min_cov) const {
        return v.find(vert.hash()) != v.end() && outDeg(vert, min_cov) == 1 && outDeg(vert.rc(), min_cov) == 1;
    }

    Edge<htype> &getOut(Vertex<htype> &vert, size_t min_cov) {
        for(Edge<htype> &edge : vert.getOutgoing()) {
            if(edge.getCoverage() >= min_cov) {
                return edge;
            }
        }
        VERIFY(false);
    }

    Path<htype> unbranching(Vertex<htype> &vert, Edge<htype> &edge, size_t minCov) {
        std::vector<Edge<htype>*> res;
        res.push_back(&edge);
        Vertex<htype> *cur = edge.end();
        while (cur != &vert && isUnbranching(*cur, minCov)) {
            res.push_back(&getOut(*cur, minCov));
            cur = res.back()->end();
        }
        return Path<htype>(vert, res);
    }

public:
    template<class I>
    Component(SparseDBG<htype> &_graph, I begin, I end) : graph(_graph), v(begin, end) {
    }

    template<class I>
    static Component<htype> neighbourhood(SparseDBG<htype> &graph, I begin, I end, size_t radius, size_t min_coverage = 0) {
        std::unordered_set<htype, alt_hasher<htype>> v;
        std::priority_queue<std::pair<size_t, htype>> queue;
        while(begin != end) {
            queue.emplace(0, *begin);
            ++begin;
        }
        while(!queue.empty()) {
            std::pair<size_t, htype> val = queue.top();
            queue.pop();
            if(v.find(val.second) != v.end())
                continue;
            v.insert(val.second);
            if(val.first > radius)
                continue;
            Vertex<htype> &vert = graph.getVertex(val.second);
            for(Edge<htype> & edge : vert.getOutgoing()) {
                if(edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
            for(Edge<htype> & edge : vert.rc().getOutgoing()) {
                if(edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
        }
        return Component<htype>(graph, v.begin(), v.end());
    }

    static Component<htype> neighbourhood(SparseDBG<htype> &graph, Contig &contig, size_t radius, size_t min_coverage = 0) {
        std::unordered_set<htype, alt_hasher<htype>> v;
        std::priority_queue<std::pair<size_t, htype>> queue;
        std::vector<PerfectAlignment<Contig, Edge<htype>>> als1 = align(graph, contig);
        Contig rc_contig = contig.RC();
        std::vector<PerfectAlignment<Contig, Edge<htype>>> als2 = align(graph, rc_contig);
//        TODO Every edge must have a full sequence stored as a composite sequence
        for(PerfectAlignment<Contig, Edge<htype>> &al : als1) {
            queue.emplace(0, al.seg_to.contig().end()->hash());
        }
        for(PerfectAlignment<Contig, Edge<htype>> &al : als2) {
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
            Vertex<htype> &vert = graph.getVertex(val.second);
            for(Edge<htype> & edge : vert.getOutgoing()) {
                if(edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
            for(Edge<htype> & edge : vert.rc().getOutgoing()) {
                if(edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
        }
        return Component<htype>(graph, v.begin(), v.end());
    }

    void printEdge(std::ostream &os, Vertex<htype> & start, Edge<htype> &edge, const std::string &extra_label = "") {
        Vertex<htype> &end = *edge.end();
        os << "\"";
        if (!start.isCanonical())
            os << "-";
        os << start.hash() % 10000 << "\" -> \"";
        if (!end.isCanonical())
            os << "-";
        os << end.hash() % 10000  << "\" [label=\"" << edge.size() << "(" << edge.getCoverage() << ")";
        if(!extra_label.empty()) {
            os << "\\n"<<extra_label;
        }
        os << "\"]\n";
    }

    void printEdge(std::ostream &os, Path<htype> & path) {
        size_t len = 0;
        size_t cov = 0;
        for(size_t i = 0; i < path.size(); i++) {
            len += path[i].size();
            cov += path[i].intCov();
        }
        Vertex<htype> &start = path.start();
        Vertex<htype> &end = *path.back().end();
        os << "\"";
        if (!start.isCanonical())
            os << "-";
        os << start.hash() % 10000 << "\" -> \"";
        if (!end.isCanonical())
            os << "-";
        os << end.hash() % 10000  << "\" [label=\"" << len << "(" << double(cov) / len << ")\"]\n";
    }

    size_t size() const {
        return v.size();
    }

    void printDot(std::ostream &os, size_t min_cov = 0) {
        const std::function<std::string(Edge<htype> &)> labeler = [](Edge<htype> &) {return "";};
        printCompressedDot(os, min_cov, labeler);
    }

    void printDot(std::ostream &os, const std::function<std::string(Edge<htype> &)> &labeler, size_t min_cov = 0) {
        os << "digraph {\nnodesep = 0.5;\n";
        for(htype vid : v) {
            Vertex<htype> &start = graph.getVertex(vid);
            for(Edge<htype> &edge : start.getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Vertex<htype> &end = *edge.end();
                printEdge(os, start, edge, labeler(edge));
                if(v.find(end.hash()) == v.end()) {
                    printEdge(os, end.rc(), start.rcEdge(edge), labeler(start.rcEdge(edge)));
                }
            }
            for(Edge<htype> &edge : start.rc().getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Vertex<htype> &end = *edge.end();
                printEdge(os, start.rc(), edge, labeler(edge));
                if(v.find(end.hash()) == v.end()) {
                    printEdge(os, end.rc(), start.rc().rcEdge(edge), labeler(start.rc().rcEdge(edge)));
                }
            }
        }
        os << "}\n";
    }

    void printCompressedDot(std::ostream &os, size_t min_cov = 0) {
        os << "digraph {\nnodesep = 0.5;\n";
        for(htype vid : v) {
            Vertex<htype> &start = graph.getVertex(vid);
            if(isUnbranching(start, min_cov))
                continue;
            for(Edge<htype> &edge : start.getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Path<htype> path = unbranching(start, edge, min_cov);
                printEdge(os, path);
                Vertex<htype> &end = *path.back().end();
                if(v.find(end.hash()) == v.end()) {
                    Path<htype> rcpath = path.RC();
                    printEdge(os, rcpath);
                }
            }
            for(Edge<htype> &edge : start.rc().getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Path<htype> path = unbranching(start.rc(), edge, min_cov);
                printEdge(os, path);
                Vertex<htype> &end = *path.back().end();
                if(v.find(end.hash()) == v.end()) {
                    Path<htype> rcpath = path.RC();
                    printEdge(os, rcpath);
                }
            }
        }
        os << "}\n";
    }
};