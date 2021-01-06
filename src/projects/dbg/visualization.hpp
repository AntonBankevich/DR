#pragma once
#include "sparse_dbg.hpp"
#include "sequences/contigs.hpp"
#include <unordered_map>


template<class U, class V>
class PerfectAlignment {
public:
    Segment<U> seg_from;
    Segment<V> seg_to;
    PerfectAlignment(const Segment<U> & seg_from_, const Segment<V> &seg_to_) : seg_from(seg_from_), seg_to(seg_to_) {
    }
};

template<typename htype>
class GraphAlignmentStorage {
private:
    std::unordered_map<const Edge<htype> *, std::vector<PerfectAlignment<Contig, Edge<htype>>>> alignments;
    std::vector<Contig*> stored_contigs;
    SparseDBG<htype> & dbg;

    void innerFill(const Contig &old_contig) {
        stored_contigs.emplace_back(new Contig(old_contig));
        Contig &contig = *stored_contigs.back();
        GraphAlignment<htype> path = dbg.align(contig.seq);
        size_t pos = 0;
        for(Segment<Edge<htype>> &al : path) {
            alignments[&al.contig()].emplace_back(Segment<Contig>(contig, pos, pos + al.size()), al);
            pos += al.size();
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