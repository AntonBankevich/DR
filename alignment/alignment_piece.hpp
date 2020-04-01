//
// Created by anton on 19.12.2019.
//

#pragma once

#include "contigs.hpp"
#include "ftree.hpp"


struct AlignmentView {
    string sfrom;
    string sto;
    string scomp;
    AlignmentView(string &&sfrom, string &&sto): sfrom(sfrom), sto(sto) {
        scomp = "";
    }

    AlignmentView &addComp() {
        std::stringstream ss;
        for(size_t i =0; i < sfrom.size(); i++) {
            if (sfrom[i] == sto[i])
                ss << "|";
            else
                ss << ".";
        }
        scomp = ss.str();
        return *this;
    }
};

inline std::ostream& operator<<(std::ostream& os, const AlignmentView& view)
{
    os << view.sfrom << std::endl;
    if (view.scomp != "")
        os << view.scomp << std::endl;
    os << view.sto << std::endl;
    return os;
}

template <class U, class V>
class AlignmentPiece;

template <class U, class V>
inline std::ostream& operator<<(std::ostream& os, const AlignmentPiece<U, V>& al);

template <class U, class V>
class AlignmentPiece {
private:
    Segment<U> seg_from;
    Segment<V> seg_to;
    FTree<size_t> ftree_from;
    FTree<size_t> ftree_to;
    friend std::ostream& operator<< <U, V>(std::ostream& os, const AlignmentPiece<U, V>& al);
public:

    AlignmentPiece(const Segment<U> & seg_from_, const Segment<V> & seg_to_, std::vector<size_t> &&positions_from, std::vector<size_t> &&positions_to):
            seg_from(seg_from_), seg_to(seg_to_), ftree_from(positions_from), ftree_to(positions_to){}

    AlignmentPiece(const Segment<U> & seg_from_, const Segment<V> & seg_to_, FTree<size_t> &&positions_from, FTree<size_t> &&positions_to):
            seg_from(seg_from_), seg_to(seg_to_), ftree_from(positions_from), ftree_to(positions_to){}


    AlignmentView view() const {
        std::stringstream ss_from, ss_to;
        std::vector<size_t> pos_from = ftree_from.get(0, ftree_from.size());
        std::vector<size_t> pos_to = ftree_from.get(0, ftree_to.size());
        for(size_t i = 0; i + 1 < pos_from.size(); i++){
            ss_from << seg_from.contig[pos_from[i]];
            ss_to << seg_to.contig[pos_to[i]];
            size_t d_from = pos_from[i + 1] - pos_from[i];
            size_t d_to = pos_to[i + 1] - pos_to[i];
            for(size_t j = 1; j < std::min(d_from, d_to); j++) {
                ss_from << seg_from.contig[pos_from[i] + j];
                ss_to << seg_to.contig[pos_to[i] + j];
            }
            for(size_t j = std::min(d_from, d_to); j < d_from; j++) {
                ss_from << seg_from.contig[pos_from[i + j + 1]];
                ss_to << "-";
            }
            for(size_t j = std::min(d_from, d_to); j < d_to; j++) {
                ss_from << "-";
                ss_to << seg_to.contig[pos_to[i + j + 1]];
            }
        }
        ss_from << seg_from.contig[pos_from.back()];
        ss_to << seg_to.contig[pos_to.back()];
        return AlignmentView(ss_from.str(), ss_to.str());
    }
};

template <class U, class V>
inline std::ostream& operator<<(std::ostream& os, const AlignmentPiece<U, V>& al)
{
    os << "(" << al.seg_from << "->" <<al.seg_to << ")";
    return os;
}