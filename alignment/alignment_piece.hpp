//
// Created by anton on 19.12.2019.
//

#pragma once

#include "contigs.hpp"
#include "ftree.hpp"

template <class U, class V>
class AlignmentPiece {
private:
    Segment<U> seg_from;
    Segment<V> seg_to;
    FTree<size_t> ftree_from;
    FTree<size_t> ftree_to;
public:
    AlignmentPiece(const Segment<U> & seg_from_, const Segment<V> & seg_to_, std::vector<size_t> &&positions_from, std::vector<size_t> &&positions_to):
            seg_from(seg_from_), seg_to(seg_to_), ftree_from(positions_from), ftree_to(positions_to){}

    AlignmentPiece(const Segment<U> & seg_from_, const Segment<V> & seg_to_, FTree<size_t> &&positions_from, FTree<size_t> &&positions_to):
            seg_from(seg_from_), seg_to(seg_to_), ftree_from(positions_from), ftree_to(positions_to){}


    std::pair<std::string, std::string> view() const {
        std::stringstream ss_from, ss_to;
        std::vector<size_t> pos_from = ftree_from.get(0, ftree_from.size());
        std::vector<size_t> pos_to = ftree_from.get(0, ftree_to.size());
        for(size_t i = 0; i + 1 < pos_from.size(); i++){
            ss_from << seg_from.contig[pos_from[i]];
            ss_to << seg_from.contig[pos_to[i]];
            size_t d_from = pos_from[i + 1] - pos_from[i] - 1;
            size_t d_to = pos_to[i + 1] - pos_to[i] - 1;
            for(size_t j = 0; j < std::min(d_from, d_to); j++) {
                ss_from << seg_from.contig[pos_from[i + j + 1]];
                ss_to << seg_from.contig[pos_to[i + j + 1]];
            }
            for(size_t j = std::min(d_from, d_to); j < d_from; j++) {
                ss_from << seg_from.contig[pos_from[i + j + 1]];
                ss_to << "-";
            }
            for(size_t j = std::min(d_from, d_to); j < d_to; j++) {
                ss_from << "-";
                ss_to << seg_from.contig[pos_to[i + j + 1]];
            }
        }
        return std::make_pair(ss_from.str(), ss_to.str());
    }
};