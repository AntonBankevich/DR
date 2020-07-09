//
// Created by anton on 19.12.2019.
//

#pragma once

#include <omp.h>
#include "sequences/contigs.hpp"
#include "ftree.hpp"
#include "minimap2/cigar.hpp"
#include "minimap2/minimap_interface.hpp"
#include "cigar_alignment.hpp"

template<class U, class V>
struct PositionalAlignment {
    Segment<U> seg_from;
    Segment<V> seg_to;
    std::vector<size_t> positions_from;
    std::vector<size_t> positions_to;

    PositionalAlignment(const Segment<U> & seg_from_, const Segment<V> & seg_to_, std::vector<size_t> &&_positions_from, std::vector<size_t> &&_positions_to):
            seg_from(seg_from_), seg_to(seg_to_), positions_from(_positions_from), positions_to(_positions_to){}

    PositionalAlignment(const U & contig_from, const V contig_to, std::vector<size_t> &&_positions_from, std::vector<size_t> &&_positions_to):
            seg_from(Segment<U>(contig_from, positions_from.front(), positions_from.back())),
            seg_to(Segment<V>(contig_to, positions_to.front(), positions_to.back())),
            positions_from(_positions_from), positions_to(_positions_to){}


    PositionalAlignment(const PositionalAlignment &other) = delete;

    PositionalAlignment(PositionalAlignment<U, V> &&other)  noexcept :
            seg_from(other.seg_from), seg_to(other.seg_to), positions_from(other.positions_from), positions_to(other.positions_to) {
    }

    explicit PositionalAlignment(const CigarAlignment<U, V> &other): seg_from(other.seg_from), seg_to(other.seg_to) {
        size_t cur_from = other.seg_from.left;
        size_t cur_to = other.seg_to.left;
        const Sequence &contig_from = other.seg_from.contig.seq;
        const Sequence &contig_to = other.seg_to.contig.seq;
        const mm_extra_t * cigar_container = other.cigar_container;
        for(size_t i = 0; i < cigar_container->n_cigar; i++){
            size_t block_len = cigar_container->cigar[i] >> 4u;
            if ((cigar_container->cigar[i] & 15u) == 1) {
                cur_from += block_len;
            } else if ((cigar_container->cigar[i] & 15u) == 2) {
                cur_to += block_len;
            } else {
                for(size_t j = 0; j < block_len; j++) {
                    if (contig_from[cur_from + j] == contig_to[cur_to + j]) {
                        positions_from.push_back(cur_from + j);
                        positions_to.push_back(cur_to + j);
                    }
                }
                cur_from += block_len;
                cur_to += block_len;
            }
        }
    }

    PositionalAlignment<U, V> subalignment(size_t from, size_t to) const {
        return PositionalAlignment<U, V>(seg_from.contig, seg_to.contig(),
                                  std::vector<size_t>(positions_from.begin() + from, positions_from.begin() + to + 1),
                                  std::vector<size_t>(positions_to.begin() + from, positions_to.begin() + to + 1));
    }

    std::vector<string> treestringForm() const {
        std::vector<char> l1;
        std::vector<char> l2;
        std::vector<char> diff;
        for(size_t i = 0; i + 1 < positions_from.size(); i++) {
            const size_t len1 = positions_from[i + 1] - positions_from[i] - 1;
            const size_t len2 = positions_to[i + 1] - positions_to[i] - 1;
            for(size_t j = 0; j < std::min(len1, len2) + 1; j++) {
                l1.push_back(seg_from.contig[positions_from[i + j]]);
                l2.push_back(seg_to.contig[positions_to[i + j]]);
            }
            for(size_t j = std::min(len1, len2) + 1; j <= len1; j++) {
                l1.push_back(seg_from.contig[positions_from[i + j]]);
                l2.push_back('-');
            }
            for(size_t j = std::min(len1, len2) + 1; j <= len2; j++) {
                l1.push_back('-');
                l2.push_back(seg_to.contig[positions_to[i + j]]);
            }
        }
        l1.push_back(seg_from.contig[positions_from.back()]);
        l2.push_back(seg_to.contig[positions_to.back()]);
        for(size_t i = 0; i < l1.size(); i++) {
            if(l1[i] == l2[i])
                diff.push_back('|');
            else
                diff.push_back('-');
        }
        return {std::string(l1.begin(), l1.end()), std::string(diff.begin(), diff.end()), std::string(l2.begin(), l2.end())};
    }
};


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
    if (!view.scomp.empty())
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

    explicit  AlignmentPiece(PositionalAlignment<U, V> &&other):
            seg_from(other.seg_from), seg_to(other.seg_to), ftree_from(other.positions_from), ftree_to(other.positions_to){}

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



template<class U, class V>
class FTreeTranslator: public AlignmentTranslator<U, V, AlignmentPiece<U, V>> {
public:
    FTreeTranslator() = default;

    virtual AlignmentPiece<U, V> translateOne(CigarAlignment<U, V> &&al) const {
        return AlignmentPiece<U, V>(PositionalAlignment<U, V>(al));
    }

    virtual ~FTreeTranslator()  = default;
};


