//
// Created by anton on 19.12.2019.
//

#pragma once

#include <omp.h>
#include "../sequences/contigs.hpp"
#include "ftree.hpp"
#include "../minimap2/cigar.hpp"
#include "../minimap2/minimap_interface.hpp"

template<class U, class V>
struct CigarAlignment {
    Segment<U> seg_from;
    Segment<V> seg_to;
    mm_extra_t *cigar_container; //For compatibility with minimap

    CigarAlignment(RawAlignment &&other, const U &contig_from, const V &contig_to):
                        seg_from(contig_from, other.seg_from.left, other.seg_from.right),
                        seg_to(contig_to, other.seg_to.left, other.seg_to.right), cigar_container(other.cigar_container) {
//        cout << "CigarAlignment From moved raw: " << cigar_container << endl;
        other.cigar_container = nullptr;
    }

    CigarAlignment(const Segment<U> &segFrom, const Segment<V> &segTo, mm_extra_t *cigarContainer) :
                        seg_from(segFrom), seg_to(segTo), cigar_container(cigarContainer) {}

    CigarAlignment(const CigarAlignment &other) = delete;

    CigarAlignment(CigarAlignment &&other) : seg_from(other.seg_from), seg_to(other.seg_to), cigar_container(other.cigar_container) {
        other.cigar_container = nullptr;
    }

    ~CigarAlignment() {
        if (cigar_container != nullptr) {
//            cout << "CigarAlignment destructor" << cigar_container << endl;
            free(cigar_container);
        }
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

template<class U, class V>
class FTreeTranslator {
public:
    AlignmentPiece<U, V> translate(CigarAlignment<U, V> &al) {
        size_t cur_from = al.seg_from.left;
        size_t cur_to = al.seg_to.left;
        const Sequence &contig_from = al.seg_from.contig.seq;
        const Sequence &contig_to = al.seg_to.contig.seq;
        const mm_extra_t * cigar_container = al.cigar_container;
        std::vector<size_t> from;
        std::vector<size_t> to;
        for(size_t i = 0; i < cigar_container->n_cigar; i++){
            size_t block_len = cigar_container->cigar[i] >> 4u;
            if ((cigar_container->cigar[i] & 15u) == 1) {
                cur_from += block_len;
            } else if ((cigar_container->cigar[i] & 15u) == 2) {
                cur_to += block_len;
            } else {
                for(size_t j = 0; j < block_len; j++) {
                    if (contig_from[cur_from + j] == contig_to[cur_to + j]) {
                        from.push_back(cur_from + j);
                        to.push_back(cur_to + j);
                    }
                }
                cur_from += block_len;
                cur_to += block_len;
            }
        }
        return AlignmentPiece<U, V>(al.seg_from, al.seg_to, FTree<size_t>(from), FTree<size_t>(to));
    }

    std::vector<AlignmentPiece<U, V>> translate(std::vector<CigarAlignment<U, V>> &als) {
        std::vector<AlignmentPiece<U, V>> res;
        for(CigarAlignment<U, V> &al: als) {
            res.emplace_back(std::move(translate(al)));
        }
        return std::move(res);
    }

    std::vector<std::vector<AlignmentPiece<U, V>>> translate(std::vector<std::vector<CigarAlignment<U, V>>> &als, size_t threads_num = 1) {
        std::vector<std::vector<AlignmentPiece<U, V>>> res;
        res.resize(als.size());
        omp_set_dynamic(0);
        omp_set_num_threads(threads_num);
#pragma omp parallel
        {
#pragma omp parallel for default(none) shared(als, res)
            for(size_t i = 0; i < als.size(); i++) {
                res[i] = translate(als[i]);
            }
        };
        return std::move(res);
    }
};

