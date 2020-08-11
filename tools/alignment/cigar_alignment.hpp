//
// Created by anton on 08.07.2020.
//

#pragma once

#include <vector>
#include <omp.h>
#include "sequences/contigs.hpp"
#include "minimap2/cigar.hpp"
#include "minimap2/minimap_interface.hpp"

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

    CigarAlignment(CigarAlignment &&other)  noexcept : seg_from(other.seg_from), seg_to(other.seg_to), cigar_container(other.cigar_container) {
        other.cigar_container = nullptr;
    }

    ~CigarAlignment() {
        if (cigar_container != nullptr) {
//            cout << "CigarAlignment destructor" << cigar_container << endl;
            free(cigar_container);
//            cout << "CigarAlignment destructor done" << cigar_container << endl;
        }
    }
};

template<class U, class V, class A>
class AlignmentTranslator {
public:
    AlignmentTranslator() = default;

    virtual A translateOne(CigarAlignment<U, V> &&al) const = 0;

    std::vector<A> translate(std::vector<CigarAlignment<U, V>> &als) {
        std::vector<A> res;
        for(CigarAlignment<U, V> &al: als) {
            res.emplace_back(std::move(translateOne(std::move(al))));
        }
        return std::move(res);
    }

    std::vector<std::vector<A>> translate(std::vector<std::vector<CigarAlignment<U, V>>> &als, size_t threads_num = 1) {
        std::vector<std::vector<A>> res;
        res.resize(als.size());
        omp_set_dynamic(0);
        omp_set_num_threads(1);
#pragma omp parallel for default(none) shared(als, res)
        for(size_t i = 0; i < als.size(); i++) {
            res[i] = translate(als[i]);
        }
        return std::move(res);
    }

    virtual ~AlignmentTranslator() = default;
};