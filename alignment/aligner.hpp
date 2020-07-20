//
// Created by anton on 18.03.2020.
//

#pragma once


#include "minimap2/minimap_interface.hpp"
#include "sequences/contigs.hpp"
#include "alignment_piece.hpp"
#include "cigar_alignment.hpp"
#include "marked_alignment.hpp"
#include <omp.h>
#include <fstream>

template <class R>
class RawAligner {
public:
    size_t default_threads;
    std::vector<mm_idx_t *> index;
    std::vector<string> ref_seqs;
    const std::vector<const R *> reference;

    explicit RawAligner(const std::vector<const R *> & ref, size_t _default_threads = 1): default_threads(_default_threads), reference(ref) {
        destroyIndex(index);
        ref_seqs.reserve(ref.size());
        for(const R * s: ref) {
            ref_seqs.emplace_back(s->str());
        }
        index = constructIndex(ref_seqs, default_threads);
    }
    ~RawAligner() {
        destroyIndex(index);
    }

    std::string itoa(size_t i) {
        std::stringstream ss;
        ss << i;
        return ss.str();
    }

    template <class T>
    std::vector<std::vector<CigarAlignment<T, R>>> align(const std::vector<const T *> & reads, size_t thread_num = size_t(-1)) {
        if (thread_num == size_t(-1))
            thread_num = default_threads;
        std::vector<string> read_seqs;
        read_seqs.reserve(reads.size());
        for(const T * read: reads) {
            read_seqs.push_back(read->str());
        }
        std::vector<std::vector<CigarAlignment<T, R>>> hits;
        hits.resize(reads.size());
        const size_t step = std::max<size_t>(reads.size() / (thread_num * 5), 1u);
        omp_set_dynamic(0);
        omp_set_num_threads(thread_num);
#pragma omp parallel for default(none) shared(hits, reads, read_seqs, step)
        for(size_t i = 0; i < (reads.size() + step - 1) / step; i++) {
            size_t from = i * step;
            size_t to = std::min(from + step, reads.size());
            VERIFY(to > from)
            std::stringstream ss;
            std::vector<std::vector<RawAlignment>> results = run_minimap(read_seqs.data() + from, read_seqs.data() + to, from, index);
            for(size_t j = 0; j < results.size(); j++) {
                std::vector<CigarAlignment<T, R>> & dump = hits[from + j];
                for(RawAlignment & rawAlignment: results[j]) {
                    dump.emplace_back(std::move(rawAlignment), *reads[from + j], *(reference[rawAlignment.seg_to.id]));
                }
            }
        }
        return std::move(hits);
    }// Should be const method but need to modify minimap for that
};

namespace alignment_recipes {
    template<class U, class V>
    std::vector<std::vector<MarkedAlignment<U, V>>> MarkAlign(const std::vector<const U *> & reads, const std::vector<const V *> & ref, const HMM &hmm, size_t thread_num = 1) {
        RawAligner<U> aligner(ref, thread_num);
        std::vector<std::vector<AlignmentPiece<U, V>>> res;
        if (thread_num == size_t(-1)) {
            thread_num = aligner.default_threads;
        }
        std::cout << "Aligning" << std::endl;
        std::vector<std::vector<CigarAlignment<U, V>>> raw_alignments = aligner.template align<V>(reads, thread_num);
        MarkingTranslator<U, V> translator(hmm);
        std::cout << "Translating" << std::endl;
        return translator.translate(raw_alignments, thread_num);
    }

    template<class U, class V>
    std::vector<std::vector<MarkedAlignment<U, V>>> MarkAlign(const std::vector<U> & reads, const std::vector<V> & ref, const HMM &hmm, size_t thread_num = 1) {
        std::vector<const U*> readlinks;
        for(auto &read : reads)
            readlinks.push_back(&read);
        std::vector<const U*> reflinks;
        for(auto &r : ref)
            reflinks.push_back(&r);
        return MarkAlign(readlinks, reflinks, hmm, thread_num);
    }

    template<class U, class V>
    std::vector<std::vector<CigarAlignment<U, V>>> CigarAlign(const std::vector<const U *> & reads, const std::vector<const V *> & ref, size_t thread_num = 1) {
        RawAligner<U> aligner(ref, thread_num);
        std::vector<std::vector<AlignmentPiece<U, V>>> res;
        if (thread_num == size_t(-1)) {
            thread_num = aligner.default_threads;
        }
        return aligner.align(reads, thread_num);
    }

    template<class U, class V>
    std::vector<std::vector<AlignmentPiece<U, V>>> SimpleAlign(const std::vector<const U *> & reads, const std::vector<const V *> & ref, size_t thread_num = 1) {
        RawAligner<U> aligner(ref, thread_num);
        SimpleAlign(reads, aligner, thread_num);
    }

    template<class U, class V>
    std::vector<std::vector<AlignmentPiece<U, V>>> SimpleAlign(const std::vector<U> & reads, const std::vector<V> & ref, size_t thread_num = 1) {
        std::vector<U*> readlinks;
        for(auto &read : reads)
            readlinks.push_back(&read);
        std::vector<U*> reflinks;
        for(auto &r : ref)
            reflinks.push_back(&r);
        RawAligner<U> aligner(reflinks, thread_num);
        SimpleAlign(readlinks, aligner, thread_num);
    }

    template<class U, class V>
    std::vector<std::vector<AlignmentPiece<U, V>>> SimpleAlign(const std::vector<const U *> & reads, RawAligner<V> &aligner, size_t thread_num = size_t(-1)) {
        std::vector<std::vector<AlignmentPiece<U, V>>> res;
        {
            if (thread_num == size_t(-1)) {
                thread_num = aligner.default_threads;
            }
            std::vector<std::vector<CigarAlignment<U, V>>> raw_alignments = aligner.align(reads, thread_num);
            FTreeTranslator<U, V> translator;
            res = translator.translate(raw_alignments, thread_num);
        }
        return std::move(res);
    }

}

