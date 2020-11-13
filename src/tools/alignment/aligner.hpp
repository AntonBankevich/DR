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
#include <common/oneline_utils.hpp>
#include <common/logging.hpp>

template <class R>
class RawAligner {
public:
    size_t default_threads;
    std::vector<mm_idx_t *> index;
    std::vector<string> ref_seqs;
    const std::vector<const R *> reference;
    const char * preset = nullptr;

    explicit RawAligner(const std::vector<const R *> & ref, size_t _default_threads, const char * _preset):
            default_threads(_default_threads), reference(ref), preset(_preset) {
        destroyIndex(index);
        ref_seqs.reserve(ref.size());
        for(const R * s: ref) {
            ref_seqs.emplace_back(s->str());
        }
        index = constructIndex(ref_seqs, default_threads, preset);
    }

    explicit RawAligner(const std::vector<R> & ref, size_t _default_threads, const char * _preset):
            RawAligner(oneline::map<R, const R*, typename std::vector<R>::const_iterator>(ref.begin(), ref.end(), [](const R & val)->const R* {return &val;}), _default_threads, _preset) {
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
    std::vector<CigarAlignment<T, R>> align(const T & read) {
        std::vector<RawAlignment> raw_results = run_minimap(read.seq.str(), 0, index, preset);
        std::vector<CigarAlignment<T, R>> final_results;
        for(RawAlignment & rawAlignment: raw_results) {
            final_results.emplace_back(std::move(rawAlignment), read, *(reference[rawAlignment.seg_to.id]));
        }
        return std::move(final_results);
    }


    template <class T>
    std::vector<std::vector<CigarAlignment<T, R>>> align1(const std::vector<T> & reads, size_t thread_num = size_t(-1)) {
        std::vector<const T *> read_refs = oneline::map<R, const R*, typename std::vector<R>::const_iterator>(reads.begin(), reads.end(), [](const R & val)->const R* {return &val;});
        return std::move(align(read_refs, thread_num));
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
        size_t step = std::max<size_t>(reads.size() / (thread_num * 5), 1u);
        omp_set_dynamic(0);
        omp_set_num_threads(thread_num);
#pragma omp parallel for default(none) shared(hits, reads, read_seqs, step)
        for(size_t i = 0; i < (reads.size() + step - 1) / step; i++) {
            size_t from = i * step;
            size_t to = std::min(from + step, reads.size());
//            VERIFY(to > from)
            std::stringstream ss;
            std::vector<std::vector<RawAlignment>> results = run_minimap(read_seqs.data() + from, read_seqs.data() + to, from, index, preset);
            for(size_t j = 0; j < results.size(); j++) {
                std::vector<CigarAlignment<T, R>> & dump = hits[from + j];
                for(RawAlignment & rawAlignment: results[j]) {
                    dump.emplace_back(std::move(rawAlignment), *reads[from + j], *(reference[rawAlignment.seg_to.id]));
                }
            }
        }
        return std::move(hits);
    } // Should be const method but need to modify minimap for that

    template <class I>
    std::vector<std::vector<CigarAlignment<typename I::value_type, R>>> alignBlock(I begin, I end) {
        typedef typename I::value_type T;
        std::vector<string> read_seqs;
        I it = begin;
        while (it!= end) {
            read_seqs.push_back((*it).str());
            ++it;
        }
        std::vector<std::vector<CigarAlignment<T, R>>> hits;
        hits.resize(read_seqs.size());
        std::vector<std::vector<RawAlignment>> results = run_minimap(read_seqs.data(), read_seqs.data() + read_seqs.size(), 0, index, preset);
        for(size_t j = 0; j < results.size(); j++) {
            std::vector<CigarAlignment<T, R>> & dump = hits[j];
            for(RawAlignment & rawAlignment: results[j]) {
                dump.emplace_back(std::move(rawAlignment), *(begin + j), *(reference[rawAlignment.seg_to.id]));
            }
        }
        return std::move(hits);
    } // Should be const method but need to modify minimap for that

};

namespace alignment_recipes {
    template<class U, class V>
    std::vector<std::vector<MarkedAlignment<U, V>>> MarkAlign(const std::vector<const U *> & reads, const std::vector<const V *> & ref, const HMM &hmm, const char * preset, size_t thread_num = 1) {
        RawAligner<U> aligner(ref, thread_num, preset);
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
    std::vector<std::vector<MarkedAlignment<U, V>>> MarkAlign(const std::vector<U> & reads, const std::vector<V> & ref, const HMM &hmm, const char * preset, size_t thread_num = 1) {
        std::vector<const U*> readlinks;
        for(auto &read : reads)
            readlinks.push_back(&read);
        std::vector<const U*> reflinks;
        for(auto &r : ref)
            reflinks.push_back(&r);
        return MarkAlign(readlinks, reflinks, hmm, preset, thread_num);
    }

    template<class U, class V>
    std::vector<std::vector<MarkedAlignment<U, V>>> FilterAlign(const std::vector<U> & reads, const std::vector<V> & ref, const HMM &hmm, const char * preset, size_t thread_num = 1) {
        std::vector<const U*> readlinks;
        for(auto &read : reads)
            readlinks.push_back(&read);
        std::vector<const U*> reflinks;
        for(auto &r : ref)
            reflinks.push_back(&r);
        return MarkAlign(readlinks, reflinks, hmm, preset, thread_num);
    }

    template<class U, class V>
    std::vector<std::vector<AlignmentPiece<U, V>>> SimpleAlign(const std::vector<const U *> & reads, const std::vector<const V *> & ref, const char * preset, size_t thread_num = 1) {
        RawAligner<U> aligner(ref, thread_num, preset);
        SimpleAlign(reads, aligner, thread_num);
    }

    template<class U, class V>
    std::vector<std::vector<AlignmentPiece<U, V>>> SimpleAlign(const std::vector<U> & reads, const std::vector<V> & ref,
            const char * preset, size_t thread_num) {
        std::vector<U*> readlinks;
        for(auto &read : reads)
            readlinks.push_back(&read);
        std::vector<U*> reflinks;
        for(auto &r : ref)
            reflinks.push_back(&r);
        RawAligner<U> aligner(reflinks, thread_num, preset);
        SimpleAlign(readlinks, aligner, thread_num);
    }

    template<class U, class V>
    std::vector<std::vector<AlignmentPiece<U, V>>> SimpleAlign(const std::vector<const U *> & reads, RawAligner<V> &aligner, size_t thread_num) {
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

    template<class U, class V>
    std::vector<std::vector<CigarAlignment<U, V>>> SplitAlign(const std::vector<const U *> & reads, const std::vector<const V *> & ref, size_t break_size, size_t min_size, const char * preset, size_t thread_num = 1) {
        RawAligner<U> aligner(ref, thread_num, preset);
        if (thread_num == size_t(-1)) {
            thread_num = aligner.default_threads;
        }
        std::vector<std::vector<CigarAlignment<U, V>>> res = aligner.align(reads, thread_num);
        for (std::vector<CigarAlignment<U, V>> &read_res : res) {
            if (read_res.empty())
                continue;
            std::vector<CigarAlignment<U, V>> splitted;
            for(CigarAlignment<U, V> &al : read_res) {
                std::vector<CigarAlignment<U, V>> tmp = CigarAlignment<U,V>::split(std::move(al), break_size, min_size);
                splitted.insert(splitted.end(), tmp.begin(), tmp.end());
            }
            if(splitted.size() == 1) {
                continue;
            }
            std::function<bool(const CigarAlignment<U, V> &, const CigarAlignment<U, V> &)> compare =
                    [] (const CigarAlignment<U, V> &a, const CigarAlignment<U, V> &b) {
                if(a.seg_to.contig().id != b.seg_to.contig.id())
                    return a.seg_to.contig().id < b.seg_to.contig.id();
                return a.seg_to.right < b.seg_to.right || (a.seg_to.right == b.seg_to.right && a.seg_to.left < b.seg_to.right);
            };
            std::sort(splitted.begin(), splitted.end(), compare);
            std::vector<bool> filter(splitted.size(), false);
            for(size_t i = 0; i < splitted.size(); i++) {
                size_t j = i;
                while(j > 0 && splitted[j - 1].seg_to.contig == splitted[i].seg_to.contig && splitted[j - 1].seg_to.right > splitted[i].seg_to.left) {
                    j--;
                    if(splitted[j])
                        continue;
                    if(splitted[i].deepInter(splitted[j])) {
                        if(splitted[i].size() > splitted[j].size()) {
                            splitted[j] = true;
                        } else {
                            splitted[i] = true;
                        }
                    }
                }
            }
            std::vector<CigarAlignment<U, V>> filtered;
            for(size_t i = 0; i < splitted.size(); i++) {
                if(!filter[i])
                    filtered.emplace_back(std::move(splitted[i]));
            }
            read_res = filtered;
        }
        return std::move(res);
    }


    //This method takes Contig generator as an input, aligns and processes generated reads.
    template<class I, class R>
    void AlignAndProcess(I begin, I end, RawAligner<Contig> &aligner,
                               const std::function<void(size_t, const Contig &, const std::vector<CigarAlignment<Contig, R>> &)>& task,
                               logging::Logger &logger, size_t threads, size_t max_records = size_t(-1)) {
        logger.info() << "Starting parallel alignment and processing" << std::endl;
        omp_set_num_threads(threads);
        size_t bucket_length = 1024 * 1024;
        size_t buffer_size = 1024 * 1024;
        size_t max_length = 1024 * 1024 * 1024;
        size_t total = 0;
        size_t total_len = 0;
        std::vector<Contig> prev_items;
        while (begin != end) {
            std::vector<Contig> items;
            items.reserve(buffer_size);
            std::vector<std::vector<CigarAlignment<Contig, R>>> als(prev_items.size());
#pragma omp parallel default(none) shared(begin, end, items, als, prev_items, total, buffer_size, max_length, bucket_length, aligner, task, threads)
            {
#pragma omp single
                {
#pragma omp task default(none) shared(aligner, begin, end, items, als, buffer_size, max_length, threads, prev_items)
                    {
                        size_t clen = 0;
                        while (begin != end && items.size() < buffer_size && clen < max_length) {
                            items.emplace_back((*begin).makeContig());
                            clen += items.back().size();
                            ++begin;
                        }
                    }
                    size_t step = std::max<size_t>(1u, std::min<size_t>(als.size() / threads / 8u, 1000u));
                    for (size_t i = 0; i < als.size(); i += step) {
#pragma omp task default(none) firstprivate(i), shared(als, aligner, prev_items, step)
                        {
                            std::vector<std::vector<CigarAlignment<Contig, R>>> chunk =
                                    aligner.alignBlock(prev_items.begin() + i,
                                                       prev_items.begin() + std::min(i + step, als.size()));
                            for (size_t j = i; j < i + chunk.size(); j++) {
                                als[j] = std::move(chunk[j - i]);
                            }
                        }
                    }
                }
            }
#pragma omp parallel for default(none) shared(prev_items, als, task, total)
            for (size_t i = 0; i < als.size(); i++) {
                task(total + i, prev_items[i], als[i]);
            }
            if (!prev_items.empty()) {
                size_t clen = 0;
                for (Contig &contig : prev_items) {
                    clen += contig.size();
                }
                total_len += clen;
                total += prev_items.size();
                logger.info() << prev_items.size() << " items of total length " << clen << " processed " << std::endl;
            }
            std::swap(items, prev_items);
            if(total > max_records)
                break;
        }
        {
            std::vector<std::vector<CigarAlignment<Contig, R>>> als(prev_items.size());
#pragma omp parallel default(none) shared(begin, end, als, prev_items, aligner, task, threads)
            {
#pragma omp single
                {
                    size_t step = std::max<size_t>(1u, std::min<size_t>(als.size() / threads / 8u, 1000u));
                    for (size_t i = 0; i < als.size(); i += step) {
#pragma omp task default(none) firstprivate(i), shared(als, aligner, prev_items, step)
                        {
                            std::vector<std::vector<CigarAlignment<Contig, R>>> chunk =
                                    aligner.alignBlock(prev_items.begin() + i,
                                                       prev_items.begin() + std::min(i + step, als.size()));
                            for (size_t j = i; j < i + chunk.size(); j++)
                                als[j] = std::move(chunk[j - i]);
                        }
                    }
                }
            }
#pragma omp parallel for default(none) shared(prev_items, als, task, total)
            for (size_t i = 0; i < als.size(); i++) {
                task(total + i, prev_items[i], als[i]);
            }
            size_t clen = 0;
            for (Contig &contig : prev_items) {
                clen += contig.size();
            }
            total_len += clen;
            total += prev_items.size();
            logger.info() << prev_items.size() << " items of total length " << clen << " processed " << std::endl;
        }
        logger.info() << "Finished parallel alignment and processing. Processed " << total <<
               " sequences with total length " << total_len << std::endl;
    }
}

