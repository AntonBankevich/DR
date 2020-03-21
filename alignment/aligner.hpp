//
// Created by anton on 18.03.2020.
//

#pragma once


#include <pthread.h>
#include "../minimap2/minimap_interface.hpp"
#include "contigs.hpp"
#include "alignment_piece.hpp"

//class NoTranslator {
//public:
//    typedef RawAlignment AType;
//    static AType translate(RawAlignment &&alignment);
//};

//class EqualTranslator {
//public:
//    typedef RawAlignment AType;
//    static AType translate(RawAlignment &&alignment) {
//
//    }
//};


class SimpleAlignmentIndex {
public:
    std::vector<mm_idx_t *> index;
    std::vector<string> seqs;
    explicit SimpleAlignmentIndex(std::vector<mm_idx_t *> &&_index, std::vector<string> &&_seqs): index(_index), seqs(_seqs) {
    }
};

template <class T>
class AlignmentIndex: public SimpleAlignmentIndex {
public:
    const std::vector<T *> & reference;

    AlignmentIndex(const std::vector<T *> &_reference, std::vector<string> &&_seqs, std::vector<mm_idx_t *> &&_index)
            : SimpleAlignmentIndex(std::move(_index), std::move(_seqs)), reference(reference) {
    }
};

class AlignmentRecord {
public:
    RawSegment from;
    RawSegment to;
    FTree<size_t> pos_from;
    FTree<size_t> pos_to;
    AlignmentRecord(RawSegment from, RawSegment to, FTree<size_t> &&posFrom, FTree<size_t> &&posTo)
            : from(from), to(to), pos_from(posFrom), pos_to(posTo) {
    }
};

class Aligner {
private:
    size_t default_threads;
public:
    explicit Aligner(size_t _default_threads = 1);

    std::vector<AlignmentRecord> doTheWork(const std::vector<string> &read_seqs, SimpleAlignmentIndex &ref_index, size_t thread_num) const;

    template <class T>
    AlignmentIndex<T> prepareIndex(const std::vector<T *> & ref, size_t threads = size_t(-1)) {
        std::vector<mm_idx_t *> index;
        destroyIndex(index);
        if (threads == size_t(-1))
            threads = default_threads;
        std::vector<string> seqs;
        seqs.reserve(ref.size());
        for(T *s: ref) {
            seqs.emplace_back(s->str());
        }
        index = constructIndex(seqs, threads);
        return AlignmentIndex<T>(ref, std::move(seqs), std::move(index));
    }

    template <class U, class V>
    std::vector<AlignmentPiece<U, V>> alignAndTranslate(const std::vector<U *> & reads, AlignmentIndex<V> & index, size_t thread_num = size_t(-1)) const {
        if (thread_num == size_t(-1))
            thread_num = default_threads;
        std::vector<string> read_seqs;
        read_seqs.reserve(reads.size());
        for(U * read: reads) {
            read_seqs.push_back(read->str());
        }
        std::vector<AlignmentRecord> res = doTheWork(read_seqs, index, thread_num);
        std::vector<AlignmentPiece<U, V>> final;
        final.reserve(res.size());
        for(AlignmentRecord &rec: res) {
            final.emplace_back(Segment<U>(*reads[rec.from.id], rec.from.left, rec.from.right),
                               Segment<U>(*reads[rec.to.id], rec.to.left, rec.to.right),
                               std::move(rec.pos_from), std::move(rec.pos_to));
        }
        return final;
    }// Should be const method but need to modify minimap for that
};

struct WorkerData {
    size_t index;
    size_t worker_cnt;
    const std::vector<string> &read_seqs;
    SimpleAlignmentIndex &ref_index;
    size_t step;
    std::vector<std::vector<AlignmentRecord>> hits;
    WorkerData(const std::vector<string> &read_seqs, SimpleAlignmentIndex &_index, size_t _step);
};

