//
// Created by anton on 7/22/20.
//

#pragma once
#include <vector>
#include <numeric>
#include <sequences/sequence.hpp>
#include <unordered_map>
#include "rolling_hash.hpp"
#include "common/omp_utils.hpp"
#include "common/logging.hpp"

template<class htype>
class ExtensionsRecord {
private:
    std::vector<std::pair<Sequence, htype>> rightExtensions;
    std::vector<std::pair<Sequence, htype>> leftExtensions;
    omp_lock_t writelock{};
public:
    Sequence seq;

    void clear() {
        rightExtensions.clear();
        leftExtensions.clear();
    }

    static void addSequence(const Sequence &s, std::vector<std::pair<Sequence, htype>> &extensions, htype next) {
        if(s.empty()) {
            return;
        }
        for(std::pair<Sequence, htype> & rightExtension : extensions) {
            if (s.size() <= rightExtension.first.size()) {
                if (s == rightExtension.first.Subseq(0, s.size())) {
                    return;
                }
            } else if (s.Subseq(0, rightExtension.first.size()) == rightExtension.first) {
                rightExtension = std::make_pair(Sequence(s.str()), next);
                return;
            }
        }
        extensions.emplace_back(Sequence(s.str()), next);
    }
public:
    ExtensionsRecord() {
        omp_init_lock(&writelock);
    }

    void setSequence(const Sequence &_seq) {
        if(seq.empty()) {
            seq = Sequence(_seq.str());
        }
    }

    void addRightExtension(const Sequence &s) {
        omp_set_lock(&writelock);
        addSequence(s, rightExtensions);
        omp_unset_lock(&writelock);
    }

    void addLeftExtension(const Sequence &s) {
        omp_set_lock(&writelock);
        addSequence(s, leftExtensions);
        omp_unset_lock(&writelock);
    }

    void addExtension(const Sequence &left, const Sequence &right, htype prev, htype next) {
        omp_set_lock(&writelock);
        addSequence(right, rightExtensions, next);
        addSequence(left, leftExtensions, prev);
        omp_unset_lock(&writelock);
    }

    std::vector<std::pair<Sequence, htype>> &getRightExtensions() {
        return rightExtensions;
    }
    std::vector<std::pair<Sequence, htype>> &getLeftExtensions() {
        return leftExtensions;
    }

    bool isJunction() const {
        return leftExtensions.size() != 1 || rightExtensions.size() != 1;
    }
};

template<typename htype>
class SparseDBG {
private:
public:
//    TODO: replace with perfect hash map? It is parallel, maybe faster and compact.
    std::unordered_map<htype, ExtensionsRecord<htype>> v;
    const RollingHash<htype> hasher;
    size_t w;
    SparseDBG(std::vector<htype> &&_hashs, RollingHash<htype> _hasher, size_t _w) : hasher(_hasher), w(_w) {
        std::vector<htype> hashs(_hashs);
        for(htype hash: hashs) {
            v[hash] = ExtensionsRecord<htype>();
        }
    }

    SparseDBG(SparseDBG<htype> &&other): hasher(other.hasher) {
        std::swap(v, other.v);
        w = other.w;
    }

    void addVertex(htype h) {
        v[h] = ExtensionsRecord<htype>();
    }

    std::vector<KWH<htype>> extractMinimizers(const Sequence &seq) const {
        std::vector<KWH<htype>> res;
        KWH<htype> kwh(hasher, seq, 0);
        while(true) {
            if(v.find(kwh.hash) != v.end()) {
                res.emplace_back(kwh);
            }
            if (!kwh.hasNext())
                break;
            kwh = kwh.next();
        }
        return std::move(res);
    }

    void processRead(const Sequence & _seq) {
        Sequence seq(_seq.str());
        std::vector<KWH<htype>> kmers = extractMinimizers(seq);
        std::vector<std::pair<size_t, htype>> pos;
        pos.push_back(std::make_pair(0, htype(-1)));
        for(size_t i = 0; i < kmers.size(); i++) {
            pos.push_back(std::make_pair(kmers[i].pos, kmers[i].hash));
        }
        pos.push_back(std::make_pair(seq.size() - hasher.k, htype(-1)));
        for(size_t i = 1; i + 1 < pos.size(); i++) {
            ExtensionsRecord<htype> &rec = v[pos[i].second];
//            TODO: if too memory heavy save only some of the labels
            rec.setSequence(seq.Subseq(pos[i].first, pos[i].first + hasher.k));
            rec.addExtension(!seq.Subseq(pos[i - 1].first, pos[i].first),
                    seq.Subseq(pos[i].first + hasher.k, pos[i + 1].first + hasher.k),
                    pos[i - 1].second, pos[i + 1].second);
        }
    }

    Sequence walkForward(htype start_hash, ExtensionsRecord<htype> &rec, size_t t) {
        std::vector<Sequence> res;
        if(rec.getLeftExtensions().size() > 0) {
            res.push_back(!(rec.getLeftExtensions()[0].first));
        }
        res.push_back(rec.seq);
        res.push_back(rec.getRightExtensions()[t].first);
        htype hash = rec.getRightExtensions()[t].second;
        while(hash != htype(-1) && hash != start_hash) {
            ExtensionsRecord<htype> &next = v.find(hash)->second;
            if (next.isJunction()) {
                if(next.getRightExtensions().size() > 0) {
                    res.push_back(next.getRightExtensions()[0].first);
                }
                break;
            }
            res.push_back(next.getRightExtensions()[0].first);
            hash = next.getRightExtensions()[0].second;
        }
        return Sequence::Concat(res);
    }

    void print() {
        for(auto &val: v) {
            ExtensionsRecord<htype> &tmp = val.second;
            std::cout << size_t(val.first) << " " << tmp.seq << " >: ";
            for(std::pair<Sequence, htype> &out: tmp.getRightExtensions()) {
                std::cout << out.first << " " << size_t(out.second) << "; ";
            }
            std::cout << " <: ";
            for(std::pair<Sequence, htype> &out: tmp.getLeftExtensions()) {
                std::cout << !out.first << " " << size_t(out.second) << "; ";
            }
            std::cout << std::endl;
        }
    }
    void printStats() {
        std::vector<size_t> arr(10);
        size_t n11 = 0;
        size_t n01 = 0;
        for(auto &val: v) {
            ExtensionsRecord<htype> &tmp = val.second;
            arr[std::min(arr.size() - 1, tmp.getLeftExtensions().size())] += 1;
            arr[std::min(arr.size() - 1, tmp.getRightExtensions().size())] += 1;
            if (tmp.getLeftExtensions().size() == 1 && tmp.getRightExtensions().size() == 1) {
                n11 += 1;
            }
            if (tmp.getLeftExtensions().size() + tmp.getRightExtensions().size() == 1) {
                n01 += 1;
            }
        }
        size_t total = std::accumulate(arr.begin(), arr.end(), size_t(0));
        std::cout << "Sparse graph statistics:" << std::endl << "Total edges: " << total << std::endl;
        std::cout << "Number of end vertices: " << n01 << std::endl;
        std::cout << "Number of unbranching vertices: " << n11 << std::endl << "Distribution of degrees:" << std::endl;
        for(size_t i = 0; i < arr.size(); i++) {
            std::cout << i << " " << arr[i] << std::endl;
        }
    }

};

template<typename htype, class Iterator>
void fillSparseDBGEdges(SparseDBG<htype> &sdbg, const Time &time, Iterator begin, Iterator end, size_t threads,
        const RollingHash<htype> &hasher, const size_t w) {
    std::cout << time.get() << "Starting to fill edges" << std::endl;
    const size_t buffer_size = 1000000000;
    while(begin != end) {
        size_t tlen = 0;
        std::cout << time.get() << "Starting new round" << std::endl;
        std::vector<Sequence> reads;
        reads.reserve(1000000);
#pragma omp parallel default(none) shared(hasher, begin, end, tlen, buffer_size, std::cout, time, reads, w, sdbg)
        {
#pragma omp single
            {
                while (begin != end and tlen < buffer_size) {
                    reads.push_back(*begin);
                    ++begin;
                    tlen += reads.back().size();
                    if(reads.back().size() < hasher.k + 1) {
                        continue;
                    }
                    size_t index = reads.size() - 1;
#pragma omp task default(none) shared(reads, index, sdbg)
                    {
                        sdbg.processRead(reads[index]);
                    }
                }
                std::cout << time.get() << tlen  << " nucleotides in " << reads.size() <<
                          " sequences were collected. Processing in progress  " << std::endl;
            }
        }
    }
    std::cout << time.get() << "Sparse graph edges filled." << std::endl;
}


template<typename htype>
SparseDBG<htype> constructSparseDBGFromReads(Time &time, const string &reads_file, size_t threads, const RollingHash<htype> &hasher,
                                       std::vector<htype> &&hash_list, const size_t w) {
    std::cout << time.get() << "Starting construction of sparse de Bruijn graph" << std::endl;
    SparseDBG<htype> sdbg(std::move(hash_list), hasher, w);
    std::cout << time.get() << "Vertex map constructed." << std::endl;
    io::SeqReader reader(reads_file);
    fillSparseDBGEdges(sdbg, time, reader.seqbegin(), reader.seqend(), threads, hasher, w);
    return std::move(sdbg);
}


template<typename htype>
void tieTips(const Time &time, SparseDBG<htype> &sdbg, size_t w, size_t threads) {
    std::cout << time.get() << " Collecting tips " << std::endl;
    ParallelRecordCollector<htype> new_minimizers(threads);
//    TODO reduce memory consumption!! A lot of duplicated k-mer storing
    ParallelRecordCollector<Sequence> old_edges(threads);
#pragma omp parallel default(none) shared(sdbg, old_edges, new_minimizers)
    {
#pragma omp single
        {
            for (auto &it: sdbg.v) {
                ExtensionsRecord<htype> &rec = it.second;
                VERIFY(!rec.seq.empty());
#pragma omp task default(none) shared(sdbg, old_edges, new_minimizers, rec)
                {
                    for (size_t i = 0; i < rec.getRightExtensions().size(); i++) {
                        const auto &ext = rec.getRightExtensions()[i];
                        Sequence seq = rec.seq + ext.first;
                        old_edges.add(seq);
                        if (ext.second == htype(-1)) {
                            htype hash = sdbg.hasher.hash(seq, ext.first.size());
                            new_minimizers.emplace_back(hash);
                        }
                    }
                    for (size_t i = 0; i < rec.getLeftExtensions().size(); i++) {
                        const auto &ext = rec.getLeftExtensions()[i];
                        Sequence seq = !ext.first + rec.seq;
                        old_edges.add(seq);
                        if (ext.second == htype(-1)) {
                            htype hash = sdbg.hasher.hash(seq, 0);
                            new_minimizers.emplace_back(hash);
                        }
                    }
                    rec.clear();
                }
            }
        }
    }
    std::cout << time.get() << "Added " << new_minimizers.size() << " artificial minimizers from tips." << std::endl;
    std::cout << time.get() << "Collected " << old_edges.size() << " old edges." << std::endl;
    for(auto it = new_minimizers.begin(); it != new_minimizers.end(); ++it) {
        sdbg.addVertex(*it);
    }
    std::cout << time.get() << "New minimizers added to sparse graph." << std::endl;
    std::cout << time.get() << "Refilling graph with edges." << std::endl;
    fillSparseDBGEdges(sdbg, time, old_edges.begin(), old_edges.end(), threads, sdbg.hasher, w);
    std::cout << time.get() << "Finished fixing sparse de Bruijn graph." << std::endl;
}

