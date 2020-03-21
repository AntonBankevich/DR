//
// Created by anton on 18.03.2020.
//

#include "aligner.hpp"
#include "contigs.hpp"
#include "../minimap2/minimap_interface.hpp"


Aligner::Aligner(size_t _default_threads) {
    default_threads = _default_threads;
}

void *minimap_worker(void *_shared) {
    auto *shared = (WorkerData *)(_shared);
    size_t tid = __sync_fetch_and_add(&shared->worker_cnt, 1u);
    while (true) {
        size_t left = __sync_fetch_and_add(&shared->index, shared->step);
        size_t right = left + shared->step;
        if (left > shared->read_seqs.size()) {
            return nullptr;
        }
        if (right > shared->read_seqs.size()) {
            right = shared->read_seqs.size();
        }
        std::vector<RawAlignment> results = run_minimap(shared->read_seqs.data() + left, shared->read_seqs.data() + right, left, shared->ref_index.index);
        size_t t = left / shared->step;
        shared->hits[t].reserve(results.size());
        for(RawAlignment &res: results) {
            std::vector<size_t> from;
            std::vector<size_t> to;
            size_t cur_from = res.seg_from.left;
            size_t cur_to = res.seg_to.left;
            const string &contig_from = shared->read_seqs[res.seg_from.id];
            const string &contig_to = shared->ref_index.seqs[res.seg_to.id];
            for(size_t i = 0; i < res.cigar_container->n_cigar; i++){
                size_t block_len = res.cigar_container->cigar[i] >> 4u;
                if ((res.cigar_container->cigar[i] & 15u) == 1) {
                    cur_from += block_len;
                } else if ((res.cigar_container->cigar[i] & 15u) == 2) {
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
            shared->hits[t].push_back(AlignmentRecord(RawSegment(res.seg_from.id, from[0], from[from.size() - 1] + 1),
                                                      RawSegment(res.seg_from.id, to[0], to[to.size() - 1] + 1),
                                                      FTree<size_t>(from), FTree<size_t>(to)));

        }
    }
}

WorkerData::WorkerData(const std::vector<string> &_read_seqs, SimpleAlignmentIndex &_ref_index, size_t _step) :
        index(0), worker_cnt(0), read_seqs(_read_seqs), ref_index(_ref_index), step(_step){
    hits.resize((read_seqs.size()) / _step + 1);//last one may be empty
}

std::vector<AlignmentRecord>
Aligner::doTheWork(const std::vector<string> &read_seqs, SimpleAlignmentIndex &ref_index, size_t thread_num) const {
    WorkerData wd(read_seqs, ref_index, 20);
    std::vector<pthread_t> threads;
    for(size_t i = 0; i < thread_num; i++) {
        threads.push_back(pthread_t());
        pthread_create(&threads.back(), nullptr, &minimap_worker, &wd);
    }
    for(pthread_t &thread: threads)
        pthread_join(thread, nullptr);
    std::vector<AlignmentRecord> result;
    for(auto & hit : wd.hits) {
        std::move(hit.begin(), hit.end(), std::back_inserter(result));
//        result.insert(result.end(), wd.hits[i].begin(), wd.hits[i].end());
    }
    return result;
}
