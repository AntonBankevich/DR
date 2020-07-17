// To compile:
//   gcc -g -O2 example.c libminimap2.a -lz

#include <cstdlib>
#include <cassert>
#include <cstdio>
#include <pthread.h>
#include <zlib.h>
#include <vector>
#include <string>
#include <iostream>
#include "minimap.h"
#include "kseq.h"
#include <ctime>
#include <sstream>
#include "minimap_interface.hpp"



char** tochar(const std::vector<std::string> & val) {
    char ** res = (char**)malloc(sizeof(char*) * val.size());
    for (size_t i = 0; i < val.size(); i++)
        res[i] = strdup(val[i].c_str());
    return res;
}


void destroyIndex(std::vector<mm_idx_t *> & ref) {
    for (mm_idx_t *mi: ref) {
        mm_idx_destroy(mi);
    }
}

RawAlignment::RawAlignment(RawSegment seg_from_, RawSegment seg_to_, bool _rc) : seg_from(seg_from_), seg_to(seg_to_), rc(_rc) {
    cigar_container = nullptr;
}

RawAlignment::RawAlignment(RawAlignment &&other) noexcept : seg_from(other.seg_from), seg_to(other.seg_to), rc(other.rc), cigar_container(other.cigar_container) {
//    std::cout << "RawAlignment move: " << cigar_container << std::endl;
    other.cigar_container = nullptr;
}

RawAlignment::~RawAlignment() {
    if (cigar_container != nullptr) {
//        std::cout << "RawAlignment destructor: " << cigar_container << std::endl;
        free(cigar_container);
//        std::cout << "RawAlignment destructor done: " << cigar_container << std::endl;
    }
}


std::vector<mm_idx_t *> constructIndex(std::vector<std::string> &ref, size_t threads) {
    if(threads < 3)
        threads = 3;
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    mm_verbose = 2; // disable message output to stderr
    mm_set_opt(nullptr, &iopt, &mopt);
    mopt.flag |= MM_F_CIGAR; // perform alignment

    char ** ref_seq = (char**)malloc(sizeof(char*) * ref.size());
    char ** ref_name = (char**)malloc(sizeof(char*) * ref.size());
    for (size_t i = 0; i < ref.size(); i++) {
        ref_seq[i] = strdup(ref[i].c_str());
        ref_name[i] = strdup(std::to_string(i).c_str());
    }
    runtime_sequences rs;
    rs.s = ref_seq;
    rs.names = ref_name;
    rs.cur = 0;
    rs.size = ref.size();
    mm_idx_reader_t *r = mm_idx_reader_open(nullptr, &rs, &iopt, 0);
    std::vector<mm_idx_t *> res;
    mm_idx_t *mi = nullptr;
    while ((mi = mm_idx_reader_read(r, int(threads))) != nullptr) {
        res.push_back(mi);
    }
    free(ref_seq);
    free(ref_name);
    mm_idx_reader_close(r); // close the index reader
    return res;
}

std::vector<std::vector<RawAlignment>> run_minimap(const std::string * reads_from, const std::string * reads_to, size_t read_id, std::vector<mm_idx_t *> & ref)
{
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    mm_verbose = 2; // disable message output to stderr
    mm_set_opt(0, &iopt, &mopt);
    mopt.flag |= MM_F_CIGAR; // perform alignment
    std::vector<std::vector<RawAlignment>> result;
    result.resize(reads_to - reads_from);
    size_t total = 0;
    for (mm_idx_t * mi: ref) {
        mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
        mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
//		gzrewind(f);
//		kseq_rewind(ks);
        size_t cnt = read_id;
//		while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
        for (size_t i = 0; i < reads_to - reads_from; i++) {
            mm_reg1_t *reg;
            int j, n_reg;
            const std::string &read = *(reads_from + i);
            char * read_seq = strdup(read.c_str());
            reg = mm_map(mi, read.size(), read_seq, &n_reg, tbuf, &mopt, nullptr); // get all hits for the query
            for (j = 0; j < n_reg; ++j) {
                mm_reg1_t *hit = &reg[j];
                if (hit->rev != 0)
                    continue;
                total += hit->p->n_cigar;
                RawSegment seg_from(cnt, hit->qs, hit->qe);
                RawSegment seg_to(hit->rid, hit->rs, hit->re);
                RawAlignment al(seg_from, seg_to, bool(hit->rev));
                result[i].emplace_back(seg_from, seg_to, bool(hit->rev));
                result[i].back().cigar_container = hit->p;
            }
            free(reg);
            cnt += 1;
        }
        mm_tbuf_destroy(tbuf);
    }
//    Control over this memory is transferred to RawSequence objects in the result container
//	kseq_destroy(ks); // close the query file
//	gzclose(f);
//    for (size_t i = 0; i < ref.size(); i++) {
//        free(ref_name[i]);
//        free(ref_seq[i]);
//    }
    return result;
}
