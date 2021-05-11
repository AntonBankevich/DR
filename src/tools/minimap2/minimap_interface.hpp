//
// Created by anton on 30.01.2020.
//

#pragma once
#include <utility>
#include <vector>
#include <string>
#include "minimap.h"
#include "sequences/contigs.hpp"
#include "ksw2.h"

void destroyIndex(std::vector<mm_idx_t *> & ref);

class RawSegment {
public:
    size_t id;
    size_t left;
    size_t right;
    RawSegment(size_t _id, size_t _left, size_t _right) : id(_id), left(_left), right(_right){}
};

class RawAlignment {
public:
    RawSegment seg_from;
    RawSegment seg_to;
    bool rc;
    mm_extra_t *cigar_container; //For compatibility with minimap
    RawAlignment(RawSegment seg_from_, RawSegment seg_to_, bool _rc);
    RawAlignment(RawAlignment &&other) noexcept ;
    RawAlignment(const RawAlignment &other) = delete;
    ~RawAlignment();
};

std::vector<RawAlignment> run_minimap(const std::string & read, size_t read_id, std::vector<mm_idx_t *> & ref, const char * preset = nullptr);
std::vector<std::vector<RawAlignment>> run_minimap(const std::string * reads_from, const std::string * reads_to, size_t read_id,
        std::vector<mm_idx_t *> & ref, const char * preset = nullptr);
std::vector<mm_idx_t *> constructIndex(std::vector<std::string> &ref, size_t threads, const char *preset = nullptr);

struct cigar_pair {
    char type;
    size_t length;
    cigar_pair(char type, size_t len):type(type), length(len) {}
};
std::vector<cigar_pair> align_example(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape, int width);
