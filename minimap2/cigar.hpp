#pragma once

// minimap2 alignment
typedef struct {
    uint32_t capacity;                  // the capacity of cigar[]
    int32_t dp_score, dp_max, dp_max2;  // DP score; score of the max-scoring segment; score of the best alternate mappings
    uint32_t n_ambi:30, trans_strand:2; // number of ambiguous bases; transcript strand: 0 for unknown, 1 for +, 2 for -
    uint32_t n_cigar;                   // number of cigar operations in cigar[]
    uint32_t cigar[];
} mm_extra_t;