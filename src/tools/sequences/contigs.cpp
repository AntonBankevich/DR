//
// Created by anton on 19.12.2019.
//

#include "contigs.hpp"

bool StringContig::homopolymer_compressing = false;
size_t StringContig::min_dimer_to_compress = 1000000000;
size_t StringContig::max_dimer_size = 1000000000;
size_t StringContig::dimer_step = 1;
