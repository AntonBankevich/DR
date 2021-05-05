//
// Created by anton on 8/3/20.
//

#pragma once

#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include "sequences/seqio.hpp"
#include "common/logging.hpp"
#include "common/omp_utils.hpp"

std::vector<htype> constructMinimizers(logging::Logger &logger, const io::Library &reads_file, size_t threads,
                                       const RollingHash &hasher, const size_t w);

