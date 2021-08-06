#pragma once

#include "dbg/sparse_dbg.hpp"
#include "common/logging.hpp"

std::experimental::filesystem::path CrudeCorrect(logging::Logger &logger, dbg::SparseDBG &dbg,
                  const std::experimental::filesystem::path &dir, const size_t w, const io::Library &reads_lib,
                  size_t threads, size_t threshold);