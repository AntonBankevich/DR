#pragma once
#include "dbg/graph_modification.hpp"
#include "dbg/compact_path.hpp"

GraphAlignment correctFromStart(const GraphAlignment &al);
size_t correct_dimers(logging::Logger &logger, RecordStorage &reads_storage, size_t k, size_t threads);