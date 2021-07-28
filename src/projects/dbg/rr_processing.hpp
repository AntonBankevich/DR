#pragma once

#include "component.hpp"
#include "sparse_dbg.hpp"
#include "multiplicity_estimation.hpp"

class RepeatResolver {
private:
    const dbg::SparseDBG &dbg;
    const AbstractUniquenessStorage &storage;
public:
    RepeatResolver(const dbg::SparseDBG &dbg, const AbstractUniquenessStorage &storage) : dbg(dbg), storage(storage) {
    }
};