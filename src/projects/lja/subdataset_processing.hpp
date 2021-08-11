#pragma once

#include "../dbg/component.hpp"
#include "../dbg/sparse_dbg.hpp"
#include "../dbg/graph_alignment_storage.hpp"
#include "../dbg/graph_printing.hpp"
#include "../dbg/visualization.hpp"
#include "common/hash_utils.hpp"
#include <utility>

class RepeatResolver {
private:
    SparseDBG &dbg;
    RecordStorage &readStorage;
    std::experimental::filesystem::path dir;

public:
    struct Subdataset {
        Subdataset(size_t id, Component component, std::experimental::filesystem::path dir) :
                            id(id), component(std::move(component)), dir(std::move(dir)) {}
        size_t id;
        Component component;
        std::vector<size_t> reads;
        std::experimental::filesystem::path dir;
        bool operator<(const Subdataset &other) const {
            if(component.size() != other.component.size())
                return component.size() > other.component.size();
            return this < &other;
        }
    };
    RepeatResolver(SparseDBG &dbg, RecordStorage &readStorage, std::experimental::filesystem::path dir) :
                dbg(dbg), readStorage(readStorage), dir(std::move(dir)) {
    }

    std::vector<Subdataset> SplitDataset(const std::function<bool(const Edge &)> &is_unique);
    void prepareDataset(const Subdataset &subdataset);
    std::vector<Contig> ResolveRepeats(logging::Logger &logger, size_t threads,
                                       const std::function<bool(const Edge &)> &is_unique = [](const Edge &){return false;});
    std::vector<Contig> CollectResults(logging::Logger &logger, size_t threads, const std::vector<Contig> &contigs,
                                       const std::experimental::filesystem::path &merging,
                                       const std::function<bool(const Edge &)> &is_unique = [](const Edge &){return false;});
};

void PrintFasta(const std::vector<Contig> &contigs, const std::experimental::filesystem::path &path);
void PrintAlignments(logging::Logger &logger, size_t threads, std::vector<Contig> &contigs,
                     const RecordStorage &readStorage, size_t k,
                     const std::experimental::filesystem::path &dir);
