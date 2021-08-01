#pragma once

#include "component.hpp"
#include "sparse_dbg.hpp"
#include "graph_alignment_storage.hpp"
#include "graph_printing.hpp"
#include "visualization.hpp"
#include "common/hash_utils.hpp"
#include <utility>

class RepeatResolver {
private:
    SparseDBG &dbg;
    RecordStorage &readStorage;
    std::experimental::filesystem::path dir;

    struct Subdataset {
        Subdataset(Component component, std::experimental::filesystem::path dir) : component(std::move(component)),
                                                                                                 dir(std::move(dir)) {}
        Component component;
        std::experimental::filesystem::path dir;
        bool operator<(const Subdataset &other) const {
            if(component.size() != other.component.size())
                return component.size() > other.component.size();
            return this < &other;
        }
    };
public:
    RepeatResolver(SparseDBG &dbg, RecordStorage &readStorage, std::experimental::filesystem::path dir) :
                dbg(dbg), readStorage(readStorage), dir(std::move(dir)) {
    }

    std::vector<Subdataset> SplitDataset(const std::function<bool(const Edge &)> &is_unique);
    std::vector<Contig> ResolveRepeats(logging::Logger &logger, size_t threads,
                                       const std::function<bool(const Edge &)> &is_unique = [](const Edge &){return false;});
    std::vector<Contig> CollectResults(logging::Logger &logger, size_t threads, const std::vector<Contig> &contigs,
                                       const std::function<bool(const Edge &)> &is_unique = [](const Edge &){return false;});
};

void PrintFasta(const std::vector<Contig> &contigs, const std::experimental::filesystem::path &path);
void PrintAlignments(logging::Logger &logger, size_t threads, std::vector<Contig> &contigs,
                     const RecordStorage &readStorage, size_t k, size_t w,
                     const std::experimental::filesystem::path &dir);
