#pragma once

#include "multiplicity_estimation.hpp"
#include "sparse_dbg.hpp"
#include "compact_path.hpp"
#include <experimental/filesystem>

void MultCorrect(SparseDBG &sdbg, logging::Logger &logger,
                    const std::experimental::filesystem::path &out_reads,
//                    const std::experimental::filesystem::path &bad_reads,
                    const std::experimental::filesystem::path &multiplicity_figures,
                    const io::Library &reads_lib, size_t unique_threshold,
                    size_t threads, const size_t min_read_size, bool dump) {
    size_t k = sdbg.hasher().k;
    ensure_dir_existance(multiplicity_figures);
    logger.info() << "Collecting info from reads" << std::endl;
//    size_t extension_size = std::max(std::min(min_read_size * 3 / 4, sdbg.hasher().k * 11 / 2), sdbg.hasher().k * 3 / 2);
    size_t min_extension = sdbg.hasher().k * 2 / 3;
    RecordStorage reads_storage(sdbg, min_extension, 100000, true);
    io::SeqReader readReader(reads_lib);
    reads_storage.fill(readReader.begin(), readReader.end(), min_read_size, logger, threads);
    size_t clow = 0;
    size_t cerr = 0;
    size_t both = 0;
    {
        UniqueClassificator classificator(sdbg);
        classificator.classify(logger, unique_threshold, multiplicity_figures);
    }
    logger.info() << "Printing reads to disk" << std::endl;
    std::ofstream ors;
    std::ofstream brs;
    ors.open(out_reads);
//    brs.open(bad_reads);
    for(auto it = reads_storage.begin(); it != reads_storage.end(); ++it) {
        AlignedRead &alignedRead = *it;
        ors << ">" << alignedRead.id << "\n" << alignedRead.path.getAlignment().Seq() << "\n";
//        for(auto edge_it : alignedRead.path.getAlignment().path()) {
//            if(edge_it->getCoverage() < threshold) {
//                brs << ">" << alignedRead.id << "\n" << alignedRead.path.getAlignment().Seq() << "\n";
//                break;
//            }
//        }
    }
    ors.close();
//    brs.close();
}