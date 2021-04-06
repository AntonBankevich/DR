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
    {
        UniqueClassificator classificator(sdbg);
        classificator.classify(logger, unique_threshold, multiplicity_figures);
        std::unordered_map<Edge *, CompactPath> unique_extensions;
        for(Edge *edge : classificator.unique_set) {
            const VertexRecord &rec = reads_storage.getRecord(*edge->start());
            Sequence seq= edge->seq.Subseq(0, 1);
            while(true) {
                Sequence best;
                size_t best_val = 0;
                for(char c = 0; c < char(4); c++) {
                    Sequence next = seq + Sequence(std::vector<char>({c}));
                    size_t val = rec.countStartsWith(next);
                    if(val > best_val) {
                        best_val = val;
                        best = next;
                    }
                }
                if(best_val == 0)
                    break;
                seq = best;
            }
            unique_extensions.emplace(edge, CompactPath(*edge->end(), seq.Subseq(1), 0, 0));
        }
        for(AlignedRead &alignedRead : reads_storage) {
            GraphAlignment al = alignedRead.path.getAlignment();
            bool corrected = false;
            bool bad = false;
            for(size_t i = 0; i < al.size(); i++) {
                if(unique_extensions.find(&al[i].contig()) == unique_extensions.end())
                    continue;
                CompactPath &compactPath = unique_extensions.find(&al[i].contig())->second;
                if(compactPath.seq().nonContradicts(alignedRead.path.seq()))
                    continue;
                corrected = true;
                GraphAlignment new_al = al.subPath(0, i + 1);
                size_t corrected_len = al.subPath(i + 1, al.size()).len();
                GraphAlignment replacement = compactPath.getAlignment();
                if(replacement.len() < corrected_len) {
                    logger.info() << "Need to correct more than known " << alignedRead.id << "\n"
                            << CompactPath(al.subPath(i + 1, al.size())) << "\n" << compactPath << std::endl;
                    bad = true;
                    break;
                }
                size_t pos = 0;
                while(corrected_len > 0) {
                    if(corrected_len < replacement[pos].size()) {
                        new_al += replacement[pos].shrinkRight(replacement[pos].size() - corrected_len);
                        corrected_len = 0;
                    } else {
                        new_al += replacement[pos];
                        corrected_len -= replacement[pos].size();
                    }
                }
                al = new_al;
            }
            if(corrected) {
                reads_storage.reroute(alignedRead, alignedRead.path.getAlignment(), al);
                logger << "Corrected read " << alignedRead.id << " " << alignedRead.path << std::endl;
            }
        }
    }
    logger.info() << "Printing corrected reads to disk" << std::endl;
    std::ofstream ors;
    std::ofstream brs;
    ors.open(out_reads);
    for(auto it = reads_storage.begin(); it != reads_storage.end(); ++it) {
        AlignedRead &alignedRead = *it;
        ors << ">" << alignedRead.id << "\n" << alignedRead.path.getAlignment().Seq() << "\n";
    }
    ors.close();
}