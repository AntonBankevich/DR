#pragma once

#include "multiplicity_estimation.hpp"
#include "sparse_dbg.hpp"
#include "compact_path.hpp"
#include <experimental/filesystem>

void correctRead(logging::Logger &logger, std::unordered_map<Edge *, CompactPath> &unique_extensions,
                 const AlignedRead &alignedRead, GraphAlignment &al, bool &corrected);

void printAl(logging::Logger &logger, std::unordered_map<Edge *, CompactPath> &unique_extensions,
             const GraphAlignment &al);

std::unordered_map<Edge *, CompactPath> constructUniqueExtensions(logging::Logger &logger,
                                                                  const RecordStorage &reads_storage, const UniqueClassificator &classificator) {
    std::unordered_map<Edge *, CompactPath> unique_extensions;
    for(Edge *edge : classificator.unique_set) {
        const VertexRecord &rec = reads_storage.getRecord(*edge->start());
        Sequence seq = edge->seq.Subseq(0, 1);
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
        logger << "Unique extension for edge " << edge->str() << " " << seq.Subseq(1) << std::endl;
    }
    return std::move(unique_extensions);
}

GraphAlignment correctRead(logging::Logger &logger, std::string &read_id,
                           std::unordered_map<Edge *, CompactPath> &unique_extensions,
                           const GraphAlignment &initial_al) {
    CompactPath initialCompactPath(initial_al);
    GraphAlignment al = initial_al;
    bool bad;
    bool corrected = false;
    for(size_t i = 0; i + 1 < al.size(); i++) {
        logger << al[i].contig().str() << std::endl;
        if(unique_extensions.find(&al[i].contig()) == unique_extensions.end())
            continue;
        CompactPath &compactPath = unique_extensions.find(&al[i].contig())->second;
        logger << compactPath << " " << CompactPath(al.subPath(i + 1, al.size())) << std::endl;
        if(compactPath.seq().nonContradicts(CompactPath(al.subPath(i + 1, al.size())).cpath()))
            continue;
        logger << "Corrected" << std::endl;
        corrected = true;
        GraphAlignment new_al = al.subPath(0, i + 1);
        size_t corrected_len = al.subPath(i + 1, al.size()).len();
        GraphAlignment replacement = compactPath.getAlignment();
        logger << replacement.len() << " " << corrected_len << std::endl;
        if(replacement.len() < corrected_len) {
            size_t deficite = corrected_len - replacement.len();
            logger.info() << "Need to correct more than known " << read_id << "\n"
                          << CompactPath(al.subPath(i + 1, al.size())) << "\n" << compactPath << std::endl;
            new_al += replacement;
            while(new_al.finish().outDeg() == 1 && deficite > 0) {
                size_t len = std::min(deficite, new_al.finish()[0].size());
                new_al += Segment<Edge>(new_al.finish()[0], 0, len);
                deficite -= len;
            }
            bad = true;
        } else {
            for (Segment<Edge> &rep_seg : replacement) {
                if (corrected_len <= rep_seg.size()) {
                    new_al += rep_seg.shrinkRight(rep_seg.size() - corrected_len);
                    corrected_len = 0;
                    break;
                } else {
                    new_al += rep_seg;
                    corrected_len -= rep_seg.size();
                }
            }
        }
        al = new_al;
    }
    if(corrected)
        return std::move(al);
    else
        return initial_al;
}

void correctReads(logging::Logger &logger, RecordStorage &reads_storage,
                  std::unordered_map<Edge *, CompactPath> &unique_extensions) {
    for(AlignedRead &alignedRead : reads_storage) {
        const GraphAlignment al = alignedRead.path.getAlignment();
        if(al.size() > 1) {
            logger << "Processing " << alignedRead.id << std::endl;
            printAl(logger, unique_extensions, al);
            GraphAlignment corrected1 = correctRead(logger, alignedRead.id, unique_extensions, al);
            printAl(logger, unique_extensions, corrected1);
            printAl(logger, unique_extensions, corrected1.RC());
            GraphAlignment corrected2 = correctRead(logger, alignedRead.id, unique_extensions, corrected1.RC()).RC();
            printAl(logger, unique_extensions, corrected2.RC());
            printAl(logger, unique_extensions, corrected2);
            printAl(logger, unique_extensions, al);
            if(al != corrected2) {
                reads_storage.reroute(alignedRead, al, corrected2);
                logger << "Corrected read " << alignedRead.id << " " << alignedRead.path << std::endl;
            } else {
                logger << "No corrections were made" << std::endl;
            }
        }
    }
}

void printAl(logging::Logger &logger, std::unordered_map<Edge *, CompactPath> &unique_extensions,
             const GraphAlignment &al) {
    for(auto &piece : al) {
        logger << piece.contig().str() << " ";
        if(unique_extensions.find(&piece.contig()) != unique_extensions.end()) {
            logger << "+ ";
        }
    }
    logger << std::endl;
}

void MultCorrect(dbg::SparseDBG &sdbg, logging::Logger &logger,
                 const std::experimental::filesystem::path &dir,
                 const io::Library &reads_lib, size_t unique_threshold,
                 size_t threads, const size_t min_read_size, bool dump) {
    const std::experimental::filesystem::path fig_before = dir / "before.dot";
    const std::experimental::filesystem::path fig_after = dir / "after.dot";
    const std::experimental::filesystem::path out_reads = dir / "corrected.fasta";
    const std::experimental::filesystem::path out_alignments = dir / "alignments.txt";
    const std::experimental::filesystem::path multiplicity_figures = dir / "figs";
    size_t k = sdbg.hasher().k;
    ensure_dir_existance(multiplicity_figures);
    logger.info() << "Collecting info from reads" << std::endl;
//    size_t extension_size = std::max(std::min(min_read_size * 3 / 4, sdbg.hasher().k * 11 / 2), sdbg.hasher().k * 3 / 2);
    size_t min_extension = sdbg.hasher().k * 2 / 3;
    RecordStorage reads_storage(sdbg, min_extension, 100000, true);
    io::SeqReader readReader(reads_lib);
    reads_storage.fill(readReader.begin(), readReader.end(), min_read_size, logger, threads);
    {
        UniqueClassificator classificator(sdbg, reads_storage);
        classificator.classify(logger, unique_threshold, multiplicity_figures);
        std::unordered_map<Edge *, CompactPath> unique_extensions =
                constructUniqueExtensions(logger, reads_storage, classificator);
        Component all(sdbg);
        std::function<std::string(Edge &)> labeler = [&reads_storage](Edge &edge){
            const VertexRecord &rec = reads_storage.getRecord(*edge.start());
            std::stringstream ss;
            for(const auto &ext : rec) {
                if(ext.first[0] == edge.seq[0])
                    ss << ext.first << "(" << ext.second << ")\n";
            }
            return ss.str();
        };
        std::function<std::string(Edge &)> colorer = [&classificator](Edge &edge){
            if(classificator.isUnique(edge))
                return "red";
            else
                return "black";
        };
        {
            std::ofstream os;
            os.open(fig_before);
            all.printDot(os, labeler, colorer);
            os.close();
        }
        correctReads(logger, reads_storage, unique_extensions);
        {
            std::ofstream os;
            os.open(fig_after);
            all.printDot(os, labeler, colorer);
            os.close();
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
    reads_storage.printAlignments(logger, out_alignments);
    ors.close();
}