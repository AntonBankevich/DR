#pragma once

#include "diploidy_analysis.hpp"
#include "multiplicity_estimation.hpp"
#include "sparse_dbg.hpp"
#include "compact_path.hpp"
#include <experimental/filesystem>

void correctRead(logging::Logger &logger, std::unordered_map<const Edge *, CompactPath> &unique_extensions,
                 const AlignedRead &alignedRead, GraphAlignment &al, bool &corrected);

void printAl(logging::Logger &logger, std::unordered_map<const Edge *, CompactPath> &unique_extensions,
             const GraphAlignment &al);

std::unordered_map<const Edge *, CompactPath> constructUniqueExtensions(logging::Logger &logger, SparseDBG &dbg,
                                                                  const RecordStorage &reads_storage, const UniqueClassificator &classificator) {
    std::unordered_map<const Edge *, CompactPath> unique_extensions;
    for(Edge &edge : dbg.edges()) {
        if(!classificator.isUnique(edge) || unique_extensions.find(&edge) != unique_extensions.end())
            continue;
        const VertexRecord &rec = reads_storage.getRecord(*edge.start());
        Sequence seq = edge.seq.Subseq(0, 1);
        dbg::Path path(*edge.start());
        path += edge;
        while(true) {
            Sequence best;
            size_t best_val = 0;
            Edge *next_edge = nullptr;
            for(Edge &next_candidate : path.finish()) {
                if(classificator.isError(next_candidate))
                    continue;
                Sequence next = seq + next_candidate.seq.Subseq(0, 1);
                size_t val = rec.countStartsWith(next);
                if(val > best_val) {
                    best_val = val;
                    best = next;
                    next_edge = &next_candidate;
                }
            }
            if(best_val == 0)
                break;
            seq = best;
            path += *next_edge;
            if(classificator.isUnique(*next_edge))
                break;
        }
        CompactPath res(*edge.end(), seq.Subseq(1), 0, 0);
        unique_extensions.emplace(&edge, CompactPath(*edge.end(), seq.Subseq(1), 0, 0));
        logger << "Unique extension for edge " << edge.str() << " " << seq.Subseq(1) << std::endl;
        if(classificator.isUnique(path.back())) {
            Edge &last_rc_edge = path.back().rc();
            CompactPath res1(path.RC().subPath(1, path.size()));
            unique_extensions.emplace(&last_rc_edge, res1);
            logger << "Unique extension for edge " << last_rc_edge.str() << " " << res1.cpath() << std::endl;
        }
    }
    return std::move(unique_extensions);
}

GraphAlignment correctRead(logging::Logger &logger, std::string &read_id,
                           std::unordered_map<const Edge *, CompactPath> &unique_extensions,
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
        while(replacement.len() < corrected_len &&
                    unique_extensions.find(&replacement.back().contig()) != unique_extensions.end()) {
            replacement += unique_extensions[&replacement.back().contig()].getAlignment();
        }
        logger << CompactPath(replacement) << " " << CompactPath(al.subPath(i + 1, al.size())) << std::endl;
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
                  std::unordered_map<const Edge *, CompactPath> &unique_extensions) {
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

void printAl(logging::Logger &logger, std::unordered_map<const Edge *, CompactPath> &unique_extensions,
             const GraphAlignment &al) {
    for(auto &piece : al) {
        logger << piece.contig().str() << " ";
        if(unique_extensions.find(&piece.contig()) != unique_extensions.end()) {
            logger << "+ ";
        }
    }
    logger << std::endl;
}

void NewMultCorrect(dbg::SparseDBG &sdbg, logging::Logger &logger,
                 const std::experimental::filesystem::path &dir,
                 RecordStorage &reads_storage, size_t unique_threshold,
                 size_t threads, bool dump) {
    const std::experimental::filesystem::path fig_before = dir / "before.dot";
    const std::experimental::filesystem::path fig_after = dir / "after.dot";
    const std::experimental::filesystem::path out_reads = dir / "corrected.fasta";
    const std::experimental::filesystem::path out_alignments = dir / "alignments.txt";
    const std::experimental::filesystem::path multiplicity_figures = dir / "figs";

    ensure_dir_existance(multiplicity_figures);
    SetUniquenessStorage initial_unique = BulgePathAnalyser(sdbg, unique_threshold).uniqueEdges();
    MultiplicityBoundsEstimator estimator(sdbg, initial_unique);
    estimator.update(logger, 3, multiplicity_figures);
}


void MultCorrect(dbg::SparseDBG &sdbg, logging::Logger &logger,
                 const std::experimental::filesystem::path &dir,
                 RecordStorage &reads_storage, size_t unique_threshold,
                 size_t threads, bool dump) {
    const std::experimental::filesystem::path fig_before = dir / "before.dot";
    const std::experimental::filesystem::path fig_after = dir / "after.dot";
    const std::experimental::filesystem::path out_reads = dir / "corrected.fasta";
    const std::experimental::filesystem::path out_alignments = dir / "alignments.txt";
    const std::experimental::filesystem::path multiplicity_figures = dir / "figs";
    size_t k = sdbg.hasher().getK();
    ensure_dir_existance(multiplicity_figures);
    logger.info() << "Collecting info from reads" << std::endl;
    {
        UniqueClassificator classificator(sdbg, reads_storage);
        classificator.classify(logger, unique_threshold, multiplicity_figures);
        std::unordered_map<const Edge *, CompactPath> unique_extensions =
                constructUniqueExtensions(logger, sdbg, reads_storage, classificator);
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
        std::function<std::string(Edge &)> colorer = classificator.colorer();
        {
            std::ofstream os;
            os.open(fig_before);
            printDot(os, all, labeler, colorer);
            os.close();
        }
        correctReads(logger, reads_storage, unique_extensions);
        {
            std::ofstream os;
            os.open(fig_after);
            printDot(os, all, labeler, colorer);
            os.close();
        }
    }
    logger.info() << "Collecting bad edges" << std::endl;
    std::unordered_set<Edge *> bad_edges;
    for(Edge & edge : sdbg.edges()) {
        if(edge.size() > k + 5000)
            continue;
        if(reads_storage.getRecord(*edge.start()).isDisconnected(edge) ||
                reads_storage.getRecord(*edge.rc().start()).isDisconnected(edge.rc())) {
            bad_edges.emplace(&edge);
        }
    }
    logger.info() << "Printing corrected reads to disk" << std::endl;
    std::ofstream ors;
    std::ofstream brs;
    ors.open(out_reads);
    for(auto it = reads_storage.begin(); it != reads_storage.end(); ++it) {
        AlignedRead &alignedRead = *it;
        GraphAlignment al = alignedRead.path.getAlignment();
        bool ok = true;
        for(Segment<dbg::Edge> &seg : al) {
            if(seg.contig().getCoverage() < 2 || bad_edges.find(&seg.contig()) != bad_edges.end()) {
                ok = false;
                break;
            }
        }
        if(ok) {
            ors << ">" << alignedRead.id << "\n" << alignedRead.path.getAlignment().Seq() << "\n";
        } else {
            logger << "Could not correct read " << alignedRead.id << ". Removing it from dataset." << std::endl;
            alignedRead.path = CompactPath();
        }
    }
    reads_storage.printAlignments(logger, out_alignments);
    ors.close();
}