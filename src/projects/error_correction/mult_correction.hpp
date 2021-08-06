#pragma once

#include "lja/diploidy_analysis.hpp"
#include "multiplicity_estimation.hpp"
#include "dbg/sparse_dbg.hpp"
#include "dbg/compact_path.hpp"
#include "dbg/graph_alignment_storage.hpp"
#include <experimental/filesystem>

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

struct UEdge {
    UEdge(Edge *from, Edge *to, const CompactPath &cpath, size_t support) : from(from), to(to), cpath(cpath),
                                                                            support(support) {}

    Edge *from;
    Edge *to;
    CompactPath cpath;
    size_t support;

    UEdge RC() const {
        return {&to->rc(), &from->rc(), cpath.RC(), support};
    }
};
//std::unordered_map<const Edge *, CompactPath> constructUniqueExtensions(logging::Logger &logger, SparseDBG &dbg,
//                                                                         const RecordStorage &reads_storage, const UniqueClassificator &classificator) {
//    std::unordered_map<Edge *, std::vector<UEdge>> bg;
//    for(Edge &edge : dbg.edges()) {
//        if(!classificator.isUnique(edge))
//            continue;
//        Vertex &start = *edge.start();
//        std::vector<Sequence> extensions;
//        for(auto & c : reads_storage.getRecord(start)) {
//            GraphAlignment al = CompactPath(start, c.first).getAlignment();
//            if(al.front().contig() != edge)
//                continue;
//            for(size_t i = 0; i < al.size(); i++) {
//                Segment<Edge> &seg = al[i];
//                if(classificator.isUnique(seg.contig())) {
//                    al = al.subalignment(0, i + 1);
//                    break;
//                }
//            }
//            if(!classificator.isUnique(al.back().contig()))
//                continue;
//            extensions.emplace_back(CompactPath(al).cpath());
//        }
//        std::sort(extensions.begin(), extensions.end());
//        extensions.erase(std::unique(extensions.begin(), extensions.end()), extensions.end());
//        for(Sequence &extension : extensions) {
//            CompactPath new_path(start, extension);
//            bg[&edge].emplace_back(&edge, &new_path.getAlignment().back().contig(), new_path,
//                                   reads_storage.getRecord(start).countStartsWith(extension));
//        }
//    }
//    std::unordered_map<Edge *, UEdge> choice;
//}

std::unordered_map<const Edge *, CompactPath> constructUniqueExtensions(logging::Logger &logger, SparseDBG &dbg,
                                        const RecordStorage &reads_storage, const UniqueClassificator &classificator) {
    std::unordered_map<const Edge *, CompactPath> unique_extensions;
    for(Edge &edge : dbg.edges()) {
        if (!classificator.isUnique(edge) || unique_extensions.find(&edge) != unique_extensions.end())
            continue;
        Vertex & start = *edge.start();
        const VertexRecord &rec = reads_storage.getRecord(start);
        Sequence seq = edge.seq.Subseq(0, 1);
        CompactPath path = rec.getFullUniqueExtension(seq, 1, 0);
        GraphAlignment al = path.getAlignment();
        for(size_t i = 1; i < al.size(); i++) {
            Segment<Edge> &seg = al[i];
            if(classificator.isUnique(seg.contig())) {
                al = al.subalignment(0, i + 1);
                break;
            }
        }
        if(!classificator.isUnique(al.back().contig()) || al.size() == 1)
            continue;
        bool ok = true;
        for(const auto &it : rec) {
            if(it.second != 0 && !path.cpath().nonContradicts(it.first)) {
                ok = false;
                break;
            }
        }
        if(!ok)
            continue;
        unique_extensions.emplace(&edge, CompactPath(*edge.end(), CompactPath(al).cpath().Subseq(1), 0, 0));
        CompactPath res1(al.RC().subalignment(1, al.size()));
        unique_extensions.emplace(&al.back().contig().rc(), res1);
    }
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
//        logger << "Unique extension for edge " << edge.str() << " " << seq.Subseq(1) << std::endl;
        if(classificator.isUnique(path.back())) {
            Edge &last_rc_edge = path.back().rc();
            CompactPath res1(path.RC().subPath(1, path.size()));
            unique_extensions.emplace(&last_rc_edge, res1);
//            logger << "Unique extension for edge " << last_rc_edge.str() << " " << res1.cpath() << std::endl;
        }
    }
    return std::move(unique_extensions);
}

GraphAlignment correctRead(std::unordered_map<const Edge *, CompactPath> &unique_extensions,
                           const GraphAlignment &initial_al) {
    CompactPath initialCompactPath(initial_al);
    GraphAlignment al = initial_al;
    bool bad;
    bool corrected = false;
    for(size_t i = 0; i + 1 < al.size(); i++) {
        if(unique_extensions.find(&al[i].contig()) == unique_extensions.end())
            continue;
        CompactPath &compactPath = unique_extensions.find(&al[i].contig())->second;
        if(compactPath.cpath().nonContradicts(CompactPath(al.subalignment(i + 1, al.size())).cpath()))
            continue;
        corrected = true;
        GraphAlignment new_al = al.subalignment(0, i + 1);
        size_t corrected_len = al.subalignment(i + 1, al.size()).len();
        GraphAlignment replacement = compactPath.getAlignment();
        while(replacement.len() < corrected_len &&
                    unique_extensions.find(&replacement.back().contig()) != unique_extensions.end()) {
            replacement += unique_extensions[&replacement.back().contig()].getAlignment();
        }
        if(replacement.len() < corrected_len) {
            size_t deficite = corrected_len - replacement.len();
//            logger.info() << "Need to correct more than known " << read_id << "\n"
//                          << CompactPath(al.subalignment(i + 1, al.size())) << "\n" << compactPath << std::endl;
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

void correctReads(logging::Logger &logger, size_t threads, RecordStorage &reads_storage,
                  std::unordered_map<const Edge *, CompactPath> &unique_extensions) {
    omp_set_num_threads(threads);
    logger.info() << "Correcting reads using unique edge extensions" << std::endl;
#pragma omp parallel for default(none) shared(reads_storage, unique_extensions)
    for(size_t i = 0; i < reads_storage.size(); i++) {
        AlignedRead &alignedRead = reads_storage[i];
        if(!alignedRead.valid())
            continue;
        const GraphAlignment al = alignedRead.path.getAlignment();
        if(al.size() > 1) {
            GraphAlignment corrected1 = correctRead(unique_extensions, al);
            GraphAlignment corrected2 = correctRead(unique_extensions, corrected1.RC()).RC();
            if(al != corrected2) {
                reads_storage.reroute(alignedRead, al, corrected2, "mult correction");
            }
        }
    }
    reads_storage.applyCorrections(logger, threads);
}

void NewMultCorrect(dbg::SparseDBG &sdbg, logging::Logger &logger,
                 const std::experimental::filesystem::path &dir,
                 RecordStorage &reads_storage, size_t unique_threshold,
                 size_t threads, bool dump) {
    const std::experimental::filesystem::path fig_before = dir / "before.dot";
    const std::experimental::filesystem::path fig_after = dir / "after.dot";
    const std::experimental::filesystem::path out_reads = dir / "corrected.fasta";
    const std::experimental::filesystem::path out_alignments = dir / "alignments.txt";
    const std::experimental::filesystem::path multiplicity_figures = dir / "mult_figs";

    recreate_dir(multiplicity_figures);
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
//    const std::experimental::filesystem::path out_reads = dir / "corrected.fasta";
//    const std::experimental::filesystem::path out_alignments = dir / "alignments.txt";
    const std::experimental::filesystem::path full_alignments = dir / "full_alignments.txt";
    const std::experimental::filesystem::path multiplicity_figures = dir / "mult_figs";
    size_t k = sdbg.hasher().getK();
    ensure_dir_existance(multiplicity_figures);
    {
        UniqueClassificator classificator(sdbg, reads_storage);
        classificator.classify(logger, unique_threshold, multiplicity_figures);
        std::unordered_map<const Edge *, CompactPath> unique_extensions =
                constructUniqueExtensions(logger, sdbg, reads_storage, classificator);
        Component all(sdbg);
        std::function<std::string(Edge &)> colorer = classificator.colorer();
        {
            std::ofstream os;
            os.open(fig_before);
            printDot(os, all, reads_storage.labeler(), colorer);
            os.close();
        }
        correctReads(logger, threads, reads_storage, unique_extensions);
        {
            std::ofstream os;
            os.open(fig_after);
            printDot(os, all, reads_storage.labeler(), colorer);
            os.close();
        }
    }
    logger.info() << "Collecting bad edges" << std::endl;
    std::unordered_set<Edge const *> bad_edges;
    for(Edge & edge : sdbg.edges()) {
        if(edge.size() > k + 5000)
            continue;
        if(reads_storage.getRecord(*edge.start()).isDisconnected(edge) ||
                reads_storage.getRecord(*edge.rc().start()).isDisconnected(edge.rc())) {
            bad_edges.emplace(&edge);
        }
    }
//    logger.info() << "Printing corrected reads to disk" << std::endl;
//    std::ofstream ors;
    std::ofstream brs;
//    ors.open(out_reads);
    std::function<bool(const Edge&)> is_bad = [&bad_edges](const Edge &edge) {
        return edge.getCoverage() < 2 || bad_edges.find(&edge) != bad_edges.end();
    };
    reads_storage.invalidateBad(logger, threads, is_bad, "after_mult");
//    for(auto & alignedRead : reads_storage) {
//        if(alignedRead.valid()) {
//            ors << ">" << alignedRead.id << "\n" << alignedRead.path.getAlignment().Seq() << "\n";
//        }
//    }
//    RemoveUncovered(logger, threads, sdbg, {&reads_storage});
//    reads_storage.printAlignments(logger, out_alignments);
//    ors.close();
}