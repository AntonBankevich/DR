#pragma once
#include "compact_path.hpp"


template <typename htype>
std::vector<Path<htype>> FindAlternatives(Vertex<htype> &start, Edge<htype> &edge, size_t max_diff) {
    size_t k = start.seq.size();
    std::vector<GraphAlignment<htype>> als = GraphAlignment<htype>(start).allExtensions(max_diff);
    max_diff = std::min(max_diff, edge.size());
    std::vector<Path<htype>> res;
    for(GraphAlignment<htype> &diff_al : als) {
        if(diff_al.size() > 0 && diff_al.front().contig() == edge) {
            continue;
        }
        for(size_t i = 0; i < max_diff; i++) {
            GraphAlignment<htype> al = diff_al;
            if(i > 0 && al.size() > 0 && al.lastNucl() == edge.seq[i - 1])
                continue;
            Sequence seq = edge.seq.Subseq(i, edge.size());
            al.extend(seq);
            if(al.valid() && al.endClosed() && al.back().contig().end() == edge.end()){
                res.emplace_back(al.path());
            }
        }
    }
    return res;
}

template <typename htype>
std::vector<Path<htype>> FilterAlternatives(const Path<htype> &initial, std::vector<Path<htype>> &als, double threshold) {
    size_t len = initial.len();
    std::vector<Path<htype>> res;
    size_t k = initial.getVertex(0).seq.size();
    size_t max_diff = std::max<size_t>(10, len / 200);
    for(Path<htype> &al : als) {
        bool ok = true;
        for(size_t i = 0; i < al.size(); i++) {
            if(al[i].getCoverage() < threshold) {
                ok = false;
                break;
            }
        }
        if(!ok) {
            continue;
        }
        size_t al_len = al.len();
        if(len > al_len + max_diff || al_len > len + max_diff)
            continue;
        res.emplace_back(al);
    }
    return res;
}

template<typename htype>
void initialCorrect(SparseDBG<htype> &sdbg, logging::Logger &logger,
                    const std::experimental::filesystem::path &out_file,
                    const io::Library &reads_lib,
                    const std::experimental::filesystem::path &ref,
                     size_t threads, const size_t min_read_size) {
    logger.info() << "Collecting info from reads" << std::endl;
    size_t extension_size = std::max(std::min(min_read_size * 2 / 3, sdbg.hasher().k * 3), sdbg.hasher().k * 3 / 2);
    size_t min_extension = sdbg.hasher().k * 2 / 3;
    RecordStorage<htype> reads_storage(sdbg, min_extension, extension_size);
    io::SeqReader readReader(reads_lib);
    reads_storage.fill(readReader.begin(), readReader.end(), min_read_size, logger, threads);
    logger.info() << "Collecting info from reference" << std::endl;
    RecordStorage<htype> ref_storage(sdbg, min_extension, extension_size);
    io::SeqReader refReader(ref);
    ref_storage.fill(refReader.begin(), refReader.end(), min_read_size, logger, threads);
    readReader.reset();
    std::ofstream out;
    out.open(out_file);
    for(auto it = reads_storage.begin(); it != reads_storage.end(); ++it) {
        AlignedRead<htype> &alignedRead = *it;
        CompactPath<htype> &initial_cpath = alignedRead.path;
        GraphAlignment<htype> path = initial_cpath.getAlignment();
        for(size_t i = 0; i < path.size(); i++) {
            Edge<htype> &edge = path[i].contig();
            Vertex<htype> &start = path.getVertex(i);
            if (edge.getCoverage() >= 2 || edge.size() > 3 * start.seq.size() || edge.end()->outDeg() == 0 || start.inDeg() == 0)
                continue;
            out << edge.size() << " " << edge.getCoverage();
            logger << "New bulge " <<   alignedRead.id << " " << i << " " << path.size() << " " <<
                                        edge.size() << " " << edge.getCoverage() << std::endl;
            logger <<   start.hash() << start.isCanonical() << " " << start.outDeg() << " " <<
                        start.inDeg() << " " << edge.end()->outDeg() << std::endl;
            std::vector<Path<htype>> alts = FindAlternatives(start, edge, 10);
            std::vector<Path<htype>> filtered = FilterAlternatives(Path<htype>(start, {&edge}), alts, 2);
            logger << alts.size() << " " << filtered.size() << std::endl;
            out << " " << filtered.size();
            std::vector<Path<htype>> read_alternatives = reads_storage.getRecord(start).getAlternatives(*edge.end(), 5);
            std::vector<Path<htype>> read_alternatives_filtered = FilterAlternatives(Path<htype>(start, {&edge}), read_alternatives, 2);;
            out << " " << read_alternatives_filtered.size();
//            if(filtered.size() > 0) {
//                size_t cut = edge.size() - std::min<size_t>(edge.size() / 2, 20);
//                logger << edge.seq.Subseq(0, edge.size() - cut) << std::endl;
//                for(size_t j = 0; j < filtered.size(); j++) {
//                    Sequence tmp = filtered[j].truncSeq();
//                    logger << tmp.Subseq(0, tmp.size() - std::min(tmp.size() / 2, cut)) << std::endl;
//                    CompactPath<htype> cpath(filtered[j]);
//                    size_t cnt_reads = reads_storage.getRecord(start).countStartsWith(cpath.cpath());
//                    size_t cnt_ref = ref_storage.getRecord(start).countStartsWith(cpath.cpath());
//                    logger << "Support " << cpath << " " << cnt_reads << " " << cnt_ref << std::endl;
//                    if(j == 0)
//                        out << " " << cnt_reads << " " << cnt_ref;
//                }
//            } else {
//                out << " 0 0";
//            }
//            out << " " << reads_storage.getRecord(start).uniqueOut(2) << std::endl;
//            logger << "Read records" << std::endl;
//            logger << reads_storage.getRecord(start) << std::endl;
//            logger << "Ref records" << std::endl;
//            logger << ref_storage.getRecord(start) << std::endl;
            out << std::endl;
        }
    }
    out.close();
}