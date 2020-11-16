#pragma once
#include "compact_path.hpp"


template <typename htype>
std::vector<Path<htype>> FindAlternatives(Path<htype> &path, size_t max_diff) {
    size_t k = path.start().seq.size();
    std::vector<GraphAlignment<htype>> als = GraphAlignment<htype>(path.start()).allExtensions(max_diff);
    max_diff = std::min(max_diff, path.len());
    std::vector<Path<htype>> res;
    Sequence path_seq = path.truncSeq();
    for(GraphAlignment<htype> &diff_al : als) {
        size_t path_pos = 0;
        size_t edge_pos = size_t (-1);
        for(size_t i = 0; i < max_diff; i++) {
            GraphAlignment<htype> al = diff_al;
            if(i > 0 && al.size() > 0 && al.lastNucl() == path[path_pos].seq[edge_pos])
                continue;
            Sequence seq = path_seq.Subseq(i, path_seq.size());
            al.extend(seq);
            if(al.valid() && al.endClosed() && al.back().contig().end() == &path.finish()){
                res.emplace_back(al.path());
            }
            edge_pos += 1;
            if(edge_pos == path[path_pos].size()) {
                path_pos += 1;
                edge_pos = 0;
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
                     double threshold, size_t threads, const size_t min_read_size) {
    size_t k = sdbg.hasher().k;
    logger.info() << "Collecting info from reads" << std::endl;
    size_t extension_size = std::max(std::min(min_read_size * 2 / 3, sdbg.hasher().k * 3), sdbg.hasher().k * 3 / 2);
    size_t min_extension = sdbg.hasher().k * 2 / 3;
    RecordStorage<htype> reads_storage(sdbg, min_extension, extension_size, true);
    io::SeqReader readReader(reads_lib);
    reads_storage.fill(readReader.begin(), readReader.end(), min_read_size, logger, threads);
    logger.info() << "Collecting info from reference" << std::endl;
    RecordStorage<htype> ref_storage(sdbg, min_extension, extension_size, false);
    io::SeqReader refReader(ref);
    ref_storage.fill(refReader.begin(), refReader.end(), min_read_size, logger, threads);
    readReader.reset();
    std::ofstream out;
    out.open(out_file);
    for(auto it = reads_storage.begin(); it != reads_storage.end(); ++it) {
        AlignedRead<htype> &alignedRead = *it;
        CompactPath<htype> &initial_cpath = alignedRead.path;
        GraphAlignment<htype> path = initial_cpath.getAlignment();
        std::vector<Segment<Edge<htype>>> corrected_path;
        bool corrected = false;
        for(size_t path_pos = 0; path_pos < path.size(); path_pos++) {
            Edge<htype> &edge = path[path_pos].contig();
            if (edge.getCoverage() >= threshold || edge.size() > 5 * k ||
                    edge.end()->outDeg() == 0 || path.getVertex(path_pos).inDeg() == 0) {
                corrected_path.emplace_back(path[path_pos]);
                continue;
            }
            size_t start_ind = path_pos;
            size_t end_ind = path_pos + 1;
            size_t size = edge.size();
            while(start_ind > 0 && path[start_ind - 1].contig().getCoverage() < 10) {
                start_ind -= 1;
                size += path[start_ind].size();
            }
            while(end_ind < path.size() && path[end_ind].contig().getCoverage() < 10) {
                size += path[end_ind].size();
                end_ind += 1;
            }
            Vertex<htype> &start = path.getVertex(start_ind);
            Vertex<htype> &end = path.getVertex(end_ind);
            if(start.inDeg() == 0 || end.outDeg() == 0 || size > 5 * k) {
                for(size_t i = path_pos; i < end_ind; i++) {
                    corrected_path.emplace_back(path[i]);
                }
                path_pos = end_ind - 1;
                continue;
            }
            out << size << " " << end_ind - start_ind << " " << edge.getCoverage();
            logger << "New bulge " <<   alignedRead.id << " " << start_ind << "-" << end_ind << " " << path.size() << " " <<
                                        size << " " << edge.getCoverage() << std::endl;
            logger <<   start.hash() << start.isCanonical() << " " << start.outDeg() << " " <<
                        start.inDeg() << " " << edge.end()->outDeg() << std::endl;
            Path<htype> bulgePath = path.path().subPath(start_ind, end_ind);
            std::vector<Path<htype>> alts = FindAlternatives(bulgePath, 10);
            std::vector<Path<htype>> filtered = FilterAlternatives(Path<htype>(start, {&edge}), alts, threshold);
            logger << alts.size() << " " << filtered.size() << std::endl;
            out << " " << filtered.size();
            std::vector<Path<htype>> read_alternatives = reads_storage.getRecord(start).getAlternatives(end, threshold);
            std::vector<Path<htype>> read_alternatives_filtered = FilterAlternatives(Path<htype>(start, {&edge}), read_alternatives, threshold);;
            out << " " << read_alternatives_filtered.size();
            out << std::endl;
            if(read_alternatives_filtered.size() == 1) {
                corrected_path.erase(corrected_path.end() - (path_pos - start_ind), corrected_path.end());
                for(size_t i = 0; i < read_alternatives_filtered[0].size(); i++) {
                    corrected_path.emplace_back(read_alternatives_filtered[0][i], 0, read_alternatives_filtered[0][i].size());
                }
                corrected = true;
            } else {
                for(size_t i = path_pos; i < end_ind; i++) {
                    corrected_path.emplace_back(path[i]);
                }
            }
            path_pos = end_ind - 1;
        }
        if(corrected) {
            GraphAlignment<htype> corrected_alignment(&initial_cpath.start(), std::move(corrected_path));
            reads_storage.reroute(alignedRead, path, corrected_alignment);
        }
    }
    out.close();
}