//
// Created by anton on 8/7/20.
//

#include "visualization.hpp"
#include "sparse_dbg.hpp"
#include "common/hash_utils.hpp"
#include "common/cl_parser.hpp"
#include "common/logging.hpp"
#include "common/omp_utils.hpp"
#include "common/zip_utils.hpp"
#include "sequences/contigs.hpp"
#include "sequences/seqio.hpp"
#include "graph_algorithms.hpp"
#include <omp.h>
#include <vector>
#include <string>
#include <iostream>

int main(int argc, char **argv) {
    CLParser parser({"dbg=", "contigs=none", "threads=8", "k-mer-size=", "output-dir=", "base=239"}, {"segment"},
                    {"k=k-mer-size", "t=threads", "o=output-dir"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    size_t threads = std::stoi(parser.getValue("threads"));
    omp_set_num_threads(threads);
    logging::Logger logger;
    logger.info() << "Reading genome" << std::endl;
    StringContig::homopolymer_compressing = true;
    std::experimental::filesystem::path dbg_file(parser.getValue("dbg"));
    std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    size_t k = std::stoull(parser.getValue("k-mer-size"));
    RollingHash hasher(k, std::stoi(parser.getValue("base")));
    io::Library contigs_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("contigs"));
    std::vector<std::string> contig_segments = parser.getListValue("segment");

    SparseDBG dbg = LoadDBGFromFasta({std::experimental::filesystem::path(dbg_file)}, hasher, logger, threads);
    dbg.fillAnchors(1000, logger, threads);

    std::vector<std::tuple<std::string, size_t, size_t, std::string>> seg_recs;
    for(const std::string& s : contig_segments) {
        std::vector<std::string> parsed = split(s, "[,]");
        seg_recs.emplace_back(parsed[0], std::stoull(parsed[1]), std::stoull(parsed[2]), s);
    }
    std::vector<Contig> segs;
    io::SeqReader reader(contigs_lib);
    for(StringContig scontig : reader) {
        Contig contig = scontig.makeContig();
        for(auto & seg_rec : seg_recs) {
            if(std::get<0>(seg_rec) == contig.id) {
                segs.emplace_back(contig.seq.Subseq(std::get<1>(seg_rec), std::get<2>(seg_rec)), std::get<3>(seg_rec));
            } else if (std::get<0>(seg_rec) == "-" + contig.id) {
                segs.emplace_back((!contig.seq).Subseq(std::get<1>(seg_rec), std::get<2>(seg_rec)), std::get<3>(seg_rec));
            }
        }
    }
    std::vector<Component> components;
    for(Contig &seg : segs) {
        const std::experimental::filesystem::path seg_file = dir / ("seg_" + seg.id + ".dot");
        logger.info() << "Printing segment " << seg.id << " to file dot file " << (seg_file) << std::endl;
        std::ofstream coordinates_dot;
        components.emplace_back(Component::neighbourhood(dbg, seg, dbg.hasher().getK() + 100));
    }

    return 0;
}