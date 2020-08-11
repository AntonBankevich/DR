//
// Created by anton on 17.07.2020.
//
#define _GLIBCXX_PARALLEL
#include "dbg_construction.hpp"
#include "dbg_disjointigs.hpp"
#include "minimizer_selection.hpp"
#include "sparse_dbg.hpp"
#include "rolling_hash.hpp"
#include "hash_utils.hpp"
#include "sequences/seqio.hpp"
#include "common/dir_utils.hpp"
#include "common/cl_parser.hpp"
#include "common/logging.hpp"
#include <iostream>
#include <queue>
#include <omp.h>
#include <unordered_set>
using logging::Logger;


typedef unsigned __int128 htype128;

template<typename htype>
std::vector<Sequence> constructDisjointigs(const CLParser &parser, const std::experimental::filesystem::path &dir,
                                           const RollingHash<htype> &hasher, const size_t w, const string &reads_file, size_t threads,
                                           logging::Logger & logger) {
    std::vector<Sequence> disjointigs;
    if (parser.getValue("disjointigs") == "none") {
        std::vector<htype> hash_list;
        if (parser.getValue("unique") != "none") {
            std::ifstream is;
            is.open(parser.getValue("unique"));
            readHashs(is, hash_list);
            is.close();
        } else {
            hash_list = constructMinimizers(logger, reads_file, threads, hasher, w);
            std::ofstream os;
            os.open(std::string(dir.c_str()) + "/unique.save");
            writeHashs(os, hash_list);
            os.close();
        }

        SparseDBG<htype> sdbg = constructSparseDBGFromReads(logger, reads_file, threads, hasher, std::move(hash_list), w);
        sdbg.printStats(logger);
        //    sdbg.print();

        tieTips(logger, sdbg, w, threads);
        sdbg.printStats(logger);
        //    sdbg.print();

        disjointigs = extractDisjointigs(logger, sdbg, threads);
        std::ofstream df;
        df.open(std::string(dir.c_str()) + "/disjointigs.fasta");
        for(size_t i = 0; i < disjointigs.size(); i++) {
            df << ">" << i << std::endl;
            df << disjointigs[i] << std::endl;
        }
        df.close();
    } else {
        io::SeqReader reader(parser.getValue("disjointigs"));
        size_t cnt = 0;
        while(!reader.eof()) {
            disjointigs.push_back(reader.read().seq);
            cnt++;
        }
    }
    return disjointigs;
}

int main(int argc, char **argv) {
    CLParser parser({"vertices=none", "reads=", "unique=none", "output-dir=", "threads=8",
                     "k-mer-size=7000", "window=3000", "base=239", "debug", "disjointigs=none"},
            {"o=output-dir", "t=threads", "k=k-mer-size","w=window"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    Logger logger;
    logger.addLogFile(ls.newLoggerFile());
    const RollingHash<htype128> hasher(std::stoi(parser.getValue("k-mer-size")), std::stoi(parser.getValue("base")));
    const size_t w = std::stoi(parser.getValue("window"));
    std::string reads_file = parser.getValue("reads");
    size_t threads = std::stoi(parser.getValue("threads"));
    omp_set_num_threads(threads);
    std::vector<Sequence> disjointigs = constructDisjointigs(parser, dir, hasher, w, reads_file, threads, logger);
    std::vector<htype128> vertices;
    if (parser.getValue("vertices") == "none") {
        vertices = findJunctions(logger, disjointigs, hasher, threads);
        std::ofstream os;
        os.open(std::string(dir.c_str()) + "/vertices.save");
        writeHashs(os, vertices);
        os.close();
    } else {
        logger << "Loading vertex hashs from file" << parser.getValue("vertices") << std::endl;
        std::ifstream is;
        is.open(parser.getValue("vertices"));
        readHashs(is, vertices);
        is.close();
    }
    SparseDBG<htype128> dbg = constructDBG(logger, vertices, disjointigs, hasher);
    logger << "Printing graph to file" << std::endl;
    std::ofstream edges;
    edges.open(std::string(dir.c_str()) + "/graph.fasta");
    dbg.printFasta(edges);
    edges.close();
    logger << "DBG construction finished" << std::endl;
    return 0;
}
