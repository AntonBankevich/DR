#include "ff.hpp"
#include "graph_algorithms.hpp"
#include "dbg_construction.hpp"
#include "crude_correct.hpp"
#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include "visualization.hpp"
#include <sequences/seqio.hpp>
#include <common/cl_parser.hpp>
#include <common/dir_utils.hpp>
#include <common/logging.hpp>
#include <iostream>

int main(int argc, char **argv) {
    CLParser parser({"output-dir=", "threads=16", "k-mer-size=5000", "window=2000", "base=239","cov-threshold=2", "compress"},
                    {"reads"},
                    {"o=output-dir", "t=threads", "k=k-mer-size","w=window"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    if(parser.getCheck("compress"))
        StringContig::homopolymer_compressing = true;
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile());
    for(size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    size_t k = std::stoi(parser.getValue("k-mer-size"));
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    const size_t w = std::stoi(parser.getValue("window"));
    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    size_t threads = std::stoi(parser.getValue("threads"));
    size_t threshold = std::__cxx11::stoull(parser.getValue("cov-threshold"));
    RollingHash hasher(k, std::stoi(parser.getValue("base")));
    std::function<void()> initial_dbg_task = [&dir, &logger, &hasher, w, &lib, threads, threshold] {
        const std::experimental::filesystem::path initial_dir = dir / "initial_graph";
        ensure_dir_existance(initial_dir);
        SparseDBG dbg = DBGPipeline(logger, hasher, w, lib, initial_dir, threads);
        dbg.fillAnchors(w, logger, threads);
        CalculateCoverage(dir, hasher, w, lib, threads, logger, dbg);

        {
            logger.info() << "Printing initial graph to fasta file " << (initial_dir / "graph.fasta") << std::endl;
            std::ofstream edges;
            edges.open(initial_dir / "graph.fasta");
            dbg.printFasta(edges);
            edges.close();
            logger.info() << "Printing graph to dot file " << (initial_dir / "graph.dot") << std::endl;
            std::ofstream dot;
            dot.open(initial_dir / "graph.dot");
            dbg.printDot(dot, true);
            dot.close();
            Component comp(dbg);
            DrawSplit(comp, initial_dir/ "split");
        }
        CrudeCorrect(logger, dbg, dir, w, lib, threads, threshold);
    };
    runInFork(initial_dbg_task);
    std::experimental::filesystem::path corrected_reads = dir / "corrected.fasta";
    io::Library corrected_lib = {corrected_reads};
    logger.info() << "Reconstructing dbg from corrected reads" << std::endl;
    SparseDBG dbg_corrected = DBGPipeline(logger, hasher, w, corrected_lib, dir, threads);
    dbg_corrected.fillAnchors(w, logger, threads);
    CalculateCoverage(dir, hasher, w, corrected_lib, threads, logger, dbg_corrected);
    {
        logger.info() << "Printing graph to fasta file " << (dir / "graph.fasta") << std::endl;
        std::ofstream edges;
        edges.open(dir / "graph.fasta");
        dbg_corrected.printFasta(edges);
        edges.close();
        logger.info() << "Printing graph to dot file " << (dir / "graph.dot") << std::endl;
        std::ofstream dot;
        dot.open(dir / "graph.dot");
        dbg_corrected.printDot(dot, true);
        dot.close();
    }

    std::experimental::filesystem::path alignments_file = alignLib(logger, dbg_corrected, corrected_lib, hasher, w, dir, threads);

    logger.info() << "DBG construction finished" << std::endl;
    logger.info() << "Corrected DBG edges can be found in " << (dir/"graph.fasta") << std::endl;
    logger.info() << "Corrected DBG in dot format can be found in " << (dir/"graph.dot") << std::endl;
    logger.info() << "Corrected reads can be found in " << corrected_reads << std::endl;
    logger.info() << "Corrected reads alignments can be found in " << alignments_file << std::endl;
    return 0;
}
