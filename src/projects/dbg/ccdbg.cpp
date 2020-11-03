
#include "dbg_construction.hpp"
#include "crude_correct.hpp"
#include "rolling_hash.hpp"
#include "hash_utils.hpp"
#include <sequences/seqio.hpp>
#include <common/cl_parser.hpp>
#include <common/dir_utils.hpp>
#include <common/logging.hpp>
#include <iostream>

int main(int argc, char **argv) {
    CLParser parser({"output-dir=", "threads=16", "k-mer-size=5000", "window=2000", "base=239","cov-threshold=2"},
                    {"reads"},
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
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile());
    for(size_t i = 0; i < argc; i++) {
        logger.noTimeSpace() << argv[i] << " ";
    }
    logger.noTimeSpace() << std::endl;
    size_t k = std::stoi(parser.getValue("k-mer-size"));
    if (k % 2 == 0) {
        logger << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    RollingHash<htype128> hasher(k, std::stoi(parser.getValue("base")));
    const size_t w = std::stoi(parser.getValue("window"));
    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    size_t threads = std::stoi(parser.getValue("threads"));
    omp_set_num_threads(threads);

    const std::experimental::filesystem::path initial_dir = dir / "initial_graph";
    ensure_dir_existance(initial_dir);
    SparseDBG<htype128> dbg = DBGPipeline(logger, hasher, w, lib, initial_dir, threads);
    dbg.fillAnchors(w, logger, threads);
    CalculateCoverage(dir, hasher, w, lib, threads, logger, dbg);

    {
        logger << "Printing initial graph to fasta file " << (initial_dir / "graph.fasta") << std::endl;
        std::ofstream edges;
        edges.open(initial_dir / "graph.fasta");
        dbg.printFasta(edges);
        edges.close();
        logger << "Printing graph to dot file " << (initial_dir / "graph.dot") << std::endl;
        std::ofstream dot;
        dot.open(initial_dir / "graph.dot");
        dbg.printDot(dot, true);
        dot.close();
    }

    size_t threshold = std::__cxx11::stoull(parser.getValue("cov-threshold"));
    std::experimental::filesystem::path corrected_reads = CrudeCorrect(logger, dbg, dir, w, lib, threads, threshold);

    io::Library corrected_lib = {corrected_reads};
    logger << "Reconstructing dbg from corrected reads" << std::endl;
    SparseDBG<htype128> dbg_corrected = DBGPipeline(logger, hasher, w, corrected_lib, dir, threads);
    dbg_corrected.fillAnchors(w, logger, threads);
    CalculateCoverage(dir, hasher, w, corrected_lib, threads, logger, dbg_corrected);
    {
        logger << "Printing graph to fasta file " << (dir / "graph.fasta") << std::endl;
        std::ofstream edges;
        edges.open(dir / "graph.fasta");
        dbg_corrected.printFasta(edges);
        edges.close();
        logger << "Printing graph to dot file " << (dir / "graph.dot") << std::endl;
        std::ofstream dot;
        dot.open(dir / "graph.dot");
        dbg_corrected.printDot(dot, true);
        dot.close();
    }

    std::experimental::filesystem::path alignments_file = alignLib(logger, dbg_corrected, corrected_lib, hasher, w, dir, threads);

    logger << "DBG construction finished" << std::endl;
    logger << "Corrected DBG edges can be found in " << (dir/"graph.fasta") << std::endl;
    logger << "Corrected DBG in dot format can be found in " << (dir/"graph.dot") << std::endl;
    logger << "Corrected reads can be found in " << corrected_reads << std::endl;
    logger << "Corrected reads alignments can be found in " << alignments_file << std::endl;
    return 0;
}
