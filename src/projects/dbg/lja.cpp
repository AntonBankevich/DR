#include "mult_correction.hpp"
#include "parameter_estimator.hpp"
#include "visualization.hpp"
#include "error_correction.hpp"
#include "graph_algorithms.hpp"
#include "dbg_construction.hpp"
#include "dbg_disjointigs.hpp"
#include "minimizer_selection.hpp"
#include "sparse_dbg.hpp"
#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include "crude_correct.hpp"
#include "common/hash_utils.hpp"
#include "initial_correction.hpp"
#include "sequences/seqio.hpp"
#include "common/dir_utils.hpp"
#include "common/cl_parser.hpp"
#include "common/logging.hpp"
#include <iostream>
#include <queue>
#include <omp.h>
#include <unordered_set>
#include <wait.h>

std::experimental::filesystem::path InitialCorrection(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                                                            const io::Library &reads_lib, size_t threads, size_t k, size_t w,
                                                            double threshold, double reliable_coverage,
                                                            bool remove_bad, bool dump) {
    logger.info() << "Performing initial correction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    RollingHash hasher(k, 239);
    std::function<void()> ic_task = [&dir, &logger, &hasher, k, w, &reads_lib, threads, threshold, reliable_coverage, dump] {
        SparseDBG dbg = DBGPipeline(logger, hasher, w, reads_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        size_t extension_size = std::max<size_t>(k * 5 / 2, 3000);
        initialCorrect(dbg, logger, dir / "correction.txt", dir / "corrected.fasta", dir / "good.fasta",
                       dir / "bad.fasta", dir / "new_reliable.fasta", reads_lib, {},
                       threshold, 2 * threshold, reliable_coverage, threads,
                       w + k - 1, extension_size, dump);
        Component comp(dbg);
        DrawSplit(comp, dir / "split");
    };
    runInFork(ic_task);
    std::experimental::filesystem::path res;
    if(remove_bad) {
        res = dir / "good.fasta";
    } else {
        res = dir / "corrected.fasta";
    }
    logger.info() << "Initial correction results with k = " << k << " printed to " << res << std::endl;
    return res;
}

std::experimental::filesystem::path CrudeCorrection(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                     const io::Library &reads_lib, size_t threads, size_t k, size_t w,
                     double threshold) {
    logger.info() << "Performing crude correction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    RollingHash hasher(k, 239);
    std::function<void()> cc_task = [&dir, &logger, &hasher, w, &reads_lib, threads, threshold] {
        SparseDBG dbg = DBGPipeline(logger, hasher, w, reads_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        CalculateCoverage(dir, hasher, w, reads_lib, threads, logger, dbg);
        CrudeCorrect(logger, dbg, dir, w, reads_lib, threads, threshold);
    };
    runInFork(cc_task);
    std::experimental::filesystem::path res = dir / "corrected.fasta";
    logger.info() << "Crude correction results with k = " << k << " printed to " << res << std::endl;
    return res;
}

std::experimental::filesystem::path MultCorrection(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                     const io::Library &reads_lib, size_t threads, size_t k, size_t w, size_t unique_threshold, bool dump) {
    logger.info() << "Performing multiplicity-based correction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    RollingHash hasher(k, 239);
    std::function<void()> mc_task = [&dir, &logger, &hasher, k, w, &reads_lib, threads, unique_threshold, dump] {
        SparseDBG dbg = DBGPipeline(logger, hasher, w, reads_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        MultCorrect(dbg, logger, dir, reads_lib, unique_threshold, threads, w + k - 1, dump);
    };
    runInFork(mc_task);
    std::experimental::filesystem::path res = dir / "corrected.fasta";
    logger.info() << "Multiplicity-based correction results with k = " << k << " printed to " << res << std::endl;
    return res;
}


int main(int argc, char **argv) {
    CLParser parser({"output-dir=", "threads=16", "k-mer-size=511", "window=2000", "K-mer-size=5001", "Window=500",
                     "cov-threshold=2", "rel-threshold=7", "Cov-threshold=2", "Rel-threshold=7", "crude-threshold=3",
                     "unique-threshold=50000", "dump"},
                    {"reads"},
                    {"o=output-dir", "t=threads", "k=k-mer-size","w=window", "K=K-mer-size","W=Window"},
                    "Error message not implemented");
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters:" << std::endl;
        std::cout << parser.check() << "\n" << std::endl;
        std::cout << parser.message() << std::endl;
        return 1;
    }
    StringContig::needs_compressing = true;
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile());
    for(size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    logger.info() << "LJA pipeline started" << std::endl;

    size_t threads = std::stoi(parser.getValue("threads"));
    bool dump = parser.getCheck("dump");

    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));

    size_t k = std::stoi(parser.getValue("k-mer-size"));
    size_t w = std::stoi(parser.getValue("window"));
    double threshold = std::stod(parser.getValue("cov-threshold"));
    double reliable_coverage = std::stod(parser.getValue("rel-threshold"));
    std::experimental::filesystem::path corrected1 =
            InitialCorrection(logger, dir / "initial1", lib, threads, k, w,
                              threshold, reliable_coverage, false, dump);

    size_t K = std::stoi(parser.getValue("K-mer-size"));
    size_t W = std::stoi(parser.getValue("Window"));
    double Threshold = std::stod(parser.getValue("Cov-threshold"));
    double Reliable_coverage = std::stod(parser.getValue("Rel-threshold"));
    std::experimental::filesystem::path corrected2 =
            InitialCorrection(logger, dir / "initial2", {corrected1}, threads, K, W,
                              Threshold, Reliable_coverage, false, dump);

    double crude_threshold = std::stod(parser.getValue("crude-threshold"));
    std::experimental::filesystem::path corrected3 =
            CrudeCorrection(logger, dir / "crude", {corrected2}, threads, K, W, crude_threshold);

    size_t unique_threshold = std::stoi(parser.getValue("unique-threshold"));
    std::experimental::filesystem::path corrected4 =
            MultCorrection(logger, dir / "mult", {corrected3}, threads, K, W, unique_threshold, dump);

    logger.info() << "Final corrected reads can be hound here: " << corrected4 << std::endl;
    logger.info() << "LJA pipeline finished" << std::endl;
    return 0;
}
