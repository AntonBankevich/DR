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
                                                            bool remove_bad, bool skip, bool dump) {
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
        RecordStorage readStorage(dbg, 0, extension_size, true);
        RecordStorage refStorage(dbg, 0, extension_size, false);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), w + k - 1, logger, threads);
        initialCorrect(dbg, logger, dir / "correction.txt", dir / "corrected.fasta", dir / "good.fasta",
                       dir / "bad.fasta", dir / "new_reliable.fasta", readStorage, refStorage,
                       threshold, 2 * threshold, reliable_coverage, threads, dump);
        Component comp(dbg);
        DrawSplit(comp, dir / "split");
    };
    if(!skip)
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
                     double threshold, bool skip) {
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
    if(!skip)
        runInFork(cc_task);
    std::experimental::filesystem::path res = dir / "corrected.fasta";
    logger.info() << "Crude correction results with k = " << k << " printed to " << res << std::endl;
    return res;
}

std::experimental::filesystem::path MultCorrection(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                     const io::Library &reads_lib, size_t threads, size_t k, size_t w, size_t unique_threshold, bool skip, bool dump) {
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
        RecordStorage readStorage(dbg, 0, 100000, true);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), w + k - 1, logger, threads);
        MultCorrect(dbg, logger, dir, readStorage, unique_threshold, threads, dump);
    };
    if(!skip)
        runInFork(mc_task);
    std::experimental::filesystem::path res = dir / "corrected.fasta";
    logger.info() << "Multiplicity-based correction results with k = " << k << " printed to " << res << std::endl;
    return res;
}

std::experimental::filesystem::path SplitDataset(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                                                   const io::Library &reads_lib, size_t threads, size_t k, size_t w, bool skip, bool dump) {
    logger.info() << "Performing dataset splitting with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    RollingHash hasher(k, 239);
    std::experimental::filesystem::path res = dir/"split";
    recreate_dir(res);
    std::function<void()> split_task = [&dir, &logger, &hasher, k, w, &reads_lib, threads, &res, dump] {
        SparseDBG dbg = DBGPipeline(logger, hasher, w, reads_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        std::vector<Component> comps = LengthSplitter(4000000000ul).split(dbg);
        std::vector<std::ofstream *> os;
        std::vector<std::experimental::filesystem::path> result;
        for(size_t i = 0; i < comps.size(); i++) {
            os.emplace_back(new std::ofstream());
            result.emplace_back(res / std::to_string(i + 1));
            ensure_dir_existance(result.back());
            os.back()->open(result.back() / "corrected.fasta");
        }
        io::SeqReader read_reader(reads_lib);
        for(StringContig scontig : read_reader) {
            string initial_seq = scontig.seq;
            Contig contig = scontig.makeContig();
            if(contig.size() < hasher.getK() + w - 1)
                return;
            GraphAlignment al = dbg.align(contig.seq);
            for(size_t j = 0; j < comps.size(); j++) {
                for(size_t i = 0; i <= al.size(); i++) {
                    if(comps[j].v.find(al.getVertex(i).hash()) != comps[j].v.end()) {
                        *os[j] << ">" << contig.id << "\n" << initial_seq << "\n";
                        break;
                    }
                }
            }
        };
        for(size_t j = 0; j < comps.size(); j++) {
            os[j]->close();
            delete os[j];
            os[j] = nullptr;
        }
    };
    if(!skip)
        runInFork(split_task);
    logger.info() << "Splitted datasets with k = " << k << " were printed to " << dir << std::endl;
    return res;
}

void ConstructSubdataset(logging::Logger &logger, const std::experimental::filesystem::path &subdir,
                                                 size_t threads, size_t k, size_t w, bool skip, bool dump) {
    logger.info() << "Constructing subdataset graphs with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    RollingHash hasher(k, 239);
    std::function<void()> split_task = [&subdir, &logger, &hasher, k, w, threads, dump] {
        logger.info() << "Constructing subdataset graphs" << std::endl;
        io::Library sublib = {subdir/"corrected.fasta"};
        SparseDBG subdbg = DBGPipeline(logger, hasher, w, sublib, subdir, threads);
        subdbg.fillAnchors(w, logger, threads);
        CalculateCoverage(subdir, hasher, w, sublib, threads, logger, subdbg);
        RecordStorage reads_storage(subdbg, 0, 100000, true);
        io::SeqReader readReader(sublib);
        reads_storage.fill(readReader.begin(), readReader.end(), hasher.getK() + w - 1, logger, threads);
        reads_storage.printAlignments(logger, subdir/"alignments.txt");
        std::ofstream edges;
        edges.open(subdir / "graph.fasta");
        subdbg.printFasta(edges);
        edges.close();
        std::ofstream gfa;
        gfa.open(subdir / "graph.gfa");
        subdbg.printGFA(gfa, true);
        gfa.close();
        std::ofstream dot;
        dot.open(subdir / "graph.dot");
        subdbg.printDot(dot, true);
        dot.close();
    };
    if (!skip)
        runInFork(split_task);
    logger.info() << "Splitted dataset with k = " << k << " were printed to " << subdir << std::endl;
}


int main(int argc, char **argv) {
    CLParser parser({"output-dir=", "threads=16", "k-mer-size=511", "window=2000", "K-mer-size=5001", "Window=500",
                     "cov-threshold=2", "rel-threshold=7", "Cov-threshold=2", "Rel-threshold=7", "crude-threshold=3",
                     "unique-threshold=50000", "dump", "dimer-compress=1000000", "restart-from=none"},
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
    StringContig::homopolymer_compressing = true;
    StringContig::dimer_compressing = std::stoull(parser.getValue("dimer-compress"));
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile());
    for(size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    std::string first_stage = parser.getValue("restart-from");
    bool skip = first_stage != "none";
    logger.info() << "LJA pipeline started" << std::endl;

    size_t threads = std::stoi(parser.getValue("threads"));
    bool dump = parser.getCheck("dump");

    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));

    size_t k = std::stoi(parser.getValue("k-mer-size"));
    size_t w = std::stoi(parser.getValue("window"));
    double threshold = std::stod(parser.getValue("cov-threshold"));
    double reliable_coverage = std::stod(parser.getValue("rel-threshold"));
    if(first_stage == "initial1")
        skip = false;
    std::experimental::filesystem::path corrected1 =
            InitialCorrection(logger, dir / "initial1", lib, threads, k, w,
                              threshold, reliable_coverage, false, skip, dump);

    size_t K = std::stoi(parser.getValue("K-mer-size"));
    size_t W = std::stoi(parser.getValue("Window"));


    double Threshold = std::stod(parser.getValue("Cov-threshold"));
    double Reliable_coverage = std::stod(parser.getValue("Rel-threshold"));
    if(first_stage == "initial2")
        skip = false;
    std::experimental::filesystem::path corrected2 =
            InitialCorrection(logger, dir / "initial2", {corrected1}, threads, K, W,
                              Threshold, Reliable_coverage, false, skip, dump);

    double crude_threshold = std::stod(parser.getValue("crude-threshold"));
    if(first_stage == "crude")
        skip = false;
    std::experimental::filesystem::path corrected3 =
            CrudeCorrection(logger, dir / "crude", {corrected2}, threads, K, W, crude_threshold, skip);

    size_t unique_threshold = std::stoi(parser.getValue("unique-threshold"));
    if(first_stage == "mult")
        skip = false;
    std::experimental::filesystem::path corrected4 =
            MultCorrection(logger, dir / "mult", {corrected3}, threads, K, W, unique_threshold, skip, dump);

    if(first_stage == "split")
        skip = false;
    std::experimental::filesystem::path split_dir =
            SplitDataset(logger, dir / "subdatasets", {corrected4}, threads, K, W, skip, dump);

    if(first_stage == "sub")
        skip = false;
    for(auto & subdir : std::experimental::filesystem::directory_iterator(split_dir)) {
        ConstructSubdataset(logger, subdir, threads, K, W, skip, dump);
    }

    logger.info() << "Final corrected reads can be hound here: " << corrected4 << std::endl;
    logger.info() << "Subdatasets for connected components can be found here: " << split_dir << std::endl;
    logger.info() << "LJA pipeline finished" << std::endl;
    return 0;
}
