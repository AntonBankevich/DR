#include "gap_closing.hpp"
#include "mult_correction.hpp"
#include "parameter_estimator.hpp"
#include "visualization.hpp"
#include "tip_correction.hpp"
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
#include "graph_printing.hpp"
#include <iostream>
#include <queue>
#include <omp.h>
#include <unordered_set>
#include <wait.h>

static size_t stage_num = 0;
void PrintPaths(logging::Logger &logger, const std::experimental::filesystem::path &dir, const std::string &stage,
                SparseDBG &dbg, RecordStorage &readStorage, const io::Library &paths_lib, bool small) {
    stage_num += 1;
    std::string stage_name = logging::itos(stage_num) + "_" + stage;
    printDot(dir / (stage_name + ".dot"), Component(dbg));
    dbg.printFastaOld(dir / (stage_name + ".fasta"));
    readStorage.printFullAlignments(logger, dir / (stage_name + ".als"));
    ensure_dir_existance(dir);
    std::vector<Contig> paths = io::SeqReader(paths_lib).readAllContigs();
    GraphAlignmentStorage storage(dbg);
    for(Contig &contig : paths) {
        storage.fill(contig);
    }
    for(Contig &contig : paths) {
        ensure_dir_existance(dir / contig.getId());
        Component comp = small ? Component::neighbourhood(dbg, contig, dbg.hasher().getK() + 500) :
                Component::longEdgeNeighbourhood(dbg, contig, 20000);
        std::function<std::string(Edge &)> labeler = readStorage.labeler() + storage.labeler();
        printDot(dir / contig.getId() / (stage_name + ".dot"), comp,labeler);
    }
}

std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path> InitialCorrection(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                                                                                                      const io::Library &reads_lib, const io::Library &pseudo_reads_lib,
                                                                                                      const io::Library &paths_lib, size_t threads, size_t k, size_t w,
                                                                                                      double threshold, double reliable_coverage,
                                                                                                      bool close_gaps, bool remove_bad, bool skip, bool dump, bool load) {
    logger.info() << "Performing initial correction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k, 239);
    std::function<void()> ic_task = [&dir, &logger, &hasher, close_gaps, load, remove_bad, k, w, &reads_lib,
                                     &pseudo_reads_lib, &paths_lib, threads, threshold, reliable_coverage, dump] {
        io::Library construction_lib = reads_lib + pseudo_reads_lib;
        SparseDBG dbg = load ? DBGPipeline(logger, hasher, w, reads_lib, dir, threads, (dir/"disjointigs.fasta").string(), (dir/"vertices.save").string()) :
                DBGPipeline(logger, hasher, w, reads_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        size_t extension_size = std::max<size_t>(k * 5 / 2, 3000);
        RecordStorage readStorage(dbg, 0, extension_size, threads, dir/"read_log.txt", true);
        RecordStorage refStorage(dbg, 0, extension_size, threads, "/dev/null", false);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        PrintPaths(logger, dir/ "paths", "paths1", dbg, readStorage, paths_lib, true);
        initialCorrect(dbg, logger, dir / "correction.txt", readStorage, refStorage,
                       threshold, 2 * threshold, reliable_coverage, threads, dump);
        PrintPaths(logger, dir/ "paths", "paths2", dbg, readStorage, paths_lib, true);
        if(close_gaps) {
            GapColserPipeline(logger, dbg, readStorage, refStorage, threads);
            PrintPaths(logger, dir/ "paths", "paths3", dbg, readStorage, paths_lib, true);
        }
        if(remove_bad) {
            readStorage.invalidateBad(logger, threshold);
            RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
            PrintPaths(logger, dir/ "paths", "paths4", dbg, readStorage, paths_lib, false);
        }
        readStorage.printFasta(logger, dir / "corrected.fasta");
        DrawSplit(Component(dbg), dir / "split");
        dbg.printFastaOld(dir / "graph.fasta");
    };
    if(!skip)
        runInFork(ic_task);
    std::experimental::filesystem::path res;
    res = dir / "corrected.fasta";
    logger.info() << "Initial correction results with k = " << k << " printed to " << res << std::endl;
    return {res, dir / "graph.fasta"};
}

void SplitDataset(SparseDBG &dbg, RecordStorage &readStorage, const std::experimental::filesystem::path &res) {
    size_t k = dbg.hasher().getK();
    std::vector<Component> comps = LengthSplitter(4000000000ul).split(dbg);
    std::vector<std::ofstream *> os;
    std::vector<std::ofstream *> alignments;
    std::vector<std::experimental::filesystem::path> result;
    for(size_t i = 0; i < comps.size(); i++) {
        os.emplace_back(new std::ofstream());
        alignments.emplace_back(new std::ofstream());
        result.emplace_back(res / std::to_string(i + 1));
        ensure_dir_existance(result.back());
        os.back()->open(result.back() / "corrected.fasta");
        alignments.back()->open(result.back() / "alignments.txt");
        std::ofstream gos;
        gos.open(result.back() / "graph.fasta");
        printFasta(gos, comps[i]);
        gos.close();
        std::ofstream log;
        log.open(result.back() / "dbg.log");
        log << "-k " << k << std::endl;
        log.close();
        std::ofstream dot;
        dot.open(result.back() / "graph.dot");
        printDot(dot, comps[i], readStorage.labeler());
        dot.close();
    }
    for(AlignedRead &read : readStorage) {
        if(!read.valid())
            continue;
        GraphAlignment al = read.path.getAlignment();
        Sequence seq = al.Seq();
        std::stringstream ss;
        ss  << read.id << " " << al.start().hash() << int(al.start().isCanonical())
            << " " << read.path.cpath().str() << "\n";
        CompactPath rc(al.RC());
        ss  << "-" << read.id << " " << rc.start().hash() << int(rc.start().isCanonical())
            << " " << rc.cpath().str() << "\n";
        std::string alignment_record = ss.str();
        for(size_t j = 0; j < comps.size(); j++) {
            for(size_t i = 0; i <= al.size(); i++) {
                if(comps[j].v.find(al.getVertex(i).hash()) != comps[j].v.end()) {
                    *os[j] << ">" << read.id << "\n" << seq << "\n";
                    *alignments[j] << alignment_record;
                    break;
                }
            }
        }
    };
    for(size_t j = 0; j < comps.size(); j++) {
        os[j]->close();
        alignments[j]->close();
        delete os[j];
        delete alignments[j];
        os[j] = nullptr;
        alignments[j] = nullptr;
    }
}

std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path> SecondPhase(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                                                                                                      const io::Library &reads_lib, const io::Library &pseudo_reads_lib,
                                                                                                      const io::Library &paths_lib, size_t threads, size_t k, size_t w,
                                                                                                      double threshold, double reliable_coverage, size_t unique_threshold,
                                                                                                      bool skip, bool dump, bool load) {
    logger.info() << "Performing initial correction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k, 239);
    std::function<void()> ic_task = [&dir, &logger, &hasher, load, k, w, &reads_lib, &pseudo_reads_lib, &paths_lib,
                                     threads, threshold, reliable_coverage, dump, unique_threshold] {
        io::Library construction_lib = reads_lib + pseudo_reads_lib;
        SparseDBG dbg = load ? DBGPipeline(logger, hasher, w, reads_lib, dir, threads, (dir/"disjointigs.fasta").string(), (dir/"vertices.save").string()) :
                        DBGPipeline(logger, hasher, w, reads_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        size_t extension_size = std::max<size_t>(k * 5 / 2, 3000);
        RecordStorage readStorage(dbg, 0, extension_size, threads, dir/"read_log.txt", true);
        RecordStorage refStorage(dbg, 0, extension_size, threads, "/dev/null", false);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        PrintPaths(logger, dir/ "paths", "initial", dbg, readStorage, paths_lib, false);
        initialCorrect(dbg, logger, dir / "correction.txt", readStorage, refStorage,
                       threshold, 2 * threshold, reliable_coverage, threads, dump);
        PrintPaths(logger, dir/ "paths", "low", dbg, readStorage, paths_lib, false);
        GapColserPipeline(logger, dbg, readStorage, refStorage, threads);
        PrintPaths(logger, dir/ "paths", "gap1", dbg, readStorage, paths_lib, false);
        readStorage.invalidateBad(logger, threshold);
        PrintPaths(logger, dir/ "paths", "bad", dbg, readStorage, paths_lib, false);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
        PrintPaths(logger, dir/ "paths", "uncovered1", dbg, readStorage, paths_lib, false);
        MultCorrect(dbg, logger, dir, readStorage, unique_threshold, threads, dump);
        PrintPaths(logger, dir/ "paths", "mult", dbg, readStorage, paths_lib, false);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
        PrintPaths(logger, dir/ "paths", "uncovered2", dbg, readStorage, paths_lib, false);
        GapColserPipeline(logger, dbg, readStorage, refStorage, threads);
        PrintPaths(logger, dir/ "paths", "gap2", dbg, readStorage, paths_lib, false);
        DrawSplit(Component(dbg), dir / "figs");
        SplitDataset(dbg, readStorage, dir / "split");
        readStorage.printFasta(logger, dir / "corrected.fasta");
        dbg.printFastaOld(dir / "graph.fasta");
    };
    if(!skip)
        runInFork(ic_task);
    std::experimental::filesystem::path res;
    res = dir / "corrected.fasta";
    logger.info() << "Second phase results with k = " << k << " printed to " << res << std::endl;
    return {res, dir / "split"};
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
    hashing::RollingHash hasher(k, 239);
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

std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path>
        MultCorrection(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                     const io::Library &reads_lib, const io::Library &pseudo_reads_lib, size_t threads, size_t k, size_t w, size_t unique_threshold, bool skip, bool dump) {
    logger.info() << "Performing multiplicity-based correction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k, 239);
    std::function<void()> mc_task = [&dir, &logger, &hasher, k, w, &reads_lib, &pseudo_reads_lib, threads, unique_threshold, dump] {
        const io::Library construction_lib = reads_lib + pseudo_reads_lib;
        SparseDBG dbg = DBGPipeline(logger, hasher, w, construction_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        RecordStorage readStorage(dbg, 0, 1000000, threads, dir/"read_log.txt", true);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        MultCorrect(dbg, logger, dir, readStorage, unique_threshold, threads, dump);
        std::ofstream edges;
        dbg.printFastaOld(dir / "graph.fasta");
        readStorage.printAlignments(logger, dir/"alignments.txt");
        readStorage.printFasta(logger, dir/"corrected.fasta");
    };
    if(!skip)
        runInFork(mc_task);
    std::experimental::filesystem::path res = dir / "corrected.fasta";
    logger.info() << "Multiplicity-based correction results with k = " << k << " printed to " << res << std::endl;
    return {res, dir / "graph.fasta"};
}

std::experimental::filesystem::path SplitDatasetStage(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                                                   const io::Library &reads_lib, const io::Library &pseudo_reads_lib, size_t threads, size_t k, size_t w, bool skip, bool dump) {
    logger.info() << "Performing dataset splitting with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k, 239);
    std::experimental::filesystem::path res = dir/"split";
    recreate_dir(res);
    std::function<void()> split_task = [&dir, &logger, &hasher, k, w, &reads_lib, &pseudo_reads_lib, threads, &res, dump] {
        io::Library construction_lib = reads_lib + pseudo_reads_lib;
        SparseDBG dbg = DBGPipeline(logger, hasher, w, construction_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        RecordStorage readStorage(dbg, 0, 1000000, threads, dir/"read_log.txt", true);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        {
            std::ofstream gos;
            gos.open(dir / "graph.dot");
            printDot(gos, Component(dbg), readStorage.labeler());
            gos.close();
        }
        SplitDataset(dbg, readStorage, res);
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
    hashing::RollingHash hasher(k, 239);
    std::function<void()> split_task = [&subdir, &logger, &hasher, k, w, threads, dump] {
        logger.info() << "Constructing subdataset graphs" << std::endl;
        io::Library sublib = {subdir/"corrected.fasta"};
        SparseDBG subdbg = DBGPipeline(logger, hasher, w, sublib, subdir, threads);
        subdbg.fillAnchors(w, logger, threads);
        CalculateCoverage(subdir, hasher, w, sublib, threads, logger, subdbg);
        RecordStorage reads_storage(subdbg, 0, 100000, threads, subdir/"read_log.txt", true);
        io::SeqReader readReader(sublib);
        reads_storage.fill(readReader.begin(), readReader.end(), subdbg, hasher.getK() + w - 1, logger, threads);
        reads_storage.printAlignments(logger, subdir/"alignments.txt");
        std::ofstream edges;
        edges.open(subdir / "graph.fasta");
        printFasta(edges, Component(subdbg));
        edges.close();
        std::ofstream gfa;
        gfa.open(subdir / "graph.gfa");
        printGFA(gfa, Component(subdbg), true);
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
                     "unique-threshold=50000", "dump", "dimer-compress=1000000000,1000000000,1", "restart-from=none", "load"},
                    {"reads", "paths"},
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
    StringContig::SetDimerParameters(parser.getValue("dimer-compress"));
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    ensure_dir_existance(dir / "paths");
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile());
    for(size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    std::string first_stage = parser.getValue("restart-from");
    bool skip = first_stage != "none";
    bool load = parser.getCheck("load");
    logger.info() << "LJA pipeline started" << std::endl;

    size_t threads = std::stoi(parser.getValue("threads"));
    bool dump = parser.getCheck("dump");

    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::Library paths = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("paths"));

    size_t k = std::stoi(parser.getValue("k-mer-size"));
    size_t w = std::stoi(parser.getValue("window"));
    double threshold = std::stod(parser.getValue("cov-threshold"));
    double reliable_coverage = std::stod(parser.getValue("rel-threshold"));
    if(first_stage == "initial1")
        skip = false;
    std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path> corrected1 =
            InitialCorrection(logger, dir / "initial1", lib, {}, paths, threads, k, w,
                              threshold, reliable_coverage, false, false, skip, dump, load);
    if(first_stage == "initial1")
        load = false;

    size_t K = std::stoi(parser.getValue("K-mer-size"));
    size_t W = std::stoi(parser.getValue("Window"));

    double Threshold = std::stod(parser.getValue("Cov-threshold"));
    double Reliable_coverage = std::stod(parser.getValue("Rel-threshold"));
    size_t unique_threshold = std::stoi(parser.getValue("unique-threshold"));

    if(first_stage == "phase2")
        skip = false;
    std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path> corrected2 =
            SecondPhase(logger, dir / "phase2", {corrected1.first}, {corrected1.second}, paths,
                        threads, K, W, Threshold, Reliable_coverage, unique_threshold, skip, dump, load);
    if(first_stage == "phase2")
        load = false;

//    double Threshold = std::stod(parser.getValue("Cov-threshold"));
//    double Reliable_coverage = std::stod(parser.getValue("Rel-threshold"));
//    if(first_stage == "initial2")
//        skip = false;
//    std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path> corrected2 =
//            InitialCorrection(logger, dir / "initial2", {corrected1.first}, {corrected1.second}, threads, K, W,
//                              Threshold, Reliable_coverage, true, true, skip, dump, load);
//    if(first_stage == "initial2")
//        load = false;
//
//    size_t unique_threshold = std::stoi(parser.getValue("unique-threshold"));
//    if(first_stage == "mult")
//        skip = false;
//    std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path> corrected4 =
//            MultCorrection(logger, dir / "mult", {corrected2.first}, {corrected2.second}, threads, K, W, unique_threshold, skip, dump);
//    if(first_stage == "mult")
//        load = false;
//
//    if(first_stage == "split")
//        skip = false;
//    std::experimental::filesystem::path split_dir =
//            SplitDatasetStage(logger, dir / "subdatasets", {corrected4.first}, {corrected4.second}, threads, K, W, skip, dump);
//    if(first_stage == "split")
//        load = false;

    logger.info() << "Final corrected reads can be found here: " << corrected2.first << std::endl;
    logger.info() << "Subdatasets for connected components can be found here: " << corrected2.second << std::endl;
    logger.info() << "LJA pipeline finished" << std::endl;
    return 0;
}
