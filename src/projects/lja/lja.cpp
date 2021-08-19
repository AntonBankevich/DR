#include "subdataset_processing.hpp"
#include "gap_closing.hpp"
#include "error_correction/mult_correction.hpp"
#include "dbg/dbg_construction.hpp"
#include "common/rolling_hash.hpp"
#include "error_correction/crude_correct.hpp"
#include "error_correction/initial_correction.hpp"
#include "sequences/seqio.hpp"
#include "common/dir_utils.hpp"
#include "common/cl_parser.hpp"
#include "common/logging.hpp"
#include "error_correction/manyk_correction.hpp"
#include <iostream>
#include <queue>
#include <unordered_set>
#include <wait.h>
#include <error_correction/dimer_correction.hpp>
#include <polishing/homopolish.hpp>

static size_t stage_num = 0;
std::vector<Contig> ref;
void PrintPaths(logging::Logger &logger, const std::experimental::filesystem::path &dir, const std::string &stage,
                SparseDBG &dbg, RecordStorage &readStorage, const io::Library &paths_lib, bool small) {
    stage_num += 1;
    std::string stage_name = itos(stage_num) + "_" + stage;
    logger.info() << "Dumping current state. Stage id: " << stage_name << std::endl;
    ensure_dir_existance(dir);
    ensure_dir_existance(dir / "paths");
    printDot(dir / (stage_name + ".dot"), Component(dbg), readStorage.labeler());
    dbg.printFastaOld(dir / (stage_name + ".fasta"));
    if(!small)
        readStorage.printFullAlignments(logger, dir / (stage_name + ".als"));
    std::vector<Contig> paths = io::SeqReader(paths_lib).readAllContigs();
    GraphAlignmentStorage storage(dbg);
    for(Contig &contig : paths) {
        storage.fill(contig);
    }
    for(Contig &contig : paths) {
        ensure_dir_existance(dir / "paths" / contig.getId());
        Component comp = small ? Component::neighbourhood(dbg, contig, dbg.hasher().getK() + 500) :
                Component::longEdgeNeighbourhood(dbg, contig, 20000);
        std::function<std::string(Edge &)> labeler = readStorage.labeler() + storage.labeler();
        printDot(dir / "paths" / contig.getId() / (stage_name + ".dot"), comp, labeler);
    }
    std::ofstream ref_os;
    ref_os.open(dir / (stage_name + ".ref"));
    for(Contig &contig : ref){
        ref_os << contig.getId() << std::endl;
        for(const PerfectAlignment<Contig, Edge> &al : GraphAligner(dbg).carefulAlign(contig)) {
            ref_os << al.seg_from << " " << al.seg_to << std::endl;
        }
    }
    ref_os.close();
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
        PrintPaths(logger, dir/ "state_dump", "initial", dbg, readStorage, paths_lib, true);
        initialCorrect(dbg, logger, dir / "correction.txt", readStorage, refStorage,
                       threshold, 2 * threshold, reliable_coverage, threads, dump);
        PrintPaths(logger, dir/ "paths", "low", dbg, readStorage, paths_lib, true);
        if(close_gaps) {
            GapColserPipeline(logger, dbg, readStorage, refStorage, threads);
            PrintPaths(logger, dir/ "paths", "gap", dbg, readStorage, paths_lib, true);
        }
        if(remove_bad) {
            readStorage.invalidateBad(logger, threads, threshold, "after_gap");
            RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
            PrintPaths(logger, dir/ "paths", "bad", dbg, readStorage, paths_lib, true);
        }
        readStorage.printFasta(logger, dir / "corrected.fasta");
        DrawSplit(Component(dbg), dir / "split");
        dbg.printFastaOld(dir / "graph.fasta");
        readStorage.printFullAlignments(logger, dir / "fullals.txt");
    };
    if(!skip)
        runInFork(ic_task);
    std::experimental::filesystem::path res;
    res = dir / "corrected.fasta";
    logger.info() << "Initial correction results with k = " << k << " printed to " << res << std::endl;
    return {res, dir / "graph.fasta"};
}

std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path> AlternativeCorrection(logging::Logger &logger, const std::experimental::filesystem::path &dir,
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
        size_t extension_size = std::max<size_t>(k * 2, 1000);
        RecordStorage readStorage(dbg, 0, extension_size, threads, dir/"read_log.txt", true);
        RecordStorage refStorage(dbg, 0, extension_size, threads, "/dev/null", false);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        PrintPaths(logger, dir/ "state_dump", "initial", dbg, readStorage, paths_lib, true);

        correct_dimers(logger, readStorage, k, threads);
//        correctAT(logger, readStorage, k, threads);
        ManyKCorrect(logger, dbg, readStorage, threshold, reliable_coverage, 800, 4, threads);
        PrintPaths(logger, dir/ "state_dump", "mk800", dbg, readStorage, paths_lib, true);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, std::max<size_t>(k * 5 / 2, 3000));
        ManyKCorrect(logger, dbg, readStorage, threshold, reliable_coverage, 2000, 4, threads);
        PrintPaths(logger, dir/ "state_dump", "mk2000", dbg, readStorage, paths_lib, true);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, std::max<size_t>(k * 7 / 2, 5000));
        correct_dimers(logger, readStorage, k, threads);
//        correctAT(logger, readStorage, k, threads);
        correctLowCoveredRegions(logger, dbg, readStorage, refStorage, "/dev/null", threshold, reliable_coverage, k, threads, dump);
        ManyKCorrect(logger, dbg, readStorage, threshold, reliable_coverage, 3500, 4, threads);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
        PrintPaths(logger, dir/ "state_dump", "mk3500", dbg, readStorage, paths_lib, false);
        readStorage.printFasta(logger, dir / "corrected.fasta");
        DrawSplit(Component(dbg), dir / "split");
        dbg.printFastaOld(dir / "graph.fasta");
//        readStorage.printFullAlignments(logger, dir / "fullals.txt");
    };
    if(!skip)
        runInFork(ic_task);
    std::experimental::filesystem::path res;
    res = dir / "corrected.fasta";
    logger.info() << "Initial correction results with k = " << k << " printed to " << res << std::endl;
    return {res, dir / "graph.fasta"};
}

std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path> SecondPhase(
        logging::Logger &logger, const std::experimental::filesystem::path &dir,
        const io::Library &reads_lib, const io::Library &pseudo_reads_lib,
        const io::Library &paths_lib, size_t threads, size_t k, size_t w,
        double threshold, double reliable_coverage, size_t unique_threshold,
        bool diploid, bool skip, bool dump, bool load) {
    logger.info() << "Performing second phase of correction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k, 239);
    std::function<void()> ic_task = [&dir, &logger, &hasher, load, k, w, &reads_lib, &pseudo_reads_lib, &paths_lib,
                                     threads, threshold, reliable_coverage, dump, unique_threshold, diploid] {
        io::Library construction_lib = reads_lib + pseudo_reads_lib;
        SparseDBG dbg = load ? DBGPipeline(logger, hasher, w, reads_lib, dir, threads, (dir/"disjointigs.fasta").string(), (dir/"vertices.save").string()) :
                        DBGPipeline(logger, hasher, w, reads_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        size_t extension_size = std::max<size_t>(k * 5 / 2, 3000);
        RecordStorage readStorage(dbg, 0, extension_size, threads, dir/"read_log.txt", true);
        RecordStorage refStorage(dbg, 0, extension_size, threads, "/dev/null", false);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        DrawSplit(Component(dbg), dir / "before_figs", readStorage.labeler(), 25000);
        PrintPaths(logger, dir/ "state_dump", "initial", dbg, readStorage, paths_lib, false);
        initialCorrect(dbg, logger, dir / "correction.txt", readStorage, refStorage,
                       threshold, 2 * threshold, reliable_coverage, threads, dump);
        PrintPaths(logger, dir/ "state_dump", "low", dbg, readStorage, paths_lib, false);
        GapColserPipeline(logger, dbg, readStorage, refStorage, threads);
        PrintPaths(logger, dir/ "state_dump", "gap1", dbg, readStorage, paths_lib, false);
        readStorage.invalidateBad(logger, threads, threshold, "after_gap1");
        PrintPaths(logger, dir/ "state_dump", "bad", dbg, readStorage, paths_lib, false);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
        PrintPaths(logger, dir/ "state_dump", "uncovered1", dbg, readStorage, paths_lib, false);
        MultCorrect(dbg, logger, dir, readStorage, unique_threshold, threads, diploid, dump);
        PrintPaths(logger, dir/ "state_dump", "mult", dbg, readStorage, paths_lib, false);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
        PrintPaths(logger, dir/ "state_dump", "uncovered2", dbg, readStorage, paths_lib, false);
        GapColserPipeline(logger, dbg, readStorage, refStorage, threads);
        PrintPaths(logger, dir/ "state_dump", "gap2", dbg, readStorage, paths_lib, false);
        DrawSplit(Component(dbg), dir / "split_figs", readStorage.labeler());
        RepeatResolver rr(dbg, readStorage, dir / "split");
        std::function<bool(const dbg::Edge &)> is_unique = [unique_threshold](const Edge &edge) {
            return edge.size() > unique_threshold;
        };
        dbg.printFastaOld(dir / "graph.fasta");
        printDot(dir / "graph.dot", Component(dbg), readStorage.labeler());
        std::vector<Contig> partial_contigs = rr.ResolveRepeats(logger, threads, is_unique);
        logger.info()<< "Printing partial repeat resolution results to " << (dir / "partial.fasta") << std::endl;
        PrintFasta(partial_contigs, dir / "partial.fasta");
        std::vector<Contig> contigs = rr.CollectResults(logger, threads, partial_contigs, dir / "merging.txt", is_unique);
        PrintAlignments(logger, threads, contigs, readStorage, k, dir / "uncompressing");
        readStorage.printFasta(logger, dir / "corrected.fasta");
    };
    if(!skip)
        runInFork(ic_task);
    std::experimental::filesystem::path res;
    res = dir / "corrected.fasta";
    logger.info() << "Second phase results with k = " << k << " printed to " << res << std::endl;
    return {res, dir / "uncompressing"};
}

std::experimental::filesystem::path PolishingPhase(
        logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &output_file,
        const std::experimental::filesystem::path &contigs_file,
        const std::experimental::filesystem::path &alignments,
        const io::Library &reads, size_t dicompress, bool skip) {
    logger.info() << "Performing polishing and homopolymer uncompression" << std::endl;
    std::function<void()> ic_task = [&logger, threads, &output_file, &contigs_file, &alignments, &reads, dicompress] {
        Polish(logger, threads, output_file, contigs_file, alignments, reads, dicompress);
    };
    if(!skip)
        runInFork(ic_task);
    return output_file;
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
                     const io::Library &reads_lib, const io::Library &pseudo_reads_lib, size_t threads, size_t k, size_t w, size_t unique_threshold, bool diploid,
                     bool skip, bool dump) {
    logger.info() << "Performing multiplicity-based correction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k, 239);
    std::function<void()> mc_task = [&dir, &logger, &hasher, k, w, &reads_lib, &pseudo_reads_lib,
                                     threads, unique_threshold, dump, diploid] {
        const io::Library construction_lib = reads_lib + pseudo_reads_lib;
        SparseDBG dbg = DBGPipeline(logger, hasher, w, construction_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        RecordStorage readStorage(dbg, 0, 1000000, threads, dir/"read_log.txt", true);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        MultCorrect(dbg, logger, dir, readStorage, unique_threshold, threads, diploid, dump);
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
//        CalculateCoverage(subdir, hasher, w, sublib, threads, logger, subdbg);
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
        printDot(subdir / "graph.dot", Component(subdbg), reads_storage.labeler());
    };
    if (!skip)
        runInFork(split_task);
    logger.info() << "Splitted dataset with k = " << k << " were printed to " << subdir << std::endl;
}


int main(int argc, char **argv) {
    CLParser parser({"output-dir=", "threads=16", "k-mer-size=511", "window=2000", "K-mer-size=5001", "Window=500",
                     "cov-threshold=2", "rel-threshold=7", "Cov-threshold=2", "Rel-threshold=7", "crude-threshold=3",
                     "unique-threshold=50000", "dump", "dimer-compress=1000000000,1000000000,1", "restart-from=none", "load",
                     "alternative", "diploid", "debug"},
                    {"reads", "paths", "ref"},
                    {"o=output-dir", "t=threads", "k=k-mer-size","w=window", "K=K-mer-size","W=Window"},
                    "Error message not implemented");
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters:" << std::endl;
        std::cout << parser.check() << "\n" << std::endl;
        std::cout << parser.message() << std::endl;
        return 1;
    }

    bool debug = parser.getCheck("debug");
    StringContig::homopolymer_compressing = true;
    StringContig::SetDimerParameters(parser.getValue("dimer-compress"));
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), debug ? logging::debug : logging::trace);
    for(size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    bool diplod = parser.getCheck("diploid");
    std::string first_stage = parser.getValue("restart-from");
    bool skip = first_stage != "none";
    bool load = parser.getCheck("load");
    logger.info() << "LJA pipeline started" << std::endl;

    size_t threads = std::stoi(parser.getValue("threads"));
    bool dump = parser.getCheck("dump");

    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::Library paths = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("paths"));
    io::Library ref_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("ref"));
    if(!io::CheckLibrary(lib + paths +ref_lib)) {
        exit(1);
    }
    ref = io::SeqReader(ref_lib).readAllContigs();
    size_t k = std::stoi(parser.getValue("k-mer-size"));
    size_t w = std::stoi(parser.getValue("window"));
    double threshold = std::stod(parser.getValue("cov-threshold"));
    double reliable_coverage = std::stod(parser.getValue("rel-threshold"));
    std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path> corrected1;
    if (first_stage == "alternative")
        skip = false;
    corrected1 = AlternativeCorrection(logger, dir / "alternative", lib, {}, paths, threads, k, w,
                                  threshold, reliable_coverage, false, false, skip, dump, load);
    if (first_stage == "alternative" || first_stage == "none")
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
                        threads, K, W, Threshold, Reliable_coverage, unique_threshold, diplod, skip, dump, load);
    if(first_stage == "phase2")
        load = false;
    if(first_stage == "polishing")
        skip = false;
    std::experimental::filesystem::path final_contigs =
            PolishingPhase(logger, threads, dir / "assembly.fasta",
                           corrected2.second / "contigs.fasta",
                           corrected2.second / "good_alignments.txt",
                            lib, StringContig::max_dimer_size / 2, skip);
    if(first_stage == "polishing")
        load = false;
    logger.info() << "Final corrected reads can be found here: " << corrected2.first << std::endl;
    logger.info() << "Final homopolymer compressed contigs can be found here: " << (corrected2.second / "contigs.fasta") << std::endl;
//    logger.info() << "Subdatasets for connected components can be found here: " << corrected2.second << std::endl;
    logger.info() << "Final assembly can be found here: " << final_contigs << std::endl;
    logger.info() << "LJA pipeline finished" << std::endl;
    return 0;
}