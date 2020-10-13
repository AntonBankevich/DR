//
// Created by anton on 17.07.2020.
//
#define _GLIBCXX_PARALLEL
#include "error_correction.hpp"
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
#include "hash_utils.hpp"
#include <iostream>
#include <queue>
#include <omp.h>
#include <unordered_set>
#include <wait.h>

using logging::Logger;

//
//class htype128 {
//    unsigned __int_128 val;
//public:
//    template<class T>
//    htype128(cosnt T &_val) : val(_val) {
//    }
//
//    htype128 operator+(const htype128& other) const {
//        return htype128(val + other.val);
//    }
//
//    htype128 operator+(const htype128& other) const {
//        return htype128(val + other.val);
//    }
//
//    htype128 operator*(const htype128& other) const {
//        return htype128(val * other.val);
//    }
//
//    htype128 operator/(const htype128& other) const {
//        return htype128(val / other.val);
//    }
//
//    bool operator<(const htype128& other) const {
//        return hval < other.val;
//    }
//    bool operator>(const htype128& other) const {
//        return hval > other.val;
//    }
//    bool operator<=(const htype128& other) const {
//        return hval <= other.val;
//    }
//    bool operator>=(const htype128& other) const {
//        return hval >= other.val;
//    }
//    bool operator==(const htype128& other) const {
//        return hval == other.val;
//    }
//    bool operator!=(const htype128& other) const {
//        return hval != other.val;
//    }
//};
//
//std::ostream &operator<<(std::ostream &os, htype128 val) {
//    while(val > 0) {
//        os << size_t(val % 10);
//        val /= 10;
//    }
//    return os;
//}

template<typename htype>
std::vector<Sequence> constructDisjointigs(const RollingHash<htype> &hasher, size_t w, const io::Library &reads_file,
                                           const std::vector<htype128> & hash_list, size_t cov_threshold, size_t threads,
                                           logging::Logger & logger) {
    std::vector<Sequence> disjointigs;
    SparseDBG<htype> sdbg = constructSparseDBGFromReads(logger, reads_file, threads, hasher, hash_list, w);
    sdbg.printStats(logger);
    sdbg.checkSeqFilled(threads, logger);

    tieTips(logger, sdbg, w, threads);
    sdbg.checkSeqFilled(threads, logger);
    sdbg.printStats(logger);
//    std::ofstream os;
//    os.open("sdbg.fasta");
//    sdbg.printFasta(os);
//    os.close();

    disjointigs = extractDisjointigs(logger, sdbg, threads);
    return disjointigs;
}

struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const
    {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

void analyseGenome(SparseDBG<htype128> &dbg, const std::string &ref_file, Logger &logger) {
    logger << "Reading reference" << std::endl;
    std::vector<StringContig> ref = io::SeqReader(ref_file).readAll();
    logger << "Finished reading reference. Starting alignment" << std::endl;
    std::vector<std::pair<Edge<htype128> const *, size_t>> path;
    for(StringContig & contig : ref) {
        auto tmp = dbg.align(contig.makeCompressedSequence());
        logger << "Aligned chromosome " << contig.id << " . Path length " << tmp.size() << std::endl;
        path.insert(path.end(), tmp.begin(), tmp.end());
    }
    logger << "Reference path consists of " << path.size() << " edges" << std::endl;
//    std::vector<size_t> fill_distr(11);
//    std::vector<size_t> fill_distr_len(11);
//    for(const std::pair<const Edge<htype128> *, size_t> & val : path) {
//        size_t ind = val.second * 10 / val.first->size();
//        fill_distr[ind] += 1;
//        fill_distr_len[ind] += 1;
//    }
//    logger << "Edge filling distribution" <<std::endl;
//    logger.noTimeSpace() << fill_distr << std::endl << fill_distr_len << std::endl;
    std::unordered_set<Edge<htype128> const *> eset;
    std::vector<size_t> cov(50);
    std::vector<size_t> cov_len(50);
    std::vector<size_t> cov_bad(50);
    std::vector<size_t> cov_bad_len(50);
    std::vector<size_t> cov_good(50);
    std::vector<size_t> cov_good_len(50);
    for(size_t i = 0; i < path.size(); i++) {
        const Edge<htype128> &edge = *path[i - 1].first;
        eset.emplace(path[i].first);
        size_t cov_val = std::min(cov.size() - 1, size_t(edge.getCoverage()));
        cov_good[cov_val] += 2;
        cov_good_len[cov_val] += edge.size() * 2;
    }
    for(auto & pair : dbg) {
        Vertex<htype128> &vert = pair.second;
        for (Edge<htype128> &edge : vert.getOutgoing()) {
            size_t cov_val = std::min(cov.size() - 1, size_t(edge.getCoverage()));
            if (eset.find(&edge) == eset.end() && eset.find(&vert.rcEdge(edge)) == eset.end()) {
                cov_bad[cov_val] += 1;
                cov_bad_len[cov_val] += edge.size();
            }
            cov[cov_val] += 1;
            cov_len[cov_val] += edge.size();
        }
    }
    logger << "All coverages" << std::endl;
    logger.noTimeSpace() << cov << std::endl << cov_len << std::endl;
    logger << "Coverages of edges in genome path" << std::endl;
    logger.noTimeSpace() << cov_good << std::endl << cov_good_len << std::endl;
    logger << "Coverages of edges outside genome path" << std::endl;
    logger.noTimeSpace() << cov_bad << std::endl << cov_bad_len << std::endl;
}

void
CalculateCoverage(const std::experimental::filesystem::path &dir, const RollingHash<htype128> &hasher, const size_t w,
                  const io::Library &lib, size_t threads, Logger &logger, SparseDBG<htype128> &dbg) {
    logger << "Calculating edge coverage." << std::endl;
    io::SeqReader reader(lib);
    fillCoverage(dbg, logger, reader.begin(), reader.end(), threads, hasher, w + hasher.k - 1);
    std::ofstream os;
    os.open(dir / "coverages.save");
    os << dbg.size() << std::endl;
    for (std::pair<const htype128, Vertex<htype128>> &pair : dbg) {
        Vertex<htype128> &v = pair.second;
        os << v.hash() << " " << v.outDeg() << " " << v.inDeg() << std::endl;
        for (const Edge<htype128> &edge : v.getOutgoing()) {
            os << size_t(edge.seq[0]) << " " << edge.intCov() << std::endl;
        }
        for (const Edge<htype128> &edge : v.rc().getOutgoing()) {
            os << size_t(edge.seq[0]) << " " << edge.intCov() << std::endl;
        }
    }
    dbg.printCoverageStats(logger);
    os.close();
}

void LoadCoverage(const CLParser &parser, Logger &logger, SparseDBG<htype128> &dbg) {
    logger << "Loading edge coverages." << std::endl;
    std::ifstream is;
    is.open(parser.getValue("coverages"));
    size_t n;
    is >> n;
    for (size_t i = 0; i < n; i++) {
        htype128 vid;
        is >> vid;
        Vertex<htype128> *v = &dbg.getVertex(vid);
        size_t inDeg, outDeg;
        is >> outDeg >> inDeg;
        for (size_t j = 0; j < inDeg + outDeg; j++) {
            if (j == outDeg)
                v = &v->rc();
            size_t next;
            is >> next;
            Edge<htype128> &edge = v->getOutgoing(char(next));
            size_t cov;
            is >> cov;
            edge.incCov(cov);
        }
    }
    is.close();
    logger << "Finished loading edge coverages." << std::endl;
}

int main(int argc, char **argv) {
    CLParser parser({"vertices=none", "unique=none", "dbg=none", "coverages=none", "segments=none", "dbg=none", "output-dir=",
                     "threads=8", "k-mer-size=5000", "window=3000", "base=239", "debug", "disjointigs=none", "reference=none",
                     "correct", "align=none", "simplify"},
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
    Logger logger;
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
    io::Library lib_for_coverage = lib;
    if (parser.getValue("reference") != "none") {
        logger << "Added reference to graph construction. Careful, some edges may have coverage 0" << std::endl;
        lib.insert(lib.begin(), std::experimental::filesystem::path(parser.getValue("reference")));
    }
    size_t threads = std::stoi(parser.getValue("threads"));
    omp_set_num_threads(threads);
    std::vector<Sequence> disjointigs;
    if (parser.getValue("disjointigs") == "none") {
        pid_t p = fork();
        if (p < 0) {
            std::cout << "Fork failed" << std::endl;
            return 1;
        }
        if(p == 0) {
            std::vector<htype128> hash_list;
            if (parser.getValue("unique") == "none") {
                hash_list = constructMinimizers(logger, lib, threads, hasher, w);
                std::ofstream os;
                os.open(std::string(dir.c_str()) + "/unique.save");
                writeHashs(os, hash_list);
                os.close();
            } else {
                logger << "Loading minimizers from file " << parser.getValue("unique") << std::endl;
                std::ifstream is;
                is.open(parser.getValue("unique"));
                readHashs(is, hash_list);
                is.close();
            }
            disjointigs = constructDisjointigs(hasher, w, lib, hash_list, 1, threads, logger);
            hash_list.clear();
            std::ofstream df;
            df.open(dir / "disjointigs.fasta");
            for (size_t i = 0; i < disjointigs.size(); i++) {
                df << ">" << i << std::endl;
                df << disjointigs[i] << std::endl;
            }
            df.close();
            return 0;
        } else {
            int status = 0;
            std::cout << "Waiting" << std::endl;
            waitpid(p, &status, 0);
            if (WEXITSTATUS(status) || WIFSIGNALED(status)) {
                std::cout << "Child process crashed" << std::endl;
                return 1;
            }
            logger << "Loading disjointigs from file " << (dir / "disjointigs.fasta") << std::endl;
            io::SeqReader reader(dir / "disjointigs.fasta");
            while(!reader.eof()) {
                disjointigs.push_back(reader.read().makeCompressedSequence());
            }
        }
    } else {
        logger << "Loading disjointigs from file " << parser.getValue("disjointigs") << std::endl;
        io::SeqReader reader(parser.getValue("disjointigs"));
        while(!reader.eof()) {
            disjointigs.push_back(reader.read().makeCompressedSequence());
        }
    }
    std::vector<htype128> vertices;
    if (parser.getValue("dbg") == "none") {
        if (parser.getValue("vertices") == "none") {
            vertices = findJunctions(logger, disjointigs, hasher, threads);
            std::ofstream os;
            os.open(std::string(dir.c_str()) + "/vertices.save");
            writeHashs(os, vertices);
            os.close();
        } else {
            logger << "Loading vertex hashs from file " << parser.getValue("vertices") << std::endl;
            std::ifstream is;
            is.open(parser.getValue("vertices"));
            readHashs(is, vertices);
            is.close();
        }
    }

    SparseDBG<htype128> dbg = parser.getValue("dbg") == "none" ?
            constructDBG(logger, vertices, disjointigs, hasher, threads) :
          SparseDBG<htype128>::loadDBGFromFasta({std::experimental::filesystem::path(parser.getValue("dbg"))},
                  hasher, logger, threads);
    dbg.checkConsistency(threads, logger);

    if(parser.getValue("dbg") == "none") {
        logger << "Printing graph to file " << (dir / "graph.fasta") << std::endl;
        std::ofstream edges;
        edges.open(dir / "graph.fasta");
        dbg.printFasta(edges);
        edges.close();
    }

    if (parser.getValue("align") != "none" || parser.getCheck("correct") || parser.getValue("segments") != "none"
                || parser.getValue("reference") != "none") {
        dbg.fillAnchors(w, logger, threads);
    }

    if (parser.getCheck("simplify") || parser.getCheck("correct") ||
            parser.getValue("segments") != "none" || parser.getValue("reference") != "none") {
        if (parser.getValue("coverages") == "none") {
            CalculateCoverage(dir, hasher, w, lib_for_coverage, threads, logger, dbg);
        } else {
            LoadCoverage(parser, logger, dbg);
        }
    }

    if (parser.getValue("align") != "none") {
        std::ofstream os(dir / "alignments.txt");
        io::SeqReader reader(parser.getValue("align"));
        std::string acgt = "ACGT";
        for(auto contig : reader) {
            Contig read = contig.makeCompressedContig();
            if(read.size() < w + hasher.k - 1)
                continue;
            Path<htype128> path = dbg.align(read.seq).path();
            os << read.id << " " << path.start().hash() << int(path.start().isCanonical()) << " ";
            for (size_t i = 0; i < path.size(); i++) {
                os << acgt[path[i].seq[0]];
            }
            os << "\n";
            Path<htype128> path1 = dbg.align(!read.seq).path();
            os << "-" << read.id << " " << path1.start().hash() << int(path1.start().isCanonical()) << " ";
            for (size_t i = 0; i < path1.size(); i++) {
                os << acgt[path1[i].seq[0]];
            }
            os << "\n";
        }
        os.close();
    }

//    findTips(logger, dbg, threads);
    if(parser.getCheck("correct")) {
        io::SeqReader reader(lib);
        error_correction::correctSequences(dbg, logger, reader.begin(), reader.end(),
                                           dir / "corrected.fasta", dir / "bad.fasta", threads, w + hasher.k - 1);
    }
    if (parser.getValue("segments") != "none") {
        logger << "Drawing components" << std::endl;
        io::SeqReader segs(parser.getValue("segments"));
        ensure_dir_existance(dir / "pictures");
        size_t cnt = 0;
        for(StringContig scontig : segs) {
            Contig s = scontig.makeCompressedContig();
            std::vector<std::pair<const Edge<htype128> *, size_t>> path = dbg.carefulAlign(s.seq);
            std::vector<htype128> hashs;
            hashs.reserve(path.size());
            for(auto & tmp : path) {
                hashs.push_back(tmp.first->end()->hash());
            }
            Component<htype128> comp = Component<htype128>::neighbourhood(dbg, hashs.begin(), hashs.end(), 600, 2);
            std::ofstream os;
            os.open(dir/ "pictures" / (logging::itos(cnt) + ".dot"));
            comp.printCompressedDot(os, 2);
            os.close();
            logger << cnt << " " << comp.size() << std::endl;
            cnt += 1;
        }
    }
    if (parser.getValue("reference") != "none") {
        analyseGenome(dbg, parser.getValue("reference"), logger);
    }
    if (parser.getCheck("simplify")) {
        logger << "Removing low covered edges" << std::endl;
        std::vector<Sequence> edges;
        std::vector<htype128> vertices_again;
        for(auto & it : dbg) {
            vertices_again.push_back(it.first);
            Vertex<htype128> &vert = it.second;
            for(Edge<htype128> & edge : vert.getOutgoing()) {
                if (edge.getCoverage() >= 4)
                    edges.push_back(vert.seq + edge.seq);
            }
            for(Edge<htype128> & edge : vert.rc().getOutgoing()) {
                if (edge.getCoverage() >= 4)
                    edges.push_back(vert.rc().seq + edge.seq);
            }
        }
        SparseDBG<htype128> simp_dbg(vertices_again.begin(), vertices_again.end(), hasher);
        simp_dbg.fillSparseDBGEdges(edges.begin(), edges.end(), logger, threads, 0);
        mergeAll(logger, simp_dbg, threads);
        std::ofstream simp_os;
        simp_os.open(dir / "simp_graph.fasta");
        simp_dbg.printFasta(simp_os);
        simp_os.close();
    }
    logger << "DBG construction finished" << std::endl;
    return 0;
}
