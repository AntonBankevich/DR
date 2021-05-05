#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <vector>

int main(int argc, char **argv) {
    CLParser parser({"ref=", "kmer-size=", "width", "output-dir"}, {"reads"}, {"k=kmer-size", "w=width", "o=output-dir"});
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
        logger << argv[i] << " ";
    }
    logger << std::endl;
    size_t k = std::stoi(parser.getValue("k-mer-size"));
    RollingHash hasher(k, 239);
    size_t w = std::stoi(parser.getValue("width"));
    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    const std::experimental::filesystem::path reff(parser.getValue("ref"));
    std::unordered_map<htype, std::vector<std::pair<Contig *, size_t>>, alt_hasher<htype>> position_map;

    StringContig::needs_compressing = true;
    std::vector<Contig> ref;
    io::SeqReader reader(reff);
    for(StringContig scontig : reader) {
        ref.emplace_back(scontig.makeContig());
        ref.emplace_back(ref.back().RC());
    }
    for(Contig &contig : ref) {
        for(size_t pos = 0; pos + k <= contig.size(); pos += w) {
            htype h = hasher.hash(contig.seq, pos);
            position_map[h].emplace_back(&contig, pos);
        }
    }
    std::ofstream os;
    os.open(dir/"alignments.txt");
    io::SeqReader readReader(lib);
    for(StringContig scontig : readReader) {
        Contig read = scontig.makeContig();
        VERIFY(read.size() >= k + w - 1);
        std::vector<std::pair<Contig *, size_t>> res;
        KWH kwh(hasher, read.seq, 0);
        while (true) {
            if (position_map.find(kwh.fHash()) != position_map.end()) {
                for(std::pair<Contig *, size_t> &pos : position_map[kwh.fHash()]) {
                    if (kwh.pos >= pos.second && kwh.pos - pos.second + read.size()<= pos.first->size()) {
                        res.emplace_back(pos.first, pos.second - kwh.pos);
                    }
                }
            }
            if (!kwh.hasNext())
                break;
            kwh = kwh.next();
        }
        std::sort(res.begin(), res.end());
        res.erase(std::unique(res.begin(), res.end()), res.end());
        VERIFY(res.size() > 0);
        bool ok = false;
        for(std::pair<Contig *, size_t> &pos : res) {
            if(read.seq == pos.first->seq.Subseq(pos.second, read.size())) {
                os << read.id << " " << pos.first->id << " " << pos.second << " " << pos.second + read.size() << std::endl;
                ok = true;
            }
        }
        VERIFY(ok);
    }
    os.close();

    return 0;
}