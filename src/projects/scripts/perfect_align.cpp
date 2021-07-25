#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <vector>
using namespace hashing;

struct ReadAndAls {
    std::string id;
    std::vector<Segment<Contig>> als;
    ReadAndAls(const std::string &id, const std::vector<Segment<Contig>> &als) : id(id), als(als) {}
};

int main(int argc, char **argv) {
    CLParser parser({"ref=", "kmer-size=", "width=", "output-dir=", "dimer-compress=1000000000,1000000000,1"}, {"reads"},
                    {"k=kmer-size", "w=width", "o=output-dir"});
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
    size_t k = std::stoi(parser.getValue("kmer-size"));
    RollingHash hasher(k, 239);
    size_t w = std::stoi(parser.getValue("width"));
    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    const std::experimental::filesystem::path reff(parser.getValue("ref"));
    std::unordered_map<htype, std::vector<std::pair<Contig *, size_t>>, alt_hasher<htype>> position_map;

    StringContig::homopolymer_compressing = true;
    StringContig::SetDimerParameters(parser.getValue("dimer-compress"));
    std::vector<Contig> ref;
    io::SeqReader reader(reff);
    std::ofstream refos;
    refos.open(dir/"contigs.fasta");
    for(StringContig scontig : reader) {
        Contig new_contig = scontig.makeContig();
        if (new_contig.seq <= !new_contig.seq) {
            ref.emplace_back(scontig.makeContig());
            ref.emplace_back(ref.back().RC());
            refos << ">" << scontig.id << "\n" << scontig.seq << std::endl;
        }
    }
    refos.close();
    for(Contig &contig : ref) {
        for(size_t pos = 1; pos + k <= contig.size(); pos += w) {
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
        std::vector<std::pair<Segment<Contig>, Segment<Contig>>> res;
        KWH kwh(hasher, read.seq, 0);
        while (true) {
            if (position_map.find(kwh.fHash()) != position_map.end()) {
                for(std::pair<Contig *, size_t> &pos : position_map[kwh.fHash()]) {
                    size_t shrink_left = 0;
                    if(kwh.pos > pos.second) {
                        shrink_left = kwh.pos - pos.second;
                    }
                    size_t shrink_right = 0;
                    if(-kwh.pos + pos.second + read.size()> pos.first->size())
                        shrink_right = -kwh.pos + pos.second + read.size() - pos.first->size();
                    Segment<Contig> seg_from(read, shrink_left, read.size() - shrink_right);
                    Segment<Contig> seg_to(*pos.first, pos.second - kwh.pos + shrink_left,
                                           pos.second - kwh.pos + read.size() - shrink_right);
                    res.emplace_back(seg_from, seg_to);
                }
            }
            if (!kwh.hasNext())
                break;
            kwh = kwh.next();
        }
        std::sort(res.begin(), res.end());
        res.erase(std::unique(res.begin(), res.end()), res.end());
        for(std::pair<Segment<Contig>, Segment<Contig>> &al : res) {
            if(al.first.seq() == al.second.seq())
                os << read.id << " " << al.second.contig().getId() << " " << al.second.left << " " << al.second.right << std::endl;
        }
    }
    os.close();

    return 0;
}
