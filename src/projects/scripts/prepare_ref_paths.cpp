#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include "common/string_utils.hpp"
#include <vector>

using namespace hashing;

struct ReadAndAls {
    std::string id;
    std::vector<Segment<Contig>> als;
    ReadAndAls(const std::string &id, const std::vector<Segment<Contig>> &als) : id(id), als(als) {}
};

int main(int argc, char **argv) {
    CLParser parser({"ref=", "coordinates=", "output-file=", "radius="}, {},
                    {"o=output-file"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    const std::experimental::filesystem::path out_file(parser.getValue("output-file"));
    const std::experimental::filesystem::path ref_file(parser.getValue("ref"));
    const std::experimental::filesystem::path coordinates_file(parser.getValue("ref"));
    for(size_t i = 0; i < argc; i++) {
        std::cout << argv[i] << " ";
    }
    std::cout << std::endl;
    std::unordered_map<std::string, Contig> ref;
    StringContig::homopolymer_compressing = false;
    for(StringContig scontig : io::SeqReader({ref_file})) {
        ref[scontig.id] = scontig.makeContig();
    }
    size_t radius = std::stoull(parser.getValue("radius"));
    std::string line;
    std::ifstream is;
    std::ofstream os;
    is.open(coordinates_file);
    os.open(out_file);
    while(getline(is, line)) {
        std::vector<std::string> s = split(line);
        std::string chr_name = s[0];
        size_t left = std::stoull(s[1]);
        size_t right = std::stoull(s[2]);
        std::string our = s[3];
        std::string their = s[4];
        Contig & chr = ref[chr_name];
        VERIFY(left < right);
        VERIFY(right <= chr.size());
        left = left - std::min(left, radius);
        right = std::min(right + radius, chr.size());
        os << ">" << chr_name << "_" << left << "_" << right << "\n" << chr.seq.Subseq(left, right) << "\n";
    }
    is.close();
    os.close();
    return 0;
}
