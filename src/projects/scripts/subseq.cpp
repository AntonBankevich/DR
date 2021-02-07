#include <alignment/aligner.hpp>
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <vector>

int main(int argc, char **argv) {
    CLParser parser({"ref=", "chr=", "from=", "to="}, {}, {});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    std::string chr_name = parser.getValue("chr");
    size_t from = std::stoull(parser.getValue("from"));
    size_t to = std::stoull(parser.getValue("to"));
    std::experimental::filesystem::path ref_file(parser.getValue("ref"));
    io::SeqReader ref_reader(ref_file);
    std::vector<StringContig> ref;
    for(auto it = ref_reader.begin(); it != ref_reader.end(); ++it) {
        ref.emplace_back(*it);
    }
    for(auto & chr : ref) {
        if(chr.id == chr_name) {
            std::cout << ">" << chr.id << "[" << from << "," << to << "]" << std::endl << chr.seq.substr(from, to - from);
        }
    }
    return 0;
}