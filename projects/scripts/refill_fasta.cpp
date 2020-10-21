#include <alignment/aligner.hpp>
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <unordered_set>
#include <vector>

int main(int argc, char **argv) {
    CLParser parser({"output="}, {"reads"},{"o=output"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    std::unordered_map<std::string, std::string> rids;
    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::SeqReader reader(lib);
    for(StringContig contig : reader) {
        rids.emplace(contig.id, contig.seq);
    }
    const std::experimental::filesystem::path out(parser.getValue("output"));
    std::ofstream os;
    os.open(out);
    for (std::pair<const std::string, std::string> it : rids) {
        os << ">" << it.first << "\n" << it.second << "\n";
    }
    os.close();
    return 0;
}