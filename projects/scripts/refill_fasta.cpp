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
    CLParser parser({"output=", "subdataset="}, {"reads"},{"o=output"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    std::unordered_set<std::string> rids;
    io::SeqReader subreader(parser.getValue("subdataset"));
    for(StringContig contig : subreader) {
        rids.emplace(contig.id);
    }
    const std::experimental::filesystem::path out(parser.getValue("output"));
    std::ofstream os;
    os.open(out);
    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::SeqReader reader(lib);
    for (StringContig read : reader) {
        if (rids.find(read.id) != rids.end()) {
            os << ">" << read.id << "\n" << read.seq << "\n";
        }
    }
    os.close();
    return 0;
}