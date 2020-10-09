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
    CLParser parser({"output=", "paf=", "chr="}, {"reads"},
                    {"o=output", "t=threads"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    std::unordered_map<std::string, Sequence> ref;
    std::ifstream is;
    is.open(parser.getValue("paf"));
    const std::experimental::filesystem::path out(parser.getValue("output"));
    std::ofstream os;
    std::string line;
    std::unordered_set<std::string> read_set;
    while (std::getline(is, line))
    {
        std::vector<std::string> ss = split(line);
        read_set.insert(ss[0]);
    }
    is.close();
    os.open(out);
    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::SeqReader reader(lib);
    for (StringContig read : reader) {
        if (read_set.find(read.id) != read_set.end()) {
            os << ">" << read.id << "\n" << read.seq << "\n";
        }
    }
    os.close();
    return 0;
}