#include <alignment/aligner.hpp>
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <vector>

int main(int argc, char **argv) {
    CLParser parser({}, {"reads"}, {});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    const std::experimental::filesystem::path out(parser.getValue("output"));
    std::vector<size_t> lens;
    io::Library libReads = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::SeqReader reads(libReads);
    for(const StringContig & tmp : reads) {
        lens.push_back(tmp.size());
    }
    std::sort(lens.begin(), lens.end());
    std::cout << lens.size() << " " << lens[lens.size() / 2] << std::endl;
    return 0;
}