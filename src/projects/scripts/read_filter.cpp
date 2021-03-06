#include <alignment/aligner.hpp>
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <vector>

int main(int argc, char **argv) {
    CLParser parser({"output="}, {"reads"},
                    {"o=output", "t=threads"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    const std::experimental::filesystem::path out(parser.getValue("output"));
    std::ofstream os;
    os.open(out);
    io::Library libReads = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::SeqReader reads(libReads);
    for(const StringContig & tmp : reads) {
        size_t score = std::stoull(split(tmp.id).back());
        if (score < 20000) {
            os << ">" << tmp.id << "\n" << tmp.seq << std::endl;
        }
    }
    os.close();
    return 0;
}