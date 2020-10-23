#include <common/cl_parser.hpp>
#include <iostream>
#include <omp.h>
#include <common/logging.hpp>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>

int main(int argc, char **argv) {
    CLParser parser({"threshold=", "out="}, {"reads"},
                    {"t=threads"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    logging::Logger logger;
    logger << "Reading and processing reads" << std::endl;
    io::Library libReads = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::SeqReader reads(libReads);
    std::ofstream os;
    os.open(parser.getValue("out"));
    size_t min = std::stoi(parser.getValue("threshold"));
    for (StringContig contig : reads) {
        if(contig.size() >= min) {
            os << ">" << contig.id << "\n" << contig.seq << "\n";
        }
    }
    os.close();
}
