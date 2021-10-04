#include <alignment/aligner.hpp>
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <vector>

int main(int argc, char **argv) {
    CLParser parser({"contigs=", "min-contig=0", "genome-size="}, {},
                    {"t=threads"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    std::experimental::filesystem::path libAssembly(parser.getValue("contigs"));
    io::SeqReader contigs(libAssembly);
    size_t min_contig = std::stoull(parser.getValue("min-contig"));
    size_t genome_size = std::stoull(parser.getValue("genome-size"));
    std::vector<size_t> lens;
    for(const StringContig & contig : contigs) {
        if (contig.size() >= min_contig)
            lens.push_back(contig.size());
    }
    std::cout << "Number of contigs " << lens.size() << std::endl;
    std::sort(lens.rbegin(), lens.rend());
    size_t sum = 0;
    for(unsigned long len : lens) {
        sum += len;
        if(sum >= genome_size / 2) {
            std::cout << "NG50 " << len << std::endl;
            break;
        }
    }
    if(sum < genome_size / 2) {
        std::cout << "NG50 NA" << std::endl;
    }
    std::cout << "Total length " << std::accumulate(lens.begin(), lens.end(), size_t(0)) << std::endl;
    return 0;
}