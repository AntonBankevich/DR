#include <alignment/aligner.hpp>
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <vector>

int main(int argc, char **argv) {
    CLParser parser({"paf="}, {}, {});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    std::ifstream is;
    is.open(parser.getValue("paf"));
    std::string line;
    std::vector<size_t> lens;
    std::string prev = "";
    size_t max = 0;
    while (std::getline(is, line))
    {
        std::vector<std::string> ss = split(line);
        size_t from = std::stoull(ss[7]);
        size_t to = std::stoull(ss[8]);
        if (ss[0] != prev) {
            lens.push_back(max);
            max = 0;
        }
        max = std::max(max, to - from);
        prev = ss[0];
    }
    lens.push_back(max);
    is.close();
    lens.erase(lens.begin());
    std::sort(lens.begin(), lens.end());
    std::cout << lens.size() << " " << lens[lens.size() / 2] << std::endl;
    return 0;
}