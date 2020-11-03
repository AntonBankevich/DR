#include <alignment/aligner.hpp>
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <vector>

int main(int argc, char **argv) {
    CLParser parser({"output=", "paf=", "left=0", "right=none", "chr=none"}, {},
                    {"o=output"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    size_t left = std::stoull(parser.getValue("left"));
    size_t right;
    if(parser.getValue("right") == "none") {
        right = size_t(-1);
    } else
        right = std::stoull(parser.getValue("right"));
    std::string chr = parser.getValue("chr");
    std::ifstream is;
    is.open(parser.getValue("paf"));
    const std::experimental::filesystem::path out(parser.getValue("output"));
    std::cout << "Reading paf" << std::endl;
    std::ofstream os;
    os.open(out);
    std::string line;
    std::vector<std::pair<size_t , size_t>> segs;
    while (std::getline(is, line)) {
        std::vector<std::string> ss = split(line);
        if(chr != ss[5] && chr != "none") {
            continue;
        }
        size_t from = std::stoull(ss[7]);
        size_t to = std::stoull(ss[8]);
        if (from > right || to < left)
            continue;
        segs.emplace_back(std::max(from, left), std::min(to, right));
    }
    is.close();
    std::cout << "Read " << segs.size() << " records" << std::endl;
    std::vector<size_t> last_to_end;
    last_to_end.resize(right - left, 0);
    std::sort(segs.begin(), segs.end());
    for (auto p : segs) {
        last_to_end[p.first] = std::max(p.second, last_to_end[p.first]);
    }
    std::cout << "Filled last_to_end 1" << std::endl;
    size_t  rb = last_to_end[0];
    for(size_t i = 0; i < last_to_end.size(); i++) {
        rb = std::max(i, std::max(rb, last_to_end[i]));
        last_to_end[i] = rb - i;
    }
    std::cout << "Filled last_to_end 2" << std::endl;
    std::sort(last_to_end.begin(), last_to_end.end());
    size_t cur = 0;
    for(size_t i = 1; i < 15000; i++) {
        while(cur < last_to_end.size() && last_to_end[cur] < i)
            cur += 1;
        os << i << " " << cur - i + 1 << " " << last_to_end.size() - i + 1 << std::endl;
    }
    os.close();
    return 0;
}