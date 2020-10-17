#include <alignment/aligner.hpp>
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <vector>

int main(int argc, char **argv) {
    CLParser parser({"output=", "paf=", "left=0", "right=none", "chr=none"}, {"reference"},
                    {"o=output", "t=threads"});
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
    io::Library libRef = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reference"));
    io::SeqReader ref_reader(libRef);
    std::unordered_map<std::string, Sequence> ref;
    for(StringContig tmp : ref_reader) {
        std::cout << tmp.id << " " << tmp.size() << std::endl;
        ref.emplace(tmp.id, Sequence(tmp.seq));
    }
    std::ifstream is;
    is.open(parser.getValue("paf"));
    const std::experimental::filesystem::path out(parser.getValue("output"));
    std::ofstream os;
    os.open(out);
    std::string line;
    while (std::getline(is, line)) {
        std::vector<std::string> ss = split(line);
        if(chr != ss[5] && chr != "none") {
            continue;
        }
        size_t from = std::stoull(ss[7]);
        size_t to = std::stoull(ss[8]);
        if (from > right || to < left)
            continue;
        Sequence seq = ref[ss[5]].Subseq(from, to);
        if(ss[4] == "-")
            seq = !seq;
        os << ">" << ss[0] << "\n" << seq << std::endl;
    }
    os.close();
    is.close();
    return 0;
}