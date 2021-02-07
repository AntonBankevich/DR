#include <alignment/aligner.hpp>
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <vector>

int main(int argc, char **argv) {
    CLParser parser({"reads=", "ref="}, {}, {});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    std::experimental::filesystem::path ref_file(parser.getValue("ref"));
    io::SeqReader ref_reader(ref_file);
    std::vector<StringContig> ref;
    for(auto it = ref_reader.begin(); it != ref_reader.end(); ++it) {
        ref.emplace_back(*it);
    }
    std::experimental::filesystem::path read_file(parser.getValue("ref"));
    io::SeqReader read_reader(ref_file);
    std::vector<StringContig> reads;
    for(auto it = read_reader.begin(); it != read_reader.end(); ++it) {
        reads.emplace_back(*it);
    }
    for(auto & read : reads) {
        std::cout << read.id << std::endl;
        std::string rc = basic::Reverse(read.seq);
        for(auto & chr : ref) {
            size_t pos = chr.seq.find(read.seq);
            if(pos != size_t(-1)) {
                std::cout << chr.id << "[" << pos <<"," << pos + read.size() << "]" << std::endl;
            }
        }
    }
    return 0;
}