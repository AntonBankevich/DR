#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <unordered_set>
#include <vector>
#include <iostream>



using namespace std;
using logging::Logger;



int main(int argc, char **argv) {
    CLParser parser({}, {"reads"},
                    {});

    parser.parseCL(argc, argv);

    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
//        io::Library lib = oneline::initialize<std::experimental::filesystem::path>("reads.fasta");
    io::SeqReader reader(lib);
    StringContig cur;
    while (!reader.eof()) {
        cur = reader.read();
        cur.seq.erase(std::unique(cur.seq.begin(), cur.seq.end()), cur.seq.end());
        cout << ">" << cur.id <<endl << cur.seq << endl;
    }

 }
