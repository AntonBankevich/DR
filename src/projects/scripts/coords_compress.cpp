#include <alignment/aligner.hpp>
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <ssw/ssw_cpp.h>

using namespace std;
using logging::Logger;

int main(int argc, char **argv) {
    CLParser parser({"ref=", "position=", "name="}, {},
                    {});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    std::vector<Contig> assembly = io::SeqReader(parser.getValue("ref")).readAllContigs();
    int pos = stoi(parser.getValue("position"));
    string name = parser.getValue("name");
    size_t pos_name = 0;
    for (; pos_name < assembly.size(); pos_name ++) {
//        cout << assembly[pos_name].id <<  " " << name << endl;
        if (name == assembly[pos_name].id){

            break;
        }
    }
    if (pos_name == assembly.size()) {
        cout << "Wrong name";
        return 0;
    }
    int comp_pos = 0;
    for (int i = 0; i < pos; i++) {
        if (assembly[pos_name][i] != assembly[pos_name][i+1])
            comp_pos ++;
    }
    cout << name << " " <<comp_pos << " " << pos << " " << assembly[0][pos];
    return 0;
}