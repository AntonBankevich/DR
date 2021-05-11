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


//compypaste!!!!
size_t compressRead(const string& read, size_t compression_length, size_t pos) {
//TODO:: is it fast?
//        logger.trace() << read << endl;
    string res = "";
    size_t current_coord = 0;
    for (size_t i = 0; i < read.length(); i++) {
        if (i == pos)
            return current_coord;
        bool compressed = false;
        if (i > 0 && read[i] == read[i - 1]) {
            compressed = true;
        } else if (current_coord > 2*compression_length) {
            bool in_repeat = (read[i] ==res[current_coord - 2] || read[i] == res[current_coord - 1]);
            if (in_repeat) {
                for (size_t j = 1; j < compression_length; j++) {
                    if (res[current_coord  - 2 * j] != res[current_coord -2]) {
                        in_repeat = false;
                        break;
                    }
                }
            }
            if (in_repeat) {
                for (size_t j = 1; j < compression_length; j++) {
                    if (res[current_coord -1 - 2 * j] != res[current_coord -1]) {
                        in_repeat = false;
                        break;
                    }
                }
            }
            compressed |= in_repeat;
        }
        if (!compressed) {
            res = std::move(res) + read[i];
            current_coord ++;
        }
    }
//        logger.trace() << res.str() << endl;
    return -1;
}

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

    int comp_pos = compressRead(assembly[pos_name].str(), 16, pos);

    cout << name << " " <<comp_pos << " " << pos << " " << assembly[0][pos];
    return 0;
}