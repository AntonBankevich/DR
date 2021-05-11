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
#include <sstream>
#include <string>
#include <ssw/ssw_cpp.h>
#include <cassert>
using namespace std;
using logging::Logger;

int main(int argc, char **argv) {
//polished assembly as reference
    CLParser parser({"ref=", "snps=", "cutoff="}, {},
                    {});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
//dimer multiplicity
    size_t cutoff = stoi(parser.getValue("cutoff"));
    std::vector<Contig> assembly = io::SeqReader(parser.getValue("ref")).readAllContigs();
    std::string line;
    ifstream infile(parser.getValue("snps"));
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        int pos;
        string name;
        string chr;
        string aux1, aux2, aux3;
        iss >> chr;
        if (chr != "chrX") {
            cout << line << endl;
            continue;
        }
        iss>>name >> aux1 >> aux2 >> aux3 >> pos;
        size_t pos_name = 0;
        for (; pos_name < assembly.size(); pos_name ++) {
//            cout << assembly[pos_name].id <<  " " << name << endl;
            if (name == assembly[pos_name].id){
                break;
            }
        }
        size_t cont_length = assembly[pos_name].size();
        size_t min_dimer = 0;
        for (int dir = -1; dir < 2; dir +=2) {
            size_t next_pos = pos;
            assert(pos < assembly[pos_name].size());
//            cout << pos <<" " << assembly[pos_name].size();
            for (; next_pos < cont_length && next_pos >= 0; next_pos+= dir) {
                if (assembly[pos_name][next_pos] != assembly[pos_name][pos])
                    break;
            }
            if (assembly[pos_name][next_pos] == assembly[pos_name][pos])
                break;
//            cout << " starting " << assembly[pos_name][pos] << " " << assembly[pos_name][next_pos] << endl;
            size_t switches = 0;
            for (int add = -1; add <= 1; add += 2) {
                size_t cur_pos = pos;
                while (cur_pos > 0 && cur_pos < cont_length - 1) {
                    if (assembly[pos_name][cur_pos] == assembly[pos_name][cur_pos + add]) {
                        cur_pos += add;
                    } else if (assembly[pos_name][cur_pos + add] != assembly[pos_name][pos] &&
                               assembly[pos_name][cur_pos + add] != assembly[pos_name][next_pos]) {
                        break;
                    } else {
                        switches++;
                        cur_pos += add;
                    }
                }
//                cout << cur_pos << endl;
            }
//            cout << "Dimers multiplicity " << switches / 2 << endl;
            min_dimer = max(min_dimer, switches/2);
        }
        cout << line;
        if (min_dimer >= cutoff) {
            cout<< " dimer " << min_dimer;
        }
        cout << endl;
    }
    return 0;
}