#include <common/cl_parser.hpp>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <unordered_map>

struct ContigRec {
    std::string name = "";
    double coverage = 0;
    size_t len = 0;
    size_t indeg = 0;
    size_t outdeg = 0;
    bool circ = false;
};
int main(int argc, char **argv) {
    std::experimental::filesystem::path gfa_file(argv[1]);
    std::ifstream is;
    is.open(gfa_file);
    std::unordered_map<std::string, ContigRec> recs;
    for( std::string line; getline(is, line); ) {
        std::vector<std::string> tokens = split(line);
        if(tokens[0] == "S") {
            ContigRec rec;
            rec.name = tokens[1];
            rec.len = tokens[2].size();
            rec.coverage = std::stod(tokens[3].substr(5));
            recs[rec.name] = std::move(rec);
        } else if(tokens[0] == "L") {
            ContigRec &rec = recs[tokens[1]];
            ContigRec &rec2 = recs[tokens[3]];
            if(tokens[1] == tokens[3] && tokens[2] == tokens[4]) {
                if(!rec.circ) {
                    rec.indeg += 1;
                    rec.outdeg += 1;
                    rec.circ = true;
                }
            } else {
                if(tokens[2] == "+")
                    rec.outdeg += 1;
                else
                    rec.indeg += 1;
                if(tokens[4] == "-")
                    rec2.outdeg += 1;
                else
                    rec2.indeg += 1;
            }
        }
    }
    for(const auto& it : recs) {
        const ContigRec &rec = it.second;
        std::cout << rec.name << " " << rec.coverage << " " << rec.len << " " << rec.circ << " " << rec.indeg << " " << rec.outdeg << "\n";
    }
    return 0;
}
