//
// Created by anton on 9/24/20.
//

#include "common/cl_parser.hpp"
#include "common/dir_utils.hpp"
#include <experimental/filesystem>
#include <unordered_map>
#include <iostream>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>

int main(int argc, char **argv) {
    CLParser parser({"dir="}, {"reads"}, {});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    const std::experimental::filesystem::path dir(parser.getValue("dir"));
    std::unordered_map<std::string, Contig> map;
    std::vector<Contig> reads = io::SeqReader(parser.getValue("reference")).readAllContigs();
    for(Contig & read : reads)
        map[read.id] = read;
    for(const std::experimental::filesystem::directory_entry& subdir: std::experimental::filesystem::directory_iterator(dir)) {
        std::experimental::filesystem::path read_list_path = subdir.path() / "reads.txt";
        std::experimental::filesystem::path output_path = subdir.path() / "reads.fasta";
        std::ifstream is;
        is.open(read_list_path);
        std::ofstream os;
        os.open(output_path);
        size_t n;
        is >> n;
        for (size_t i = 0; i < n; i++) {
            size_t m;
            is >> m;
            std::string s;
            for(size_t j = 0; j < m; j++) {
                is >> s;
                Contig &read = map[s];
                os << ">" << read.id << "\n" << read.seq << "\n";
            }
        }
        is.close();
        os.close();
    }
}