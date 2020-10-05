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
    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    std::cout << "Reading reads from " << parser.getListValue("reads")  << std::endl;
    std::vector<Contig> reads = io::SeqReader(lib).readAllContigs();
    std::cout << "Putting reads into map" << std::endl;
    for(Contig & read : reads)
        map[split(read.id)[0]] = read;
    std::cout << "Writing reads to disk" << std::endl;
    for(const std::experimental::filesystem::directory_entry& subdir: std::experimental::filesystem::directory_iterator(dir)) {
        if (!isNonnegativeNumber(subdir.path().filename().string()))
            continue;
        std::cout << subdir << std::endl;
        std::experimental::filesystem::path read_list_path = subdir.path() / "reads.txt";
        std::experimental::filesystem::path output_path = subdir.path() / "reads.fasta";
        std::ifstream is;
        is.open(read_list_path);
        std::ofstream os;
        os.open(output_path);
        size_t n;
        is >> n;
        for (size_t i = 0; i < n; i++) {
            std::string s;
            is >> s;
            VERIFY(map.find(s) != map.end());
            Contig &read = map[s];
            os << ">" << s << "\n" << read.seq << "\n";
        }
        is.close();
        os.close();
    }
}