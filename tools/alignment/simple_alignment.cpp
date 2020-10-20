//
// Created by anton on 08.07.2020.
//
#include <iostream>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>
#include <common/dir_utils.hpp>
#include "common/cl_parser.hpp"
#include "aligner.hpp"

int main(int argc, char **argv) {
    CLParser parser({"reads=", "ref=", "output-dir=", "threads=8"}, {}, {"o=output-dir", "t=threads"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    std::cout << "Reading reads" <<std::endl;
    std::vector<Contig> reads = io::SeqReader(parser.getValue("reads")).readAllContigs();
    std::cout << "Read " << reads.size() << " reads" <<std::endl;
    std::cout << "Reading reference" <<std::endl;
    std::vector<Contig> reference = io::SeqReader(parser.getValue("ref")).readAllContigs();
    std::cout << "Read " << reference.size() << " reference contigs" <<std::endl;
    std::cout << "Aligning reads" <<std::endl;
    std::ifstream is;
    is.open("parameters/hmm.txt");
    HMM hmm(HMM::load(is));
    is.close();
    std::vector<std::vector<MarkedAlignment<Contig, Contig>>> alignments =
            alignment_recipes::MarkAlign(reads, reference, hmm, "map-pb", std::stoi(parser.getValue("threads")));
    std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    cout << std::experimental::filesystem::absolute(dir) << endl;
    std::ofstream os;
    os.open(dir / "log.info");
    std::cout << "Printing results" <<std::endl;
    for(auto & v: alignments) {
        for(auto & al: v) {
            os << al << std::endl;
            for(const string & s : al.fourstringForm()) {
                os << s << std::endl;
            }
        }
    }
    os.close();
    std::cout << "Finished" <<std::endl;
    return 0;
}