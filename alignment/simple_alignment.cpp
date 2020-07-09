//
// Created by anton on 08.07.2020.
//
#include <iostream>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>
#include "common/cl_parser.hpp"
#include "aligner.hpp"

int main(int argc, char **argv) {
    CLParser parser({"reads=", "ref=", "output-dir=", "threads=8"}, {"o=output-dir", "t=threads"});
    if (!parser.check()) {
        std::cout << "Incorrect parameters" << std::endl;
        return 1;
    }
    std::vector<Contig> reads = io::SeqReader(parser.getValue("reads")).readAll();
    std::vector<Contig> reference = io::SeqReader(parser.getValue("ref")).readAll();
    std::vector<std::vector<RawAlignment>> alignments = alignment_recipes::SimpleAlign(reads, reference, std::stoi(parser.getValue("threads")));
    return 0;
}