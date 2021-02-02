//
// Created by anton on 8/7/20.
//

#include "common/hash_utils.hpp"
#include "common/cl_parser.hpp"
#include "common/logging.hpp"
#include "common/omp_utils.hpp"
#include "common/zip_utils.hpp"
#include "sequences/contigs.hpp"
#include "sequences/seqio.hpp"
#include <omp.h>
#include <vector>
#include <string>
#include <iostream>

int main(int argc, char **argv) {
    CLParser parser({"dbg=", "contigs=none", "threads=8", "k-mer-size=", "output-dir="}, {"segment"},
                    {"k=k-mer-size", "t=threads", "o=output-dir"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    size_t threads = std::stoi(parser.getValue("threads"));
    omp_set_num_threads(threads);
    logging::Logger logger;
    logger.info() << "Reading genome" << std::endl;
    StringContig::needs_compressing = true;
    SparseDBG<htype128> dbg = SparseDBG<htype128>::loadDBGFromFasta({std::experimental::filesystem::path(dbg_file)}, hasher, logger, threads);

    return 0;
}