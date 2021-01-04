//
// Created by Andrey Bzikadze on 12/17/20.
//

#include <common/cl_parser.hpp>
#include "exact_matcher.hpp"

using namespace scripts::exact_matcher;

int main(int argc, char **argv) {
    CLParser parser({"output=", "reference=", "query=", "base=239"}, {},
                    {"o=output", "r=reference", "q=query"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }

    PositionsVector positions = exact_match<htype128>(parser.getValue("query"),
                                                      parser.getValue("reference"));
    output_matches(positions, parser.getValue("output"));
}