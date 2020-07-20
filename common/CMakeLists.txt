project(common)
set(CMAKE_CXX_STANDARD 14)

include_directories(.)
add_library(common STATIC
        cl_parser.cpp oneline_utils.hpp)
find_package(OpenMP)
target_link_libraries (common m OpenMP::OpenMP_CXX stdc++fs)
