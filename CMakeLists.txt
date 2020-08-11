cmake_minimum_required(VERSION 3.1)
project(DR)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "-std=c++14 -lstdc++fs -ggdb3" )
SET(CMAKE_SHARED_LINKER_FLAGS "-static-libgcc static-libstdc++ -Wall -Wc++-compat -O2 -msse4.1 -DHAVE_KALLOC -DKSW_CPU_DISPATCH -D_FILE_OFFSET_BITS=64 -ltbb -fsigned-char -fsanitize=address")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tools)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
find_package(OpenMP)
add_subdirectory(tests)
add_subdirectory(tools)
add_subdirectory(projects)


