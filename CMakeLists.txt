cmake_minimum_required(VERSION 3.15)
project(DR)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "-std=c++14 -lstdc++fs -ggdb3" )
SET(CMAKE_SHARED_LINKER_FLAGS "-static-libgcc static-libstdc++ -Wall -Wc++-compat -O2 -msse4.1 -DHAVE_KALLOC -DKSW_CPU_DISPATCH -D_FILE_OFFSET_BITS=64 -ltbb -fsigned-char -fsanitize=address")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
find_package(OpenMP)
add_subdirectory(common)
add_subdirectory(tests)
add_subdirectory(alignment)
add_subdirectory(sequences)
add_subdirectory(dbg)
add_subdirectory(minimap2)
include_directories(minimap2)

add_executable(DR main.cpp)
target_link_libraries(DR common minimap sequence alignment)

