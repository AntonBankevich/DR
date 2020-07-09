cmake_minimum_required(VERSION 3.15)
project(DR)
set(CMAKE_CXX_STANDARD 14)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
find_package(OpenMP)
add_subdirectory(common)
add_subdirectory(tests)
add_subdirectory(alignment)
add_subdirectory(sequences)


add_subdirectory(minimap2)
include_directories(minimap2)
SET(CMAKE_SHARED_LINKER_FLAGS "-Wall -Wc++-compat -O2 -msse4.1 -DHAVE_KALLOC -DKSW_CPU_DISPATCH -D_FILE_OFFSET_BITS=64 -fsigned-char -fsanitize=address")

add_executable(DR main.cpp)
target_link_libraries(DR common minimap sequence alignment)

