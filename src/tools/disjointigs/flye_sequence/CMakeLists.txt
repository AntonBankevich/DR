project(flye-sequences)
set(CMAKE_CXX_STANDARD 14)

add_library(flye-sequence STATIC alignment.cpp consensus_generator.cpp edlib.cpp overlap.cpp sequence.cpp sequence_container.cpp)
