project(alignments)
set(CMAKE_CXX_STANDARD 14)

include_directories(.)
add_library(alignment STATIC
        alignment_piece.cpp aligner.cpp aligner.hpp hmm.cpp simple_alignment.cpp)

find_package(OpenMP)
target_link_libraries (alignment common sequence minimap ${CMAKE_THREAD_LIBS_INIT} ${ZLIB_LIBRARIES} m OpenMP::OpenMP_CXX)
add_executable(align simple_alignment.cpp)
target_link_libraries(align alignment common sequence)