project(alignments)
set(CMAKE_CXX_STANDARD 14)

include_directories(.)
add_library(alignment STATIC
        alignment_piece.cpp contigs.cpp aligner.cpp aligner.hpp)

target_link_libraries (alignment minimap ${CMAKE_THREAD_LIBS_INIT} ${ZLIB_LIBRARIES} m)
