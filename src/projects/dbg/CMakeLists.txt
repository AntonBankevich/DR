project(debruijn)
set(CMAKE_CXX_STANDARD 14)


add_library(dbg_basic STATIC sparse_dbg.cpp graph_algorithms.cpp dbg_disjointigs.cpp dbg_construction.cpp minimizer_selection.cpp ../error_correction/multiplicity_estimation.cpp paths.cpp graph_alignment_storage.cpp component.cpp)
target_link_libraries (dbg_basic m ${OpenMP_CXX_FLAGS} stdc++fs)

