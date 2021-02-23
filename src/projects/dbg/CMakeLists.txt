project(debruijn)
set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP)
add_executable(dbg dbg.cpp crude_correct.cpp sparse_dbg.cpp graph_algorithms.cpp graph_algorithms.hpp dbg_disjointigs.cpp dbg_construction.cpp minimizer_selection.cpp)
target_link_libraries(dbg common sequence m)

add_executable(ccdbg ccdbg.cpp crude_correct.cpp sparse_dbg.cpp graph_algorithms.cpp graph_algorithms.hpp dbg_disjointigs.cpp dbg_construction.cpp minimizer_selection.cpp)
target_link_libraries(ccdbg common sequence m)

add_executable(extract extract.cpp)
target_link_libraries(extract common sequence m)
