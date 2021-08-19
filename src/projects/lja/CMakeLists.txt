project(debruijn)
set(CMAKE_CXX_STANDARD 14)

add_executable(lja lja.cpp subdataset_processing.cpp)
target_link_libraries(lja ec dbg_basic homopolish common sequence m)

add_executable(dbg dbg.cpp subdataset_processing.cpp)
target_link_libraries(dbg ec dbg_basic common sequence m)
