project(debruijn)
set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP)
add_executable(dbg dbg.cpp crude_correct.cpp)
add_executable(ccdbg ccdbg.cpp crude_correct.cpp)
target_link_libraries(dbg common sequence m)
target_link_libraries(ccdbg common sequence m)
