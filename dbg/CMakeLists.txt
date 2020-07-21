project(debruijn)
set(CMAKE_CXX_STANDARD 17)

include_directories(.)
find_package(OpenMP)
add_executable(dbg de_bruijn.cpp)
target_link_libraries(dbg common sequence)
