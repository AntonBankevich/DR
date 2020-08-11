project(debruijn)
set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP)
add_executable(dbg de_bruijn.cpp)
target_link_libraries(dbg -static-libgcc -static-libstdc++ common sequence m)
