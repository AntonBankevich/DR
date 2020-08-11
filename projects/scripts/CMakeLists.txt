project(scripts)
set(CMAKE_CXX_STANDARD 14)
add_executable(analyse_hifi main.cpp)
target_link_libraries(analyse_hifi common minimap sequence alignment)
