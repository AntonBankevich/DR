project(mosaic)
set(CMAKE_CXX_STANDARD 14)
add_executable(mosaic main.cpp)
target_link_libraries(mosaic common minimap sequence alignment)
