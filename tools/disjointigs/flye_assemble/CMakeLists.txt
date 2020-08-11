project(flye-assemble)
set(CMAKE_CXX_STANDARD 14)

add_library(flye-assemble STATIC chimera.cpp extender.cpp main_assemble.cpp parameters_estimator.cpp)
target_link_libraries(flye-assemble flye-sequences)
