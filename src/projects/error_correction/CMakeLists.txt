project(debruijn)
set(CMAKE_CXX_STANDARD 14)

add_library(ec STATIC correction_utils.cpp crude_correct.cpp manyk_correction.cpp multiplicity_estimation.cpp initial_correction.cpp dimer_correction.cpp)
target_link_libraries (ec dbg_basic m)
