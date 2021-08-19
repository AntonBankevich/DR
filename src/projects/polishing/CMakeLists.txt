project(polishing)
add_library(homopolish STATIC homopolish.cpp)
target_link_libraries(homopolish common sequence spoa ksw2)

add_executable(run_polishing run_polishing.cpp)
target_link_libraries(run_polishing homopolish)
