cmake_minimum_required(VERSION 3.15)
project(minimap2_master)
set(CMAKE_CXX_STANDARD 14)

include_directories(.)
#include_directories(sse2neon)

#add_executable(minimap2_master
#        kthread.c kalloc.c misc.c bseq.c sketch.c sdust.c options.c index.c chain.c align.c hit.c map.c format.c pe.c esterr.c splitidx.c ksw2_ll_sse.c)
#set( CMAKE_CXX_FLAGS "-Wall -Wc++-compat -O2 -msse4.1 -DHAVE_KALLOC -DKSW_CPU_DISPATCH -D_FILE_OFFSET_BITS=64 -fsigned-char -fsanitize=address" )
#add_executable(minimap2_master
#        main.c kthread.c kalloc.c misc.c bseq.c sketch.c sdust.c options.c index.c chain.c align.c hit.c map.c format.c pe.c esterr.c splitidx.c ksw2_ll_sse.c
#        ksw2_extz2_sse.c ksw2_extd2_sse.c ksw2_exts2_sse.c ksw2_dispatch.c)
add_library(minimap STATIC
        minimap_interface.cpp kthread.c kalloc.c misc.c bseq.c sketch.c sdust.c options.c index.c chain.c align.c hit.c map.c format.c pe.c esterr.c splitidx.c ksw2_ll_sse.c
        ksw2_extz2_sse.c ksw2_extd2_sse.c ksw2_exts2_sse.c ksw2_dispatch.c kstring.h kseq.c)
find_package (Threads)
find_package (ZLIB)
target_link_libraries (minimap ${CMAKE_THREAD_LIBS_INIT} ${ZLIB_LIBRARIES} m)

