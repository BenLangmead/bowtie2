cmake_minimum_required(VERSION 2.8) 

cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0005 NEW)

project(bowtie2 LANGUAGES CXX VERSION "2.3.4.3")

enable_testing()

include(CTest)
include(ProcessorCount)
include(CheckSymbolExists)

ProcessorCount(NUM_CORES)

option(BOWTIE_MM "enable bowtie2 memory mapping" ON)
option(BOWITE_SHARED_MM "enable shared memory mapping" ON)

set(CMAKE_CXX_STANDARD 98)
set(NO_TBB ${NO_TBB})
set(NO_SPINLOCK, ${NO_SPINLOCK})
set(WITH_THREAD_PROFILING ${WITH_THREAD_PROFILING})
set(POPCNT_CAPABILITY ${POPCNT_CAPABILITY})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -g3 -Wall -msse2")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -funroll-loops")
set(INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Installation directory for executables")

set(CMAKE_EXPORT_COMPILE_COMMANDS on)

if (MINGW)
    option(BOWTIE_MM "Memory mapped files not supported on Windows" OFF)
    option(BOWTIE_SHARED_MEM "Shared memory not supported on Windows" OFF)
endif(MINGW)

if (APPLE)
    set(CMAKE_XCODE_ATTRIBUTE_DEBUG_INFORMATION_FORMAT "dwarf-with-dsym")
endif(APPLE)

set(BOWTIE2_BIN_LIST
    bowtie2-build-s
    bowtie2-build-l
    bowtie2-align-s
    bowtie2-align-l
    bowtie2-inspect-s
    bowtie2-inspect-l
    )

set(SHARED_CPPS
    ccnt_lut.cpp
    ref_read.cpp
    alphabet.cpp
    shmem.cpp
    edit.cpp
    bt2_idx.cpp
    bt2_io.cpp
    bt2_util.cpp
    reference.cpp
    ds.cpp
    multikey_qsort.cpp
    limit.cpp
    random_source.cpp
    )

set(SEARCH_CPPS
    qual.cpp pat.cpp sam.cpp
    read_qseq.cpp aligner_seed_policy.cpp 
    aligner_seed.cpp 
    aligner_seed2.cpp 
    aligner_sw.cpp 
    aligner_sw_driver.cpp aligner_cache.cpp 
    aligner_result.cpp ref_coord.cpp mask.cpp 
    pe.cpp aln_sink.cpp dp_framer.cpp 
    scoring.cpp presets.cpp unique.cpp 
    simple_func.cpp 
    random_util.cpp 
    aligner_bt.cpp sse_util.cpp 
    aligner_swsse.cpp outq.cpp 
    aligner_swsse_loc_i16.cpp 
    aligner_swsse_ee_i16.cpp 
    aligner_swsse_loc_u8.cpp 
    aligner_swsse_ee_u8.cpp 
    aligner_driver.cpp
    bowtie_main.cpp
    bt2_search.cpp
    )

set(BUILD_CPPS
    bt2_build.cpp
    diff_sample.cpp
    bowtie_build_main.cpp)

set(INSPECT_CPPS
    bt2_inspect.cpp 
    )

string(TIMESTAMP BUILD_DATE)
execute_process(COMMAND ${CMAKE_CXX_COMPILER} -v OUTPUT_VARIABLE COMPILER_VERSION)

add_definitions(
    -DBOWTIE2
    -DBUILD_HOST="${CMAKE_HOST_SYSTEM}"
    -DBUILD_TIME="${BUILD_DATE}"
    -DCOMPILER_VERSION="${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_VERSION}"
    -DBOWTIE2_VERSION="${PROJECT_VERSION}"
    -D_LARGEFILE_SOURCE
    -D_FILE_OFFSET_BITS=64
    -D_GNU_SOURCE
    )

if (WITH_AFFINITY)
    add_definitions(-DWITH_AFFINITY=1)
endif()

if (NO_SPINLOCK)
    add_definitions(-DNO_SPINLOCK)
endif()

if (POPCNT_CAPABILITY)
    add_definitions(-DPOPCNT_CAPABILITY)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I third_party")
endif()

if (WITH_THREAD_PROFILING)
    add_definitions(-DPER_THREAD_TIMING=1)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    add_definitions(-DNDEBUG)
endif()

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
find_path(TBB_INCLUDE_PATH tbb)
find_library(TBB_LIBRARY_PATH tbb)
find_library(TBB_MALLOC_LIBRARY_PATH tbbmalloc)

if (ZLIB_FOUND)
    link_libraries(${ZLIB_LIBRARIES})
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

if (Threads_FOUND)
    link_libraries(Threads::Threads)
endif()

if (NOT TBB_LIBRARY_PATH)
    set(SHARED_CPPS ${SHARED_CPPS} tinythread.cpp)
else()
    if (NOT NO_QUEUELOCK)
        add_definitions(-DWITH_TBB -DWITH_QUEUELOCK=1 -DNO_SPINLOCK)
    endif()
    include_directories(${TBB_INCLUDE_PATH})
    link_libraries(${TBB_LIBRARY_PATH} ${TBB_MALLOC_LIBRARY_PATH})
endif()

include_directories(${PROJECT_SOURCE_DIR})
get_directory_property(COMPILER_DEFS COMPILE_DEFINITIONS)
add_definitions(-DCOMPILER_OPTIONS="${CMAKE_CXX_FLAGS}")

add_executable(bowtie2-align-s ${SEARCH_CPPS} ${SHARED_CPPS})
add_executable(bowtie2-align-l ${SEARCH_CPPS} ${SHARED_CPPS})
add_executable(bowtie2-build-s ${BUILD_CPPS} ${SHARED_CPPS})
add_executable(bowtie2-build-l ${BUILD_CPPS} ${SHARED_CPPS})
add_executable(bowtie2-inspect-s ${INSPECT_CPPS} ${SHARED_CPPS})
add_executable(bowtie2-inspect-l ${INSPECT_CPPS} ${SHARED_CPPS})

set_target_properties(bowtie2-align-l bowtie2-build-l bowtie2-inspect-l PROPERTIES COMPILE_FLAGS "-DBOWTIE2_64BIT_INDEX")
set_target_properties(bowtie2-inspect-s bowtie2-inspect-l PROPERTIES COMPILE_FLAGS "-DBOWTIE_INSPECT_MAIN")
set_target_properties(${BOWTIE2_BIN_LIST} PROPERTIES DEBUG_POSTFIX "-debug")

add_test(NAME simple-align COMMAND ${PROJECT_SOURCE_DIR}/bowtie2 -x example/index/lambda_virus example/reads/longreads.fq)
