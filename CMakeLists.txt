cmake_minimum_required(VERSION 2.8)

cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0005 NEW)

# set(CMAKE_OSX_DEPLOYMENT_TARGET "10.6")
project(bowtie2 LANGUAGES CXX VERSION "2.3.4.3")

enable_testing()

include(CTest)
include(ExternalProject)
include(ProcessorCount)
include(CheckSymbolExists)

ProcessorCount(NUM_CORES)

option(BOWTIE_MM "enable bowtie2 memory mapping" ON)
option(BOWITE_SHARED_MM "enable shared memory mapping" ON)

set(CMAKE_CXX_STANDARD 98)
set(NO_TBB ${NO_TBB})
set(NO_SPINLOCK ${NO_SPINLOCK})
set(USE_SRA 0)
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
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
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

if (USE_SRA)
  add_definitions(-DUSE_SRA)
  link_directories(${CMAKE_SOURCE_DIR}/ngs-2.9.2/install/lib64)
  link_directories(${CMAKE_SOURCE_DIR}/ncbi-vdb-2.9.2-1/install/lib64)
  include_directories(${CMAKE_SOURCE_DIR}/ngs-2.9.2/install/include)
  include_directories(${CMAKE_SOURCE_DIR}/ncbi-vdb-2.9.2-1/install/include)
  link_libraries(ncbi-ngs-c++-static ngs-c++-static ncbi-vdb-static dl)
endif()

if (WITH_THREAD_PROFILING)
  add_definitions(-DPER_THREAD_TIMING=1)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  add_definitions(-DNDEBUG)
endif()

find_package(ZLIB REQUIRED)

if (NO_TBB)
  find_package(Threads REQUIRED)
  link_libraries(Threads::Threads)
  set(SHARED_CPPS ${SHARED_CPPS} tinythread.cpp)
else()
  find_path(TBB_INCLUDE_PATH tbb)
  find_library(TBB_LIB_PATH tbb)
  find_library(TBBMALLOC_PROXY_LIB_PATH tbbmalloc_proxy)

  if (TBB_LIB_PATH)
    include_directories(${TBB_INCLUDE_PATH})
    link_libraries(${TBB_LIB_PATH} ${TBBMALLOC_PROXY_LIB_PATH})
  else()
    set(TBB_VER 2019_U4)

    ExternalProject_Add(tbb_project
      URL https://github.com/01org/tbb/archive/2019_U4.tar.gz
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/tbb-${TBB_VER}
      BUILD_IN_SOURCE 1
      CONFIGURE_COMMAND ""
      BUILD_COMMAND make tbb_build_prefix=${TBB_VER}
      INSTALL_COMMAND ""
      PREFIX=${CMAKE_CURRENT_BINARY_DIR}/tbb_project
      )

    ExternalProject_Get_Property(tbb_project SOURCE_DIR)

    set(TBB_SOURCE_DIR ${SOURCE_DIR})
    set(TBB_DEBUG_DIR ${TBB_SOURCE_DIR}/build/${TBB_VER}_debug)
    set(TBB_RELEASE_DIR ${TBB_SOURCE_DIR}/build/${TBB_VER}_release)
    include_directories(${TBB_SOURCE_DIR}/include)

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(TBB_LIBS ${TBB_DEBUG_DIR}/tbb_debug ${TBB_DEBUG_DIR}/tbbmalloc_proxy_debug)
      add_library(tbb_debug SHARED IMPORTED)
      add_library(tbbmalloc_proxy_debug SHARED IMPORTED)
      add_dependencies(tbb_debug tbb_project)
      add_dependencies(tbbmalloc_proxy_debug tbb_project)
      # set_property(TARGET tbb_debug PROPERTY IMPORTED_LOCATION ${TBB_DEBUG_DIR}/)
    else()
      set(TBB_LIBS ${TBB_RELEASE_DIR}/tbb ${TBB_RELEASE_DIR}/tbbmalloc_proxy)
      add_library(tbb SHARED IMPORTED)
      add_library(tbbmalloc_proxy SHARED IMPORTED)
      set_property(TARGET tbb PROPERTY IMPORTED_LOCATION ${TBB_RELEASE_DIR}/libtbb.dylib)
      set_property(TARGET tbbmalloc_proxy PROPERTY IMPORTED_LOCATION ${TBB_RELEASE_DIR}/libtbbmalloc_proxy.dylib)
      add_dependencies(tbb tbb_project)
      add_dependencies(tbbmalloc_proxy tbb_project)
    endif()
    link_libraries(tbb tbbmalloc_proxy)
  endif()
  add_definitions(-DWITH_TBB -DWITH_QUEUELOCK=1 -DNO_SPINLOCK)
endif()

if (ZLIB_FOUND)
  link_libraries(${ZLIB_LIBRARIES})
  include_directories(${ZLIB_INCLUDE_DIRS})
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

# target_link_libraries(bowtie2-align-s tbb tbbmalloc_proxy)

set_target_properties(bowtie2-align-l bowtie2-build-l bowtie2-inspect-l PROPERTIES COMPILE_FLAGS "-DBOWTIE2_64BIT_INDEX")
set_target_properties(bowtie2-inspect-s bowtie2-inspect-l PROPERTIES COMPILE_FLAGS "-DBOWTIE_INSPECT_MAIN")
set_target_properties(${BOWTIE2_BIN_LIST} PROPERTIES DEBUG_POSTFIX "-debug")

add_test(NAME simple-align COMMAND ${PROJECT_SOURCE_DIR}/bowtie2 -x example/index/lambda_virus example/reads/longreads.fq)
