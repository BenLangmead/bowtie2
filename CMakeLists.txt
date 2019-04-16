cmake_minimum_required(VERSION 2.8)

cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0005 NEW)

project(bowtie2 LANGUAGES CXX VERSION "2.3.5.1")

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
  set(NGS_VER 2.9.2)
  set(NCBI_VDB_VER 2.9.2-1)

  find_package(Java COMPONENTS Development REQUIRED)
  find_package(PythonInterp REQUIRED)
  find_package(Perl REQUIRED)

  find_program(MAKE_EXE NAMES gmake nmake make)
  ExternalProject_add(ngs_project
    URL https://github.com/ncbi/ngs/archive/${NGS_VER}.tar.gz
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/ngs-${NGS_VER}
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ./configure --prefix=install
    BUILD_COMMAND ${MAKE_EXE}
    INSTALL_COMMAND ${MAKE_EXE} install
    )

  ExternalProject_Get_Property(ngs_project SOURCE_DIR)
  set(NGS_INSTALL_DIR ${SOURCE_DIR}/ngs-sdk/install)

  link_directories(${NGS_INSTALL_DIR}/lib64)
  include_directories(${NGS_INSTALL_DIR}/include)

  add_library(ngs-c++-static STATIC IMPORTED)
  add_dependencies(ngs-c++-static ngs_project)
  find_library(LIBNGS_CPP ngs-c++-static PATHS ${NGS_INSTALL_DIR}/lib64 NO_DEFAULT_PATH)
  set_property(TARGET ngs-c++-static PROPERTY IMPORTED_LOCATION ${LIBNGS_CPP})

  link_libraries(ngs-c++-static)

  ExternalProject_add(ncbi_vdb_project
    URL https://github.com/ncbi/ncbi-vdb/archive/${NCBI_VDB_VER}.tar.gz
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/ncbi-vdb-${NCBI_VDB_VER}
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ./configure --prefix=install --with-ngs-sdk=${NGS_INSTALL_DIR}
    BUILD_COMMAND ${MAKE_EXE}
    INSTALL_COMMAND ${MAKE_EXE} install
    DEPENDS ngs_project
    )

  ExternalProject_Get_Property(ncbi_vdb_project SOURCE_DIR)
  set(NCBI_VDB_INSTALL_DIR ${SOURCE_DIR}/install)

  link_directories(${NCBI_VDB_INSTALL_DIR}/lib64)
  include_directories(${NCBI_VDB_INSTALL_DIR}/include)

  add_library(ncbi-ngs-c++-static STATIC IMPORTED)
  add_dependencies(ncbi-ngs-c++-static ncbi_vdb_project)
  find_library(LIBNCBI_NGS_CPP ncbi-ngs-c++-static PATHS ${NCBI_VDB_INSTALL_DIR}/lib64 NO_DEFAULT_ATH)
  set_property(TARGET ncbi-ngs-c++-static PROPERTY IMPORTED_LOCATION ${LIBNCBI_NGS_CPP})

  add_library(ncbi-vdb-static STATIC IMPORTED)
  add_dependencies(ncbi-vdb-static ncbi_vdb_project)
  find_library(LIBNCBI_VDB ncbi-vdb-static PATHS ${NCBI_VDB_INSTALL_DIR}/lib64 NO_DEFAULT_PATH)
  set_property(TARGET ncbi-vdb-static PROPERTY IMPORTED_LOCATION ${LIBNCBI_VDB})

  link_libraries(ncbi-ngs-c++-static ncbi-vdb-static)
  add_definitions(-DUSE_SRA)
endif()

if (WITH_THREAD_PROFILING)
  add_definitions(-DPER_THREAD_TIMING=1)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  add_definitions(-DNDEBUG)
endif()

find_package(ZLIB REQUIRED)


find_package(Threads REQUIRED)
link_libraries(Threads::Threads)

if (NO_TBB)
  set(SHARED_CPPS ${SHARED_CPPS} tinythread.cpp)
else()
  find_library(TBB_LIB_PATH tbb)

  if (TBB_LIB_PATH)
    find_path(TBB_INCLUDE_PATH tbb)
    link_libraries(${TBB_LIB_PATH})
    find_library(TBBMALLOC_PROXY_LIB_PATH tbbmalloc_proxy)
    if (NOT TBBMALLOC_PROXY_LIB_PATH)
      find_library(TBBMALLOC_LIB_PATH tbbmalloc)
      link_libraries(${TBBMALLOC_LIB_PATH})
    else()
      link_libraries(${TBBMALLOC_PROXY_LIB_PATH})
    endif()
    include_directories(${TBB_INCLUDE_PATH})
  else()
    set(TBB_VER 2019_U4)

    if (CMAKE_SYSTEM_NAME STREQUAL Linux)
      set(tbb_os linux)
    elseif (CMAKE_SYSTEM_NAME STREQUAL Darwin)
      set(tbb_os macos)
    elseif (CMAKE_SYSTEM_NAME STREQUAL Windows)
      set(tbb_os windows)
    else()
      set(tbb_os ${CMAKE_SYSTEM_NAME})
    endif()

    find_program(MAKE_EXE NAMES gmake nmake make)
    ExternalProject_Add(tbb_project
      URL https://github.com/01org/tbb/archive/2019_U4.tar.gz
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/tbb-${TBB_VER}
      BUILD_IN_SOURCE 1
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ${MAKE_EXE} tbb_build_prefix=${TBB_VER} tbb_os=${tbb_os}
      INSTALL_COMMAND ""
      PREFIX=${CMAKE_CURRENT_BINARY_DIR}/tbb_project
      )

    ExternalProject_Get_Property(tbb_project SOURCE_DIR)

    set(TBB_SOURCE_DIR ${SOURCE_DIR})
    set(TBB_RELEASE_DIR ${TBB_SOURCE_DIR}/build/${TBB_VER}_release)
    include_directories(${TBB_SOURCE_DIR}/include)

    link_directories(${TBB_RELEASE_DIR})
    add_library(tbb SHARED IMPORTED)
    add_dependencies(tbb tbb_project)
    find_library(LIBTBB tbb PATHS ${TBB_RELEASE_DIR} NO_DEFAULT_PATH)
    set_property(TARGET tbb PROPERTY IMPORTED_LOCATION ${LIBTBB})

    add_library(tbbmalloc_proxy SHARED IMPORTED)
    add_dependencies(tbbmalloc_proxy tbb_project)
    find_library(LIBTBBMALLOC_PROXY tbbmalloc_proxy PATHS ${TBB_RELEASE_DIR} NO_DEFAULT_PATH)
    if (NOT LIBTBBMALLOC_PROXY)
      add_library(tbbmalloc SHARED IMPORTED)
      add_dependencies(tbbmalloc tbb_project)
      find_library(LIBTBBMALLOC tbbmalloc PATHS ${TBB_RELEASE_DIR} NO_DEFAULT_PATH)
      set_property(TARGET tbbmalloc PROPERTY IMPORTED_LOCATION ${LIBTBBMALLOC})
      link_libraries(tbb tbbmalloc)
    else()
      set_property(TARGET tbbmalloc_proxy PROPERTY IMPORTED_LOCATION ${LIBTBBMALLOC_PROXY})
      link_libraries(tbb tbbmalloc_proxy)
    endif()
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

set_target_properties(bowtie2-align-l bowtie2-build-l bowtie2-inspect-l PROPERTIES COMPILE_FLAGS "-DBOWTIE2_64BIT_INDEX")
set_target_properties(bowtie2-inspect-s bowtie2-inspect-l PROPERTIES COMPILE_FLAGS "-DBOWTIE_INSPECT_MAIN")
set_target_properties(${BOWTIE2_BIN_LIST} PROPERTIES DEBUG_POSTFIX -debug)

add_test(NAME simple-align COMMAND ${PROJECT_SOURCE_DIR}/bowtie2 -x example/index/lambda_virus example/reads/longreads.fq)
