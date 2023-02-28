cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0005 NEW)

file(STRINGS ${CMAKE_SOURCE_DIR}/BOWTIE2_VERSION PROJECT_VERSION)
project(bowtie2 LANGUAGES CXX VERSION ${PROJECT_VERSION})

enable_testing()

include(CTest)
include(ExternalProject)
include(ProcessorCount)
include(CheckSymbolExists)

ProcessorCount(NUM_CORES)

option(BOWTIE_MM "enable bowtie2 memory mapping" ON)
option(BOWTIE_SHARED_MEM "enable shared memory mapping" OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(USE_SRA ${USE_SRA})
set(WITH_THREAD_PROFILING ${WITH_THREAD_PROFILING})
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -g3 -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -funroll-loops")
set(INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Installation directory for executables")

include_directories(${SOURCE_DIR}/third_party)

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

set(BOWTIE2_WRAPPER_SCRIPTS
  bowtie2
  bowtie2-build
  bowtie2-inspect
  )

set(SHARED_CPPS
  ccnt_lut.cpp
  ref_read.cpp
  alphabet.cpp
  shmem.cpp
  edit.cpp
  bt2_idx.cpp
  bt2_io.cpp
  bt2_locks.cpp
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

if (BOWTIE_MM)
  add_definitions(-DBOWTIE_MM)
endif()

if (BOWTIE_SHARED_MEM)
  add_definitions(-DBOWTIE_SHARED_MEM)
endif()

if (WITH_AFFINITY)
  add_definitions(-DWITH_AFFINITY=1)
endif()

if (NO_SPINLOCK)
  add_definitions(-DNO_SPINLOCK)
endif()

if (NOT NO_QUEUELOCK)
  add_definitions(-DNO_SPINLOCK)
  add_definitions(-DWITH_QUEUELOCK=1)
endif()

 execute_process(COMMAND uname -m OUTPUT_VARIABLE ARCH)
if (${ARCH} MATCHES aarch64|arm64|390x|ppc64|ppc64le)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-simd")
elseif (${ARCH} MATCHES amd64|x64_64)
  add_definitions(-DPOPCNT_CAPABILITY)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
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
else()
  message("--------------------------------------------------------------------------")
  message("ATTN: Targets built in debug mode will have `-debug' appended to its name.")
  message("See https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html")
  message("--------------------------------------------------------------------------")
  set(CMAKE_EXECUTABLE_SUFFIX "-debug")
endif()

find_package(ZLIB REQUIRED)

find_package(Threads REQUIRED)
link_libraries(Threads::Threads)

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

install(TARGETS ${BOWTIE2_BIN_LIST} DESTINATION bin)
install(FILES ${BOWTIE2_WRAPPER_SCRIPTS}
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  DESTINATION bin)

add_test(NAME simple-align
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMAND bowtie2 -x example/index/lambda_virus example/reads/longreads.fq -u 10)
add_test(NAME simple-build
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMAND bowtie2-build-s -c GGGCGGCGACCTCGCGGGTTTTCGCTA out)
add_test(NAME simple-inspect
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMAND bowtie2-inspect-s out && rm out*)
