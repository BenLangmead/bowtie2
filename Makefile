#
# Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
#
# This file is part of Bowtie 2.
#
# Bowtie 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Bowtie 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
#

#
# Makefile for bowtie, bowtie2-build, bowtie2-inspect
#

prefix := /usr/local
bindir := $(prefix)/bin

LDLIBS := -lz
GCC_PREFIX := $(shell dirname `which gcc`)
GCC_SUFFIX :=
CC ?= $(GCC_PREFIX)/gcc$(GCC_SUFFIX)
CPP ?= $(GCC_PREFIX)/g++$(GCC_SUFFIX)
CXX ?= $(CPP)
CXXFLAGS += -std=c++98
ifeq (aarch64,$(shell uname -m))
	CXXFLAGS += -fopenmp-simd
	CPPFLAGS += -Ithird_party/simde
endif

HEADERS := $(wildcard *.h)
BOWTIE_MM := 1
BOWTIE_SHARED_MEM :=

NGS_VER ?= 2.9.2
VDB_VER ?= 2.9.2-1

# Detect Cygwin or MinGW
WINDOWS :=
MINGW :=
ifneq (,$(findstring MINGW,$(shell uname)))
	WINDOWS := 1
	MINGW := 1
	# POSIX memory-mapped files not currently supported on Windows
	BOWTIE_MM :=
	BOWTIE_SHARED_MEM :=
endif

MACOS :=
ifneq (,$(findstring Darwin,$(shell uname)))
	MACOS := 1
	ifeq (1,$(shell uname -r | awk -F. '{ if ($$1 > 12 && $$1 < 16) print 1; }'))
		CXXFLAGS += -stdlib=libstdc++
	endif
	ifdef STATIC_BUILD
		CXXFLAGS += -mmacosx-version-min=10.9
	endif
endif

ifdef STATIC_BUILD
	LDFLAGS += -L$(CURDIR)/.tmp/lib
	CPPFLAGS += -I$(CURDIR)/.tmp/include
endif

POPCNT_CAPABILITY ?= 1
ifeq (1, $(POPCNT_CAPABILITY))
	CXXFLAGS += -DPOPCNT_CAPABILITY
	CPPFLAGS += -I third_party
endif

MM_DEF :=

ifeq (1,$(BOWTIE_MM))
	MM_DEF := -DBOWTIE_MM
endif

SHMEM_DEF :=

ifdef BOWTIE_SHARED_MEM
	SHMEM_DEF := -DBOWTIE_SHARED_MEM
endif

PTHREAD_PKG :=
PTHREAD_LIB :=

#if we're not using TBB, then we can't use queuing locks
ifeq (1,$(NO_TBB))
	NO_QUEUELOCK := 1
endif

ifeq (1,$(MINGW))
	PTHREAD_LIB :=
else
	PTHREAD_LIB := -lpthread
endif

ifeq (1,$(NO_SPINLOCK))
	CXXFLAGS += -DNO_SPINLOCK
endif

#default is to use Intel TBB
ifneq (1,$(NO_TBB))
	LDLIBS += $(PTHREAD_LIB) -ltbb
	ifdef STATIC_BUILD
		LDLIBS += -ltbbmalloc
	ifndef MINGW
		LDLIBS += -ldl
	endif
	else
		LDLIBS += -ltbbmalloc_proxy
	endif
	CXXFLAGS += -DWITH_TBB
else
	LDLIBS += $(PTHREAD_LIB)
endif

USE_SRA ?=
ifeq (1, $(USE_SRA))
ifdef MINGW
	$(error "SRA binaries cannot be built on MINGW")
else
	LDFLAGS += -L$(CURDIR)/.tmp/lib64

	ifndef ($(STATIC_BUILD))
		CPPFLAGS += -I$(CURDIR)/.tmp/include
	endif

	LDLIBS += -lncbi-ngs-c++-static
	LDLIBS += -lngs-c++-static
	LDLIBS += -lncbi-vdb-static
	LDLIBS += -ldl

	CXXFLAGS += -DUSE_SRA
endif
endif

ifeq (1,$(WITH_THREAD_PROFILING))
	CXXFLAGS += -DPER_THREAD_TIMING=1
endif

ifeq (1,$(WITH_AFFINITY))
	CXXFLAGS += -DWITH_AFFINITY=1
endif

#default is to use Intel TBB's queuing lock for better thread scaling performance
ifneq (1,$(NO_QUEUELOCK))
	CXXFLAGS += -DNO_SPINLOCK
	CXXFLAGS += -DWITH_QUEUELOCK=1
endif

SHARED_CPPS := ccnt_lut.cpp ref_read.cpp alphabet.cpp shmem.cpp \
               edit.cpp bt2_idx.cpp bt2_io.cpp bt2_util.cpp \
               reference.cpp ds.cpp multikey_qsort.cpp limit.cpp \
			   random_source.cpp

ifeq (1,$(NO_TBB))
	SHARED_CPPS += tinythread.cpp
endif

SEARCH_CPPS := qual.cpp pat.cpp sam.cpp \
               read_qseq.cpp aligner_seed_policy.cpp \
               aligner_seed.cpp \
			   aligner_seed2.cpp \
			   aligner_sw.cpp \
			   aligner_sw_driver.cpp aligner_cache.cpp \
			   aligner_result.cpp ref_coord.cpp mask.cpp \
			   pe.cpp aln_sink.cpp dp_framer.cpp \
			   scoring.cpp presets.cpp unique.cpp \
			   simple_func.cpp \
			   random_util.cpp \
			   aligner_bt.cpp sse_util.cpp \
			   aligner_swsse.cpp outq.cpp \
			   aligner_swsse_loc_i16.cpp \
			   aligner_swsse_ee_i16.cpp \
			   aligner_swsse_loc_u8.cpp \
			   aligner_swsse_ee_u8.cpp \
			   aligner_driver.cpp

SEARCH_CPPS_MAIN := $(SEARCH_CPPS) bowtie_main.cpp

DP_CPPS := qual.cpp aligner_sw.cpp aligner_result.cpp ref_coord.cpp mask.cpp \
          simple_func.cpp sse_util.cpp aligner_bt.cpp aligner_swsse.cpp \
		  aligner_swsse_loc_i16.cpp aligner_swsse_ee_i16.cpp \
		  aligner_swsse_loc_u8.cpp aligner_swsse_ee_u8.cpp scoring.cpp

BUILD_CPPS := diff_sample.cpp
BUILD_CPPS_MAIN := $(BUILD_CPPS) bowtie_build_main.cpp

SEARCH_FRAGMENTS := $(wildcard search_*_phase*.c)
VERSION := $(shell cat VERSION)

BITS := 32
ifeq (x86_64,$(shell uname -m))
	BITS := 64
endif
ifeq (amd64,$(shell uname -m))
	BITS := 64
endif
ifeq (aarch64,$(shell uname -m))
	BITS := 64
endif
# msys will always be 32 bit so look at the cpu arch instead.
ifneq (,$(findstring AMD64,$(PROCESSOR_ARCHITEW6432)))
	ifeq (1,$(MINGW))
		BITS := 64
	endif
endif
ifeq (32,$(BITS))
  $(error bowtie2 compilation requires a 64-bit platform )
endif

SSE_FLAG := -msse2
M64_FLAG := -m64
ifeq (aarch64,$(shell uname -m))
	SSE_FLAG =
	M64_FLAG =
endif

DEBUG_FLAGS    := -O0 -g3 $(M64_FLAG) $(SSE_FLAG)
RELEASE_FLAGS  := -O3 $(M64_FLAG) $(SSE_FLAG) -funroll-loops -g3
NOASSERT_FLAGS := -DNDEBUG
FILE_FLAGS     := -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE
DEBUG_DEFS     = -DCOMPILER_OPTIONS="\"$(DEBUG_FLAGS) $(CXXFLAGS)\""
RELEASE_DEFS   = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS) $(CXXFLAGS)\""

ifeq (0,$(shell $(CPP) -E -fsanitize=address,undefined btypes.h 2>&1 > /dev/null; echo $$?))
	SANITIZER_FLAGS := -fsanitize=address,undefined
endif
ifndef SANITIZER_FLAGS
	ifeq (0,$(shell $(CPP) -E -fsanitize=address btypes.h 2>&1 > /dev/null; echo $$?))
		SANITIZER_FLAGS := -fsanitize=address
	endif
endif
ifndef SANITIZER_FLAGS
	ifeq (0,$(shell $(CPP) -E -fsanitize=undefined btypes.h 2>&1 > /dev/null; echo $$?))
		SANITIZER_FLAGS := -fsanitize=undefined
	endif
endif

BOWTIE2_BIN_LIST :=     bowtie2-build-s \
                        bowtie2-build-l \
                        bowtie2-align-s \
                        bowtie2-align-l \
                        bowtie2-inspect-s \
                        bowtie2-inspect-l
BOWTIE2_BIN_LIST_DBG := bowtie2-build-s-debug \
                        bowtie2-build-l-debug \
                        bowtie2-align-s-debug \
                        bowtie2-align-l-debug \
                        bowtie2-inspect-s-debug \
                        bowtie2-inspect-l-debug
BOWTIE2_BIN_LIST_SAN := bowtie2-build-s-sanitized \
                        bowtie2-build-l-sanitized \
                        bowtie2-align-s-sanitized \
                        bowtie2-align-l-sanitized \
                        bowtie2-inspect-s-sanitized \
                        bowtie2-inspect-l-sanitized

GENERAL_LIST := $(wildcard scripts/*.sh) \
                $(wildcard scripts/*.pl) \
                doc/manual.html \
                doc/README \
		Dockerfile-aarch64 \
                doc/style.css \
		$(wildcard example/index/*.bt2) \
		$(wildcard example/reads/*.fq) \
		$(wildcard example/reads/*.pl) \
		example/reads/combined_reads.bam \
		example/reference/lambda_virus.fa \
                $(PTHREAD_PKG) \
		bowtie2 \
		bowtie2-build \
		bowtie2-inspect \
                AUTHORS \
                LICENSE \
                NEWS \
                MANUAL \
                MANUAL.markdown \
		README.md \
                TUTORIAL \
                VERSION

ifeq (1,$(WINDOWS))
	BOWTIE2_BIN_LIST := $(BOWTIE2_BIN_LIST) bowtie2.bat bowtie2-build.bat bowtie2-inspect.bat
    ifneq (1,$(NO_TBB))
	    CXXFLAGS += -static-libgcc -static-libstdc++
	else
	    CXXFLAGS += -static -static-libgcc -static-libstdc++
	endif
endif

# This is helpful on Windows under MinGW/MSYS, where Make might go for
# the Windows FIND tool instead.
FIND := $(shell which find)

SRC_PKG_LIST = $(wildcard *.h) \
               $(wildcard *.hh) \
               $(wildcard *.c) \
               $(wildcard *.cpp) \
               $(wildcard third_party/*) \
               Makefile \
               CMakeLists.txt \
               $(GENERAL_LIST)

ifeq (1,$(WINDOWS))
	BIN_PKG_LIST := $(GENERAL_LIST) bowtie2.bat bowtie2-build.bat bowtie2-inspect.bat
else
	BIN_PKG_LIST := $(GENERAL_LIST)
endif

.PHONY: all allall both both-debug

all: $(BOWTIE2_BIN_LIST) ;
allall: $(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_DBG) $(BOWTIE2_BIN_LIST_SAN) ;
both: bowtie2-align-s bowtie2-build-s bowtie2-align-l bowtie2-build-l ;
both-debug: bowtie2-align-s-debug bowtie2-build-s-debug bowtie2-align-l-debug bowtie2-build-l-debug ;
both-sanitized: bowtie2-align-s-sanitized bowtie2-build-s-sanitized bowtie2-align-l-sanitized bowtie2-build-l-sanitized ;

DEFS := -fno-strict-aliasing \
        -DBOWTIE2_VERSION="\"`cat VERSION`\"" \
        -DBUILD_HOST="\"${HOSTNAME:-`hostname`}\"" \
        -DBUILD_TIME="\"`date -u -r NEWS`\"" \
        -DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\"" \
        $(FILE_FLAGS) \
        $(PREF_DEF) \
        $(MM_DEF) \
        $(SHMEM_DEF)

# set compiler flags for all sanitized builds
$(BOWTIE2_BIN_LIST_SAN): CXXFLAGS += $(SANITIZER_FLAGS)

#
# bowtie2-build targets
#

bowtie2-build-s-sanitized bowtie2-build-s: bt2_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(CXXFLAGS) \
		$(DEFS) -DBOWTIE2 $(NOASSERT_FLAGS) -Wall \
		$(CPPFLAGS) \
		-o $@ $< \
		$(SHARED_CPPS) $(BUILD_CPPS_MAIN) \
		$(LDFLAGS) $(LDLIBS)

bowtie2-build-l-sanitized bowtie2-build-l: bt2_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(CXXFLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) -Wall \
		$(CPPFLAGS) \
		-o $@ $< \
		$(SHARED_CPPS) $(BUILD_CPPS_MAIN) \
		$(LDFLAGS) $(LDLIBS)

bowtie2-build-s-debug: bt2_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(CXXFLAGS) \
		$(DEFS) -DBOWTIE2 -Wall \
		$(CPPFLAGS) \
		-o $@ $< \
		$(SHARED_CPPS) $(BUILD_CPPS_MAIN) \
		$(LDFLAGS) $(LDLIBS)

bowtie2-build-l-debug: bt2_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(CXXFLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
		$(CPPFLAGS) \
		-o $@ $< \
		$(SHARED_CPPS) $(BUILD_CPPS_MAIN) \
		$(LDFLAGS) $(LDLIBS)

#
# bowtie2-align targets
#

bowtie2-align-s-sanitized bowtie2-align-s: bt2_search.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(CXXFLAGS) \
		$(DEFS) -DBOWTIE2 $(NOASSERT_FLAGS) -Wall \
		$(CPPFLAGS) \
		-o $@ $< \
		$(SHARED_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LDFLAGS) $(LDLIBS)

bowtie2-align-l-sanitized bowtie2-align-l: bt2_search.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(CXXFLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) -Wall \
		$(CPPFLAGS) \
		-o $@ $< \
		$(SHARED_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LDFLAGS) $(LDLIBS)

bowtie2-align-s-debug: bt2_search.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(CXXFLAGS) \
		$(DEFS) -DBOWTIE2 -Wall \
		$(CPPFLAGS) \
		-o $@ $< \
		$(SHARED_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LDFLAGS) $(LDLIBS)

bowtie2-align-l-debug: bt2_search.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(CXXFLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
		$(CPPFLAGS) \
		-o $@ $< \
		$(SHARED_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LDFLAGS) $(LDLIBS)

#
# bowtie2-inspect targets
#

bowtie2-inspect-s-sanitized bowtie2-inspect-s: bt2_inspect.cpp $(HEADERS) $(SHARED_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
		$(RELEASE_DEFS) $(CXXFLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_INSPECT_MAIN -Wall \
		$(CPPFLAGS) -I . \
		-o $@ $< \
		$(SHARED_CPPS) \
		$(LDFLAGS) $(LDLIBS)

bowtie2-inspect-l-sanitized bowtie2-inspect-l: bt2_inspect.cpp $(HEADERS) $(SHARED_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
		$(RELEASE_DEFS) $(CXXFLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_INSPECT_MAIN  -DBOWTIE_64BIT_INDEX -Wall \
		$(CPPFLAGS) -I . \
		-o $@ $< \
		$(SHARED_CPPS) \
		$(LDFLAGS) $(LDLIBS)

bowtie2-inspect-s-debug: bt2_inspect.cpp $(HEADERS) $(SHARED_CPPS)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(CXXFLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_INSPECT_MAIN -Wall \
		$(CPPFLAGS) -I . \
		-o $@ $< \
		$(SHARED_CPPS) \
		$(LDFLAGS) $(LDLIBS)

bowtie2-inspect-l-debug: bt2_inspect.cpp $(HEADERS) $(SHARED_CPPS)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(CXXFLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -DBOWTIE_INSPECT_MAIN -Wall \
		$(CPPFLAGS) -I . \
		-o $@ $< \
		$(SHARED_CPPS) \
		$(LDFLAGS) $(LDLIBS)

#
# bowtie2-dp targets
#

bowtie2-dp: bt2_dp.cpp $(HEADERS) $(SHARED_CPPS) $(DP_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
		$(RELEASE_DEFS) $(CXXFLAGS) $(NOASSERT_FLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_DP_MAIN -Wall \
		$(CPPFLAGS) -I . \
		-o $@ $< \
		$(DP_CPPS) $(SHARED_CPPS) \
		$(LDFLAGS) $(LDLIBS)

bowtie2-dp-debug: bt2_dp.cpp $(HEADERS) $(SHARED_CPPS) $(DP_CPPS)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(CXXFLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_DP_MAIN -Wall \
		$(CPPFLAGS) -I . \
		-o $@ $< \
		$(DP_CPPS) $(SHARED_CPPS) \
		$(LDFLAGS) $(LDLIBS)

bowtie2.bat:
	echo "@echo off" > bowtie2.bat
	echo "perl %~dp0/bowtie2 %*" >> bowtie2.bat

bowtie2-build.bat:
	echo "@echo off" > bowtie2-build.bat
	echo "python %~dp0/bowtie2-build %*" >> bowtie2-build.bat

bowtie2-inspect.bat:
	echo "@echo off" > bowtie2-inspect.bat
	echo "python %~dp0/bowtie2-inspect %*" >> bowtie2-inspect.bat

.PHONY: bowtie2-src-pkg
bowtie2-src-pkg: $(SRC_PKG_LIST)
	chmod a+x scripts/*.sh scripts/*.pl
	mkdir .src.tmp
	mkdir .src.tmp/bowtie2-$(VERSION)
	zip tmp.zip $(SRC_PKG_LIST)
	mv tmp.zip .src.tmp/bowtie2-$(VERSION)
	cd .src.tmp/bowtie2-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .src.tmp ; zip -r bowtie2-$(VERSION)-source.zip bowtie2-$(VERSION)
	cp .src.tmp/bowtie2-$(VERSION)-source.zip .
	rm -rf .src.tmp

.PHONY: bowtie2-bin-pkg
bowtie2-bin-pkg: PKG_DIR := bowtie2-$(VERSION)-$(if $(USE_SRA),sra-)$(if $(MACOS),macos,$(if $(MINGW),mingw,linux))-x86_64
bowtie2-bin-pkg: $(BIN_PKG_LIST) $(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_DBG)
	chmod a+x scripts/*.sh scripts/*.pl
	rm -rf .bin.tmp
	mkdir -p .bin.tmp/$(PKG_DIR)
	if [ -f bowtie2-align-s.exe ] ; then \
		zip tmp.zip $(BIN_PKG_LIST) $(addsuffix .exe,$(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_DBG)) ; \
	else \
		zip tmp.zip $(BIN_PKG_LIST) $(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_DBG) ; \
	fi
	mv tmp.zip .bin.tmp/$(PKG_DIR)
	cd .bin.tmp/$(PKG_DIR) ; unzip tmp.zip ; rm -f tmp.zip
	cd .bin.tmp ; zip -r $(PKG_DIR).zip $(PKG_DIR)
	cp .bin.tmp/$(PKG_DIR).zip .
	rm -rf .bin.tmp

bowtie2-seeds-debug: aligner_seed.cpp ccnt_lut.cpp alphabet.cpp aligner_seed.h bt2_idx.cpp bt2_io.cpp
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(CXXFLAGS) \
		-DSCAN_MAIN \
		$(DEFS) -Wall \
		$(CPPFLAGS) -I . \
		-o $@ $< \
		aligner_seed.cpp bt2_idx.cpp ccnt_lut.cpp alphabet.cpp bt2_io.cpp \
		$(LDFLAGS) $(LDLIBS)

.PHONY: doc
doc: doc/manual.html MANUAL

doc/manual.html: MANUAL.markdown
	echo "<h1>Table of Contents</h1>" > .tmp.head
	pandoc -B .tmp.head \
	       --metadata title:"Bowtie 2 Manual"\
	       --css doc/style.css -o $@ \
	       --from markdown --to HTML \
	       --table-of-contents $^
	rm -f .tmp.head

MANUAL: MANUAL.markdown
	pandoc -f markdown -t plain $^ -o $@

.PHONY: install
install: all
	mkdir -p $(DESTDIR)$(bindir)
	for file in $(BOWTIE2_BIN_LIST) bowtie2-inspect bowtie2-build bowtie2 ; do \
		cp -f $$file $(DESTDIR)$(bindir) ; \
	done

.PHONY: simple-test
simple-test: perl-deps both both-debug both-sanitized
	eval `perl -I $(CURDIR)/.tmp/lib/perl5 -Mlocal::lib=$(CURDIR)/.tmp` ; \
	sh ./scripts/test/simple_tests.sh

.PHONY: random-test
random-test: all perl-deps
	eval `perl -I $(CURDIR)/.tmp/lib/perl5 -Mlocal::lib=$(CURDIR)/.tmp` ; \
	sh ./scripts/sim/run.sh $(if $(NUM_CORES), $(NUM_CORES), 2)

.PHONY: perl-deps
perl-deps:
	if [ ! -d "$(CURDIR)/.tmp" ]; then \
		mkdir -p $(CURDIR)/.tmp/include $(CURDIR)/.tmp/lib ; \
	fi
	DL=$$( ( which wget >/dev/null 2>&1 && echo "wget --no-check-certificate -O-" ) || echo "curl -L") ; \
	$$DL http://cpanmin.us | perl - -l $(CURDIR)/.tmp App::cpanminus local::lib ; \
	eval `perl -I $(CURDIR)/.tmp/lib/perl5 -Mlocal::lib=$(CURDIR)/.tmp` ; \
	$(CURDIR)/.tmp/bin/cpanm --force File::Which Math::Random Clone Test::Deep Sys::Info ; \

.PHONY: static-libs
static-libs:
	if [ ! -d "$(CURDIR)/.tmp" ] ; then \
		mkdir -p $(CURDIR)/.tmp/include $(CURDIR)/.tmp/lib ; \
	fi ; \
	if [ `uname` = "Darwin" ]; then \
		export CFLAGS=-mmacosx-version-min=10.9 ; \
		export CXXFLAGS=-mmacosx-version-min=10.9 ; \
	fi ; \
	cd $(CURDIR)/.tmp ; \
	DL=$$( ( which wget >/dev/null 2>&1 && echo "wget --no-check-certificate" ) || echo "curl -LOk") ; \
	if [ ! -f "$(CURDIR)/.tmp/include/zlib.h" ] ; then \
		$$DL https://zlib.net/zlib-1.2.11.tar.gz && tar xzf zlib-1.2.11.tar.gz && cd zlib-1.2.11 ; \
		$(if $(MINGW), mingw32-make -f win32/Makefile.gcc, ./configure --static && make) ; \
		cp zlib.h zconf.h $(CURDIR)/.tmp/include && cp libz.a $(CURDIR)/.tmp/lib ; \
		rm -f zlib-1.2.11 ; \
	fi ; \
	if [ ! -d "$(CURDIR)/.tmp/include/tbb" ] ; then \
		cd $(CURDIR)/.tmp ; \
		$$DL https://github.com/01org/tbb/archive/2019_U4.tar.gz && tar xzf 2019_U4.tar.gz && cd tbb-2019_U4; \
		$(if $(MINGW), mingw32-make compiler=gcc arch=ia64 runtime=mingw, make) extra_inc=big_iron.inc -j4 \
		&& cp -r include/tbb $(CURDIR)/.tmp/include && cp build/*_release/*.a $(CURDIR)/.tmp/lib ; \
		rm -f 2019_U4.tar.gz ; \
	fi

.PHONY: sra-deps
sra-deps:
	DL=$$( ( which wget >/dev/null 2>&1 && echo "wget --no-check-certificate" ) || echo "curl -LOk") ; \
	if [ ! -d "$(CURDIR)/.tmp" ] ; then \
		mkdir -p $(CURDIR)/.tmp/include $(CURDIR)/.tmp/lib ; \
	fi ; \
	if [ `uname` = "Darwin" ]; then \
		export CFLAGS=-mmacosx-version-min=10.9 ; \
		export CXXFLAGS=-mmacosx-version-min=10.9 ; \
	fi ; \
	if [ ! -f "$(CURDIR)/.tmp/include/ngs/Alignment.hpp" ] ; then \
		if [ ! -d "$(CURDIR)/.tmp/ngs-$(NGS_VER)/ngs-sdk" ] ; then \
			cd $(CURDIR)/.tmp ; \
			$$DL https://github.com/ncbi/ngs/archive/$(NGS_VER).tar.gz ; \
			tar xzvf $(NGS_VER).tar.gz ; \
			rm -f $(NGS_VER).tar.gz ; \
		fi ; \
		cd $(CURDIR)/.tmp/ngs-$(NGS_VER) && ./configure --prefix=$(CURDIR)/.tmp --build-prefix=`pwd`/build ; \
		make && make install ; \
	fi ; \
	if [ ! -f "$(CURDIR)/.tmp/include/ncbi-vdb/NGS.hpp" ] ; then \
		if [ ! -d "$(CURDIR)/.tmp/ncbi-vdb-$(VDB_VER)/vdb3" ] ; then \
			cd $(CURDIR)/.tmp ; \
	 		$$DL https://github.com/ncbi/ncbi-vdb/archive/$(VDB_VER).tar.gz ; \
	 		tar zxvf $(VDB_VER).tar.gz ; \
	 		rm -f $(VDB_VER).tar.gz ; \
	 	fi ; \
	 	cd $(CURDIR)/.tmp/ncbi-vdb-$(VDB_VER) && ./configure --prefix=$(CURDIR)/.tmp --build-prefix=`pwd`/build --with-ngs-sdk=$(CURDIR)/.tmp && make && make install ; \
	 fi ;

.PHONY: test
test: simple-test random-test

.PHONY: clean
clean:
	rm -f $(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_DBG) $(BOWTIE2_BIN_LIST_SAN) \
	$(addsuffix .exe,$(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_DBG)) \
	bowtie2-*.zip
	rm -f core.* .tmp.head
	rm -rf *.dSYM
	rm -rf .tmp
