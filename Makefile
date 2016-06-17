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

prefix = /usr/local
bindir = $(prefix)/bin

INC =
GCC_PREFIX = $(shell dirname `which gcc`)
GCC_SUFFIX =
CC ?= $(GCC_PREFIX)/gcc$(GCC_SUFFIX)
CPP ?= $(GCC_PREFIX)/g++$(GCC_SUFFIX)
CXX ?= $(CPP)
HEADERS = $(wildcard *.h)
BOWTIE_MM = 1
BOWTIE_SHARED_MEM = 0

# Detect Cygwin or MinGW
WINDOWS = 0
MINGW = 0
ifneq (,$(findstring MINGW,$(shell uname)))
	WINDOWS = 1
	MINGW = 1
	# POSIX memory-mapped files not currently supported on Windows
	BOWTIE_MM = 0
	BOWTIE_SHARED_MEM = 0
endif

MACOS = 0
ifneq (,$(findstring Darwin,$(shell uname)))
	MACOS = 1
	ifneq (,$(findstring 13,$(shell uname -r)))
		CPP = clang++
		CC = clang
		override EXTRA_FLAGS += -stdlib=libstdc++
	endif
endif

POPCNT_CAPABILITY ?= 1
ifeq (1, $(POPCNT_CAPABILITY))
    override EXTRA_FLAGS += -DPOPCNT_CAPABILITY
    INC += -I third_party
endif

MM_DEF = 

ifeq (1,$(BOWTIE_MM))
	MM_DEF = -DBOWTIE_MM
endif

SHMEM_DEF = 

ifeq (1,$(BOWTIE_SHARED_MEM))
	SHMEM_DEF = -DBOWTIE_SHARED_MEM
endif

PTHREAD_PKG =
PTHREAD_LIB = 

ifeq (1,$(MINGW))
	PTHREAD_LIB = 
else
	PTHREAD_LIB = -lpthread
endif

ifeq (1,$(NO_SPINLOCK))
	override EXTRA_FLAGS += -DNO_SPINLOCK
endif

ifeq (1,$(WITH_TBB))
	LIBS = $(PTHREAD_LIB) -ltbb -ltbbmalloc_proxy
	override EXTRA_FLAGS += -DWITH_TBB
else
	LIBS = $(PTHREAD_LIB)
endif
SEARCH_LIBS = 
BUILD_LIBS = 
INSPECT_LIBS =

ifeq (1,$(MINGW))
	BUILD_LIBS = 
	INSPECT_LIBS = 
endif

ifeq (1,$(WITH_THREAD_PROFILING))
	override EXTRA_FLAGS += -DPER_THREAD_TIMING=1
endif

ifeq (1,$(WITH_AFFINITY))
	override EXTRA_FLAGS += -DWITH_AFFINITY=1
endif

SHARED_CPPS = ccnt_lut.cpp ref_read.cpp alphabet.cpp shmem.cpp \
              edit.cpp bt2_idx.cpp bt2_io.cpp bt2_util.cpp \
              reference.cpp ds.cpp multikey_qsort.cpp limit.cpp \
			  random_source.cpp
ifneq (1,$(WITH_TBB))
	SHARED_CPPS += tinythread.cpp
endif

SEARCH_CPPS = qual.cpp pat.cpp sam.cpp \
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
SEARCH_CPPS_MAIN = $(SEARCH_CPPS) bowtie_main.cpp

DP_CPPS = qual.cpp aligner_sw.cpp aligner_result.cpp ref_coord.cpp mask.cpp \
          simple_func.cpp sse_util.cpp aligner_bt.cpp aligner_swsse.cpp \
		  aligner_swsse_loc_i16.cpp aligner_swsse_ee_i16.cpp \
		  aligner_swsse_loc_u8.cpp aligner_swsse_ee_u8.cpp scoring.cpp

BUILD_CPPS = diff_sample.cpp
BUILD_CPPS_MAIN = $(BUILD_CPPS) bowtie_build_main.cpp

SEARCH_FRAGMENTS = $(wildcard search_*_phase*.c)
VERSION = $(shell cat VERSION)

BITS=32
ifeq (x86_64,$(shell uname -m))
	BITS=64
endif
ifeq (amd64,$(shell uname -m))
	BITS=64
endif
# msys will always be 32 bit so look at the cpu arch instead.
ifneq (,$(findstring AMD64,$(PROCESSOR_ARCHITEW6432)))
	ifeq (1,$(MINGW))
		BITS=64
	endif
endif
ifeq (32,$(BITS))
  $(error bowtie2 compilation requires a 64-bit platform )
endif

SSE_FLAG=-msse2 

DEBUG_FLAGS    = -O0 -g3 -m64 $(SSE_FLAG)
DEBUG_DEFS     = -DCOMPILER_OPTIONS="\"$(DEBUG_FLAGS) $(EXTRA_FLAGS)\""
RELEASE_FLAGS  = -O3 -m64 $(SSE_FLAG) -funroll-loops -g3
RELEASE_DEFS   = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS) $(EXTRA_FLAGS)\""
NOASSERT_FLAGS = -DNDEBUG
FILE_FLAGS     = -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

BOWTIE2_BIN_LIST =     bowtie2-build-s \
                       bowtie2-build-l \
                       bowtie2-align-s \
                       bowtie2-align-l \
                       bowtie2-inspect-s \
                       bowtie2-inspect-l
BOWTIE2_BIN_LIST_AUX = bowtie2-build-s-debug \
                       bowtie2-build-l-debug \
                       bowtie2-align-s-debug \
                       bowtie2-align-l-debug \
                       bowtie2-inspect-s-debug \
                       bowtie2-inspect-l-debug

GENERAL_LIST = $(wildcard scripts/*.sh) \
               $(wildcard scripts/*.pl) \
               doc/manual.html \
               doc/README \
               doc/style.css \
			   $(wildcard example/index/*.bt2) \
			   $(wildcard example/reads/*.fq) \
			   $(wildcard example/reads/*.pl) \
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
               TUTORIAL \
               VERSION

ifeq (1,$(WINDOWS))
	BOWTIE2_BIN_LIST := $(BOWTIE2_BIN_LIST) bowtie2.bat bowtie2-build.bat bowtie2-inspect.bat 
    ifeq (1,$(WITH_TBB)) 
	    override EXTRA_FLAGS += -static-libgcc -static-libstdc++
	else
	    override EXTRA_FLAGS += -static -static-libgcc -static-libstdc++
	endif
endif

# This is helpful on Windows under MinGW/MSYS, where Make might go for
# the Windows FIND tool instead.
FIND=$(shell which find)

SRC_PKG_LIST = $(wildcard *.h) \
               $(wildcard *.hh) \
               $(wildcard *.c) \
               $(wildcard *.cpp) \
               $(wildcard third_party/*) \
               doc/strip_markdown.pl \
               Makefile \
               $(GENERAL_LIST)

ifeq (1,$(WINDOWS))
	BIN_PKG_LIST = $(GENERAL_LIST) bowtie2.bat bowtie2-build.bat bowtie2-inspect.bat 
else
	BIN_PKG_LIST = $(GENERAL_LIST)
endif

.PHONY: all allall both both-debug

all: $(BOWTIE2_BIN_LIST)

allall: $(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_AUX)

both: bowtie2-align-s bowtie2-build-s bowtie2-align-l bowtie2-build-l

both-debug: bowtie2-align-s-debug bowtie2-build-s-debug bowtie2-align-l-debug bowtie2-build-l-debug

DEFS=-fno-strict-aliasing \
     -DBOWTIE2_VERSION="\"`cat VERSION`\"" \
     -DBUILD_HOST="\"`hostname`\"" \
     -DBUILD_TIME="\"`date`\"" \
     -DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\"" \
     $(FILE_FLAGS) \
     $(PREF_DEF) \
     $(MM_DEF) \
     $(SHMEM_DEF)

#
# bowtie2-build targets
#

bowtie2-build-s: bt2_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
		$(DEFS) -DBOWTIE2 $(NOASSERT_FLAGS) -Wall \
		$(INC) \
		-o $@ $< \
		$(SHARED_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

bowtie2-build-l: bt2_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) -Wall \
		$(INC) \
		-o $@ $< \
		$(SHARED_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

bowtie2-build-s-debug: bt2_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) \
		$(DEFS) -DBOWTIE2 -Wall \
		$(INC) \
		-o $@ $< \
		$(SHARED_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

bowtie2-build-l-debug: bt2_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
		$(INC) \
		-o $@ $< \
		$(SHARED_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

#
# bowtie2-align targets
#

bowtie2-align-s: bt2_search.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
		$(DEFS) -DBOWTIE2 $(NOASSERT_FLAGS) -Wall \
		$(INC) \
		-o $@ $< \
		$(SHARED_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

bowtie2-align-l: bt2_search.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) -Wall \
		$(INC) \
		-o $@ $< \
		$(SHARED_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

bowtie2-align-s-debug: bt2_search.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(EXTRA_FLAGS) \
		$(DEFS) -DBOWTIE2 -Wall \
		$(INC) \
		-o $@ $< \
		$(SHARED_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

bowtie2-align-l-debug: bt2_search.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(EXTRA_FLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
		$(INC) \
		-o $@ $< \
		$(SHARED_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

#
# bowtie2-inspect targets
#

bowtie2-inspect-s: bt2_inspect.cpp $(HEADERS) $(SHARED_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
		$(RELEASE_DEFS) $(EXTRA_FLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_INSPECT_MAIN -Wall \
		$(INC) -I . \
		-o $@ $< \
		$(SHARED_CPPS) \
		$(LIBS) $(INSPECT_LIBS)

bowtie2-inspect-l: bt2_inspect.cpp $(HEADERS) $(SHARED_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
		$(RELEASE_DEFS) $(EXTRA_FLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_INSPECT_MAIN  -DBOWTIE_64BIT_INDEX -Wall \
		$(INC) -I . \
		-o $@ $< \
		$(SHARED_CPPS) \
		$(LIBS) $(INSPECT_LIBS)

bowtie2-inspect-s-debug: bt2_inspect.cpp $(HEADERS) $(SHARED_CPPS) 
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(EXTRA_FLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_INSPECT_MAIN -Wall \
		$(INC) -I . \
		-o $@ $< \
		$(SHARED_CPPS) \
		$(LIBS) $(INSPECT_LIBS)

bowtie2-inspect-l-debug: bt2_inspect.cpp $(HEADERS) $(SHARED_CPPS)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(EXTRA_FLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -DBOWTIE_INSPECT_MAIN -Wall \
		$(INC) -I . \
		-o $@ $< \
		$(SHARED_CPPS) \
		$(LIBS) $(INSPECT_LIBS)

#
# bowtie2-dp targets
#

bowtie2-dp: bt2_dp.cpp $(HEADERS) $(SHARED_CPPS) $(DP_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
		$(RELEASE_DEFS) $(EXTRA_FLAGS) $(NOASSERT_FLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_DP_MAIN -Wall \
		$(INC) -I . \
		-o $@ $< \
		$(DP_CPPS) $(SHARED_CPPS) \
		$(LIBS) $(SEARCH_LIBS)

bowtie2-dp-debug: bt2_dp.cpp $(HEADERS) $(SHARED_CPPS) $(DP_CPPS)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(EXTRA_FLAGS) \
		$(DEFS) -DBOWTIE2 -DBOWTIE_DP_MAIN -Wall \
		$(INC) -I . \
		-o $@ $< \
		$(DP_CPPS) $(SHARED_CPPS) \
		$(LIBS) $(SEARCH_LIBS)

bowtie2.bat:
	echo "@echo off" > bowtie2.bat
	echo "perl %~dp0/bowtie2 %*" >> bowtie2.bat

bowtie2-build.bat:
	echo "@echo off" > bowtie2-build.bat
	echo "python %~dp0/bowtie2-build %*" >> bowtie2-build.bat

bowtie2-inspect.bat:
	echo "@echo off" > bowtie2-inspect.bat
	echo "python %~dp0/bowtie2-inspect %*" >> bowtie2-inspect.bat

.PHONY: bowtie2-src
bowtie2-src: $(SRC_PKG_LIST)
	chmod a+x scripts/*.sh scripts/*.pl
	mkdir .src.tmp
	mkdir .src.tmp/bowtie2-$(VERSION)
	zip tmp.zip $(SRC_PKG_LIST)
	mv tmp.zip .src.tmp/bowtie2-$(VERSION)
	cd .src.tmp/bowtie2-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .src.tmp ; zip -r bowtie2-$(VERSION)-source.zip bowtie2-$(VERSION)
	cp .src.tmp/bowtie2-$(VERSION)-source.zip .
	rm -rf .src.tmp

.PHONY: bowtie2-bin
bowtie2-bin: $(BIN_PKG_LIST) $(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_AUX) 
	chmod a+x scripts/*.sh scripts/*.pl
	rm -rf .bin.tmp
	mkdir .bin.tmp
	mkdir .bin.tmp/bowtie2-$(VERSION)
	if [ -f bowtie2-align-s.exe ] ; then \
		zip tmp.zip $(BIN_PKG_LIST) $(addsuffix .exe,$(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_AUX)) ; \
	else \
		zip tmp.zip $(BIN_PKG_LIST) $(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_AUX) ; \
	fi
	mv tmp.zip .bin.tmp/bowtie2-$(VERSION)
	cd .bin.tmp/bowtie2-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .bin.tmp ; zip -r bowtie2-$(VERSION).zip bowtie2-$(VERSION)
	cp .bin.tmp/bowtie2-$(VERSION).zip .
	rm -rf .bin.tmp

bowtie2-seeds-debug: aligner_seed.cpp ccnt_lut.cpp alphabet.cpp aligner_seed.h bt2_idx.cpp bt2_io.cpp
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(EXTRA_FLAGS) \
		-DSCAN_MAIN \
		$(DEFS) -Wall \
		$(INC) -I . \
		-o $@ $< \
		aligner_seed.cpp bt2_idx.cpp ccnt_lut.cpp alphabet.cpp bt2_io.cpp \
		$(LIBS)

.PHONY: doc
doc: doc/manual.html MANUAL

doc/manual.html: MANUAL.markdown
	echo "<h1>Table of Contents</h1>" > .tmp.head
	pandoc -T "Bowtie 2 Manual" -B .tmp.head \
	       --css style.css -o $@ \
	       --from markdown --to HTML \
	       --table-of-contents $^
	rm -f .tmp.head

MANUAL: MANUAL.markdown
	perl doc/strip_markdown.pl < $^ > $@

.PHONY: install
install: all
	mkdir -p $(DESTDIR)$(bindir)
	for file in $(BOWTIE2_BIN_LIST) bowtie2-inspect bowtie2-build bowtie2 ; do \
		cp -f $$file $(DESTDIR)$(bindir) ; \
	done

.PHONY: clean
clean:
	rm -f $(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_AUX) \
	$(addsuffix .exe,$(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_AUX)) \
	bowtie2-src.zip bowtie2-bin.zip
	rm -f core.* .tmp.head
	rm -rf *.dSYM
