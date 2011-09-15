#
# Makefile for bowtie, bowtie-build, bowtie-inspect
#

INC =
GCC_PREFIX = $(shell dirname `which gcc`)
GCC_SUFFIX =
CC = $(GCC_PREFIX)/gcc$(GCC_SUFFIX)
CPP = $(GCC_PREFIX)/g++$(GCC_SUFFIX)
CXX = $(CPP)
HEADERS = $(wildcard *.h)
BOWTIE_PTHREADS = 1
BOWTIE_MM = 1
BOWTIE_SHARED_MEM = 0

# Detect Cygwin or MinGW
WINDOWS = 0
ifneq (,$(findstring CYGWIN,$(shell uname)))
WINDOWS = 1
# POSIX memory-mapped files not currently supported on Windows
BOWTIE_MM = 0
BOWTIE_SHARED_MEM = 0
else
ifneq (,$(findstring MINGW,$(shell uname)))
WINDOWS = 1
# POSIX memory-mapped files not currently supported on Windows
BOWTIE_MM = 0
BOWTIE_SHARED_MEM = 0
endif
endif

MACOS = 0
ifneq (,$(findstring Darwin,$(shell uname)))
MACOS = 1
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
PTHREAD_DEF =
ifeq (1,$(BOWTIE_PTHREADS))
PTHREAD_DEF = -DBOWTIE_PTHREADS
ifeq (1,$(WINDOWS))
# pthreads for windows forces us to be specific about the library
PTHREAD_LIB = -lpthreadGC2
PTHREAD_PKG = pthreadGC2.dll
else
# There's also -pthread, but that only seems to work on Linux
PTHREAD_LIB = -lpthread
endif
endif

PREFETCH_LOCALITY = 2
PREF_DEF = -DPREFETCH_LOCALITY=$(PREFETCH_LOCALITY)

LIBS = 
SEARCH_LIBS = $(PTHREAD_LIB)
BUILD_LIBS =

SHARED_CPPS = ccnt_lut.cpp ref_read.cpp alphabet.cpp shmem.cpp \
              edit.cpp ebwt.cpp ebwt_io.cpp ebwt_util.cpp \
              reference.cpp ds.cpp multikey_qsort.cpp
SEARCH_CPPS = qual.cpp pat.cpp \
              log.cpp sam.cpp \
              color.cpp color_dec.cpp \
              read_qseq.cpp aligner_seed_policy.cpp \
              aligner_seed.cpp aligner_sw.cpp aligner_sw_col.cpp \
			  aligner_sw_driver.cpp aligner_cache.cpp \
			  aligner_result.cpp ref_coord.cpp mask.cpp \
			  pe.cpp aln_sink.cpp read_sink.cpp dp_framer.cpp \
			  scoring.cpp sa_rescomb.cpp \
			  seed_scan.cpp presets.cpp \
			  aligner_swsse.cpp \
			  aligner_swsse_loc_i16.cpp \
			  aligner_swsse_ee_i16.cpp \
			  aligner_swsse_loc_u8.cpp \
			  aligner_swsse_ee_u8.cpp
SEARCH_CPPS_MAIN = $(SEARCH_CPPS) bowtie_main.cpp

BUILD_CPPS = diff_sample.cpp
BUILD_CPPS_MAIN = $(BUILD_CPPS) bowtie_build_main.cpp

SEARCH_FRAGMENTS = $(wildcard search_*_phase*.c)
VERSION = $(shell cat VERSION)
EXTRA_FLAGS =

# Convert BITS=?? to a -m flag
BITS_FLAG =
ifeq (32,$(BITS))
BITS_FLAG = -m32
endif
ifeq (64,$(BITS))
BITS_FLAG = -m64
endif

DEBUG_FLAGS    = -O0 -g3 $(BITS_FLAG)
DEBUG_DEFS     = -DCOMPILER_OPTIONS="\"$(DEBUG_FLAGS) $(EXTRA_FLAGS)\""
RELEASE_FLAGS  = -O3 $(BITS_FLAG) -funroll-loops
RELEASE_DEFS   = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS) $(EXTRA_FLAGS)\""
NOASSERT_FLAGS = -DNDEBUG
FILE_FLAGS     = -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

BOWTIE2_BIN_LIST =     bowtie2-build \
                       bowtie2-align \
                       bowtie2-inspect
BOWTIE2_BIN_LIST_AUX = bowtie2-build-debug \
                       bowtie2-align-debug \
                       bowtie2-inspect-debug

GENERAL_LIST = $(wildcard scripts/*.sh) \
               $(wildcard scripts/*.pl) \
               $(wildcard indexes/e_coli*) \
               $(wildcard genomes/NC_008253.fna) \
               $(wildcard reads/e_coli_1000.*) \
               $(wildcard reads/e_coli_1000_*) \
               doc/manual.html \
               doc/README \
               doc/style.css \
               reads/e_coli_10000snp.fa \
               reads/e_coli_10000snp.fq \
               $(PTHREAD_PKG) \
               AUTHORS \
               COPYING \
               NEWS \
               MANUAL \
               MANUAL.markdown \
               TUTORIAL \
               VERSION

# This is helpful on Windows under MinGW/MSYS, where Make might go for
# the Windows FIND tool instead.
FIND=$(shell which find)

SRC_PKG_LIST = $(wildcard *.h) \
               $(wildcard *.hh) \
               $(wildcard *.c) \
               $(wildcard *.cpp) \
               doc/strip_markdown.pl \
               Makefile \
               $(GENERAL_LIST)

BIN_PKG_LIST = $(GENERAL_LIST)

.PHONY: all allall both both-debug

all: $(BOWTIE2_BIN_LIST)

allall: $(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_AUX)

both: bowtie2 bowtie2-build

both-debug: bowtie2-align-debug bowtie2-build-debug

DEFS=-fno-strict-aliasing \
     -DBOWTIE_VERSION="\"`cat VERSION`\"" \
     -DBUILD_HOST="\"`hostname`\"" \
     -DBUILD_TIME="\"`date`\"" \
     -DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\"" \
     $(FILE_FLAGS) \
     $(PTHREAD_DEF) \
     $(PREF_DEF) \
     $(MM_DEF) \
     $(SHMEM_DEF)

define checksum
  cat $^ | cksum | sed 's/[01-9][01-9] .*//' | awk '{print $$1}' > .$@.cksum
endef

#
# bowtie-build targets
#

bowtie2-build: ebwt_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(checksum)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
		-DEBWT_BUILD_HASH=`cat .$@.cksum` \
		$(DEFS) -DBOWTIE2 $(NOASSERT_FLAGS) -Wall \
		$(INC) \
		-o $@ $< \
		$(SHARED_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

bowtie2-build-debug: ebwt_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(checksum)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) \
		-DEBWT_BUILD_HASH=`cat .$@.cksum` \
		$(DEFS) -DBOWTIE2 -Wall \
		$(INC) \
		-o $@ $< \
		$(SHARED_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

#
# bowtie targets
#

bowtie2-align: ebwt_search.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(checksum)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
		-DEBWT_SEARCH_HASH=`cat .$@.cksum` \
		$(DEFS) -DBOWTIE2 $(NOASSERT_FLAGS) -Wall \
		$(INC) \
		-o $@ $< \
		$(SHARED_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

bowtie2-align-debug: ebwt_search.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(checksum)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(EXTRA_FLAGS) \
		-DEBWT_SEARCH_HASH=`cat .$@.cksum` \
		$(DEFS) -DBOWTIE2 -Wall \
		$(INC) \
		-o $@ $< \
		$(SHARED_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

#
# bowtie-inspect targets
#

bowtie2-inspect: bowtie_inspect.cpp $(HEADERS) $(SHARED_CPPS)
	$(checksum)
	$(CXX) $(RELEASE_FLAGS) \
		$(RELEASE_DEFS) $(EXTRA_FLAGS) \
		-DEBWT_INSPECT_HASH=`cat .$@.cksum` \
		$(DEFS) -DBOWTIE2 -DBOWTIE_INSPECT_MAIN -Wall \
		$(INC) -I . \
		-o $@ $< \
		$(SHARED_CPPS) \
		$(LIBS)

bowtie2-inspect-debug: bowtie_inspect.cpp $(HEADERS) $(SHARED_CPPS) 
	$(checksum)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(EXTRA_FLAGS) \
		-DEBWT_INSPECT_HASH=`cat .$@.cksum` \
		$(DEFS) -DBOWTIE2 -DBOWTIE_INSPECT_MAIN -Wall \
		$(INC) -I . \
		-o $@ $< \
		$(SHARED_CPPS) \
		$(LIBS)

bowtie2-src.zip: $(SRC_PKG_LIST)
	chmod a+x scripts/*.sh scripts/*.pl
	mkdir .src.tmp
	mkdir .src.tmp/bowtie-$(VERSION)
	zip tmp.zip $(SRC_PKG_LIST)
	mv tmp.zip .src.tmp/bowtie-$(VERSION)
	cd .src.tmp/bowtie-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .src.tmp ; zip -r $@ bowtie-$(VERSION)
	cp .src.tmp/$@ .
	rm -rf .src.tmp

bowtie2-bin.zip: $(BIN_PKG_LIST) $(BOWTIE_BIN_LIST) $(BOWTIE_BIN_LIST_AUX) 
	chmod a+x scripts/*.sh scripts/*.pl
	rm -rf .bin.tmp
	mkdir .bin.tmp
	mkdir .bin.tmp/bowtie-$(VERSION)
	if [ -f bowtie.exe ] ; then \
		zip tmp.zip $(BIN_PKG_LIST) $(addsuffix .exe,$(BOWTIE_BIN_LIST) $(BOWTIE_BIN_LIST_AUX)) ; \
	else \
		zip tmp.zip $(BIN_PKG_LIST) $(BOWTIE_BIN_LIST) $(BOWTIE_BIN_LIST_AUX) ; \
	fi
	mv tmp.zip .bin.tmp/bowtie-$(VERSION)
	cd .bin.tmp/bowtie-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .bin.tmp ; zip -r $@ bowtie-$(VERSION)
	cp .bin.tmp/$@ .
	rm -rf .bin.tmp

bowtie2-seeds-debug: aligner_seed.cpp ccnt_lut.cpp alphabet.cpp aligner_seed.h ebwt.cpp ebwt_io.cpp
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(EXTRA_FLAGS) \
		-DSCAN_MAIN \
		$(DEFS) -Wall \
		$(INC) -I . \
		-o $@ $< \
		aligner_seed.cpp ebwt.cpp ccnt_lut.cpp alphabet.cpp ebwt.cpp ebwt_io.cpp \
		$(LIBS)

.PHONY: doc
doc: doc/manual.html MANUAL

doc/manual.html: MANUAL.markdown
	echo "<h1>Table of Contents</h1>" > .tmp.head
	pandoc -T "Bowtie Manual" -B .tmp.head \
	       --css style.css -o $@ \
	       --from markdown --to HTML \
	       --table-of-contents $^

MANUAL: MANUAL.markdown
	perl doc/strip_markdown.pl < $^ > $@

.PHONY: clean
clean:
	rm -f $(BOWTIE_BIN_LIST)  $(BOWTIE_BIN_LIST_AUX) \
	rm -f $(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_AUX) \
	$(addsuffix .exe,$(BOWTIE_BIN_LIST)  $(BOWTIE_BIN_LIST_AUX)) \
	$(addsuffix .exe,$(BOWTIE2_BIN_LIST) $(BOWTIE2_BIN_LIST_AUX)) \
	bowtie-src.zip bowtie-bin.zip
	rm -f core.*
