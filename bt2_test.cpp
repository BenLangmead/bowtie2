/*
 * Copyright 2011, Ben Langmead <blangmea@jhsph.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <utility>
#include "group_walk.h"
#include "random_source.h"
#include "reference.h"
#include "bt2_idx.h"

using namespace std;

/**
 * A reference.  A reference string and the corresponding Ebwts.
 */
struct BowtieIndex {
	
	BowtieIndex() : refstrs(MISC_CAT) {
		refstrs.clear();
		ebwts.first = ebwts.second = NULL;
		cebwts.first = cebwts.second = NULL;
		ref = NULL;
	}

	~BowtieIndex() { reset(); }
	
	/**
	 *
	 */
	void init(
		const char *ref,
		bool noref = false,
		bool nocolor = false,
		const char *fn = "bowtie-ebwt-test")
	{
		EList<string> refs;
		refs.push_back(string(ref));
		init(refs, noref, nocolor, fn);
	}
	
	/**
	 * Initialzie this bowtie index 
	 */
	void init(
		const EList<string>& refs,
		bool noref = false,
		bool nocolor = false,
		const char *fn = "bowtie-ebwt-test")
	{
		reset();
		refstrs = refs;
		ebwts = Ebwt::fromStrings<SString<char> >(
			refstrs,
			false, // colorspace? no
			false, // packed? no
			REF_READ_REVERSE,           // how to reverse
			Ebwt::default_bigEndian,    // store big endinan?
			Ebwt::default_lineRate,     // line rate (default)
			Ebwt::default_offRate,      // SA sample rate (default)
			Ebwt::default_ftabChars,    // # of chars consumed in ftab lookup (default)
			fn,                         // basename
			Ebwt::default_useBlockwise, // blockwise builder?
			Ebwt::default_bmax,         // bmax (default)
			Ebwt::default_bmaxMultSqrt, // bmaxSqrtMult (default)
			Ebwt::default_bmaxDivN,     // bmaxDivN (default)
			Ebwt::default_dcv,          // dcv (default)
			Ebwt::default_seed,         // pseudo-random seed
			false,  // verbose
			false,  // pass along memory exceptions?
			true);  // sanity check?
		ebwts.first->loadIntoMemory(
			false,  // colorspace? no
			-1,     // fw
			true,   // load SA sample?
			true,   // load ftab?
			true,   // load rstarts?
			true,   // load reference names?
			false); // verbose?
		ebwts.second->loadIntoMemory(
			false,  // colorspace? no
			1,      // yes, must be entire reverse
			false,  // load SA sample?
			true,   // load ftab?
			false,  // load rstarts?
			true,   // load reference names?
			false); // verbose?
		if(!nocolor) {
			cebwts = Ebwt::fromStrings<SString<char> >(
				refstrs,
				true,  // colorspace? yes
				false, // packed? no
				REF_READ_REVERSE,           // how to reverse
				Ebwt::default_bigEndian,    // store big endinan?
				Ebwt::default_lineRate,     // line rate (default)
				Ebwt::default_offRate,      // SA sample rate (default)
				Ebwt::default_ftabChars,    // # of chars consumed in ftab lookup (default)
				fn,                         // basename
				Ebwt::default_useBlockwise, // blockwise builder?
				Ebwt::default_bmax,         // bmax (default)
				Ebwt::default_bmaxMultSqrt, // bmaxSqrtMult (default)
				Ebwt::default_bmaxDivN,     // bmaxDivN (default)
				Ebwt::default_dcv,          // dcv (default)
				Ebwt::default_seed,         // pseudo-random seed
				false,  // verbose
				false,  // pass along memory exceptions?
				true);  // sanity check?
			cebwts.first->loadIntoMemory(
				true,   // colorspace? yes
				-1,     // fw
				true,   // load SA sample?
				true,   // load ftab?
				true,   // load rstarts?
				true,   // load reference names?
				false); // verbose?
			cebwts.second->loadIntoMemory(
				true,   // colorspace? yes
				1,      // yes, must be entire reverse
				false,  // load SA sample?
				true,   // load ftab?
				false,  // load rstarts?
				true,   // load reference names?
				false); // verbose?
		} else {
			if(cebwts.first != NULL) {
				delete cebwts.first;
				cebwts.first = NULL;
			}
			if(cebwts.second != NULL) {
				delete cebwts.second;
				cebwts.second = NULL;
			}
		}
		if(!noref) {
			ref = new BitPairReference(
				fn,
				true,    // colorspace?
				true,    // sanity check?
				&refstrs,// infiles
				NULL,    // origs
				true,    // infilesSeq
				false,   // useMm
				false,   // useShmem
				false,   // mmSweep
				false,   // verbose
				false);  // startVerbose
		} else {
			if(ref != NULL) {
				delete ref;
				ref = NULL;
			}
		}
	}
	
	/**
	 * Free all memory and set all pointers to NULL.
	 */
	void reset() {
		if(ebwts.first   != NULL) delete ebwts.first;
		if(ebwts.second  != NULL) delete ebwts.second;
		if(cebwts.first  != NULL) delete cebwts.first;
		if(cebwts.second != NULL) delete cebwts.second;
		if(ref != NULL) delete ref;
		ebwts.first = ebwts.second = NULL;
		cebwts.first = cebwts.second = NULL;
		ref = NULL;
		refstrs.clear();
	}

	EList<string> refstrs;     // reference strings
	pair<Ebwt*, Ebwt*> ebwts;  // nucleotide-space Ebwt objects
	pair<Ebwt*, Ebwt*> cebwts; // colorspace Ebwt objects
	BitPairReference *ref;     // reference sequences
};

/**
 * Unit tests for components of the SeedAligner.
 */
void ebwt_tests() {
	BowtieIndex bi;
	bi.init("");
	return;
}

#ifdef EBWT_TEST_MAIN

#include <getopt.h>
#include <string>

enum {
	ARG_NOFW = 256,
	ARG_NORC,
	ARG_MM,
	ARG_SHMEM,
	ARG_TESTS
};

static const char *short_opts = "vCt";
static struct option long_opts[] = {
	{(char*)"verbose",  no_argument, 0, 'v'},
	{(char*)"color",    no_argument, 0, 'C'},
	{(char*)"timing",   no_argument, 0, 't'},
	{(char*)"nofw",     no_argument, 0, ARG_NOFW},
	{(char*)"norc",     no_argument, 0, ARG_NORC},
	{(char*)"mm",       no_argument, 0, ARG_MM},
	{(char*)"shmem",    no_argument, 0, ARG_SHMEM},
	{(char*)"tests",    no_argument, 0, ARG_TESTS},
};

static void printUsage(ostream& os) {
	os << "Usage: ebwt-test [options]* <index> <patterns>" << endl;
	os << "Options:" << endl;
	os << "  --mm                memory-mapped mode" << endl;
	os << "  --shmem             shared memory mode" << endl;
	os << "  --nofw              don't align forward-oriented read" << endl;
	os << "  --norc              don't align reverse-complemented read" << endl;
	os << "  -t/--timing         show timing information" << endl;
	os << "  -C/--color          colorspace mode" << endl;
	os << "  -v/--verbose        talkative mode" << endl;
}

bool gNorc = false;
bool gNofw = false;
bool gColor = false;
int gVerbose = 0;
int gGapBarrier = 1;
bool gColorExEnds = true;
int gSnpPhred = 30;

extern void ebwt_tests();

/**
 * A way of feeding simply tests to the seed alignment infrastructure.
 */
int main(int argc, char **argv) {
	bool useMm = false;
	bool useShmem = false;
	bool mmSweep = false;
	bool noRefNames = false;
	bool sanity = false;
	bool timing = false;
	int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(argc, argv, short_opts, long_opts, &option_index);
		switch (next_option) {
			case 'v':       gVerbose = true; break;
			case 'C':       gColor   = true; break;
			case 't':       timing   = true; break;
			case ARG_NOFW:  gNofw    = true; break;
			case ARG_NORC:  gNorc    = true; break;
			case ARG_MM:    useMm    = true; break;
			case ARG_SHMEM: useShmem = true; break;
			case ARG_TESTS: {
				ebwt_tests();
				return 0;
			}
			case -1: break;
			default: {
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
			}
		}
	} while(next_option != -1);
	char *reffn;
	if(optind >= argc) {
		cerr << "No reference; quitting..." << endl;
		return 1;
	}
	reffn = argv[optind++];
	if(optind >= argc) {
		cerr << "No reads; quitting..." << endl;
		return 1;
	}
	string ebwtBase(reffn);
	BitPairReference ref(
		ebwtBase,    // base path
		gColor,      // whether we expect it to be colorspace
		sanity,      // whether to sanity-check reference as it's loaded
		NULL,        // fasta files to sanity check reference against
		NULL,        // another way of specifying original sequences
		false,       // true -> infiles (2 args ago) contains raw seqs
		useMm,       // use memory mapping to load index?
		useShmem,    // use shared memory (not memory mapping)
		mmSweep,     // touch all the pages after memory-mapping the index
		gVerbose,    // verbose
		gVerbose);   // verbose but just for startup messages
	Timer *t = new Timer(cerr, "Time loading fw index: ", timing);
	Ebwt ebwtFw(
		ebwtBase,
		gColor,      // index is colorspace
		0,           // don't need entireReverse for fw index
		true,        // index is for the forward direction
		5,           // offrate (irrelevant)
		0,           // offrate-plus (irrelevant)
		useMm,       // whether to use memory-mapped files
		useShmem,    // whether to use shared memory
		mmSweep,     // sweep memory-mapped files
		!noRefNames, // load names?
		false,       // load SA sample?
		true,        // load ftab?
		true,        // load rstarts?
		NULL,        // reference map, or NULL if none is needed
		gVerbose,    // whether to be talkative
		gVerbose,    // talkative during initialization
		false,       // handle memory exceptions, don't pass them up
		sanity);
	delete t;
	t = new Timer(cerr, "Time loading bw index: ", timing);
	Ebwt ebwtBw(
		ebwtBase + ".rev",
		gColor,      // index is colorspace
		1,           // need entireReverse
		false,       // index is for the backward direction
		5,           // offrate (irrelevant)
		0,           // offrate-plus (irrelevant)
		useMm,       // whether to use memory-mapped files
		useShmem,    // whether to use shared memory
		mmSweep,     // sweep memory-mapped files
		!noRefNames, // load names?
		false,       // load SA sample?
		true,        // load ftab?
		false,       // load rstarts?
		NULL,        // reference map, or NULL if none is needed
		gVerbose,    // whether to be talkative
		gVerbose,    // talkative during initialization
		false,       // handle memory exceptions, don't pass them up
		sanity);
	delete t;
	for(int i = optind; i < argc; i++) {
	}
}
#endif

