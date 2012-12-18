/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
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

#include "aligner_cache.h"

#ifndef NDEBUG
/**
 * Check that this QVal is internally consistent and consistent
 * with the contents of the given cache.
 */
bool QVal::repOk(const AlignmentCache& ac) const {
	if(rangen_ > 0) {
		assert_lt(i_, ac.qSize());
		assert_leq(i_ + rangen_, ac.qSize());
	}
	assert_geq(eltn_, rangen_);
	return true;
}
#endif

#ifndef NDEBUG
/**
 * Check that this SAVal is internally consistent and consistent
 * with the contents of the given cache.
 */
bool SAVal::repOk(const AlignmentCache& ac) const {
	assert(len == 0 || i < ac.saSize());
	assert_leq(i + len, ac.saSize());
	return true;
}
#endif

/**
 * Add a new association between a read sequnce ('seq') and a
 * reference sequence ('')
 */
bool AlignmentCache::addOnTheFly(
	QVal& qv,         // qval that points to the range of reference substrings
	const SAKey& sak, // the key holding the reference substring
	uint32_t topf,    // top range elt in BWT index
	uint32_t botf,    // bottom range elt in BWT index
	uint32_t topb,    // top range elt in BWT' index
	uint32_t botb,    // bottom range elt in BWT' index
	bool getLock)
{
	ThreadSafe ts(lockPtr(), shared_ && getLock);
	bool added = true;
	// If this is the first reference sequence we're associating with
	// the query sequence, initialize the QVal.
	if(!qv.valid()) {
		qv.init((uint32_t)qlist_.size(), 0, 0);
	}
	qv.addRange(botf-topf); // update tally for # ranges and # elts
	if(!qlist_.add(pool(), sak)) {
		return false; // Exhausted pool memory
	}
#ifndef NDEBUG
	for(size_t i = qv.offset(); i < qlist_.size(); i++) {
		if(i > qv.offset()) {
			assert(qlist_.get(i) != qlist_.get(i-1));
		}
	}
#endif
	assert_eq(qv.offset() + qv.numRanges(), qlist_.size());
	SANode *s = samap_.add(pool(), sak, &added);
	if(s == NULL) {
		return false; // Exhausted pool memory
	}
	assert(s->key.repOk());
	if(added) {
		s->payload.i = (uint32_t)salist_.size();
		s->payload.len = botf - topf;
		s->payload.topf = topf;
		s->payload.topb = topb;
		for(size_t j = 0; j < (botf-topf); j++) {
			if(!salist_.add(pool(), 0xffffffff)) {
				// Change the payload's len field
				s->payload.len = (uint32_t)j;
				return false; // Exhausted pool memory
			}
		}
		assert(s->payload.repOk(*this));
	}
	// Now that we know all allocations have succeeded, we can do a few final
	// updates
	
	return true; 
}

#ifdef ALIGNER_CACHE_MAIN

#include <iostream>
#include <getopt.h>
#include <string>
#include "random_source.h"

using namespace std;

enum {
	ARG_TESTS = 256
};

static const char *short_opts = "vCt";
static struct option long_opts[] = {
	{(char*)"verbose",  no_argument, 0, 'v'},
	{(char*)"tests",    no_argument, 0, ARG_TESTS},
};

static void printUsage(ostream& os) {
	os << "Usage: bowtie2-cache [options]*" << endl;
	os << "Options:" << endl;
	os << "  --tests       run unit tests" << endl;
	os << "  -v/--verbose  talkative mode" << endl;
}

int gVerbose = 0;

static void add(
	RedBlack<QKey, QVal>& t,
	Pool& p,
	const char *dna)
{
	QKey qk;
	qk.init(BTDnaString(dna, true));
	t.add(p, qk, NULL);
}

/**
 * Small tests for the AlignmentCache.
 */
static void aligner_cache_tests() {
	RedBlack<QKey, QVal> rb(1024);
	Pool p(64 * 1024, 1024);
	// Small test
	add(rb, p, "ACGTCGATCGT");
	add(rb, p, "ACATCGATCGT");
	add(rb, p, "ACGACGATCGT");
	add(rb, p, "ACGTAGATCGT");
	add(rb, p, "ACGTCAATCGT");
	add(rb, p, "ACGTCGCTCGT");
	add(rb, p, "ACGTCGAACGT");
	assert_eq(7, rb.size());
	rb.clear();
	p.clear();
	// Another small test
	add(rb, p, "ACGTCGATCGT");
	add(rb, p, "CCGTCGATCGT");
	add(rb, p, "TCGTCGATCGT");
	add(rb, p, "GCGTCGATCGT");
	add(rb, p, "AAGTCGATCGT");
	assert_eq(5, rb.size());
	rb.clear();
	p.clear();
	// Regression test (attempt to make it smaller)
	add(rb, p, "CCTA");
	add(rb, p, "AGAA");
	add(rb, p, "TCTA");
	add(rb, p, "GATC");
	add(rb, p, "CTGC");
	add(rb, p, "TTGC");
	add(rb, p, "GCCG");
	add(rb, p, "GGAT");
	rb.clear();
	p.clear();
	// Regression test
	add(rb, p, "CCTA");
	add(rb, p, "AGAA");
	add(rb, p, "TCTA");
	add(rb, p, "GATC");
	add(rb, p, "CTGC");
	add(rb, p, "CATC");
	add(rb, p, "CAAA");
	add(rb, p, "CTAT");
	add(rb, p, "CTCA");
	add(rb, p, "TTGC");
	add(rb, p, "GCCG");
	add(rb, p, "GGAT");
	assert_eq(12, rb.size());
	rb.clear();
	p.clear();
	// Larger random test
	EList<BTDnaString> strs;
	char buf[5];
	for(int i = 0; i < 4; i++) {
		for(int j = 0; j < 4; j++) {
			for(int k = 0; k < 4; k++) {
				for(int m = 0; m < 4; m++) {
					buf[0] = "ACGT"[i];
					buf[1] = "ACGT"[j];
					buf[2] = "ACGT"[k];
					buf[3] = "ACGT"[m];
					buf[4] = '\0';
					strs.push_back(BTDnaString(buf, true));
				}
			}
		}
	}
	// Add all of the 4-mers in several different random orders
	RandomSource rand;
	for(uint32_t runs = 0; runs < 100; runs++) {
		rb.clear();
		p.clear();
		assert_eq(0, rb.size());
		rand.init(runs);
		EList<bool> used;
		used.resize(256);
		for(int i = 0; i < 256; i++) used[i] = false;
		for(int i = 0; i < 256; i++) {
			int r = rand.nextU32() % (256-i);
			int unused = 0;
			bool added = false;
			for(int j = 0; j < 256; j++) {
				if(!used[j] && unused == r) {
					used[j] = true;
					QKey qk;
					qk.init(strs[j]);
					rb.add(p, qk, NULL);
					added = true;
					break;
				}
				if(!used[j]) unused++;
			}
			assert(added);
		}
	}
}

/**
 * A way of feeding simply tests to the seed alignment infrastructure.
 */
int main(int argc, char **argv) {
	int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(argc, argv, short_opts, long_opts, &option_index);
		switch (next_option) {
			case 'v':       gVerbose = true; break;
			case ARG_TESTS: aligner_cache_tests(); return 0;
			case -1: break;
			default: {
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
			}
		}
	} while(next_option != -1);
}
#endif
