/*
 *  sim.cpp
 *  bowtie-beta1
 *
 *  Created by Benjamin Langmead on 1/22/11.
 *  Copyright 2011 Johns Hopkins University. All rights reserved.
 *
 */

#include "sim.h"
#include "util.h"

/**
 * Generate a random DNA sequence with length chosen uniformly randomly
 * from the given range.  Return the length.
 */
uint32_t SimRefBuilder::buildSeq(
	RandomSource& rnd,
	BTDnaMask& seq,
	TSzPair lenrange)
{
	float probs[4] = { 0.25f, 0.25f, 0.25f, 0.25f };
	uint32_t len = rnd.nextU32Range(lenrange.first, lenrange.second);
	for(size_t i = 0; i < len; i++) {
		seq.appendChar("ACGT"[rnd.nextFromProbs(probs, 4)]);
	}
	return len;
}

/**
 * Generates a random set of reference with number of references and
 * length of references chosen uniformly randomly from the given
 * ranges.
 */
void SimRefBuilder::build(
	RandomSource& rnd,
	SimRef& sr,
	TSzPair numref,
	TSzPair reflen)
{
	char buf[1024];
	// Determine number of reference sequences
	uint32_t nref = rnd.nextU32Range(numref.first, numref.second);
	for(size_t i = 0; i < nref; i++) {
		// Build one reference sequence
		sr.refs_.expand();
		size_t sz = SimRefBuilder::buildSeq(rnd, sr.refs_.back(), reflen);
		sr.totlen_ += sz;
		itoa10(i+1, buf);
		sr.names_.push_back(string(buf));
	}
}

/**
 * Generates a random set of reference with number of references and
 * length of references chosen uniformly randomly from the given
 * ranges.
 */
void SimReadBuilder::build(
	RandomSource& rnd,       // pseudo-random generator
	const SimRef& sr,        // reference
	EList<BTDnaMask>& reads, // append reads here
	TSzPair numreads,        // # reads
	bool    samelen,         // reads all the same length?
	TSzPair readlen)         // range of read lengths
{
}

#ifdef SIM_MAIN
int main(int argc, char**argv) {
	SimDriver sd;
	sd.run(5);
	return 0;
}
#endif
