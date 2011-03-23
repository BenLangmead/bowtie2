/*
 *  mask.h
 */

#ifndef MASK_H_
#define MASK_H_

#include <iostream>
#include "random_source.h"

// 5-bit pop count
extern int alts5[32];

// Index of lowest set bit
extern int firsts5[32];

/**
 * Return 1 if a 2-bit-encoded base ('i') matches any bit in the mask ('j') and
 * the mask < 16.  Returns -1 if either the reference or the read character was
 * ambiguous.  Returns 0 if the characters unambiguously mismatch.
 */
static inline int matchesEx(int i, int j) {
	if(j >= 16 || i > 3) {
		// read and/or ref was ambiguous
		return -1;
	}
	return (((1 << i) & j) != 0) ? 1 : 0;
}

/**
 * Return 1 if a 2-bit-encoded base ('i') matches any bit in the mask ('j').
 */
static inline bool matches(int i, int j) {
	return ((1 << i) & j) != 0;
}

/**
 * Given a mask with up to 5 bits, return an index corresponding to a
 * set bit in the mask, randomly chosen from among all set bits.
 */
static inline int randFromMask(RandomSource& rand, int mask) {
	assert_gt(mask, 0);
	if(alts5[mask] == 1) {
		// only one to pick from, pick it via lookup table
		return firsts5[mask];
	}
	assert_gt(mask, 0);
	assert_lt(mask, 32);
	int r = rand.nextU32() % alts5[mask];
	assert_geq(r, 0);
	assert_lt(r, alts5[mask]);
	// could do the following via lookup table too
	for(int i = 0; i < 5; i++) {
		if((mask & (1 << i)) != 0) {
			if(r == 0) return i;
			r--;
		}
	}
	std::cerr << "Shouldn't get here" << std::endl;
	throw 1;
	return -1;
}

#endif /*ndef MASK_H_*/
