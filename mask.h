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
static inline int randFromMask(RandomSource& rnd, int mask) {
	assert_gt(mask, 0);
	if(alts5[mask] == 1) {
		// only one to pick from, pick it via lookup table
		return firsts5[mask];
	}
	assert_gt(mask, 0);
	assert_lt(mask, 32);
	int r = rnd.nextU32() % alts5[mask];
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
