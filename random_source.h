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

#ifndef RANDOM_GEN_H_
#define RANDOM_GEN_H_

#include <stdint.h>
#include "assert_helpers.h"

//#define MERSENNE_TWISTER

#ifndef MERSENNE_TWISTER

/**
 * Simple pseudo-random linear congruential generator, a la Numerical
 * Recipes.
 */
class RandomSource {
public:
	static const uint32_t DEFUALT_A = 1664525;
	static const uint32_t DEFUALT_C = 1013904223;

	RandomSource() :
		a(DEFUALT_A), c(DEFUALT_C), inited_(false) { }
	RandomSource(uint32_t _last) :
		a(DEFUALT_A), c(DEFUALT_C), last(_last), inited_(true) { }
	RandomSource(uint32_t _a, uint32_t _c) :
		a(_a), c(_c), inited_(false) { }

	void init(uint32_t seed = 0) {
		last = seed;
		inited_ = true;
		lastOff = 30;
	}

	uint32_t nextU32() {
		assert(inited_);
		uint32_t ret;
		last = a * last + c;
		ret = last >> 16;
		last = a * last + c;
		ret ^= last;
		lastOff = 0;
		return ret;
	}

	/**
	 * Return a pseudo-random unsigned 32-bit integer sampled uniformly
	 * from [lo, hi].
	 */
	uint32_t nextU32Range(uint32_t lo, uint32_t hi) {
		uint32_t ret = lo;
		if(hi > lo) {
			ret += (nextU32() % (hi-lo+1));
		}
		return ret;
	}

	/**
	 * Get next 2-bit unsigned integer.
	 */
	uint32_t nextU2() {
		assert(inited_);
		if(lastOff > 30) {
			nextU32();
		}
		uint32_t ret = (last >> lastOff) & 3;
		lastOff += 2;
		return ret;
	}

	/**
	 * Get next boolean.
	 */
	bool nextBool() {
		assert(inited_);
		if(lastOff > 31) {
			nextU32();
		}
		uint32_t ret = (last >> lastOff) & 1;
		lastOff++;
		return ret;
	}
	
	/**
	 * Return an unsigned int chosen by picking randomly from among
	 * options weighted by probabilies supplied as the elements of the
	 * 'weights' array of length 'numWeights'.  The weights should add
	 * to 1.
	 */
	uint32_t nextFromProbs(
		const float* weights,
		size_t numWeights)
	{
		float f = nextFloat();
		float tot = 0.0f; // total weight seen so far
		for(uint32_t i = 0; i < numWeights; i++) {
			tot += weights[i];
			if(f < tot) return i;
		}
		return (uint32_t)(numWeights-1);
	}

	float nextFloat() {
		assert(inited_);
		return (float)nextU32() / (float)0xffffffff;
	}

	static uint32_t nextU32(uint32_t last,
	                        uint32_t a = DEFUALT_A,
	                        uint32_t c = DEFUALT_C)
	{
		return (a * last) + c;
	}
	
	uint32_t currentA() const { return a; }
	uint32_t currentC() const { return c; }
	uint32_t currentLast() const { return last; }

private:
	uint32_t a;
	uint32_t c;
	uint32_t last;
	uint32_t lastOff;
	bool inited_;
};

#else

class RandomSource { // Mersenne Twister random number generator

public:

	// default constructor: uses default seed only if this is the first instance
	RandomSource() {
		reset();
	}
	
	// constructor with 32 bit int as seed
	RandomSource(uint32_t s) {
		init(s);
	}
	
	// constructor with array of size 32 bit ints as seed
	RandomSource(const uint32_t* array, int size) {
		init(array, size);
	}
	
	void reset() {
		state_[0] = 0;
		p_ = 0;
		inited_ = false;
	}
	
	virtual ~RandomSource() { }
	
	// the two seed functions
	void init(uint32_t); // seed with 32 bit integer
	void init(const uint32_t*, int size); // seed with array

	/**
	 * Return next 1-bit unsigned integer.
	 */
	bool nextBool() {
		return (nextU32() & 1) == 0;
	}
	
	/**
	 * Get next unsigned 32-bit integer.
	 */
	inline uint32_t nextU32() {
		assert(inited_);
		if(p_ == n) {
			gen_state(); // new state vector needed
		}
		// gen_state() is split off to be non-inline, because it is only called once
		// in every 624 calls and otherwise irand() would become too big to get inlined
		uint32_t x = state_[p_++];
		x ^= (x >> 11);
		x ^= (x << 7) & 0x9D2C5680UL;
		x ^= (x << 15) & 0xEFC60000UL;
		x ^= (x >> 18);
		return x;
	}
	
	/**
	 * Return next float between 0 and 1.
	 */
	float nextFloat() {
		assert(inited_);
		return (float)nextU32() / (float)0xffffffff;
	}
	
protected: // used by derived classes, otherwise not accessible; use the ()-operator

	static const int n = 624, m = 397; // compile time constants

	// the variables below are static (no duplicates can exist)
	uint32_t state_[n]; // state vector array
	int p_; // position in state array
	
	bool inited_; // true if init function has been called
	
	// private functions used to generate the pseudo random numbers
	uint32_t twiddle(uint32_t u, uint32_t v) {
		return (((u & 0x80000000UL) | (v & 0x7FFFFFFFUL)) >> 1) ^ ((v & 1UL) ? 0x9908B0DFUL : 0x0UL);
	}
	
	void gen_state(); // generate new state
	
};

#endif

#endif /*RANDOM_GEN_H_*/
