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

#ifndef RANDOM_GEN_H_
#define RANDOM_GEN_H_

#include <stdint.h>
#include "assert_helpers.h"

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

#endif /*RANDOM_GEN_H_*/
