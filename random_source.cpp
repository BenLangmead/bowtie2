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

#include "random_source.h"
#include "random_util.h"

#ifdef MERSENNE_TWISTER

void RandomSource::gen_state() {
	for(int i = 0; i < (n - m); ++i) {
		state_[i] = state_[i + m] ^ twiddle(state_[i], state_[i + 1]);
	}
	for(int i = n - m; i < (n - 1); ++i) {
		state_[i] = state_[i + m - n] ^ twiddle(state_[i], state_[i + 1]);
	}
	state_[n - 1] = state_[m - 1] ^ twiddle(state_[n - 1], state_[0]);
	p_ = 0; // reset position
}

void RandomSource::init(uint32_t s) {  // init by 32 bit seed
	reset();
	state_[0] = s;
	for(int i = 1; i < n; ++i) {
		state_[i] = 1812433253UL * (state_[i - 1] ^ (state_[i - 1] >> 30)) + i;
	}
	p_ = n; // force gen_state() to be called for next random number
	inited_ = true;
}

void RandomSource::init(const uint32_t* array, int size) { // init by array
	init(19650218UL);
	int i = 1, j = 0;
	for(int k = ((n > size) ? n : size); k; --k) {
		state_[i] = (state_[i] ^ ((state_[i - 1] ^ (state_[i - 1] >> 30)) * 1664525UL)) + array[j] + j; // non linear
		++j; j %= size;
		if((++i) == n) { state_[0] = state_[n - 1]; i = 1; }
	}
	for(int k = n - 1; k; --k) {
		state_[i] = (state_[i] ^ ((state_[i - 1] ^ (state_[i - 1] >> 30)) * 1566083941UL)) - i;
		if((++i) == n) { state_[0] = state_[n - 1]; i = 1; }
	}
	state_[0] = 0x80000000UL; // MSB is 1; assuring non-zero initial array
	p_ = n; // force gen_state() to be called for next random number
	inited_ = true;
}

#endif

#ifdef MAIN_RANDOM_SOURCE

using namespace std;

int main(void) {
	cerr << "Test 1" << endl;
	{
		RandomSource rnd;
		int cnts[32];
		for(size_t i = 0; i < 32; i++) {
			cnts[i] = 0;
		}
		for(uint32_t j = 0; j < 10; j++) {
			rnd.init(j);
			for(size_t i = 0; i < 10000; i++) {
				uint32_t rndi = rnd.nextU32();
				for(size_t i = 0; i < 32; i++) {
					if((rndi & 1) != 0) {
						cnts[i]++;
					}
					rndi >>= 1;
				}
			}
			for(size_t i = 0; i < 32; i++) {
				cerr << i << ": " << cnts[i] << endl;
			}
		}
	}

	cerr << "Test 2" << endl;
	{
		int cnts[4][4];
		for(size_t i = 0; i < 4; i++) {
			for(size_t j = 0; j < 4; j++) {
				cnts[i][j] = 0;
			}
		}
		RandomSource rnd;
		Random1toN rn1n;
		for(size_t i = 0; i < 100; i++) {
			rnd.init((uint32_t)i);
			rn1n.init(4, true);
			uint32_t ri = rn1n.next(rnd);
			cnts[ri][0]++;
			ri = rn1n.next(rnd);
			cnts[ri][1]++;
			ri = rn1n.next(rnd);
			cnts[ri][2]++;
			ri = rn1n.next(rnd);
			cnts[ri][3]++;
		}
		for(size_t i = 0; i < 4; i++) {
			for(size_t j = 0; j < 4; j++) {
				cerr << cnts[i][j];
				if(j < 3) {
					cerr << ", ";
				}
			}
			cerr << endl;
		}
	}
}

#endif
