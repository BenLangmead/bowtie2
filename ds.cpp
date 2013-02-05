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

#include "ds.h"

MemoryTally gMemTally;

/**
 * Tally a memory allocation of size amt bytes.
 */
void MemoryTally::add(int cat, uint64_t amt) {
	ThreadSafe ts(&mutex_m);
	tots_[cat] += amt;
	tot_ += amt;
	if(tots_[cat] > peaks_[cat]) {
		peaks_[cat] = tots_[cat];
	}
	if(tot_ > peak_) {
		peak_ = tot_;
	}
}

/**
 * Tally a memory free of size amt bytes.
 */
void MemoryTally::del(int cat, uint64_t amt) {
	ThreadSafe ts(&mutex_m);
	assert_geq(tots_[cat], amt);
	assert_geq(tot_, amt);
	tots_[cat] -= amt;
	tot_ -= amt;
}
	
#ifdef MAIN_DS

#include <limits>
#include "random_source.h"

using namespace std;

int main(void) {
	cerr << "Test EHeap 1...";
	{
		EHeap<float> h;
		h.insert(0.5f);  // 1
		h.insert(0.6f);  // 2
		h.insert(0.25f); // 3
		h.insert(0.75f); // 4
		h.insert(0.1f);  // 5
		h.insert(0.9f);  // 6
		h.insert(0.4f);  // 7
		assert_eq(7, h.size());
		if(h.pop() != 0.1f) {
			throw 1;
		}
		assert_eq(6, h.size());
		if(h.pop() != 0.25f) {
			throw 1;
		}
		assert_eq(5, h.size());
		if(h.pop() != 0.4f) {
			throw 1;
		}
		assert_eq(4, h.size());
		if(h.pop() != 0.5f) {
			throw 1;
		}
		assert_eq(3, h.size());
		if(h.pop() != 0.6f) {
			throw 1;
		}
		assert_eq(2, h.size());
		if(h.pop() != 0.75f) {
			throw 1;
		}
		assert_eq(1, h.size());
		if(h.pop() != 0.9f) {
			throw 1;
		}
		assert_eq(0, h.size());
		assert(h.empty());
	}
	cerr << "PASSED" << endl;

	cerr << "Test EHeap 2...";
	{
		EHeap<size_t> h;
		RandomSource rnd(12);
		size_t lim = 2000;
		while(h.size() < lim) {
			h.insert(rnd.nextU32());
		}
		size_t last = std::numeric_limits<size_t>::max();
		bool first = true;
		while(!h.empty()) {
			size_t p = h.pop();
			assert(first || p >= last);
			last = p;
			first = false;
		}
	}
	cerr << "PASSED" << endl;

	cerr << "Test EBitList 1...";
	{
		EBitList<128> l;
		assert_eq(0, l.size());
		assert_eq(std::numeric_limits<size_t>::max(), l.max());
		
		assert(!l.test(0));
		assert(!l.test(1));
		assert(!l.test(10));
		
		for(int i = 0; i < 3; i++) {
			l.set(10);
			assert(!l.test(0));
			assert(!l.test(1));
			assert(!l.test(9));
			assert(l.test(10));
			assert(!l.test(11));
		}
		
		assert_eq(10, l.max());
		l.clear();
		assert(!l.test(10));
		assert_eq(std::numeric_limits<size_t>::max(), l.max());
		
		RandomSource rnd(12);
		size_t lim = 2000;
		for(size_t i = 0; i < lim; i++) {
			uint32_t ri = rnd.nextU32() % 10000;
			l.set(ri);
			assert(l.test(ri));
		}
	}
	cerr << "PASSED" << endl;
}

#endif /*def MAIN_SSTRING*/
