/**
 * seed_scan.cpp
 *
 * Somewhat generic routines for searching for any of a set of up to 64-bit
 * words in a longer string of 2-bit characters.  Uses a simple 1-byte filter
 */

#include "seed_scan.h"

#ifdef MAIN_SEED_SCAN
#include <iostream>
#include "random_source.h"

using namespace std;

int main(void) {
	cerr << "Case 1 ... ";
	{
		SeedScanTable tab;
		tab.init(1);
		tab.add(make_pair(0, 0), 7, 2);
		SeedScanner sc;
		sc.init(tab);
		sc.nextChar(0);
		assert(sc.hits().empty());
		sc.nextChar(0);
		assert(sc.hits().empty());
		for(size_t i = 0; i < 100; i++) {
			sc.nextChar(0);
			assert(sc.hits().empty());
		}
		sc.nextChar(1);
		assert(sc.hits().empty());
		sc.nextChar(3);
		assert(!sc.hits().empty());
		assert_eq(1, sc.hits().size());
		assert_eq(0, sc.hits()[0].ns());
		for(size_t i = 0; i < 100; i++) {
			sc.nextChar(0);
			assert_eq(1, sc.hits().size());
		}
		assert_eq(1, sc.hits().size());
		sc.nextChar(1);
		assert_eq(1, sc.hits().size());
		sc.nextChar(1);
		assert_eq(1, sc.hits().size());
		sc.nextChar(1);
		assert_eq(1, sc.hits().size());
		sc.nextChar(2);
		assert_eq(1, sc.hits().size());
		sc.nextChar(3);
		assert_eq(1, sc.hits().size());
		sc.nextChar(1);
		assert_eq(1, sc.hits().size());
		sc.nextChar(3);
		assert_eq(2, sc.hits().size());
		assert_eq(0, sc.hits()[1].ns());
	}
	cerr << "PASSED" << endl;

	cerr << "Case 2 (with unmatchable) ... ";
	{
		SeedScanTable tab;
		tab.init(1);
		tab.add(make_pair(0, 0), 255, 4);
		SeedScanner sc;
		sc.init(tab);
		sc.nextChar(0);
		assert_eq(0, sc.hits().size());
		sc.nextChar(1);
		assert_eq(0, sc.hits().size());
		sc.nextChar(2);
		assert_eq(0, sc.hits().size());
		sc.nextChar(3);
		assert_eq(0, sc.hits().size());
		sc.nextChar(3);
		assert_eq(0, sc.hits().size());
		sc.nextChar(3);
		assert_eq(0, sc.hits().size());
		sc.nextChar(3);
		assert_eq(1, sc.hits().size());
		assert_eq(0, sc.hits()[0].ns());
		assert(make_pair((uint32_t)0, (uint32_t)0) == sc.hits()[0].id());
		assert_eq(3, sc.hits()[0].off());
		sc.nextChar(3);
		assert_eq(2, sc.hits().size());
		assert_eq(0, sc.hits()[1].ns());
		assert(make_pair((uint32_t)0, (uint32_t)0) == sc.hits()[0].id());
		assert_eq(3, sc.hits()[0].off());
		assert(make_pair((uint32_t)0, (uint32_t)0) == sc.hits()[1].id());
		assert_eq(4, sc.hits()[1].off());
		sc.nextChar(4);
		assert_eq(2, sc.hits().size());
		sc.nextChar(3);
		assert_eq(2, sc.hits().size());
		sc.nextChar(3);
		assert_eq(2, sc.hits().size());
		sc.nextChar(3);
		assert_eq(2, sc.hits().size());
		sc.nextChar(3);
		assert_eq(3, sc.hits().size());
		assert_eq(1, sc.hits()[2].ns());
	}
	cerr << "PASSED" << endl;

	cerr << "Case 3 (with unmatchable, As) ... ";
	{
		SeedScanTable tab;
		tab.init(1);
		tab.add(make_pair(0, 0), 0, 4);
		SeedScanner sc;
		sc.init(tab);
		sc.nextChar(3);
		assert_eq(0, sc.hits().size());
		sc.nextChar(2);
		assert_eq(0, sc.hits().size());
		sc.nextChar(1);
		assert_eq(0, sc.hits().size());
		sc.nextChar(0);
		assert_eq(0, sc.hits().size());
		sc.nextChar(0);
		assert_eq(0, sc.hits().size());
		sc.nextChar(0);
		assert_eq(0, sc.hits().size());
		sc.nextChar(0);
		assert_eq(1, sc.hits().size());
		assert_eq(0, sc.hits()[0].ns());
		assert(make_pair((uint32_t)0, (uint32_t)0) == sc.hits()[0].id());
		assert_eq(3, sc.hits()[0].off());
		sc.nextChar(0);
		assert_eq(2, sc.hits().size());
		assert_eq(0, sc.hits()[1].ns());
		assert(make_pair((uint32_t)0, (uint32_t)0) == sc.hits()[0].id());
		assert_eq(3, sc.hits()[0].off());
		assert(make_pair((uint32_t)0, (uint32_t)0) == sc.hits()[1].id());
		assert_eq(4, sc.hits()[1].off());
		sc.nextChar(4);
		assert_eq(2, sc.hits().size());
		sc.nextChar(0);
		assert_eq(2, sc.hits().size());
		sc.nextChar(0);
		assert_eq(2, sc.hits().size());
		sc.nextChar(0);
		assert_eq(2, sc.hits().size());
		sc.nextChar(0);
		assert_eq(3, sc.hits().size());
		assert_eq(1, sc.hits()[2].ns());
	}
	cerr << "PASSED" << endl;

	cerr << "Case 4 (randomness) ... ";
	{
		RandomSource rnd(2345);
		SeedScanTable tab;
		uint64_t seq1 = 0x333;
		uint64_t seq2 = 0x17f;
		tab.init(1);         // key bps
		tab.add(
			make_pair(0, 0), // id
			seq1,            // seq
			5);              // len
		tab.add(
			make_pair(0, 0), // id
			seq2,            // seq
			5);              // len
		SeedScanner sc;
		sc.init(tab);
		uint64_t seq = 0;
		size_t prevHits = 0;
		size_t ns = 0;
		for(size_t i = 0; i < 100000; i++) {
			int c = (int)(rnd.nextU32() % 5);
			if(c == 4) {
				ns++;
				seq = 0;
			} else {
				seq <<= 2;
				seq |= c;
				seq &= 0x3ff;
			}
			sc.nextChar(c);
			if(sc.hits().size() > prevHits) {
				assert_eq(sc.hits().size(), prevHits+1);
				for(size_t j = prevHits; j < sc.hits().size(); j++) {
					assert_eq(ns, sc.hits()[j].ns());
				}
				prevHits = sc.hits().size();
				assert(seq == seq1 || seq == seq2);
			}
		}
	}
	cerr << "PASSED" << endl;
}
#endif /*def MAIN_SEED_SCAN*/
