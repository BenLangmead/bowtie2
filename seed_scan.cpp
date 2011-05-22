/**
 * seed_scan.cpp
 *
 * Somewhat generic routines for searching for any of a set of up to 64-bit
 * words in a longer string of 2-bit characters.  Uses a simple 1-byte filter
 */

#include "seed_scan.h"

#ifdef MAIN_SEED_SCAN
#include <iostream>

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
		assert(make_pair((uint32_t)0, (uint32_t)0) == sc.hits()[0].id());
		assert_eq(6, sc.hits()[0].off());
		sc.nextChar(3);
		assert_eq(2, sc.hits().size());
		assert(make_pair((uint32_t)0, (uint32_t)0) == sc.hits()[0].id());
		assert_eq(6, sc.hits()[0].off());
		assert(make_pair((uint32_t)0, (uint32_t)0) == sc.hits()[1].id());
		assert_eq(7, sc.hits()[1].off());
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
		assert(make_pair((uint32_t)0, (uint32_t)0) == sc.hits()[0].id());
		assert_eq(6, sc.hits()[0].off());
		sc.nextChar(0);
		assert_eq(2, sc.hits().size());
		assert(make_pair((uint32_t)0, (uint32_t)0) == sc.hits()[0].id());
		assert_eq(6, sc.hits()[0].off());
		assert(make_pair((uint32_t)0, (uint32_t)0) == sc.hits()[1].id());
		assert_eq(7, sc.hits()[1].off());
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
	}
	cerr << "PASSED" << endl;
}
#endif /*def MAIN_SEED_SCAN*/
