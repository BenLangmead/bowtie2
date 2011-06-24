/**
 * ival_list.h
 */

#include "ival_list.h"

#ifdef MAIN_IVAL_DS

#include <iostream>
#include "random_source.h"

using namespace std;

int main(void) {
	cerr << "Case 1 ... ";
	{
		EIvalMergeList list((size_t)5);
		list.add(Interval(0, 10, true, 10));
		list.add(Interval(0, 30, true, 10));
		list.add(Interval(0, 20, true, 10));
		assert(!list.locusPresent(Coord(0, 5, true)));
		assert(!list.locusPresent(Coord(0, 9, true)));
		assert(list.locusPresent(Coord(0, 10, true)));
		assert(list.locusPresent(Coord(0, 11, true)));
		assert(list.locusPresent(Coord(0, 19, true)));
		assert(list.locusPresent(Coord(0, 20, true)));
		assert(list.locusPresent(Coord(0, 21, true)));
		assert(list.locusPresent(Coord(0, 29, true)));
		assert(list.locusPresent(Coord(0, 30, true)));
		assert(list.locusPresent(Coord(0, 31, true)));
		assert(list.locusPresent(Coord(0, 39, true)));
		assert(!list.locusPresent(Coord(0, 40, true)));
		assert(!list.locusPresent(Coord(0, 41, true)));
	}
	cerr << " PASSED" << endl;

	cerr << "Case 2 ... ";
	{
		EIvalMergeList list((size_t)5);
		list.add(Interval(0, 10, true, 10));
		for(size_t i = 5; i < 45; i++) {
			assert(list.locusPresent(Coord(0, i, true)) == (i >= 10 && i < 20));
		}
		list.clear();
		list.add(Interval(0, 15, true, 10));
		for(size_t i = 5; i < 45; i++) {
			assert(list.locusPresent(Coord(0, i, true)) == (i >= 15 && i < 25));
		}
	}
	cerr << " PASSED" << endl;

	cerr << "Case 3 ... ";
	{
		EIvalMergeList list((size_t)5);
		for(size_t i = 0; i < 20; i++) {
			list.add(Interval(0, 10*i, true, 9));
		}
		for(size_t i = 0; i < 200; i++) {
			assert(list.locusPresent(Coord(0, i, true)) == ((i % 10) != 9));
			assert(!list.locusPresent(Coord(0, i, false)));
			assert(!list.locusPresent(Coord(1, i, true)));
		}
	}
	cerr << " PASSED" << endl;

	cerr << "Case 4 ... ";
	{
		EIvalMergeList list((size_t)5);
		for(int i = 19; i >= 0; i--) {
			list.add(Interval(0, 10*i, true, 9));
		}
		for(size_t i = 0; i < 200; i++) {
			assert(list.locusPresent(Coord(0, i, true)) == ((i % 10) != 9));
			assert(!list.locusPresent(Coord(0, i, false)));
			assert(!list.locusPresent(Coord(1, i, true)));
		}
	}
	cerr << " PASSED" << endl;

	cerr << "Random testing (1 ref) ... ";
	{
		RandomSource rnd(34523);
		for(size_t c = 0; c < 10; c++) {
			EIvalMergeList list1((size_t)16);
			EIvalMergeList list2((size_t)2000);
			size_t num_intervals = 20;
			uint32_t max_width = 100;
			for(size_t i = 0; i < num_intervals; i++) {
				uint32_t start = rnd.nextU32() % max_width/2;
				uint32_t end = (rnd.nextU32() % (max_width - start - 1) + start)+1;
				assert_lt(end, max_width);
				assert_gt(end, start);
				list1.add(Interval(0, start, false, end-start));
				list2.add(Interval(0, start, false, end-start));
			}
			assert_geq(num_intervals, list1.size());
			assert_geq(num_intervals, list2.size());
			assert(list1.repOk());
			assert(list2.repOk());
			for(uint32_t i = 0; i < max_width+1; i++) {
				assert(list1.repOk());
				assert(list2.repOk());
				ASSERT_ONLY(bool l1 = list1.locusPresent(Coord(0, i, true)));
				ASSERT_ONLY(bool l2 = list2.locusPresent(Coord(0, i, true)));
				assert_eq(l1, l2);
			}
		}
	}
	cerr << " PASSED" << endl;

	cerr << "Random testing (few refs) ... ";
	{
		RandomSource rnd(34523);
		for(size_t c = 0; c < 10; c++) {
			EIvalMergeList list1((size_t)16);
			EIvalMergeList list2((size_t)2000);
			size_t num_intervals = 20;
			uint32_t max_width = 100;
			for(size_t i = 0; i < num_intervals; i++) {
				uint32_t start = rnd.nextU32() % max_width/2;
				uint32_t end = (rnd.nextU32() % (max_width - start - 1) + start)+1;
				assert_lt(end, max_width);
				assert_gt(end, start);
				bool orient = (rnd.nextU2() == 0);
				TRefId ref = (TRefId)(rnd.nextU32() % 5);
				list1.add(Interval(ref, start, orient, end-start));
				list2.add(Interval(ref, start, orient, end-start));
			}
			assert_geq(num_intervals, list1.size());
			assert_geq(num_intervals, list2.size());
			assert(list1.repOk());
			assert(list2.repOk());
			for(uint32_t i = 0; i < max_width+1; i++) {
				assert(list1.repOk());
				assert(list2.repOk());
				for(int fwi = 0; fwi < 2; fwi++) {
					bool fw = (fwi == 0);
					for(TRefId refi = 0; refi < 5; refi++) {
						ASSERT_ONLY(bool l1 = list1.locusPresent(Coord(refi, i, fw)));
						ASSERT_ONLY(bool l2 = list2.locusPresent(Coord(refi, i, fw)));
						assert_eq(l1, l2);
					}
				}
			}
		}
	}
	cerr << " PASSED" << endl;
}

#endif /*def MAIN_IVAL_DS*/
