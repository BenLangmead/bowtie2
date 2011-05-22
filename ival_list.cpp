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
		list.add(10, 20);
		list.add(30, 40);
		list.add(20, 30);
		assert(!list.locusPresent(5));
		assert(!list.locusPresent(9));
		assert(list.locusPresent(10));
		assert(list.locusPresent(11));
		assert(list.locusPresent(19));
		assert(list.locusPresent(20));
		assert(list.locusPresent(21));
		assert(list.locusPresent(29));
		assert(list.locusPresent(30));
		assert(list.locusPresent(31));
		assert(list.locusPresent(39));
		assert(!list.locusPresent(40));
		assert(!list.locusPresent(41));
	}
	cerr << " PASSED" << endl;

	cerr << "Case 2 ... ";
	{
		EIvalMergeList list((size_t)5);
		list.add(10, 20);
		for(size_t i = 5; i < 45; i++) {
			assert(list.locusPresent(i) == (i >= 10 && i < 20));
		}
		list.clear();
		list.add(15, 25);
		for(size_t i = 5; i < 45; i++) {
			assert(list.locusPresent(i) == (i >= 15 && i < 25));
		}
	}
	cerr << " PASSED" << endl;

	cerr << "Case 3 ... ";
	{
		EIvalMergeList list((size_t)5);
		for(size_t i = 0; i < 20; i++) {
			list.add(10*i, 10*i+9);
		}
		for(size_t i = 0; i < 200; i++) {
			assert(list.locusPresent(i) == ((i % 10) != 9));
		}
	}
	cerr << " PASSED" << endl;

	cerr << "Case 4 ... ";
	{
		EIvalMergeList list((size_t)5);
		for(int i = 19; i >= 0; i--) {
			list.add(10*i, 10*i+9);
		}
		for(size_t i = 0; i < 200; i++) {
			assert(list.locusPresent(i) == ((i % 10) != 9));
		}
	}
	cerr << " PASSED" << endl;

	cerr << "Random testing ... ";
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
				list1.add(start, end);
				list2.add(start, end);
			}
			assert_geq(num_intervals, list1.size());
			assert_geq(num_intervals, list2.size());
			assert(list1.repOk());
			assert(list2.repOk());
			for(uint32_t i = 0; i < max_width+1; i++) {
				assert(list1.repOk());
				assert(list2.repOk());
				ASSERT_ONLY(bool l1 = list1.locusPresent(i));
				ASSERT_ONLY(bool l2 = list2.locusPresent(i));
				assert_eq(l1, l2);
			}
		}
	}
	cerr << " PASSED" << endl;
}

#endif /*def MAIN_IVAL_DS*/
