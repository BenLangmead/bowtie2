/*
 * random_test.cpp
 *
 * Pseudo-random tests for Bowite2.  Desiderata for the pseudo-random
 * tester were:
 *
 * 1. It should test the pseudo-random number generator itself.
 * 2. It should avoid the overhead of invoking the bowtie2 executables
 *    repeatedly, by invoking the bowtie() function directly.
 * 3. Where sets of options should yield the same results, it should
 *    confirm that the results are the same
 * 4. 
 */

#include <iostream>
#include <string.h>
#include "random_source.h"

using namespace std;

int main(void) {
	RandomSource rand;
	rand.init(0);
	uint32_t ts[32];
	memset(ts, 0, 32*sizeof(uint32_t));
	uint32_t r = rand.nextU32();
	cout << "Without reseeding:" << endl;
	for(int i = 0; i < 10000; i++) {
		uint32_t nr = rand.nextU32();
		for(int j = 0; j < 32; j++) {
			if(((r >> j) & 1) != ((nr >> j) & 1)) {
				ts[j]++;
			}
		}
	}
	for(int j = 0; j < 32; j++) {
		cout << ts[j] << endl;
	}
	memset(ts, 0, 32*sizeof(uint32_t));
	rand.init(0);
	r = rand.nextU32();
	cout << "With reseeding:" << endl;
	for(int i = 0; i < 10000; i++) {
		rand.init(i+1);
		uint32_t nr = rand.nextU32();
		for(int j = 0; j < 32; j++) {
			if(((r >> j) & 1) != ((nr >> j) & 1)) {
				ts[j]++;
			}
		}
	}
	for(int j = 0; j < 32; j++) {
		cout << ts[j] << endl;
	}
}
