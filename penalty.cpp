//
//  penalty.cpp
//

#include <iostream>
#include "penalty.h"

using namespace std;

/**
 * Init an array that maps quality to penalty according to 'type'
 * and 'cons'
 */
void Penalties::initPens(
	int *pens,  // array to fill
	int type,   // penalty type; qual | rounded qual | constant
	int cons)   // constant for when penalty type is constant
{
	if(type == COST_MODEL_ROUNDED_QUAL) {
		for(int i = 0; i < 256; i++) {
			pens[i] = qualRounds[i];
		}
	} else if(type == COST_MODEL_QUAL) {
		for(int i = 0; i < 256; i++) {
			pens[i] = i;
		}
	} else if(type == COST_MODEL_CONSTANT) {
		assert_gt(mmcost, 0);
		for(int i = 0; i < 256; i++) {
			pens[i] = cons;
		}
	} else {
		throw 1;
	}
#ifndef NDEBUG
	for(int i = 0; i < 256; i++) {
		assert_geq(pens[i], 0);
	}
#endif
}

/**
 * Given the score floor for valid alignments and the length of the read,
 * calculate the maximum possible number of read gaps that could occur in a
 * valid alignment.
 */
int Penalties::maxReadGaps(
	int64_t minsc,
	size_t rdlen) const
{
	// Score if all characters match
	int64_t sc = (int64_t)(rdlen * match);
	assert_geq(sc, minsc);
	// Now convert matches to read gaps until sc calls below minsc
	bool first = true;
	int num = 0;
	while(sc >= minsc) {
		sc -= match;
		if(first) {
			first = false;
			// Subtract both penalties
			sc -= readGapOpen();
		} else {
			// Subtract just the extension penalty
			sc -= readGapExtend();
		}
		num++;
	}
	assert_gt(num, 0);
	return num-1;
}

/**
 * Given the score floor for valid alignments and the length of the read,
 * calculate the maximum possible number of reference gaps that could occur
 * in a valid alignment.
 */
int Penalties::maxRefGaps(
	int64_t minsc,
	size_t rdlen) const
{
	// Score if all characters match
	int64_t sc = (int64_t)(rdlen * match);
	assert_geq(sc, minsc);
	// Now convert matches to read gaps until sc calls below minsc
	bool first = true;
	int num = 0;
	while(sc >= minsc) {
		if(first) {
			first = false;
			// Subtract both penalties
			sc -= refGapOpen();
		} else {
			// Subtract just the extension penalty
			sc -= refGapExtend();
		}
		num++;
	}
	assert_gt(num, 0);
	return num-1;
}

/**
 * Given a read sequence, return true iff the read passes the N filter.
 * The N filter rejects reads with more than the number of Ns calculated by
 * taking nCeilConst + nCeilLinear * read length.
 */
bool Penalties::nFilter(const BTDnaString& rd) const {
	float ns = 0.0f;
	size_t rdlen = rd.length();
	float maxns = nCeilConst + nCeilLinear * rdlen;
	if(maxns < 0.0f) maxns = 0.0f;
	assert_gt(rd.length(), 0);
	for(size_t i = 0; i < rdlen; i++) {
		if(rd[i] == 4) {
			ns += 1.0f;
			if(ns > maxns) return false; // doesn't pass
		}
	}
	return true; // passes
}

/**
 * Given a read sequence, return true iff the read passes the N filter.
 * The N filter rejects reads with more than the number of Ns calculated by
 * taking nCeilConst + nCeilLinear * read length.
 *
 * For paired-end reads, there is a	question of how to apply the filter.
 * The filter could be applied to both mates separately, which might then
 * prevent paired-end alignment.  Or the filter could be applied to the
 * reads as though they're concatenated together.  The latter approach has
 * pros and cons.  The pro is that we can use paired-end information to
 * recover alignments for mates that would not have passed the N filter on
 * their own.  The con is that we might not want to do that, since the
 * non-N portion of the bad mate might contain particularly unreliable
 * information.
 */
void Penalties::nFilterPair(
	const BTDnaString* rd1, // mate 1
	const BTDnaString* rd2, // mate 2
	bool& filt1,            // true -> mate 1 rejected by filter
	bool& filt2)            // true -> mate 2 rejected by filter
	const
{
	// Both fail to pass by default
	filt1 = filt2 = false;
	if(rd1 != NULL && rd2 != NULL && ncatpair) {
		size_t rdlen1 = rd1->length();
		size_t rdlen2 = rd2->length();
		float maxns = nCeilConst + nCeilLinear * (rdlen1 + rdlen2);
		size_t ns = 0;
		for(size_t i = 0; i < rdlen1; i++) {
			if((*rd1)[i] == 4) ns++;
			if((float)ns > maxns) {
				// doesn't pass
				return;
			}
		}
		for(size_t i = 0; i < rdlen2; i++) {
			if((*rd2)[i] == 4) ns++;
			if((float)ns > maxns) {
				// doesn't pass
				return;
			}
		}
		// Both pass
		filt1 = filt2 = true;
	} else {
		if(rd1 != NULL) filt1 = nFilter(*rd1);
		if(rd2 != NULL) filt2 = nFilter(*rd2);
	}
}

#ifdef PENALTY_MAIN

int main() {
	{
		cout << "Case 1: Simple 1 ... ";
		Penalties pen = Penalties::bwaLike();
		assert_eq(1, pen.match);
		
		assert_eq(0, pen.maxReadGaps(0, 10));  // 10 - 1 - 15 = -6
		assert_eq(0, pen.maxReadGaps(0, 11));  // 11 - 1 - 15 = -5
		assert_eq(0, pen.maxReadGaps(0, 12));  // 12 - 1 - 15 = -4
		assert_eq(0, pen.maxReadGaps(0, 13));  // 13 - 1 - 15 = -3
		assert_eq(0, pen.maxReadGaps(0, 14));  // 14 - 1 - 15 = -2
		assert_eq(0, pen.maxReadGaps(0, 15));  // 15 - 1 - 15 = -1
		assert_eq(1, pen.maxReadGaps(0, 16));  // 16 - 1 - 15 =  0
		assert_eq(1, pen.maxReadGaps(0, 17));  // 17 - 2 - 19 = -4
		assert_eq(1, pen.maxReadGaps(0, 18));  // 18 - 2 - 19 = -3
		assert_eq(1, pen.maxReadGaps(0, 19));  // 19 - 2 - 19 = -2
		assert_eq(1, pen.maxReadGaps(0, 20));  // 20 - 2 - 19 = -1
		assert_eq(2, pen.maxReadGaps(0, 21));  // 21 - 2 - 19 =  0
		
		assert_eq(0, pen.maxRefGaps(0, 10));   // 10 - 0 - 15 = -5
		assert_eq(0, pen.maxRefGaps(0, 11));   // 11 - 0 - 15 = -4
		assert_eq(0, pen.maxRefGaps(0, 12));   // 12 - 0 - 15 = -3
		assert_eq(0, pen.maxRefGaps(0, 13));   // 13 - 0 - 15 = -2
		assert_eq(0, pen.maxRefGaps(0, 14));   // 14 - 0 - 15 = -1
		assert_eq(1, pen.maxRefGaps(0, 15));   // 15 - 0 - 15 =  0
		assert_eq(1, pen.maxRefGaps(0, 16));   // 16 - 0 - 19 = -3
		assert_eq(1, pen.maxRefGaps(0, 17));   // 17 - 0 - 19 = -2
		assert_eq(1, pen.maxRefGaps(0, 18));   // 18 - 0 - 19 = -1
		assert_eq(2, pen.maxRefGaps(0, 19));   // 19 - 0 - 19 =  0
		assert_eq(2, pen.maxRefGaps(0, 20));   // 20 - 0 - 23 = -3
		assert_eq(2, pen.maxRefGaps(0, 21));   // 21 - 0 - 23 = -2
		
		// N ceiling: const=2, linear=0.1
		assert_eq(2, pen.nCeil(1));
		assert_eq(2, pen.nCeil(3));
		assert_eq(2, pen.nCeil(5));
		assert_eq(2, pen.nCeil(7));
		assert_eq(2, pen.nCeil(9));
		assert_eq(3, pen.nCeil(10));
		for(int i = 0; i < 30; i++) {
			assert_eq(3, pen.n(i));
			assert_eq(3, pen.mm(i));
		}
		assert_eq(3, pen.snp);
		cout << "PASSED" << endl;
	}
	{
		cout << "Case 2: Simple 2 ... ";
		Penalties pen(
			4,               // reward for a match
			COST_MODEL_QUAL, // how to penalize mismatches
			0,               // constant if mm pelanty is a constant
			30,              // penalty for nuc mm in decoded colorspace als
			3.0f,            // max # ref Ns allowed in alignment; const coeff
			0.4f,            // max # ref Ns allowed in alignment; linear coeff
			COST_MODEL_QUAL, // how to penalize Ns in the read
			0,               // constant if N pelanty is a constant
			true,            // whether to concatenate mates before N filtering
			25,              // constant coeff for cost of gap in the read
			25,              // constant coeff for cost of gap in the ref
			10,              // coeff of linear term for cost of gap in read
			10               // coeff of linear term for cost of gap in ref
		);

		assert_eq(4, pen.match);
		assert_eq(COST_MODEL_QUAL, pen.mmcostType);
		assert_eq(COST_MODEL_QUAL, pen.npenType);
		
		assert_eq(0, pen.maxReadGaps(0, 8));  // 32 - 4 - 35 = -7
		assert_eq(0, pen.maxReadGaps(0, 9));  // 36 - 4 - 35 = -3
		assert_eq(1, pen.maxReadGaps(0, 10)); // 40 - 4 - 35 =  1
		assert_eq(1, pen.maxReadGaps(0, 11)); // 44 - 8 - 45 = -9
		assert_eq(1, pen.maxReadGaps(0, 12)); // 48 - 8 - 45 = -5
		assert_eq(1, pen.maxReadGaps(0, 13)); // 52 - 8 - 45 = -1
		assert_eq(2, pen.maxReadGaps(0, 14)); // 56 - 8 - 45 =  3
		
		assert_eq(0, pen.maxRefGaps(0, 8));   // 32 - 0 - 35 = -3
		assert_eq(1, pen.maxRefGaps(0, 9));   // 36 - 0 - 35 =  1
		assert_eq(1, pen.maxRefGaps(0, 10));  // 40 - 0 - 45 = -5
		assert_eq(1, pen.maxRefGaps(0, 11));  // 44 - 0 - 45 = -1
		assert_eq(2, pen.maxRefGaps(0, 12));  // 48 - 0 - 45 =  3
		assert_eq(2, pen.maxRefGaps(0, 13));  // 52 - 0 - 55 = -3
		assert_eq(3, pen.maxRefGaps(0, 14));  // 56 - 0 - 55 =  1

		// N ceiling: const=3, linear=0.4
		assert_eq(3, pen.nCeil(1));
		assert_eq(3, pen.nCeil(2));
		assert_eq(4, pen.nCeil(3));
		assert_eq(4, pen.nCeil(4));
		assert_eq(5, pen.nCeil(5));
		assert_eq(5, pen.nCeil(6));
		assert_eq(5, pen.nCeil(7));
		assert_eq(6, pen.nCeil(8));
		assert_eq(6, pen.nCeil(9));

		for(int i = 0; i < 256; i++) {
			assert_eq(i, pen.n(i));
			assert_eq(i, pen.mm(i));
		}
		assert_eq(30, pen.snp);
		cout << "PASSED" << endl;
	}
}

#endif /*def PENALTY_MAIN*/
