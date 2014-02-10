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

#ifndef MULTIKEY_QSORT_H_
#define MULTIKEY_QSORT_H_

#include <iostream>
#include "sequence_io.h"
#include "alphabet.h"
#include "assert_helpers.h"
#include "diff_sample.h"
#include "sstring.h"
#include "btypes.h"

using namespace std;

/**
 * Swap elements a and b in s
 */
template <typename TStr, typename TPos>
static inline void swap(TStr& s, size_t slen, TPos a, TPos b) {
	assert_lt(a, slen);
	assert_lt(b, slen);
	swap(s[a], s[b]);
}

/**
 * Swap elements a and b in array s
 */
template <typename TVal, typename TPos>
static inline void swap(TVal* s, size_t slen, TPos a, TPos b) {
	assert_lt(a, slen);
	assert_lt(b, slen);
	swap(s[a], s[b]);
}

/**
 * Helper macro for swapping elements a and b in s.  Does some additional
 * sainty checking w/r/t begin and end (which are parameters to the sorting
 * routines below).
 */
#define SWAP(s, a, b) { \
	assert_geq(a, begin); \
	assert_geq(b, begin); \
	assert_lt(a, end); \
	assert_lt(b, end); \
	swap(s, slen, a, b); \
}

/**
 * Helper macro for swapping the same pair of elements a and b in two different
 * strings s and s2.  This is a helpful variant if, for example, the caller
 * would like to see how their input was permuted by the sort routine (in that
 * case, the caller would let s2 be an array s2[] where s2 is the same length
 * as s and s2[i] = i).
 */
#define SWAP2(s, s2, a, b) { \
	SWAP(s, a, b); \
	swap(s2, slen, a, b); \
}

#define SWAP1(s, s2, a, b) { \
	SWAP(s, a, b); \
}

/**
 * Helper macro that swaps a range of elements [i, i+n) with another
 * range [j, j+n) in s.
 */
#define VECSWAP(s, i, j, n) { \
	if(n > 0) { vecswap(s, slen, i, j, n, begin, end); } \
}

/**
 * Helper macro that swaps a range of elements [i, i+n) with another
 * range [j, j+n) both in s and s2.
 */
#define VECSWAP2(s, s2, i, j, n) { \
	if(n > 0) { vecswap2(s, slen, s2, i, j, n, begin, end); } \
}

/**
 * Helper function that swaps a range of elements [i, i+n) with another
 * range [j, j+n) in s.  begin and end represent the current range under
 * consideration by the caller (one of the recursive multikey_quicksort
 * routines below).
 */
template <typename TStr, typename TPos>
static inline void vecswap(TStr& s, size_t slen, TPos i, TPos j, TPos n, TPos begin, TPos end) {
	assert_geq(i, begin);
	assert_geq(j, begin);
	assert_lt(i, end);
	assert_lt(j, end);
	while(n-- > 0) {
		assert_geq(n, 0);
		TPos a = i+n;
		TPos b = j+n;
		assert_geq(a, begin);
		assert_geq(b, begin);
		assert_lt(a, end);
		assert_lt(b, end);
		swap(s, slen, a, b);
	}
}

template <typename TVal, typename TPos>
static inline void vecswap(TVal *s, size_t slen, TPos i, TPos j, TPos n, TPos begin, TPos end) {
	assert_geq(i, begin);
	assert_geq(j, begin);
	assert_lt(i, end);
	assert_lt(j, end);
	while(n-- > 0) {
		assert_geq(n, 0);
		TPos a = i+n;
		TPos b = j+n;
		assert_geq(a, begin);
		assert_geq(b, begin);
		assert_lt(a, end);
		assert_lt(b, end);
		swap(s, slen, a, b);
	}
}

/**
 * Helper function that swaps a range of elements [i, i+n) with another range
 * [j, j+n) both in s and s2.  begin and end represent the current range under
 * consideration by the caller (one of the recursive multikey_quicksort
 * routines below).
 */
template <typename TStr, typename TPos>
static inline void vecswap2(
	TStr& s,
	size_t slen,
	TStr& s2,
	TPos i,
	TPos j,
	TPos n,
	TPos begin,
	TPos end)
{
	assert_geq(i, begin);
	assert_geq(j, begin);
	assert_lt(i, end);
	assert_lt(j, end);
	while(n-- > 0) {
		assert_geq(n, 0);
		TPos a = i+n;
		TPos b = j+n;
		assert_geq(a, begin);
		assert_geq(b, begin);
		assert_lt(a, end);
		assert_lt(b, end);
		swap(s, slen, a, b);
		swap(s2, slen, a, b);
	}
}

template <typename TVal, typename TPos>
static inline void vecswap2(TVal* s, size_t slen, TVal* s2, TPos i, TPos j, TPos n, TPos begin, TPos end) {
	assert_geq(i, begin);
	assert_geq(j, begin);
	assert_lt(i, end);
	assert_lt(j, end);
	while(n-- > 0) {
		assert_geq(n, 0);
		TPos a = i+n;
		TPos b = j+n;
		assert_geq(a, begin);
		assert_geq(b, begin);
		assert_lt(a, end);
		assert_lt(b, end);
		swap(s, slen, a, b);
		swap(s2, slen, a, b);
	}
}

/// Retrieve an int-ized version of the ath character of string s, or,
/// if a goes off the end of s, return a (user-specified) int greater
/// than any TAlphabet character - 'hi'.
#define CHAR_AT(ss, aa) ((length(s[ss]) > aa) ? (int)(s[ss][aa]) : hi)

/// Retrieve an int-ized version of the ath character of string s, or,
/// if a goes off the end of s, return a (user-specified) int greater
/// than any TAlphabet character - 'hi'.
#define CHAR_AT_SUF(si, off) \
	(((off + s[si]) < hlen) ? ((int)(host[off + s[si]])) : (hi))

/// Retrieve an int-ized version of the ath character of string s, or,
/// if a goes off the end of s, return a (user-specified) int greater
/// than any TAlphabet character - 'hi'.

#define CHAR_AT_SUF_U8(si, off) char_at_suf_u8(host, hlen, s, si, off, hi)

// Note that CHOOSE_AND_SWAP_RANDOM_PIVOT is unused
#define CHOOSE_AND_SWAP_RANDOM_PIVOT(sw, ch) {                            \
	/* Note: rand() didn't really cut it here; it seemed to run out of */ \
	/* randomness and, after a time, returned the same thing over and */  \
	/* over again */                                                      \
	a = (rand() % n) + begin; /* choose pivot between begin and end */  \
	assert_lt(a, end); assert_geq(a, begin);                              \
	sw(s, s2, begin, a); /* move pivot to beginning */                    \
}

/**
 * Ad-hoc DNA-centric way of choose a pretty good pivot without using
 * the pseudo-random number generator.  We try to get a 1 or 2 if
 * possible, since they'll split things more evenly than a 0 or 4.  We
 * also avoid swapping in the event that we choose the first element.
 */
#define CHOOSE_AND_SWAP_SMART_PIVOT(sw, ch) {                                    \
	a = begin; /* choose first elt */                                            \
	/* now try to find a better elt */                                           \
	if(n >= 5) { /* n is the difference between begin and end */                 \
		if     (ch(begin+1, depth) == 1 || ch(begin+1, depth) == 2) a = begin+1; \
		else if(ch(begin+2, depth) == 1 || ch(begin+2, depth) == 2) a = begin+2; \
		else if(ch(begin+3, depth) == 1 || ch(begin+3, depth) == 2) a = begin+3; \
		else if(ch(begin+4, depth) == 1 || ch(begin+4, depth) == 2) a = begin+4; \
		if(a != begin) sw(s, s2, begin, a); /* move pivot to beginning */        \
	}                                                                            \
	/* the element at [begin] now holds the pivot value */                       \
}

#define CHOOSE_AND_SWAP_PIVOT CHOOSE_AND_SWAP_SMART_PIVOT

#ifndef NDEBUG

/**
 * Assert that the range of chars at depth 'depth' in strings 'begin'
 * to 'end' in string-of-suffix-offsets s is parititioned properly
 * according to the ternary paritioning strategy of Bentley and McIlroy
 * (*prior to* swapping the = regions to the center)
 */
template<typename THost>
bool assertPartitionedSuf(
	const THost& host,
	TIndexOffU *s,
	size_t slen,
	int hi,
	int pivot,
	size_t begin,
	size_t end,
	size_t depth)
{
	size_t hlen = host.length();
	int state = 0; // 0 -> 1st = section, 1 -> < section, 2 -> > section, 3 -> 2nd = section
	for(size_t i = begin; i < end; i++) {
		switch(state) {
		case 0:
			if       (CHAR_AT_SUF(i, depth) < pivot)  { state = 1; break; }
			else if  (CHAR_AT_SUF(i, depth) > pivot)  { state = 2; break; }
			assert_eq(CHAR_AT_SUF(i, depth), pivot);  break;
		case 1:
			if       (CHAR_AT_SUF(i, depth) > pivot)  { state = 2; break; }
			else if  (CHAR_AT_SUF(i, depth) == pivot) { state = 3; break; }
			assert_lt(CHAR_AT_SUF(i, depth), pivot);  break;
		case 2:
			if       (CHAR_AT_SUF(i, depth) == pivot) { state = 3; break; }
			assert_gt(CHAR_AT_SUF(i, depth), pivot);	 break;
		case 3:
			assert_eq(CHAR_AT_SUF(i, depth), pivot);	 break;
		}
	}
	return true;
}

/**
 * Assert that the range of chars at depth 'depth' in strings 'begin'
 * to 'end' in string-of-suffix-offsets s is parititioned properly
 * according to the ternary paritioning strategy of Bentley and McIlroy
 * (*after* swapping the = regions to the center)
 */
template<typename THost>
bool assertPartitionedSuf2(
	const THost& host,
	TIndexOffU *s,
	size_t slen,
	int hi,
	int pivot,
	size_t begin,
	size_t end,
	size_t depth)
{
	size_t hlen = host.length();
	int state = 0; // 0 -> < section, 1 -> = section, 2 -> > section
	for(size_t i = begin; i < end; i++) {
		switch(state) {
		case 0:
			if       (CHAR_AT_SUF(i, depth) == pivot) { state = 1; break; }
			else if  (CHAR_AT_SUF(i, depth) > pivot)  { state = 2; break; }
			assert_lt(CHAR_AT_SUF(i, depth), pivot);  break;
		case 1:
			if       (CHAR_AT_SUF(i, depth) > pivot)  { state = 2; break; }
			assert_eq(CHAR_AT_SUF(i, depth), pivot);  break;
		case 2:
			assert_gt(CHAR_AT_SUF(i, depth), pivot);  break;
		}
	}
	return true;
}
#endif

/**
 * Assert that string s of suffix offsets into string 'host' is a seemingly
 * legitimate suffix-offset list (at this time, we just check that it doesn't
 * list any suffix twice).
 */
static inline void sanityCheckInputSufs(TIndexOffU *s, size_t slen) {
	assert_gt(slen, 0);
	for(size_t i = 0; i < slen; i++) {
		// Actually, it's convenient to allow the caller to provide
		// suffix offsets thare are off the end of the host string.
		// See, e.g., build() in diff_sample.cpp.
		//assert_lt(s[i], length(host));
		for(size_t j = i+1; j < slen; j++) {
			assert_neq(s[i], s[j]);
		}
	}
}

/**
 * Assert that the string s of suffix offsets into  'host' really are in
 * lexicographical order up to depth 'upto'.
 */
template <typename T>
void sanityCheckOrderedSufs(
	const T& host,
	size_t hlen,
	TIndexOffU *s,
	size_t slen,
	size_t upto,
	size_t lower = 0,
	size_t upper = OFF_MASK)
{
	assert_lt(s[0], hlen);
	upper = min<size_t>(upper, slen-1);
	for(size_t i = lower; i < upper; i++) {
		// Allow s[i+t] to point off the end of the string; this is
		// convenient for some callers
		if(s[i+1] >= hlen) continue;
#ifndef NDEBUG
		if(upto == OFF_MASK) {
			assert(sstr_suf_lt(host, s[i], hlen, host, s[i+1], hlen, false));
		} else {
			if(sstr_suf_upto_lt(host, s[i], host, s[i+1], upto, false)) {
				// operator > treats shorter strings as
				// lexicographically smaller, but we want to opposite
				//assert(isPrefix(suffix(host, s[i+1]), suffix(host, s[i])));
			}
		}
#endif
	}
}

/**
 * Main multikey quicksort function for suffixes.  Based on Bentley &
 * Sedgewick's algorithm on p.5 of their paper "Fast Algorithms for
 * Sorting and Searching Strings".  That algorithm has been extended in
 * three ways:
 *
 *  1. Deal with keys of different lengths by checking bounds and
 *     considering off-the-end values to be 'hi' (b/c our goal is the
 *     BWT transform, we're biased toward considring prefixes as
 *     lexicographically *greater* than their extensions).
 *  2. The multikey_qsort_suffixes version takes a single host string
 *     and a list of suffix offsets as input.  This reduces memory
 *     footprint compared to an approach that treats its input
 *     generically as a set of strings (not necessarily suffixes), thus
 *     requiring that we store at least two integers worth of
 *     information for each string.
 *  3. Sorting functions take an extra "upto" parameter that upper-
 *     bounds the depth to which the function sorts.
 *
 * TODO: Consult a tie-breaker (like a difference cover sample) if two
 * keys share a long prefix.
 */
template<typename T>
void mkeyQSortSuf(
	const T& host,
	size_t hlen,
	TIndexOffU *s,
	size_t slen,
	int hi,
	size_t begin,
	size_t end,
	size_t depth,
	size_t upto = OFF_MASK)
{
	// Helper for making the recursive call; sanity-checks arguments to
	// make sure that the problem actually got smaller.
	#define MQS_RECURSE_SUF(nbegin, nend, ndepth) { \
		assert(nbegin > begin || nend < end || ndepth > depth); \
		if(ndepth < upto) { /* don't exceed depth of 'upto' */ \
			mkeyQSortSuf(host, hlen, s, slen, hi, nbegin, nend, ndepth, upto); \
		} \
	}
	assert_leq(begin, slen);
	assert_leq(end, slen);
	size_t a, b, c, d, /*e,*/ r;
	size_t n = end - begin;
	if(n <= 1) return;                 // 1-element list already sorted
	CHOOSE_AND_SWAP_PIVOT(SWAP1, CHAR_AT_SUF); // pick pivot, swap it into [begin]
	int v = CHAR_AT_SUF(begin, depth); // v <- randomly-selected pivot value
	#ifndef NDEBUG
	{
		bool stillInBounds = false;
		for(size_t i = begin; i < end; i++) {
			if(depth < (hlen-s[i])) {
				stillInBounds = true;
				break;
			} else { /* already fell off this suffix */ }
		}
		assert(stillInBounds); // >=1 suffix must still be in bounds
	}
	#endif
	a = b = begin;
	c = d = end-1;
	while(true) {
		// Invariant: everything before a is = pivot, everything
		// between a and b is <
		int bc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v >= (bc = CHAR_AT_SUF(b, depth))) {
			if(v == bc) {
				SWAP(s, a, b); a++;
			}
			b++;
		}
		// Invariant: everything after d is = pivot, everything
		// between c and d is >
		int cc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v <= (cc = CHAR_AT_SUF(c, depth))) {
			if(v == cc) {
				SWAP(s, c, d); d--;
			}
			c--;
		}
		if(b > c) break;
		SWAP(s, b, c);
		b++;
		c--;
	}
	assert(a > begin || c < end-1);                      // there was at least one =s
	assert_lt(d-c, n); // they can't all have been > pivot
	assert_lt(b-a, n); // they can't all have been < pivot
	assert(assertPartitionedSuf(host, s, slen, hi, v, begin, end, depth));  // check pre-=-swap invariant
	r = min(a-begin, b-a); VECSWAP(s, begin, b-r,   r);  // swap left = to center
	r = min(d-c, end-d-1); VECSWAP(s, b,     end-r, r);  // swap right = to center
	assert(assertPartitionedSuf2(host, s, slen, hi, v, begin, end, depth)); // check post-=-swap invariant
	r = b-a; // r <- # of <'s
	if(r > 0) {
		MQS_RECURSE_SUF(begin, begin + r, depth); // recurse on <'s
	}
	// Do not recurse on ='s if the pivot was the off-the-end value;
	// they're already fully sorted
	if(v != hi) {
		MQS_RECURSE_SUF(begin + r, begin + r + (a-begin) + (end-d-1), depth+1); // recurse on ='s
	}
	r = d-c; // r <- # of >'s excluding those exhausted
	if(r > 0 && v < hi-1) {
		MQS_RECURSE_SUF(end-r, end, depth); // recurse on >'s
	}
}

/**
 * Toplevel function for multikey quicksort over suffixes.
 */
template<typename T>
void mkeyQSortSuf(
	const T& host,
	TIndexOffU *s,
	size_t slen,
	int hi,
	bool verbose = false,
	bool sanityCheck = false,
	size_t upto = OFF_MASK)
{
	size_t hlen = host.length();
	assert_gt(slen, 0);
	if(sanityCheck) sanityCheckInputSufs(s, slen);
	mkeyQSortSuf(host, hlen, s, slen, hi, (size_t)0, slen, (size_t)0, upto);
	if(sanityCheck) sanityCheckOrderedSufs(host, hlen, s, slen, upto);
}

/**
 * Just like mkeyQSortSuf but all swaps are applied to s2 as well as s.
 * This is a helpful variant if, for example, the caller would like to
 * see how their input was permuted by the sort routine (in that case,
 * the caller would let s2 be an array s2[] where s2 is the same length
 * as s and s2[i] = i).
 */
template<typename T>
void mkeyQSortSuf2(
	const T& host,
	size_t hlen,
	TIndexOffU *s,
	size_t slen,
	TIndexOffU *s2,
	int hi,
	size_t begin,
	size_t end,
	size_t depth,
	size_t upto = OFF_MASK)
{
	// Helper for making the recursive call; sanity-checks arguments to
	// make sure that the problem actually got smaller.
	#define MQS_RECURSE_SUF_DS(nbegin, nend, ndepth) { \
		assert(nbegin > begin || nend < end || ndepth > depth); \
		if(ndepth < upto) { /* don't exceed depth of 'upto' */ \
			mkeyQSortSuf2(host, hlen, s, slen, s2, hi, nbegin, nend, ndepth, upto); \
		} \
	}
	assert_leq(begin, slen);
	assert_leq(end, slen);
	size_t a, b, c, d, /*e,*/ r;
	size_t n = end - begin;
	if(n <= 1) return;                 // 1-element list already sorted
	CHOOSE_AND_SWAP_PIVOT(SWAP2, CHAR_AT_SUF); // pick pivot, swap it into [begin]
	int v = CHAR_AT_SUF(begin, depth); // v <- randomly-selected pivot value
	#ifndef NDEBUG
	{
		bool stillInBounds = false;
		for(size_t i = begin; i < end; i++) {
			if(depth < (hlen-s[i])) {
				stillInBounds = true;
				break;
			} else { /* already fell off this suffix */ }
		}
		assert(stillInBounds); // >=1 suffix must still be in bounds
	}
	#endif
	a = b = begin;
	c = d = /*e =*/ end-1;
	while(true) {
		// Invariant: everything before a is = pivot, everything
		// between a and b is <
		int bc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v >= (bc = CHAR_AT_SUF(b, depth))) {
			if(v == bc) {
				SWAP2(s, s2, a, b); a++;
			}
			b++;
		}
		// Invariant: everything after d is = pivot, everything
		// between c and d is >
		int cc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v <= (cc = CHAR_AT_SUF(c, depth))) {
			if(v == cc) {
				SWAP2(s, s2, c, d); d--; /*e--;*/
			}
			//else if(c == e && v == hi) e--;
			c--;
		}
		if(b > c) break;
		SWAP2(s, s2, b, c);
		b++;
		c--;
	}
	assert(a > begin || c < end-1);                      // there was at least one =s
	assert_lt(/*e*/d-c, n); // they can't all have been > pivot
	assert_lt(b-a, n); // they can't all have been < pivot
	assert(assertPartitionedSuf(host, s, slen, hi, v, begin, end, depth));  // check pre-=-swap invariant
	r = min(a-begin, b-a); VECSWAP2(s, s2, begin, b-r,   r);  // swap left = to center
	r = min(d-c, end-d-1); VECSWAP2(s, s2, b,     end-r, r);  // swap right = to center
	assert(assertPartitionedSuf2(host, s, slen, hi, v, begin, end, depth)); // check post-=-swap invariant
	r = b-a; // r <- # of <'s
	if(r > 0) {
		MQS_RECURSE_SUF_DS(begin, begin + r, depth); // recurse on <'s
	}
	// Do not recurse on ='s if the pivot was the off-the-end value;
	// they're already fully sorted
	if(v != hi) {
		MQS_RECURSE_SUF_DS(begin + r, begin + r + (a-begin) + (end-d-1), depth+1); // recurse on ='s
	}
	r = d-c;   // r <- # of >'s excluding those exhausted
	if(r > 0 && v < hi-1) {
		MQS_RECURSE_SUF_DS(end-r, end, depth); // recurse on >'s
	}
}

/**
 * Toplevel function for multikey quicksort over suffixes with double
 * swapping.
 */
template<typename T>
void mkeyQSortSuf2(
	const T& host,
	TIndexOffU *s,
	size_t slen,
	TIndexOffU *s2,
	int hi,
	bool verbose = false,
	bool sanityCheck = false,
	size_t upto = OFF_MASK)
{
	size_t hlen = host.length();
	if(sanityCheck) sanityCheckInputSufs(s, slen);
	TIndexOffU *sOrig = NULL;
	if(sanityCheck) {
		sOrig = new TIndexOffU[slen];
		memcpy(sOrig, s, OFF_SIZE * slen);
	}
	mkeyQSortSuf2(host, hlen, s, slen, s2, hi, (size_t)0, slen, (size_t)0, upto);
	if(sanityCheck) {
		sanityCheckOrderedSufs(host, hlen, s, slen, upto);
		for(size_t i = 0; i < slen; i++) {
			assert_eq(s[i], sOrig[s2[i]]);
		}
		delete[] sOrig;
	}
}

// Ugly but necessary; otherwise the compiler chokes dramatically on
// the DifferenceCoverSample<> template args to the next few functions
template <typename T>
class DifferenceCoverSample;

/**
 * Constant time
 */
template<typename T1, typename T2> inline
bool sufDcLt(
	const T1& host,
	const T2& s1,
	const T2& s2,
	const DifferenceCoverSample<T1>& dc,
	bool sanityCheck = false)
{
	size_t diff = dc.tieBreakOff(s1, s2);
	ASSERT_ONLY(size_t hlen = host.length());
	assert_lt(diff, dc.v());
	assert_lt(diff, hlen-s1);
	assert_lt(diff, hlen-s2);
	if(sanityCheck) {
		for(size_t i = 0; i < diff; i++) {
			assert_eq(host[s1+i], host[s2+i]);
		}
	}
	bool ret = dc.breakTie(s1+diff, s2+diff) < 0;
#ifndef NDEBUG
	if(sanityCheck && ret != sstr_suf_lt(host, s1, hlen, host, s2, hlen, false)) {
		assert(false);
	}
#endif
	return ret;
}

/**
 * k log(k)
 */
template<typename T> inline
void qsortSufDc(
	const T& host,
	size_t hlen,
	TIndexOffU* s,
	size_t slen,
	const DifferenceCoverSample<T>& dc,
	size_t begin,
	size_t end,
	bool sanityCheck = false)
{
	assert_leq(end, slen);
	assert_lt(begin, slen);
	assert_gt(end, begin);
	size_t n = end - begin;
	if(n <= 1) return;                 // 1-element list already sorted
	// Note: rand() didn't really cut it here; it seemed to run out of
	// randomness and, after a time, returned the same thing over and
	// over again
	size_t a = (rand() % n) + begin; // choose pivot between begin and end
	assert_lt(a, end);
	assert_geq(a, begin);
	SWAP(s, end-1, a); // move pivot to end
	size_t cur = 0;
	for(size_t i = begin; i < end-1; i++) {
		if(sufDcLt(host, s[i], s[end-1], dc, sanityCheck)) {
			if(sanityCheck)
				assert(dollarLt(suffix(host, s[i]), suffix(host, s[end-1])));
			assert_lt(begin + cur, end-1);
			SWAP(s, i, begin + cur);
			cur++;
		}
	}
	// Put pivot into place
	assert_lt(cur, end-begin);
	SWAP(s, end-1, begin+cur);
	if(begin+cur > begin) qsortSufDc(host, hlen, s, slen, dc, begin, begin+cur);
	if(end > begin+cur+1) qsortSufDc(host, hlen, s, slen, dc, begin+cur+1, end);
}

/**
 * Toplevel function for multikey quicksort over suffixes.
 */
template<typename T1, typename T2>
void mkeyQSortSufDcU8(
	const T1& host1,
	const T2& host,
	size_t hlen,
	TIndexOffU* s,
	size_t slen,
	const DifferenceCoverSample<T1>& dc,
	int hi,
	bool verbose = false,
	bool sanityCheck = false)
{
	if(sanityCheck) sanityCheckInputSufs(s, slen);
	mkeyQSortSufDcU8(host1, host, hlen, s, slen, dc, hi, 0, slen, 0, sanityCheck);
	if(sanityCheck) sanityCheckOrderedSufs(host1, hlen, s, slen, OFF_MASK);
}

/**
 * Return a boolean indicating whether s1 < s2 using the difference
 * cover to break the tie.
 */
template<typename T1, typename T2> inline
bool sufDcLtU8(
	const T1& host1,
	const T2& host,
	size_t hlen,
	size_t s1,
	size_t s2,
	const DifferenceCoverSample<T1>& dc,
	bool sanityCheck = false)
{
	hlen += 0;
	size_t diff = dc.tieBreakOff((TIndexOffU)s1, (TIndexOffU)s2);
	assert_lt(diff, dc.v());
	assert_lt(diff, hlen-s1);
	assert_lt(diff, hlen-s2);
	if(sanityCheck) {
		for(size_t i = 0; i < diff; i++) {
			assert_eq(host[s1+i], host1[s2+i]);
		}
	}
	bool ret = dc.breakTie((TIndexOffU)(s1+diff), (TIndexOffU)(s2+diff)) < 0;
	// Sanity-check return value using dollarLt
#ifndef NDEBUG
	bool ret2 = sstr_suf_lt(host1, s1, hlen, host, s2, hlen, false);
	assert(!sanityCheck || ret == ret2);
#endif
	return ret;
}

/**
 * k log(k)
 */
template<typename T1, typename T2> inline
void qsortSufDcU8(
	const T1& host1,
	const T2& host,
	size_t hlen,
	TIndexOffU* s,
	size_t slen,
	const DifferenceCoverSample<T1>& dc,
	size_t begin,
	size_t end,
	bool sanityCheck = false)
{
	assert_leq(end, slen);
	assert_lt(begin, slen);
	assert_gt(end, begin);
	size_t n = end - begin;
	if(n <= 1) return;                 // 1-element list already sorted
	// Note: rand() didn't really cut it here; it seemed to run out of
	// randomness and, after a time, returned the same thing over and
	// over again
	size_t a = (rand() % n) + begin; // choose pivot between begin and end
	assert_lt(a, end);
	assert_geq(a, begin);
	SWAP(s, end-1, a); // move pivot to end
	size_t cur = 0;
	for(size_t i = begin; i < end-1; i++) {
		if(sufDcLtU8(host1, host, hlen, s[i], s[end-1], dc, sanityCheck)) {
#ifndef NDEBUG
			if(sanityCheck) {
				assert(sstr_suf_lt(host1, s[i], hlen, host1, s[end-1], hlen, false));
			}
			assert_lt(begin + cur, end-1);
#endif
			SWAP(s, i, begin + cur);
			cur++;
		}
	}
	// Put pivot into place
	assert_lt(cur, end-begin);
	SWAP(s, end-1, begin+cur);
	if(begin+cur > begin) qsortSufDcU8(host1, host, hlen, s, slen, dc, begin, begin+cur);
	if(end > begin+cur+1) qsortSufDcU8(host1, host, hlen, s, slen, dc, begin+cur+1, end);
}

#define BUCKET_SORT_CUTOFF (4 * 1024 * 1024)
#define SELECTION_SORT_CUTOFF 6

// 5 64-element buckets for bucket-sorting A, C, G, T, $
extern TIndexOffU bkts[4][4 * 1024 * 1024];

/**
 * Straightforwardly obtain a uint8_t-ized version of t[off].  This
 * works fine as long as TStr is not packed.
 */
template<typename TStr>
inline uint8_t get_uint8(const TStr& t, size_t off) {
	return t[off];
}

/**
 * For incomprehensible generic-programming reasons, getting a uint8_t
 * version of a character in a packed String<> requires casting first
 * to Dna then to uint8_t.
 */
template<>
inline uint8_t get_uint8<S2bDnaString>(const S2bDnaString& t, size_t off) {
	return (uint8_t)t[off];
}

/**
 * Return character at offset 'off' from the 'si'th suffix in the array
 * 's' of suffixes.  If the character is out-of-bounds, return hi.
 */
template<typename TStr>
static inline int char_at_suf_u8(
	const TStr& host,
	size_t hlen,
	TIndexOffU* s,
	size_t si,
	size_t off,
	uint8_t hi)
{
	return ((off+s[si]) < hlen) ? get_uint8(host, off+s[si]) : (hi);
}

template<typename T1, typename T2>
static void selectionSortSufDcU8(
		const T1& host1,
		const T2& host,
        size_t hlen,
        TIndexOffU* s,
        size_t slen,
        const DifferenceCoverSample<T1>& dc,
        uint8_t hi,
        size_t begin,
        size_t end,
        size_t depth,
        bool sanityCheck = false)
{
#define ASSERT_SUF_LT(l, r) \
	if(sanityCheck && \
	   !sstr_suf_lt(host1, s[l], hlen, host1, s[r], hlen, false)) { \
		assert(false); \
	}

	assert_gt(end, begin+1);
	assert_leq(end-begin, SELECTION_SORT_CUTOFF);
	assert_eq(hi, 4);
	size_t v = dc.v();
	if(end == begin+2) {
		size_t off = dc.tieBreakOff(s[begin], s[begin+1]);
		if(off + s[begin] >= hlen ||
		   off + s[begin+1] >= hlen)
		{
			off = OFF_MASK;
		}
		if(off != OFF_MASK) {
			if(off < depth) {
				qsortSufDcU8<T1,T2>(host1, host, hlen, s, slen, dc,
				                    begin, end, sanityCheck);
				// It's helpful for debugging if we call this here
				if(sanityCheck) {
					sanityCheckOrderedSufs(host1, hlen, s, slen,
					                       OFF_MASK, begin, end);
				}
				return;
			}
			v = off - depth + 1;
		}
	}
	assert_leq(v, dc.v());
	size_t lim = v;
	assert_geq(lim, 0);
	for(size_t i = begin; i < end-1; i++) {
		size_t targ = i;
		size_t targoff = depth + s[i];
		for(size_t j = i+1; j < end; j++) {
			assert_neq(j, targ);
			size_t joff = depth + s[j];
			size_t k;
			for(k = 0; k <= lim; k++) {
				assert_neq(j, targ);
				uint8_t jc = (k + joff < hlen)    ? get_uint8(host, k + joff)    : hi;
				uint8_t tc = (k + targoff < hlen) ? get_uint8(host, k + targoff) : hi;
				assert(jc != hi || tc != hi);
				if(jc > tc) {
					// the jth suffix is greater than the current
					// smallest suffix
					ASSERT_SUF_LT(targ, j);
					break;
				} else if(jc < tc) {
					// the jth suffix is less than the current smallest
					// suffix, so update smallest to be j
					ASSERT_SUF_LT(j, targ);
					targ = j;
					targoff = joff;
					break;
				} else if(k == lim) {
					// Check whether either string ends immediately
					// after this character
					assert_leq(k + joff + 1, hlen);
					assert_leq(k + targoff + 1, hlen);
					if(k + joff + 1 == hlen) {
						// targ < j
						assert_neq(k + targoff + 1, hlen);
						ASSERT_SUF_LT(targ, j);
						break;
					} else if(k + targoff + 1 == hlen) {
						// j < targ
						ASSERT_SUF_LT(j, targ);
						targ = j;
						targoff = joff;
						break;
					}
				} else {
					// They're equal so far, keep going
				}
			}
			// The jth suffix was equal to the current smallest suffix
			// up to the difference-cover period, so disambiguate with
			// difference cover
			if(k == lim+1) {
				assert_neq(j, targ);
				if(sufDcLtU8(host1, host, hlen, s[j], s[targ], dc, sanityCheck)) {
					// j < targ
					assert(!sufDcLtU8(host1, host, hlen, s[targ], s[j], dc, sanityCheck));
					ASSERT_SUF_LT(j, targ);
					targ = j;
					targoff = joff;
				} else {
					assert(sufDcLtU8(host1, host, hlen, s[targ], s[j], dc, sanityCheck));
					ASSERT_SUF_LT(targ, j); // !
				}
			}
		}
		if(i != targ) {
			ASSERT_SUF_LT(targ, i);
			// swap i and targ
			TIndexOffU tmp = s[i];
			s[i] = s[targ];
			s[targ] = tmp;
		}
		for(size_t j = i+1; j < end; j++) {
			ASSERT_SUF_LT(i, j);
		}
	}
	if(sanityCheck) {
		sanityCheckOrderedSufs(host1, hlen, s, slen, OFF_MASK, begin, end);
	}
}

template<typename T1, typename T2>
static void bucketSortSufDcU8(
		const T1& host1,
		const T2& host,
        size_t hlen,
        TIndexOffU* s,
        size_t slen,
        const DifferenceCoverSample<T1>& dc,
        uint8_t hi,
        size_t begin,
        size_t end,
        size_t depth,
        bool sanityCheck = false)
{
	size_t cnts[] = { 0, 0, 0, 0, 0 };
	#define BKT_RECURSE_SUF_DC_U8(nbegin, nend) { \
		bucketSortSufDcU8<T1,T2>(host1, host, hlen, s, slen, dc, hi, \
		                         (nbegin), (nend), depth+1, sanityCheck); \
	}
	assert_gt(end, begin);
	assert_leq(end-begin, BUCKET_SORT_CUTOFF);
	assert_eq(hi, 4);
	if(end == begin+1) return; // 1-element list already sorted
	if(depth > dc.v()) {
		// Quicksort the remaining suffixes using difference cover
		// for constant-time comparisons; this is O(k*log(k)) where
		// k=(end-begin)
		qsortSufDcU8<T1,T2>(host1, host, hlen, s, slen, dc, begin, end, sanityCheck);
		return;
	}
	if(end-begin <= SELECTION_SORT_CUTOFF) {
		// Bucket sort remaining items
		selectionSortSufDcU8(host1, host, hlen, s, slen, dc, hi,
		                     begin, end, depth, sanityCheck);
		if(sanityCheck) {
			sanityCheckOrderedSufs(host1, hlen, s, slen,
			                       OFF_MASK, begin, end);
		}
		return;
	}
	for(size_t i = begin; i < end; i++) {
		size_t off = depth + s[i];
		uint8_t c = (off < hlen) ? get_uint8(host, off) : hi;
		assert_leq(c, 4);
		if(c == 0) {
			s[begin + cnts[0]++] = s[i];
		} else {
			bkts[c-1][cnts[c]++] = s[i];
		}
	}
	assert_eq(cnts[0] + cnts[1] + cnts[2] + cnts[3] + cnts[4], end - begin);
	size_t cur = begin + cnts[0];
	if(cnts[1] > 0) { memcpy(&s[cur], bkts[0], cnts[1] << (OFF_SIZE/4 + 1)); cur += cnts[1]; }
	if(cnts[2] > 0) { memcpy(&s[cur], bkts[1], cnts[2] << (OFF_SIZE/4 + 1)); cur += cnts[2]; }
	if(cnts[3] > 0) { memcpy(&s[cur], bkts[2], cnts[3] << (OFF_SIZE/4 + 1)); cur += cnts[3]; }
	if(cnts[4] > 0) { memcpy(&s[cur], bkts[3], cnts[4] << (OFF_SIZE/4 + 1)); }
	// This frame is now totally finished with bkts[][], so recursive
	// callees can safely clobber it; we're not done with cnts[], but
	// that's local to the stack frame.
	cur = begin;
	if(cnts[0] > 0) {
		BKT_RECURSE_SUF_DC_U8(cur, cur + cnts[0]); cur += cnts[0];
	}
	if(cnts[1] > 0) {
		BKT_RECURSE_SUF_DC_U8(cur, cur + cnts[1]); cur += cnts[1];
	}
	if(cnts[2] > 0) {
		BKT_RECURSE_SUF_DC_U8(cur, cur + cnts[2]); cur += cnts[2];
	}
	if(cnts[3] > 0) {
		BKT_RECURSE_SUF_DC_U8(cur, cur + cnts[3]);
	}
	// Done
}

/**
 * Main multikey quicksort function for suffixes.  Based on Bentley &
 * Sedgewick's algorithm on p.5 of their paper "Fast Algorithms for
 * Sorting and Searching Strings".  That algorithm has been extended in
 * three ways:
 *
 *  1. Deal with keys of different lengths by checking bounds and
 *     considering off-the-end values to be 'hi' (b/c our goal is the
 *     BWT transform, we're biased toward considring prefixes as
 *     lexicographically *greater* than their extensions).
 *  2. The multikey_qsort_suffixes version takes a single host string
 *     and a list of suffix offsets as input.  This reduces memory
 *     footprint compared to an approach that treats its input
 *     generically as a set of strings (not necessarily suffixes), thus
 *     requiring that we store at least two integers worth of
 *     information for each string.
 *  3. Sorting functions take an extra "upto" parameter that upper-
 *     bounds the depth to which the function sorts.
 */
template<typename T1, typename T2>
void mkeyQSortSufDcU8(
	const T1& host1,
	const T2& host,
	size_t hlen,
	TIndexOffU* s,
	size_t slen,
	const DifferenceCoverSample<T1>& dc,
	int hi,
	size_t begin,
	size_t end,
	size_t depth,
	bool sanityCheck = false)
{
	// Helper for making the recursive call; sanity-checks arguments to
	// make sure that the problem actually got smaller.
	#define MQS_RECURSE_SUF_DC_U8(nbegin, nend, ndepth) { \
		assert(nbegin > begin || nend < end || ndepth > depth); \
		mkeyQSortSufDcU8(host1, host, hlen, s, slen, dc, hi, nbegin, nend, ndepth, sanityCheck); \
	}
	assert_leq(begin, slen);
	assert_leq(end, slen);
	size_t n = end - begin;
	if(n <= 1) return; // 1-element list already sorted
	if(depth > dc.v()) {
		// Quicksort the remaining suffixes using difference cover
		// for constant-time comparisons; this is O(k*log(k)) where
		// k=(end-begin)
		qsortSufDcU8<T1,T2>(host1, host, hlen, s, slen, dc, begin, end, sanityCheck);
		if(sanityCheck) {
			sanityCheckOrderedSufs(host1, hlen, s, slen, OFF_MASK, begin, end);
		}
		return;
	}
	if(n <= BUCKET_SORT_CUTOFF) {
		// Bucket sort remaining items
		bucketSortSufDcU8(host1, host, hlen, s, slen, dc,
		                  (uint8_t)hi, begin, end, depth, sanityCheck);
		if(sanityCheck) {
			sanityCheckOrderedSufs(host1, hlen, s, slen, OFF_MASK, begin, end);
		}
		return;
	}
	size_t a, b, c, d, r;
	CHOOSE_AND_SWAP_PIVOT(SWAP1, CHAR_AT_SUF_U8); // choose pivot, swap to begin
	int v = CHAR_AT_SUF_U8(begin, depth); // v <- pivot value
	#ifndef NDEBUG
	{
		bool stillInBounds = false;
		for(size_t i = begin; i < end; i++) {
			if(depth < (hlen-s[i])) {
				stillInBounds = true;
				break;
			} else { /* already fell off this suffix */ }
		}
		assert(stillInBounds); // >=1 suffix must still be in bounds
	}
	#endif
	a = b = begin;
	c = d = end-1;
	while(true) {
		// Invariant: everything before a is = pivot, everything
		// between a and b is <
		int bc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v >= (bc = CHAR_AT_SUF_U8(b, depth))) {
			if(v == bc) {
				SWAP(s, a, b); a++;
			}
			b++;
		}
		// Invariant: everything after d is = pivot, everything
		// between c and d is >
		int cc = 0; // shouldn't have to init but gcc on Mac complains
		//bool hiLatch = true;
		while(b <= c && v <= (cc = CHAR_AT_SUF_U8(c, depth))) {
			if(v == cc) {
				SWAP(s, c, d); d--;
			}
			//else if(hiLatch && cc == hi) { }
			c--;
		}
		if(b > c) break;
		SWAP(s, b, c);
		b++;
		c--;
	}
	assert(a > begin || c < end-1);                      // there was at least one =s
	assert_lt(d-c, n); // they can't all have been > pivot
	assert_lt(b-a, n); // they can't all have been < pivot
	r = min(a-begin, b-a); VECSWAP(s, begin, b-r,   r);  // swap left = to center
	r = min(d-c, end-d-1); VECSWAP(s, b,     end-r, r);  // swap right = to center
	r = b-a; // r <- # of <'s
	if(r > 0) {
		MQS_RECURSE_SUF_DC_U8(begin, begin + r, depth); // recurse on <'s
	}
	// Do not recurse on ='s if the pivot was the off-the-end value;
	// they're already fully sorted
	if(v != hi) {
		MQS_RECURSE_SUF_DC_U8(begin + r, begin + r + (a-begin) + (end-d-1), depth+1); // recurse on ='s
	}
	r = d-c; // r <- # of >'s excluding those exhausted
	if(r > 0 && v < hi-1) {
		MQS_RECURSE_SUF_DC_U8(end-r, end, depth); // recurse on >'s
	}
}


#endif /*MULTIKEY_QSORT_H_*/
