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

#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "bt2_idx.h"

using namespace std;

///////////////////////////////////////////////////////////////////////
//
// Functions for searching Ebwts
//
///////////////////////////////////////////////////////////////////////

/**
 * Return the final character in row i (i.e. the i'th character in the
 * BWT transform).  Note that the 'L' in the name of the function
 * stands for 'last', as in the literature.
 */
int Ebwt::rowL(const SideLocus& l) const {
	// Extract and return appropriate bit-pair
#ifdef SIXTY4_FORMAT
	return (((uint64_t*)l.side(this->ebwt()))[l._by >> 3] >> ((((l._by & 7) << 2) + l._bp) << 1)) & 3;
#else
	return unpack_2b_from_8b(l.side(this->ebwt())[l._by], l._bp);
#endif
}

/**
 * Return the final character in row i (i.e. the i'th character in the
 * BWT transform).  Note that the 'L' in the name of the function
 * stands for 'last', as in the literature.
 */
int Ebwt::rowL(uint32_t i) const {
	// Extract and return appropriate bit-pair
	SideLocus l;
	l.initFromRow(i, _eh, ebwt());
	return rowL(l);
}

/**
 * Inline-function version of the above.  This does not always seem to
 * be inlined
 */
#if 0
// Use gcc's intrinsic popcountll.  I don't recommend it because it
// seems to be somewhat slower than the bit-bashing pop64 routine both
// on an AMD server and on an Intel workstation.  On the other hand,
// perhaps when the builtin is used GCC is smart enough to insert a
// pop-count instruction on architectures that have one (e.g. Itanium).
// For now, it's disabled.
#define pop64(x) __builtin_popcountll(x)
#elif 0
__declspec naked int __stdcall pop64
(uint64_t v)
{
	static const uint64_t C55 = 0x5555555555555555ll;
	static const uint64_t C33 = 0x3333333333333333ll;
	static const uint64_t C0F = 0x0F0F0F0F0F0F0F0Fll;
	__asm {
		MOVD      MM0, [ESP+4] ;v_low
		PUNPCKLDQ MM0, [ESP+8] ;v
		MOVQ      MM1, MM0     ;v
		PSRLD     MM0, 1       ;v >> 1
		PAND      MM0, [C55]   ;(v >> 1) & 0x55555555
		PSUBD     MM1, MM0     ;w = v - ((v >> 1) & 0x55555555)
		MOVQ      MM0, MM1     ;w
		PSRLD     MM1, 2       ;w >> 2
		PAND      MM0, [C33]   ;w & 0x33333333
		PAND      MM1, [C33]   ;(w >> 2)  & 0x33333333
		PADDD     MM0, MM1     ;x = (w & 0x33333333) +
		; ((w >> 2) & 0x33333333)
		MOVQ      MM1, MM0     ;x
		PSRLD     MM0, 4       ;x >> 4
		PADDD     MM0, MM1     ;x + (x >> 4)
		PAND      MM0, [C0F]   ;y = (x + (x >> 4) & 0x0F0F0F0F)
		PXOR      MM1, MM1     ;0
		PSADBW    (MM0, MM1)   ;sum across all 8 bytes
		MOVD      EAX, MM0     ;result in EAX per calling
		; convention
		EMMS                   ;clear MMX state
		RET  8                 ;pop 8-byte argument off stack
		; and return
	}
}
#elif 0
// Use a bytewise LUT version of popcount.  This is slower than the
// bit-bashing pop64 routine both on an AMD server and on an Intel
// workstation.  It seems to be about the same speed as the GCC builtin
// on Intel, and a bit faster than it on AMD.  For now, it's disabled.
const int popcntU8Table[256] = {
    0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};

// Use this bytewise population count table
inline static int pop64(uint64_t x) {
	const unsigned char * p = (const unsigned char *) &x;
	return popcntU8Table[p[0]] +
	popcntU8Table[p[1]] +
	popcntU8Table[p[2]] +
	popcntU8Table[p[3]] +
	popcntU8Table[p[4]] +
	popcntU8Table[p[5]] +
	popcntU8Table[p[6]] +
	popcntU8Table[p[7]];
}
#else
// Use this standard bit-bashing population count
inline static int pop64(uint64_t x) {
	// Lots of cache misses on following lines (>10K)
	x = x - ((x >> 1) & 0x5555555555555555llu);
	x = (x & 0x3333333333333333llu) + ((x >> 2) & 0x3333333333333333llu);
	x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0Fllu;
	x = x + (x >> 8);
	x = x + (x >> 16);
	x = x + (x >> 32);
	return (int)(x & 0x3Fllu);
}
#endif

/**
 * Tricky-bit-bashing bitpair counting for given two-bit value (0-3)
 * within a 64-bit argument.
 */
inline static int countInU64(int c, uint64_t dw) {
	uint64_t dwA  = dw &  0xAAAAAAAAAAAAAAAAllu;
	uint64_t dwNA = dw & ~0xAAAAAAAAAAAAAAAAllu;
	uint64_t tmp;
	switch(c) {
		case 0:
			tmp = (dwA >> 1) | dwNA;
			break;
		case 1:
			tmp = ~(dwA >> 1) & dwNA;
			break;
		case 2:
			tmp = (dwA >> 1) & ~dwNA;
			break;
		case 3:
			tmp = (dwA >> 1) & dwNA;
			break;
		default:
			throw;
	}
	tmp = pop64(tmp); // Gets 7.62% in profile
	if(c == 0) {
		tmp = 32 - tmp;
	}
	assert_leq(tmp, 32);
	assert_geq(tmp, 0);
	return (int)tmp;
}

/**
 * Tricky-bit-bashing bitpair counting for given two-bit value (0-3)
 * within a 64-bit argument.
 *
 * Function gets 2.32% in profile
 */
inline static void countInU64Ex(uint64_t dw, uint32_t* arrs) {
	// Cache misses here (~9K)
	uint64_t dwA  = dw &  0xAAAAAAAAAAAAAAAAllu;
	uint64_t dwNA = dw & ~0xAAAAAAAAAAAAAAAAllu;
	arrs[0] += (32 - pop64((dwA >> 1) | dwNA));
	arrs[1] += pop64(~(dwA >> 1) & dwNA);
	arrs[2] += pop64((dwA >> 1) & ~dwNA);
	arrs[3] += pop64((dwA >> 1) & dwNA);
}

#define WITHIN_FCHR(x) \
	assert_leq(x[0], this->fchr()[1]); \
	assert_leq(x[1], this->fchr()[2]); \
	assert_leq(x[2], this->fchr()[3]); \
	assert_leq(x[3], this->fchr()[4])

#define WITHIN_FCHR_DOLLARA(x) \
	assert_leq(x[0], this->fchr()[1]+1); \
	assert_leq(x[1], this->fchr()[2]); \
	assert_leq(x[2], this->fchr()[3]); \
	assert_leq(x[3], this->fchr()[4])

/**
 * Given top and bot loci, calculate counts of all four DNA chars up to
 * those loci.  Also, update a set of tops and bots for the reverse
 * index/direction using the idea from the bi-directional BWT paper.
 */
void Ebwt::mapBiLFEx(
	const SideLocus& ltop,
	const SideLocus& lbot,
	uint32_t *tops,
	uint32_t *bots,
	uint32_t *topsP, // topsP[0] = top
	uint32_t *botsP
	ASSERT_ONLY(, bool overrideSanity)
	) const
{
	// TODO: Where there's overlap, reuse the count for the overlapping
	// portion
#ifdef EBWT_STATS
	const_cast<Ebwt*>(this)->mapBiLFEx_++;
#endif
#ifndef NDEBUG
	for(int i = 0; i < 4; i++) {
		assert_eq(0, tops[0]);  assert_eq(0, bots[0]);
	}
#endif
#ifdef BOWTIE2
	countBt2SideEx(ltop, tops);
	countBt2SideEx(lbot, bots);
#else
	if(ltop.fw_) countFwSideEx(ltop, tops); // Forward side
	else         countBwSideEx(ltop, tops); // Backward side
	if(lbot.fw_) countFwSideEx(lbot, bots); // Forward side
	else         countBwSideEx(lbot, bots); // Backward side
#endif
#ifndef NDEBUG
	if(_sanity && !overrideSanity) {
		// Make sure results match up with individual calls to mapLF;
		// be sure to override sanity-checking in the callee, or we'll
		// have infinite recursion
		assert_eq(mapLF(ltop, 0, true), tops[0]);
		assert_eq(mapLF(ltop, 1, true), tops[1]);
		assert_eq(mapLF(ltop, 2, true), tops[2]);
		assert_eq(mapLF(ltop, 3, true), tops[3]);
		assert_eq(mapLF(lbot, 0, true), bots[0]);
		assert_eq(mapLF(lbot, 1, true), bots[1]);
		assert_eq(mapLF(lbot, 2, true), bots[2]);
		assert_eq(mapLF(lbot, 3, true), bots[3]);
	}
#endif
	// bots[0..3] - tops[0..3] = # of ways to extend the suffix with an
	// A, C, G, T
	botsP[0] = topsP[0] + (bots[0] - tops[0]);
	topsP[1] = botsP[0];
	botsP[1] = topsP[1] + (bots[1] - tops[1]);
	topsP[2] = botsP[1];
	botsP[2] = topsP[2] + (bots[2] - tops[2]);
	topsP[3] = botsP[2];
	botsP[3] = topsP[3] + (bots[3] - tops[3]);
}

/**
 * Given top and bot rows, calculate counts of all four DNA chars up to
 * those loci.
 */
void Ebwt::mapLFEx(
	uint32_t top,
	uint32_t bot,
	uint32_t *tops,
	uint32_t *bots
	ASSERT_ONLY(, bool overrideSanity)
	) const
{
	SideLocus ltop, lbot;
	SideLocus::initFromTopBot(top, bot, _eh, ebwt(), ltop, lbot);
	mapLFEx(ltop, lbot, tops, bots ASSERT_ONLY(, overrideSanity));
}

/**
 * Given top and bot loci, calculate counts of all four DNA chars up to
 * those loci.  Used for more advanced backtracking-search.
 */
void Ebwt::mapLFEx(
	const SideLocus& ltop,
	const SideLocus& lbot,
	uint32_t *tops,
	uint32_t *bots
	ASSERT_ONLY(, bool overrideSanity)
	) const
{
	assert(ltop.repOk(*this));
	assert(lbot.repOk(*this));
#ifdef EBWT_STATS
	const_cast<Ebwt*>(this)->mapLFExs_++;
#endif
	assert_eq(0, tops[0]); assert_eq(0, bots[0]);
	assert_eq(0, tops[1]); assert_eq(0, bots[1]);
	assert_eq(0, tops[2]); assert_eq(0, bots[2]);
	assert_eq(0, tops[3]); assert_eq(0, bots[3]);
#ifdef BOWTIE2
	countBt2SideEx(ltop, tops);
	countBt2SideEx(lbot, bots);
#endif
#ifndef NDEBUG
	if(_sanity && !overrideSanity) {
		// Make sure results match up with individual calls to mapLF;
		// be sure to override sanity-checking in the callee, or we'll
		// have infinite recursion
		assert_eq(mapLF(ltop, 0, true), tops[0]);
		assert_eq(mapLF(ltop, 1, true), tops[1]);
		assert_eq(mapLF(ltop, 2, true), tops[2]);
		assert_eq(mapLF(ltop, 3, true), tops[3]);
		assert_eq(mapLF(lbot, 0, true), bots[0]);
		assert_eq(mapLF(lbot, 1, true), bots[1]);
		assert_eq(mapLF(lbot, 2, true), bots[2]);
		assert_eq(mapLF(lbot, 3, true), bots[3]);
	}
#endif
}

/**
 * Given top and bot loci, calculate counts of all four DNA chars up to
 * those loci.  Used for more advanced backtracking-search.
 */
void Ebwt::mapLFRange(
	SideLocus& ltop,
	SideLocus& lbot,
	uint32_t num,        // Number of elts
	uint32_t* cntsUpto,  // A/C/G/T counts up to top
	uint32_t* cntsIn,    // A/C/G/T counts within range
	EList<bool> *masks
	ASSERT_ONLY(, bool overrideSanity)
	) const
{
	assert(ltop.repOk(*this));
	assert(lbot.repOk(*this));
	assert_eq(num, lbot.toBWRow() - ltop.toBWRow());
#ifdef EBWT_STATS
	const_cast<Ebwt*>(this)->mapLFRanges_++;
#endif
	assert_eq(0, cntsUpto[0]); assert_eq(0, cntsIn[0]);
	assert_eq(0, cntsUpto[1]); assert_eq(0, cntsIn[1]);
	assert_eq(0, cntsUpto[2]); assert_eq(0, cntsIn[2]);
	assert_eq(0, cntsUpto[3]); assert_eq(0, cntsIn[3]);
#ifdef BOWTIE2
	countBt2SideRange(ltop, num, cntsUpto, cntsIn, masks);
#else
	if(ltop.fw_) countFwSideRange(ltop, num, cntsUpto, cntsIn, masks); // Forward side
	else         countBwSideRange(ltop, num, cntsUpto, cntsIn, masks); // Backward side
#endif
	assert_eq(num, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
#ifndef NDEBUG
	if(_sanity && !overrideSanity) {
		// Make sure results match up with individual calls to mapLF;
		// be sure to override sanity-checking in the callee, or we'll
		// have infinite recursion
		uint32_t tops[4] = {0, 0, 0, 0};
		uint32_t bots[4] = {0, 0, 0, 0};
		assert(ltop.repOk(*this));
		assert(lbot.repOk(*this));
		mapLFEx(ltop, lbot, tops, bots);
		for(int i = 0; i < 4; i++) {
			assert(cntsUpto[i] == tops[i] || tops[i] == bots[i]);
			if(i == 0) {
				assert(cntsIn[i] == bots[i]-tops[i] ||
				       cntsIn[i] == bots[i]-tops[i]+1);
			} else {
				assert_eq(cntsIn[i], bots[i]-tops[i]);
			}
		}
	}
#endif
}

#ifndef NDEBUG
/**
 * Given top and bot loci, calculate counts of all four DNA chars up to
 * those loci.  Used for more advanced backtracking-search.
 */
void Ebwt::mapLFEx(const SideLocus& l,
                   uint32_t *arrs
                   ASSERT_ONLY(, bool overrideSanity)
                   ) const
{
	assert_eq(0, arrs[0]);
	assert_eq(0, arrs[1]);
	assert_eq(0, arrs[2]);
	assert_eq(0, arrs[3]);
#ifdef BOWTIE2
	countBt2SideEx(l, arrs);
#else
	if(l.fw_) countFwSideEx(l, arrs); // Forward side
	else      countBwSideEx(l, arrs); // Backward side
#endif
	if(_sanity && !overrideSanity) {
		// Make sure results match up with individual calls to mapLF;
		// be sure to override sanity-checking in the callee, or we'll
		// have infinite recursion
		assert_eq(mapLF(l, 0, true), arrs[0]);
		assert_eq(mapLF(l, 1, true), arrs[1]);
		assert_eq(mapLF(l, 2, true), arrs[2]);
		assert_eq(mapLF(l, 3, true), arrs[3]);
	}
}
#endif

/**
 * Given row i, return the row that the LF mapping maps i to.
 */
uint32_t Ebwt::mapLF(const SideLocus& l
                     ASSERT_ONLY(, bool overrideSanity)
                     ) const
{
#ifdef EBWT_STATS
	const_cast<Ebwt*>(this)->mapLFs_++;
#endif
	ASSERT_ONLY(uint32_t srcrow = l.toBWRow());
	uint32_t ret;
	assert(l.side(this->ebwt()) != NULL);
	int c = rowL(l);
	assert_lt(c, 4);
	assert_geq(c, 0);
#ifdef BOWTIE2
	ret = countBt2Side(l, c);
#else
	if(l.fw_) ret = countFwSide(l, c); // Forward side
	else      ret = countBwSide(l, c); // Backward side
#endif
	assert_lt(ret, this->_eh._bwtLen);
	assert_neq(srcrow, ret);
#ifndef NDEBUG
	if(_sanity && !overrideSanity) {
		// Make sure results match up with results from mapLFEx;
		// be sure to override sanity-checking in the callee, or we'll
		// have infinite recursion
		uint32_t arrs[] = { 0, 0, 0, 0 };
		mapLFEx(l, arrs, true);
		assert_eq(arrs[c], ret);
	}
#endif
	return ret;
}

/**
 * Given row i and character c, return the row that the LF mapping maps
 * i to on character c.
 */
uint32_t Ebwt::mapLF(const SideLocus& l, int c
                     ASSERT_ONLY(, bool overrideSanity)
                     ) const
{
#ifdef EBWT_STATS
	const_cast<Ebwt*>(this)->mapLFcs_++;
#endif
	uint32_t ret;
	assert_lt(c, 4);
	assert_geq(c, 0);
#ifdef BOWTIE2
	ret = countBt2Side(l, c);
#else
	if(l.fw_) ret = countFwSide(l, c); // Forward side
	else      ret = countBwSide(l, c); // Backward side
#endif
	assert_lt(ret, this->_eh._bwtLen);
#ifndef NDEBUG
	if(_sanity && !overrideSanity) {
		// Make sure results match up with results from mapLFEx;
		// be sure to override sanity-checking in the callee, or we'll
		// have infinite recursion
		uint32_t arrs[] = { 0, 0, 0, 0 };
		mapLFEx(l, arrs, true);
		assert_eq(arrs[c], ret);
	}
#endif
	return ret;
}

/**
 * Given row and its locus information, proceed on the given character
 * and return the next row, or all-fs if we can't proceed on that
 * character.  Returns 0xffffffff if this row ends in $.
 */
uint32_t Ebwt::mapLF1(
	uint32_t row,       // starting row
	const SideLocus& l, // locus for starting row
	int c               // character to proceed on
    ASSERT_ONLY(, bool overrideSanity)) const
{
#ifdef EBWT_STATS
	const_cast<Ebwt*>(this)->mapLF1cs_++;
#endif
	if(rowL(l) != c || row == _zOff) return 0xffffffff;
	uint32_t ret;
	assert_lt(c, 4);
	assert_geq(c, 0);
#ifdef BOWTIE2
	ret = countBt2Side(l, c);
#else
	if(l.fw_) ret = countFwSide(l, c); // Forward side
	else      ret = countBwSide(l, c); // Backward side
#endif
	assert_lt(ret, this->_eh._bwtLen);
#ifndef NDEBUG
	if(_sanity && !overrideSanity) {
		// Make sure results match up with results from mapLFEx;
		// be sure to override sanity-checking in the callee, or we'll
		// have infinite recursion
		uint32_t arrs[] = { 0, 0, 0, 0 };
		mapLFEx(l, arrs, true);
		assert_eq(arrs[c], ret);
	}
#endif
	return ret;
}

/**
 * Given row and its locus information, set the row to LF(row) and
 * return the character that was in the final column.
 */
int Ebwt::mapLF1(
	uint32_t& row,      // starting row
	const SideLocus& l  // locus for starting row
    ASSERT_ONLY(, bool overrideSanity)) const
{
#ifdef EBWT_STATS
	const_cast<Ebwt*>(this)->mapLF1s_++;
#endif
	if(row == _zOff) return -1;
	int c = rowL(l);
	assert_lt(c, 4);
	assert_geq(c, 0);
#ifdef BOWTIE2
	row = countBt2Side(l, c);
#else
	if(l.fw_) row = countFwSide(l, c); // Forward side
	else      row = countBwSide(l, c); // Backward side
#endif
	assert_lt(row, this->_eh._bwtLen);
#ifndef NDEBUG
	if(_sanity && !overrideSanity) {
		// Make sure results match up with results from mapLFEx;
		// be sure to override sanity-checking in the callee, or we'll
		// have infinite recursion
		uint32_t arrs[] = { 0, 0, 0, 0 };
		mapLFEx(l, arrs, true);
		assert_eq(arrs[c], row);
	}
#endif
	return c;
}

/**
 * Take an offset into the joined text and translate it into the
 * reference of the index it falls on, the offset into the reference,
 * and the length of the reference.  Use a binary search through the
 * sorted list of reference fragment ranges t
 */
void Ebwt::joinedToTextOff(
	uint32_t qlen,
	uint32_t off,
	uint32_t& tidx,
    uint32_t& textoff,
    uint32_t& tlen,
	bool rejectStraddle,
	bool& straddled) const
{
	assert(rstarts() != NULL); // must have loaded rstarts
	uint32_t top = 0;
	uint32_t bot = _nFrag; // 1 greater than largest addressable element
	uint32_t elt = 0xffffffff;
	// Begin binary search
	while(true) {
		ASSERT_ONLY(uint32_t oldelt = elt);
		elt = top + ((bot - top) >> 1);
		assert_neq(oldelt, elt); // must have made progress
		uint32_t lower = rstarts()[elt*3];
		uint32_t upper;
		if(elt == _nFrag-1) {
			upper = _eh._len;
		} else {
			upper = rstarts()[((elt+1)*3)];
		}
		assert_gt(upper, lower);
		uint32_t fraglen = upper - lower;
		if(lower <= off) {
			if(upper > off) { // not last element, but it's within
				// off is in this range; check if it falls off
				if(off + qlen > upper) {
					straddled = true;
					if(rejectStraddle) {
						// it falls off; signal no-go and return
						tidx = 0xffffffff;
						assert_lt(elt, _nFrag-1);
						return;
					}
				}
				// This is the correct text idx whether the index is
				// forward or reverse
				tidx = rstarts()[(elt*3)+1];
				assert_lt(tidx, this->_nPat);
				assert_leq(fraglen, this->plen()[tidx]);
				// it doesn't fall off; now calculate textoff.
				// Initially it's the number of characters that precede
				// the alignment in the fragment
				uint32_t fragoff = off - rstarts()[(elt*3)];
				if(!this->fw_) {
					fragoff = fraglen - fragoff - 1;
					fragoff -= (qlen-1);
				}
				// Add the alignment's offset into the fragment
				// ('fragoff') to the fragment's offset within the text
				textoff = fragoff + rstarts()[(elt*3)+2];
				assert_lt(textoff, this->plen()[tidx]);
				break; // done with binary search
			} else {
				// 'off' belongs somewhere in the region between elt
				// and bot
				top = elt;
			}
		} else {
			// 'off' belongs somewhere in the region between top and
			// elt
			bot = elt;
		}
		// continue with binary search
	}
	tlen = this->plen()[tidx];
}

/**
 * Walk 'steps' steps to the left and return the row arrived at.  If we
 * walk through the dollar sign, return 0xffffffff.
 */
uint32_t Ebwt::walkLeft(uint32_t row, uint32_t steps) const {
	assert(offs() != NULL);
	assert_neq(0xffffffff, row);
	SideLocus l;
	if(steps > 0) l.initFromRow(row, _eh, ebwt());
	while(steps > 0) {
		if(row == _zOff) return 0xffffffff;
		uint32_t newrow = this->mapLF(l);
		assert_neq(0xffffffff, newrow);
		assert_neq(newrow, row);
		row = newrow;
		steps--;
		if(steps > 0) l.initFromRow(row, _eh, ebwt());
	}
	return row;
}

/**
 * Resolve the reference offset of the BW element 'elt'.
 */
uint32_t Ebwt::getOffset(uint32_t row) const {
	assert(offs() != NULL);
	assert_neq(0xffffffff, row);
	if(row == _zOff) return 0;
	if((row & _eh._offMask) == row) return this->offs()[row >> _eh._offRate];
	int jumps = 0;
	SideLocus l;
	l.initFromRow(row, _eh, ebwt());
	while(true) {
		uint32_t newrow = this->mapLF(l);
		jumps++;
		assert_neq(0xffffffff, newrow);
		assert_neq(newrow, row);
		row = newrow;
		if(row == _zOff) {
			return jumps;
		} else if((row & _eh._offMask) == row) {
			return jumps + this->offs()[row >> _eh._offRate];
		}
		l.initFromRow(row, _eh, ebwt());
	}
}

/**
 * Resolve the reference offset of the BW element 'elt' such that
 * the offset returned is at the right-hand side of the forward
 * reference substring involved in the hit.
 */
uint32_t Ebwt::getOffset(
	uint32_t elt,
	bool fw,
	uint32_t hitlen) const
{
	uint32_t off = getOffset(elt);
	assert_neq(0xffffffff, off);
	if(!fw) {
		assert_lt(off, _eh._len);
		off = _eh._len - off - 1;
		assert_geq(off, hitlen-1);
		off -= (hitlen-1);
		assert_lt(off, _eh._len);
	}
	return off;
}

#define WITHIN_BWT_LEN(x) \
	assert_leq(x[0], this->_eh._sideBwtLen); \
	assert_leq(x[1], this->_eh._sideBwtLen); \
	assert_leq(x[2], this->_eh._sideBwtLen); \
	assert_leq(x[3], this->_eh._sideBwtLen)

#ifdef BOWTIE2
/**
 * Count all occurrences of character c from the beginning of the
 * forward side to <by,bp> and add in the occ[] count up to the side
 * break just prior to the side.
 *
 * A Bowtie 2 side is shaped like:
 *
 * XXXXXXXXXXXXXXXX [A] [C] [G] [T]
 * --------48------ -4- -4- -4- -4-  (numbers in bytes)
 */
inline uint32_t Ebwt::countBt2Side(const SideLocus& l, int c) const {
	assert_range(0, 3, c);
	assert_range(0, (int)this->_eh._sideBwtSz-1, (int)l._by);
	assert_range(0, 3, (int)l._bp);
	const uint8_t *side = l.side(this->ebwt());
	uint32_t cCnt = countUpTo(l, c);
	assert_leq(cCnt, this->_eh._sideBwtLen);
	if(c == 0 && l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff) {
		// Adjust for the fact that we represented $ with an 'A', but
		// shouldn't count it as an 'A' here
		if((l._sideByteOff + l._by > _zEbwtByteOff) ||
		   (l._sideByteOff + l._by == _zEbwtByteOff && l._bp > _zEbwtBpOff))
		{
			cCnt--; // Adjust for '$' looking like an 'A'
		}
	}
	uint32_t ret;
	// Now factor in the occ[] count at the side break
	const uint8_t *acgt8 = side + _eh._sideBwtSz;
	const uint32_t *acgt = reinterpret_cast<const uint32_t*>(acgt8);
	assert_leq(acgt[0], this->_eh._numSides * this->_eh._sideBwtLen); // b/c it's used as padding
	assert_leq(acgt[1], this->_eh._len);
	assert_leq(acgt[2], this->_eh._len);
	assert_leq(acgt[3], this->_eh._len);
	ret = acgt[c] + cCnt + this->fchr()[c];
#ifndef NDEBUG
	assert_leq(ret, this->fchr()[c+1]); // can't have jumpded into next char's section
	if(c == 0) {
		assert_leq(cCnt, this->_eh._sideBwtLen);
	} else {
		assert_leq(ret, this->_eh._bwtLen);
	}
#endif
	return ret;
}
#else
/**
 * Count all occurrences of character c from the beginning of the
 * forward side to <by,bp> and add in the occ[] count up to the side
 * break just prior to the side.
 *
 * A forward side is shaped like:
 *
 * [A] [C] XXXXXXXXXXXXXXXX
 * -4- -4- --------56------ (numbers in bytes)
 *         ^
 *         Side ptr (result from SideLocus.side())
 *
 * And following it is a reverse side shaped like:
 * 
 * [G] [T] XXXXXXXXXXXXXXXX
 * -4- -4- --------56------ (numbers in bytes)
 *         ^
 *         Side ptr (result from SideLocus.side())
 */
inline uint32_t Ebwt::countFwSide(const SideLocus& l, int c) const {
	assert_range(0, 3, c);
	assert_range(0, (int)this->_eh._sideBwtSz-1, (int)l._by);
	assert_range(0, 3, (int)l._bp);
	const uint8_t *side = l.side(this->ebwt());
	uint32_t cCnt = countUpTo(l, c);
	assert_leq(cCnt, this->_eh._sideBwtLen);
	if(c == 0 && l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff) {
		// Adjust for the fact that we represented $ with an 'A', but
		// shouldn't count it as an 'A' here
		if((l._sideByteOff + l._by > _zEbwtByteOff) ||
		   (l._sideByteOff + l._by == _zEbwtByteOff && l._bp > _zEbwtBpOff))
		{
			cCnt--; // Adjust for '$' looking like an 'A'
		}
	}
	uint32_t ret;
	// Now factor in the occ[] count at the side break
	if(c < 2) {
		const uint8_t *ac8 = side - 8;
		const uint32_t *ac = reinterpret_cast<const uint32_t*>(ac8);
		assert_leq(ac[0], this->_eh._numSides * this->_eh._sideBwtLen); // b/c it's used as padding
		assert_leq(ac[1], this->_eh._len);
		ret = ac[c] + cCnt + this->fchr()[c];
	} else {
		const uint8_t *gt8 = side + this->_eh._sideSz - 8;
		const uint32_t *gt = reinterpret_cast<const uint32_t*>(gt8); // next
		assert_leq(gt[0], this->_eh._len); assert_leq(gt[1], this->_eh._len);
		ret = gt[c-2] + cCnt + this->fchr()[c];
	}
#ifndef NDEBUG
	assert_leq(ret, this->fchr()[c+1]); // can't have jumpded into next char's section
	if(c == 0) {
		assert_leq(cCnt, this->_eh._sideBwtLen);
	} else {
		assert_leq(ret, this->_eh._bwtLen);
	}
#endif
	return ret;
}

/**
 * Count all instances of character c from <by,bp> to the logical end
 * (actual beginning) of the backward side, and subtract that from the
 * occ[] count up to the side break.
 *
 * A reverse side is shaped like:
 *
 * [G] [T] XXXXXXXXXXXXXXXX
 * -4- -4- --------56------ (numbers in bytes)
 *         ^
 *         Side ptr (result from SideLocus.side())
 *
 * And following it is a reverse side shaped like:
 * 
 * [G] [T] XXXXXXXXXXXXXXXX
 * -4- -4- --------56------ (numbers in bytes)
 *         ^
 *         Side ptr (result from SideLocus.side())
 */
inline uint32_t Ebwt::countBwSide(const SideLocus& l, int c) const {
	assert_range(0, 3, c);
	assert_range(0, (int)this->_eh._sideBwtSz-1, (int)l._by);
	assert_range(0, 3, (int)l._bp);
	const uint8_t *side = l.side(this->ebwt());
	uint32_t cCnt = countUpTo(l, c);
	if(rowL(l) == c) cCnt++;
	assert_leq(cCnt, this->_eh._sideBwtLen);
	if(c == 0 && l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff) {
		// Adjust for the fact that we represented $ with an 'A', but
		// shouldn't count it as an 'A' here
		if((l._sideByteOff + l._by > _zEbwtByteOff) ||
		   (l._sideByteOff + l._by == _zEbwtByteOff && l._bp >= _zEbwtBpOff))
		{
			cCnt--;
		}
	}
	uint32_t ret;
	// Now factor in the occ[] count at the side break
	if(c < 2) {
		const uint32_t *ac = reinterpret_cast<const uint32_t*>(side + this->_eh._sideSz - 8);
		assert_leq(ac[0], this->_eh._numSides * this->_eh._sideBwtLen); // b/c it's used as padding
		assert_leq(ac[1], this->_eh._len);
		ret = ac[c] - cCnt + this->fchr()[c];
	} else {
		const uint32_t *gt = reinterpret_cast<const uint32_t*>(side + (2*this->_eh._sideSz) - 8); // next
		assert_leq(gt[0], this->_eh._len); assert_leq(gt[1], this->_eh._len);
		ret = gt[c-2] - cCnt + this->fchr()[c];
	}
#ifndef NDEBUG
	assert_leq(ret, this->fchr()[c+1]); // can't have jumped into next char's section
	if(c == 0) {
		assert_leq(cCnt, this->_eh._sideBwtLen);
	} else {
		assert_lt(ret, this->_eh._bwtLen);
	}
#endif
	return ret;
}
#endif

#ifdef BOWTIE2
/**
 * Count all occurrences of all four nucleotides up to the starting
 * point (which must be in a forward side) given by 'l' storing the
 * result in 'cntsUpto', then count nucleotide occurrences within the
 * range of length 'num' storing the result in 'cntsIn'.  Also, keep
 * track of the characters occurring within the range by setting
 * 'masks' accordingly (masks[1][10] == true -> 11th character is a
 * 'C', and masks[0][10] == masks[2][10] == masks[3][10] == false.
 */
inline void Ebwt::countBt2SideRange(
	SideLocus& l,        // top locus
	uint32_t num,        // number of elts in range to tall
	uint32_t* cntsUpto,  // A/C/G/T counts up to top
	uint32_t* cntsIn,    // A/C/G/T counts within range
	EList<bool> *masks) const // masks indicating which range elts = A/C/G/T
{
	assert_gt(num, 0);
	assert_range(0, (int)this->_eh._sideBwtSz-1, (int)l._by);
	assert_range(0, 3, (int)l._bp);
	countUpToEx(l, cntsUpto);
	WITHIN_FCHR_DOLLARA(cntsUpto);
	WITHIN_BWT_LEN(cntsUpto);
	const uint8_t *side = l.side(this->ebwt());
	bool adjustedForDollar = false;
	if(l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff) {
		// Adjust for the fact that we represented $ with an 'A', but
		// shouldn't count it as an 'A' here
		if((l._sideByteOff + l._by > _zEbwtByteOff) ||
		   (l._sideByteOff + l._by == _zEbwtByteOff && l._bp > _zEbwtBpOff))
		{
			adjustedForDollar = true;
			cntsUpto[0]--; // Adjust for '$' looking like an 'A'
		}
	}
	// Now factor in the occ[] count at the side break
	const uint32_t *acgt = reinterpret_cast<const uint32_t*>(side + _eh._sideBwtSz);
	assert_leq(acgt[0], this->fchr()[1] + this->_eh.sideBwtLen());
	assert_leq(acgt[1], this->fchr()[2]-this->fchr()[1]);
	assert_leq(acgt[2], this->fchr()[3]-this->fchr()[2]);
	assert_leq(acgt[3], this->fchr()[4]-this->fchr()[3]);
	assert_leq(acgt[0], this->_eh._len + this->_eh.sideBwtLen());
	assert_leq(acgt[1], this->_eh._len);
	assert_leq(acgt[2], this->_eh._len);
	assert_leq(acgt[3], this->_eh._len);
	cntsUpto[0] += (acgt[0] + this->fchr()[0]);
	cntsUpto[1] += (acgt[1] + this->fchr()[1]);
	cntsUpto[2] += (acgt[2] + this->fchr()[2]);
	cntsUpto[3] += (acgt[3] + this->fchr()[3]);
	masks[0].resize(num);
	masks[1].resize(num);
	masks[2].resize(num);
	masks[3].resize(num);
	WITHIN_FCHR_DOLLARA(cntsUpto);
	WITHIN_FCHR_DOLLARA(cntsIn);
	// 'cntsUpto' is complete now.
	// Walk forward until we've tallied the entire 'In' range
	uint32_t nm = 0;
	// Rest of this side
	nm += countBt2SideRange2(l, true, num - nm, cntsIn, masks, nm);
	assert_eq(nm, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
	assert_leq(nm, num);
	SideLocus lcopy = l;
	while(nm < num) {
		// Subsequent sides, if necessary
		lcopy.nextSide(this->_eh);
		nm += countBt2SideRange2(lcopy, false, num - nm, cntsIn, masks, nm);
		WITHIN_FCHR_DOLLARA(cntsIn);
		assert_leq(nm, num);
		assert_eq(nm, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
	}
	assert_eq(num, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
	WITHIN_FCHR_DOLLARA(cntsIn);
}
#else
/**
 * Count all occurrences of all four nucleotides up to the starting
 * point (which must be in a forward side) given by 'l' storing the
 * result in 'cntsUpto', then count nucleotide occurrences within the
 * range of length 'num' storing the result in 'cntsIn'.  Also, keep
 * track of the characters occurring within the range by setting
 * 'masks' accordingly (masks[1][10] == true -> 11th character is a
 * 'C', and masks[0][10] == masks[2][10] == masks[3][10] == false.
 */
inline void Ebwt::countFwSideRange(
	SideLocus& l,        // top locus
	uint32_t num,        // number of elts in range to tall
	uint32_t* cntsUpto,  // A/C/G/T counts up to top
	uint32_t* cntsIn,    // A/C/G/T counts within range
	EList<bool> *masks) const // masks indicating which range elts = A/C/G/T
{
	assert_gt(num, 0);
	assert_range(0, (int)this->_eh._sideBwtSz-1, (int)l._by);
	assert_range(0, 3, (int)l._bp);
	countUpToEx(l, cntsUpto);
	WITHIN_FCHR_DOLLARA(cntsUpto);
	WITHIN_BWT_LEN(cntsUpto);
	const uint8_t *side = l.side(this->ebwt());
	bool adjustedForDollar = false;
	if(l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff) {
		// Adjust for the fact that we represented $ with an 'A', but
		// shouldn't count it as an 'A' here
		if((l._sideByteOff + l._by > _zEbwtByteOff) ||
		   (l._sideByteOff + l._by == _zEbwtByteOff && l._bp > _zEbwtBpOff))
		{
			adjustedForDollar = true;
			cntsUpto[0]--; // Adjust for '$' looking like an 'A'
		}
	}
	// Now factor in the occ[] count at the side break
	const uint32_t *ac = reinterpret_cast<const uint32_t*>(side - 8);
	const uint32_t *gt = reinterpret_cast<const uint32_t*>(side + this->_eh._sideSz - 8);
	assert_leq(ac[0], this->fchr()[1] + this->_eh.sideBwtLen());
	assert_leq(ac[1], this->fchr()[2]-this->fchr()[1]);
	assert_leq(gt[0], this->fchr()[3]-this->fchr()[2]);
	assert_leq(gt[1], this->fchr()[4]-this->fchr()[3]);
	assert_leq(ac[0], this->_eh._len + this->_eh.sideBwtLen()); assert_leq(ac[1], this->_eh._len);
	assert_leq(gt[0], this->_eh._len); assert_leq(gt[1], this->_eh._len);
	cntsUpto[0] += (ac[0] + this->fchr()[0]);
	cntsUpto[1] += (ac[1] + this->fchr()[1]);
	cntsUpto[2] += (gt[0] + this->fchr()[2]);
	cntsUpto[3] += (gt[1] + this->fchr()[3]);
	masks[0].resize(num);
	masks[1].resize(num);
	masks[2].resize(num);
	masks[3].resize(num);
	WITHIN_FCHR_DOLLARA(cntsUpto);
	WITHIN_FCHR_DOLLARA(cntsIn);
	// 'cntsUpto' is complete now.
	// Walk forward until we've tallied the entire 'In' range
	uint32_t nm = 0;
	// Rest of this side
	nm += countFwSideRange2(l, true, num - nm, cntsIn, masks, nm);
	assert_eq(nm, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
	assert_leq(nm, num);
	SideLocus lcopy = l;
	while(nm < num) {
		// Subsequent sides, if necessary
		lcopy.nextSide(this->_eh);
		nm += countBwSideRange2(lcopy, false, num - nm, cntsIn, masks, nm);
		WITHIN_FCHR_DOLLARA(cntsIn);
		assert_leq(nm, num);
		assert_eq(nm, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
		if(nm == num) break;
		lcopy.nextSide(this->_eh);
		nm += countFwSideRange2(lcopy, false, num - nm, cntsIn, masks, nm);
		WITHIN_FCHR_DOLLARA(cntsIn);
		assert_leq(nm, num);
		assert_eq(nm, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
	}
	assert_eq(num, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
	WITHIN_FCHR_DOLLARA(cntsIn);
}

/**
 * Count all occurrences of all four nucleotides up to the starting
 * point (which must be in a forward side) given by 'l' storing the
 * result in 'cntsUpto', then count nucleotide occurrences within the
 * range of length 'num' storing the result in 'cntsIn'.  Also, keep
 * track of the characters occurring within the range by setting
 * 'masks' accordingly (masks[1][10] == true -> 11th character is a
 * 'C', and masks[0][10] == masks[2][10] == masks[3][10] == false.
 */
inline void Ebwt::countBwSideRange(
	SideLocus& l,        // top locus
	uint32_t num,        // number of elts in range to tall
	uint32_t* cntsUpto,  // A/C/G/T counts up to top
	uint32_t* cntsIn,    // A/C/G/T counts within range
	EList<bool> *masks) const // masks indicating which range elts = A/C/G/T
{
	assert_gt(num, 0);
	assert_range(0, (int)this->_eh._sideBwtSz-1, (int)l._by);
	assert_range(0, 3, (int)l._bp);
	const uint8_t *side = l.side(this->ebwt());
	countUpToEx(l, cntsUpto);
	cntsUpto[rowL(l)]++;
	WITHIN_BWT_LEN(cntsUpto);
	if(l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff) {
		// Adjust for the fact that we represented $ with an 'A', but
		// shouldn't count it as an 'A' here
		if((l._sideByteOff + l._by > _zEbwtByteOff) ||
		   (l._sideByteOff + l._by == _zEbwtByteOff && l._bp >= _zEbwtBpOff))
		{
			cntsUpto[0]--; // Adjust for '$' looking like an 'A'
		}
	}
	// Now factor in the occ[] count at the side break
	const uint32_t *ac = reinterpret_cast<const uint32_t*>(side + this->_eh._sideSz - 8);
	const uint32_t *gt = reinterpret_cast<const uint32_t*>(side + (2*this->_eh._sideSz) - 8);
	assert_leq(ac[0], this->fchr()[1] + this->_eh.sideBwtLen());
	assert_leq(ac[1], this->fchr()[2]-this->fchr()[1]);
	assert_leq(gt[0], this->fchr()[3]-this->fchr()[2]);
	assert_leq(gt[1], this->fchr()[4]-this->fchr()[3]);
	assert_leq(ac[0], this->_eh._len + this->_eh.sideBwtLen()); assert_leq(ac[1], this->_eh._len);
	assert_leq(gt[0], this->_eh._len); assert_leq(gt[1], this->_eh._len);
	cntsUpto[0] = (ac[0] - cntsUpto[0] + this->fchr()[0]);
	cntsUpto[1] = (ac[1] - cntsUpto[1] + this->fchr()[1]);
	cntsUpto[2] = (gt[0] - cntsUpto[2] + this->fchr()[2]);
	cntsUpto[3] = (gt[1] - cntsUpto[3] + this->fchr()[3]);
	masks[0].resize(num);
	masks[1].resize(num);
	masks[2].resize(num);
	masks[3].resize(num);
	WITHIN_FCHR_DOLLARA(cntsUpto);
	WITHIN_FCHR_DOLLARA(cntsIn);
	// 'cntsUpto' is complete now.
	// Walk forward until we've tallied the entire range
	uint32_t nm = 0;
	// Rest of this side
	nm += countBwSideRange2(l, true, num - nm, cntsIn, masks, nm);
	assert_leq(nm, num);
	assert_eq(nm, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
	SideLocus lcopy = l;
	while(nm < num) {
		// Subsequent sides, if necessary
		lcopy.nextSide(this->_eh);
		nm += countFwSideRange2(lcopy, false, num - nm, cntsIn, masks, nm);
		WITHIN_FCHR_DOLLARA(cntsIn);
		assert_leq(nm, num);
		assert_eq(nm, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
		if(nm == num) break;
		lcopy.nextSide(this->_eh);
		nm += countBwSideRange2(lcopy, false, num - nm, cntsIn, masks, nm);
		WITHIN_FCHR_DOLLARA(cntsIn);
		assert_leq(nm, num);
		assert_eq(nm, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
	}
	assert_eq(num, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
	WITHIN_FCHR_DOLLARA(cntsIn);
}
#endif

#ifdef BOWTIE2
/**
 * Count all occurrences of character c from the beginning of the
 * forward side to <by,bp> and add in the occ[] count up to the side
 * break just prior to the side.
 *
 * A forward side is shaped like:
 *
 * [A] [C] XXXXXXXXXXXXXXXX
 * -4- -4- --------56------ (numbers in bytes)
 *         ^
 *         Side ptr (result from SideLocus.side())
 *
 * And following it is a reverse side shaped like:
 * 
 * [G] [T] XXXXXXXXXXXXXXXX
 * -4- -4- --------56------ (numbers in bytes)
 *         ^
 *         Side ptr (result from SideLocus.side())
 *
 */
inline void Ebwt::countBt2SideEx(const SideLocus& l, uint32_t* arrs) const
{
	assert_range(0, (int)this->_eh._sideBwtSz-1, (int)l._by);
	assert_range(0, 3, (int)l._bp);
	countUpToEx(l, arrs);
	if(l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff) {
		// Adjust for the fact that we represented $ with an 'A', but
		// shouldn't count it as an 'A' here
		if((l._sideByteOff + l._by > _zEbwtByteOff) ||
		   (l._sideByteOff + l._by == _zEbwtByteOff && l._bp > _zEbwtBpOff))
		{
			arrs[0]--; // Adjust for '$' looking like an 'A'
		}
	}
	WITHIN_FCHR(arrs);
	WITHIN_BWT_LEN(arrs);
	// Now factor in the occ[] count at the side break
	const uint8_t *side = l.side(this->ebwt());
	const uint8_t *acgt16 = side + this->_eh._sideSz - 16;
	const uint32_t *acgt = reinterpret_cast<const uint32_t*>(acgt16);
	assert_leq(acgt[0], this->fchr()[1] + this->_eh.sideBwtLen());
	assert_leq(acgt[1], this->fchr()[2]-this->fchr()[1]);
	assert_leq(acgt[2], this->fchr()[3]-this->fchr()[2]);
	assert_leq(acgt[3], this->fchr()[4]-this->fchr()[3]);
	assert_leq(acgt[0], this->_eh._len + this->_eh.sideBwtLen());
	assert_leq(acgt[1], this->_eh._len);
	assert_leq(acgt[2], this->_eh._len);
	assert_leq(acgt[3], this->_eh._len);
	arrs[0] += (acgt[0] + this->fchr()[0]);
	arrs[1] += (acgt[1] + this->fchr()[1]);
	arrs[2] += (acgt[2] + this->fchr()[2]);
	arrs[3] += (acgt[3] + this->fchr()[3]);
	WITHIN_FCHR(arrs);
}
#else
/**
 * Count all occurrences of character c from the beginning of the
 * forward side to <by,bp> and add in the occ[] count up to the side
 * break just prior to the side.
 *
 * A forward side is shaped like:
 *
 * [A] [C] XXXXXXXXXXXXXXXX
 * -4- -4- --------56------ (numbers in bytes)
 *         ^
 *         Side ptr (result from SideLocus.side())
 *
 * And following it is a reverse side shaped like:
 * 
 * [G] [T] XXXXXXXXXXXXXXXX
 * -4- -4- --------56------ (numbers in bytes)
 *         ^
 *         Side ptr (result from SideLocus.side())
 *
 */
inline void Ebwt::countFwSideEx(const SideLocus& l, uint32_t* arrs) const
{
	assert_range(0, (int)this->_eh._sideBwtSz-1, (int)l._by);
	assert_range(0, 3, (int)l._bp);
	countUpToEx(l, arrs);
	WITHIN_FCHR(arrs);
	WITHIN_BWT_LEN(arrs);
	const uint8_t *side = l.side(this->ebwt());
	if(l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff) {
		// Adjust for the fact that we represented $ with an 'A', but
		// shouldn't count it as an 'A' here
		if((l._sideByteOff + l._by > _zEbwtByteOff) ||
		   (l._sideByteOff + l._by == _zEbwtByteOff && l._bp > _zEbwtBpOff))
		{
			arrs[0]--; // Adjust for '$' looking like an 'A'
		}
	}
	// Now factor in the occ[] count at the side break
	const uint8_t *ac8 = side - 8;
	const uint8_t *gt8 = side + this->_eh._sideSz - 8;
	const uint32_t *ac = reinterpret_cast<const uint32_t*>(ac8);
	const uint32_t *gt = reinterpret_cast<const uint32_t*>(gt8);
	assert_leq(ac[0], this->fchr()[1] + this->_eh.sideBwtLen());
	assert_leq(ac[1], this->fchr()[2]-this->fchr()[1]);
	assert_leq(gt[0], this->fchr()[3]-this->fchr()[2]);
	assert_leq(gt[1], this->fchr()[4]-this->fchr()[3]);
	assert_leq(ac[0], this->_eh._len + this->_eh.sideBwtLen()); assert_leq(ac[1], this->_eh._len);
	assert_leq(gt[0], this->_eh._len); assert_leq(gt[1], this->_eh._len);
	arrs[0] += (ac[0] + this->fchr()[0]);
	arrs[1] += (ac[1] + this->fchr()[1]);
	arrs[2] += (gt[0] + this->fchr()[2]);
	arrs[3] += (gt[1] + this->fchr()[3]);
	WITHIN_FCHR(arrs);
}

/**
 * Count all instances of character c from <by,bp> to the logical end
 * (actual beginning) of the backward side, and subtract that from the
 * occ[] count up to the side break.
 */
inline void Ebwt::countBwSideEx(const SideLocus& l, uint32_t* arrs) const {
	assert_range(0, (int)this->_eh._sideBwtSz-1, (int)l._by);
	assert_range(0, 3, (int)l._bp);
	const uint8_t *side = l.side(this->ebwt());
	countUpToEx(l, arrs);
	arrs[rowL(l)]++;
	WITHIN_BWT_LEN(arrs);
	if(l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff) {
		// Adjust for the fact that we represented $ with an 'A', but
		// shouldn't count it as an 'A' here
		if((l._sideByteOff + l._by > _zEbwtByteOff) ||
		   (l._sideByteOff + l._by == _zEbwtByteOff && l._bp >= _zEbwtBpOff))
		{
			arrs[0]--; // Adjust for '$' looking like an 'A'
		}
	}
	// Now factor in the occ[] count at the side break
	const uint32_t *ac = reinterpret_cast<const uint32_t*>(side + this->_eh._sideSz - 8);
	const uint32_t *gt = reinterpret_cast<const uint32_t*>(side + (2*this->_eh._sideSz) - 8);
	assert_leq(ac[0], this->fchr()[1] + this->_eh.sideBwtLen());
	assert_leq(ac[1], this->fchr()[2]-this->fchr()[1]);
	assert_leq(gt[0], this->fchr()[3]-this->fchr()[2]);
	assert_leq(gt[1], this->fchr()[4]-this->fchr()[3]);
	assert_leq(ac[0], this->_eh._len + this->_eh.sideBwtLen()); assert_leq(ac[1], this->_eh._len);
	assert_leq(gt[0], this->_eh._len); assert_leq(gt[1], this->_eh._len);
	arrs[0] = (ac[0] - arrs[0] + this->fchr()[0]);
	arrs[1] = (ac[1] - arrs[1] + this->fchr()[1]);
	arrs[2] = (gt[0] - arrs[2] + this->fchr()[2]);
	arrs[3] = (gt[1] - arrs[3] + this->fchr()[3]);
	WITHIN_FCHR(arrs);
}
#endif

/**
 * Counts the number of occurrences of character 'c' in the given Ebwt
 * side up to (but not including) the given byte/bitpair (by/bp).
 *
 * This is a performance-critical function.  This is the top search-
 * related hit in the time profile.
 *
 * Function gets 11.09% in profile
 */
inline uint32_t Ebwt::countUpTo(const SideLocus& l, int c) const {
	// Count occurrences of c in each 64-bit (using bit trickery);
	// Someday countInU64() and pop() functions should be
	// vectorized/SSE-ized in case that helps.
	uint32_t cCnt = 0;
	const uint8_t *side = l.side(this->ebwt());
	int i = 0;
#if 1
	for(; i + 7 < l._by; i += 8) {
		cCnt += countInU64(c, *(uint64_t*)&side[i]);
	}
#else
	for(; i + 2 < l._by; i += 2) {
		cCnt += cCntLUT_16b_4[c][*(uint16_t*)&side[i]];
	}
#endif
#ifdef SIXTY4_FORMAT
	// Calculate number of bit pairs to shift off the end
	const int bpShiftoff = 32 - (((l._by & 7) << 2) + l._bp);
	if(bpShiftoff < 32) {
		assert_lt(bpShiftoff, 32);
		const uint64_t sw = (*(uint64_t*)&side[i]) << (bpShiftoff << 1);
		cCnt += countInU64(c, sw);
		if(c == 0) cCnt -= bpShiftoff; // we turned these into As
	}
#else
	// Count occurences of c in the rest of the side (using LUT)
	for(; i < l._by; i++) {
		cCnt += cCntLUT_4[0][c][side[i]];
	}
	// Count occurences of c in the rest of the byte
	if(l._bp > 0) {
		cCnt += cCntLUT_4[(int)l._bp][c][side[i]];
	}
#endif
	return cCnt;
}

/**
 * Counts the number of occurrences of all four nucleotides in the
 * given side up to (but not including) the given byte/bitpair (by/bp).
 * Count for 'a' goes in arrs[0], 'c' in arrs[1], etc.
 */
inline void Ebwt::countUpToEx(const SideLocus& l, uint32_t* arrs) const {
	int i = 0;
	// Count occurrences of each nucleotide in each 64-bit word using
	// bit trickery; note: this seems does not seem to lend a
	// significant boost to performance in practice.  If you comment
	// out this whole loop (which won't affect correctness - it will
	// just cause the following loop to take up the slack) then runtime
	// does not change noticeably. Someday the countInU64() and pop()
	// functions should be vectorized/SSE-ized in case that helps.
	const uint8_t *side = l.side(this->ebwt());
	for(; i+7 < l._by; i += 8) {
		countInU64Ex(*(uint64_t*)&side[i], arrs);
	}
#ifdef SIXTY4_FORMAT
	// Calculate number of bit pairs to shift off the end
	const int bpShiftoff = 32 - (((l._by & 7) << 2) + l._bp);
	assert_leq(bpShiftoff, 32);
	if(bpShiftoff < 32) {
		const uint64_t sw = (*(uint64_t*)&l.side(this->ebwt())[i]) << (bpShiftoff << 1);
		countInU64Ex(sw, arrs);
		arrs[0] -= bpShiftoff;
	}
#else
	// Count occurences of nucleotides in the rest of the side (using LUT)
	// Many cache misses on following lines (~20K)
	for(; i < l._by; i++) {
		arrs[0] += cCntLUT_4[0][0][side[i]];
		arrs[1] += cCntLUT_4[0][1][side[i]];
		arrs[2] += cCntLUT_4[0][2][side[i]];
		arrs[3] += cCntLUT_4[0][3][side[i]];
	}
	// Count occurences of c in the rest of the byte
	if(l._bp > 0) {
		arrs[0] += cCntLUT_4[(int)l._bp][0][side[i]];
		arrs[1] += cCntLUT_4[(int)l._bp][1][side[i]];
		arrs[2] += cCntLUT_4[(int)l._bp][2][side[i]];
		arrs[3] += cCntLUT_4[(int)l._bp][3][side[i]];
	}
#endif
}

#ifdef BOWTIE2
/**
 * Counts the number of occurrences of all four nucleotides in the
 * given side from the given byte/bitpair (l->_by/l->_bp) (or the
 * beginning of the side if l == 0).  Count for 'a' goes in arrs[0],
 * 'c' in arrs[1], etc.
 *
 * Note: must account for $.
 *
 * Must fill in masks
 */
inline uint32_t Ebwt::countBt2SideRange2(
	const SideLocus& l,
	bool startAtLocus,
	uint32_t num,
	uint32_t* arrs,
	EList<bool> *masks,
	uint32_t maskOff) const
{
	assert(!masks[0].empty());
	assert_eq(masks[0].size(), masks[1].size());
	assert_eq(masks[0].size(), masks[2].size());
	assert_eq(masks[0].size(), masks[3].size());
	ASSERT_ONLY(uint32_t myarrs[4] = {0, 0, 0, 0});
	uint32_t nm = 0; // number of nucleotides tallied so far
	int iby = 0;      // initial byte offset
	int ibp = 0;      // initial base-pair offset
	if(startAtLocus) {
		iby = l._by;
		ibp = l._bp;
	} else {
		// Start at beginning
	}
	int by = iby, bp = ibp;
#ifdef SIXTY4_FORMAT
	throw 1; // Unsupported for now
#else
	assert_lt(bp, 4);
	assert_lt(by, (int)this->_eh._sideBwtSz);
	const uint8_t *side = l.side(this->ebwt());
	while(nm < num) {
		int c = (side[by] >> (bp * 2)) & 3;
		assert_lt(maskOff + nm, masks[c].size());
		masks[0][maskOff + nm] = masks[1][maskOff + nm] =
		masks[2][maskOff + nm] = masks[3][maskOff + nm] = false;
		assert_range(0, 3, c);
		// Note: we tally $ just like an A
		arrs[c]++; // tally it
		ASSERT_ONLY(myarrs[c]++);
		masks[c][maskOff + nm] = true; // not dead
		nm++;
		if(++bp == 4) {
			bp = 0;
			by++;
			assert_leq(by, (int)this->_eh._sideBwtSz);
			if(by == (int)this->_eh._sideBwtSz) {
				// Fell off the end of the side
				break;
			}
		}
	}
	WITHIN_FCHR_DOLLARA(arrs);
#endif
#ifndef NDEBUG
	if(_sanity) {
		// Make sure results match up with a call to mapLFEx.
		uint32_t tops[4] = {0, 0, 0, 0};
		uint32_t bots[4] = {0, 0, 0, 0};
		uint32_t top = l.toBWRow();
		uint32_t bot = top + nm;
		mapLFEx(top, bot, tops, bots);
		assert(myarrs[0] == (bots[0] - tops[0]) || myarrs[0] == (bots[0] - tops[0])+1);
		assert_eq(myarrs[1], bots[1] - tops[1]);
		assert_eq(myarrs[2], bots[2] - tops[2]);
		assert_eq(myarrs[3], bots[3] - tops[3]);
	}
#endif
	return nm;
}
#else
/**
 * Counts the number of occurrences of all four nucleotides in the
 * given side from the given byte/bitpair (l->_by/l->_bp) (or the
 * beginning of the side if l == 0).  Count for 'a' goes in arrs[0],
 * 'c' in arrs[1], etc.
 *
 * Note: must account for $.
 *
 * Must fill in masks
 */
inline uint32_t Ebwt::countFwSideRange2(
	const SideLocus& l,
	bool startAtLocus,
	uint32_t num,
	uint32_t* arrs,
	EList<bool> *masks,
	uint32_t maskOff) const
{
	assert(!masks[0].empty());
	assert_eq(masks[0].size(), masks[1].size());
	assert_eq(masks[0].size(), masks[2].size());
	assert_eq(masks[0].size(), masks[3].size());
	ASSERT_ONLY(uint32_t myarrs[4] = {0, 0, 0, 0});
	uint32_t nm = 0; // number of nucleotides tallied so far
	int iby = 0;      // initial byte offset
	int ibp = 0;      // initial base-pair offset
	if(startAtLocus) {
		iby = l._by;
		ibp = l._bp;
	} else {
		// Start at beginning
	}
	int by = iby, bp = ibp;
#ifdef SIXTY4_FORMAT
	throw 1; // Unsupported for now
#else
	assert_lt(bp, 4);
	assert_lt(by, (int)this->_eh._sideBwtSz);
	const uint8_t *side = l.side(this->ebwt());
	while(nm < num) {
		int c = (side[by] >> (bp * 2)) & 3;
		assert_lt(maskOff + nm, masks[c].size());
		masks[0][maskOff + nm] = masks[1][maskOff + nm] =
		masks[2][maskOff + nm] = masks[3][maskOff + nm] = false;
		assert_range(0, 3, c);
		// Note: we tally $ just like an A
		arrs[c]++; // tally it
		ASSERT_ONLY(myarrs[c]++);
		masks[c][maskOff + nm] = true; // not dead
		nm++;
		if(++bp == 4) {
			bp = 0;
			by++;
			assert_leq(by, (int)this->_eh._sideBwtSz);
			if(by == (int)this->_eh._sideBwtSz) {
				// Fell off the end of the side
				break;
			}
		}
	}
	WITHIN_FCHR_DOLLARA(arrs);
#endif
#ifndef NDEBUG
	if(_sanity) {
		// Make sure results match up with a call to mapLFEx.
		uint32_t tops[4] = {0, 0, 0, 0};
		uint32_t bots[4] = {0, 0, 0, 0};
		uint32_t top = l.toBWRow();
		uint32_t bot = top + nm;
		mapLFEx(top, bot, tops, bots);
		assert(myarrs[0] == (bots[0] - tops[0]) || myarrs[0] == (bots[0] - tops[0])+1);
		assert_eq(myarrs[1], bots[1] - tops[1]);
		assert_eq(myarrs[2], bots[2] - tops[2]);
		assert_eq(myarrs[3], bots[3] - tops[3]);
	}
#endif
	return nm;
}

/**
 * Counts the number of occurrences of all four nucleotides in the
 * given side from the given byte/bitpair (l->_by/l->_bp) (or the
 * beginning of the side if l == 0).  Count for 'a' goes in arrs[0],
 * 'c' in arrs[1], etc.
 *
 * Note: must account for $.
 */
inline uint32_t Ebwt::countBwSideRange2(
	const SideLocus& l,
	bool startAtLocus,
	uint32_t num,
	uint32_t* arrs,
	EList<bool> *masks,
	uint32_t maskOff) const
{
	assert(!masks[0].empty());
	assert_eq(masks[0].size(), masks[1].size());
	assert_eq(masks[0].size(), masks[2].size());
	assert_eq(masks[0].size(), masks[3].size());
	ASSERT_ONLY(uint32_t myarrs[4] = {0, 0, 0, 0});
	uint32_t nm = 0; // number of nucleotides tallied so far
	int iby = this->_eh._sideBwtSz-1; // initial byte offset
	int ibp = 3;      // initial base-pair offset
	if(startAtLocus) {
		iby = l._by;
		ibp = l._bp;
	} else {
		// Start at beginning
	}
	int by = iby, bp = ibp;
#ifdef SIXTY4_FORMAT
	throw 1; // Unsupported for now
#else
	assert_lt(bp, 4);
	assert_lt(by, (int)this->_eh._sideBwtSz);
	const uint8_t *side = l.side(this->ebwt());
	while(nm < num) {
		int c = (side[by] >> (bp * 2)) & 3;
		assert_lt(maskOff + nm, masks[c].size());
		masks[0][maskOff + nm] = masks[1][maskOff + nm] =
		masks[2][maskOff + nm] = masks[3][maskOff + nm] = false;
		assert_range(0, 3, c);
		// Note: we tally $ just like an A
		arrs[c]++; // tally it
		ASSERT_ONLY(myarrs[c]++);
		masks[c][maskOff + nm] = true;
		nm++;
		if(--bp == -1) {
			bp = 3;
			by--;
			assert_geq(by, -1);
			if(by == -1) {
				// Fell off the end of the side
				break;
			}
		}
	}
	WITHIN_FCHR_DOLLARA(arrs);
#endif
#ifndef NDEBUG
	if(_sanity) {
		// Make sure results match up with a call to mapLFEx.
		uint32_t tops[4] = {0, 0, 0, 0};
		uint32_t bots[4] = {0, 0, 0, 0};
		uint32_t top = l.toBWRow();
		uint32_t bot = top + nm;
		mapLFEx(top, bot, tops, bots);
		assert(myarrs[0] == (bots[0] - tops[0]) || myarrs[0] == (bots[0] - tops[0])+1);
		assert_eq(myarrs[1], bots[1] - tops[1]);
		assert_eq(myarrs[2], bots[2] - tops[2]);
		assert_eq(myarrs[3], bots[3] - tops[3]);
	}
#endif
	return nm;
}
#endif

/**
 * Returns true iff the index contains the given string (exactly).  The given
 * string must contain only unambiguous characters.  TODO: support ambiguous
 * characters in 'str'.
 */
bool Ebwt::contains(
	const BTDnaString& str,
	uint32_t *otop,
	uint32_t *obot) const
{
	assert(isInMemory());
	SideLocus tloc, bloc;
	if(str.empty()) {
		if(otop != NULL && obot != NULL) *otop = *obot = 0;
		return true;
	}
	int c = str[str.length()-1];
	assert_range(0, 4, c);
	uint32_t top = 0, bot = 0;
	if(c < 4) {
		top = fchr()[c];
		bot = fchr()[c+1];
	} else {
		bool set = false;
		for(int i = 0; i < 4; i++) {
			if(fchr()[c] < fchr()[c+1]) {
				if(set) {
					return false;
				} else {
					set = true;
					top = fchr()[c];
					bot = fchr()[c+1];
				}
			}
		}
	}
	assert_geq(bot, top);
	tloc.initFromRow(top, eh(), ebwt());
	bloc.initFromRow(bot, eh(), ebwt());
	ASSERT_ONLY(uint32_t lastDiff = bot - top);
	for(int i = (int)str.length()-2; i >= 0; i--) {
		c = str[i];
		assert_range(0, 4, c);
		if(c <= 3) {
			top = mapLF(tloc, c);
			bot = mapLF(bloc, c);
		} else {
			size_t sz = bot - top;
			int c1 = mapLF1(top, tloc);
			bot = mapLF(bloc, c1);
			assert_leq(bot - top, sz);
			if(bot - top < sz) {
				// Encountered an N and could not proceed through it because
				// there was more than one possible nucleotide we could replace
				// it with
				return false;
			}
		}
		assert_geq(bot, top);
		assert_leq(bot-top, lastDiff);
		ASSERT_ONLY(lastDiff = bot-top);
		if(i > 0) {
			tloc.initFromRow(top, eh(), ebwt());
			bloc.initFromRow(bot, eh(), ebwt());
		}
	}
	if(otop != NULL && obot != NULL) {
		*otop = top; *obot = bot;
	}
	return bot > top;
}

/**
 * Try to find the Bowtie index specified by the user.  First try the
 * exact path given by the user.  Then try the user-provided string
 * appended onto the path of the "indexes" subdirectory below this
 * executable, then try the provided string appended onto
 * "$BOWTIE2_INDEXES/".
 */
string adjustEbwtBase(const string& cmdline,
					  const string& ebwtFileBase,
					  bool verbose = false)
{
	string str = ebwtFileBase;
	ifstream in;
	if(verbose) cout << "Trying " << str << endl;
	in.open((str + ".1.bt2").c_str(), ios_base::in | ios::binary);
	if(!in.is_open()) {
		if(verbose) cout << "  didn't work" << endl;
		in.close();
		if(getenv("BOWTIE2_INDEXES") != NULL) {
			str = string(getenv("BOWTIE2_INDEXES")) + "/" + ebwtFileBase;
			if(verbose) cout << "Trying " << str << endl;
			in.open((str + ".1.bt2").c_str(), ios_base::in | ios::binary);
			if(!in.is_open()) {
				if(verbose) cout << "  didn't work" << endl;
				in.close();
			} else {
				if(verbose) cout << "  worked" << endl;
			}
		}
	}
	if(!in.is_open()) {
		cerr << "Could not locate a Bowtie index corresponding to basename \"" << ebwtFileBase << "\"" << endl;
		throw 1;
	}
	return str;
}
