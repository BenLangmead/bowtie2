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

#include "dp_framer.h"

using namespace std;

/**
 * Set up variables that describe the shape of a dynamic programming matrix to
 * be filled in.  The matrix is built around the diagonal containing the seed
 * hit: the "seed diagonal".  The N diagonals to the right of the seed diagonal
 * are the "RHS gap" diagonals, where N is the maximum number of read or
 * reference gaps permitted (whichever is larger).  The N diagonals to the left
 * of the seed diagonal are the "LHS gap" diagonals.
 *
 * The way the rectangle is currently formulated, there are another N diagonals
 * to the left of the "LHS gap" diagonals called the "LHS extra diagonals".  It
 * might also be possible to split the "extra diagonals" into two subsets and
 * place them both to the left of the LHS gap diagonals and to the right of the
 * RHS gap diagonals.
 *
 * The purpose of arranging and these groupings of diagonals is that a subset
 * of them, the "core diagonals", can now be considered "covered."  By
 * "covered" I mean that any alignment that overlaps a cell in any of the core
 * diagonals cannot possibly overlap another, higher-scoring alignment that
 * falls partially outside the rectangle.
 *
 * Say the read is 5 characters long, the maximum number of read or ref gaps is
 * 2, and the seed hit puts the main diagonal at offset 10 in the reference.
 * The larger rectangle explored looks like this:
 *
 *  off=10, maxgap=2
 *
 * Ref      1
 * off: 67890123456   0: seed diagonal
 *      **OO0oo++----   o: "RHS gap" diagonals
 *      -**OO0oo++---   O: "LHS gap" diagonals
 *      --**OO0oo++--   *: "LHS extra" diagonals
 *      ---**OO0oo++-   +: "RHS extra" diagonals
 *      ----**OO0oo++   -: cells that can't possibly be involved in a valid    
 *                         alignment that overlaps one of the core diagonals
 *
 * The "core diagonals" are marked with 0's, O's or o's.
 *
 * A caveat is that, for performance reasons, we place an upper limit on N -
 * the maximum number of read or reference gaps.  It is constrained to be no
 * greater than 'maxgap'.  This means that in some situations, we may report an
 * alignment that spuriously trumps a better alignment that falls partially
 * outside the rectangle.  Also, we may fail to find a valid alignment with
 * more than 'maxgap' gaps.
 *
 * Another issue is trimming: if the seed hit is sufficiently close to one or
 * both ends of the reference sequence, and either (a) overhang is not
 * permitted, or (b) the number of Ns permitted is less than the number of
 * columns that overhang the reference, then we want to exclude the trimmed
 * columns from the rectangle.
 *
 * We need to return enough information so that downstream routines can fully
 * understand the shape of the rectangle, which diagonals are which (esp. which
 * are the "core" diagonals, since we needn't examine any more seed hits from
 * those columns in the future), and how the rectangle is trimmed.  The
 * information returned should be compatible with the sort of information
 * returned by the routines that set up rectangles for mate finding.
 */
bool DynProgFramer::frameSeedExtensionRect(
	int64_t  off,      // ref offset implied by seed hit assuming no gaps
	size_t   rdlen,    // length of read sequence used in DP table (so len
	                   // of +1 nucleotide sequence for colorspace reads)
	int64_t  reflen,   // length of reference sequence aligned to
	size_t   maxrdgap, // max # of read gaps permitted in opp mate alignment
	size_t   maxrfgap, // max # of ref gaps permitted in opp mate alignment
	int64_t  maxns,    // # Ns permitted
	size_t   maxhalf,  // max width in either direction
	DPRect&  rect)     // out: DP rectangle
{
	assert_gt(rdlen, 0);
	assert_gt(reflen, 0);
	// Set N, the maximum number of reference or read gaps permitted, whichever
	// is larger.  Also, enforce ceiling: can't be larger than 'maxhalf'.
	size_t maxgap = max(maxrdgap, maxrfgap);
	maxgap = min(maxgap, maxhalf);
	// Leave room for "LHS gap" and "LHS extra" diagonals
	int64_t refl = off - 2 * maxgap;               // inclusive
	// Leave room for "RHS gap" and "RHS extra" diagonals
	int64_t refr = off + (rdlen - 1) + 2 * maxgap; // inclusive
	size_t triml = 0, trimr = 0;
	// Check if we have to trim to fit the extents of the reference
	if(trimToRef_) {
		maxns = 0; // no leeway
	} else if(maxns == (int64_t)rdlen) {
		maxns--;
	}
	// Trim from RHS of rectangle
	if(refr >= reflen + maxns) {
		trimr = (size_t)(refr - (reflen + maxns - 1));
	}
	// Trim from LHS of rectangle
	if(refl < -maxns) {
		triml = (size_t)(-refl) - (size_t)maxns;
	}
	rect.refl_pretrim = refl;
	rect.refr_pretrim = refr;
	rect.refl  = refl + triml;
	rect.refr  = refr - trimr;
	rect.triml = triml;
	rect.trimr = trimr;
	rect.maxgap = maxgap;
	// Remember which diagonals are "core" as offsets from the LHS of the
	// untrimmed rectangle
	rect.corel = maxgap;
	rect.corer = rect.corel + 2 * maxgap; // inclusive
	assert(rect.repOk());
	return !rect.entirelyTrimmed();
}

/**
 * Set up variables that describe the shape of a dynamic programming matrix to
 * be filled in.  The matrix is built around the diagonals that terminate in
 * the range of columns where the RHS of the opposite mate must fall in order
 * to satisfy the fragment-length constraint.  These are the "mate" diagonals
 * and they also happen to be the "core" diagonals in this case.
 *
 * The N diagonals to the right of the mate diagonals are the "RHS gap"
 * diagonals, where N is the maximum number of read or reference gaps permitted
 * (whichever is larger).  The N diagonals to the left of the mate diagonals
 * are the "LHS gap" diagonals.
 *
 * The purpose of arranging and these groupings of diagonals is that a subset
 * of them, the "core diagonals", can now be considered "covered."  By
 * "covered" I mean that any alignment that overlaps a cell in any of the core
 * diagonals cannot possibly overlap another, higher-scoring alignment that
 * falls partially outside the rectangle.
 *
 *   |Anchor| 
 *   o---------OO0000000000000oo------  0: mate diagonal (also core diags!)
 *   -o---------OO0000000000000oo-----  o: "RHS gap" diagonals
 *   --o---------OO0000000000000oo----  O: "LHS gap" diagonals
 *   ---oo--------OO0000000000000oo---  *: "LHS extra" diagonals
 *   -----o--------OO0000000000000oo--  -: cells that can't possibly be
 *   ------o--------OO0000000000000oo-     involved in a valid alignment that
 *   -------o--------OO0000000000000oo     overlaps one of the core diagonals
 *                     XXXXXXXXXXXXX
 *                     | RHS Range |
 *                     ^           ^
 *                     rl          rr
 *
 * The "core diagonals" are marked with 0s.
 *
 * A caveat is that, for performance reasons, we place an upper limit on N -
 * the maximum number of read or reference gaps.  It is constrained to be no
 * greater than 'maxgap'.  This means that in some situations, we may report an
 * alignment that spuriously trumps a better alignment that falls partially
 * outside the rectangle.  Also, we may fail to find a valid alignment with
 * more than 'maxgap' gaps.
 *
 * Another issue is trimming: if the seed hit is sufficiently close to one or
 * both ends of the reference sequence, and either (a) overhang is not
 * permitted, or (b) the number of Ns permitted is less than the number of
 * columns that overhang the reference, then we want to exclude the trimmed
 * columns from the rectangle.
 */
bool DynProgFramer::frameFindMateAnchorLeftRect(
	int64_t ll,       // leftmost Watson off for LHS of opp alignment
	int64_t lr,       // rightmost Watson off for LHS of opp alignment
	int64_t rl,       // leftmost Watson off for RHS of opp alignment
	int64_t rr,       // rightmost Watson off for RHS of opp alignment
	size_t  rdlen,    // length of opposite mate
	int64_t reflen,   // length of reference sequence aligned to
	size_t  maxrdgap, // max # of read gaps permitted in opp mate alignment
	size_t  maxrfgap, // max # of ref gaps permitted in opp mate alignment
	int64_t maxns,    // max # ns permitted in the alignment
	size_t  maxhalf,  // max width in either direction
	DPRect& rect)     // out: DP rectangle
	const
{
	assert_geq(lr, ll);  // LHS rightmost must be >= LHS leftmost
	assert_geq(rr, rl);  // RHS rightmost must be >= RHS leftmost
	assert_geq(rr, lr);  // RHS rightmost must be >= LHS rightmost
	assert_geq(rl, ll);  // RHS leftmost must be >= LHS leftmost
	assert_gt(rdlen, 0);
	assert_gt(reflen, 0);
	size_t triml = 0, trimr = 0;
	size_t maxgap = max(maxrdgap, maxrfgap);
	maxgap = max(maxgap, maxhalf);
	// Amount of padding we have to add to account for the fact that alignments
	// ending between en_left/en_right might start in various columns in the
	// first row
	int64_t pad_left = maxgap;
	int64_t pad_right = maxgap;
	int64_t en_left  = rl;
	int64_t en_right = rr;
	int64_t st_left  = en_left - (rdlen-1);
	ASSERT_ONLY(int64_t st_right = en_right - (rdlen-1));
	int64_t en_right_pad = en_right + pad_right;
	ASSERT_ONLY(int64_t en_left_pad  = en_left  - pad_left);
	ASSERT_ONLY(int64_t st_right_pad = st_right + pad_right);
	int64_t st_left_pad  = st_left  - pad_left;
	assert_leq(st_left, en_left);
	assert_geq(en_right, st_right);
	assert_leq(st_left_pad, en_left_pad);
	assert_geq(en_right_pad, st_right_pad);
	int64_t refl = st_left_pad;
	int64_t refr = en_right_pad;
	if(trimToRef_) {
		maxns = 0;
	} else if(maxns == (int64_t)rdlen) {
		maxns--;
	}
	// Trim from the RHS of the rectangle?
	if(refr >= reflen + maxns) {
		trimr = (size_t)(refr - (reflen + maxns - 1));
	}
	// Trim from the LHS of the rectangle?
	if(refl < -maxns) {
		triml = (size_t)(-refl) - (size_t)maxns;
	}
	size_t width = (size_t)(refr - refl + 1);
	rect.refl_pretrim = refl;
	rect.refr_pretrim = refr;
	rect.refl  = refl + triml;
	rect.refr  = refr - trimr;
	rect.triml = triml;
	rect.trimr = trimr;
	rect.maxgap = maxgap;
	rect.corel = maxgap;
	rect.corer = width - maxgap - 1; // inclusive
	assert(rect.repOk());
	return !rect.entirelyTrimmed();
}

/**
 * Set up variables that describe the shape of a dynamic programming matrix to
 * be filled in.  The matrix is built around the diagonals that begin in the
 * range of columns where the LHS of the opposite mate must fall in order to
 * satisfy the fragment-length constraint.  These are the "mate" diagonals and
 * they also happen to be the "core" diagonals in this case.
 *
 * The N diagonals to the right of the mate diagonals are the "RHS gap"
 * diagonals, where N is the maximum number of read or reference gaps permitted
 * (whichever is larger).  The N diagonals to the left of the mate diagonals
 * are the "LHS gap" diagonals.
 *
 * The purpose of arranging and these groupings of diagonals is that a subset
 * of them, the "core diagonals", can now be considered "covered."  By
 * "covered" I mean that any alignment that overlaps a cell in any of the core
 * diagonals cannot possibly overlap another, higher-scoring alignment that
 * falls partially outside the rectangle.
 *
 *    ll          lr
 *    v           v
 *    | LHS Range |
 *    XXXXXXXXXXXXX          |Anchor|
 *  OO0000000000000oo--------o--------  0: mate diagonal (also core diags!)
 *  -OO0000000000000oo--------o-------  o: "RHS gap" diagonals
 *  --OO0000000000000oo--------o------  O: "LHS gap" diagonals
 *  ---OO0000000000000oo--------oo----  *: "LHS extra" diagonals
 *  ----OO0000000000000oo---------o---  -: cells that can't possibly be
 *  -----OO0000000000000oo---------o--     involved in a valid alignment that
 *  ------OO0000000000000oo---------o-     overlaps one of the core diagonals
 *
 * The "core diagonals" are marked with 0s.
 *
 * A caveat is that, for performance reasons, we place an upper limit on N -
 * the maximum number of read or reference gaps.  It is constrained to be no
 * greater than 'maxgap'.  This means that in some situations, we may report an
 * alignment that spuriously trumps a better alignment that falls partially
 * outside the rectangle.  Also, we may fail to find a valid alignment with
 * more than 'maxgap' gaps.
 *
 * Another issue is trimming: if the seed hit is sufficiently close to one or
 * both ends of the reference sequence, and either (a) overhang is not
 * permitted, or (b) the number of Ns permitted is less than the number of
 * columns that overhang the reference, then we want to exclude the trimmed
 * columns from the rectangle.
 */
bool DynProgFramer::frameFindMateAnchorRightRect(
	int64_t ll,       // leftmost Watson off for LHS of opp alignment
	int64_t lr,       // rightmost Watson off for LHS of opp alignment
	int64_t rl,       // leftmost Watson off for RHS of opp alignment
	int64_t rr,       // rightmost Watson off for RHS of opp alignment
	size_t rdlen,     // length of opposite mate
	int64_t reflen,   // length of reference sequence aligned to
	size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
	size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
	int64_t maxns,    // max # ns permitted in the alignment
	size_t maxhalf,   // max width in either direction
	DPRect& rect)     // out: DP rectangle
	const
{
	assert_geq(lr, ll);
	assert_geq(rr, rl);
	assert_geq(rr, lr);
	assert_geq(rl, ll);
	assert_gt(rdlen, 0);
	assert_gt(reflen, 0);
	size_t triml = 0, trimr = 0;
	size_t maxgap = max(maxrdgap, maxrfgap);
	maxgap = max(maxgap, maxhalf);
	int64_t pad_left = maxgap;
	int64_t pad_right = maxgap;
	int64_t st_left = ll;
	int64_t st_right = lr;
	ASSERT_ONLY(int64_t en_left = st_left + (rdlen-1));
	int64_t en_right = st_right + (rdlen-1);
	int64_t en_right_pad = en_right + pad_right;
	ASSERT_ONLY(int64_t en_left_pad  = en_left  - pad_left);
	ASSERT_ONLY(int64_t st_right_pad = st_right + pad_right);
	int64_t st_left_pad  = st_left  - pad_left;
	assert_leq(st_left, en_left);
	assert_geq(en_right, st_right);
	assert_leq(st_left_pad, en_left_pad);
	assert_geq(en_right_pad, st_right_pad);
	// We have enough info to deduce where the boundaries of our rectangle
	// should be.  Finalize the boundaries, ignoring reference trimming for now
	int64_t refl = st_left_pad;
	int64_t refr = en_right_pad;
	if(trimToRef_) {
		maxns = 0;
	} else if(maxns == (int64_t)rdlen) {
		maxns--;
	}
	// Trim from the RHS of the rectangle?
	if(refr >= reflen + maxns) {
		trimr = (size_t)(refr - (reflen + maxns - 1));
	}
	// Trim from the LHS of the rectangle?
	if(refl < -maxns) {
		triml = (size_t)(-refl) - (size_t)maxns;
	}
	size_t width = (size_t)(refr - refl + 1);
	rect.refl_pretrim = refl;
	rect.refr_pretrim = refr;
	rect.refl  = refl + triml;
	rect.refr  = refr - trimr;
	rect.triml = triml;
	rect.trimr = trimr;
	rect.maxgap = maxgap;
	rect.corel = maxgap;
	rect.corer = width - maxgap - 1; // inclusive
	assert(rect.repOk());
	return !rect.entirelyTrimmed();
}

#ifdef MAIN_DP_FRAMER

#include <iostream>

static void testCaseFindMateAnchorLeft(
	const char *testName,
	bool trimToRef,
	int64_t ll,
	int64_t lr,
	int64_t rl,
	int64_t rr,
	size_t rdlen,
	size_t reflen,
	size_t maxrdgap,
	size_t maxrfgap,
	size_t ex_width,
	size_t ex_solwidth,
	size_t ex_trimup,
	size_t ex_trimdn,
	int64_t ex_refl,
	int64_t ex_refr,
	const char *ex_st,    // string of '0'/'1' chars
	const char *ex_en)    // string of '0'/'1' chars
{
	cerr << testName << "...";
	DynProgFramer fr(trimToRef);
	size_t width, solwidth;
	int64_t refl, refr;
	EList<bool> st, en;
	size_t trimup, trimdn;
	size_t maxhalf = 500;
	size_t maxgaps = 0;
	fr.frameFindMateAnchorLeft(
		ll,       // leftmost Watson off for LHS of opp alignment
		lr,       // rightmost Watson off for LHS of opp alignment
		rl,       // leftmost Watson off for RHS of opp alignment
		rr,       // rightmost Watson off for RHS of opp alignment
		rdlen,    // length of opposite mate
		reflen,   // length of reference sequence aligned to
		maxrdgap, // max # of read gaps permitted in opp mate alignment
		maxrfgap, // max # of ref gaps permitted in opp mate alignment
		maxns,    // max # Ns permitted
		maxhalf,  // max width in either direction
		width,    // out: calculated width stored here
		maxgaps,  // out: max # gaps
		trimup,   // out: number of bases trimmed from upstream end
		trimdn,   // out: number of bases trimmed from downstream end
		refl,     // out: ref pos of upper LHS of parallelogram
		refr,     // out: ref pos of lower RHS of parallelogram
		st,       // out: legal starting columns stored here
		en);      // out: legal ending columns stored here
	assert_eq(ex_width, width);
	assert_eq(ex_solwidth, solwidth);
	assert_eq(ex_trimup, trimup);
	assert_eq(ex_trimdn, trimdn);
	assert_eq(ex_refl, refl);
	assert_eq(ex_refr, refr);
	for(size_t i = 0; i < width; i++) {
		assert_eq((ex_st[i] == '1'), st[i]);
		assert_eq((ex_en[i] == '1'), en[i]);
	}
	cerr << "PASSED" << endl;
}

static void testCaseFindMateAnchorRight(
	const char *testName,
	bool trimToRef,
	int64_t ll,
	int64_t lr,
	int64_t rl,
	int64_t rr,
	size_t rdlen,
	size_t reflen,
	size_t maxrdgap,
	size_t maxrfgap,
	size_t ex_width,
	size_t ex_solwidth,
	size_t ex_trimup,
	size_t ex_trimdn,
	int64_t ex_refl,
	int64_t ex_refr,
	const char *ex_st,    // string of '0'/'1' chars
	const char *ex_en)    // string of '0'/'1' chars
{
	cerr << testName << "...";
	DynProgFramer fr(trimToRef);
	size_t width, solwidth;
	size_t maxgaps;
	int64_t refl, refr;
	EList<bool> st, en;
	size_t trimup, trimdn;
	size_t maxhalf = 500;
	fr.frameFindMateAnchorRight(
		ll,       // leftmost Watson off for LHS of opp alignment
		lr,       // rightmost Watson off for LHS of opp alignment
		rl,       // leftmost Watson off for RHS of opp alignment
		rr,       // rightmost Watson off for RHS of opp alignment
		rdlen,    // length of opposite mate
		reflen,   // length of reference sequence aligned to
		maxrdgap, // max # of read gaps permitted in opp mate alignment
		maxrfgap, // max # of ref gaps permitted in opp mate alignment
		maxns,    // max # Ns permitted
		maxhalf,  // max width in either direction
		width,    // out: calculated width stored here
		maxgaps,  // out: calcualted max # gaps
		trimup,   // out: number of bases trimmed from upstream end
		trimdn,   // out: number of bases trimmed from downstream end
		refl,     // out: ref pos of upper LHS of parallelogram
		refr,     // out: ref pos of lower RHS of parallelogram
		st,       // out: legal starting columns stored here
		en);      // out: legal ending columns stored here
	assert_eq(ex_width, width);
	assert_eq(ex_trimup, trimup);
	assert_eq(ex_trimdn, trimdn);
	assert_eq(ex_refl, refl);
	assert_eq(ex_refr, refr);
	for(size_t i = 0; i < width; i++) {
		assert_eq((ex_st[i] == '1'), st[i]);
		assert_eq((ex_en[i] == '1'), en[i]);
	}
	cerr << "PASSED" << endl;
}

int main(void) {
	
	///////////////////////////
	//
	// ANCHOR ON THE LEFT
	//
	///////////////////////////

	//    -------------
	//       o     o
	//        o     o
	//         o     o
	//          o     o
	//        <<<------->>>
	// 012345678901234567890
	// 0         1         2
	testCaseFindMateAnchorLeft(
		"FindMateAnchorLeft1",
		false,            // trim to reference
		3,                // left offset of upper parallelogram extent
		15,               // right offset of upper parallelogram extent
		10,               // left offset of lower parallelogram extent
		16,               // right offset of lower parallelogram extent
		5,                // length of opposite mate
		30,               // length of reference sequence aligned to
		3,                // max # of read gaps permitted in opp mate alignment
		3,                // max # of ref gaps permitted in opp mate alignment
		13,               // expected width
		0,                // expected # bases trimmed from upstream end
		0,                // expected # bases trimmed from downstream end
		3,                // ref offset of upstream column
		19,               // ref offset of downstream column
		"1111111111111",  // expected starting bools
		"0001111111000"); // expected ending bools

	//        *******
	//     <<===-----
	//       o    o
	//        o    o
	//         o    o
	//          o    o
	//         <<=----->>
	//            *******
	// 012345678901234567890
	// 0         1         2
	testCaseFindMateAnchorLeft(
		"FindMateAnchorLeft2",
		false,            // trim to reference
		9,                // left offset of left upper parallelogram extent
		14,               // right offset of left upper parallelogram extent
		10,               // left offset of left lower parallelogram extent
		15,               // right offset of left lower parallelogram extent
		5,                // length of opposite mate
		30,               // length of reference sequence aligned to
		2,                // max # of read gaps permitted in opp mate alignment
		2,                // max # of ref gaps permitted in opp mate alignment
		7,                // expected width
		3,                // expected # bases trimmed from upstream end
		0,                // expected # bases trimmed from downstream end
		7,                // ref offset of upstream column
		17,               // ref offset of downstream column
		"0011111",        // expected starting bools
		"1111100");       // expected ending bools

	//        *******
	//     <<===--->>
	//       o    o
	//        o    o
	//         o    o
	//          o    o
	//           o    o
	//         <<=----->>
	//            *******
	// 01234567890123456xxxx
	// 0         1         2
	testCaseFindMateAnchorLeft(
		"FindMateAnchorLeft3",
		true,             // trim to reference
		9,                // left offset of left upper parallelogram extent
		14,               // right offset of left upper parallelogram extent
		10,               // left offset of left lower parallelogram extent
		15,               // right offset of left lower parallelogram extent
		5,                // length of opposite mate
		17,               // length of reference sequence aligned to
		2,                // max # of read gaps permitted in opp mate alignment
		2,                // max # of ref gaps permitted in opp mate alignment
		7,                // expected width
		3,                // expected # bases trimmed from upstream end
		0,                // expected # bases trimmed from downstream end
		7,                // ref offset of upstream column
		17,               // ref offset of downstream column
		"0011111",        // expected starting bools
		"1111100");       // expected ending bools

	//        ******
	//     <<===-----
	//       o    o
	//        o    o
	//         o    o
	//          o    o
	//         <<=----=>>
	//            ******
	// 012345678901234xxxxxx
	// 0         1         2
	testCaseFindMateAnchorLeft(
		"FindMateAnchorLeft4",
		true,             // trim to reference
		9,                // left offset of left upper parallelogram extent
		14,               // right offset of left upper parallelogram extent
		10,               // left offset of left lower parallelogram extent
		15,               // right offset of left lower parallelogram extent
		5,                // length of opposite mate
		15,               // length of reference sequence aligned to
		2,                // max # of read gaps permitted in opp mate alignment
		2,                // max # of ref gaps permitted in opp mate alignment
		6,                // expected width
		3,                // expected # bases trimmed from upstream end
		1,                // expected # bases trimmed from downstream end
		7,                // ref offset of upstream column
		16,               // ref offset of downstream column
		"001111",         // expected starting bools
		"111100");        // expected ending bools

	// -1         0         2
	//  xxxxxxxxxx012345678xx
	//
	//           *******
	//        <<===-----
	//          o    o
	//           o    o
	//            o    o
	//             o    o
	//              o    o
	//            <<=----->>
	//               *******
	//                
	//  xxxxxxxxxx012345678xx
	// -1         0         2
	testCaseFindMateAnchorLeft(
		"FindMateAnchorLeft5",
		true,             // trim to reference
		1,                // left offset of left upper parallelogram extent
		7,                // right offset of left upper parallelogram extent
		2,                // left offset of left lower parallelogram extent
		7,                // right offset of left lower parallelogram extent
		5,                // length of opposite mate
		9,                // length of reference sequence aligned to
		2,                // max # of read gaps permitted in opp mate alignment
		2,                // max # of ref gaps permitted in opp mate alignment
		7,                // expected width
		3,                // expected # bases trimmed from upstream end
		0,                // expected # bases trimmed from downstream end
		-1,               // ref offset of upstream column
		9,                // ref offset of downstream column
		"0011111",        // expected starting bools
		"1111100");       // expected ending bools

	//   <<<<==-===>>
	//       o    o
	//        o    o
	//         o    o
	//          o    o
	//       <<<<------>>
	//           ******
	// 012345678901234567890
	// 0         1         2
	testCaseFindMateAnchorLeft(
		"FindMateAnchorLeft6",
		false,            // trim to reference
		8,                // left offset of left upper parallelogram extent
		8,                // right offset of left upper parallelogram extent
		10,               // left offset of left lower parallelogram extent
		15,               // right offset of left lower parallelogram extent
		5,                // length of opposite mate
		30,               // length of reference sequence aligned to
		4,                // max # of read gaps permitted in opp mate alignment
		2,                // max # of ref gaps permitted in opp mate alignment
		6,                // expected width
		4,                // expected # bases trimmed from upstream end
		2,                // expected # bases trimmed from downstream end
		6,                // ref offset of upstream column
		15,               // ref offset of downstream column
		"001000",         // expected starting bools
		"111111");        // expected ending bools

	///////////////////////////
	//
	// ANCHOR ON THE RIGHT
	//
	///////////////////////////

	//        <<<------->>>
	//           o     o
	//            o     o
	//             o     o
	//              o     o
	//            <<<------->>>
	// 012345678901234567890123456789
	// 0         1         2
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight1",
		false,            // trim to reference
		10,               // left offset of left upper parallelogram extent
		16,               // right offset of left upper parallelogram extent
		11,               // left offset of left lower parallelogram extent
		23,               // right offset of left lower parallelogram extent
		5,                // length of opposite mate
		30,               // length of reference sequence aligned to
		3,                // max # of read gaps permitted in opp mate alignment
		3,                // max # of ref gaps permitted in opp mate alignment
		13,               // expected width
		0,                // expected # bases trimmed from upstream end
		0,                // expected # bases trimmed from downstream end
		7,                // ref offset of upstream column
		23,               // ref offset of downstream column
		"0001111111000",  // expected starting bools
		"1111111111111"); // expected ending bools

	// 0         1         2
	// 012345678901234567890
	//        *******
	//     <<------>>
	//        o    o
	//         o    o
	//          o    o
	//           o    o
	//         <<===--->>
	//            *******
	// 012345678901234567890
	// 0         1         2
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight2",
		false,            // trim to reference
		6,                // left offset of left upper parallelogram extent
		11,               // right offset of left upper parallelogram extent
		13,               // left offset of left lower parallelogram extent
		18,               // right offset of left lower parallelogram extent
		5,                // length of opposite mate
		30,               // length of reference sequence aligned to
		2,                // max # of read gaps permitted in opp mate alignment
		2,                // max # of ref gaps permitted in opp mate alignment
		7,                // expected width
		3,                // expected # bases trimmed from upstream end
		0,                // expected # bases trimmed from downstream end
		7,                // ref offset of upstream column
		17,               // ref offset of downstream column
		"1111100",        // expected starting bools
		"0011111");       // expected ending bools

	// Reference trimming takes off the left_pad of the left mate
	//
	//             *******
	//          <<------>>
	//            o    o
	//             o    o
	//              o    o
	//               o    o
	//                o    o
	//              <<===--->>
	//                 *******
	//  0123456789012345678901234567890
	// -1         0         1         2
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight3",
		true,             // trim to reference
		0,                // left offset of left upper parallelogram extent
		5,                // right offset of left upper parallelogram extent
		7,                // left offset of left lower parallelogram extent
		11,               // right offset of left lower parallelogram extent
		5,                // length of opposite mate
		30,               // length of reference sequence aligned to
		2,                // max # of read gaps permitted in opp mate alignment
		2,                // max # of ref gaps permitted in opp mate alignment
		7,                // expected width
		3,                // expected # bases trimmed from upstream end
		0,                // expected # bases trimmed from downstream end
		1,                // ref offset of upstream column
		11,               // ref offset of downstream column
		"1111100",        // expected starting bools
		"0011111");       // expected ending bools

	// Reference trimming takes off the leftmost 5 positions of the left mate,
	// and takes 1 from the right mate
	//
	//            *****
	//       <<------>>
	//         o    o
	//          o    o
	//           o    o
	//            o    o
	//             o    o
	//           <<===--->>
	//                *****
	//  0987654321012345678901234567890
	// -1         0         1         2
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight4",
		true,             // trim to reference
		-3,               // left offset of left upper parallelogram extent
		2,                // right offset of left upper parallelogram extent
		4,                // left offset of left lower parallelogram extent
		10,               // right offset of left lower parallelogram extent
		5,                // length of opposite mate
		30,               // length of reference sequence aligned to
		2,                // max # of read gaps permitted in opp mate alignment
		2,                // max # of ref gaps permitted in opp mate alignment
		5,                // expected width
		5,                // expected # bases trimmed from upstream end
		0,                // expected # bases trimmed from downstream end
		0,                // ref offset of upstream column
		8,                // ref offset of downstream column
		"11100",          // expected starting bools
		"11111");         // expected ending bools

	// Reference trimming takes off the leftmost 5 positions of the left mate,
	// and takes 1 from the left of the right mate.  Also, it takes 2 from the
	// right of the right mate.
	//
	//            ***
	//       <<------>>
	//         o    o
	//          o    o
	//           o    o
	//            o    o
	//             o    o
	//           <<===--->>
	//                ***
	//  0987654321012345678901234567890
	// -1         0         1         2
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight5",
		true,             // trim to reference
		-3,               // left offset of left upper parallelogram extent
		2,                // right offset of left upper parallelogram extent
		4,                // left offset of left lower parallelogram extent
		10,               // right offset of left lower parallelogram extent
		5,                // length of opposite mate
		7,                // length of reference sequence aligned to
		2,                // max # of read gaps permitted in opp mate alignment
		2,                // max # of ref gaps permitted in opp mate alignment
		3,                // expected width
		5,                // expected # bases trimmed from upstream end
		2,                // expected # bases trimmed from downstream end
		0,                // ref offset of upstream column
		6,                // ref offset of downstream column
		"111",            // expected starting bools
		"111");           // expected ending bools

	//       ******
	//     <<------>>>>
	//        o    o
	//         o    o
	//          o    o
	//           o    o
	//         <<====-=>>>>
	//           ******
	// 012345678901234567890
	// 0         1         2
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight6",
		false,            // trim to reference
		6,                // left offset of left upper parallelogram extent
		11,               // right offset of left upper parallelogram extent
		14,               // left offset of left lower parallelogram extent
		14,               // right offset of left lower parallelogram extent
		5,                // length of opposite mate
		30,               // length of reference sequence aligned to
		4,                // max # of read gaps permitted in opp mate alignment
		2,                // max # of ref gaps permitted in opp mate alignment
		6,                // expected width
		2,                // expected # bases trimmed from upstream end
		4,                // expected # bases trimmed from downstream end
		6,                // ref offset of upstream column
		15,               // ref offset of downstream column
		"111111",         // expected starting bools
		"000010");        // expected ending bools

	//         ****
	//   <<<<==---->>
	//       o    o
	//        o    o
	//         o    o
	//          o    o
	//           o    o
	//       <<<<====-=>>
	//             ****
	// 012345678901234567890
	// 0         1         2
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight7",
		false,            // trim to reference
		6,                // left offset of left upper parallelogram extent
		11,               // right offset of left upper parallelogram extent
		14,               // left offset of left lower parallelogram extent
		14,               // right offset of left lower parallelogram extent
		5,                // length of opposite mate
		30,               // length of reference sequence aligned to
		2,                // max # of read gaps permitted in opp mate alignment
		4,                // max # of ref gaps permitted in opp mate alignment
		4,                // expected width
		6,                // expected # bases trimmed from upstream end
		2,                // expected # bases trimmed from downstream end
		8,                // ref offset of upstream column
		15,               // ref offset of downstream column
		"1111",           // expected starting bools
		"0010");          // expected ending bools
	
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight8",
		true,             // trim to reference
		-37,              // left offset of left upper parallelogram extent
		13,               // right offset of left upper parallelogram extent
		-37,              // left offset of left lower parallelogram extent
		52,               // right offset of left lower parallelogram extent
		10,               // length of opposite mate
		53,               // length of reference sequence aligned to
		0,                // max # of read gaps permitted in opp mate alignment
		0,                // max # of ref gaps permitted in opp mate alignment
		14,               // expected width
		37,               // expected # bases trimmed from upstream end
		0,                // expected # bases trimmed from downstream end
		0,                // ref offset of upstream column
		22,               // ref offset of downstream column
		"11111111111111", // expected starting bools
		"11111111111111");// expected ending bools
}

#endif /*def MAIN_DP_FRAMER*/
