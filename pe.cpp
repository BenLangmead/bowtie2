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

#include "assert_helpers.h"
#include "pe.h"

using namespace std;

/**
 * Return a PE_TYPE flag indicating, given a PE_POLICY and coordinates
 * for a paired-end alignment, what type of alignment it is, i.e.,
 * whether it's:
 *
 * 1. Straightforwardly concordant
 * 2. Mates dovetail (one extends beyond the end of the other)
 * 3. One mate contains the other but they don't dovetail
 * 4. One mate overlaps the other but neither contains the other and
 *    they don't dovetail
 * 5. Discordant
 */
int PairedEndPolicy::peClassifyPair(
	int64_t  off1,   // offset of mate 1
	size_t   len1,   // length of mate 1
	bool     fw1,    // whether mate 1 aligned to Watson
	int64_t  off2,   // offset of mate 2
	size_t   len2,   // length of mate 2
	bool     fw2)    // whether mate 2 aligned to Watson
	const
{
	assert_gt(len1, 0);
	assert_gt(len2, 0);
	// Expand the maximum fragment length if necessary to accomodate
	// the longer mate
	size_t maxfrag = maxfrag_;
	if(len1 > maxfrag && expandToFit_) maxfrag = len1;
	if(len2 > maxfrag && expandToFit_) maxfrag = len2;
	size_t minfrag = minfrag_;
	if(minfrag < 1) {
		minfrag = 1;
	}
	bool oneLeft = false;
	if(pol_ == PE_POLICY_FF) {
		if(fw1 != fw2) {
			// Bad combination of orientations
			return PE_ALS_DISCORD;
		}
		oneLeft = fw1;
	} else if(pol_ == PE_POLICY_RR) {
		if(fw1 != fw2) {
			// Bad combination of orientations
			return PE_ALS_DISCORD;
		}
		oneLeft = !fw1;
	} else if(pol_ == PE_POLICY_FR) {
		if(fw1 == fw2) {
			// Bad combination of orientations
			return PE_ALS_DISCORD;
		}
		oneLeft = fw1;
	} else if(pol_ == PE_POLICY_RF) {
		if(fw1 == fw2) {
			// Bad combination of orientations
			return PE_ALS_DISCORD;
		}
		oneLeft = !fw1;
	}
	// Calc implied fragment size
	int64_t fraglo = min<int64_t>(off1, off2);
	int64_t fraghi = max<int64_t>(off1+len1, off2+len2);
	assert_gt(fraghi, fraglo);
	size_t frag = (size_t)(fraghi - fraglo);
	if(frag > maxfrag || frag < minfrag) {
		// Pair is discordant by virtue of the extents
		return PE_ALS_DISCORD;
	}
	int64_t lo1 = off1;
	int64_t hi1 = off1 + len1 - 1;
	int64_t lo2 = off2;
	int64_t hi2 = off2 + len2 - 1;
	bool containment = false;
	// Check whether one mate entirely contains the other
	if((lo1 >= lo2 && hi1 <= hi2) ||
	   (lo2 >= lo1 && hi2 <= hi1))
	{
		containment = true;
	}
	int type = PE_ALS_NORMAL;
	// Check whether one mate overlaps the other
	bool olap = false;
	if((lo1 <= lo2 && hi1 >= lo2) ||
	   (lo1 <= hi2 && hi1 >= hi2) ||
	   containment)
	{
		// The mates overlap
		olap = true;
		if(!olapOk_) return PE_ALS_DISCORD;
		type = PE_ALS_OVERLAP;
	}
	// Check if the mates are in the wrong relative orientation,
	// without any overlap
	if(!olap) {
		if((oneLeft && lo2 < lo1) || (!oneLeft && lo1 < lo2)) {
			return PE_ALS_DISCORD;
		}
	}
	// If one mate contained the other, report that
	if(containment) {
		if(!containOk_) return PE_ALS_DISCORD;
		type = PE_ALS_CONTAIN;
	}
	// Check whether there's dovetailing; i.e. does the left mate
	// extend past the right end of the right mate, or vice versa
	bool dovetailing = false;
	if(( oneLeft && (hi1 > hi2 || lo2 < lo1)) ||
	   (!oneLeft && (hi2 > hi1 || lo1 < lo2)))
	{
		dovetailing = true;
		if(!dovetailOk_) return PE_ALS_DISCORD;
		type = PE_ALS_DOVETAIL;
	}
	return type;
}

/**
 * Given details about how one mate aligns, and some details about the
 * reference sequence it aligned to, calculate a window and orientation s.t.
 * a paired-end alignment is concordant iff the opposite mate aligns in the
 * calculated window with the calculated orientation.  The "window" is really a
 * cosntraint on which positions the extreme end of the opposite mate can fall.
 * This is a different type of constraint from the one placed on seed-extend
 * dynamic programming problems.  That constraints requires that alignments at
 * one point pass through one of a set of "core" columns.
 *
 * When the opposite mate is to the left, we're constraining where its
 * left-hand extreme can fall, i.e., which cells in the top row of the matrix
 * it can end in.  When the opposite mate is to the right, we're cosntraining
 * where its right-hand extreme can fall, i.e., which cells in the bottom row
 * of the matrix it can end in.  However, in practice we can only constrain
 * where we start the backtrace, i.e. where the RHS of the alignment falls.
 * See frameFindMateRect for details.
 *
 * This calculaton does not consider gaps - the dynamic programming framer will
 * take gaps into account.
 *
 * Returns false if no concordant alignments are possible, true otherwise.
 */
bool PairedEndPolicy::otherMate(
	bool     is1,       // true -> mate 1 aligned and we're looking
					    // for 2, false -> vice versa
	bool     fw,        // orientation of aligned mate
	int64_t  off,       // offset into the reference sequence
	int64_t  maxalcols, // max # columns spanned by alignment
	size_t   reflen,    // length of reference sequence aligned to
	size_t   len1,      // length of mate 1
	size_t   len2,      // length of mate 2
	bool&    oleft,     // out: true iff opp mate must be to right of anchor
	int64_t& oll,       // out: leftmost Watson off for LHS of opp alignment
	int64_t& olr,       // out: rightmost Watson off for LHS of opp alignment
	int64_t& orl,       // out: leftmost Watson off for RHS of opp alignment
	int64_t& orr,       // out: rightmost Watson off for RHS of opp alignment
	bool&    ofw)       // out: true iff opp mate must be on Watson strand
	const
{
	assert_gt(len1, 0);
	assert_gt(len2, 0);
	assert_gt(maxfrag_, 0);
	assert_geq(minfrag_, 0);
	assert_geq(maxfrag_, minfrag_);
	assert(maxalcols == -1 || maxalcols > 0);
	
	// Calculate whether opposite mate should align to left or to right
	// of given mate, and what strand it should align to
	pePolicyMateDir(pol_, is1, fw, oleft, ofw);
	
	size_t alen = is1 ? len1 : len2; // length of opposite mate
	
	// Expand the maximum fragment length if necessary to accomodate
	// the longer mate
	size_t maxfrag = maxfrag_;
	size_t minfrag = minfrag_;
	if(minfrag < 1) {
		minfrag = 1;
	}
	if(len1 > maxfrag && expandToFit_) maxfrag = len1;
	if(len2 > maxfrag && expandToFit_) maxfrag = len2;
	if(!expandToFit_ && (len1 > maxfrag || len2 > maxfrag)) {
		// Not possible to find a concordant alignment; one of the
		// mates is too long
		return false;
	}
	
	// Now calculate bounds within which a dynamic programming
	// algorithm should search for an alignment for the opposite mate
	if(oleft) {
		//    -----------FRAG MAX----------------
		//                 -------FRAG MIN-------
		//                               |-alen-|
		//                             Anchor mate
		//                               ^off
		//                  |------|
		//       Not concordant: LHS not outside min
		//                 |------|
		//                Concordant
		//      |------|
		//     Concordant
		//  |------|
		// Not concordant: LHS outside max
		
		//    -----------FRAG MAX----------------
		//                 -------FRAG MIN-------
		//                               |-alen-|
		//                             Anchor mate
		//                               ^off
		//    |------------|
		// LHS can't be outside this range
		//                               -----------FRAG MAX----------------
		//    |------------------------------------------------------------|
		// LHS can't be outside this range, assuming no restrictions on
		// flipping, dovetailing, containment, overlap, etc.
		//                                      |-------|
		//                                      maxalcols
		//    |-----------------------------------------|
		// LHS can't be outside this range, assuming no flipping
		//    |---------------------------------|
		// LHS can't be outside this range, assuming no dovetailing
		//    |-------------------------|
		// LHS can't be outside this range, assuming no overlap

		oll = off + alen - maxfrag;
		olr = off + alen - minfrag;
		assert_geq(olr, oll);
		
		orl = oll;
		orr = off + maxfrag - 1;
		assert_geq(olr, oll);

		// What if overlapping alignments are not allowed?
		if(!olapOk_) {
			// RHS can't be flush with or to the right of off
			orr = min<int64_t>(orr, off-1);
			if(orr < olr) olr = orr;
			assert_leq(oll, olr);
			assert_leq(orl, orr);
			assert_geq(orr, olr);
		}
		// What if dovetail alignments are not allowed?
		else if(!dovetailOk_) {
			// RHS can't be past off+alen-1
			orr = min<int64_t>(orr, off + alen - 1);
			assert_leq(oll, olr);
			assert_leq(orl, orr);
		}
		// What if flipped alignments are not allowed?
		else if(!flippingOk_ && maxalcols != -1) {
			// RHS can't be right of ???
			orr = min<int64_t>(orr, off + alen - 1 + (maxalcols-1));
			assert_leq(oll, olr);
			assert_leq(orl, orr);
		}
		assert_geq(olr, oll);
		assert_geq(orr, orl);
		assert_geq(orr, olr);
		assert_geq(orl, oll);
	} else {
		//                             -----------FRAG MAX----------------
		//                             -------FRAG MIN-------
		//  -----------FRAG MAX----------------
		//                             |-alen-|
		//                           Anchor mate
		//                             ^off
 		//                                          |------|
		//                            Not concordant: RHS not outside min
		//                                           |------|
		//                                          Concordant
		//                                                      |------|
		//                                                     Concordant
		//                                                          |------|
		//                                      Not concordant: RHS outside max
		//

		//                             -----------FRAG MAX----------------
		//                             -------FRAG MIN-------
		//  -----------FRAG MAX----------------
		//                             |-alen-|
		//                           Anchor mate
		//                             ^off
		//                                                  |------------|
		//                                      RHS can't be outside this range
		//  |------------------------------------------------------------|
		// LHS can't be outside this range, assuming no restrictions on
		// dovetailing, containment, overlap, etc.
		//                     |-------|
		//                     maxalcols
		//                     |-----------------------------------------|
		//             LHS can't be outside this range, assuming no flipping
		//                             |---------------------------------|
		//          LHS can't be outside this range, assuming no dovetailing
		//                                     |-------------------------|
		//              LHS can't be outside this range, assuming no overlap
		
		orr = off + (maxfrag - 1);
		orl  = off + (minfrag - 1);
		assert_geq(orr, orl);
		
		oll = off + alen - maxfrag;
		olr = orr;
		assert_geq(olr, oll);
		
		// What if overlapping alignments are not allowed?
		if(!olapOk_) {
			// LHS can't be left of off+alen
			oll = max<int64_t>(oll, off+alen);
			if(oll > orl) orl = oll;
			assert_leq(oll, olr);
			assert_leq(orl, orr);
			assert_geq(orl, oll);
		}
		// What if dovetail alignments are not allowed?
		else if(!dovetailOk_) {
			// LHS can't be left of off
			oll = max<int64_t>(oll, off);
			assert_leq(oll, olr);
			assert_leq(orl, orr);
		}
		// What if flipped alignments are not allowed?
		else if(!flippingOk_ && maxalcols != -1) {
			// LHS can't be left of off - maxalcols + 1
			oll = max<int64_t>(oll, off - maxalcols + 1);
			assert_leq(oll, olr);
			assert_leq(orl, orr);
		}
		assert_geq(olr, oll);
		assert_geq(orr, orl);
		assert_geq(orr, olr);
		assert_geq(orl, oll);
	}

	// Boundaries and orientation determined
	return true;
}

#ifdef MAIN_PE

#include <string>
#include <sstream>

void testCaseClassify(
	const string& name,
	int      pol,
	size_t   maxfrag,
	size_t   minfrag,
	bool     local,
	bool     flip,
	bool     dove,
	bool     cont,
	bool     olap,
	bool     expand,
	int64_t  off1,
	size_t   len1,
	bool     fw1,
	int64_t  off2,
	size_t   len2,
	bool     fw2,
	int      expect_class)
{
	PairedEndPolicy pepol(
		pol,
		maxfrag,
		minfrag,
		local,
		flip,
		dove,
		cont,
		olap,
		expand);
	int ret = pepol.peClassifyPair(
		off1,   // offset of mate 1
		len1,   // length of mate 1
		fw1,    // whether mate 1 aligned to Watson
		off2,   // offset of mate 2
		len2,   // length of mate 2
		fw2);   // whether mate 2 aligned to Watson
	assert_eq(expect_class, ret);
	cout << "peClassifyPair: " << name << "...PASSED" << endl;
}

void testCaseOtherMate(
	const string& name,
	int      pol,
	size_t   maxfrag,
	size_t   minfrag,
	bool     local,
	bool     flip,
	bool     dove,
	bool     cont,
	bool     olap,
	bool     expand,
	bool     is1,
	bool     fw,
	int64_t  off,
	int64_t  maxalcols,
	size_t   reflen,
	size_t   len1,
	size_t   len2,
	bool     expect_ret,
	bool     expect_oleft,
	int64_t  expect_oll,
	int64_t  expect_olr,
	int64_t  expect_orl,
	int64_t  expect_orr,
	bool     expect_ofw)
{
	PairedEndPolicy pepol(
		pol,
		maxfrag,
		minfrag,
		local,
		flip,
		dove,
		cont,
		olap,
		expand);
	int64_t oll = 0, olr = 0;
	int64_t orl = 0, orr = 0;
	bool oleft = false, ofw = false;
	bool ret = pepol.otherMate(
		is1,
		fw,
		off,
		maxalcols,
		reflen,
		len1,
		len2,
		oleft,
		oll,
		olr,
		orl,
		orr,
		ofw);
	assert(ret == expect_ret);
	if(ret) {
		assert_eq(expect_oleft, oleft);
		assert_eq(expect_oll, oll);
		assert_eq(expect_olr, olr);
		assert_eq(expect_orl, orl);
		assert_eq(expect_orr, orr);
		assert_eq(expect_ofw, ofw);
	}
	cout << "otherMate: " << name << "...PASSED" << endl;
}

int main(int argc, char **argv) {

	// Set of 8 cases where we look for the opposite mate to the right
	// of the anchor mate, with various combinations of policies and
	// anchor-mate orientations.

	// |--------|
	//           |--------|
	//           ^110     ^119
	// |------------------|
	//      min frag
	//                     |--------|
	//                     ^120     ^129
	// |----------------------------|
	//           max frag
	// ^
	// 100

	{
	int  policies[] = { PE_POLICY_FF, PE_POLICY_RR, PE_POLICY_FR, PE_POLICY_RF, PE_POLICY_FF, PE_POLICY_RR, PE_POLICY_FR, PE_POLICY_RF };
	bool is1[]      = { true,  true,   true,  true, false, false, false, false };
	bool fw[]       = { true,  false,  true, false, false,  true,  true, false };
	bool oleft[]    = { false, false, false, false, false, false, false, false };
	bool ofw[]      = { true,  false, false,  true, false,  true, false,  true };

	for(int i = 0; i < 8; i++) {
		ostringstream oss;
		oss << "Simple";
		oss << i;
		testCaseOtherMate(
			oss.str(),
			policies[i],  // policy
			30,           // maxfrag
			20,           // minfrag
			false,        // local
			true,         // flipping OK
			true,         // dovetail OK
			true,         // containment OK
			true,         // overlap OK
			true,         // expand-to-fit
			is1[i],       // mate 1 is anchor
			fw[i],        // anchor aligned to Watson
			100,          // anchor's offset into ref
			-1,           // max # alignment cols
			200,          // ref length
			10,           // mate 1 length
			10,           // mate 2 length
			true,         // expected return val from otherMate
			oleft[i],     // wheter to look for opposite to left
			80,           // expected leftmost pos for opp mate LHS
			129,          // expected rightmost pos for opp mate LHS
			119,          // expected leftmost pos for opp mate RHS
			129,          // expected rightmost pos for opp mate RHS
			ofw[i]);      // expected orientation in which opposite mate must align
	}
	}

	// Set of 8 cases where we look for the opposite mate to the left
	// of the anchor mate, with various combinations of policies and
	// anchor-mate orientations.

	// |--------|
	// ^100     ^109
	//           |--------|
	//           ^110     ^119
	//           |------------------|
	//                 min frag
	//                     |-Anchor-|
	//                     ^120     ^129
	// |----------------------------|
	//           max frag
	// ^
	// 100

	{
	int  policies[] = { PE_POLICY_FF, PE_POLICY_RR, PE_POLICY_FR, PE_POLICY_RF, PE_POLICY_FF, PE_POLICY_RR, PE_POLICY_FR, PE_POLICY_RF };
	bool is1[]      = { false, false, false, false,  true,  true,  true,  true };
	bool fw[]       = {  true, false, false,  true, false,  true, false,  true };
	bool oleft[]    = {  true,  true,  true,  true,  true,  true,  true,  true };
	bool ofw[]      = {  true, false,  true, false, false,  true,  true, false };
	
	for(int i = 0; i < 8; i++) {
		ostringstream oss;
		oss << "Simple";
		oss << (i+8);
		testCaseOtherMate(
			oss.str(),
			policies[i],  // policy
			30,           // maxfrag
			20,           // minfrag
			false,        // local
			true,         // flipping OK
			true,         // dovetail OK
			true,         // containment OK
			true,         // overlap OK
			true,         // expand-to-fit
			is1[i],       // mate 1 is anchor
			fw[i],        // anchor aligned to Watson
			120,          // anchor's offset into ref
			-1,           // max # alignment cols
			200,          // ref length
			10,           // mate 1 length
			10,           // mate 2 length
			true,         // expected return val from otherMate
			oleft[i],     // wheter to look for opposite to left
			100,          // expected leftmost pos for opp mate LHS
			110,          // expected rightmost pos for opp mate LHS
			100,          // expected leftmost pos for opp mate RHS
			149,          // expected rightmost pos for opp mate RHS
			ofw[i]);      // expected orientation in which opposite mate must align
	}
	}

	// Case where min frag == max frag and opposite is to the right

	// |----------------------------|
	//      min frag
	//                     |--------|
	//                     ^120     ^129
	// |----------------------------|
	//           max frag
	// ^
	// 100
	testCaseOtherMate(
		"MinFragEqMax1",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		30,           // minfrag
		false,        // local
		true,         // flipping OK
		true,         // dovetail OK
		true,         // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		false,        // mate 1 is anchor
		false,        // anchor aligned to Watson
		120,          // anchor's offset into ref
		-1,           // max # alignment cols
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		true,         // expected return val from otherMate
		true,         // wheter to look for opposite to left
		100,          // expected leftmost pos for opp mate LHS
		100,          // expected rightmost pos for opp mate LHS
		100,          // expected leftmost pos for opp mate RHS
		149,          // expected rightmost pos for opp mate RHS
		true);        // expected orientation in which opposite mate must align

	// Case where min frag == max frag and opposite is to the right

	// |----------------------------|
	//      min frag                ^129
	// |--------|
	// ^100     ^109
	// |----------------------------|
	//           max frag
	testCaseOtherMate(
		"MinFragEqMax2",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		30,           // minfrag
		false,        // local
		true,         // flipping OK
		true,         // dovetail OK
		true,         // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		true,         // mate 1 is anchor
		true,         // anchor aligned to Watson
		100,          // anchor's offset into ref
		-1,           // max # alignment cols
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		true,         // expected return val from otherMate
		false,        // wheter to look for opposite to left
		80,           // expected leftmost pos for opp mate LHS
		129,          // expected rightmost pos for opp mate LHS
		129,          // expected leftmost pos for opp mate RHS
		129,          // expected rightmost pos for opp mate RHS
		false);       // expected orientation in which opposite mate must align

	testCaseOtherMate(
		"MinFragEqMax4NoDove1",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		25,           // minfrag
		false,        // local
		true,         // flipping OK
		false,        // dovetail OK
		true,         // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		true,         // mate 1 is anchor
		true,         // anchor aligned to Watson
		100,          // anchor's offset into ref
		-1,           // max # alignment cols
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		true,         // expected return val from otherMate
		false,        // wheter to look for opposite to left
		100,          // expected leftmost pos for opp mate LHS
		129,          // expected rightmost pos for opp mate LHS
		124,          // expected leftmost pos for opp mate RHS
		129,          // expected rightmost pos for opp mate RHS
		false);       // expected orientation in which opposite mate must align

	testCaseOtherMate(
		"MinFragEqMax4NoCont1",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		25,           // minfrag
		false,        // local
		true,         // flipping OK
		false,        // dovetail OK
		false,        // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		true,         // mate 1 is anchor
		true,         // anchor aligned to Watson
		100,          // anchor's offset into ref
		-1,           // max # alignment cols
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		true,         // expected return val from otherMate
		false,        // wheter to look for opposite to left
		100,          // expected leftmost pos for opp mate LHS
		129,          // expected rightmost pos for opp mate LHS
		124,          // expected leftmost pos for opp mate RHS
		129,          // expected rightmost pos for opp mate RHS
		false);       // expected orientation in which opposite mate must align

	testCaseOtherMate(
		"MinFragEqMax4NoOlap1",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		25,           // minfrag
		false,        // local
		true,         // flipping OK
		false,        // dovetail OK
		false,        // containment OK
		false,        // overlap OK
		true,         // expand-to-fit
		true,         // mate 1 is anchor
		true,         // anchor aligned to Watson
		100,          // anchor's offset into ref
		-1,           // max # alignment cols
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		true,         // expected return val from otherMate
		false,        // wheter to look for opposite to left
		110,          // expected leftmost pos for opp mate LHS
		129,          // expected rightmost pos for opp mate LHS
		124,          // expected leftmost pos for opp mate RHS
		129,          // expected rightmost pos for opp mate RHS
		false);       // expected orientation in which opposite mate must align

	testCaseOtherMate(
		"MinFragEqMax4NoDove2",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		25,           // minfrag
		false,        // local
		true,         // flipping OK
		false,        // dovetail OK
		true,         // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		false,        // mate 1 is anchor
		false,        // anchor aligned to Watson
		120,          // anchor's offset into ref
		-1,           // max # alignment cols
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		true,         // expected return val from otherMate
		true,         // whether to look for opposite to left
		100,          // expected leftmost pos for opp mate LHS
		105,          // expected rightmost pos for opp mate LHS
		100,          // expected leftmost pos for opp mate RHS
		129,          // expected rightmost pos for opp mate RHS
		true);        // expected orientation in which opposite mate must align

	testCaseOtherMate(
		"MinFragEqMax4NoOlap2",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		25,           // minfrag
		false,        // local
		true,         // flipping OK
		false,        // dovetail OK
		false,        // containment OK
		false,        // overlap OK
		true,         // expand-to-fit
		false,        // mate 1 is anchor
		false,        // anchor aligned to Watson
		120,          // anchor's offset into ref
		-1,           // max # alignment cols
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		true,         // expected return val from otherMate
		true,         // whether to look for opposite to left
		100,          // expected leftmost pos for opp mate LHS
		105,          // expected rightmost pos for opp mate LHS
		100,          // expected leftmost pos for opp mate RHS
		119,          // expected rightmost pos for opp mate RHS
		true);        // expected orientation in which opposite mate must align

	{
	int olls[] = { 110 };
	int olrs[] = { 299 };
	int orls[] = { 149 };
	int orrs[] = { 299 };
	for(int i = 0; i < 1; i++) {
		ostringstream oss;
		oss << "Overhang1_";
		oss << (i+1);
		testCaseOtherMate(
			oss.str(),
			PE_POLICY_FR, // policy
			200,          // maxfrag
			50,           // minfrag
			false,        // local
			true,         // flipping OK
			true,         // dovetail OK
			true,         // containment OK
			false,        // overlap OK
			true,         // expand-to-fit
			true,         // mate 1 is anchor
			true,         // anchor aligned to Watson
			100,          // anchor's offset into ref
			-1,           // max # alignment cols
			200,          // ref length
			10,           // mate 1 length
			10,           // mate 2 length
			true,         // expected return val from otherMate
			false,        // whether to look for opposite to left
			olls[i],      // expected leftmost pos for opp mate LHS
			olrs[i],      // expected rightmost pos for opp mate LHS
			orls[i],      // expected leftmost pos for opp mate RHS
			orrs[i],      // expected rightmost pos for opp mate RHS
			false);       // expected orientation in which opposite mate must align
	}
	}

	{
	int olls[] = { -100 };
	int olrs[] = {   50 };
	int orls[] = { -100 };
	int orrs[] = {   89 };
	for(int i = 0; i < 1; i++) {
		ostringstream oss;
		oss << "Overhang2_";
		oss << (i+1);
		testCaseOtherMate(
			oss.str(),
			PE_POLICY_FR, // policy
			200,          // maxfrag
			50,           // minfrag
			false,        // local
			true,         // flipping OK
			true,         // dovetail OK
			true,         // containment OK
			false,        // overlap OK
			true,         // expand-to-fit
			true,         // mate 1 is anchor
			false,        // anchor aligned to Watson
			90,           // anchor's offset into ref
			-1,           // max # alignment cols
			200,          // ref length
			10,           // mate 1 length
			10,           // mate 2 length
			true,         // expected return val from otherMate
			true,         // whether to look for opposite to left
			olls[i],      // expected leftmost pos for opp mate LHS
			olrs[i],      // expected rightmost pos for opp mate LHS
			orls[i],      // expected leftmost pos for opp mate RHS
			orrs[i],      // expected rightmost pos for opp mate RHS
			true);        // expected orientation in which opposite mate must align
	}
	}

	{
	int mate2offs[] = {           150,            149,            149,            100,              99,           299,              1,            250,            250 };
	int mate2lens[] = {            50,             50,             51,            100,             101,             1,             50,             50,             51 };
	int peExpects[] = { PE_ALS_NORMAL, PE_ALS_DISCORD, PE_ALS_OVERLAP, PE_ALS_CONTAIN, PE_ALS_DOVETAIL, PE_ALS_NORMAL, PE_ALS_DISCORD,  PE_ALS_NORMAL, PE_ALS_DISCORD };

	for(int i = 0; i < 9; i++) {
		ostringstream oss;
		oss << "Simple1_";
		oss << (i);
		testCaseClassify(
			oss.str(),
			PE_POLICY_FR, // policy
			200,          // maxfrag
			100,          // minfrag
			false,        // local
			true,         // flipping OK
			true,         // dovetail OK
			true,         // containment OK
			true,         // overlap OK
			true,         // expand-to-fit
			100,          // offset of mate 1
			50,           // length of mate 1
			true,         // whether mate 1 aligned to Watson
			mate2offs[i], // offset of mate 2
			mate2lens[i], // length of mate 2
			false,        // whether mate 2 aligned to Watson
			peExpects[i]);// expectation for PE_ALS flag returned
	}
	}

	{
	int mate1offs[] = {           200,            201,            200,            200,             200,           100,            400,            100,             99 };
	int mate1lens[] = {            50,             49,             51,            100,             101,             1,             50,             50,             51 };
	int peExpects[] = { PE_ALS_NORMAL, PE_ALS_DISCORD, PE_ALS_OVERLAP, PE_ALS_CONTAIN, PE_ALS_DOVETAIL, PE_ALS_NORMAL, PE_ALS_DISCORD,  PE_ALS_NORMAL, PE_ALS_DISCORD };

	for(int i = 0; i < 9; i++) {
		ostringstream oss;
		oss << "Simple2_";
		oss << (i);
		testCaseClassify(
			oss.str(),
			PE_POLICY_FR, // policy
			200,          // maxfrag
			100,          // minfrag
			false,        // local
			true,         // flipping OK
			true,         // dovetail OK
			true,         // containment OK
			true,         // overlap OK
			true,         // expand-to-fit
			mate1offs[i], // offset of mate 1
			mate1lens[i], // length of mate 1
			true,         // whether mate 1 aligned to Watson
			250,          // offset of mate 2
			50,           // length of mate 2
			false,        // whether mate 2 aligned to Watson
			peExpects[i]);// expectation for PE_ALS flag returned
	}
	}

	testCaseOtherMate(
		"Regression1",
		PE_POLICY_FF, // policy
		50,           // maxfrag
		0,            // minfrag
		false,        // local
		true,         // flipping OK
		true,         // dovetail OK
		true,         // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		true,         // mate 1 is anchor
		false,        // anchor aligned to Watson
		3,            // anchor's offset into ref
		-1,           // max # alignment cols
		53,           // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		true,         // expected return val from otherMate
		true,         // whether to look for opposite to left
		-37,          // expected leftmost pos for opp mate LHS
		13,           // expected rightmost pos for opp mate LHS
		-37,          // expected leftmost pos for opp mate RHS
		52,           // expected rightmost pos for opp mate RHS
		false);       // expected orientation in which opposite mate must align
}

#endif /*def MAIN_PE*/
