/*
 *  pe.cpp
 */

#include <algorithm>
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
	uint32_t len1,   // length of mate 1
	bool     fw1,    // whether mate 1 aligned to Watson
	int64_t  off2,   // offset of mate 2
	uint32_t len2,   // length of mate 2
	bool     fw2)    // whether mate 2 aligned to Watson
	const
{
	assert_gt(len1, 0);
	assert_gt(len2, 0);
	// Expand the maximum fragment length if necessary to accomodate
	// the longer mate
	uint32_t maxfrag = maxfrag_;
	if(len1 > maxfrag && expandToFit_) maxfrag = len1;
	if(len2 > maxfrag && expandToFit_) maxfrag = len2;
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
	uint32_t frag = (uint32_t)(fraghi - fraglo);
	if(frag > maxfrag || frag < minfrag_) {
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
	if(lo1 <= lo2 && hi1 >= lo2 ||
	   lo1 <= hi2 && hi1 >= hi2 ||
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
 * Given details about how one mate aligns, and some details about
 * the reference sequence it aligned to, calculate a window and
 * orientation s.t. the alignment for the pair will be concordant
 * if the other mate aligns with that orientation in that window.
 *
 * Returns false if it is clearly not possible to find a concordant
 * alignment involving these mates, true otherwise.
 */
bool PairedEndPolicy::otherMate(
	bool     is1,       // true -> mate 1 aligned and we're looking
					    // for 2, false -> vice versa
	bool     fw,        // orientation of aligned mate
	int64_t  off,       // offset into the reference sequence
	uint32_t reflen,    // length of reference sequence aligned to
	uint32_t len1,      // length of mate 1
	uint32_t len2,      // length of mate 2
	int      maxgaps,   // maximum number of gaps permitted in the
	                    // alignment of the opposite mate
	int      maxohang,  // maximum overhang of dynamic programming
	                    // region off end of reference
	bool&    oleft,     // whether to look to the left for opposite mate
	int64_t& oleftoff,  // offset of leftmost character to include in
	                    // dyn prog problem looking for opposite mate
	int64_t& orightoff, // offset of rightmost character to include in
	                    // dyn prog problem looking for opposite mate
	bool&    ofw)       // whether to look for opposite mate's forward
	const               // or reverse-comp representation
{
	assert_gt(maxfrag_, 0);
	assert_geq(minfrag_, 0);
	assert_geq(maxfrag_, minfrag_);
	
	// Calculate whether opposite mate should align to left or to right
	// of given mate, and what strand it should align to
	pePolicyMateDir(pol_, is1, fw, oleft, ofw);
	
	uint32_t olen = is1 ? len2 : len1; // length of anchor mate
	uint32_t alen = is1 ? len1 : len2; // length of opposite mate
	
	// Expand the maximum fragment length if necessary to accomodate
	// the longer mate
	uint32_t maxfrag = maxfrag_;
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
		// Opposite mate is to the left.  Left-hand side of opposite
		// mate has to be inside Frag Max and outside Frag Min
		
		//    --------------FRAG MAX-------------------
		//                       -------FRAG MIN-------
		//                                     |------|
		//                                    Anchor mate
		//                        |------|
		//            Not concordant: LHS not outside min
		//             |------|
		//            Concordant
		//      |------|
		//     Concordant
		//  |------|
		// Not concordant: RHS outside max

		oleftoff  = off + alen - maxfrag;
		orightoff = off + alen - minfrag_;
		assert_geq(orightoff, oleftoff);

		//    --------------FRAG MAX-------------------
		//                       -------FRAG MIN-------
		//    ^oleftoff          ^orightoff
		//                                     |------|
		//                                    Anchor mate
		
		orightoff += (olen - 1);

		//    --------------FRAG MAX-------------------
		//                       -------FRAG MIN-------
		//    ^oleftoff                    ^orightoff
		//                                     |------|
		//                                    Anchor mate
		//                       |---------|
		//                     Opposite mate len

		// Add 'maxgaps' to orightoff to account for additional
		// alignment length owing to gaps
		orightoff += maxgaps;

		//    --------------FRAG MAX-------------------
		//                       -------FRAG MIN-------
		//    ^oleftoff                         ^orightoff
		//                                     |------|
		//                                    Anchor mate
		//                       |---------|
		//                     Opposite mate len
		//                                  |---|
		//                                 Max gaps

		// What if orightoff is now so far off the right-hand side of
		// the anchor mate that the frag-max would be exceeded if an
		// alignment really started there?
		int64_t extent = orightoff - off;
		if(extent > maxfrag) {
			orightoff -= (extent - maxfrag);
			assert_leq(oleftoff, orightoff);
		}

		// What if overlapping alignments are not allowed?
		if(!olapOk_) {
			// oleftoff can't be less than off
			if(orightoff >= off) {
				orightoff = off - 1;
				assert_leq(oleftoff, orightoff);
			}
		}
		// What if contained alignments are not allowed?
		else if(!containOk_) {
			// oleftoff can't be less than or equal to off
			if(orightoff >= off+alen-1) {
				orightoff = off+alen-2;
				assert_leq(oleftoff, orightoff);
			}
			// can we move oleftoff even further to the right?
		}
		// What if dovetail alignments are not allowed?
		else if(!dovetailOk_) {
			// orightoff can't be off RHS of anchor
			if(orightoff >= off + alen) {
				orightoff = off + alen - 1;
				assert_leq(oleftoff, orightoff);
			}
		}

		// Enforce overhang limit
		assert_lt (oleftoff,  reflen);
		assert_geq(orightoff, 0);
		if(oleftoff < -maxohang)                oleftoff  = -maxohang;
		if((orightoff - reflen + 1) > maxohang) orightoff = reflen + maxohang - 1;

		// What if the window that remains is too small to contain an
		// alignment for the opposite mate, even with a bunch of
		// reference gaps?
		if(!local_ && (orightoff - oleftoff + 1) < ((int64_t)olen - maxgaps)) {
			// Not possible to find a concordant alignment; the window
			// is too small for the opposite mate to fit in
			return false;
		}
	} else {
		// Opposite mate is to the right.  Right-hand side of opposite
		// mate has to be inside Frag Max and outside Frag Min
		
		//  --------------FRAG MAX-------------------
		//  -------FRAG MIN-------
		//  |------|
		// Anchor mate
		//               |------|
		//    Not concordant: RHS not outside min
		//                  |------|
		//                 Concordant
		//                                   |------|
		//                                  Concordant
		//                                     |------|
		//                          Not concordant: RHS outside max
		//
		
		orightoff = off + (maxfrag - 1);
		oleftoff  = off + (minfrag_ - 1);
		assert_geq(orightoff, oleftoff);
		
		//  --------------FRAG MAX-------------------
		//  -------FRAG MIN-------
		//                       ^oleftoff          ^orightoff
		//  |------|
		// Anchor mate
		
		// It is possible that oleftoff is within the anchor mate, but
		// it can't be to the left of the anchor mate
		assert_geq(oleftoff, off);

		// Subtract from oleftoff to account for the length of the
		// opposite mate
		oleftoff -= (olen - 1);

		//  --------------FRAG MAX-------------------
		//  -------FRAG MIN-------
		//             ^oleftoff                    ^orightoff
		//  |------|
		// Anchor mate
		//             |---------|
		//           Opposite mate len
		
		// Subtract 'maxgaps' from oleftoff to account for additional
		// alignment length owing to gaps
		oleftoff -= maxgaps;

		//  --------------FRAG MAX-------------------
		//  -------FRAG MIN-------
		//        ^oleftoff                         ^orightoff
		//  |------|
		// Anchor mate
		//             |---------|
		//           Opposite mate len
		//        |---|
		//       Max gaps
		
		// What if oleftoff is now so far off the left-hand side of the
		// anchor mate that the frag-max would be exceeded if an
		// alignment really started there?
		int64_t extent = off + alen - oleftoff;
		if(extent > maxfrag) {
			oleftoff += (extent - maxfrag);
			assert_leq(oleftoff, orightoff);
		}
		
		// What if overlapping alignments are not allowed?
		if(!olapOk_) {
			// oleftoff can't be less than off
			if(oleftoff < off+alen) {
				oleftoff = off+alen;
				assert_leq(oleftoff, orightoff);
			}
		}
		// What if contained alignments are not allowed?
		else if(!containOk_) {
			// oleftoff can't be less than or equal to off
			if(oleftoff <= off) {
				oleftoff = off+1;
				assert_leq(oleftoff, orightoff);
			}
			// can we move oleftoff even further to the right?
		}
		// What if dovetail alignments are not allowed?
		else if(!dovetailOk_) {
			// oleftoff can't be less than off
			if(oleftoff < off) {
				oleftoff = off;
				assert_leq(oleftoff, orightoff);
			}
		}

		// Enforce overhang limit
		assert_lt (oleftoff,  reflen);
		assert_geq(orightoff, 0);
		if(oleftoff < -maxohang)                oleftoff  = -maxohang;
		if((orightoff - reflen + 1) > maxohang) orightoff = reflen + maxohang - 1;

		// What if the window that remains is too small to contain an
		// alignment for the opposite mate, even with a bunch of
		// reference gaps?
		if(!local_ && (orightoff - oleftoff + 1) < ((int64_t)olen - maxgaps)) {
			// Not possible to find a concordant alignment; the window
			// is too small for the opposite mate to fit in
			return false;
		}
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
	uint32_t maxfrag,
	uint32_t minfrag,
	bool     local,
	bool     dove,
	bool     cont,
	bool     olap,
	bool     expand,
	int64_t  off1,
	uint32_t len1,
	bool     fw1,
	int64_t  off2,
	uint32_t len2,
	bool     fw2,
	int      expect_class)
{
	PairedEndPolicy pepol(
		pol,
		maxfrag,
		minfrag,
		local,
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
	uint32_t maxfrag,
	uint32_t minfrag,
	bool     local,
	bool     dove,
	bool     cont,
	bool     olap,
	bool     expand,
	bool     is1,
	bool     fw,
	uint32_t off,
	uint32_t reflen,
	uint32_t len1,
	uint32_t len2,
	int      maxgaps,
	int      maxohang,
	bool     expect_ret,
	bool     expect_oleft,
	int64_t  expect_oleftoff,
	int64_t  expect_orightoff,
	bool     expect_ofw)
{
	PairedEndPolicy pepol(
		pol,
		maxfrag,
		minfrag,
		local,
		dove,
		cont,
		olap,
		expand);
	int64_t oleftoff = 0, orightoff = 0;
	bool oleft = false, ofw = false;
	bool ret = pepol.otherMate(
		is1,
		fw,
		off,
		reflen,
		len1,
		len2,
		maxgaps,
		maxohang,
		oleft,
		oleftoff,
		orightoff,
		ofw);
	assert(ret == expect_ret);
	if(ret) {
		assert_eq(expect_oleft, oleft);
		assert_eq(expect_oleftoff, oleftoff);
		assert_eq(expect_orightoff, orightoff);
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
			true,         // dovetail OK
			true,         // containment OK
			true,         // overlap OK
			true,         // expand-to-fit
			is1[i],       // mate 1 is anchor
			fw[i],        // anchor aligned to Watson
			100,          // anchor's offset into ref
			200,          // ref length
			10,           // mate 1 length
			10,           // mate 2 length
			0,            // maximum # gaps in alignment for opposite mate
			0,            // maximum overhang off end of reference
			true,         // expected return val from otherMate
			oleft[i],     // wheter to look for opposite to left
			110,          // expected left extent of dynamic prog. problem
			129,          // expected right extent of dynamic prog. problem
			ofw[i]);      // expected orientation in which opposite mate must align
	}
	}

	// Set of 8 cases where we look for the opposite mate to the right
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
			true,         // dovetail OK
			true,         // containment OK
			true,         // overlap OK
			true,         // expand-to-fit
			is1[i],       // mate 1 is anchor
			fw[i],        // anchor aligned to Watson
			120,          // anchor's offset into ref
			200,          // ref length
			10,           // mate 1 length
			10,           // mate 2 length
			0,            // maximum # gaps in alignment for opposite mate
			0,            // maximum overhang off end of reference
			true,         // expected return val from otherMate
			oleft[i],     // wheter to look for opposite to left
			100,          // expected left extent of dynamic prog. problem
			119,          // expected right extent of dynamic prog. problem
			ofw[i]);      // expected orientation in which opposite mate must align
	}
	}

	// Case where min frag == max frag and opposite is to the left

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
		true,         // dovetail OK
		true,         // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		false,        // mate 1 is anchor
		false,        // anchor aligned to Watson
		120,          // anchor's offset into ref
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		0,            // maximum # gaps in alignment for opposite mate
		0,            // maximum overhang off end of reference
		true,         // expected return val from otherMate
		true,         // wheter to look for opposite to left
		100,          // expected left extent of dynamic prog. problem
		109,          // expected right extent of dynamic prog. problem
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
		true,         // dovetail OK
		true,         // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		true,         // mate 1 is anchor
		true,         // anchor aligned to Watson
		100,          // anchor's offset into ref
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		0,            // maximum # gaps in alignment for opposite mate
		0,            // maximum overhang off end of reference
		true,         // expected return val from otherMate
		false,        // wheter to look for opposite to left
		120,          // expected left extent of dynamic prog. problem
		129,          // expected right extent of dynamic prog. problem
		false);       // expected orientation in which opposite mate must align

	// Case where min frag == max frag, opposite is to the right and gap len is 2

	// |----------------------------|
	//      min frag                ^129
	// |--------|
	// ^100     ^109
	// |----------------------------|
	//           max frag
	testCaseOtherMate(
		"MinFragEqMax3",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		30,           // minfrag
		false,        // local
		true,         // dovetail OK
		true,         // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		true,         // mate 1 is anchor
		true,         // anchor aligned to Watson
		100,          // anchor's offset into ref
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		2,            // maximum # gaps in alignment for opposite mate
		0,            // maximum overhang off end of reference
		true,         // expected return val from otherMate
		false,        // wheter to look for opposite to left
		118,          // expected left extent of dynamic prog. problem
		129,          // expected right extent of dynamic prog. problem
		false);       // expected orientation in which opposite mate must align


	// Case where min frag == max frag, opposite is to the right and gap len is huge

	// |----------------------------|
	//      min frag                ^129
	// |--------|
	// ^100     ^109
	// |----------------------------|
	//           max frag
	testCaseOtherMate(
		"MinFragEqMax4",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		30,           // minfrag
		false,        // local
		true,         // dovetail OK
		true,         // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		true,         // mate 1 is anchor
		true,         // anchor aligned to Watson
		100,          // anchor's offset into ref
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		300,          // maximum # gaps in alignment for opposite mate
		0,            // maximum overhang off end of reference
		true,         // expected return val from otherMate
		false,        // wheter to look for opposite to left
		80,           // expected left extent of dynamic prog. problem
		129,          // expected right extent of dynamic prog. problem
		false);       // expected orientation in which opposite mate must align

	testCaseOtherMate(
		"MinFragEqMax4NoDove1",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		25,           // minfrag
		false,        // local
		false,        // dovetail OK
		true,         // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		true,         // mate 1 is anchor
		true,         // anchor aligned to Watson
		100,          // anchor's offset into ref
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		300,          // maximum # gaps in alignment for opposite mate
		0,            // maximum overhang off end of reference
		true,         // expected return val from otherMate
		false,        // wheter to look for opposite to left
		100,          // expected left extent of dynamic prog. problem
		129,          // expected right extent of dynamic prog. problem
		false);       // expected orientation in which opposite mate must align

	testCaseOtherMate(
		"MinFragEqMax4NoCont1",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		25,           // minfrag
		false,        // local
		false,        // dovetail OK
		false,        // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		true,         // mate 1 is anchor
		true,         // anchor aligned to Watson
		100,          // anchor's offset into ref
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		300,          // maximum # gaps in alignment for opposite mate
		0,            // maximum overhang off end of reference
		true,         // expected return val from otherMate
		false,        // wheter to look for opposite to left
		101,          // expected left extent of dynamic prog. problem
		129,          // expected right extent of dynamic prog. problem
		false);       // expected orientation in which opposite mate must align

	testCaseOtherMate(
		"MinFragEqMax4NoOlap1",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		25,           // minfrag
		false,        // local
		false,        // dovetail OK
		false,        // containment OK
		false,        // overlap OK
		true,         // expand-to-fit
		true,         // mate 1 is anchor
		true,         // anchor aligned to Watson
		100,          // anchor's offset into ref
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		300,          // maximum # gaps in alignment for opposite mate
		0,            // maximum overhang off end of reference
		true,         // expected return val from otherMate
		false,        // wheter to look for opposite to left
		110,          // expected left extent of dynamic prog. problem
		129,          // expected right extent of dynamic prog. problem
		false);       // expected orientation in which opposite mate must align

	testCaseOtherMate(
		"MinFragEqMax4NoDove2",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		25,           // minfrag
		false,        // local
		false,        // dovetail OK
		true,         // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		false,        // mate 1 is anchor
		false,        // anchor aligned to Watson
		120,          // anchor's offset into ref
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		300,          // maximum # gaps in alignment for opposite mate
		0,            // maximum overhang off end of reference
		true,         // expected return val from otherMate
		true,         // whether to look for opposite to left
		100,          // expected left extent of dynamic prog. problem
		129,          // expected right extent of dynamic prog. problem
		true);        // expected orientation in which opposite mate must align

	testCaseOtherMate(
		"MinFragEqMax4NoCont2",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		25,           // minfrag
		false,        // local
		false,        // dovetail OK
		false,        // containment OK
		true,         // overlap OK
		true,         // expand-to-fit
		false,        // mate 1 is anchor
		false,        // anchor aligned to Watson
		120,          // anchor's offset into ref
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		300,          // maximum # gaps in alignment for opposite mate
		0,            // maximum overhang off end of reference
		true,         // expected return val from otherMate
		true,         // whether to look for opposite to left
		100,          // expected left extent of dynamic prog. problem
		128,          // expected right extent of dynamic prog. problem
		true);        // expected orientation in which opposite mate must align

	testCaseOtherMate(
		"MinFragEqMax4NoOlap2",
		PE_POLICY_FR, // policy
		30,           // maxfrag
		25,           // minfrag
		false,        // local
		false,        // dovetail OK
		false,        // containment OK
		false,        // overlap OK
		true,         // expand-to-fit
		false,        // mate 1 is anchor
		false,        // anchor aligned to Watson
		120,          // anchor's offset into ref
		200,          // ref length
		10,           // mate 1 length
		10,           // mate 2 length
		300,          // maximum # gaps in alignment for opposite mate
		0,            // maximum overhang off end of reference
		true,         // expected return val from otherMate
		true,         // whether to look for opposite to left
		100,          // expected left extent of dynamic prog. problem
		119,          // expected right extent of dynamic prog. problem
		true);        // expected orientation in which opposite mate must align

	{
	int ohang[]    = {   0,   1,   0,   1,   2 };
	int maxgaps[]  = {   0,   0,   1,   1,   1 };
	int leftExt[]  = { 140, 140, 139, 139, 139 };
	int rightExt[] = { 199, 200, 199, 200, 201 };
	for(int i = 0; i < 5; i++) {
		ostringstream oss;
		oss << "Overhang1_";
		oss << (i+1);
		testCaseOtherMate(
			oss.str(),
			PE_POLICY_FR, // policy
			200,          // maxfrag
			50,           // minfrag
			false,        // local
			true,         // dovetail OK
			true,         // containment OK
			false,        // overlap OK
			true,         // expand-to-fit
			true,         // mate 1 is anchor
			true,         // anchor aligned to Watson
			100,          // anchor's offset into ref
			200,          // ref length
			10,           // mate 1 length
			10,           // mate 2 length
			maxgaps[i],   // maximum # gaps in alignment for opposite mate
			ohang[i],     // maximum overhang off end of reference
			true,         // expected return val from otherMate
			false,        // whether to look for opposite to left
			leftExt[i],   // expected left extent of dynamic prog. problem
			rightExt[i],  // expected right extent of dynamic prog. problem
			false);       // expected orientation in which opposite mate must align
	}
	}

	{
	int ohang[]    = {   0,   1,   0,   1,   2 };
	int maxgaps[]  = {   0,   0,   1,   1,   1 };
	int leftExt[]  = {   0,  -1,   0,  -1,  -2 };
	int rightExt[] = {  59,  59,  60,  60,  60 };
	for(int i = 0; i < 5; i++) {
		ostringstream oss;
		oss << "Overhang2_";
		oss << (i+1);
		testCaseOtherMate(
			oss.str(),
			PE_POLICY_FR, // policy
			200,          // maxfrag
			50,           // minfrag
			false,        // local
			true,         // dovetail OK
			true,         // containment OK
			false,        // overlap OK
			true,         // expand-to-fit
			true,         // mate 1 is anchor
			false,        // anchor aligned to Watson
			90,           // anchor's offset into ref
			200,          // ref length
			10,           // mate 1 length
			10,           // mate 2 length
			maxgaps[i],   // maximum # gaps in alignment for opposite mate
			ohang[i],     // maximum overhang off end of reference
			true,         // expected return val from otherMate
			true,         // whether to look for opposite to left
			leftExt[i],   // expected left extent of dynamic prog. problem
			rightExt[i],  // expected right extent of dynamic prog. problem
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
}

#endif /*def MAIN_PE*/
