//
//  dp_framer.cpp
//

#include <algorithm>
#include "dp_framer.h"

using namespace std;

/**
 * Given information about a seed hit and the read that the seed came from,
 * return parameters for the dynamic programming problem to solve.
 */
bool DynProgFramer::frameSeedExtension(
	int64_t off,      // ref offset implied by seed hit assuming no gaps
	size_t rdlen,     // length of read sequence used in DP table (so len
	                  // of +1 nucleotide sequence for colorspace reads)
	size_t reflen,    // length of reference sequence aligned to
	size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
	size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
	size_t& width,    // out: calculated width stored here
	size_t& trimup,   // out: number of bases trimmed from upstream end
	size_t& trimdn,   // out: number of bases trimmed from downstream end
	int64_t& refl,    // out: ref pos of upper LHS of parallelogram
	int64_t& refr,    // out: ref pos of lower RHS of parallelogram
	EList<bool>& st,  // out: legal starting columns stored here
	EList<bool>& en)  // out: legal ending columns stored here
{
	assert_gt(rdlen, 0);
	assert_gt(reflen, 0);
	size_t maxgap = max(maxrdgap, maxrfgap);
	width = 1 + (2 * maxgap);
	refl = off - maxgap;
	refr = off + (rdlen - 1) + maxgap;
	trimup = trimdn = 0;
	// Check if we have to trim to fit the extents of the reference
	if(trimToRef_) {
		trimToRef(reflen, refl, refr, trimup, trimdn);
		// Occassionally, trimToRef trims the whole problem away (because the
		// reference is actually too short to support a valid alignment).
		if(trimup >= width || trimdn >= width) {
			return false;
		}
		// Trimming affects which cells are set to false, not the width.  TODO:
		// it could also affect the width, but we're not checking for that.
	}
	assert_gt(width, 0);
	st.resize(width);
	en.resize(width);
	// Most common & simplest case is maxrdgap == maxrfgap; handle that first
	for(size_t i = 0; i < width; i++) {
		st[i] = en[i] = true;
	}
	if(maxrdgap < maxrfgap && width > 1) {
		// More read gaps than ref gaps; some cells at RHS of 'st' and LHS of
		// 'en' aren't valid
		size_t diff = maxrfgap - maxrdgap;
		for(size_t i = 0; i < diff; i++) {
			assert_geq(width, i+1);
			st[i] = en[width - i - 1] = false;
		}
	} else if(maxrfgap < maxrdgap && width > 1) {
		// More ref gaps than read gaps; some cells at LHS of 'st' and RHS of
		// 'en' aren't valid
		size_t diff = maxrdgap - maxrfgap;
		for(size_t i = 0; i < diff; i++) {
			assert_geq(width, i+1);
			st[width - i - 1] = en[i] = false;
		}
	}
	for(size_t i = 0; i < trimup; i++) {
		st[i] = false;
	}
	for(size_t i = 0; i < trimdn; i++) {
		en[width - i - 1] = false;
	}
	return true;
}

/**
 * Given information about an anchor mate hit, and information deduced by
 * PairedEndPolicy about where the opposite mate can begin and start given
 * the fragment length range, return parameters for the dynamic programming
 * problem to solve.
 */
bool DynProgFramer::frameFindMateAnchorLeft(
	int64_t ll,       // leftmost Watson off for LHS of opp alignment
	int64_t lr,       // rightmost Watson off for LHS of opp alignment
	int64_t rl,       // leftmost Watson off for RHS of opp alignment
	int64_t rr,       // rightmost Watson off for RHS of opp alignment
	size_t rdlen,     // length of opposite mate
	size_t reflen,    // length of reference sequence aligned to
	size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
	size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
	size_t& width,    // out: calculated width stored here
	size_t& trimup,   // out: number of bases trimmed from upstream end
	size_t& trimdn,   // out: number of bases trimmed from downstream end
	int64_t& refl,    // out: ref pos of upper LHS of parallelogram
	int64_t& refr,    // out: ref pos of lower RHS of parallelogram
	EList<bool>& st,  // out: legal starting columns stored here
	EList<bool>& en)  // out: legal ending columns stored here
{
	assert_geq(lr, ll);
	assert_geq(rr, rl);
	assert_geq(rr, lr);
	assert_geq(rl, ll);
	assert_gt(rdlen, 0);
	assert_gt(reflen, 0);
	trimup = trimdn = 0;
	// Amount of padding we have to add to account for the fact that alignments
	// ending between en_left/en_right might start in various columns in the
	// first row
	int64_t pad_left = maxrdgap;
	int64_t pad_right = maxrfgap;
	int64_t en_left  = rl;
	int64_t en_right = rr;
	int64_t st_left  = en_left - (rdlen-1);
	int64_t st_right = en_right - (rdlen-1);
	int64_t ltrimr = 0, ltriml = 0;
	int64_t rtrimr = 0, rtriml = 0;
	int64_t en_right_pad = en_right + pad_right;
	int64_t en_left_pad  = en_left  - pad_left;
	int64_t st_right_pad = st_right + pad_right;
	int64_t st_left_pad  = st_left  - pad_left;
	// Trim both mates to reference, if required
	if(trimToRef_) {
		if(en_right_pad >= (int64_t)reflen) rtrimr = en_right_pad - reflen + 1;
		if(en_left_pad < 0)                 rtriml = -en_left_pad;
		if(st_right_pad >= (int64_t)reflen) ltrimr = st_right_pad - reflen + 1;
		if(st_left_pad < 0)                 ltriml = -st_left_pad;
	}
	// Trim left mate to ll/lr
	if(st_right_pad > lr && (st_right_pad - lr) > ltrimr) {
		ltrimr = st_right_pad - lr;
	}
	if(st_left_pad < ll && (ll - st_left_pad) > ltriml) {
		ltriml = ll - st_left_pad;
	}
	rtrimr = max<size_t>(rtrimr, pad_right);
	rtriml = max<size_t>(rtriml, pad_left);
	// Are we trimming so much of the left (top) interval that we can also trim
	// more from the right (bottom) interval?
	{
		int64_t trimr = 0, triml = 0;
		int64_t mxrdgap = (int64_t)maxrdgap;
		if(ltrimr > rtrimr + mxrdgap) {
			trimr = ltrimr - mxrdgap;
		} else if(rtrimr > ltrimr + mxrdgap) {
			trimr = rtrimr - mxrdgap;
		} else {
			trimr = min(rtrimr, ltrimr);
		}
		int64_t mxrfgap = (int64_t)maxrfgap;
		if(ltriml > rtriml + mxrfgap) {
			triml = ltriml - mxrfgap;
		} else if(rtriml > ltriml + mxrfgap) {
			triml = rtriml - mxrfgap;
		} else {
			triml = min(rtriml, ltriml);
		}
		if(trimr > ltrimr) ltrimr = 0;
		else ltrimr -= trimr;
		if(triml > ltriml) ltriml = 0;
		else ltriml -= triml;
		if(trimr > rtrimr) rtrimr = 0;
		else rtrimr -= trimr;
		if(triml > rtriml) rtriml = 0;
		else rtriml -= triml;
		if(trimr > 0) {
			trimdn = trimr;
			if(pad_right < trimr) {
				trimr -= pad_right;
				pad_right = 0;
				en_right -= trimr;
				st_right -= trimr;
			} else pad_right -= trimr;
		}
		if(triml > 0) {
			trimup = triml;
			if(pad_left < triml) {
				triml -= pad_left;
				pad_left = 0;
				en_left += triml;
				st_left += triml;
			} else pad_left -= triml;
		}
	}
	if(en_left > en_right) {
		// Trimmed everything
		return false;
	}
	assert(!trimToRef_ ||
		((st_left+ltriml) >= 0 && en_left >= 0));
	assert(!trimToRef_ ||
		(st_right < (int64_t)reflen &&
		 (en_right-rtrimr) < (int64_t)reflen));
	// Calculate width taking gaps into account
	width = (size_t)(en_right - en_left + 1 + pad_left + pad_right);
	assert_gt(width, 0);
	st.resize(width);
	en.resize(width);
	// Apply padding
	st_left -= pad_left;
	en_left -= pad_left;
	st_right += pad_right;
	en_right += pad_right;
	assert_eq(en_right - en_left, st_right - st_left);
	for(size_t i = 0; i < width; i++) {
		st[i] = en[i] = true;
	}
	for(int64_t i = 0; i < rtriml; i++) {
		en[i] = false;
	}
	for(int64_t i = 0; i < rtrimr; i++) {
		en[width - i - 1] = false;
	}
	for(int64_t i = 0; i < ltrimr; i++) {
		st[width - i - 1] = false;
	}
	for(int64_t i = 0; i < ltriml; i++) {
		st[i] = false;
	}
	refl = st_left;
	refr = en_right;
	return true;
}

/**
 * Given information about an anchor mate hit, and information deduced by
 * PairedEndPolicy about where the opposite mate can begin and start given
 * the fragment length range, return parameters for the dynamic programming
 * problem to solve.
 */
bool DynProgFramer::frameFindMateAnchorRight(
	int64_t ll,       // leftmost Watson off for LHS of opp alignment
	int64_t lr,       // rightmost Watson off for LHS of opp alignment
	int64_t rl,       // leftmost Watson off for RHS of opp alignment
	int64_t rr,       // rightmost Watson off for RHS of opp alignment
	size_t rdlen,     // length of opposite mate
	size_t reflen,    // length of reference sequence aligned to
	size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
	size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
	size_t& width,    // out: calculated width stored here
	size_t& trimup,   // out: number of bases trimmed from upstream end
	size_t& trimdn,   // out: number of bases trimmed from downstream end
	int64_t& refl,    // out: ref pos of upper LHS of parallelogram
	int64_t& refr,    // out: ref pos of lower RHS of parallelogram
	EList<bool>& st,  // out: legal starting columns stored here
	EList<bool>& en)  // out: legal ending columns stored here
{
	assert_geq(lr, ll);
	assert_geq(rr, rl);
	assert_geq(rr, lr);
	assert_geq(rl, ll);
	assert_gt(rdlen, 0);
	assert_gt(reflen, 0);
	trimup = trimdn = 0;
	// Amount of padding we have to add to account for the fact that alignments
	// ending between en_left/en_right might start in various columns in the
	// first row
	int64_t pad_left = maxrfgap;
	int64_t pad_right = maxrdgap;
	int64_t st_left = ll;
	int64_t st_right = lr;
	int64_t en_left = st_left + (rdlen-1);
	int64_t en_right = st_right + (rdlen-1);
	int64_t ltrimr = 0, ltriml = 0;
	int64_t rtrimr = 0, rtriml = 0;
	int64_t en_right_pad = en_right + pad_right;
	int64_t en_left_pad  = en_left  - pad_left;
	int64_t st_right_pad = st_right + pad_right;
	int64_t st_left_pad  = st_left  - pad_left;
	// Trim both mates to reference, if required
	if(trimToRef_) {
		if(en_right_pad >= (int64_t)reflen) rtrimr = en_right_pad - reflen + 1;
		if(en_left_pad < 0)                 rtriml = -en_left_pad;
		if(st_right_pad >= (int64_t)reflen) ltrimr = st_right_pad - reflen + 1;
		if(st_left_pad < 0)                 ltriml = -st_left_pad;
	}
	// Trim left mate to ll/lr
	if(en_right_pad > rr && (en_right_pad - rr) > rtrimr) {
		rtrimr = en_right_pad - rr;
	}
	if(en_left_pad < rl && (rl - en_left_pad) > rtriml) {
		rtriml = rl - en_left_pad;
	}
	ltrimr = max<size_t>(ltrimr, pad_right);
	ltriml = max<size_t>(ltriml, pad_left);
	// Are we trimming so much of the left (top) interval that we can also trim
	// more from the right (bottom) interval, or vice versa?
	{
		int64_t trimr = 0, triml = 0;
		int64_t mxrfgap = (int64_t)maxrfgap;
		if(ltrimr > rtrimr + mxrfgap) {
			trimr = ltrimr - mxrfgap;
		} else if(rtrimr > ltrimr + mxrfgap) {
			trimr = rtrimr - mxrfgap;
		} else {
			trimr = min(rtrimr, ltrimr);
		}
		int64_t mxrdgap = (int64_t)maxrdgap;
		if(ltriml > rtriml + mxrdgap) {
			triml = ltriml - mxrdgap;
		} else if(rtriml > ltriml + mxrdgap) {
			triml = rtriml - mxrdgap;
		} else {
			triml = min(rtriml, ltriml);
		}
		if(trimr > ltrimr) ltrimr = 0;
		else ltrimr -= trimr;
		if(triml > ltriml) ltriml = 0;
		else ltriml -= triml;
		if(trimr > rtrimr) rtrimr = 0;
		else rtrimr -= trimr;
		if(triml > rtriml) rtriml = 0;
		else rtriml -= triml;
		if(trimr > 0) {
			trimdn = trimr;
			if(pad_right < trimr) {
				trimr -= pad_right;
				pad_right = 0;
				en_right -= trimr;
				st_right -= trimr;
			} else pad_right -= trimr;
		}
		if(triml > 0) {
			trimup = triml;
			if(pad_left < triml) {
				triml -= pad_left;
				pad_left = 0;
				en_left += triml;
				st_left += triml;
			} else pad_left -= triml;
		}
	}
	if(en_left > en_right) {
		// Trimmed everything
		return false;
	}
	assert(!trimToRef_ ||
		((st_left+ltriml) >= 0 && en_left >= 0));
	assert(!trimToRef_ ||
		(st_right < (int64_t)reflen &&
		 (en_right-rtrimr) < (int64_t)reflen));
	// Calculate width taking gaps into account
	width = (size_t)(en_right - en_left + 1 + pad_left + pad_right);
	assert_gt(width, 0);
	st.resize(width);
	en.resize(width);
	// Apply padding
	st_left -= pad_left;
	en_left -= pad_left;
	st_right += pad_right;
	en_right += pad_right;
	assert_eq(en_right - en_left, st_right - st_left);
	for(size_t i = 0; i < width; i++) {
		st[i] = en[i] = true;
	}
	for(int64_t i = 0; i < ltriml; i++) {
		st[i] = false;
	}
	for(int64_t i = 0; i < ltrimr; i++) {
		st[width - i - 1] = false;
	}
	for(int64_t i = 0; i < rtrimr; i++) {
		en[width - i - 1] = false;
	}
	for(int64_t i = 0; i < rtriml; i++) {
		en[i] = false;
	}
	refl = st_left;
	refr = en_right;
	return true;
}

#ifdef MAIN_DP_FRAMER

#include <iostream>

static void testCaseSeedExtension(
	const char *testName,
	bool trimToRef,
	int64_t off,
	size_t rdlen,
	size_t reflen,
	size_t maxrdgap,
	size_t maxrfgap,
	size_t ex_width,
	size_t ex_trimup,
	size_t ex_trimdn,
	int64_t ex_refl,
	int64_t ex_refr,
	const char *ex_st,    // string of '0'/'1' chars
	const char *ex_en)    // string of '0'/'1' chars
{
	cerr << testName << "...";
	DynProgFramer fr(trimToRef);
	size_t width;
	int64_t refl, refr;
	EList<bool> st, en;
	size_t trimup, trimdn;
	fr.frameSeedExtension(
		off,      // ref offset implied by seed hit assuming no gaps
		rdlen,    // length of read sequence used in DP table (so len
		          // of +1 nucleotide sequence for colorspace reads)
		reflen,   // length of reference sequence aligned to
		maxrdgap, // max # of read gaps permitted in opp mate alignment
		maxrfgap, // max # of ref gaps permitted in opp mate alignment
		width,    // out: calculated width stored here
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
	size_t ex_trimup,
	size_t ex_trimdn,
	int64_t ex_refl,
	int64_t ex_refr,
	const char *ex_st,    // string of '0'/'1' chars
	const char *ex_en)    // string of '0'/'1' chars
{
	cerr << testName << "...";
	DynProgFramer fr(trimToRef);
	size_t width;
	int64_t refl, refr;
	EList<bool> st, en;
	size_t trimup, trimdn;
	fr.frameFindMateAnchorLeft(
		ll,       // leftmost Watson off for LHS of opp alignment
		lr,       // rightmost Watson off for LHS of opp alignment
		rl,       // leftmost Watson off for RHS of opp alignment
		rr,       // rightmost Watson off for RHS of opp alignment
		rdlen,    // length of opposite mate
		reflen,   // length of reference sequence aligned to
		maxrdgap, // max # of read gaps permitted in opp mate alignment
		maxrfgap, // max # of ref gaps permitted in opp mate alignment
		width,    // out: calculated width stored here
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
	size_t ex_trimup,
	size_t ex_trimdn,
	int64_t ex_refl,
	int64_t ex_refr,
	const char *ex_st,    // string of '0'/'1' chars
	const char *ex_en)    // string of '0'/'1' chars
{
	cerr << testName << "...";
	DynProgFramer fr(trimToRef);
	size_t width;
	int64_t refl, refr;
	EList<bool> st, en;
	size_t trimup, trimdn;
	fr.frameFindMateAnchorRight(
		ll,       // leftmost Watson off for LHS of opp alignment
		lr,       // rightmost Watson off for LHS of opp alignment
		rl,       // leftmost Watson off for RHS of opp alignment
		rr,       // rightmost Watson off for RHS of opp alignment
		rdlen,    // length of opposite mate
		reflen,   // length of reference sequence aligned to
		maxrdgap, // max # of read gaps permitted in opp mate alignment
		maxrfgap, // max # of ref gaps permitted in opp mate alignment
		width,    // out: calculated width stored here
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
	//           v off
	//           *
	//          ooo
	//           ooo
	//            ooo
	//             ooo
	//              ooo
	// 012345678901234567890
	// 0         1         2
	// Note: length of read isn't known here
	testCaseSeedExtension(
		"SeedExtension1", // name
		false,            // trim to reference
		10,               // offset
		5,                // read length
		30,               // reference length
		1,                // max read gap
		1,                // max ref gap
		3,                // width of parallelogram
		0,                // # bases trimmed from upstream end
		0,                // # bases trimmed from downstream end
		9,                // ref offset of upstream column
		15,               // ref offset of downstream column
		"111",            // expected starting bools
		"111");           // expected ending bools
	testCaseSeedExtension(
		"SeedExtension2", // name
		false,            // trim to reference
		10,               // offset
		5,                // read length
		30,               // reference length
		0,                // max read gap
		1,                // max ref gap
		3,                // width of parallelogram
		0,                // # bases trimmed from upstream end
		0,                // # bases trimmed from downstream end
		9,                // ref offset of upstream column
		15,               // ref offset of downstream column
		"011",            // expected starting bools
		"110");           // expected ending bools
	
	testCaseSeedExtension(
		"SeedExtension3", // name
		false,            // trim to reference
		10,               // offset
		5,                // read length
		30,               // reference length
		1,                // max read gap
		0,                // max ref gap
		3,                // width of parallelogram
		0,                // expected # bases trimmed from upstream end
		0,                // expected # bases trimmed from downstream end
		9,                // ref offset of upstream column
		15,               // ref offset of downstream column
		"110",            // expected starting bools
		"011");           // expected ending bools

	testCaseSeedExtension(
		"SeedExtension4", // name
		false,            // trim to reference
		0,                // offset
		5,                // read length
		3,                // reference length
		2,                // max read gap
		2,                // max ref gap
		5,                // width of parallelogram
		0,                // expected # bases trimmed from upstream end
		0,                // expected # bases trimmed from downstream end
		-2,               // ref offset of upstream column
		6,                // ref offset of downstream column
		"11111",          // expected starting bools
		"11111");         // expected ending bools
	
	testCaseSeedExtension(
		"SeedExtension5", // name
		true,             // trim to reference
		0,                // offset
		5,                // read length
		3,                // reference length
		2,                // max read gap
		2,                // max ref gap
		5,                // width of parallelogram
		2,                // expected # bases trimmed from upstream end
		4,                // expected # bases trimmed from downstream end
		-2,               // ref offset of upstream column
		6,                // ref offset of downstream column
		"00111",          // expected starting bools
		"10000");         // expected ending bools
	
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
}

#endif /*def MAIN_DP_FRAMER*/
