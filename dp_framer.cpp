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
void DynProgFramer::frameSeedExtension(
	int64_t off,      // ref offset implied by seed hit assuming no gaps
	size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
	size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
	size_t& width,    // out: calculated width stored here
	int64_t& refi,    // out: ref pos associated with first/leftmost column
	EList<bool>& st,  // out: legal starting columns stored here
	EList<bool>& en)  // out: legal ending columns stored here
{
	size_t maxgap = max(maxrdgap, maxrfgap);
	width = 1 + (2 * maxgap);
	refi = off - maxgap;
	st.resize(width);
	en.resize(width);
	// Most common & simplest case is maxrdgap == maxrfgap; handle that first
	for(size_t i = 0; i < width; i++) {
		st[i] = en[i] = true;
	}
	if(maxrdgap < maxrfgap) {
		// More read gaps than ref gaps; some cells at RHS of 'st' and LHS of
		// 'en' aren't valid
		size_t diff = maxrfgap - maxrdgap;
		for(size_t i = 0; i < diff; i++) {
			st[i] = en[width - i - 1] = false;
		}
	} else if(maxrfgap < maxrdgap) {
		// More ref gaps than read gaps; some cells at LHS of 'st' and RHS of
		// 'en' aren't valid
		size_t diff = maxrdgap - maxrfgap;
		for(size_t i = 0; i < diff; i++) {
			st[width - i - 1] = en[i] = false;
		}
	}
	// All set for dynamic programming
}

/**
 * Given information about an anchor mate hit, and information deduced by
 * PairedEndPolicy about where the opposite mate can begin and start given
 * the fragment length range, return parameters for the dynamic programming
 * problem to solve.
 */
void DynProgFramer::frameFindMateAnchorLeft(
	int64_t ll,       // leftmost Watson off for LHS of opp alignment
	int64_t lr,       // rightmost Watson off for LHS of opp alignment
	int64_t rl,       // leftmost Watson off for RHS of opp alignment
	int64_t rr,       // rightmost Watson off for RHS of opp alignment
	size_t rdlen,     // length of opposite mate
	size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
	size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
	size_t& width,    // out: calculated width stored here
	int64_t& refi,    // out: ref pos associated with first/leftmost column
	EList<bool>& st,  // out: legal starting columns stored here
	EList<bool>& en)  // out: legal ending columns stored here
{
	assert_geq(lr, ll);
	assert_geq(rr, rl);
	assert_gt(rdlen, 0);
	// en_left / en_right are the ref offsets of the terminal cells where
	// alignments that are acceptable according to the fragment length bounds
	// can end
	int64_t en_left = rl;
	int64_t en_right = rr;
	
	size_t left_pad = maxrdgap;
	size_t right_pad = maxrfgap;
	width = (size_t)(en_right - en_left + 1) + left_pad + right_pad;
	int64_t st_left = rl - (rdlen-1) - maxrdgap;
	int64_t st_right = st_left + width - 1;
	refi = st_left;
	st.resize(width);
	en.resize(width);
	// Most entries will be true
	for(size_t i = 0; i < width; i++) {
		st[i] = en[i] = true;
	}
	// Set entries in 'en' that aren't in the acceptable range to false
	for(size_t i = 0; i < left_pad; i++) {
		en[i] = false;
	}
	for(size_t i = 0; i < right_pad; i++) {
		en[width - i - 1] = false;
	}
	// Set entries in 'st' that fall outside the ll/lr bounds calculated by
	// PairedEndPolicy to false
	if(st_left < ll) {
		size_t diff = (size_t)(ll - st_left);
		for(size_t i = 0; i < diff; i++) {
			st[i] = false;
		}
	}
	if(st_right > lr) {
		size_t diff = (size_t)(st_right - lr);
		for(size_t i = 0; i < diff; i++) {
			st[width - i - 1] = false;
		}
	}
}

/**
 * Given information about an anchor mate hit, and information deduced by
 * PairedEndPolicy about where the opposite mate can begin and start given
 * the fragment length range, return parameters for the dynamic programming
 * problem to solve.
 */
void DynProgFramer::frameFindMateAnchorRight(
	int64_t ll,       // leftmost Watson off for LHS of opp alignment
	int64_t lr,       // rightmost Watson off for LHS of opp alignment
	int64_t rl,       // leftmost Watson off for RHS of opp alignment
	int64_t rr,       // rightmost Watson off for RHS of opp alignment
	size_t rdlen,     // length of opposite mate
	size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
	size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
	size_t& width,    // out: calculated width stored here
	int64_t& refi,    // out: ref pos associated with first/leftmost column
	EList<bool>& st,  // out: legal starting columns stored here
	EList<bool>& en)  // out: legal ending columns stored here
{
	assert_geq(lr, ll);
	assert_geq(rr, rl);
	assert_gt(rdlen, 0);
	// en_left / en_right are the ref offsets of the terminal cells where
	// alignments that are acceptable according to the fragment length bounds
	// can end
	int64_t st_left = ll;
	int64_t st_right = lr;
	
	size_t left_pad = maxrfgap;
	size_t right_pad = maxrdgap;
	width = (size_t)(st_right - st_left + 1) + left_pad + right_pad;
	int64_t en_left = rl + (rdlen-1) - maxrfgap;
	int64_t en_right = en_left + width - 1;
	refi = st_left - left_pad;
	st.resize(width);
	en.resize(width);
	// Most entries will be true
	for(size_t i = 0; i < width; i++) {
		st[i] = en[i] = true;
	}
	// Set entries in 'st' that aren't in the acceptable range to false
	for(size_t i = 0; i < left_pad; i++) {
		st[i] = false;
	}
	for(size_t i = 0; i < right_pad; i++) {
		st[width - i - 1] = false;
	}
	// Set entries in 'en' that fall outside the rl/rr bounds calculated by
	// PairedEndPolicy to false
	if(en_left < rl) {
		size_t diff = (size_t)(rl - en_left);
		for(size_t i = 0; i < diff; i++) {
			en[i] = false;
		}
	}
	if(en_right > rr) {
		size_t diff = (size_t)(en_right - rr);
		for(size_t i = 0; i < diff; i++) {
			en[width - i - 1] = false;
		}
	}
}

#ifdef MAIN_DP_FRAMER

static void testCaseSeedExtension(
	int64_t     off,
	size_t      maxrdgap,
	size_t      maxrfgap,
	size_t      ex_width,
	int64_t     ex_refi,
	const char* ex_st,    // string of '0'/'1' chars
	const char* ex_en)    // string of '0'/'1' chars
{
	DynProgFramer fr;
	size_t width;
	int64_t refi;
	EList<bool> st, en;
	fr.frameSeedExtension(
		off,      // ref offset implied by seed hit assuming no gaps
		maxrdgap, // max # of read gaps permitted in opp mate alignment
		maxrfgap, // max # of ref gaps permitted in opp mate alignment
		width,    // out: calculated width stored here
		refi,     // out: ref pos associated with first/leftmost column
		st,       // out: legal starting columns stored here
		en);      // out: legal ending columns stored here
	assert_eq(ex_width, width);
	assert_eq(ex_refi, refi);
	for(size_t i = 0; i < width; i++) {
		assert_eq((ex_st[i] == '1'), st[i]);
		assert_eq((ex_en[i] == '1'), en[i]);
	}
}

int main(void) {
	testCaseSeedExtension(10, 1, 1, 3, 9, "111", "111");
	testCaseSeedExtension(10, 0, 1, 3, 9, "011", "110");
	testCaseSeedExtension(10, 1, 0, 3, 9, "110", "011");
}

#endif /*def MAIN_DP_FRAMER*/
