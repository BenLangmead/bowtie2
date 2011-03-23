/*
 *  dp_framer.h
 *
 * Classes and routines for framing dynamic programming problems.  There are 2
 * basic types of dynamic programming problems solved in Bowtie 2:
 *
 * 1. Seed extension: we found a seed hit using Burrows-Wheeler techniques and
 *    now we would like to extend it into a full alignment by doing dynamic
 *    programming in the vicinity of the seed hit.
 *
 * 2. Mate finding: we would a full alignment for one mate in a pair and now we
 *    would like to extend it into a full alignment by doing dynamic
 *    programming in the area prescribed by the maximum and minimum fragment
 *    lengths.
 *
 * By "framing" the dynamic programming problem, we mean that all of the
 * following DP inputs are calculated:
 *
 * 1. The width of the parallelogram/rectangle to explore.
 * 2. The 0-based offset of the reference position associated with the leftmost
 *    diagnomal/column in the parallelogram/rectangle to explore
 * 3. An EList<bool> of length=width encoding which columns the alignment may
 *    start in
 * 4. An EList<bool> of length=width encoding which columns the alignment may
 *    end in
 */

#ifndef DP_FRAMER_H_
#define DP_FRAMER_H_

#include <stdint.h>
#include "ds.h"

/**
 * Encapsulates routines for calculating parameters for the various types of
 * dynamic programming problems solved in Bowtie2.
 */
class DynProgFramer {

public:

	DynProgFramer() { }

	/**
	 * Given information about a seed hit and the read that the seed came from,
	 * return parameters for the dynamic programming problem to solve.
	 */
	void frameSeedExtension(
		int64_t off,      // ref offset implied by seed hit assuming no gaps
		size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
		size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
		size_t& width,    // out: calculated width stored here
		int64_t& refi,    // out: ref pos associated with first/leftmost column
		EList<bool>& st,  // out: legal starting columns stored here
		EList<bool>& en); // out: legal ending columns stored here

	/**
	 * Given information about an anchor mate hit, and information deduced by
	 * PairedEndPolicy about where the opposite mate can begin and start given
	 * the fragment length range, return parameters for the dynamic programming
	 * problem to solve.
	 */
	void frameFindMateAnchorLeft(
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
		EList<bool>& en); // out: legal ending columns stored here

	/**
	 * Given information about an anchor mate hit, and information deduced by
	 * PairedEndPolicy about where the opposite mate can begin and start given
	 * the fragment length range, return parameters for the dynamic programming
	 * problem to solve.
	 */
	void frameFindMateAnchorRight(
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
		EList<bool>& en); // out: legal ending columns stored here

};

#endif /*ndef DP_FRAMER_H_*/
