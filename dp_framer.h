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

	DynProgFramer(bool trimToRef) : trimToRef_(trimToRef) { }

	/**
	 * Given information about a seed hit and the read that the seed came from,
	 * return parameters for the dynamic programming problem to solve.
	 */
	bool frameSeedExtensionParallelogram(
		int64_t off,      // ref offset implied by seed hit assuming no gaps
		size_t rdlen,     // length of read sequence used in DP table (so len
		                  // of +1 nucleotide sequence for colorspace reads)
		size_t reflen,    // length of reference sequence aligned to
		size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
		size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
		size_t maxhalf,   // max width in either direction
		size_t& width,    // out: calculated width stored here
		size_t& maxgap,   // out: calculated width stored here
		size_t& trimup,   // out: number of bases trimmed from upstream end
		size_t& trimdn,   // out: number of bases trimmed from downstream end
		int64_t& refl,    // out: ref pos of upper LHS of parallelogram
		int64_t& refr,    // out: ref pos of lower RHS of parallelogram
		//EList<bool>& st,  // out: legal starting columns stored here
		EList<bool>& en); // out: legal ending columns stored here

	/**
	 * Similar to frameSeedExtensionParallelogram but we're being somewhat more
	 * inclusive in order to ensure all characters aling the "width" in the last
	 * row are exhaustively scored.
	 */
	bool frameSeedExtensionRect(
		int64_t off,      // ref offset implied by seed hit assuming no gaps
		size_t rdlen,     // length of read sequence used in DP table (so len
						  // of +1 nucleotide sequence for colorspace reads)
		size_t reflen,    // length of reference sequence aligned to
		size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
		size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
		size_t maxns,     // # Ns permitted
		size_t maxhalf,   // max width in either direction
		size_t& width,    // out: calculated width stored here
		size_t& solwidth, // out: # rightmost cols where soln can end
		size_t& maxgap,   // out: calculated width stored here
		size_t& trimup,   // out: number of bases trimmed from upstream end
		size_t& trimdn,   // out: number of bases trimmed from downstream end
		int64_t& refl,    // out: ref pos of upper LHS of parallelogram
		int64_t& refr,    // out: ref pos of lower RHS of parallelogram
		EList<bool>& en); // out: legal ending columns stored here

	/**
	 * Given information about an anchor mate hit, and information deduced by
	 * PairedEndPolicy about where the opposite mate can begin and start given
	 * the fragment length range, return parameters for the dynamic programming
	 * problem to solve.
	 */
	bool frameFindMate(
		bool anchorLeft,  // true iff anchor alignment is to the left
		int64_t ll,       // leftmost Watson off for LHS of opp alignment
		int64_t lr,       // rightmost Watson off for LHS of opp alignment
		int64_t rl,       // leftmost Watson off for RHS of opp alignment
		int64_t rr,       // rightmost Watson off for RHS of opp alignment
		size_t rdlen,     // length of opposite mate
		size_t reflen,    // length of reference sequence aligned to
		size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
		size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
		size_t maxhalf,   // max width in either direction
		size_t& width,    // out: calculated width stored here
		size_t& omaxgaps, // out: calculated max # gaps stored here
		size_t& trimup,   // out: number of bases trimmed from upstream end
		size_t& trimdn,   // out: number of bases trimmed from downstream end
		int64_t& refl,    // out: ref pos of upper LHS of parallelogram
		int64_t& refr,    // out: ref pos of lower RHS of parallelogram
		//EList<bool>& st,  // out: legal starting columns stored here
		EList<bool>& en)  // out: legal ending columns stored here
		const
	{
		if(anchorLeft) {
			return frameFindMateAnchorLeft(
				ll,
				lr,
				rl,
				rr,
				rdlen,
				reflen,
				maxrdgap,
				maxrfgap,
				maxhalf,
				width,
				omaxgaps,
				trimup,
				trimdn,
				refl,
				refr,
				//st,
				en);
		} else {
			return frameFindMateAnchorRight(
				ll,
				lr,
				rl,
				rr,
				rdlen,
				reflen,
				maxrdgap,
				maxrfgap,
				maxhalf,
				width,
				omaxgaps,
				trimup,
				trimdn,
				refl,
				refr,
				//st,
				en);
		}
	}

	/**
	 * Given information about an anchor mate hit, and information deduced by
	 * PairedEndPolicy about where the opposite mate can begin and start given
	 * the fragment length range, return parameters for the dynamic programming
	 * problem to solve.
	 */
	bool frameFindMateRect(
		bool anchorLeft,  // true iff anchor alignment is to the left
		int64_t ll,       // leftmost Watson off for LHS of opp alignment
		int64_t lr,       // rightmost Watson off for LHS of opp alignment
		int64_t rl,       // leftmost Watson off for RHS of opp alignment
		int64_t rr,       // rightmost Watson off for RHS of opp alignment
		size_t rdlen,     // length of opposite mate
		size_t reflen,    // length of reference sequence aligned to
		size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
		size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
		size_t maxns,     // max # Ns permitted
		size_t maxhalf,   // max width in either direction
		size_t& width,    // out: calculated width stored here
		size_t& solwidth, // out: # rightmost cols where solution can end
		size_t& omaxgaps, // out: calculated max # gaps stored here
		size_t& trimup,   // out: number of bases trimmed from upstream end
		size_t& trimdn,   // out: number of bases trimmed from downstream end
		int64_t& refl,    // out: ref pos of upper LHS of parallelogram
		int64_t& refr,    // out: ref pos of lower RHS of parallelogram
		EList<bool>& en)  // out: legal ending columns stored here
		const
	{
		if(anchorLeft) {
			return frameFindMateAnchorLeftRect(
				ll,
				lr,
				rl,
				rr,
				rdlen,
				reflen,
				maxrdgap,
				maxrfgap,
				maxns,
				maxhalf,
				width,
				solwidth,
				omaxgaps,
				trimup,
				trimdn,
				refl,
				refr,
				en);
		} else {
			return frameFindMateAnchorRightRect(
				ll,
				lr,
				rl,
				rr,
				rdlen,
				reflen,
				maxrdgap,
				maxrfgap,
				maxns,
				maxhalf,
				width,
				solwidth,
				omaxgaps,
				trimup,
				trimdn,
				refl,
				refr,
				en);
		}
	}

	/**
	 * Given information about an anchor mate hit, and information deduced by
	 * PairedEndPolicy about where the opposite mate can begin and start given
	 * the fragment length range, return parameters for the dynamic programming
	 * problem to solve.
	 */
	bool frameFindMateAnchorLeft(
		int64_t ll,       // leftmost Watson off for LHS of opp alignment
		int64_t lr,       // rightmost Watson off for LHS of opp alignment
		int64_t rl,       // leftmost Watson off for RHS of opp alignment
		int64_t rr,       // rightmost Watson off for RHS of opp alignment
		size_t rdlen,     // length of opposite mate
		size_t reflen,    // length of reference sequence aligned to
		size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
		size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
		size_t maxhalf,   // max width in either direction
		size_t& width,    // out: calculated width stored here
		size_t& omaxgaps, // out: max # gaps
		size_t& trimup,   // out: number of bases trimmed from upstream end
		size_t& trimdn,   // out: number of bases trimmed from downstream end
		int64_t& refl,    // out: ref pos of upper LHS of parallelogram
		int64_t& refr,    // out: ref pos of lower RHS of parallelogram
		//EList<bool>& st,  // out: legal starting columns stored here
		EList<bool>& en)  // out: legal ending columns stored here
		const;

	/**
	 * Given information about an anchor mate hit, and information deduced by
	 * PairedEndPolicy about where the opposite mate can begin and start given
	 * the fragment length range, return parameters for the dynamic programming
	 * problem to solve.
	 */
	bool frameFindMateAnchorLeftRect(
		int64_t ll,       // leftmost Watson off for LHS of opp alignment
		int64_t lr,       // rightmost Watson off for LHS of opp alignment
		int64_t rl,       // leftmost Watson off for RHS of opp alignment
		int64_t rr,       // rightmost Watson off for RHS of opp alignment
		size_t rdlen,     // length of opposite mate
		size_t reflen,    // length of reference sequence aligned to
		size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
		size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
		size_t maxns,     // max # Ns permitted in alignment
		size_t maxhalf,   // max width in either direction
		size_t& width,    // out: calculated width stored here
		size_t& solwidth, // out: # rightmost cols where solution can end
		size_t& omaxgaps, // out: max # gaps
		size_t& trimup,   // out: number of bases trimmed from upstream end
		size_t& trimdn,   // out: number of bases trimmed from downstream end
		int64_t& refl,    // out: ref pos of upper LHS of parallelogram
		int64_t& refr,    // out: ref pos of lower RHS of parallelogram
		EList<bool>& en)  // out: legal ending columns stored here
		const;

	/**
	 * Given information about an anchor mate hit, and information deduced by
	 * PairedEndPolicy about where the opposite mate can begin and start given
	 * the fragment length range, return parameters for the dynamic programming
	 * problem to solve.
	 */
	bool frameFindMateAnchorRight(
		int64_t ll,       // leftmost Watson off for LHS of opp alignment
		int64_t lr,       // rightmost Watson off for LHS of opp alignment
		int64_t rl,       // leftmost Watson off for RHS of opp alignment
		int64_t rr,       // rightmost Watson off for RHS of opp alignment
		size_t rdlen,     // length of opposite mate
		size_t reflen,    // length of reference sequence aligned to
		size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
		size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
		size_t maxhalf,   // max width in either direction
		size_t& width,    // out: calculated width stored here
		size_t& omaxgaps, // out: max # gaps
		size_t& trimup,   // out: number of bases trimmed from upstream end
		size_t& trimdn,   // out: number of bases trimmed from downstream end
		int64_t& refl,    // out: ref pos of upper LHS of parallelogram
		int64_t& refr,    // out: ref pos of lower RHS of parallelogram
		EList<bool>& en)  // out: legal ending columns stored here
		const;

	/**
	 * Given information about an anchor mate hit, and information deduced by
	 * PairedEndPolicy about where the opposite mate can begin and start given
	 * the fragment length range, return parameters for the dynamic programming
	 * problem to solve.
	 */
	bool frameFindMateAnchorRightRect(
		int64_t ll,       // leftmost Watson off for LHS of opp alignment
		int64_t lr,       // rightmost Watson off for LHS of opp alignment
		int64_t rl,       // leftmost Watson off for RHS of opp alignment
		int64_t rr,       // rightmost Watson off for RHS of opp alignment
		size_t rdlen,     // length of opposite mate
		size_t reflen,    // length of reference sequence aligned to
		size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
		size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
		size_t maxns,     // max # Ns permitted in alignment
		size_t maxhalf,   // max width in either direction
		size_t& width,    // out: calculated width stored here
		size_t& solwidth, // out: # rightmost cols where solution can end
		size_t& omaxgaps, // out: max # gaps
		size_t& trimup,   // out: number of bases trimmed from upstream end
		size_t& trimdn,   // out: number of bases trimmed from downstream end
		int64_t& refl,    // out: ref pos of upper LHS of parallelogram
		int64_t& refr,    // out: ref pos of lower RHS of parallelogram
		EList<bool>& en)  // out: legal ending columns stored here
		const;

protected:

	/**
	 * Trim the given parallelogram width and reference window so that neither
	 * overhangs the beginning or end of the reference.  Return true if width
	 * is still > 0 after trimming, otherwise return false.
	 */
	void trimToRef(
		size_t   reflen,  // in: length of reference sequence aligned to
		int64_t& refl,    // in/out: ref pos of upper LHS of parallelogram
		int64_t& refr,    // in/out: ref pos of lower RHS of parallelogram
		size_t&  trimup,  // out: number of bases trimmed from upstream end
		size_t&  trimdn)  // out: number of bases trimmed from downstream end
	{
		if(refl < 0) {
			trimup = (size_t)(-refl);
			//refl = 0;
		}
		if(refr >= (int64_t)reflen) {
			trimdn = (size_t)(refr - reflen + 1);
			//refr = (int64_t)reflen-1;
		}
	}

	bool trimToRef_;
};

#endif /*ndef DP_FRAMER_H_*/
