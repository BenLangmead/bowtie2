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
#include "ref_coord.h"

/**
 * Describes a dynamic programming rectangle.
 *
 * Only knows about reference offsets, not reference sequences.
 */
struct DPRect {

	DPRect(int cat = 0) /*: st(cat), en(cat)*/ {
		refl = refr = triml = trimr = corel = corer = 0;
	}

	int64_t refl;         // leftmost ref offset involved post trimming (incl)
	int64_t refr;         // rightmost ref offset involved post trimming (incl)

	int64_t refl_pretrim; // leftmost ref offset involved pre trimming (incl)
	int64_t refr_pretrim; // rightmost ref offset involved pre trimming (incl)
	
	size_t  triml;        // positions trimmed from LHS
	size_t  trimr;        // positions trimmed from RHS
	
	// If "core" diagonals are specified, then any alignment reported has to
	// overlap one of the core diagonals.  This is to avoid the situation where
	// an alignment is reported that overlaps a better-scoring alignment that
	// falls partially outside the rectangle.  This is used in both seed
	// extensions and in mate finding.  Filtering based on the core diagonals
	// should happen in the backtrace routine.  I.e. it should simply never
	// return an alignment that doesn't overlap a core diagonal, even if there
	// is such an alignment and it's valid.
	
	size_t  corel; // offset of column where leftmost "core" diagonal starts
	size_t  corer; // offset of column where rightmost "core" diagonal starts
	// [corel, corer] is an inclusive range and offsets are with respect to the
	// original, untrimmed rectangle.
	
	size_t  maxgap; // max # gaps - width of the gap bands
	
	void write(std::ostream& os) const {
		os << refl << ',' << refr << ',' << refl_pretrim << ','
		   << refr_pretrim << ',' << triml << ',' << trimr << ','
		   << corel << ',' << corer << ',' << maxgap;
	}
	
	/**
	 * Return true iff the combined effect of triml and trimr is to trim away
	 * the entire rectangle.
	 */
	bool entirelyTrimmed() const {
		bool tr = refr < refl;
		ASSERT_ONLY(size_t width = (size_t)(refr_pretrim - refl_pretrim + 1));
		assert(tr == (width <= triml + trimr));
		return tr;
	}
	
#ifndef NDEBUG
	bool repOk() const {
		assert_geq(corer, corel);
		return true;
	}
#endif
	
	/**
	 * Set the given interval to the range of diagonals that are "covered" by
	 * this dynamic programming problem.
	 */
	void initIval(Interval& iv) {
		iv.setOff(refl_pretrim + (int64_t)corel);
		iv.setLen(corer - corel + 1);
	}
};

/**
 * Encapsulates routines for calculating parameters for the various types of
 * dynamic programming problems solved in Bowtie2.
 */
class DynProgFramer {

public:

	DynProgFramer(bool trimToRef) : trimToRef_(trimToRef) { }

	/**
	 * Similar to frameSeedExtensionParallelogram but we're being somewhat more
	 * inclusive in order to ensure all characters aling the "width" in the last
	 * row are exhaustively scored.
	 */
	bool frameSeedExtensionRect(
		int64_t off,      // ref offset implied by seed hit assuming no gaps
		size_t rdlen,     // length of read sequence used in DP table (so len
						  // of +1 nucleotide sequence for colorspace reads)
		int64_t reflen,   // length of reference sequence aligned to
		size_t maxrdgap,  // max # of read gaps permitted in opp mate alignment
		size_t maxrfgap,  // max # of ref gaps permitted in opp mate alignment
		int64_t maxns,    // # Ns permitted
		size_t maxhalf,   // max width in either direction
		DPRect& rect);    // out: DP rectangle

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
		size_t  rdlen,    // length of opposite mate
		int64_t reflen,   // length of reference sequence aligned to
		size_t  maxrdgap, // max # of read gaps permitted in opp mate alignment
		size_t  maxrfgap, // max # of ref gaps permitted in opp mate alignment
		int64_t maxns,    // max # Ns permitted
		size_t  maxhalf,  // max width in either direction
		DPRect& rect)     // out: DP rectangle
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
				rect);
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
				rect);
		}
	}

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
		size_t  rdlen,    // length of opposite mate
		int64_t reflen,   // length of reference sequence aligned to
		size_t  maxrdgap, // max # of read gaps permitted in opp mate alignment
		size_t  maxrfgap, // max # of ref gaps permitted in opp mate alignment
		int64_t maxns,    // max # Ns permitted in alignment
		size_t  maxhalf,  // max width in either direction
		DPRect& rect)     // out: DP rectangle
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
		size_t  rdlen,    // length of opposite mate
		int64_t reflen,   // length of reference sequence aligned to
		size_t  maxrdgap, // max # of read gaps permitted in opp mate alignment
		size_t  maxrfgap, // max # of ref gaps permitted in opp mate alignment
		int64_t maxns,    // max # Ns permitted in alignment
		size_t  maxhalf,  // max width in either direction
		DPRect& rect)     // out: DP rectangle
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
