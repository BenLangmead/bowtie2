/*
 *  aligner_sw_col.h
 */

#ifndef ALIGNER_SW_COL_H_
#define ALIGNER_SW_COL_H_

#include <stdint.h>
#include "ds.h"
#include "threading.h"
#include "aligner_seed.h"
#include "reference.h"
#include "random_source.h"
#include "mem_ids.h"
#include "aligner_result.h"

/**
 * A bitmask encoding which backtracking paths out of a particular cell
 * correspond to optimal subpaths.  This struct is tailored to the
 * colorspace case where each transition down in the dynamic
 * programming matrix is also associated with a nucleotide assignment
 * to the source (upper) row.
 */
struct SwColorCellMask {

	/**
	 * Set all flags to 0, indicating there is no way to backtrack from
	 * this cell to an optimal answer.
	 */
	void clear() {
		*((uint16_t*)this) = 0;
	}

	/**
	 * Return true iff there are no backward paths recorded in this
	 * mask.
	 */
	inline bool empty() const {
		return *((uint16_t*)this) == 0;
	}

	/**
	 * Return true iff it's possible to extend a gap in the reference
	 * in the cell below this one.
	 */
	inline bool refExtendPossible() const {
		return rfop || rfex;
	}

	/**
	 * Return true iff it's possible to open a gap in the reference
	 * in the cell below this one (false implies that only extension
	 * is possible).
	 */
	inline bool refOpenPossible() const {
		return diag || rfop || rfex;
	}

	/**
	 * Return true iff it's possible to extend a gap in the read
	 * in the cell to the right of this one.
	 */
	inline bool readExtendPossible() const {
		return rdop || rdex;
	}

	/**
	 * Return true iff it's possible to open a gap in the read in the
	 * cell to the right of this one (false implies that only extension
	 * is possible).
	 */
	inline bool readOpenPossible() const {
		return diag || rdop || rdex;
	}
	
	/**
	 * Return true iff there is >0 possible way to backtrack from this
	 * cell.
	 */
	inline int numPossible() const {
		int num = 0;
		num += mask2popcnt[diag];
		num += mask2popcnt[rfop];
		num += mask2popcnt[rfex];
		num += rdop;
		num += rdex;
		return num;
	}

	/**
	 * Select a path for backtracking from this cell.  If there is a
	 * tie among eligible paths, break it randomly.  Return value is
	 * a pair where first = a flag indicating the backtrack type (see
	 * enum defining SW_BT_* above), and second = a selection for
	 * the read character for the next row up.  second should be
	 * ignored if the backtrack type is a gap in the read.
	 */
	std::pair<int, int> randBacktrack(RandomSource& rand, bool& branch);

	uint16_t diag     : 4;
	uint16_t rfop     : 4;
	uint16_t rfex     : 4;
	uint16_t rdop     : 1;
	uint16_t rdex     : 1;
	uint16_t reserved : 2;
};

/**
 * Encapsulates a backtrace stack frame.  Includes enough information that we
 * can "pop" back up to this frame and choose to make a different backtracking
 * decision.  The information included is:
 *
 * 1. The mask at the decision point.  When we first move through the mask and
 *    when we backtrack to it, we're careful to mask out the bit corresponding
 *    to the path we're taking.  When we move through it after removing the
 *    last bit from the mask, we're careful to pop it from the stack.
 * 2. The sizes of the edit lists.  When we backtrack, we resize the lists back
 *    down to these sizes to get rid of any edits introduced since the branch
 *    point.
 */
struct DpColFrame {

	/**
	 * Initialize a new DpNucFrame stack frame.
	 */
	void init(
		size_t   nedsz_,
		size_t   aedsz_,
		size_t   cedsz_,
		size_t   celsz_,
		size_t   row_,
		size_t   col_,
		int      curC_,
		int      gaps_,
		AlnScore score_)
	{
		nedsz = nedsz_;
		aedsz = aedsz_;
		cedsz = cedsz_;
		celsz = celsz_;
		row   = row_;
		col   = col_;
		curC  = curC_;
		gaps  = gaps_;
		score = score_;
	}

	size_t   nedsz; // size of the nucleotide edit list at branch (before
	                // adding the branch edit)
	size_t   aedsz; // size of ambiguous nucleotide edit list at branch
	size_t   cedsz; // size of color edit list at branch
	size_t   celsz; // size of cell-traversed list at branch
	size_t   row;   // row of cell where branch occurred
	size_t   col;   // column of cell where branch occurred
	int      curC;  // character cell we're in
	int      gaps;  // gaps before branch occurred
	AlnScore score; // score where branch occurred
};

/**
 * Encapsulates all information needed to encode the optimal subproblem
 * at a cell in a colorspace SW matrix.
 */
struct SwColorCell {

	/**
	 * Clear this cell so that it's ready for updates.
	 */
	void clear() {
		// Initially, best scores are all invalid
		best[0] = best[1] = best[2] = best[3] = AlnScore::INVALID();
		// Initially, there's no way to backtrack from this cell
		mask[0].clear();
		mask[1].clear();
		mask[2].clear();
		mask[3].clear();
		empty = true;
		reportedThru_ = false;
		reportedFrom_[0] = reportedFrom_[1] =
		reportedFrom_[2] = reportedFrom_[3] = false;
		ASSERT_ONLY(finalized = false);
	}
	
	/**
	 * Return true if any bests are valid.
	 */
	bool valid() const {
		for(int i = 0; i < 4; i++) {
			if(VALID_AL_SCORE(best[i])) {
				assert(!mask[i].empty());
				return true;
			} else {
				assert(mask[i].empty());
			}
		}
		return false;
	}

	/**
	 * We finished updating the cell; set empty and finalized
	 * appropriately.
	 */
	inline bool finalize(TAlScore minsc) {
		ASSERT_ONLY(finalized = true);
		assert(VALID_SCORE(minsc));
		for(int i = 0; i < 4; i++) {
			if(!mask[i].empty()) {
				assert(VALID_AL_SCORE(best[i]));
				assert_geq(best[i].score(), minsc);
				empty = false;
#ifdef NDEBUG
				break;
#endif
			}
		}
		return !empty;
	}

	/**
	 * Determine whether this cell has a solution with the given score.  If so,
	 * return true.  Otherwise, return false.
	 */
	inline bool bestSolutionGeq(
		const TAlScore& min,
		AlnScore& bst)
	{
		if(reportedThru_ || empty) {
			return false;
		}
		AlnScore bestSoFar;
		bestSoFar.invalidate();
		for(int i = 0; i < 4; i++) {
			if(!reportedFrom_[i] &&
			   best[i].score() >= min &&
			   best[i] > bestSoFar)
			{
				bestSoFar = best[i];
			}
		}
		if(bestSoFar.valid()) {
			bst = bestSoFar;
			return true;
		}
		return false;
	}

	/**
	 * Determine whether this cell has a solution with the given score.  If so,
	 * return true.  Otherwise, return false.
	 */
	inline bool hasSolutionEq(const TAlScore& eq) {
		if(reportedThru_ || empty) {
			return false;
		}
		for(int i = 0; i < 4; i++) {
			if(!reportedFrom_[i] && best[i].score() == eq) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Determine whether this cell has a solution with a score greater than or
	 * equal to the given score.  If so, return true.  Otherwise, return false.
	 */
	inline bool hasSolution(const TAlScore& eq) {
		if(reportedThru_ || empty) {
			return false;
		}
		for(int i = 0; i < 4; i++) {
			if(!reportedFrom_[i] && best[i].score() >= eq) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Determine whether this cell has a solution with the given score.  If so,
	 * return true.  Otherwise, return false.
	 */
	inline bool nextSolutionEq(
		const TAlScore& eq,
		int& nuc)
	{
		if(reportedThru_ || empty) {
			return false;
		}
		for(int i = 0; i < 4; i++) {
			if(!reportedFrom_[i] && best[i].score() == eq) {
				reportedFrom_[i] = true;
				nuc = i;     // set decoded nucleotide for last alignment pos
				return true; // found a cell to potentially backtrack from
			}
		}
		return false;
	}
	
	/**
	 * Mark this cell as "reported through," meaning that an already-reported
	 * alignment goes through the cell.  Future alignments that go through this
	 * cell are usually filtered out and not reported, on the theory that
	 * they are redundant with the previously-reported alignment.
	 */
	inline void setReportedThrough() {
		reportedThru_ = true;
	}

	// Best incoming score for each 'to' character
	AlnScore best[4];
	// Mask for tied-for-best incoming paths for each 'to' character
	SwColorCellMask mask[4];
	
	// True iff there are no ways to backtrack through this cell as part of a
	// valid alignment.
	bool empty;

	// Initialized to false, set to true once an alignment that moves through
	// the cell is reported.
	bool reportedThru_;
	
	// Initialized to false, set to true once an alignment for which the
	// backtrace begins at this cell is reported.
	bool reportedFrom_[4];

	ASSERT_ONLY(bool finalized);
};

#endif /*ndef ALIGNER_SW_COL_H_*/
