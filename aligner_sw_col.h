/*
 *  aligner_sw_col.h
 */

#ifndef ALIGNER_SW_COL_H_
#define ALIGNER_SW_COL_H_

#include <stdint.h>
#include "aligner_sw_common.h"
#include "aligner_result.h"
#include "mask.h"

/**
 * A bitmask encoding which backtracking paths out of a particular cell
 * correspond to optimal subpaths.  This struct is tailored to the
 * colorspace case where each transition down in the dynamic
 * programming matrix is also associated with a nucleotide assignment
 * to the source (upper) row.
 */
struct SwColorCellMask {

	typedef uint32_t TMask;

	/**
	 * Set all flags to 0 (meaning either: there's no way to backtrack from
	 * this cell to an optimal answer, or we haven't set the mask yet)
	 */
	void clear() {
		*((TMask*)this) = 0;
	}

	/**
	 * Return true iff the mask is empty.
	 */
	inline bool empty() const {
		return *((TMask*)this) == 0;
	}

	/**
	 * Return the number of equally good ways to backtrack from the oall table
	 * version of this cell.
	 */
	inline int numOverallPossible() const {
		return alts5[oall_diag] +
		       alts5[oall_rfop] +
			   alts5[oall_rfex] +
			   oall_rdop +
			   oall_rdex;
	}

	/**
	 * Return the number of equally good ways to backtrack from the rdgap table
	 * version of this cell.
	 */
	inline int numReadGapPossible() const {
		return rdgap_op + rdgap_ex;
	}

	/**
	 * Return the number of equally good ways to backtrack from the rfgap table
	 * version of this cell.
	 */
	inline int numRefGapPossible() const {
		return alts5[rfgap_op] + alts5[rfgap_ex];
	}
	
	/**
	 * Return true iff there is >0 possible way to backtrack from this
	 * cell.
	 */
	inline int numPossible(int ct) const {
		if(ct == SW_BT_CELL_OALL) {
			return numOverallPossible();
		} else if(ct == SW_BT_CELL_RDGAP) {
			return numReadGapPossible();
		} else {
			assert_eq(SW_BT_CELL_RFGAP, ct);
			return numRefGapPossible();
		}
	}

	/**
	 * Select a path for backtracking from this cell.  If there is a
	 * tie among eligible paths, break it randomly.  Return value is
	 * a pair where first = a flag indicating the backtrack type (see
	 * enum defining SW_BT_* above), and second = a selection for
	 * the read character for the next row up.  second should be
	 * ignored if the backtrack type is a gap in the read.
	 */
	std::pair<int, int> randOverallBacktrack(
		RandomSource& rand,
		bool& branch,
		bool clear);

	/**
	 * Select a path for backtracking from this cell.  If there is a
	 * tie among eligible paths, break it randomly.  Return value is
	 * a pair where first = a flag indicating the backtrack type (see
	 * enum defining SW_BT_* above), and second = a selection for
	 * the read character for the next row up.  second should be
	 * ignored if the backtrack type is a gap in the read.
	 */
	std::pair<int, int> randReadGapBacktrack(
		RandomSource& rand,
		bool& branch,
		bool clear);

	/**
	 * Select a path for backtracking from this cell.  If there is a
	 * tie among eligible paths, break it randomly.  Return value is
	 * a pair where first = a flag indicating the backtrack type (see
	 * enum defining SW_BT_* above), and second = a selection for
	 * the read character for the next row up.  second should be
	 * ignored if the backtrack type is a gap in the read.
	 */
	std::pair<int, int> randRefGapBacktrack(
		RandomSource& rand,
		bool& branch,
		bool clear);

	/**
	 * Select a path for backtracking from this cell.  If there is a
	 * tie among eligible paths, break it randomly.  Return value is
	 * a pair where first = a flag indicating the backtrack type (see
	 * enum defining SW_BT_* above), and second = a selection for
	 * the read character for the next row up.  second should be
	 * ignored if the backtrack type is a gap in the read.
	 */
	std::pair<int, int> randBacktrack(
		int ct,
		RandomSource& rand,
		bool& branch,
		bool clear)
	{
		if(ct == SW_BT_CELL_OALL) {
			return randOverallBacktrack(rand, branch, clear);
		} else if(ct == SW_BT_CELL_RDGAP) {
			return randReadGapBacktrack(rand, branch, clear);
		} else {
			assert_eq(SW_BT_CELL_RFGAP, ct);
			return randRefGapBacktrack(rand, branch, clear);
		}
	}

	/**
	 * Clear mask for overall-table cell.
	 */
	inline void clearOverallMask() {
		oall_diag = oall_rfop = oall_rfex = oall_rdop = oall_rdex = 0;
	}

	/**
	 * Clear mask for read-gap-table cell.
	 */
	inline void clearReadGapMask() {
		rdgap_op = rdgap_ex = 0;
	}

	/**
	 * Clear mask for ref-gap-table cell.
	 */
	inline void clearRefGapMask() {
		rfgap_op = rfgap_ex = 0;
	}

	// Overall (oall) table
	TMask oall_diag : 4;
	TMask oall_rfop : 4;
	TMask oall_rfex : 4;
	TMask oall_rdop : 1;
	TMask oall_rdex : 1;

	// Read gap (rdgap) table
	TMask rdgap_op  : 1;
	TMask rdgap_ex  : 1;

	// Reference gap (rfgap) table
	TMask rfgap_op  : 4;
	TMask rfgap_ex  : 4;
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
		size_t   gaps_,
		size_t   readGaps_,
		size_t   refGaps_,
		AlnScore score_,
		int      ct_
		ASSERT_ONLY(, size_t drdsz_)
		)
	{
		nedsz    = nedsz_;
		aedsz    = aedsz_;
		cedsz    = cedsz_;
		celsz    = celsz_;
		row      = row_;
		col      = col_;
		curC     = curC_;
		gaps     = gaps_;
		readGaps = readGaps_;
		refGaps  = refGaps_;
		score    = score_;
		ct       = ct_;
		ASSERT_ONLY(drdsz = drdsz_);
	}

	size_t   nedsz;    // size of the nucleotide edit list at branch (before
	                   // adding the branch edit)
	size_t   aedsz;    // size of ambiguous nucleotide edit list at branch
	size_t   cedsz;    // size of color edit list at branch
	size_t   celsz;    // size of cell-traversed list at branch
	size_t   row;      // row of cell where branch occurred
	size_t   col;      // column of cell where branch occurred
	int      curC;     // character cell we're in
	size_t   gaps;     // gaps before branch occurred
	size_t   readGaps; // read gaps before branch occurred
	size_t   refGaps;  // ref gaps before branch occurred
	AlnScore score;    // score where branch occurred
	int      ct;       // table type (oall, rdgap or rfgap)
	ASSERT_ONLY(size_t drdsz); // size of decoded string
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
		oallBest [0] = oallBest [1] = oallBest [2] = oallBest [3] = AlnScore::INVALID();
		rdgapBest[0] = rdgapBest[1] = rdgapBest[2] = rdgapBest[3] = AlnScore::INVALID();
		rfgapBest[0] = rfgapBest[1] = rfgapBest[2] = rfgapBest[3] = AlnScore::INVALID();
		// Initially, there's no way to backtrack from this cell
		mask[0].clear();
		mask[1].clear();
		mask[2].clear();
		mask[3].clear();
		empty = true;
		terminal[0] = terminal[1] = terminal[2] = terminal[3] = false;
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
			if(VALID_AL_SCORE(oallBest[i])) {
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
	inline bool finalize(TAlScore floorsc);

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
			   oallBest[i].score() >= min &&
			   oallBest[i] > bestSoFar)
			{
				bestSoFar = oallBest[i];
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
			if(!reportedFrom_[i] && oallBest[i].score() == eq) {
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
			if(!reportedFrom_[i] && oallBest[i].score() >= eq) {
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
			if(!reportedFrom_[i] && oallBest[i].score() == eq) {
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
	
	/**
	 * Return true iff we can backtrace through this cell.
	 */
	inline bool canMoveThrough(int curC, int ct) const {
		bool empty = mask[curC].numPossible(ct) == 0;
		bool backtrack = false;
		if(empty) {
			// There were legitimate ways to backtrace from this cell, but
			// they are now foreclosed because they are redundant with
			// alignments already reported
			backtrack = !terminal[curC];
		}
		if(reportedThru_) {
			backtrack = true;
		}
		return !backtrack;
	}

	/**
	 * Return the best AlnScore field for the given cell type.
	 */
	const AlnScore& best(int ch, int ct) const {
		if(ct == SW_BT_CELL_OALL) {
			return oallBest[ch];
		} else if(ct == SW_BT_CELL_RDGAP) {
			return rdgapBest[ch];
		} else {
			assert_eq(SW_BT_CELL_RFGAP, ct);
			return rfgapBest[ch];
		}
	}

	// Best incoming score for each 'to' character...
	AlnScore oallBest[4];  // ...assuming nothing about the incoming transition
	AlnScore rdgapBest[4]; // ...assuming incoming transition is a read gap
	AlnScore rfgapBest[4]; // ...assuming incoming transition is a ref gap
	
	// Mask for tied-for-best incoming paths for each 'to' character
	SwColorCellMask mask[4];
	
	// True iff there are no ways to backtrack through this cell as part of a
	// valid alignment.
	bool empty;

	// True iff the cell's best score is >= the floor but there are no cells to
	// backtrack to.
	bool terminal[4];

	// Initialized to false, set to true once an alignment that moves through
	// the cell is reported.
	bool reportedThru_;
	
	// Initialized to false, set to true once an alignment for which the
	// backtrace begins at this cell is reported.
	bool reportedFrom_[4];

	ASSERT_ONLY(bool finalized);
};

#endif /*ndef ALIGNER_SW_COL_H_*/
