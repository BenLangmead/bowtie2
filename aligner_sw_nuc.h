/*
 * aligner_sw_nuc.h
 */

#ifndef ALIGNER_SW_NUC_H_
#define ALIGNER_SW_NUC_H_

#include <stdint.h>
#include "aligner_sw_common.h"
#include "aligner_result.h"

/**
 * Encapsulates a bitmask.  The bitmask encodes which backtracking paths out of
 * a cell lie on optimal subpaths.
 */
struct SwNucCellMask {

	typedef uint16_t TMask;

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
		return oall_diag + oall_rfop + oall_rfex + oall_rdop + oall_rdex;
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
		return rfgap_op + rfgap_ex;
	}
	
	/**
	 * Return the number of equally good ways to backtrack from the given table
	 * version of this cell.
	 */
	inline int numPossible(int ct) const {
		if(ct == SW_BT_CELL_OALL) {
			return numOverallPossible();
		} else if(ct == SW_BT_CELL_RDGAP) {
			return numReadGapPossible();
		} else {
			return numRefGapPossible();
		}
	}

	/**
	 * Select a path for backtracking from the oall table version of this cell.
	 * If there is a tie among eligible paths, break it randomly.  Return value
	 * is a flag indicating the backtrack type (see enum defining SW_BT_*
	 * above).
	 */
	int randOverallBacktrack(RandomSource& rand, bool& branch, bool clear);

	/**
	 * Select a path for backtracking from the rdgap table version of this cell.
	 * If there is a tie among eligible paths, break it randomly.  Return value
	 * is a flag indicating the backtrack type (see enum defining SW_BT_*
	 * above).
	 */
	int randReadGapBacktrack(RandomSource& rand, bool& branch, bool clear);

	/**
	 * Select a path for backtracking from the rfgap table version of this cell.
	 * If there is a tie among eligible paths, break it randomly.  Return value
	 * is a flag indicating the backtrack type (see enum defining SW_BT_*
	 * above).
	 */
	int randRefGapBacktrack(RandomSource& rand, bool& branch, bool clear);
	
	/**
	 * Select a path for backtracking from the given table version of this cell.
	 * If there is a tie among eligible paths, break it randomly.  Return value
	 * is a flag indicating the backtrack type (see enum defining SW_BT_*
	 * above).
	 */
	int randBacktrack(int ct, RandomSource& rand, bool& branch, bool clear) {
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
	TMask oall_diag : 1;
	TMask oall_rfop : 1;
	TMask oall_rfex : 1;
	TMask oall_rdop : 1;
	TMask oall_rdex : 1;

	// Read gap (rdgap) table
	TMask rdgap_op  : 1;
	TMask rdgap_ex  : 1;

	// Reference gap (rfgap) table
	TMask rfgap_op  : 1;
	TMask rfgap_ex  : 1;
	
	TMask reserved  : 7;
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
struct DpNucFrame {

	/**
	 * Initialize a new DpNucFrame stack frame.
	 */
	void init(
		size_t   nedsz_,
		size_t   aedsz_,
		size_t   celsz_,
		size_t   row_,
		size_t   col_,
		size_t   gaps_,
		size_t   readGaps_,
		size_t   refGaps_,
		AlnScore score_,
		int      ct_)
	{
		nedsz    = nedsz_;
		aedsz    = aedsz_;
		celsz    = celsz_;
		row      = row_;
		col      = col_;
		gaps     = gaps_;
		readGaps = readGaps_;
		refGaps  = refGaps_;
		score    = score_;
		ct       = ct_;
	}

	size_t   nedsz;    // size of the nucleotide edit list at branch (before
	                   // adding the branch edit)
	size_t   aedsz;    // size of ambiguous nucleotide edit list at branch
	size_t   celsz;    // size of cell-traversed list at branch
	size_t   row;      // row of cell where branch occurred
	size_t   col;      // column of cell where branch occurred
	size_t   gaps;     // number of gaps before branch occurred
	size_t   readGaps; // number of read gaps before branch occurred
	size_t   refGaps;  // number of ref gaps before branch occurred
	AlnScore score;    // score where branch occurred
	int      ct;       // table type (oall, rdgap or rfgap)
};

/**
 * Encapsulates a cell that we might want to backtrace from.
 */
struct DpNucBtCandidate {

	DpNucBtCandidate() { reset(); }
	
	DpNucBtCandidate(size_t row_, size_t col_, TAlScore score_) {
		init(row_, col_, score_);
	}
	
	void reset() { init(0, 0, 0); }
	
	void init(size_t row_, size_t col_, TAlScore score_) {
		row = row_;
		col = col_;
		score = score_;
	}

	/**
	 * Return true if this candidate is "greater than" (should be considered
	 * later than) the given candidate.
	 */
	bool operator>(const DpNucBtCandidate& o) const {
		if(score < o.score) return true;
		if(score > o.score) return false;
		if(row < o.row) return true;
		if(row > o.row) return false;
		if(col < o.col) return true;
		if(col > o.col) return false;
		return false;
	}

	/**
	 * Return true if this candidate is "less than" (should be considered
	 * sooner than) the given candidate.
	 */
	bool operator<(const DpNucBtCandidate& o) const {
		if(score > o.score) return true;
		if(score < o.score) return false;
		if(row > o.row) return true;
		if(row < o.row) return false;
		if(col > o.col) return true;
		if(col < o.col) return false;
		return false;
	}
	
	/**
	 * Return true if this candidate equals the given candidate.
	 */
	bool operator==(const DpNucBtCandidate& o) const {
		return row   == o.row &&
		       col   == o.col &&
			   score == o.score;
	}

	bool operator>=(const DpNucBtCandidate& o) const { return !((*this) < o); }
	bool operator<=(const DpNucBtCandidate& o) const { return !((*this) > o); }
	
	/**
	 * Check internal consistency.
	 */
	bool repOk() const {
		assert(VALID_SCORE(score));
		return true;
	}

	size_t   row;   // cell row
	size_t   col;   // cell column w/r/t LHS of rectangle
	TAlScore score; // score fo alignment
};

/**
 * Encapsulates all information needed to encode the optimal subproblem
 * at a cell in a colorspace SW matrix.
 *
 * Besides scores and masks, we also need to keep track of a few variables that
 * determine how we should backtrack through and into the cell:
 *
 * empty:
 *
 *   Set to true upon reset()/clear().  Set to false in finalize() if both (a)
 *   the mask IS NOT empty, (b) the best score is not less than the score floor.
 *   See 'terminal' below.
 *
 * terminal:
 *
 *   Set to false upon reset()/clear().  Set to true in finalize() if both (a)
 *   the mask IS empty, (b) the best score is not less than the score floor.
 *   See 'empty' above.
 *
 * reportedThru_:
 *
 *   Set to false upon reset()/clear().  Once backtraceNucleotides() has picked
 *   a valid path to backtrace through, it sets reportedThru_ to true for all
 *   cells on the path.
 *
 * backtraceCandidate:
 *
 *   Set to false upon reset()/clear().  In the UPDATE_SOLS macro, it is set to
 *   true iff the cell is one that we might want to backtrace from.  In
 *   end-to-end alignment mode, this means that it must (a) have a score not
 *   less than the minimum, (b) be in the last row, (c) not be in a diagonal
 *   disallowed by the en_ vector.  In local alignment mode, this means that it
 *   must (a) have a score not less than the minimum, (b) be an improvement
 *   over the score diagonally before it, (c) the score diagnoally after cannot
 *   be a further improvement, and (d) cannot be in a diagonal disallowed by
 *   the en_ vector.
 *
 * finalized (debug only):
 *
 *   Set to false upon reset()/clear().  Set to true by user when the mask and
 *   best scores are finished being updated and the 'empty' and 'termina;'
 *   fields have been set appropriately.
 */
struct SwNucCell {

	/**
	 * Clear this cell so that it's ready for updates.
	 */
	void clear() {
		// Initially, best score is invalid
		oallBest = rdgapBest = rfgapBest = AlnScore::INVALID();
		// Initially, there's no way to backtrack from this cell
		mask.clear();
		terminal = false;
		empty = true;
		reportedThru_ = false;
		backtraceCandidate = false;
		assert(mask.empty());
		assert(!oallBest.valid());
		assert(!rdgapBest.valid());
		assert(!rfgapBest.valid());
		ASSERT_ONLY(finalized = false);
		assert(!valid());
		assert(repOk());
	}

	/**
	 * Return true if best is valid.
	 */
	bool valid() const {
		bool val = VALID_AL_SCORE(oallBest);
		//assert(!val || VALID_AL_SCORE(rdgapBest));
		//assert(!val || VALID_AL_SCORE(rfgapBest));
		return val;
	}

	/**
	 * Determine whether this cell has a solution with the given score.  If so,
	 * return true.  Otherwise, return false.
	 */
	inline bool bestSolutionGeq(const TAlScore& min, AlnScore& bst) {
		if(hasSolutionGeq(min)) {
			bst = oallBest;
			return true;
		}
		return false;
	}

	/**
	 * Determine whether this cell has a solution with the given score.  If so,
	 * return true.  Otherwise, return false.
	 */
	inline bool hasSolutionEq(const TAlScore& eq) {
		return !empty && !reportedThru_ && oallBest.score() == eq;
	}

	/**
	 * Determine whether this cell has a solution with a score greater than or
	 * equal to the given score.  If so, return true.  Otherwise, return false.
	 */
	inline bool hasSolutionGeq(const TAlScore& eq) {
		return !empty && !reportedThru_ && oallBest.score() >= eq;
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
	 * Return true iff we can backtrace through this cell.  Called by
	 * backtraceNucleotides() to see if we should abort our current backtrace
	 * without reporting the alignment.
	 */
	inline bool canMoveThrough(int ct) const {
		bool empty = mask.numPossible(ct) == 0;
		bool backtrack = false;
		if(empty) {
			// There were legitimate ways to backtrace from this cell, but
			// they are now foreclosed because they are redundant with
			// alignments already reported
			backtrack = !terminal;
		}
		if(reportedThru_) {
			backtrack = true;
		}
		return !backtrack;
	}

	/**
	 * We finished updating the cell; set empty and finalized
	 * appropriately.
	 */
	inline bool finalize(TAlScore floorsc);
	
	/**
	 * Check that cell is internally consistent
	 */
	bool repOk() const;
	
	/**
	 * Return the best AlnScore field for the given cell type.
	 */
	const AlnScore& best(int ct) const {
		if(ct == SW_BT_CELL_OALL) {
			return oallBest;
		} else if(ct == SW_BT_CELL_RDGAP) {
			return rdgapBest;
		} else {
			assert_eq(SW_BT_CELL_RFGAP, ct);
			return rfgapBest;
		}
	}

	// Best incoming score...
	AlnScore oallBest;  // ...assuming nothing about the incoming transition
	AlnScore rdgapBest; // ...assuming incoming transition is a read gap
	AlnScore rfgapBest; // ...assuming incoming transition is a ref gap
	
	// Mask for tied-for-best incoming paths
	SwNucCellMask mask;
	
	// True iff there are no ways to backtrack through this cell as part of a
	// valid alignment.
	bool empty;

	// True iff the cell's best score is >= the floor but there are no cells to
	// backtrack to.
	bool terminal;
	
	// Initialized to false, set to true once an alignment that moves through
	// the cell is reported.
	bool reportedThru_;

	// Initialized to false, set to true iff it is found to be a candidate cell
	// for starting a bactrace.
	bool backtraceCandidate;

	ASSERT_ONLY(bool finalized);
};

#endif /*def ALIGNER_SW_NUC_H_*/
