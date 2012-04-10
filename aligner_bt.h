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

#ifndef ALIGNER_BT_H_
#define ALIGNER_BT_H_

#include <utility>
#include <stdint.h>
#include "aligner_sw_common.h"
#include "aligner_result.h"
#include "scoring.h"
#include "edit.h"
#include "limit.h"
#include "dp_framer.h"
#include "sse_util.h"

/* Say we've filled in a DP matrix in a cost-only manner, not saving the scores
 * for each of the cells.  At the end, we obtain a list of candidate cells and
 * we'd like to backtrace from them.  The per-cell scores are gone, but we have
 * to re-create the correct path somehow.  Hopefully we can do this without
 * recreating most or al of the score matrix, since this takes too much memory.
 *
 * Approach 1: Naively refill the matrix.
 *
 *  Just refill the matrix, perhaps backwards starting from the backtrace cell.
 *  Since this involves recreating all or most of the score matrix, this is not
 *  a good approach.
 *
 * Approach 2: Naive backtracking.
 *
 *  Conduct a search through the space of possible backtraces, rooted at the
 *  candidate cell.  To speed things along, we can prioritize paths that have a
 *  high score and that align more characters from the read.
 *
 *  The approach is simple, but it's neither fast nor memory-efficient in
 *  general.
 *
 * Approach 3: Refilling with checkpoints.
 *
 *  Refill the matrix "backwards" starting from the candidate cell, but use
 *  checkpoints to ensure that only a series of relatively small triangles or
 *  rectangles need to be refilled.  The checkpoints must include elements from
 *  the H, E and F matrices; not just H.  After each refill, we backtrace
 *  through the refilled area, then discard/reuse the fill memory.  I call each
 *  such fill/backtrace a mini-fill/backtrace.
 *
 *  If there's only one path to be found, then this is O(m+n).  But what if
 *  there are many?  And what if we would like to avoid paths that overlap in
 *  one or more cells?  There are two ways we can make this more efficient:
 *
 *   1. Remember the re-calculated E/F/H values and try to retrieve them
 *   2. Keep a record of cells that have already been traversed
 *
 *  Legend:
 *
 *  1: Candidate cell
 *  2: Final cell from first mini-fill/backtrace
 *  3: Final cell from second mini-fill/backtrace (third not shown)
 *  +: Checkpointed cell
 *  *: Cell filled from first or second mini-fill/backtrace
 *  -: Unfilled cell
 *
 *        ---++--------++--------++----
 *        --++--------++*-------++-----
 *        -++--(etc)-++**------++------
 *        ++--------+3***-----++-------
 *        +--------++****----++--------
 *        --------++*****---++--------+
 *        -------++******--++--------++
 *        ------++*******-++*-------++-
 *        -----++********++**------++--
 *        ----++********2+***-----++---
 *        ---++--------++****----++----
 *        --++--------++*****---++-----
 *        -++--------++*****1--++------
 *        ++--------++--------++-------
 *
 * Approach 4: Backtracking with checkpoints.
 *
 *  Conduct a search through the space of possible backtraces, rooted at the
 *  candidate cell.  Use "checkpoints" to prune.  That is, when a backtrace
 *  moves through a cell with a checkpointed score, consider the score
 *  accumulated so far and the cell's saved score; abort if those two scores
 *  add to something less than a valid score.  Note we're only checkpointing H
 *  in this case (possibly; see "subtle point"), not E or F.
 *
 *  Subtle point: checkpoint scores are a result of moving forward through
 *  the matrix whereas backtracking scores result from moving backward.  This
 *  matters becuase the two paths that meet up at a cell might have both
 *  factored in a gap open penalty for the same gap, in which case we will
 *  underestimate the overall score and prune a good path.  Here are two ideas
 *  for how to resolve this:
 *
 *   Idea 1: when we combine the forward and backward scores to find an overall
 *   score, and our backtrack procedure *just* made a horizontal or vertical
 *   move, add in a "bonus" equal to the gap open penalty of the appropraite
 *   type (read gap open for horizontal, ref gap open for vertical). This might
 *   overcompensate, since
 *
 *   Idea 2: keep the E and F values for the checkpoints around, in addition to
 *   the H values.  When it comes time to combine the score from the forward
 *   and backward paths, we consider the last move we made in the backward
 *   backtrace.  If it's a read gap (horizontal move), then we calculate the
 *   overall score as:
 *
 *     max(Score-backward + H-forward, Score-backward + E-forward + read-open)
 *
 *   If it's a reference gap (vertical move), then we calculate the overall
 *   score as:
 *
 *     max(Score-backward + H-forward, Score-backward + F-forward + ref-open)
 *
 *   What does it mean to abort a backtrack?  If we're starting a new branch
 *   and there is a checkpoing in the bottommost cell of the branch, and the
 *   overall score is less than the target, then we can simply ignore the
 *   branch.  If the checkpoint occurs in the middle of a string of matches, we
 *   need to curtail the branch such that it doesn't include the checkpointed
 *   cell and we won't ever try to enter the checkpointed cell, e.g., on a
 *   mismatch.
 *
 * Approaches 3 and 4 seem reasonable, and could be combined.  For simplicity,
 * we implement only approach 4 for now.
 *
 * Checkpoint information is propagated from the fill process to the backtracer
 * via a 
 */

enum {
	BT_NOT_FOUND = 1,      // could not obtain the backtrace because it
	                       // overlapped a previous solution
	BT_FOUND,              // obtained a valid backtrace
	BT_REJECTED_N,         // backtrace rejected because it had too many Ns
	BT_REJECTED_CORE_DIAG  // backtrace rejected because it failed to overlap a
	                       // core diagonal
};

/**
 * Parameters for a matrix of potential backtrace problems to solve.
 * Encapsulates information about:
 *
 * The problem given a particular reference substring:
 *
 * - The query string (nucleotides and qualities)
 * - The reference substring (incl. orientation, offset into overall sequence)
 * - Checkpoints (i.e. values of matrix cells)
 * - Scoring scheme and other thresholds
 *
 * The problem given a particular reference substring AND a particular row and
 * column from which to backtrace:
 *
 * - The row and column
 * - The target score
 */
class BtBranchProblem {

public:

	/**
	 * Create new uninitialized problem.
	 */
	BtBranchProblem() { reset(); }

	/**
	 * Initialize a new problem.
	 */
	void initRef(
		const char          *qry,    // query string (along rows)
		const char          *qual,   // query quality string (along rows)
		size_t               qrylen, // query string (along rows) length
		const char          *ref,    // reference string (along columns)
		size_t               reflen, // reference string (along columns) length
		TRefId               refid,  // reference id
		TRefOff              refoff, // reference offset
		bool                 fw,     // orientation of problem
		const DPRect*        rect,   // dynamic programming rectangle filled out
		const Checkpointer*  cper,   // checkpointer
		const Scoring       *sc,     // scoring scheme
		size_t               nceil)  // max # Ns allowed in alignment
	{
		qry_    = qry;
		qual_   = qual;
		qrylen_ = qrylen;
		ref_    = ref;
		reflen_ = reflen;
		refid_  = refid;
		refoff_ = refoff;
		fw_     = fw;
		rect_   = rect;
		cper_   = cper;
		sc_     = sc;
		nceil_  = nceil;
	}

	/**
	 * Initialize a new problem.
	 */
	void initBt(
		size_t   row,   // row
		size_t   col,   // column
		bool     fill,  // use a filling rather than a backtracking strategy
		bool     usecp, // use checkpoints to short-circuit while backtracking
		TAlScore targ)  // target score
	{
		row_    = row;
		col_    = col;
		targ_   = targ;
		fill_   = fill;
		usecp_  = usecp;
		if(fill) {
			assert(usecp_);
		}
	}

	/**
	 * Reset to uninitialized state.
	 */
	void reset() {
		qry_ = qual_ = ref_ = NULL;
		cper_ = NULL;
		rect_ = NULL;
		sc_ = NULL;
		qrylen_ = reflen_ = refid_ = refoff_ = row_ = col_ = targ_ = nceil_ = 0;
		fill_ = fw_ = usecp_ = false;
	}
	
	/**
	 * Return true iff the BtBranchProblem has been initialized.
	 */
	bool inited() const {
		return qry_ != NULL;
	}
	
	/**
	 * Sanity-check the problem.
	 */
	bool repOk() {
		assert_gt(qrylen_, 0);
		assert_gt(reflen_, 0);
		assert_lt(row_, qrylen_);
		assert_lt(col_, reflen_);
		return true;
	}

protected:

	const char         *qry_;    // query string (along rows)
	const char         *qual_;   // query quality string (along rows)
	size_t              qrylen_; // query string (along rows) length
	const char         *ref_;    // reference string (along columns)
	size_t              reflen_; // reference string (along columns) length
	TRefId              refid_;  // reference id
	TRefOff             refoff_; // reference offset
	bool                fw_;     // orientation of problem
	const DPRect*       rect_;   // dynamic programming rectangle filled out
	size_t              row_;    // starting row
	size_t              col_;    // starting column
	TAlScore            targ_;   // target score
	const Checkpointer *cper_;   // checkpointer
	bool                fill_;   // use mini-fills
	bool                usecp_;  // use checkpointing?
	const Scoring      *sc_;     // scoring scheme
	size_t              nceil_;  // max # Ns allowed in alignment
	
	friend class BtBranch;
	friend class BtBranchQ;
	friend class BtBranchTracer;
};

/**
 * Encapsulates a "branch" which is a diagonal of cells (possibly of length 0)
 * in the matrix where all the cells are matches.  These stretches are linked
 * together by edits to form a full backtrace path through the matrix.  Lengths
 * are measured w/r/t to the number of rows traversed by the path, so a branch
 * that represents a read gap extension could have length = 0.
 *
 * At the end of the day, the full backtrace path is represented as a list of
 * BtBranch's where each BtBranch represents a stretch of matching cells (and
 * up to one mismatching cell at its bottom extreme) ending in an edit (or in
 * the bottommost row, in which case the edit is uninitialized).  Each
 * BtBranch's row and col fields indicate the bottommost cell involved in the
 * diagonal stretch of matches, and the len_ field indicates the length of the
 * stretch of matches.  Note that the edits themselves also correspond to
 * movement through the matrix.
 *
 * A related issue is how we record which cells have been visited so that we
 * never report a pair of paths both traversing the same (row, col) of the
 * overall DP matrix.  This gets a little tricky because we have to take into
 * account the cells covered by *edits* in addition to the cells covered by the
 * stretches of matches.  For instance: imagine a mismatch.  That takes up a
 * cell of the DP matrix, but it may or may not be preceded by a string of
 * matches.  It's hard to imagine how to represent this unless we let the
 * mismatch "count toward" the len_ of the branch and let (row, col) refer to
 * the cell where the mismatch occurs.
 *
 * We need BtBranches to "live forever" so that we can make some BtBranches
 * parents of others using parent pointers.  For this reason, BtBranch's are
 * stored in an EFactory object in the BtBranchTracer class.
 */
class BtBranch {

public:

	BtBranch() { reset(); }

	BtBranch(
		const BtBranchProblem& prob,
		size_t parentId,
		TAlScore penalty,
		TAlScore score_en,
		int64_t row,
		int64_t col,
		Edit e,
		int hef,
		bool root,
		bool extend)
	{
		init(prob, parentId, penalty, score_en, row, col, e, hef, root, extend);
	}
	
	/**
	 * Reset to uninitialized state.
	 */
	void reset() {
		parentId_ = 0;
		score_st_ = score_en_ = len_ = row_ = col_ = 0;
		curtailed_ = false;
		e_.reset();
	}
	
	/**
	 * Caller gives us score_en, row and col.  We figure out score_st and len_
	 * by comparing characters from the strings.
	 */
	void init(
		const BtBranchProblem& prob,
		size_t parentId,
		TAlScore penalty,
		TAlScore score_en,
		int64_t row,
		int64_t col,
		Edit e,
		int hef,
		bool root,
		bool extend);
	
	/**
	 * Return true iff this branch ends in a solution to the backtrace problem.
	 */
	bool isSolution(const BtBranchProblem& prob) const {
		const bool end2end = prob.sc_->monotone;
		return score_st_ == prob.targ_ && (!end2end || endsInFirstRow());
	}
	
	/**
	 * Return true iff this branch could potentially lead to a valid alignment.
	 */
	bool isValid(const BtBranchProblem& prob) const {
		int64_t scoreFloor = prob.sc_->monotone ? MIN_I64 : 0;
		if(score_st_ < scoreFloor) {
			// Dipped below the score floor
			return false;
		}
		if(isSolution(prob)) {
			// It's a solution, so it's also valid
			return true;
		}
		if((int64_t)len_ > row_) {
			// Went all the way to the top row
			//assert_leq(score_st_, prob.targ_);
			return score_st_ == prob.targ_;
		} else {
			int64_t match = prob.sc_->match();
			int64_t bonusLeft = (row_ + 1 - len_) * match;
			return score_st_ + bonusLeft >= prob.targ_;
		}
	}
	
	/**
	 * Return true iff this branch overlaps with the given branch.
	 */
	bool overlap(const BtBranchProblem& prob, const BtBranch& bt) const {
		// Calculate this branch's diagonal
		assert_lt(row_, (int64_t)prob.qrylen_);
		size_t fromend = prob.qrylen_ - row_ - 1;
		size_t diag = fromend + col_;
		int64_t lo = 0, hi = row_ + 1;
		if(len_ == 0) {
			lo = row_;
		} else {
			lo = row_ - (len_ - 1);
		}
		// Calculate other branch's diagonal
		assert_lt(bt.row_, (int64_t)prob.qrylen_);
		size_t ofromend = prob.qrylen_ - bt.row_ - 1;
		size_t odiag = ofromend + bt.col_;
		if(diag != odiag) {
			return false;
		}
		int64_t olo = 0, ohi = bt.row_ + 1;
		if(bt.len_ == 0) {
			olo = bt.row_;
		} else {
			olo = bt.row_ - (bt.len_ - 1);
		}
		int64_t losm = olo, hism = ohi;
		if(hi - lo < ohi - olo) {
			swap(lo, losm);
			swap(hi, hism);
		}
		if((lo <= losm && hi > losm) || (lo <  hism && hi >= hism)) {
			return true;
		}
		return false;
	}
	
	/**
	 * Return true iff this branch is higher priority than the branch 'o'.
	 */
	bool operator<(const BtBranch& o) const {
		// Prioritize uppermost above score
		if(uppermostRow() != o.uppermostRow()) {
			return uppermostRow() < o.uppermostRow();
		}
		if(score_st_ != o.score_st_) return score_st_ > o.score_st_;
		if(row_      != o.row_)      return row_ < o.row_;
		if(col_      != o.col_)      return col_ > o.col_;
		if(parentId_ != o.parentId_) return parentId_ > o.parentId_;
		assert(false);
		return false;
	}
	
	/**
	 * Return true iff the topmost cell involved in this branch is in the top
	 * row.
	 */
	bool endsInFirstRow() const {
		assert_leq((int64_t)len_, row_ + 1);
		return (int64_t)len_ == row_+1;
	}
	
	/**
	 * Return the uppermost row covered by this branch.
	 */
	size_t uppermostRow() const {
		assert_geq(row_ + 1, (int64_t)len_);
		return row_ + 1 - (int64_t)len_;
	}

	/**
	 * Return the leftmost column covered by this branch.
	 */
	size_t leftmostCol() const {
		assert_geq(col_ + 1, (int64_t)len_);
		return col_ + 1 - (int64_t)len_;
	}
	
	/**
	 * Sanity-check this BtBranch.
	 */
	bool repOk() const {
		assert(root_ || e_.inited());
		assert_gt(len_, 0);
		assert_geq(col_ + 1, (int64_t)len_);
		assert_geq(row_ + 1, (int64_t)len_);
		return true;
	}

protected:

	// ID of the parent branch.
	size_t   parentId_;

	// Penalty associated with the edit at the bottom of this branch (0 if
	// there is no edit)
	TAlScore penalty_;
	
	// Score at the beginning of the branch
	TAlScore score_st_;
	
	// Score at the end of the branch (taking the edit into account)
	TAlScore score_en_;
	
	// Length of the branch.  That is, the total number of diagonal cells
	// involved in all the matches and in the edit (if any).  Should always be
	// > 0.
	size_t   len_;
	
	// The row of the final (bottommost) cell in the branch.  This might be the
	// bottommost match if the branch has no associated edit.  Otherwise, it's
	// the cell occupied by the edit.
	int64_t  row_;
	
	// The column of the final (bottommost) cell in the branch.
	int64_t  col_;
	
	// The edit at the bottom of the branch.  If this is the bottommost branch
	// in the alignment and it does not end in an edit, then this remains
	// uninitialized.
	Edit     e_;
	
	// True iff this is the bottommost branch in the alignment.  We can't just
	// use row_ to tell us this because local alignments don't necessarily end
	// in the last row.
	bool     root_;
	
	bool     curtailed_;  // true -> pruned at a checkpoint where we otherwise
	                      // would have had a match

friend class BtBranchQ;
friend class BtBranchTracer;

};

/**
 * Instantiate and solve best-first branch-based backtraces.
 */
class BtBranchTracer {

public:

	explicit BtBranchTracer() :
		prob_(), bs_(), seenPaths_(DP_CAT), sawcell_(DP_CAT), doTri_() { }

	/**
	 * Add a branch to the queue.
	 */
	void add(size_t id) {
		assert(!bs_[id].isSolution(prob_));
		unsorted_.push_back(make_pair(bs_[id].score_st_, id));
	}
	
	/**
	 * Add a branch to the list of solutions.
	 */
	void addSolution(size_t id) {
		assert(bs_[id].isSolution(prob_));
		solutions_.push_back(id);
	}

	/**
	 * Given a potential branch to add to the queue, see if we can follow the
	 * branch a little further first.  If it's still valid, or if we reach a
	 * choice between valid outgoing paths, go ahead and add it to the queue.
	 */
	void examineBranch(
		int64_t row,
		int64_t col,
		const Edit& e,
		TAlScore pen,
		TAlScore sc,
		size_t parentId);

	/**
	 * Take all possible ways of leaving the given branch and add them to the
	 * branch queue.
	 */
	void addOffshoots(size_t bid);
	
	/**
	 * Get the best branch and remove it from the priority queue.
	 */
	size_t best(RandomSource& rnd) {
		assert(!empty());
		flushUnsorted();
		assert_gt(sortedSel_ ? sorted1_.size() : sorted2_.size(), cur_);
		// Perhaps shuffle everyone who's tied for first?
		size_t id = sortedSel_ ? sorted1_[cur_] : sorted2_[cur_];
		cur_++;
		return id;
	}
	
	/**
	 * Return true iff there are no branches left to try.
	 */
	bool empty() const {
		return size() == 0;
	}
	
	/**
	 * Return the size, i.e. the total number of branches contained.
	 */
	size_t size() const {
		return unsorted_.size() +
		       (sortedSel_ ? sorted1_.size() : sorted2_.size()) - cur_;
	}

	/**
	 * Return true iff there are no solutions left to try.
	 */
	bool emptySolution() const {
		return sizeSolution() == 0;
	}
	
	/**
	 * Return the size of the solution set so far.
	 */
	size_t sizeSolution() const {
		return solutions_.size();
	}
	
	/**
	 * Sort unsorted branches, merge them with master sorted list.
	 */
	void flushUnsorted();
	
	/**
	 * Sanity-check the queue.
	 */
	bool repOk() const {
		assert_lt(cur_, (sortedSel_ ? sorted1_.size() : sorted2_.size()));
		return true;
	}
	
	/**
	 * Initialize the tracer with respect to a new read.  This involves
	 * resetting all the state relating to the set of cells already visited
	 */
	void initRef(
		const char*         rd,     // in: read sequence
		const char*         qu,     // in: quality sequence
		size_t              rdlen,  // in: read sequence length
		const char*         rf,     // in: reference sequence
		size_t              rflen,  // in: reference sequence length
		TRefId              refid,  // in: reference id
		TRefOff             refoff, // in: reference offset
		bool                fw,     // in: orientation
		const DPRect       *rect,   // in: DP rectangle
		const Checkpointer *cper,   // in: checkpointer
		const Scoring&      sc,     // in: scoring scheme
		size_t              nceil)  // in: N ceiling
	{
		prob_.initRef(rd, qu, rdlen, rf, rflen, refid, refoff, fw, rect, cper, &sc, nceil);
		const size_t ndiag = rflen + rdlen - 1;
		seenPaths_.resize(ndiag);
		for(size_t i = 0; i < ndiag; i++) {
			seenPaths_[i].clear();
		}
		// clear each of the per-column sets
		if(sawcell_.size() < rflen) {
			size_t isz = sawcell_.size();
			sawcell_.resize(rflen);
			for(size_t i = isz; i < rflen; i++) {
				sawcell_[i].setCat(DP_CAT);
			}
		}
		for(size_t i = 0; i < rflen; i++) {
			sawcell_[i].setCat(DP_CAT);
			sawcell_[i].clear(); // clear the set
		}
	}
	
	/**
	 * Initialize with a new backtrace.
	 */
	void initBt(
		TAlScore       escore, // in: alignment score
		size_t         row,    // in: start in this row
		size_t         col,    // in: start in this column
		bool           fill,   // in: use mini-filling?
		bool           usecp,  // in: use checkpointing?
		bool           doTri,  // in: triangle-shaped mini-fills?
		RandomSource&  rnd)    // in: random gen, to choose among equal paths
	{
		prob_.initBt(row, col, fill, usecp, escore);
		Edit e; e.reset();
		unsorted_.clear();
		solutions_.clear();
		sorted1_.clear();
		sorted2_.clear();
		cur_ = 0;
		nmm_ = 0;         // number of mismatches attempted
		nnmm_ = 0;        // number of mismatches involving N attempted
		nrdop_ = 0;       // number of read gap opens attempted
		nrfop_ = 0;       // number of ref gap opens attempted
		nrdex_ = 0;       // number of read gap extensions attempted
		nrfex_ = 0;       // number of ref gap extensions attempted
		nmmPrune_ = 0;    // number of mismatches attempted
		nnmmPrune_ = 0;   // number of mismatches involving N attempted
		nrdopPrune_ = 0;  // number of read gap opens attempted
		nrfopPrune_ = 0;  // number of ref gap opens attempted
		nrdexPrune_ = 0;  // number of read gap extensions attempted
		nrfexPrune_ = 0;  // number of ref gap extensions attempted
		row_ = row;
		col_ = col;
		doTri_ = doTri;
		bs_.clear();
		if(!prob_.fill_) {
			size_t id = bs_.alloc();
			bs_[id].init(
				prob_,
				0,     // parent id
				0,     // penalty
				0,     // starting score
				row,   // row
				col,   // column
				e,
				0,
			    true,  // this is the root
				true); // this should be extend with exact matches
			if(bs_[id].isSolution(prob_)) {
				addSolution(id);
			} else {
				add(id);
			}
		} else {
			int64_t row = row_, col = col_;
			TAlScore targsc = prob_.targ_;
			int hef = 0;
			bool done = false, abort = false;
			size_t depth = 0;
			while(!done && !abort) {
				// Accumulate edits as we go.  We can do this by adding
				// BtBranches to the bs_ structure.  Each step of the backtrace
				// either involves an edit (thereby starting a new branch) or
				// extends the previous branch by one more position.
				//
				// Note: if the BtBranches are in line, then trySolution can be
				// used to populate the SwResult and check for various
				// situations where we might reject the alignment (i.e. due to
				// a cell having been visited previously).
				if(doTri_) {
					triangleFill(
						row,          // row of cell to backtrace from
						col,          // column of cell to backtrace from
						hef,          // cell to bt from: H (0), E (1), or F (2)
						targsc,       // score of cell to backtrace from
						prob_.targ_,  // score of alignment we're looking for
						rnd,          // pseudo-random generator
						row,          // out: row we ended up in after bt
						col,          // out: column we ended up in after bt
						hef,          // out: H/E/F after backtrace
						targsc,       // out: score up to cell we ended up in
						done,         // out: finished tracing out an alignment?
						abort);       // out: aborted b/c cell was seen before?
				} else {
					squareFill(
						row,          // row of cell to backtrace from
						col,          // column of cell to backtrace from
						hef,          // cell to bt from: H (0), E (1), or F (2)
						targsc,       // score of cell to backtrace from
						prob_.targ_,  // score of alignment we're looking for
						rnd,          // pseudo-random generator
						row,          // out: row we ended up in after bt
						col,          // out: column we ended up in after bt
						hef,          // out: H/E/F after backtrace
						targsc,       // out: score up to cell we ended up in
						done,         // out: finished tracing out an alignment?
						abort);       // out: aborted b/c cell was seen before?
				}
				if(depth >= ndep_.size()) {
					ndep_.resize(depth+1);
					ndep_[depth] = 1;
				} else {
					ndep_[depth]++;
				}
				depth++;
				assert((row >= 0 && col >= 0) || done);
			}
		}
		ASSERT_ONLY(seen_.clear());
	}
	
	/**
	 * Get the next valid alignment given the backtrace problem.  Return false
	 * if there is no valid solution, e.g., if 
	 */
	bool nextAlignment(
		size_t maxiter,
		SwResult& res,
		size_t& off,
		size_t& nrej,
		size_t& niter,
		RandomSource& rnd);
	
	/**
	 * Return true iff this tracer has been initialized
	 */
	bool inited() const {
		return prob_.inited();
	}
	
	/**
	 * Return true iff the mini-fills are triangle-shaped.
	 */
	bool doTri() const { return doTri_; }

	/**
	 * Fill in a triangle of the DP table and backtrace from the given cell to
	 * a cell in the previous checkpoint, or to the terminal cell.
	 */
	void triangleFill(
		int64_t rw,          // row of cell to backtrace from
		int64_t cl,          // column of cell to backtrace from
		int hef,             // cell to backtrace from is H (0), E (1), or F (2)
		TAlScore targ,       // score of cell to backtrace from
		TAlScore targ_final, // score of alignment we're looking for
		RandomSource& rnd,   // pseudo-random generator
		int64_t& row_new,    // out: row we ended up in after backtrace
		int64_t& col_new,    // out: column we ended up in after backtrace
		int& hef_new,        // out: H/E/F after backtrace
		TAlScore& targ_new,  // out: score up to cell we ended up in
		bool& done,          // out: finished tracing out an alignment?
		bool& abort);        // out: aborted b/c cell was seen before?

	/**
	 * Fill in a square of the DP table and backtrace from the given cell to
	 * a cell in the previous checkpoint, or to the terminal cell.
	 */
	void squareFill(
		int64_t rw,          // row of cell to backtrace from
		int64_t cl,          // column of cell to backtrace from
		int hef,             // cell to backtrace from is H (0), E (1), or F (2)
		TAlScore targ,       // score of cell to backtrace from
		TAlScore targ_final, // score of alignment we're looking for
		RandomSource& rnd,   // pseudo-random generator
		int64_t& row_new,    // out: row we ended up in after backtrace
		int64_t& col_new,    // out: column we ended up in after backtrace
		int& hef_new,        // out: H/E/F after backtrace
		TAlScore& targ_new,  // out: score up to cell we ended up in
		bool& done,          // out: finished tracing out an alignment?
		bool& abort);        // out: aborted b/c cell was seen before?

protected:

	/**
	 * Get the next valid alignment given a backtrace problem.  Return false
	 * if there is no valid solution.  Use a backtracking search to find the
	 * solution.  This can be very slow.
	 */
	bool nextAlignmentBacktrace(
		size_t maxiter,
		SwResult& res,
		size_t& off,
		size_t& nrej,
		size_t& niter,
		RandomSource& rnd);

	/**
	 * Get the next valid alignment given a backtrace problem.  Return false
	 * if there is no valid solution.  Use a triangle-fill backtrace to find
	 * the solution.  This is usually fast (it's O(m + n)).
	 */
	bool nextAlignmentFill(
		size_t maxiter,
		SwResult& res,
		size_t& off,
		size_t& nrej,
		size_t& niter,
		RandomSource& rnd);

	/**
	 * Try all the solutions accumulated so far.  Solutions might be rejected
	 * if they, for instance, overlap a previous solution, have too many Ns,
	 * fail to overlap a core diagonal, etc.
	 */
	bool trySolutions(
		bool lookForOlap,
		SwResult& res,
		size_t& off,
		size_t& nrej,
		RandomSource& rnd,
		bool& success);
	
	/**
	 * See if a given solution branch works as a solution (i.e. doesn't overlap
	 * another one, have too many Ns, fail to overlap a core diagonal, etc.)
	 */
	int trySolution(
		size_t id,
		bool lookForOlap,
		SwResult& res,
		size_t& off,
		size_t& nrej,
		RandomSource& rnd);

	BtBranchProblem    prob_; // problem configuration
	EFactory<BtBranch> bs_;   // global BtBranch factory
	
	// already reported alignments going through these diagonal segments
	ELList<std::pair<size_t, size_t> > seenPaths_;
	ELSet<size_t> sawcell_; // cells already backtraced through
	
	EList<std::pair<TAlScore, size_t> > unsorted_;  // unsorted list of as-yet-unflished BtBranches
	EList<size_t> sorted1_;   // list of BtBranch, sorted by score
	EList<size_t> sorted2_;   // list of BtBranch, sorted by score
	EList<size_t> solutions_; // list of solution branches
	bool          sortedSel_; // true -> 1, false -> 2
	size_t        cur_;       // cursor into sorted list to start from
	
	size_t        nmm_;         // number of mismatches attempted
	size_t        nnmm_;        // number of mismatches involving N attempted
	size_t        nrdop_;       // number of read gap opens attempted
	size_t        nrfop_;       // number of ref gap opens attempted
	size_t        nrdex_;       // number of read gap extensions attempted
	size_t        nrfex_;       // number of ref gap extensions attempted
	
	size_t        nmmPrune_;    // 
	size_t        nnmmPrune_;   // 
	size_t        nrdopPrune_;  // 
	size_t        nrfopPrune_;  // 
	size_t        nrdexPrune_;  // 
	size_t        nrfexPrune_;  // 
	
	size_t        row_;         // row
	size_t        col_;         // column

	bool           doTri_;      // true -> fill in triangles; false -> squares
	EList<CpQuad>  sq_;         // square to fill when doing mini-fills
	ELList<CpQuad> tri_;        // triangle to fill when doing mini-fills
	EList<size_t>  ndep_;       // # triangles mini-filled at various depths

#ifndef NDEBUG
	ESet<size_t>  seen_;        // seedn branch ids; should never see same twice
#endif
};

#endif /*ndef ALIGNER_BT_H_*/
