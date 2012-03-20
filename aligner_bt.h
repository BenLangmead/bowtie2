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
 *  checkpoints to ensure that only relatively small triangles need to be
 *  refilled.  Note that checkpoints must include elements from the H, E and F
 *  arrays; not just H. After each refill, we backtrace through the refilled
 *  area, then discard/reuse the fill memory.  I call each such fill/backtrace
 *  a mini-fill/backtrace  See following diagram.
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
			assert(cper_->hasEF());
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
 * Encapsulates a "branch" which is a diagonal stretch of the backtrace where
 * all the cells in the stretch are matches.
 *
 * We need BtBranches to "live forever" so that we can make some BtBranches
 * parents of others using parent pointers.
 */
class BtBranch {

public:

	BtBranch() { reset(); }

	BtBranch(
		const BtBranchProblem& prob,
		size_t parentId,
		TAlScore score_en,
		int64_t row,
		int64_t col,
		Edit e)
	{
		init(prob, parentId, score_en, row, col, e);
	}
	
	/**
	 * Reset to uninitialized state.
	 */
	void reset() {
		parentId_ = 0;
		score_st_ = score_en_ = score_best_ = len_ = row_ = col_ = 0;
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
		TAlScore score_en,
		int64_t row,
		int64_t col,
		Edit e);
	
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
		int64_t scoreFloor = prob.sc_->monotone ? std::numeric_limits<int64_t>::min() : 0;
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
		assert_geq(col_ + 1, (int64_t)len_);
		assert_geq(row_ + 1, (int64_t)len_);
		return true;
	}

protected:

	size_t   parentId_;   // pointer to parent branch
	TAlScore score_best_; // best possible score leading through this branch
	TAlScore score_st_;   // score at beginning of branch
	TAlScore score_en_;   // score at end of branch
	size_t   len_;        // length of branch
	int64_t  row_;        // row of lower-right hand cell
	int64_t  col_;        // col of lower-right hand cell
	Edit     e_;          // edit that separates this branch from parent
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

	BtBranchTracer() : prob_(), bs_() { }

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
		bs_.clear();
		if(!prob_.fill_) {
			size_t id = bs_.alloc();
			bs_[id].init(prob_, NULL, 0 /* starting score */, row, col, e);
			if(bs_[id].isSolution(prob_)) {
				addSolution(id);
			} else {
				add(id);
			}
		} else {
			assert(prob_.cper_->doCheckpoints());
			assert(prob_.cper_->hasEF());
			size_t row = row_, col = col_;
			TAlScore targsc = prob_.targ_;
			int hef = 0;
			bool done = false;
			while(!done) {
				// Accumulate edits as we go.  We can do this by adding
				// BtBranches to the bs_ structure.  Each step of the backtrace
				// either involves an edit (thereby starting a new branch) or
				// extends the previous branch by one more position.
				//
				// Note: if the BtBranches are in line, then trySolution can be
				// used to populate the SwResult and check for various
				// situations where we might reject the alignment (i.e. due to
				// a cell having been visited previously).
				cerr << "triangleFill(" << row << ", " << col << ", " << hef << ", " << targsc << ")" << endl;
				triangleFill(
					row,          // row of cell to backtrace from
					col,          // column of cell to backtrace from
					hef,          // cell to backtrace from is H (0), E (1), or F (2)
					targsc,       // score of cell to backtrace from
					rnd,          // pseudo-random generator
					row,          // out: row we ended up in after backtrace
					col,          // out: column we ended up in after backtrace
					hef,          // out: H/E/F after backtrace
					targsc,       // out: score up to cell we ended up in
					done);        // whether we finished tracing out an alignment
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
	 * Fill in a triangle of the DP table and backtrace from the given cell to
	 * a cell in the previous checkpoint, or to the terminal cell.
	 */
	void triangleFill(
		size_t rw,          // row of cell to backtrace from
		size_t cl,          // column of cell to backtrace from
		int hef,            // cell to backtrace from is H (0), E (1), or F (2)
		TAlScore targ,      // score of cell to backtrace from
		RandomSource& rnd,  // pseudo-random generator
		size_t& row_new,    // out: row we ended up in after backtrace
		size_t& col_new,    // out: column we ended up in after backtrace
		int& hef_new,       // out: H/E/F after backtrace
		TAlScore& targ_new, // out: score up to cell we ended up in
		bool& done);        // whether we finished tracing out an alignment

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
		SwResult& res,
		size_t& off,
		size_t& nrej,
		RandomSource& rnd,
		bool& success)
	{
		if(solutions_.size() > 0) {
			for(size_t i = 0; i < solutions_.size(); i++) {
				int ret = trySolution(solutions_[i], res, off, nrej, rnd);
				if(ret == BT_FOUND) {
					success = true;
					return true; // there were solutions and one was good
				}
			}
			solutions_.clear();
			success = false;
			return true; // there were solutions but none were good
		}
		return false; // there were no solutions to check
	}
	
	/**
	 * See if a given solution branch works as a solution (i.e. doesn't overlap
	 * another one, have too many Ns, fail to overlap a core diagonal, etc.)
	 */
	int trySolution(
		size_t id,
		SwResult& res,
		size_t& off,
		size_t& nrej,
		RandomSource& rnd);

	BtBranchProblem    prob_; // problem configuration
	EFactory<BtBranch> bs_;   // global BtBranch factory
	
	// already reported alignments going through these diagonal segments
	ELList<std::pair<size_t, size_t> > seenPaths_;
	
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

	ELList<_CpQuad> tri_;       // triangle to fill when doing mini-fills

#ifndef NDEBUG
	ESet<size_t>  seen_;      // seedn branch ids; should never see same twice
#endif
};

#endif /*ndef ALIGNER_BT_H_*/
