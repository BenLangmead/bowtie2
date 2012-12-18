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

#ifndef ALIGNER_SW_NUC_H_
#define ALIGNER_SW_NUC_H_

#include <stdint.h>
#include "aligner_sw_common.h"
#include "aligner_result.h"

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

enum {
	BT_CAND_FATE_SUCCEEDED = 1,
	BT_CAND_FATE_FAILED,
	BT_CAND_FATE_FILT_START,     // skipped b/c starting cell already explored
	BT_CAND_FATE_FILT_DOMINATED, // skipped b/c it was dominated
	BT_CAND_FATE_FILT_SCORE      // skipped b/c score not interesting anymore
};

/**
 * Encapsulates a cell that we might want to backtrace from.
 */
struct DpBtCandidate {

	DpBtCandidate() { reset(); }
	
	DpBtCandidate(size_t row_, size_t col_, TAlScore score_) {
		init(row_, col_, score_);
	}
	
	void reset() { init(0, 0, 0); }
	
	void init(size_t row_, size_t col_, TAlScore score_) {
		row = row_;
		col = col_;
		score = score_;
		// 0 = invalid; this should be set later according to what happens
		// before / during the backtrace
		fate = 0; 
	}
	
	/** 
	 * Return true iff this candidate is (heuristically) dominated by the given
	 * candidate.  We say that candidate A dominates candidate B if (a) B is
	 * somewhere in the N x N square that extends up and to the left of A,
	 * where N is an arbitrary number like 20, and (b) B's score is <= than
	 * A's.
	 */
	inline bool dominatedBy(const DpBtCandidate& o) {
		const size_t SQ = 40;
		size_t rowhi = row;
		size_t rowlo = o.row;
		if(rowhi < rowlo) swap(rowhi, rowlo);
		size_t colhi = col;
		size_t collo = o.col;
		if(colhi < collo) swap(colhi, collo);
		return (colhi - collo) <= SQ &&
		       (rowhi - rowlo) <= SQ;
	}

	/**
	 * Return true if this candidate is "greater than" (should be considered
	 * later than) the given candidate.
	 */
	bool operator>(const DpBtCandidate& o) const {
		if(score < o.score) return true;
		if(score > o.score) return false;
		if(row   < o.row  ) return true;
		if(row   > o.row  ) return false;
		if(col   < o.col  ) return true;
		if(col   > o.col  ) return false;
		return false;
	}

	/**
	 * Return true if this candidate is "less than" (should be considered
	 * sooner than) the given candidate.
	 */
	bool operator<(const DpBtCandidate& o) const {
		if(score > o.score) return true;
		if(score < o.score) return false;
		if(row   > o.row  ) return true;
		if(row   < o.row  ) return false;
		if(col   > o.col  ) return true;
		if(col   < o.col  ) return false;
		return false;
	}
	
	/**
	 * Return true if this candidate equals the given candidate.
	 */
	bool operator==(const DpBtCandidate& o) const {
		return row   == o.row &&
		       col   == o.col &&
			   score == o.score;
	}
	bool operator>=(const DpBtCandidate& o) const { return !((*this) < o); }
	bool operator<=(const DpBtCandidate& o) const { return !((*this) > o); }
	
#ifndef NDEBUG
	/**
	 * Check internal consistency.
	 */
	bool repOk() const {
		assert(VALID_SCORE(score));
		return true;
	}
#endif

	size_t   row;   // cell row
	size_t   col;   // cell column w/r/t LHS of rectangle
	TAlScore score; // score fo alignment
	int      fate;  // flag indicating whether we succeeded, failed, skipped
};

template <typename T>
class NBest {

public:

	NBest<T>() { nelt_ = nbest_ = n_ = 0; }
	
	bool inited() const { return nelt_ > 0; }
	
	void init(size_t nelt, size_t nbest) {
		nelt_ = nelt;
		nbest_ = nbest;
		elts_.resize(nelt * nbest);
		ncur_.resize(nelt);
		ncur_.fill(0);
		n_ = 0;
	}
	
	/**
	 * Add a new result to bin 'elt'.  Where it gets prioritized in the list of
	 * results in that bin depends on the result of operator>.
	 */
	bool add(size_t elt, const T& o) {
		assert_lt(elt, nelt_);
		const size_t ncur = ncur_[elt];
		assert_leq(ncur, nbest_);
		n_++;
		for(size_t i = 0; i < nbest_ && i <= ncur; i++) {
			if(o > elts_[nbest_ * elt + i] || i >= ncur) {
				// Insert it here
				// Move everyone from here on down by one slot
				for(int j = (int)ncur; j > (int)i; j--) {
					if(j < (int)nbest_) {
						elts_[nbest_ * elt + j] = elts_[nbest_ * elt + j - 1];
					}
				}
				elts_[nbest_ * elt + i] = o;
				if(ncur < nbest_) {
					ncur_[elt]++;
				}
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Return true iff there are no solutions.
	 */
	bool empty() const {
		return n_ == 0;
	}
	
	/**
	 * Dump all the items in our payload into the given EList.
	 */
	template<typename TList>
	void dump(TList& l) const {
		if(empty()) return;
		for(size_t i = 0; i < nelt_; i++) {
			assert_leq(ncur_[i], nbest_);
			for(size_t j = 0; j < ncur_[i]; j++) {
				l.push_back(elts_[i * nbest_ + j]);
			}
		}
	}

protected:

	size_t        nelt_;
	size_t        nbest_;
	EList<T>      elts_;
	EList<size_t> ncur_;
	size_t        n_;     // total # results added
};

#endif /*def ALIGNER_SW_NUC_H_*/
