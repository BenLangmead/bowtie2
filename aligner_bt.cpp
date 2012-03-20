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

#include "aligner_bt.h"
#include "mask.h"
#include <algorithm>

using namespace std;

/**
 * Fill in a triangle of the DP table and backtrace from the given cell to
 * a cell in the previous checkpoint, or to the terminal cell.
 */
void BtBranchTracer::triangleFill(
	size_t rw,          // row of cell to backtrace from
	size_t cl,          // column of cell to backtrace from
	int hef,            // cell to backtrace from is H (0), E (1), or F (2)
	TAlScore targ,      // score of cell to backtrace from
	RandomSource& rnd,  // pseudo-random generator
	size_t& row_new,    // out: row we ended up in after backtrace
	size_t& col_new,    // out: column we ended up in after backtrace
	int& hef_new,       // out: H/E/F after backtrace
	TAlScore& targ_new, // out: score up to cell we ended up in
	bool& done)         // whether we finished tracing out an alignment
{
	assert_range(0, 2, hef);
	assert_lt(rw, prob_.qrylen_);
	assert_lt(cl, prob_.reflen_);
	assert(prob_.usecp_ && prob_.fill_ && prob_.cper_->doCheckpoints());
	int64_t row = rw, col = cl;
	const int64_t colmin = 0;
	const int64_t rowmin = 0;
	const int64_t colmax = prob_.reflen_ - 1;
	const int64_t rowmax = prob_.qrylen_ - 1;
	assert_leq(col, prob_.cper_->hicol());
	assert_geq(col, prob_.cper_->locol());
	assert_geq(prob_.cper_->per(), 2);
	size_t mod = (row + col) & prob_.cper_->lomask();
	assert_lt(mod, prob_.cper_->per());
	// Allocate room for diags
	size_t depth = mod+1;
	assert_leq(depth, prob_.cper_->per());
	size_t breadth = depth;
	tri_.resize(depth);
	// Allocate room for each diag
	for(size_t i = 0; i < depth; i++) {
		tri_[i].resize(breadth - i);
	}
	bool upperleft = false;
	size_t off = (row + col) >> prob_.cper_->perpow2();
	if(off == 0) {
		upperleft = true;
	} else {
		off--;
	}
	const TAlScore sc_rdo = prob_.sc_->readGapOpen();
	const TAlScore sc_rde = prob_.sc_->readGapExtend();
	const TAlScore sc_rfo = prob_.sc_->refGapOpen();
	const TAlScore sc_rfe = prob_.sc_->refGapExtend();
	const bool local = !prob_.sc_->monotone;
	int64_t row_lo = row - (int64_t)mod;
	const _CpQuad *prev2 = NULL, *prev1 = NULL;
	ASSERT_ONLY(TAlScore bestFromCp = std::numeric_limits<TAlScore>::min());
	if(!upperleft) {
		// Read-only pointer to cells in diagonal -2.  Start one row above the
		// target row.
		prev2 = prob_.cper_->qdiag1sPtr() + (off * prob_.cper_->nrow() + row_lo - 1);
		// Read-only pointer to cells in diagonal -1.  Start one row above the
		// target row
		prev1 = prob_.cper_->qdiag2sPtr() + (off * prob_.cper_->nrow() + row_lo - 1);
#ifndef NDEBUG
		for(size_t i = 0; i < depth+1; i++) {
			bestFromCp = max<TAlScore>(bestFromCp, prev1[i].sc[0]);
			if(i < depth) {
				bestFromCp = max<TAlScore>(bestFromCp, prev2[i].sc[0]);
			}
		}
		size_t rowc = row - mod, colc = col;
		if(prob_.cper_->isCheckpointed(rowc-1, colc)) {
			assert_eq(prob_.cper_->score(rowc-1, colc, 0), prev1[0].sc[0]);
		}
		if(prob_.cper_->isCheckpointed(rowc-1, colc-1)) {
			assert_eq(prob_.cper_->score(rowc-1, colc-1, 0), prev2[0].sc[0]);
		}
#endif
	}
	// Pointer to cells in current diagonal
	// For each diagonal we need to fill in
	for(size_t i = 0; i < depth; i++) {
		_CpQuad * cur = tri_[i].ptr();
		_CpQuad * curc = cur;
		size_t doff = mod - i; // # diagonals we are away from target diag
		assert_geq(row, doff);
		int64_t rowc = row - doff;
		int64_t colc = col;
		size_t neval = 0; // # cells evaluated in this diag
		const _CpQuad *last = NULL;
		// Fill this diagonal from upper right to lower left
		for(size_t j = 0; j < breadth; j++) {
			if(rowc >= rowmin && rowc <= rowmax &&
			   colc >= colmin && colc <= colmax)
			{
				neval++;
				int64_t fromend = prob_.qrylen_ - rowc - 1;
				bool allowGaps = fromend >= prob_.sc_->gapbar && rowc >= prob_.sc_->gapbar;
				// Fill this cell
				// Some things we might want to calculate about this cell up front:
				// 1. How many matches are possible from this cell to the cell in
				//    row, col, in case this allows us to prune
				// Get character from read
				int qc = prob_.qry_[rowc];
				// Get quality value from read
				int qq = prob_.qual_[rowc];
				assert_geq(qq, 33);
				// Get character from reference
				int rc = prob_.ref_[colc];
				assert_range(0, 16, rc);
				TAlScore sc_diag = prob_.sc_->score(qc, rc, qq - 33);
				TAlScore sc_h_up = std::numeric_limits<TAlScore>::min();
				TAlScore sc_f_up = std::numeric_limits<TAlScore>::min();
				TAlScore sc_h_lf = std::numeric_limits<TAlScore>::min();
				TAlScore sc_e_lf = std::numeric_limits<TAlScore>::min();
				if(allowGaps) {
					if(rowc > 0) {
						sc_h_up = prev1[j+0].sc[0] - sc_rfo;
						sc_f_up = prev1[j+0].sc[2] - sc_rfe;
					}
					if(colc > 0) {
						sc_h_lf = prev1[j+1].sc[0] - sc_rdo;
						sc_e_lf = prev1[j+1].sc[1] - sc_rde;
					}
				}
				TAlScore sc_h_dg =
					(rowc > 0 && colc > 0) ? (prev2[j+0].sc[0] + sc_diag) :
					                         std::numeric_limits<TAlScore>::min();
				if(local) {
					// Local-alignment clamping
					sc_h_up = max<TAlScore>(sc_h_up, 0);
					sc_f_up = max<TAlScore>(sc_f_up, 0);
					sc_h_lf = max<TAlScore>(sc_h_lf, 0);
					sc_e_lf = max<TAlScore>(sc_e_lf, 0);
					sc_h_dg = max<TAlScore>(sc_h_dg, 0);
				}
				int mask = 1;
				// Calculate best ways into H, E, F cells starting with H.
				// Mask bits:
				// H: 1=diag, 2=hhoriz, 4=ehoriz, 8=hvert, 16=fvert
				// E: 32=hhoriz, 64=ehoriz
				// F: 128=hvert, 256=fvert
				TAlScore sc_best = sc_h_dg;
				if(sc_h_lf >= sc_best) {
					if(sc_h_lf > sc_best) mask = 0;
					mask |= 2;
					sc_best = sc_h_lf;
				}
				if(sc_e_lf >= sc_best) {
					if(sc_e_lf > sc_best) mask = 0;
					mask |= 4;
					sc_best = sc_e_lf;
				}
				if(sc_h_up >= sc_best) {
					if(sc_h_up > sc_best) mask = 0;
					mask |= 8;
					sc_best = sc_h_up;
				}
				if(sc_f_up >= sc_best) {
					if(sc_f_up > sc_best) mask = 0;
					mask |= 16;
					sc_best = sc_f_up;
				}
				assert(local || sc_best <= bestFromCp);
				// Calculate best way into E cell
				TAlScore sc_e_best = sc_h_lf;
				if(sc_h_lf >= sc_e_lf) {
					if(sc_h_lf == sc_e_lf) {
						mask |= 64;
					}
					mask |= 32;
				} else {
					sc_e_best = sc_e_lf;
					mask |= 64;
				}
				assert(local || sc_e_best <= bestFromCp);
				// Calculate best way into F cell
				TAlScore sc_f_best = sc_h_up;
				if(sc_h_up >= sc_f_up) {
					if(sc_h_up == sc_f_up) {
						mask |= 256;
					}
					mask |= 128;
				} else {
					sc_f_best = sc_f_up;
					mask |= 256;
				}
				assert(local || sc_f_best <= bestFromCp);
				// Install results in cur
				assert(!prob_.sc_->monotone || sc_best <= 0);
				assert(!prob_.sc_->monotone || sc_e_best <= 0);
				assert(!prob_.sc_->monotone || sc_f_best <= 0);
				assert( prob_.sc_->monotone || sc_best >= 0);
				assert( prob_.sc_->monotone || sc_e_best >= 0);
				assert( prob_.sc_->monotone || sc_f_best >= 0);
				curc->sc[0] = sc_best;
				curc->sc[1] = sc_e_best;
				curc->sc[2] = sc_f_best;
#ifndef NDEBUG
				if(prob_.cper_->isCheckpointed(rowc, colc)) {
					assert_eq(prob_.cper_->score(rowc, colc, 0), sc_best);
					assert_eq(prob_.cper_->score(rowc, colc, 1), sc_e_best);
					assert_eq(prob_.cper_->score(rowc, colc, 2), sc_f_best);
				}
#endif
				curc->sc[3] = mask;
				last = curc;
			}
			// Update row, col
			assert_gt(colc, 0);
			assert_lt(rowc, prob_.qrylen_);
			rowc++;
			colc--;
			curc++;
		}
		if(i == depth-1) {
			// Final iteration
			assert(last != NULL);
			assert_eq(1, neval);
			assert_eq(targ, last->sc[hef]);
			assert_neq(0, last->sc[3]);
		} else {
			breadth--;
			prev2 = prev1 + 1;
			prev1 = cur;
		}
	}
	
	//if(bs_.empty()) {
	//	// Start an initial branch
	//	size_t id = bs_.alloc();
	//}
	// Now backtrack through the triangle
	size_t rowc = row, colc = col;
	int hefc = hef;
	size_t idx_orig = (row + col) >> prob_.cper_->perpow2();
	while(true) {
		// What depth are we?
		size_t mod = (rowc + colc) & prob_.cper_->lomask();
		assert_lt(mod, prob_.cper_->per());
		_CpQuad * cur = tri_[mod].ptr();
		assert(!local || cur->sc[0] > 0);
		int64_t row_off = rowc - row_lo - mod;
		assert_geq(row_off, 0);
		int mask = cur[row_off].sc[3];
		assert_gt(mask, 0);
		int sel = -1;
		if(hefc == 0) {
			mask &= 31;
			sel = randFromMask(rnd, mask);
		} else if(hefc == 1) {
			mask >>= 5;
			mask &= 3;
			sel = randFromMask(rnd, mask);
		} else {
			assert_eq(2, hefc);
			mask >>= 7;
			mask &= 3;
			sel = randFromMask(rnd, mask);
		}
		assert_geq(sel, 0);
		if(sel == 0) {
			hefc = 0;
			assert_gt(rowc, 0);
			assert_gt(colc, 0);
			rowc--;
			colc--;
		} else if((sel >= 1 && sel <= 2) || (sel >= 5 && sel <= 6)) {
			hefc = 1;
			assert_gt(colc, 0);
			colc--;
		} else {
			hefc = 2;
			assert_gt(rowc, 0);
			rowc--;
		}
		size_t mod_new = (rowc + colc) & prob_.cper_->lomask();
		size_t idx = (rowc + colc) >> prob_.cper_->perpow2();
		assert_lt(mod_new, prob_.cper_->per());
		_CpQuad * cur_new = tri_[mod_new].ptr();
		int64_t row_off_new = rowc - row_lo - mod_new;
		// Check whether we made it to the top row or to 
		if(rowc == 0 || (local && cur->sc[0] == 0)) {
			row_new = rowc; col_new = colc;
			done = true;
			targ_new = cur_new[row_off_new].sc[hefc];
			hef_new = hefc;
			return;
		}
		if(idx < idx_orig) {
			assert(prob_.cper_->isCheckpointed(rowc, colc));
			row_new = rowc; col_new = colc;
			done = false;
			targ_new = prob_.cper_->score(rowc, colc, hefc);
			hef_new = hefc;
			return;
		}
	}
	assert(false);
}

/**
 * Caller gives us score_en, row and col.  We figure out score_st and len_
 * by comparing characters from the strings.
 */
void BtBranch::init(
	const BtBranchProblem& prob,
	size_t parentId,
	TAlScore score_en,
	int64_t row,
	int64_t col,
	Edit e)
{
	// TODO: if the first cell is no good w/r/t checkpoint score, we should
	// have aborted by now.  In this function we're only checking if cells that
	// we would otherwise "match through" can be pruned.
	score_en_ = score_en;
	score_st_ = score_en_;
	row_ = row;
	col_ = col;
	parentId_ = parentId;
	e_ = e;
	assert_lt(row, (int64_t)prob.qrylen_);
	assert_lt(col, (int64_t)prob.reflen_);
	int64_t rowc = row, colc = col;
	len_ = 0;
	int64_t match = prob.sc_->match();
	bool cp = prob.usecp_ && prob.cper_->doCheckpoints(); // Are there are any checkpoints?
	size_t iters = 0;
	curtailed_ = false;
	while(rowc >= 0 && colc >= 0) {
		int rfm = prob.ref_[colc];
		assert_range(0, 16, rfm);
		int rdc = prob.qry_[rowc];
		//int rdq = prob.qual_[rowc];
		bool matches = (rfm & (1 << rdc)) != 0;
		if(!matches) {
			// What's the mismatch penalty?
			//int sc = prob.sc_->score(rdc, rfm, rdq - 33);
			break;
		}
		// Get score from checkpointer
		score_st_ += match;
		if(cp && rowc - 1 >= 0 && colc - 1 >= 0 &&
		   prob.cper_->isCheckpointed(rowc - 1, colc - 1))
		{
			// Possibly prune
			int16_t cpsc = prob.cper_->hScore(rowc - 1, colc - 1);
			if(cpsc + score_st_ < prob.targ_) {
				curtailed_ = true;
				break;
			}
		}
		iters++;
		rowc--; colc--;
	}
	assert_geq(rowc, -1);
	assert_geq(colc, -1);
	len_ = (int64_t)row - rowc;
	assert_leq((int64_t)len_, row_+1);
	assert_leq((int64_t)len_, col_+1);
	assert_leq((int64_t)score_st_, (int64_t)prob.qrylen_ * match);
	//if(len_ > row_) {
	//	assert_leq(score_st_, prob.targ_);
	//}
	if(uppermostRow() > 0) {
		// Also, check if we're doing so well score-wise that no possible gap
		// or mismatch could drive us 
	
		// If root edit is a gap, perhaps explore extensions?
		
		// Calculate a somewhat tighter upper bound on the score we could get
		// by including this branch on our path. We do this by taking the
		// smallest penalty of any of the penalties that could be incurred upon
		// leaving this branch.
	}
}

/**
 * Given a potential branch to add to the queue, see if we can follow the
 * branch a little further first.  If it's still valid, or if we reach a
 * choice between valid outgoing paths, go ahead and add it to the queue.
 */
void BtBranchTracer::examineBranch(
	int64_t row,
	int64_t col,
	const Edit& e,
	TAlScore sc,
	size_t parentId)
{
	size_t id = bs_.alloc();
	bs_[id].init(prob_, parentId, sc, row, col, e);
	if(bs_[id].isSolution(prob_)) {
		assert(bs_[id].isValid(prob_));
		addSolution(id);
	} else {
		// Check if this branch is legit
		if(bs_[id].isValid(prob_)) {
			add(id);
		} else {
			bs_.pop();
		}
	}
}

/**
 * Take all possible ways of leaving the given branch and add them to the
 * branch queue.
 */
void BtBranchTracer::addOffshoots(size_t bid) {
	BtBranch& b = bs_[bid];
	TAlScore sc = b.score_en_;
	int64_t match = prob_.sc_->match();
	int64_t scoreFloor = prob_.sc_->monotone ? std::numeric_limits<int64_t>::min() : 0;
	bool cp = prob_.usecp_ && prob_.cper_->doCheckpoints(); // Are there are any checkpoints?
	ASSERT_ONLY(TAlScore perfectScore = prob_.sc_->perfectScore(prob_.qrylen_));
	assert_leq(prob_.targ_, perfectScore);
	// For each cell in the branch
	for(size_t i = 0 ; i <= b.len_; i++) {
		assert_leq((int64_t)i, b.row_+1);
		assert_leq((int64_t)i, b.col_+1);
		int64_t row = b.row_ - i, col = b.col_ - i;
		int64_t bonusLeft = (row + 1) * match;
		int64_t fromend = prob_.qrylen_ - row - 1;
		bool allowGaps = fromend >= prob_.sc_->gapbar && row >= prob_.sc_->gapbar;
		if(allowGaps && row >= 0 && col >= 0) {
			if(col > 0) {
				// Try a read gap - it's either an extension or an open
				bool extend = b.e_.inited() && b.e_.isReadGap() && i == 0;
				TAlScore rdgapPen = extend ?
					prob_.sc_->readGapExtend() : prob_.sc_->readGapOpen();
				bool prune = false;
				assert_gt(rdgapPen, 0);
				if(cp && prob_.cper_->isCheckpointed(row, col - 1)) {
					// Possibly prune
					int16_t cpsc = prob_.cper_->hScore(row, col - 1);
					assert_leq(cpsc, perfectScore);
					assert_geq(prob_.sc_->readGapOpen(), prob_.sc_->readGapExtend());
					TAlScore bonus = prob_.sc_->readGapOpen() - prob_.sc_->readGapExtend();
					assert_geq(bonus, 0);
					if(cpsc + bonus + sc - rdgapPen < prob_.targ_) {
						prune = true;
					}
				}
				if(prune) {
					if(extend) { nrdexPrune_++; } else { nrdopPrune_++; }
				} else if(sc - rdgapPen >= scoreFloor && sc - rdgapPen + bonusLeft >= prob_.targ_) {
					// Yes, we can introduce a read gap here
					Edit e((int)row, mask2dna[(int)prob_.ref_[col]], '-', EDIT_TYPE_READ_GAP);
					assert(e.isReadGap());
					examineBranch(row, col - 1, e, sc - rdgapPen, bid);
					if(extend) { nrdex_++; } else { nrdop_++; }
				}
			}
			if(row > 0) {
				// Try a reference gap - it's either an extension or an open
				bool extend = b.e_.inited() && b.e_.isRefGap() && i == 0;
				TAlScore rfgapPen = (b.e_.inited() && b.e_.isRefGap()) ?
					prob_.sc_->refGapExtend() : prob_.sc_->refGapOpen();
				bool prune = false;
				assert_gt(rfgapPen, 0);
				if(cp && prob_.cper_->isCheckpointed(row - 1, col)) {
					// Possibly prune
					int16_t cpsc = prob_.cper_->hScore(row - 1, col);
					assert_leq(cpsc, perfectScore);
					assert_geq(prob_.sc_->refGapOpen(), prob_.sc_->refGapExtend());
					TAlScore bonus = prob_.sc_->refGapOpen() - prob_.sc_->refGapExtend();
					assert_geq(bonus, 0);
					if(cpsc + bonus + sc - rfgapPen < prob_.targ_) {
						prune = true;
					}
				}
				if(prune) {
					if(extend) { nrfexPrune_++; } else { nrfopPrune_++; }
				} else if(sc - rfgapPen >= scoreFloor && sc - rfgapPen + bonusLeft >= prob_.targ_) {
					// Yes, we can introduce a ref gap here
					Edit e((int)row, '-', "ACGTN"[(int)prob_.qry_[row]], EDIT_TYPE_REF_GAP);
					assert(e.isRefGap());
					examineBranch(row - 1, col, e, sc - rfgapPen, bid);
					if(extend) { nrfex_++; } else { nrfop_++; }
				}
			}
		}
		// If we're at the top of the branch but not yet at the top of
		// the DP table, a mismatch branch is also possible.
		if(i == b.len_ && !b.curtailed_ && row >= 0 && col >= 0) {
			int rfm = prob_.ref_[col];
			assert_lt(row, (int64_t)prob_.qrylen_);
			int rdc = prob_.qry_[row];
			int rdq = prob_.qual_[row];
			int scdiff = prob_.sc_->score(rdc, rfm, rdq - 33);
			assert_lt(scdiff, 0); // at end of branch, so can't match
			bool prune = false;
			if(cp && row > 0 && col > 0 && prob_.cper_->isCheckpointed(row - 1, col - 1)) {
				// Possibly prune
				int16_t cpsc = prob_.cper_->hScore(row - 1, col - 1);
				assert_leq(cpsc, perfectScore);
				assert_leq(cpsc + scdiff + sc, perfectScore);
				if(cpsc + scdiff + sc < prob_.targ_) {
					prune = true;
				}
			}
			if(prune) {
				nmm_++;
			} else  {
				// Yes, we can introduce a mismatch here
				if(sc + scdiff >= scoreFloor && sc + scdiff + bonusLeft >= prob_.targ_) {
					Edit e((int)row, mask2dna[rfm], "ACGTN"[rdc], EDIT_TYPE_MM);
					bool nmm = (mask2dna[rfm] == 'N' || rdc > 4);
					assert_neq(e.chr, e.qchr);
					examineBranch(row - 1, col - 1, e, sc + scdiff, bid);
					if(nmm) { nnmm_++; } else { nmm_++; }
				}
			}
		}
		sc += match;
	}
}

/**
 * Sort unsorted branches, merge them with master sorted list.
 */
void BtBranchTracer::flushUnsorted() {
	if(unsorted_.empty()) {
		return;
	}
	unsorted_.sort();
	unsorted_.reverse();
#ifndef NDEBUG
	for(size_t i = 1; i < unsorted_.size(); i++) {
		assert_leq(bs_[unsorted_[i].second].score_st_, bs_[unsorted_[i-1].second].score_st_);
	}
#endif
	EList<size_t> *src2 = sortedSel_ ? &sorted1_ : &sorted2_;
	EList<size_t> *dest = sortedSel_ ? &sorted2_ : &sorted1_;
	// Merge src1 and src2 into dest
	dest->clear();
	size_t cur1 = 0, cur2 = cur_;
	while(cur1 < unsorted_.size() || cur2 < src2->size()) {
		// Take from 1 or 2 next?
		bool take1 = true;
		if(cur1 == unsorted_.size()) {
			take1 = false;
		} else if(cur2 == src2->size()) {
			take1 = true;
		} else {
			assert_neq(unsorted_[cur1].second, (*src2)[cur2]);
			take1 = bs_[unsorted_[cur1].second] < bs_[(*src2)[cur2]];
		}
		if(take1) {
			dest->push_back(unsorted_[cur1++].second); // Take from list 1
		} else {
			dest->push_back((*src2)[cur2++]); // Take from list 2
		}
	}
	assert_eq(cur1, unsorted_.size());
	assert_eq(cur2, src2->size());
	sortedSel_ = !sortedSel_;
	cur_ = 0;
	unsorted_.clear();
}

/**
 * Given the id of a branch that completes a successful backtrace, turn the
 * chain of branches into 
 */
int BtBranchTracer::trySolution(
	size_t id,
	SwResult& res,
	size_t& off,
	size_t& nrej,
	RandomSource& rnd)
{
	AlnScore score;
	BtBranch *br = &bs_[id];
	// 'br' corresponds to the leftmost edit in a right-to-left
	// chain of edits.  
	EList<Edit>& ned = res.alres.ned();
	const BtBranch *cur = br;
	size_t ns = 0, nrefns = 0;
	size_t ngap = 0;
	while(cur->e_.inited()) {
		if(cur->e_.isMismatch()) {
			if(cur->e_.qchr == 'N' || cur->e_.chr == 'N') {
				if(cur->e_.chr == 'N') {
					nrefns++;
				}
				ns++;
			}
		} else if(cur->e_.isGap()) {
			ngap++;
		}
		ned.push_back(cur->e_);
		cur = &bs_[cur->parentId_];
	}
	if(ns > prob_.nceil_) {
		// Alignment has too many Ns in it!
		res.reset();
		nrej++;
		return BT_REJECTED_N;
	}
	// Update 'seenPaths_'
	cur = br;
	bool rejSeen = false; // set =true if we overlap prev path
	bool rejCore = true; // set =true if we don't touch core diag
	while(true) {
		// Consider row, col, len, then do something
		int64_t row = cur->row_, col = cur->col_;
		assert_lt(row, (int64_t)prob_.qrylen_);
		size_t fromend = prob_.qrylen_ - row - 1;
		size_t diag = fromend + col;
		// Calculate the diagonal within the *trimmed* rectangle,
		// i.e. the rectangle we dealt with in align, gather and
		// backtrack.
		int64_t diagi = col - row;
		// Now adjust to the diagonal within the *untrimmed*
		// rectangle by adding on the amount trimmed from the left.
		diagi += prob_.rect_->triml;
		assert_lt(diag, seenPaths_.size());
		int64_t newlo = (int64_t)cur->len_ > row ? 0 : (row - cur->len_);
		int64_t newhi = row+1;
		assert_geq(newlo, 0);
		assert_geq(newhi, 0);
		assert((int64_t)cur->len_ > row || newhi > newlo);
		if(newhi > newlo) {
			bool added = false;
			// Does it overlap a core diagonal?
			if(diagi >= 0) {
				size_t diag = (size_t)diagi;
				if(diag >= prob_.rect_->corel &&
				   diag <= prob_.rect_->corer)
				{
					rejCore = false;
				}
			}
			const size_t sz = seenPaths_[diag].size();
			for(size_t i = 0; i < sz; i++) {
				// Does the new interval overlap this already-seen
				// interval?  Also of interest: does it abut this
				// already-seen interval?  If so, we should merge them.
				size_t lo = seenPaths_[diag][i].first;
				size_t hi = seenPaths_[diag][i].second;
				assert_lt(lo, hi);
				size_t lo_sm = newlo, hi_sm = newhi;
				if(hi - lo < hi_sm - lo_sm) {
					swap(lo, lo_sm);
					swap(hi, hi_sm);
				}
				if((lo <= lo_sm && hi > lo_sm) ||
				   (lo <  hi_sm && hi >= hi_sm))
				{
					// One or both of the shorter interval's end points
					// are contained in the longer interval - so they
					// overlap.
					rejSeen = true;
					// Merge them into one longer interval
					seenPaths_[diag][i].first = min(lo, lo_sm);
					seenPaths_[diag][i].second = max(hi, hi_sm);
#ifndef NDEBUG
					for(int64_t ii = seenPaths_[diag][i].first;
					    ii < seenPaths_[diag][i].second;
						ii++)
					{
						cerr << "trySolution rejected (" << ii << ", " << (ii + col - row) << ")" << endl;
					}
#endif
					added = true;
					break;
				} else if(hi == lo_sm || lo == hi_sm) {
					// Merge them into one longer interval
					seenPaths_[diag][i].first = min(lo, lo_sm);
					seenPaths_[diag][i].second = max(hi, hi_sm);
#ifndef NDEBUG
					for(int64_t ii = seenPaths_[diag][i].first;
					    ii < seenPaths_[diag][i].second;
						ii++)
					{
						cerr << "trySolution rejected (" << ii << ", " << (ii + col - row) << ")" << endl;
					}
#endif
					added = true;
					// Keep going in case it overlaps one of the other
					// intervals
				}
			}
			if(!added) {
				seenPaths_[diag].push_back(make_pair(newlo, newhi));
			}
		}
		// After the merging that may have occurred above, it's no
		// longer guarnateed that all the overlapping intervals in
		// the list have been merged.  That's OK though.  We'll
		// still get correct answers to overlap queries.
		if(!cur->e_.inited()) {
			assert_eq(0, cur->parentId_);
			break;
		}
		cur = &bs_[cur->parentId_];
	} // while(cur->e_.inited())
	if(rejSeen) {
		return BT_NOT_FOUND;
	}
	if(rejCore) {
		return BT_REJECTED_CORE_DIAG;
	}
	off = br->leftmostCol();
	score.score_ = prob_.targ_;
	score.ns_    = ns;
	score.gaps_  = ngap;
	res.alres.setScore(score);
	res.alres.setRefNs(nrefns);
	size_t trimBeg = br->uppermostRow();
	size_t trimEnd = prob_.qrylen_ - prob_.row_ - 1;
	assert_leq(trimBeg, prob_.qrylen_);
	assert_leq(trimEnd, prob_.qrylen_);
	TRefOff refoff = off + prob_.refoff_ + prob_.rect_->refl;
	res.alres.setShape(
		prob_.refid_,                   // ref id
		refoff,                         // 0-based ref offset
		prob_.fw_,                      // aligned to Watson?
		prob_.qrylen_,                  // read length
		true,                           // pretrim soft?
		0,                              // pretrim 5' end
		0,                              // pretrim 3' end
		true,                           // alignment trim soft?
		prob_.fw_ ? trimBeg : trimEnd,  // alignment trim 5' end
		prob_.fw_ ? trimEnd : trimBeg); // alignment trim 3' end
	return BT_FOUND;
}

/**
 * Get the next valid alignment given a backtrace problem.  Return false
 * if there is no valid solution.  Use a backtracking search to find the
 * solution.  This can be very slow.
 */
bool BtBranchTracer::nextAlignmentBacktrace(
	size_t maxiter,
	SwResult& res,
	size_t& off,
	size_t& nrej,
	size_t& niter,
	RandomSource& rnd)
{
	assert(!empty() || !emptySolution());
	assert(prob_.inited());
	//ASSERT_ONLY(TAlScore lastScore = std::numeric_limits<TAlScore>::max());
	// There's a subtle case where we might fail to backtracing in
	// local-alignment mode.  The basic fact to remember is that when we're
	// backtracing from the highest-scoring cell in the table, we're guaranteed
	// to be able to backtrace without ever dipping below 0.  But if we're
	// backtracing from a cell other than the highest-scoring cell in the
	// table, we might dip below 0.  Dipping below 0 implies that there's a
	// shorted local alignment with a better score.  In which case, it's
	// perfectly fair for us to abandon any path that dips below the floor, and
	// this might result in the queue becoming empty before we finish.
	bool result = false;
	niter = 0;
	while(!empty()) {
		if(trySolutions(res, off, nrej, rnd, result)) {
			return result;
		}
		if(niter++ >= maxiter) {
			break;
		}
		size_t brid = best(rnd); // put best branch in 'br'
		assert(!seen_.contains(brid));
		ASSERT_ONLY(seen_.insert(brid));
		BtBranch *br = &bs_[brid];
		cerr << brid
		     << ": targ:" << prob_.targ_
		     << ", sc:" << br->score_st_
		     << ", row:" << br->uppermostRow()
			 << ", nmm:" << nmm_
			 << ", nnmm:" << nnmm_
			 << ", nrdop:" << nrdop_
			 << ", nrfop:" << nrfop_
			 << ", nrdex:" << nrdex_
			 << ", nrfex:" << nrfex_
			 << ", nrdop_pr: " << nrdopPrune_
			 << ", nrfop_pr: " << nrfopPrune_
			 << ", nrdex_pr: " << nrdexPrune_
			 << ", nrfex_pr: " << nrfexPrune_
			 << endl;
		addOffshoots(brid);
	}
	if(trySolutions(res, off, nrej, rnd, result)) {
		return result;
	}
	return false;
}

/**
 * Get the next valid alignment given a backtrace problem.  Return false
 * if there is no valid solution.  Use a triangle-fill backtrace to find
 * the solution.  This is usually fast (it's O(m + n)).
 */
bool BtBranchTracer::nextAlignmentFill(
	size_t maxiter,
	SwResult& res,
	size_t& off,
	size_t& nrej,
	size_t& niter,
	RandomSource& rnd)
{
	assert(prob_.cper_->doCheckpoints());
	assert(prob_.cper_->hasEF());
	return false;
}

/**
 * Get the next valid alignment given the backtrace problem.  Return false
 * if there is no valid solution, e.g., if 
 */
bool BtBranchTracer::nextAlignment(
	size_t maxiter,
	SwResult& res,
	size_t& off,
	size_t& nrej,
	size_t& niter,
	RandomSource& rnd)
{
	if(prob_.fill_) {
		return nextAlignmentFill(
			maxiter,
			res,
			off,
			nrej,
			niter,
			rnd);
	} else {
		return nextAlignmentBacktrace(
			maxiter,
			res,
			off,
			nrej,
			niter,
			rnd);
	}
}

#ifdef MAIN_ALIGNER_BT

#include <iostream>

int main(int argc, char **argv) {
	size_t off = 0;
	RandomSource rnd(77);
	BtBranchTracer tr;
	Scoring sc = Scoring::base1();
	SwResult res;
	tr.init(
		"ACGTACGT", // in: read sequence
		"IIIIIIII", // in: quality sequence
		8,          // in: read sequence length
		"ACGTACGT", // in: reference sequence
		8,          // in: reference sequence length
		0,          // in: reference id
		0,          // in: reference offset
		true,       // in: orientation
		sc,         // in: scoring scheme
		0,          // in: N ceiling
		8,          // in: alignment score
		7,          // start in this row
		7,          // start in this column
		rnd);       // random gen, to choose among equal paths
	size_t nrej = 0;
	tr.nextAlignment(
		res,
		off,
		nrej,
		rnd);
}

#endif /*def MAIN_ALIGNER_BT*/
