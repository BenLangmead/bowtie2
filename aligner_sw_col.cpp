/*
 *  aligner_sw_col.cpp
 */

#include "aligner_sw.h"
#include "search_globals.h"
#include "scoring.h"
#include "mask.h"

/**
 * Select a path for backtracking from this cell.  If there is a
 * tie among eligible paths, break it randomly.  Return value is
 * a pair where first = a flag indicating the backtrack type (see
 * enum defining SW_BT_* above), and second = a selection for
 * the read character for the next row up.  second should be
 * ignored if the backtrack type is a gap in the read.
 */
pair<int, int> SwColorCellMask::randBacktrack(
	RandomSource& rand,
	bool& branch)
{
	ASSERT_ONLY(int num = numPossible());
	std::pair<int, int> ret;
	ret.second = -1;
	int i =
		((diag != 0) << 0) |
		((rfop != 0) << 1) |
		((rfex != 0) << 2) |
		((rdop != 0) << 3) |
		((rdex != 0) << 4);
	// Count the total number of different choices we could make
	int totChoices =
		alts5[diag] +
		alts5[rfop] +
		alts5[rfex] +
		alts5[rdop] +
		alts5[rdex];
	// If we're choosing from among >1 possibilities, inform caller so that
	// caller can add a frame to the backtrack stack.
	branch = totChoices > 1;
	ret.first = randFromMask(rand, i);
	assert_lt(ret.first, 5);
	assert(ret.first == SW_BT_DIAG ||
		   ret.first == SW_BT_REF_OPEN ||
		   ret.first == SW_BT_REF_EXTEND ||
		   ret.first == SW_BT_READ_OPEN ||
		   ret.first == SW_BT_READ_EXTEND);
	// Clear the bit associated with the path chosen
	if(ret.first < 3) {
		// Must choose character for next row
		if(ret.first == SW_BT_DIAG) {
			assert(diag != 0);
			ret.second = randFromMask(rand, diag);
			diag &= ~(1 << ret.second);
		} else if(ret.first == SW_BT_REF_OPEN) {
			assert(rfop != 0);
			ret.second = randFromMask(rand, rfop);
			rfop &= ~(1 << ret.second);
		} else if(ret.first == SW_BT_REF_EXTEND) {
			assert(rfex != 0);
			ret.second = randFromMask(rand, rfex);
			rfex &= ~(1 << ret.second);
		}
	} else {
		if(ret.first == SW_BT_READ_OPEN) {
			rdop = 0;
		}
		else if(ret.first == SW_BT_READ_EXTEND) {
			rdex = 0;
		}
	}
	assert_eq(num-1, numPossible());
	return ret;
}

/**
 * We finished updating the cell; set empty and finalized
 * appropriately.
 */
inline bool SwColorCell::finalize(TAlScore floorsc) {
	ASSERT_ONLY(finalized = true);
	for(int i = 0; i < 4; i++) {
		if(!mask[i].empty() && best[i].score() >= floorsc) {
			assert(VALID_AL_SCORE(best[i]));
			assert_geq(best[i].score(), floorsc);
			empty = false;
#ifdef NDEBUG
			break;
#endif
		}
	}
	return !empty;
}

/**
 * Given the dynamic programming table and a cell (both the table offset
 * and the reference character), trace backwards from the cell and install
 * the edits and score/penalty in the appropriate fields of res.  The
 * RandomSource is used to break ties among equally good ways of tracing
 * back.
 *
 * Note that the subject nucleotide sequence is decoded at the same time as
 * the alignment is constructed.  So the traceback reveals both the
 * nucleotide decoding (in ned) and the colorspace error pattern (in ced).
 * The approach is very similar to the one described in the SHRiMP
 * paper:
 *
 * Rumble SM, Lacroute P, Dalca AV, Fiume M, Sidow A, Brudno M. SHRiMP:
 * accurate mapping of short color-space reads. PLoS Comput Biol. 2009
 * May;5(5)
 *
 * Whenever we enter a cell, we check whether the read/ref coordinates of
 * that cell correspond to a cell we traversed constructing a previous
 * alignment.  If so, we backtrack to the last decision point, mask out the
 * path that led to the previously observed cell, and continue along a
 * different path; or, if there are no more paths to try, we give up.
 *
 * If an alignment is found, 'off' is set to the alignment's upstream-most
 * reference character's offset into the chromosome and true is returned.
 * Otherwise, false is returned.
 */
bool SwAligner::backtrackColors(
	TAlScore      escore,  // score we expect to get over backtrack
	SwResult&     res,     // out: store results (edits and scores) here
	size_t&       off,     // out: leftmost ref char involved in aln
	size_t        row,     // start in this row (w/r/t full matrix)
	size_t        col,     // start in this column (w/r/t full matrix)
	int           lastC,   // cell to backtrace from in lower-right corner
	RandomSource& rand)    // pseudo-random generator
{
	ELList<SwColorCell>& tab = ctab_;
	assert_lt(row, tab.size());
	assert_geq(escore, minsc_);
	btcstack_.clear();
	btcells_.clear();
	size_t tabcol = col - row;
	ASSERT_ONLY(BTDnaString drd);
	ASSERT_ONLY(drd.resize(rdf_-rdi_+1));
	int curC = lastC;
	AlnScore score; score.score_ = 0;
	score.gaps_ = score.ns_ = 0;
	bool refExtend = false, readExtend = false;
	ASSERT_ONLY(size_t origCol = col);
	int gaps = 0;
	res.alres.reset();
	EList<Edit>& ned = res.alres.ned();
	//EList<Edit>& aed = res.alres.aed();
	EList<Edit>& ced = res.alres.ced();
	res.ndn = lastC;
	assert_geq(score.score(), floorsc_);
	assert(!sc_->monotone || escore <= 0);
	assert(!sc_->monotone || score.score() >= escore);
	while((int)row >= 0) {
		res.swbts++;
		assert_lt(row, tab.size());
		assert_leq(col, origCol);
		assert_range(0, 3, curC);
		ASSERT_ONLY(drd.set(curC, row));
		assert_geq(col, row);
		tabcol = col - row;
		assert_geq(tabcol, 0);
		// Cell was involved in a previously-reported alignment?
		if(tab[row][tabcol].reportedThru_) {
			if(!btcstack_.empty()) {
				// Pop record off the top of the stack
				ned.resize(btcstack_.back().nedsz);
				btcells_.resize(btnstack_.back().celsz);
				//aed.resize(btcstack_.back().aedsz);
				row   = btcstack_.back().row;
				col   = btcstack_.back().col;
				gaps  = btcstack_.back().gaps;
				score = btcstack_.back().score;
				curC  = btcstack_.back().curC;
				btcstack_.pop_back();
				assert(!sc_->monotone || score.score() >= escore);
				continue;
			} else {
				// No branch points to revisit; just give up
				return false;
			}
		}
		assert(!tab[row][tabcol].reportedThru_);
		assert(!sc_->monotone || score.score() >= escore);
		if(row == 0) {
			btcells_.expand();
			btcells_.back().first = row;
			btcells_.back().second = tabcol;
			break;
		}
		assert_gt(tab[row][tabcol].mask[curC].numPossible(), 0);
		bool branch = false;
		pair<int, int> cur =
			tab[row][tabcol].mask[curC].randBacktrack(rand, branch);
		if(branch) {
			assert_gt(tab[row][tabcol].mask[curC].numPossible(), 0);
			// Add a frame to the backtrack stack
			btcstack_.expand();
			btcstack_.back().init(
				ned.size(),
				0,          // aed.size()
				ced.size(),
				btcells_.size(),
				row,
				col,
				curC,
				gaps,
				score);
		}
		btcells_.expand();
		btcells_.back().first = row;
		btcells_.back().second = tabcol;
		switch(cur.first) {
			// Move up and to the left.  If the reference nucleotide in the
			// source row mismatches the decoded nucleotide, penalize
			// it and add a nucleotide mismatch.  If the color being
			// traversed is a miscall, penalize that.
			case SW_BT_DIAG: {
				assert_neq(-1, cur.second);
				refExtend = readExtend = false;
				assert_gt(row, 0); assert_gt(col, 0);
				// Check for base mismatch at source (lower-right) cell
				int refNmask = (int)rf_[rfi_+col];
				assert_gt(refNmask, 0);
				int m = matchesEx(curC, refNmask);
				if(m != 1) {
					Edit e((int)row, mask2dna[refNmask], "ACGTN"[curC], EDIT_TYPE_MM);
					assert(e.repOk());
					ned.push_back(e);
					score.score_ -= ((m == 0) ? sc_->snp : sc_->n(30));
					assert_geq(score.score(), floorsc_);
					assert(!sc_->monotone || score.score() >= escore);
				} else {
					score.score_ += sc_->match(30);
					assert_geq(score.score(), floorsc_);
					assert(!sc_->monotone || score.score() >= escore);
				}
				if(m == -1) {
					score.ns_++;
				}
				row--; col--;
				assert_lt(row, tab.size());
				assert(VALID_AL_SCORE(score));
				// Check for color mismatch
				int readC = (*rd_)[rdi_ + row];
				int decC = dinuc2color[cur.second][curC];
				assert_range(0, 3, decC);
				if(decC != readC) {
					Edit e((int)row, "ACGT"[decC], "ACGTN"[readC], EDIT_TYPE_MM);
					score.score_ -= QUAL(row); // color mismatch
					assert_geq(score.score(), floorsc_);
					assert(!sc_->monotone || score.score() >= escore);
					assert(e.repOk());
					ced.push_back(e);
				}
				if(readC > 3) {
					score.ns_++;
				}
				assert_lt(tabcol, tab[row].size());
				assert(VALID_AL_SCORE(score));
				lastC = curC;
				curC = cur.second;
				break;
			}
			// Move up.  Add an edit encoding the ref gap.  If the color
			// being traversed is a miscall, penalize that.
			case SW_BT_REF_OPEN: {
				refExtend = true; readExtend = false;
				assert_gt(row, 0);
				Edit e((int)row, '-', "ACGTN"[curC], EDIT_TYPE_DEL);
				assert(e.repOk());
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				row--;
				// Check for color miscall
				int decC = dinuc2color[cur.second][curC];
				assert_range(0, 3, decC);
				int rdC = (*rd_)[rdi_+row];
				if(decC != rdC) {
					score.score_ -= QUAL(row);
					assert_geq(score.score(), floorsc_);
					assert(!sc_->monotone || score.score() >= escore);
					Edit e((int)row, "ACGT"[decC], "ACGTN"[rdC], EDIT_TYPE_MM);
					assert(e.repOk());
					ced.push_back(e);
				}
				if(rdC > 3) {
					score.ns_++;
				}
				score.score_ -= sc_->refGapOpen();
				score.gaps_++;
				assert_geq(score.score(), floorsc_);
				assert(!sc_->monotone || score.score() >= escore);
				gaps++;
#ifndef NDEBUG
				assert_leq(score.gaps_, rdgap_ + rfgap_);
				assert_lt(row, tab.size());
				tabcol = col - row;
				assert_lt(tabcol, tab[row].size());
#endif
				lastC = curC;
				curC = cur.second;
				break;
			}
			// Move up.  Add an edit encoding the ref gap.  If the color
			// being traversed is a miscall, penalize that.
			case SW_BT_REF_EXTEND: {
				refExtend = true; readExtend = false;
				assert_gt(row, 1);
				Edit e((int)row, '-', "ACGTN"[curC], EDIT_TYPE_DEL);
				assert(e.repOk());
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				row--;
				// Check for color mismatch
				int decC = dinuc2color[cur.second][curC];
				assert_range(0, 3, decC);
				int rdC = (*rd_)[rdi_+row];
				if(decC != rdC) {
					score.score_ -= QUAL(row);
					assert_geq(score.score(), floorsc_);
					assert(!sc_->monotone || score.score() >= escore);
					Edit e((int)row, "ACGT"[decC], "ACGTN"[rdC], EDIT_TYPE_MM);
					assert(e.repOk());
					ced.push_back(e);
				}
				if(rdC > 3) {
					score.ns_++;
				}
				score.score_ -= sc_->refGapExtend();
				assert_geq(score.score(), floorsc_);
				assert(!sc_->monotone || score.score() >= escore);
				score.gaps_++;
				gaps++;
#ifndef NDEBUG
				assert_leq(score.gaps_, rdgap_ + rfgap_);
				assert_lt(row, tab.size());
				tabcol = col - row;
				assert_lt(tabcol, tab[row].size());
#endif
				lastC = curC;
				curC = cur.second;
				break;
			}
			case SW_BT_READ_OPEN: {
				refExtend = false; readExtend = true;
				assert_gt(col, 0);
				Edit e((int)row+1, "ACGTN"[firsts5[(int)rf_[rfi_+col]]], '-', EDIT_TYPE_INS);
				assert(e.repOk());
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				col--;
				score.score_ -= sc_->readGapOpen();
				score.gaps_++;
				assert_geq(score.score(), floorsc_);
				assert(!sc_->monotone || score.score() >= escore);
				gaps++;
#ifndef NDEBUG
				assert_leq(score.gaps_, rdgap_ + rfgap_);
				assert_lt(row, tab.size());
				tabcol = col - row;
				assert_lt(tabcol, tab[row].size());
#endif
				lastC = curC;
				break;
			}
			case SW_BT_READ_EXTEND: {
				refExtend = false; readExtend = true;
				assert_gt(col, 1);
				Edit e((int)row+1, "ACGTN"[firsts5[(int)rf_[rfi_+col]]], '-', EDIT_TYPE_INS);
				assert(e.repOk());
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				col--;
				score.score_ -= sc_->readGapExtend();
				score.gaps_++;
				assert_geq(score.score(), floorsc_);
				assert(!sc_->monotone || score.score() >= escore);
				gaps++;
#ifndef NDEBUG
				assert_leq(score.gaps_, rdgap_ + rfgap_);
				assert_lt(row, tab.size());
				tabcol = col - row;
				assert_lt(tabcol, tab[row].size());
#endif
				lastC = curC;
				break;
			}
			default: throw 1;
		}
	} // while((int)row > 0)
	for(size_t i = 0; i < btcells_.size(); i++) {
		size_t rw = btcells_[i].first;
		size_t cl = btcells_[i].second;
		assert(!tab[rw][cl].reportedThru_);
		tab[rw][cl].setReportedThrough();
	}
	assert_eq(0, row);
	int refNmask = (int)rf_[rfi_+col]; // get last char in ref involved in alignment
	assert_gt(refNmask, 0);
	int m = matchesEx(curC, refNmask);
	if(m != 1) {
		Edit e((int)row, mask2dna[refNmask], "ACGTN"[curC], EDIT_TYPE_MM);
		assert(e.repOk());
		ned.push_back(e);
		score.score_ -= ((m == 0) ? sc_->snp : sc_->n(30));
		assert_geq(score.score(), floorsc_);
		assert(!sc_->monotone || score.score() >= escore);
	}
	if(m == -1) {
		score.ns_++;
	}
	ASSERT_ONLY(drd.set(curC, 0));
	res.nup = curC;
	res.reverse();
	assert(Edit::repOk(ced, (*rd_)));
#ifndef NDEBUG
	BTDnaString refstr;
	for(size_t i = col; i <= origCol; i++) {
		refstr.append(firsts5[(int)rf_[rfi_+i]]);
	}
	BTDnaString editstr;
	Edit::toRef(drd, ned, editstr);
	if(refstr != editstr) {
		cerr << "Decoded nucleotides and edits don't match reference:" << endl;
		cerr << "           score: " << score.score() << " (" << gaps << " gaps)" << endl;
		cerr << "           edits: ";
		Edit::print(cerr, ned);
		cerr << endl;
		cerr << "    decoded nucs: " << drd << endl;
		cerr << "     edited nucs: " << editstr << endl;
		cerr << "  reference nucs: " << refstr << endl;
		assert(0);
	}
#endif
	// done
	assert_eq(score.score(), escore);
	assert_leq(gaps, rdgap_ + rfgap_);
	// Dummy values for refid and fw
	res.alres.setScore(score);
	off = (col - row);
	return true;
}

/**
 * Update a SwCell's best[] and mask[] arrays with respect to its
 * neighbor on the left.
 */
inline void SwAligner::updateColorHoriz(
	const SwColorCell& lc,
	SwColorCell&       dstc,
	int                rfm)
{
	assert(lc.finalized);
	if(lc.empty) return;
	bool ninvolved = rfm > 15;
	for(int to = 0; to < 4; to++) {
		if(lc.empty) {
			// There's no way to get to the cell to our left within
			// budget if the read character assigned to this row is
			// 'to', so try next character assignment.
			continue;
		}
		AlnScore leftBest = lc.best[to];
		assert(!sc_->monotone || leftBest.score() <= 0);
		const SwColorCellMask& fromMask = lc.mask[to];
		AlnScore& myBest = dstc.best[to];
		SwColorCellMask& myMask = dstc.mask[to];
		if(ninvolved) leftBest.incNs(nceil_);
		if(!VALID_AL_SCORE(leftBest)) continue;
		// *Don't* penalize for a nucleotide mismatch because we must
		// have already done that in a previous vertical or diagonal
		// step.
		if(fromMask.readExtendPossible()) {
			// Read gap extension possible?
			AlnScore ex = leftBest;
			assert_leq(ex.score(), 0);
			assert(VALID_AL_SCORE(ex));
			ex.score_ -= sc_->readGapExtend();
			assert_leq(ex.score(), 0);
			if(ex.score_ >= floorsc_ && ex >= myBest) {
				if(ex > myBest) {
					myMask.clear();
					myBest = ex;
				}
				myMask.rdex = 1;
				assert(VALID_AL_SCORE(myBest));
			}
		}
		if(fromMask.readOpenPossible()) {
			// Read gap open possible?
			AlnScore ex = leftBest;
			assert(VALID_AL_SCORE(ex));
			ex.score_ -= sc_->readGapOpen();
			if(ex.score_ >= floorsc_ && ex >= myBest) {
				if(ex > myBest) {
					myMask.clear();
					myBest = ex;
				}
				myMask.rdop = 1;
				assert(VALID_AL_SCORE(myBest));
			}
		}
	}
}

/**
 * Update a SwCell's best[] and mask[] arrays with respect to its
 * neighbor on the left.
 */
inline void SwAligner::updateColorVert(
	const SwColorCell& uc,        // cell above
	SwColorCell&       dstc,      // destination cell
	int                c,         // color b/t this row, one above
	int                penmm)     // penalty to incur for color miscall
{
	assert(uc.finalized);
	if(uc.empty) return;
	for(int to = 0; to < 4; to++) {
		// Assuming that the read character in this row is 'to'...
		// No SNP penalty because destination read char aligns to a
		// gap in the reference.
		AlnScore from[] = {
			uc.best[0], uc.best[1],
			uc.best[2], uc.best[3] };
		// Calculate which 'from' character isn't going to
		// cause a color mismatch
		const int goodfrom = nuccol2nuc[to][c];
		// Reward the 'from' that corroborates the color
		if(goodfrom < 4) from[goodfrom].score_ += penmm;
		for(int fr = 0; fr < 4; fr++) {
			AlnScore frBest = from[fr];
			const SwColorCellMask& frMask = uc.mask[fr];
			AlnScore& myBest = dstc.best[to];
			SwColorCellMask& myMask = dstc.mask[to];
			if(c > 3) {
				frBest.incNs(nceil_);
			}
			if(!VALID_AL_SCORE(frBest)) continue;
			frBest.score_ -= penmm;
			assert(!sc_->monotone || frBest.score() <= 0);
			if(frMask.refExtendPossible()) {
				// Extend is possible
				frBest.score_ -= sc_->refGapExtend();
				assert_leq(frBest.score_, 0);
				if(frBest.score_ >= floorsc_ && frBest >= myBest) {
					if(frBest > myBest) {
						myBest = frBest;
						myMask.clear();
					}
					myMask.rfex |= (1 << fr);
					assert(VALID_AL_SCORE(myBest));
				}
				// put it back
				frBest.score_ += sc_->refGapExtend();
			}
			if(frMask.refOpenPossible()){
				// Open is possible
				frBest.score_ -= sc_->refGapOpen();
				if(frBest.score_ >= floorsc_ && frBest >= myBest) {
					if(frBest > myBest) {
						myBest = frBest;
						myMask.clear();
					}
					myMask.rfop |= (1 << fr);
					assert(VALID_AL_SCORE(myBest));
				}
			}
		}
	}
}

/**
 * Update a SwCell's best[] and mask[] arrays with respect to its
 * neighbor up and to the left.  SNPs are charged
 */
inline void SwAligner::updateColorDiag(
	const SwColorCell& uc,        // cell above and to the left
	SwColorCell&       dstc,      // destination cell
	int                refMask,   // ref mask associated with destination cell
	int                c,         // color being traversed
	int                penmm)     // penalty to incur for color miscall
{
	assert(uc.finalized);
	if(uc.empty) return;
	bool ninvolved = (c > 3 || refMask > 15);
	for(int to = 0; to < 4; to++) {
		// Assuming that the read character in this row is 'to'...
		int add = (matches(to, refMask) ? sc_->match(30) : -((refMask > 15) ? sc_->n(30) : sc_->snp));
		AlnScore from[] = {
			uc.best[0] + add, uc.best[1] + add,
			uc.best[2] + add, uc.best[3] + add };
		// Calculate which 'from' character isn't going to
		// cause a color mismatch
		const int goodfrom = nuccol2nuc[to][c];
		if(goodfrom < 4) {
			from[goodfrom].score_ += penmm;
		}
		for(int fr = 0; fr < 4; fr++) {
			AlnScore frBest = from[fr];
			AlnScore& myBest = dstc.best[to];
			SwColorCellMask& myMask = dstc.mask[to];
			if(ninvolved) frBest.incNs(nceil_);
			if(!VALID_AL_SCORE(frBest)) continue;
			frBest.score_ -= penmm;
			assert(!sc_->monotone || frBest.score() <= 0);
			if(frBest.score_ >= floorsc_ && frBest >= myBest) {
				if(frBest > myBest) {
					myBest = frBest;
					myMask.clear();
					assert_eq(0, myMask.diag);
				}
				myMask.diag |= (1 << fr);
				assert(VALID_AL_SCORE(myBest));
			}
		}
	}
}

// Assumes:
// 1. 'row' = 0-based row of DP table
// 2. 'col' = 0-based col of DP table
// 3. 'minsc_' = minimum AlignmentScore required for a cell to have a sol
// 4. 'solrows_' is a pairs encoding the range of rows s.t. all rows with
//    solutions are in the range
// 5. 'solcols_' is a list of pairs, where each pair cor
// 6. 'solbest_'
// 7. 'solrowbest_'
#define UPDATE_SOLS(cur, row, col) { \
	if(en_ == NULL || (*en_)[col]) { \
		for(int I = 0; I < 4; I++) { \
			if(cur.best[I].score() >= minsc_) { \
				if((int64_t)row >= solrowlo_) { \
					/* Score is acceptable */ \
					if(row < solrows_.first)           solrows_.first = row; \
					if(row > solrows_.second)          solrows_.second = row; \
					if(col < solcols_[row].first)      solcols_[row].first  = col; \
					if(col > solcols_[row].second)     solcols_[row].second = col; \
					if(cur.best[I].score() > solbest_) \
						solbest_ = cur.best[I].score(); \
					assert(!sc_->monotone || solbest_ <= 0); \
					if(cur.best[I].score() > solrowbest_[row]) \
						solrowbest_[row] = cur.best[I].score(); \
					assert(!sc_->monotone || solrowbest_[row] <= 0); \
					nsols_++; \
					soldone_ = false; \
					assert(repOk()); \
				} \
			} \
		} \
	} \
}

#define FINALIZE_CELL(r, c) { \
	assert(!tab[r][c].finalized); \
	if(tab[r][c].finalize(floorsc_)) { \
		UPDATE_SOLS(tab[r][c], r, c); \
		validInRow = true; \
	} \
	assert(tab[r][c].finalized); \
}

/**
 * Align the colorspace read 'read' as aligned against the reference
 * string 'rf' while simultaneously decoding the nucleotide sequence.
 * The approach is very similar to what's described in the SHRiMP
 * paper:
 *
 * Rumble SM, Lacroute P, Dalca AV, Fiume M, Sidow A, Brudno M. SHRiMP:
 * accurate mapping of short color-space reads. PLoS Comput Biol. 2009
 * May;5(5)
 *
 * If an alignment is found, its offset relative to rdi_ is returned.
 * E.g. if an alignment is found that occurs starting at rdi, 0 is
 * returned.  If no alignment is found, -1 is returned.
 */
bool SwAligner::alignColors(RandomSource& rnd) {
	typedef SwColorCell TCell;
	assert_leq(rdf_, rd_->length());
	assert_leq(rdf_, qu_->length());
	assert_lt(rfi_, rff_);
	assert_lt(rdi_, rdf_);
	assert_eq(rd_->length(), qu_->length());
	assert_geq(sc_->gapbar, 1);
	assert(repOk());
#ifndef NDEBUG
	for(size_t i = rfi_; i < rff_; i++) {
		assert_range(0, 16, (int)rf_[i]);
	}
#endif
	//
	// Initialize the first row
	//
	ELList<TCell>& tab = ctab_;
	tab.resize(1); // add first row to row list
	solrows_ = SwAligner::EXTREMES;        // at first, no rows have sols
	solcols_.resize(rdf_ - rdi_ + 1);
	solcols_.fill(SwAligner::EXTREMES);    // at first, no cols have sols
	solrowbest_.resize(rdf_ - rdi_ + 1);
	solrowbest_.fill(std::numeric_limits<TAlScore>::min()); // no row has sol
	const size_t wlo = 0;
	const size_t whi = (int)(width_ - 1);
	tab[0].resize(whi-wlo+1); // add columns to first row
	bool validInRow = !sc_->monotone;;
	// Calculate starting values for the rest of the columns in the
	// first row.  No need to consider colors yet since we won't have
	// consecutive pairs of read nucleotides until we get to the second
	// row.  We just consider reference mutations.
	for(size_t col = 0; col <= whi; col++) {
		tab[0][col].clear(); // clear the cell; masks and scores
		int rfm = rf_[rfi_+col];
		// Can we start from here?
		bool canStart = (st_ == NULL || (*st_)[col]);
		if(canStart) {
			for(int to = 0; to < 4; to++) {
				int m = matchesEx(to, rfm);
				if(m == 1) {
					// The assigned subject nucleotide matches the reference;
					// no penalty
					assert_lt(rfm, 16);
					tab[0][col].best[to].gaps_ = 0;
					tab[0][col].best[to].ns_ = 0;
					tab[0][col].best[to].score_ = 0;
					tab[0][col].mask[to].diag = 0xf;
				} else if(m == 0 && sc_->snp <= -minsc_) {
					// The assigned subject nucleotide does not match the
					// reference nucleotide, so we add a SNP penalty
					tab[0][col].best[to].gaps_ = 0;
					tab[0][col].best[to].ns_ = 0;
					tab[0][col].best[to].score_ = -sc_->snp;
					tab[0][col].mask[to].diag = 0xf;
				} else if(m == -1) {
					// The mask or query character is an N
					tab[0][col].best[to].gaps_ = 0;
					tab[0][col].best[to].ns_ = 1;
					tab[0][col].best[to].score_ = -sc_->n(30);
					tab[0][col].mask[to].diag = 0xf;
				} else {
					// Leave mask[to] cleared
				}
			}
		}
		// Calculate horizontals if barrier allows
		if(sc_->gapbar <= 1 && col > 0) {
			updateColorHoriz(tab[0][col-1], tab[0][col], rfm);
			ncups_++;
		}
		FINALIZE_CELL(0, col);
	}
	nrowups_++;
	if(!validInRow) {
		nrowskips_ += (rdf_ - rdi_);
		return !soldone_;
	}

	//
	// Calculate all subsequent rows
	//

	// Do rest of table
	for(size_t row = 1; row <= rdf_ - rdi_; row++) {
		nrowups_++;
		tab.expand(); // add another row
		bool onlyDiagInto =
			(row+1             <= (size_t)sc_->gapbar ||
			 rdf_ - rdi_ - row <= (size_t)sc_->gapbar);
		tab.back().resize(whi-wlo+1); // add enough space for columns
		assert_gt(row, 0);
		assert_leq(row, qu_->length());
		assert_leq(row, rd_->length());
		int c = (*rd_)[row-1];   // read character in this row
		int q = QUAL(row-1); // quality for the read character; should be Phred
		assert_geq(q, 0);
		const int mmpen = ((c > 3) ? sc_->n(30) : sc_->mm(q));
		validInRow = !sc_->monotone;;		
		//
		// Handle col == wlo case before (and the col == whi case
		// after) the upcoming inner loop to reduce the number of
		// guards needed inside it.
		//
		assert_leq(wlo, whi);
		size_t col = wlo;
		TCell& cur = tab[row][0];
		cur.clear();
		if(!tab[row-1][0].empty) {
			const size_t fullcol = col + row;
			updateColorDiag(
				tab[row-1][0],     // cell diagonally above and to the left
				cur,               // destination cell
				rf_[rfi_ + fullcol], // ref mask associated with destination cell
				c,                 // color being traversed
				mmpen);            // penalty to incur for color miscall
		}
		if(!onlyDiagInto && col < whi && !tab[row-1][1].empty) {
			updateColorVert(
				tab[row-1][1],     // cell above
				cur,               // destination cell
				c,                 // color being traversed
				mmpen);            // penalty to incur for color miscall
		}
		ncups_++;
		// 'cur' is now initialized
		FINALIZE_CELL(row, col);
		// Iterate from leftmost to rightmost inner diagonals
		for(col = wlo+1; col < whi; col++) {
			const size_t fullcol = col + row;
			int r = rf_[rfi_ + fullcol];
			TCell& cur = tab[row][col-wlo];
			cur.clear();
			TCell& dg = tab[row-1][col-wlo];
			// The mismatch penalty is a function of the read character
			// in this row & the reference character in this column
			// (specifically: whether they match and whether either is
			// an N) as well as the quality value of the read
			// character.
			updateColorDiag(
				dg,            // cell diagonally above and to the left
				cur,           // destination cell
				r,             // ref mask associated with destination column
				c,             // nucleotide in destination row
				mmpen);        // penalty to incur for color miscall
			if(!onlyDiagInto) {
				TCell& up = tab[row-1][col-wlo+1];
				updateColorVert(
					up,        // cell above
					cur,       // destination cell
					c,         // nucleotide in destination row
					mmpen);    // penalty to incur for color miscall
				// Can do horizontal
				TCell& lf = tab[row][col-wlo-1];
				updateColorHoriz(
					lf,        // cell to the left
					cur,       // destination cell
					r);        // ref mask associated with destination column
			}
			ncups_++;
			// 'cur' is now initialized
			FINALIZE_CELL(row, col);
		} // end loop over inner diagonals
		//
		// Handle the col == whi case (provided wlo != whi) after the
		// the prior inner loop to reduce the number of guards needed
		// inside it.
		//
		if(whi > wlo) {
			col = whi;
			TCell& cur = tab[row][col-wlo];
			cur.clear();
			const size_t fullcol = col + row;
			const int r = rf_[rfi_ + fullcol];
			TCell& dg = tab[row-1][col-wlo];
			if(!dg.empty) {
				updateColorDiag(
					dg,        // cell diagonally above and to the left
					cur,       // destination cell
					r,         // ref mask associated with destination column
					c,         // nucleotide in destination row
					mmpen);    // penalty to incur for color miscall
			}
			TCell& lf = tab[row][col-wlo-1];
			if(!onlyDiagInto && !lf.empty) {
				updateColorHoriz(
					lf,        // cell to the left
					cur,       // destination cell
					r);        // ref mask associated with destination column
			}
			ncups_++;
			// 'cur' is now initialized
			FINALIZE_CELL(row, col);
		}
		if(!validInRow) {
			assert_geq(rdf_ - rdi_, row);
			nrowskips_ += (rdf_ - rdi_ - row);
			return !soldone_;
		}
	}
	assert_eq(tab.size(), rd_->length()+1);
	return !soldone_;
}

