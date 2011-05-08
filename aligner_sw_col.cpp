/*
 *  aligner_sw_col.cpp
 */

#include "aligner_sw.h"
#include "aligner_sw_col.h"
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
pair<int, int> SwColorCellMask::randOverallBacktrack(
	RandomSource& rand,
	bool& branch,
	bool clear)
{
	ASSERT_ONLY(int num = numOverallPossible());
	std::pair<int, int> ret;
	ret.second = -1;
	int i =
		((oall_diag != 0) << 0) |
		((oall_rfop != 0) << 1) |
		((oall_rfex != 0) << 2) |
		((oall_rdop != 0) << 3) |
		((oall_rdex != 0) << 4);
	// Count the total number of different choices we could make
	int totChoices =
		alts5[oall_diag] +
		alts5[oall_rfop] +
		alts5[oall_rfex] +
		oall_rdop +
		oall_rdex;
	// If we're choosing from among >1 possibilities, inform caller so that
	// caller can add a frame to the backtrack stack.
	branch = totChoices > 1;
	ret.first = randFromMask(rand, i);
	assert_lt(ret.first, 5);
	assert(ret.first == SW_BT_OALL_DIAG ||
		   ret.first == SW_BT_OALL_REF_OPEN ||
		   ret.first == SW_BT_OALL_REF_EXTEND ||
		   ret.first == SW_BT_OALL_READ_OPEN ||
		   ret.first == SW_BT_OALL_READ_EXTEND);
	// Clear the bit associated with the path chosen
	if(ret.first < 3) {
		// Must choose character for next row
		if(ret.first == SW_BT_OALL_DIAG) {
			assert(oall_diag != 0);
			ret.second = randFromMask(rand, oall_diag);
			if(clear) oall_diag &= ~(1 << ret.second);
		} else if(ret.first == SW_BT_OALL_REF_OPEN) {
			assert(oall_rfop != 0);
			ret.second = randFromMask(rand, oall_rfop);
			if(clear) oall_rfop &= ~(1 << ret.second);
		} else if(ret.first == SW_BT_OALL_REF_EXTEND) {
			assert(oall_rfex != 0);
			ret.second = randFromMask(rand, oall_rfex);
			if(clear) oall_rfex &= ~(1 << ret.second);
		}
	} else if(clear) {
		if(ret.first == SW_BT_OALL_READ_OPEN) {
			oall_rdop = 0;
		}
		else if(ret.first == SW_BT_OALL_READ_EXTEND) {
			oall_rdex = 0;
		}
	}
	assert_eq(num - (clear ? 1 : 0), numOverallPossible());
	return ret;
}

/**
 * Select a path for backtracking from this cell.  If there is a
 * tie among eligible paths, break it randomly.  Return value is
 * a pair where first = a flag indicating the backtrack type (see
 * enum defining SW_BT_* above), and second = a selection for
 * the read character for the next row up.  second should be
 * ignored if the backtrack type is a gap in the read.
 */
pair<int, int> SwColorCellMask::randReadGapBacktrack(
	RandomSource& rand,
	bool& branch,
	bool clear)
{
	ASSERT_ONLY(int num = numReadGapPossible());
	std::pair<int, int> ret;
	ret.second = -1;
	int i =
		((rdgap_op != 0) << 0) |
		((rdgap_ex != 0) << 1);
	// Count the total number of different choices we could make
	int totChoices = rdgap_op + rdgap_ex;
	// If we're choosing from among >1 possibilities, inform caller so that
	// caller can add a frame to the backtrack stack.
	branch = totChoices > 1;
	ret.first = randFromMask(rand, i) + SW_BT_RDGAP_OPEN;
	assert(ret.first == SW_BT_RDGAP_OPEN ||
		   ret.first == SW_BT_RDGAP_EXTEND);
	if(clear) {
		if(ret.first == SW_BT_RDGAP_OPEN) {
			rdgap_op = 0;
		} else {
			rdgap_ex = 0;
		}
	}
	assert_eq(num - (clear ? 1 : 0), numReadGapPossible());
	return ret;
}

/**
 * Select a path for backtracking from this cell.  If there is a
 * tie among eligible paths, break it randomly.  Return value is
 * a pair where first = a flag indicating the backtrack type (see
 * enum defining SW_BT_* above), and second = a selection for
 * the read character for the next row up.  second should be
 * ignored if the backtrack type is a gap in the read.
 */
pair<int, int> SwColorCellMask::randRefGapBacktrack(
	RandomSource& rand,
	bool& branch,
	bool clear)
{
	ASSERT_ONLY(int num = numRefGapPossible());
	std::pair<int, int> ret;
	ret.second = -1;
	int i =
		((rfgap_op != 0) << 0) |
		((rfgap_ex != 0) << 1);
	// Count the total number of different choices we could make
	int totChoices =
		alts5[rfgap_op] +
		alts5[rfgap_ex];
	// If we're choosing from among >1 possibilities, inform caller so that
	// caller can add a frame to the backtrack stack.
	branch = totChoices > 1;
	ret.first = randFromMask(rand, i) + SW_BT_RFGAP_OPEN;
	assert(ret.first == SW_BT_RFGAP_OPEN ||
		   ret.first == SW_BT_RFGAP_EXTEND);
	if(ret.first == SW_BT_RFGAP_OPEN) {
		assert(rfgap_op != 0);
		ret.second = randFromMask(rand, rfgap_op);
		if(clear) rfgap_op &= ~(1 << ret.second);
	} else if(ret.first == SW_BT_RFGAP_EXTEND) {
		assert(rfgap_ex != 0);
		ret.second = randFromMask(rand, rfgap_ex);
		if(clear) rfgap_ex &= ~(1 << ret.second);
	}
	assert_eq(num - (clear ? 1 : 0), numRefGapPossible());
	return ret;
}

/**
 * We finished updating the cell; set empty and finalized appropriately.
 */
inline bool SwColorCell::finalize(TAlScore floorsc) {
	ASSERT_ONLY(finalized = true);
	for(int i = 0; i < 4; i++) {
		bool aboveFloor = oallBest[i].score() >= floorsc;
		if(!mask[i].empty() && aboveFloor) {
			assert(VALID_AL_SCORE(oallBest[i]));
			assert_geq(oallBest[i].score(), floorsc);
			empty = false;
#ifdef NDEBUG
			break;
#endif
		} else if(aboveFloor) {
			// This cell & nucleotide could be at the top end of an alignment
			terminal[i] = true;
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
	typedef SwColorCell TCell;
	ELList<TCell>& tab = ctab_;
	assert_lt(row, tab.size());
	assert_geq(escore, minsc_);
	btcstack_.clear();
	btcells_.clear();
	size_t tabcol = col - row;
	ASSERT_ONLY(res.alres.drd.clear());
	int curC = lastC;
	AlnScore score; score.score_ = 0;
	score.gaps_ = score.ns_ = 0;
	bool refExtend = false, readExtend = false;
	size_t origCol = col;
	size_t gaps = 0, readGaps = 0, refGaps = 0;
	res.alres.reset();
	EList<Edit>& ned = res.alres.ned();
	//EList<Edit>& aed = res.alres.aed();
	EList<Edit>& ced = res.alres.ced();
	res.ndn = lastC;
	size_t trimEnd = dpRows() - row - 1; 
	size_t trimBeg = 0;
	assert(!sc_->monotone || escore <= 0);
	assert(!sc_->monotone || score.score() >= escore);
	int ct = SW_BT_CELL_OALL; // cell type
	ASSERT_ONLY(res.alres.drd.append(curC));
	while((int)row >= 0) {
		nbts_++;
		assert_lt(row, tab.size());
		assert_leq(col, origCol);
		assert_range(0, 3, curC);
		assert_geq(col, row);
		tabcol = col - row;
		assert_geq(tabcol, 0);
		TCell& curc = tab[row][tabcol];
		bool empty = curc.mask[curC].numPossible(ct) == 0;
		assert_eq(gaps, Edit::numGaps(ned));
		assert_leq(gaps, rdgap_ + rfgap_);
		// Cell was involved in a previously-reported alignment?
		if(!curc.canMoveThrough(curC, ct)) {
			if(!btcstack_.empty()) {
				// Pop record off the top of the stack
				ned.resize(btcstack_.back().nedsz);
				ced.resize(btcstack_.back().cedsz);
				btcells_.resize(btcstack_.back().celsz);
				//aed.resize(btcstack_.back().aedsz);
				row      = btcstack_.back().row;
				col      = btcstack_.back().col;
				gaps     = btcstack_.back().gaps;
				readGaps = btcstack_.back().readGaps;
				refGaps  = btcstack_.back().refGaps;
				score    = btcstack_.back().score;
				curC     = btcstack_.back().curC;
				ct       = btcstack_.back().ct;
				ASSERT_ONLY(res.alres.drd.resize(btcstack_.back().drdsz));
				btcstack_.pop_back();
				assert(!sc_->monotone || score.score() >= escore);
				continue;
			} else {
				// No branch points to revisit; just give up
				return false;
			}
		}
		assert(!curc.reportedThru_);
		assert(!sc_->monotone || score.score() >= escore);
		if(row == 0) {
			btcells_.expand();
			btcells_.back().first = row;
			btcells_.back().second = tabcol;
			break;
		}
		if(empty) {
			assert_eq(SW_BT_CELL_OALL, ct);
			// This cell is at the end of a legitimate alignment
			trimBeg = row;
			assert_eq(btcells_.size(), dpRows() - trimBeg - trimEnd + readGaps);
			break;
		}
		bool branch = false;
		pair<int, int> cur =
			curc.mask[curC].randBacktrack(ct, rand, branch, true);
		if(branch) {
			assert_gt(curc.mask[curC].numPossible(ct), 0);
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
				readGaps,
				refGaps,
				score,
				ct
				ASSERT_ONLY(, res.alres.drd.length())
				);
		}
		btcells_.expand();
		btcells_.back().first = row;
		btcells_.back().second = tabcol;
		ASSERT_ONLY(TAlScore origScore = score.score_);
		switch(cur.first) {
			// Move up and to the left.  If the reference nucleotide in the
			// source row mismatches the decoded nucleotide, penalize
			// it and add a nucleotide mismatch.  If the color being
			// traversed is a miscall, penalize that.
			case SW_BT_OALL_DIAG: {
				assert_neq(-1, cur.second);
				refExtend = readExtend = false;
				assert_gt(row, 0); assert_gt(col, 0);
				// Check for base mismatch at source (lower-right) cell
				int refNmask = (int)rf_[rfi_+col];
				assert_gt(refNmask, 0);
				int m = matchesEx(curC, refNmask);
				ASSERT_ONLY(int lastct = ct);
				ct = SW_BT_CELL_OALL;
				if(m != 1) {
					Edit e((int)row, mask2dna[refNmask], "ACGTN"[curC], EDIT_TYPE_MM);
					assert(e.repOk());
					assert(ned.empty() || ned.back().pos >= row);
					ned.push_back(e);
					int pen = ((m == 0) ? sc_->snp : sc_->n(30));
					score.score_ -= pen;
					assert(!sc_->monotone || score.score() >= escore);
				} else {
					int64_t bonus = sc_->match(30);
					score.score_ += bonus;
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
					int pen = QUAL(row);
					score.score_ -= pen; // color mismatch
					assert(!sc_->monotone || score.score() >= escore);
					assert(e.repOk());
					assert(ced.empty() || ced.back().pos >= row);
					ced.push_back(e);
				}
				if(readC > 3) {
					score.ns_++;
				}
				assert_lt(tabcol, tab[row].size());
				assert(VALID_AL_SCORE(score));
				lastC = curC;
				curC = cur.second;
				ASSERT_ONLY(int64_t scoreDiff = origScore - score.score_);
				assert_eq(scoreDiff, tab[row][col-row].best(curC, ct).score() - tab[row+1][col-row].best(lastC, lastct).score());
				ASSERT_ONLY(res.alres.drd.append(curC));
				break;
			}
			// Move up.  Add an edit encoding the ref gap.  If the color
			// being traversed is a miscall, penalize that.
			case SW_BT_OALL_REF_OPEN:
			case SW_BT_RFGAP_OPEN:
			{
				refExtend = true; readExtend = false;
				assert_gt(row, 0);
				Edit e((int)row, '-', "ACGTN"[curC], EDIT_TYPE_REF_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				row--;
				ASSERT_ONLY(int lastct = ct);
				ct = SW_BT_CELL_OALL;
				// Check for color miscall
				int decC = dinuc2color[cur.second][curC];
				assert_range(0, 3, decC);
				int rdC = (*rd_)[rdi_+row];
				if(decC != rdC) {
					score.score_ -= QUAL(row);
					assert(!sc_->monotone || score.score() >= escore);
					Edit e((int)row, "ACGT"[decC], "ACGTN"[rdC], EDIT_TYPE_MM);
					assert(e.repOk());
					assert(ced.empty() || ced.back().pos >= row);
					ced.push_back(e);
				}
				if(rdC > 3) {
					score.ns_++;
				}
				score.score_ -= sc_->refGapOpen();
				assert(!sc_->monotone || score.score() >= escore);
				gaps++; refGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				assert_lt(row, tab.size());
				assert_lt(col - row, tab[row].size());
				lastC = curC;
				curC = cur.second;
				ASSERT_ONLY(int64_t scoreDiff = origScore - score.score_);
				assert_eq(scoreDiff, tab[row][col-row].best(curC, ct).score() - tab[row+1][col-(row+1)].best(lastC, lastct).score());
				ASSERT_ONLY(res.alres.drd.append(curC));
				break;
			}
			// Move up.  Add an edit encoding the ref gap.  If the color
			// being traversed is a miscall, penalize that.
			case SW_BT_OALL_REF_EXTEND:
			case SW_BT_RFGAP_EXTEND:
			{
				refExtend = true; readExtend = false;
				assert_gt(row, 1);
				Edit e((int)row, '-', "ACGTN"[curC], EDIT_TYPE_REF_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				row--;
				ASSERT_ONLY(int lastct = ct);
				ct = SW_BT_CELL_RFGAP;
				// Check for color mismatch
				int decC = dinuc2color[cur.second][curC];
				assert_range(0, 3, decC);
				int rdC = (*rd_)[rdi_+row];
				if(decC != rdC) {
					score.score_ -= QUAL(row);
					assert(!sc_->monotone || score.score() >= escore);
					Edit e((int)row, "ACGT"[decC], "ACGTN"[rdC], EDIT_TYPE_MM);
					assert(e.repOk());
					assert(ced.empty() || ced.back().pos >= row);
					ced.push_back(e);
				}
				if(rdC > 3) {
					score.ns_++;
				}
				score.score_ -= sc_->refGapExtend();
				assert(!sc_->monotone || score.score() >= escore);
				gaps++; refGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				assert_lt(row, tab.size());
				assert_lt(col - row, tab[row].size());
				lastC = curC;
				curC = cur.second;
				ASSERT_ONLY(int64_t scoreDiff = origScore - score.score_);
				assert_eq(scoreDiff, tab[row][col-row].best(curC, ct).score() - tab[row+1][col-(row+1)].best(lastC, lastct).score());
				ASSERT_ONLY(res.alres.drd.append(curC));
				break;
			}
			case SW_BT_OALL_READ_OPEN:
			case SW_BT_RDGAP_OPEN:
			{
				refExtend = false; readExtend = true;
				assert_gt(col, 0);
				Edit e((int)row+1, "ACGTN"[firsts5[(int)rf_[rfi_+col]]], '-', EDIT_TYPE_READ_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				col--;
				ASSERT_ONLY(int lastct = ct);
				ct = SW_BT_CELL_OALL;
				score.score_ -= sc_->readGapOpen();
				assert(!sc_->monotone || score.score() >= escore);
				gaps++; readGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				assert_lt(row, tab.size());
				assert_lt(col - row, tab[row].size());
				lastC = curC;
				ASSERT_ONLY(int64_t scoreDiff = origScore - score.score_);
				assert_eq(scoreDiff, tab[row][col-row].best(curC, ct).score() - tab[row][(col+1)-row].best(lastC, lastct).score());
				break;
			}
			case SW_BT_OALL_READ_EXTEND:
			case SW_BT_RDGAP_EXTEND:
			{
				refExtend = false; readExtend = true;
				assert_gt(col, 1);
				Edit e((int)row+1, "ACGTN"[firsts5[(int)rf_[rfi_+col]]], '-', EDIT_TYPE_READ_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				col--;
				ASSERT_ONLY(int lastct = ct);
				ct = SW_BT_CELL_RDGAP;
				score.score_ -= sc_->readGapExtend();
				assert(!sc_->monotone || score.score() >= escore);
				gaps++; readGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				assert_lt(row, tab.size());
				assert_lt(col - row, tab[row].size());
				lastC = curC;
				ASSERT_ONLY(int64_t scoreDiff = origScore - score.score_);
				assert_eq(scoreDiff, tab[row][col-row].best(curC, ct).score() - tab[row][(col+1)-row].best(lastC, lastct).score());
				break;
			}
			default: throw 1;
		}
	} // while((int)row > 0)
	assert_geq(col, 0);
	assert_eq(SW_BT_CELL_OALL, ct);
	// The number of cells in the backtracs should equal the number of read
	// bases after trimming plus the number of gaps
	assert_eq(btcells_.size(), dpRows() - trimBeg - trimEnd + readGaps);
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
	int pen = 0;
	if(m != 1) {
		Edit e((int)row, mask2dna[refNmask], "ACGTN"[curC], EDIT_TYPE_MM);
		assert(e.repOk());
		assert(ned.empty() || ned.back().pos >= row);
		ned.push_back(e);
		pen = ((m == 0) ? sc_->snp : sc_->n(30));
		score.score_ -= pen;
		assert(!sc_->monotone || score.score() >= escore);
	}
	// TODO: check that pen matches score in final cell?
	if(m == -1) {
		score.ns_++;
	}
	ASSERT_ONLY(res.alres.drd.reverse());
	res.nup = curC;
	res.reverse();
	assert(Edit::repOk(ced, (*rd_)));
#ifndef NDEBUG
	size_t gapsCheck = 0;
	for(size_t i = 0; i < ned.size(); i++) {
		if(ned[i].isGap()) gapsCheck++;
	}
	assert_eq(gaps, gapsCheck);
	static BTDnaString refstr;
	refstr.clear();
	for(size_t i = col; i <= origCol; i++) {
		refstr.append(firsts5[(int)rf_[rfi_+i]]);
	}
	static BTDnaString editstr;
	editstr.clear();
	Edit::toRef(res.alres.drd, ned, editstr);
	if(refstr != editstr) {
		cerr << "Decoded nucleotides and edits don't match reference:" << endl;
		cerr << "           score: " << score.score() << " (" << gaps << " gaps)" << endl;
		cerr << "           edits: ";
		Edit::print(cerr, ned);
		cerr << endl;
		cerr << "    decoded nucs: " << res.alres.drd << endl;
		cerr << "     edited nucs: " << editstr << endl;
		cerr << "  reference nucs: " << refstr << endl;
		assert(0);
	}
#endif
	// done
	assert_eq(score.score(), escore);
	assert_leq(gaps, rdgap_ + rfgap_);
	// Dummy values for refid and fw
	off = (col - row);
	score.gaps_ = gaps;
	res.alres.setScore(score);
	res.alres.setShape(
		refidx_,                  // ref id
		off + rfi_ + refoff_,     // 0-based ref offset
		fw_,                      // aligned to Watson?
		rdf_ - rdi_,              // read length
		color_,                   // read was colorspace?
		true,                     // pretrim soft?
		0,                        // pretrim 5' end
		0,                        // pretrim 3' end
		false,                    // alignment trim soft?
		fw_ ? trimBeg : trimEnd,  // alignment trim 5' end
		fw_ ? trimEnd : trimBeg); // alignment trim 3' end
	size_t refns = 0;
	for(size_t i = col; i <= origCol; i++) {
		if((int)rf_[rfi_+i] > 15) {
			refns++;
		}
	}
	res.alres.setRefNs(refns);
	assert(res.repOk());
	return true;
}

#if 0
// Count it
#define IF_COUNT_N(x...) x
#define COUNT_N(x) { \
	if(ninvolved) { \
		x.incNs(nceil_); \
	} \
}
#else
// Do nothing
#define IF_COUNT_N(x...)
#define COUNT_N(x)
#endif

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
	IF_COUNT_N(bool ninvolved = rfm > 15);
	for(int to = 0; to < 4; to++) {
		AlnScore leftOallBest  = lc.oallBest[to];
		AlnScore leftRdgapBest = lc.rdgapBest[to];
		AlnScore& myOallBest   = dstc.oallBest[to];
		AlnScore& myRdgapBest  = dstc.rdgapBest[to];
		SwColorCellMask& myMask  = dstc.mask[to];
		assert(!sc_->monotone || leftOallBest.score()  <= 0);
		assert(!sc_->monotone || leftRdgapBest.score() <= 0);
		assert_leq(leftRdgapBest, leftOallBest);
		if(!VALID_AL_SCORE(leftOallBest)) {
			assert(!VALID_AL_SCORE(leftRdgapBest));
			return;
		}
		COUNT_N(myOallBest);
		COUNT_N(myRdgapBest);
		// Consider read gap extension
		AlnScore ex = leftRdgapBest - sc_->readGapExtend();
		if(ex.score_ >= floorsc_) {
			if(ex >= myOallBest) {
				if(ex > myOallBest) {
					myMask.clearOverallMask();
					myOallBest = ex;
				}
				myMask.oall_rdex = 1;
			}
			if(ex >= myRdgapBest) {
				if(ex > myRdgapBest) {
					myMask.clearReadGapMask();
					myRdgapBest = ex;
				}
				myMask.rdgap_ex = 1;
			}
		}
		// Consider read gap open
		ex = leftOallBest - sc_->readGapOpen();
		if(ex.score_ >= floorsc_) {
			if(ex >= myOallBest) {
				if(ex > myOallBest) {
					myMask.clearOverallMask();
					myOallBest = ex;
				}
				myMask.oall_rdop = 1;
			}
			if(ex >= myRdgapBest) {
				if(ex > myRdgapBest) {
					myMask.clearReadGapMask();
					myRdgapBest = ex;
				}
				myMask.rdgap_op = 1;
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
		AlnScore fromOall[] = {
			uc.oallBest[0], uc.oallBest[1],
			uc.oallBest[2], uc.oallBest[3] };
		AlnScore fromRfgap[] = {
			uc.rfgapBest[0], uc.rfgapBest[1],
			uc.rfgapBest[2], uc.rfgapBest[3] };
		// Calculate which 'from' character isn't going to cause a color mm
		const int goodfrom = nuccol2nuc[to][c];
		// Neutralize penalties for 'from' chars that corroborate color
		if(goodfrom < 4) {
			if(fromOall [goodfrom].valid()) fromOall [goodfrom].score_ += penmm;
			if(fromRfgap[goodfrom].valid()) fromRfgap[goodfrom].score_ += penmm;
		}
		for(int fr = 0; fr < 4; fr++) {
			AlnScore  upOallBest  = fromOall[fr] - penmm;
			AlnScore  upRfgapBest = fromRfgap[fr] - penmm;
			AlnScore& myOallBest  = dstc.oallBest[to];
			AlnScore& myRfgapBest = dstc.rfgapBest[to];
			SwColorCellMask& myMask = dstc.mask[to];
			assert(!sc_->monotone || upOallBest.score()  <= 0);
			assert(!sc_->monotone || upRfgapBest.score() <= 0);
			assert_leq(upRfgapBest, upOallBest);
			if(!VALID_AL_SCORE(upOallBest)) {
				assert(!VALID_AL_SCORE(upRfgapBest));
				return;
			}
			COUNT_N(myOallBest);
			COUNT_N(myRfgapBest);
			// Consider reference gap extension
			AlnScore ex = upRfgapBest - sc_->refGapExtend();
			if(ex.score_ >= floorsc_) {
				if(ex >= myOallBest) {
					if(ex > myOallBest) {
						myMask.clearOverallMask();
						myOallBest = ex;
					}
					myMask.oall_rfex |= (1 << fr);
				}
				if(ex >= myRfgapBest) {
					if(ex > myRfgapBest) {
						myMask.clearRefGapMask();
						myRfgapBest = ex;
					}
					myMask.rfgap_ex |= (1 << fr);
				}
			}
			// Consider reference gap open
			ex = upOallBest - sc_->refGapOpen();
			if(ex.score_ >= floorsc_) {
				if(ex >= myOallBest) {
					if(ex > myOallBest) {
						myMask.clearOverallMask();
						myOallBest = ex;
					}
					myMask.oall_rfop |= (1 << fr);
				}
				if(ex >= myRfgapBest) {
					if(ex > myRfgapBest) {
						myMask.clearRefGapMask();
						myRfgapBest = ex;
					}
					myMask.rfgap_op |= (1 << fr);
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
	int                penmm,     // penalty to incur for color miscall
	bool&              improved)
{
	assert(uc.finalized);
	if(uc.empty) return;
	IF_COUNT_N(bool ninvolved = (c > 3 || refMask > 15));
	
	// TODO!!!
	improved = true;
	
	for(int to = 0; to < 4; to++) {
		// Assuming that the read character in this row is 'to'...
		TAlScore add =
			(matches(to, refMask) ?
				sc_->match(30) :
				-((refMask > 15) ? sc_->n(30) : sc_->snp));
		AlnScore fromOall[] = {
			uc.oallBest[0], uc.oallBest[1],
			uc.oallBest[2], uc.oallBest[3] };
		const int goodfrom = nuccol2nuc[to][c];
		if(goodfrom < 4) {
			if(fromOall[goodfrom].valid()) {
				// Neutralize penalty for 'from' char that corroborates color
				fromOall[goodfrom].score_ += penmm;
			}
		}
		for(int fr = 0; fr < 4; fr++) {
			AlnScore  dcOallBest  = fromOall[fr] - penmm;
			AlnScore& myOallBest  = dstc.oallBest[to];
			SwColorCellMask& myMask = dstc.mask[to];
			assert(!sc_->monotone || dcOallBest.score() <= 0);
			if(!VALID_AL_SCORE(dcOallBest)) {
				continue;
			}
			COUNT_N(dcOallBest);
			AlnScore ex = dcOallBest + (int)add;
			if(ex.score_ >= floorsc_ && ex >= myOallBest) {
				if(ex > myOallBest) {
					myOallBest = ex;
					myMask.clear();
				}
				myMask.oall_diag |= (1 << fr);
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
#define UPDATE_SOLS(cur, row, col, improved) { \
	assert_lt(col, width_); \
	if(en_ == NULL || (*en_)[col]) { \
		/* Column is acceptable */ \
		for(int I = 0; I < 4; I++) { \
			if(cur.oallBest[I].score() >= minsc_) { \
				/* Score and column are acceptable */ \
				const bool local = !sc_->monotone; \
				/* For local alignment, a cell is only a solution candidate */ \
				/* if the score was improved when we moved into the cell. */ \
				if(!local || improved) { \
					/* Score is acceptable */ \
					if((int64_t)row >= solrowlo_) { \
						/* Row is acceptable */ \
						/* This cell is now a good candidate for backtrace BUT */ \
						/* in local alignment mode we still need to know */ \
						/* whether this solution is improved upon in the next */ \
						/* row.  If it is improved upon later, this flag is set */ \
						/* to false later. */ \
						tab[row][col].backtraceCandidate = true; \
						if(row > 0) tab[row-1][col].backtraceCandidate = false; \
						assert(repOk()); \
					} \
				} \
			} \
		} \
	} \
}

// c is the column offset with respect to the LHS of the rectangle; for offset
// w/r/t LHS of parallelogram, use r+c
#define FINALIZE_CELL(r, c, improved) { \
	assert_lt(col, width_); \
	assert(!tab[r][c].finalized); \
	if(tab[r][c].finalize(floorsc_)) { \
		UPDATE_SOLS(tab[r][c], r, c, improved); \
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
void SwAligner::alignColors(RandomSource& rnd) {
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
	const size_t wlo = 0;
	const size_t whi = (int)(width_ - 1);
	tab[0].resize(whi-wlo+1); // add columns to first row
	//TAlScore loBound = sc_->monotone ? minsc_ : floorsc_;
	bool validInRow = !sc_->monotone;;
	// Calculate starting values for the rest of the columns in the
	// first row.  No need to consider colors yet since we won't have
	// consecutive pairs of read nucleotides until we get to the second
	// row.  We just consider reference mutations.
	for(size_t col = 0; col <= whi; col++) {
		TCell& curc = tab[0][col];
		curc.clear(); // clear the cell; masks and scores
		int rfm = rf_[rfi_+col];
		// Can we start from here?
		bool canStart = (st_ == NULL || (*st_)[col]);
		bool improved = false;
		if(canStart) {
			for(int to = 0; to < 4; to++) {
				curc.oallBest[to].invalidate();
				curc.rdgapBest[to].invalidate();
				curc.rfgapBest[to].invalidate();
				int m = matchesEx(to, rfm);
				if(m == 1) {
					// The assigned subject nucleotide matches the reference;
					// no penalty
					assert_lt(rfm, 16);
					curc.oallBest[to].gaps_ = 0;
					curc.oallBest[to].ns_ = 0;
					curc.oallBest[to].score_ = 0;
					curc.mask[to].oall_diag = 0xf;
					improved = true;
				} else if(m == 0 && sc_->snp <= -minsc_) {
					// The assigned subject nucleotide does not match the
					// reference nucleotide, so we add a SNP penalty
					curc.oallBest[to].gaps_ = 0;
					curc.oallBest[to].ns_ = 0;
					curc.oallBest[to].score_ = -sc_->snp;
					curc.mask[to].oall_diag = 0xf;
				} else if(m == -1) {
					// The mask or query character is an N
					curc.oallBest[to].gaps_ = 0;
					curc.oallBest[to].ns_ = 1;
					curc.oallBest[to].score_ = -sc_->n(30);
					curc.mask[to].oall_diag = 0xf;
				} else {
					// Leave mask[to] cleared
				}
			}
		}
		// Calculate horizontals if barrier allows
		if(sc_->gapbar < 1 && col > 0) {
			updateColorHoriz(tab[0][col-1], curc, rfm);
			ncups_++;
		}
		FINALIZE_CELL(0, col, improved);
	}
	nrowups_++;
	if(!validInRow) {
		nrowskips_ += (rdf_ - rdi_);
		return;
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
		bool improved = false;
		if(!tab[row-1][0].empty) {
			const size_t fullcol = col + row;
			updateColorDiag(
				tab[row-1][0],     // cell diagonally above and to the left
				cur,               // destination cell
				rf_[rfi_ + fullcol], // ref mask associated with destination cell
				c,                 // color being traversed
				mmpen,             // penalty to incur for color miscall
				improved);
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
		FINALIZE_CELL(row, col, improved);
		// Iterate from leftmost to rightmost inner diagonals
		for(col = wlo+1; col < whi; col++) {
			const size_t fullcol = col + row;
			int r = rf_[rfi_ + fullcol];
			TCell& cur = tab[row][col-wlo];
			cur.clear();
			TCell& dg = tab[row-1][col-wlo];
			improved = false;
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
				mmpen,         // penalty to incur for color miscall
				improved);
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
			FINALIZE_CELL(row, col, improved);
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
			improved = false;
			if(!dg.empty) {
				updateColorDiag(
					dg,        // cell diagonally above and to the left
					cur,       // destination cell
					r,         // ref mask associated with destination column
					c,         // nucleotide in destination row
					mmpen,     // penalty to incur for color miscall
					improved);
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
			FINALIZE_CELL(row, col, improved);
		}
		if(!validInRow) {
			assert_geq(rdf_ - rdi_, row);
			nrowskips_ += (rdf_ - rdi_ - row);
			return;
		}
	}
	assert_eq(tab.size(), rd_->length()+1);
}

