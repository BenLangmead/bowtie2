/*
 *  aligner_sw_col.cpp
 */

#include "aligner_sw.h"
#include "search_globals.h"
#include "penalty.h"
#include "mask.h"

/**
 * Select a path for backtracking from this cell.  If there is a
 * tie among eligible paths, break it randomly.  Return value is
 * a pair where first = a flag indicating the backtrack type (see
 * enum defining SW_BT_* above), and second = a selection for
 * the read character for the next row up.  second should be
 * ignored if the backtrack type is a gap in the read.
 */
pair<int, int> SwColorCellMask::randBacktrack(RandomSource& rand) {
	ASSERT_ONLY(int num = numPossible());
	std::pair<int, int> ret;
	ret.second = -1;
	int i = ((diag != 0) << 0) |
			((rfop != 0) << 1) |
			((rfex != 0) << 2) |
			((rdop != 0) << 3) |
			((rdex != 0) << 4);
	ret.first = randFromMask(rand, i);
	assert_lt(ret.first, 5);
	assert(ret.first == SW_BT_DIAG ||
		   ret.first == SW_BT_REF_OPEN ||
		   ret.first == SW_BT_REF_EXTEND ||
		   ret.first == SW_BT_READ_OPEN ||
		   ret.first == SW_BT_READ_EXTEND);
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
 * Given the dynamic programming table, trace backwards from the lower
 * right-hand corner and populate 'decoded' 'nedits', and 'cedits'
 * accordingly.
 *
 * The approach is very similar to the one described in the SHRiMP
 * paper:
 *
 * Rumble SM, Lacroute P, Dalca AV, Fiume M, Sidow A, Brudno M. SHRiMP:
 * accurate mapping of short color-space reads. PLoS Comput Biol. 2009
 * May;5(5)
 *
 * The return value is the alignment's upstream-most reference
 * character's offset with respect to rfi.
 */
int SwAligner::backtrackColors(
	const BTDnaString& rd, // read sequence
	const BTString& qu,    // read qualities
	size_t rdi,            // offset of first read char to align
	size_t rdf,            // offset of last read char to align
	BTString& rf,          // reference sequence
	size_t rfi,            // offset of first reference char to align to
	size_t rff,            // offset of last reference char to align to
	AlignmentScore expectScore,   // score we expect to get over backtrack
	int readGaps,          // max # gaps in read
	int refGaps,           // max # gaps in ref
	const SwParams& pa,    // params for SW alignment
	const Penalties& pen,  // penalties for edit types
	int& nup,              // upstream decoded nucleotide
	int& ndn,              // downstream decoded nucleotide
	SwResult& res,         // store results (edits and scores) here
	int col,               // start in this column (w/r/t the full matrix)
	int lastC,             // character to backtrace from in lower-right corner
	RandomSource& rand)    // pseudo-random generator
{
	ELList<SwColorCell>& tab = ctab_;
	int row = rd.length();
	assert_eq(row, (int)tab.size()-1);
	ASSERT_ONLY(BTDnaString drd);
	ASSERT_ONLY(drd.resize(rdf-rdi+1));
	int curC = lastC;
	AlignmentScore score; score.score_ = 0;
	score.gaps_ = score.ns_ = 0;
	bool refExtend = false, readExtend = false;
	ASSERT_ONLY(int origCol = col);
	ASSERT_ONLY(int gaps = 0);
	res.alres.reset();
	EList<Edit>& ned = res.alres.ned();
	//EList<Edit>& aed = res.alres.aed();
	EList<Edit>& ced = res.alres.ced();
	ndn = lastC;
	while(row > 0) {
		res.swbts++;
		assert_range(0, (int)tab.size()-1, row);
		assert_range(0, origCol, col);
		assert_range(0, 3, curC);
		ASSERT_ONLY(drd.set(curC, row));
		int tabcol = col - row;
		assert_geq(tabcol, 0);
		ASSERT_ONLY(AlignmentScore scoreThisRound = score);
		ASSERT_ONLY(AlignmentScore bestThisRound = tab[row][tabcol].best[curC]);
		assert_gt(tab[row][tabcol].mask[curC].numPossible(), 0);
		pair<int, int> cur =
			tab[row][tabcol].mask[curC].randBacktrack(rand);
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
				int refNmask = (int)rf[rfi+col];
				int m = matches(curC, refNmask);
				if(m != 1) {
					Edit e(row, mask2dna[refNmask], "ACGTN"[curC], EDIT_TYPE_MM);
					assert(e.repOk());
					ned.push_back(e);
					score.score_ -= ((m == 0) ? pen.snp : pen.n(30));
				}
				if(m == -1) {
					score.ns_++;
				}
				row--; col--;
				assert_range(0, (int)tab.size()-1, row);
				assert(VALID_AL_SCORE(score));
				// Check for color mismatch
				int readC = rd[rdi + row];
				int decC = dinuc2color[cur.second][curC];
				assert_range(0, 3, decC);
				if(decC != readC) {
					Edit e(row, "ACGT"[decC], "ACGTN"[readC], EDIT_TYPE_MM);
					score.score_ -= QUAL(row); // color mismatch
					assert(e.repOk());
					ced.push_back(e);
				}
				if(readC > 3) {
					score.ns_++;
				}
				assert_range(0, (int)tab[row].size()-1, tabcol);
				assert(VALID_AL_SCORE(score));
				ASSERT_ONLY(scoreThisRound = score - scoreThisRound);
				ASSERT_ONLY(bestThisRound -= tab[row][tabcol].best[cur.second]);
				// Make sure that both changed in the same way
				assert_eq(scoreThisRound, bestThisRound);
				lastC = curC;
				curC = cur.second;
				break;
			}
			// Move up.  Add an edit encoding the ref gap.  If the color
			// being traversed is a miscall, penalize that.
			case SW_BT_REF_OPEN: {
				refExtend = true; readExtend = false;
				assert_gt(row, 0);
				Edit e(row, '-', "ACGTN"[curC], EDIT_TYPE_DEL);
				assert(e.repOk());
				ned.push_back(e);
				assert_geq(row, pa.gapBar);
				assert_geq((int)(rdf-rdi-row-1), pa.gapBar-1);
				row--;
				// Check for color miscall
				int decC = dinuc2color[cur.second][curC];
				assert_range(0, 3, decC);
				int rdC = rd[rdi+row];
				if(decC != rdC) {
					score.score_ -= QUAL(row);
					Edit e(row, "ACGT"[decC], "ACGTN"[rdC], EDIT_TYPE_MM);
					assert(e.repOk());
					ced.push_back(e);
				}
				if(rdC > 3) {
					score.ns_++;
				}
				score.score_ -= pen.refOpen;
				score.gaps_++;
#ifndef NDEBUG
				gaps++;
				assert_leq(score.gaps_, readGaps + refGaps);
				assert_range(0, (int)tab.size()-1, row);
				tabcol = col - row;
				assert_range(0, (int)tab[row].size()-1, tabcol);
				scoreThisRound = score - scoreThisRound;
				bestThisRound -= tab[row][tabcol].best[cur.second];
				assert_eq(scoreThisRound, bestThisRound);
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
				Edit e(row, '-', "ACGTN"[curC], EDIT_TYPE_DEL);
				assert(e.repOk());
				ned.push_back(e);
				assert_geq(row, pa.gapBar);
				assert_geq((int)(rdf-rdi-row-1), pa.gapBar-1);
				row--;
				// Check for color mismatch
				int decC = dinuc2color[cur.second][curC];
				assert_range(0, 3, decC);
				int rdC = rd[rdi+row];
				if(decC != rdC) {
					score.score_ -= QUAL(row);
					Edit e(row, "ACGT"[decC], "ACGTN"[rdC], EDIT_TYPE_MM);
					assert(e.repOk());
					ced.push_back(e);
				}
				if(rdC > 3) {
					score.ns_++;
				}
				score.score_ -= pen.refExConst;
				score.gaps_++;
#ifndef NDEBUG
				gaps++;
				assert_leq(score.gaps_, readGaps + refGaps);
				assert_range(0, (int)tab.size()-1, row);
				tabcol = col - row;
				assert_range(0, (int)tab[row].size()-1, tabcol);
				scoreThisRound = score - scoreThisRound;
				bestThisRound -= tab[row][tabcol].best[cur.second];
				assert_eq(scoreThisRound, bestThisRound);
#endif
				lastC = curC;
				curC = cur.second;
				break;
			}
			case SW_BT_READ_OPEN: {
				refExtend = false; readExtend = true;
				assert_gt(col, 0);
				Edit e(row+1, "ACGTN"[firsts5[(int)rf[rfi+col]]], '-', EDIT_TYPE_INS);
				assert(e.repOk());
				ned.push_back(e);
				assert_geq(row, pa.gapBar);
				assert_geq((int)(rdf-rdi-row-1), pa.gapBar-1);
				col--;
				score.score_ -= pen.readOpen;
				score.gaps_++;
#ifndef NDEBUG
				gaps++;
				assert_leq(score.gaps_, readGaps + refGaps);
				assert_range(0, (int)tab.size()-1, row);
				tabcol = col - row;
				assert_range(0, (int)tab[row].size()-1, tabcol);
				scoreThisRound = score - scoreThisRound;
				bestThisRound -= tab[row][tabcol].best[curC];
				assert_eq(scoreThisRound, bestThisRound);
#endif
				lastC = curC;
				break;
			}
			case SW_BT_READ_EXTEND: {
				refExtend = false; readExtend = true;
				assert_gt(col, 1);
				Edit e(row+1, "ACGTN"[firsts5[(int)rf[rfi+col]]], '-', EDIT_TYPE_INS);
				assert(e.repOk());
				ned.push_back(e);
				assert_geq(row, pa.gapBar);
				assert_geq((int)(rdf-rdi-row-1), pa.gapBar-1);
				col--;
				score.score_ -= pen.readExConst;
				score.gaps_++;
#ifndef NDEBUG
				gaps++;
				assert_leq(score.gaps_, readGaps + refGaps);
				assert_range(0, (int)tab.size()-1, row);
				tabcol = col - row;
				assert_range(0, (int)tab[row].size()-1, tabcol);
				scoreThisRound = score - scoreThisRound;
				bestThisRound -= tab[row][tabcol].best[curC];
				assert_eq(scoreThisRound, bestThisRound);
#endif
				lastC = curC;
				break;
			}
			default: throw 1;
		}
	}
	assert_eq(0, row);
	int refNmask = (int)rf[rfi+col]; // get last char in ref involved in alignment
	int m = matches(curC, refNmask);
	if(m != 1) {
		Edit e(row, mask2dna[refNmask], "ACGTN"[curC], EDIT_TYPE_MM);
		assert(e.repOk());
		ned.push_back(e);
		score.score_ -= ((m == 0) ? pen.snp : pen.n(30));
	}
	if(m == -1) {
		score.ns_++;
	}
	ASSERT_ONLY(drd.set(curC, 0));
	nup = curC;
	res.reverse();
	assert(Edit::repOk(ced, rd));
#ifndef NDEBUG
	BTDnaString refstr;
	for(int i = col; i <= origCol; i++) {
		refstr.append(firsts5[(int)rf[rfi+i]]);
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
	assert_eq(score.score(), expectScore.score());
	//assert_leq(score.gaps(), expectScore.gaps());
	assert_leq(gaps, readGaps + refGaps);
	// Dummy values for refid and fw
	res.alres.setScore(score);
	return col;
}

/**
 * Update a SwCell's best[] and mask[] arrays with respect to its
 * neighbor on the left.
 */
inline void SwColorCell::updateHoriz(
	const SwColorCell& lc,
	int                rfm,
	const Penalties&   pen,
	size_t             nceil,
	int                penceil)
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
		AlignmentScore leftBest = lc.best[to];
		const SwColorCellMask& fromMask = lc.mask[to];
		AlignmentScore& myBest = best[to];
		SwColorCellMask& myMask = mask[to];
		assert_leq(leftBest.score(), 0);
		if(ninvolved) leftBest.incNs(nceil);
		if(!VALID_AL_SCORE(leftBest)) continue;
		// *Don't* penalize for a nucleotide mismatch because we must
		// have already done that in a previous vertical or diagonal
		// step.
		if(fromMask.readExtendPossible()) {
			// Read gap extension possible?
			AlignmentScore ex = leftBest;
			assert_leq(ex.score(), 0);
			assert(VALID_AL_SCORE(ex));
			ex.score_ -= pen.readExConst;
			assert_leq(ex.score(), 0);
			if(-ex.score_ <= penceil && ex >= myBest) {
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
			AlignmentScore ex = leftBest;
			assert_leq(ex.score_, 0);
			assert(VALID_AL_SCORE(ex));
			ex.score_ -= pen.readOpen;
			assert_leq(ex.score_, 0);
			if(-ex.score_ <= penceil && ex >= myBest) {
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
inline void SwColorCell::updateVert(
	const SwColorCell& uc,        // cell above
	int                c,         // color b/t this row, one above
	int                penmm,     // penalty to incur for color miscall
	const Penalties&   pen,
	size_t             nceil,
	int                penceil)
{
	assert(uc.finalized);
	if(uc.empty) return;
	for(int to = 0; to < 4; to++) {
		// Assuming that the read character in this row is 'to'...
		// No SNP penalty because destination read char aligns to a
		// gap in the reference.
		AlignmentScore from[] = {
			uc.best[0], uc.best[1],
			uc.best[2], uc.best[3] };
		assert(!VALID_AL_SCORE(from[0]) || from[0].score_ <= 0);
		assert(!VALID_AL_SCORE(from[1]) || from[1].score_ <= 0);
		assert(!VALID_AL_SCORE(from[2]) || from[2].score_ <= 0);
		assert(!VALID_AL_SCORE(from[3]) || from[3].score_ <= 0);
		// Calculate which 'from' character isn't going to
		// cause a color mismatch
		const int goodfrom = nuccol2nuc[to][c];
		// Reward the 'from' that corroborates the color
		if(goodfrom < 4) from[goodfrom].score_ += penmm;
		for(int fr = 0; fr < 4; fr++) {
			AlignmentScore frBest = from[fr];
			const SwColorCellMask& frMask = uc.mask[fr];
			AlignmentScore& myBest = best[to];
			SwColorCellMask& myMask = mask[to];
			if(c > 3) {
				frBest.incNs(nceil);
			}
			if(!VALID_AL_SCORE(frBest)) continue;
			frBest.score_ -= penmm;
			assert_leq(frBest.score_, 0);
			if(frMask.refExtendPossible()) {
				// Extend is possible
				frBest.score_ -= pen.refExConst;
				assert_leq(frBest.score_, 0);
				if(-frBest.score_ <= penceil && frBest >= myBest) {
					if(frBest > myBest) {
						myBest = frBest;
						myMask.clear();
					}
					myMask.rfex |= (1 << fr);
					assert(VALID_AL_SCORE(myBest));
				}
				// put it back
				frBest.score_ += pen.refExConst;
			}
			if(frMask.refOpenPossible()){
				// Open is possible
				frBest.score_ -= pen.refOpen;
				if(-frBest.score_ <= penceil && frBest >= myBest) {
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
inline void SwColorCell::updateDiag(
	const SwColorCell& uc,        // cell above and to the left
	int                refMask,   // ref mask associated with destination cell
	int                c,         // color being traversed
	int                penmm,     // penalty to incur for color miscall
	const Penalties&   pens,      // Penalties
	size_t             nceil,     // max # Ns allowed
	int                penceil)   // penalty ceiling
{
	assert(uc.finalized);
	if(uc.empty) return;
	bool ninvolved = (c > 3 || refMask > 15);
	for(int to = 0; to < 4; to++) {
		// Assuming that the read character in this row is 'to'...
		int add = ((matches(to, refMask) == 1) ? 0 : ((refMask > 15) ? pens.n(30) : pens.snp));
		AlignmentScore from[] = {
			uc.best[0] - add, uc.best[1] - add,
			uc.best[2] - add, uc.best[3] - add };
		assert(!VALID_AL_SCORE(from[0]) || from[0].score_ <= 0);
		assert(!VALID_AL_SCORE(from[1]) || from[1].score_ <= 0);
		assert(!VALID_AL_SCORE(from[2]) || from[2].score_ <= 0);
		assert(!VALID_AL_SCORE(from[3]) || from[3].score_ <= 0);
		// Calculate which 'from' character isn't going to
		// cause a color mismatch
		const int goodfrom = nuccol2nuc[to][c];
		if(goodfrom < 4) {
			from[goodfrom].score_ += penmm;
		}
		for(int fr = 0; fr < 4; fr++) {
			AlignmentScore frBest = from[fr];
			AlignmentScore& myBest = best[to];
			SwColorCellMask& myMask = mask[to];
			if(ninvolved) frBest.incNs(nceil);
			if(!VALID_AL_SCORE(frBest)) continue;
			frBest.score_ -= penmm;
			assert_leq(frBest.score_, 0);
			if(-frBest.score_ <= penceil && frBest >= myBest) {
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
 * If an alignment is found, its offset relative to rdi is returned.
 * E.g. if an alignment is found that occurs starting at rdi, 0 is
 * returned.  If no alignment is found, -1 is returned.
 */
int SwAligner::alignColors(
	const BTDnaString& rd, // read color sequence
	const BTString& qu,    // read qualities
	size_t rdi,            // offset of first character within 'read' to consider
	size_t rdf,            // offset of last char (exclusive) in 'read' to consider
	BTString& rf,          // reference sequence, as masks
	size_t rfi,            // offset of first character within 'rf' to consider
	size_t rff,            // offset of last char (exclusive) in 'rf' to consider
	const SwParams& pa,    // parameters governing Smith-Waterman problem
	const Penalties& pen,  // penalties for various edits
	int penceil,           // penalty ceiling for valid alignments
	int& nup,              // upstream decoded nucleotide
	int& ndn,              // downstream decoded nucleotide
	SwResult& res,         // edits and scores
	RandomSource& rnd)     // pseudo-random generator
{
	typedef SwColorCell TCell;
	assert_leq(rdf, rd.length());
	assert_leq(rdf, qu.length());
	assert_leq(rff, rf.length());
	assert_geq(rf.length(), rd.length()+1);
	assert_lt(rfi, rff);
	assert_lt(rdi, rdf);
	assert_eq(rd.length(), qu.length());
	assert_geq(pa.gapBar, 1);
	res.sws++;
#ifndef NDEBUG
	for(size_t i = rfi; i < rff; i++) {
		assert_range(0, 16, (int)rf[i]);
	}
#endif
	
	// Calculate the largest possible number of read and reference gaps
	// given 'penceil' and 'pens'
	int readGaps = pen.maxReadGaps(penceil);
	int refGaps  = pen.maxRefGaps(penceil);
	assert_geq(readGaps, 0);
	assert_geq(refGaps, 0);
	size_t nceil = pen.nCeil(rd.length());

	//
	// Initialize the first row
	//
	
	ELList<TCell>& tab = ctab_;
	tab.resize(1); // add first row to row list
	int maxGaps = max(readGaps, refGaps);
	const int wlo = 0;
	const int whi = maxGaps * 2;
	tab[0].resize(whi-wlo+1); // add columns to first row
	bool validInRow = false;
	// Calculate starting values for the rest of the columns in the
	// first row.  No need to consider colors yet since we won't have
	// consecutive pairs of read nucleotides until we get to the second
	// row.  We just consider reference mutations.
	for(int col = 0; col <= whi; col++) {
		tab[0][col].clear(); // clear the cell; masks and scores
		int fromEnd = whi - col;
		int rfm = rf[rfi+col];
		// Can we start from here?
		if(col >= refGaps - readGaps && fromEnd >= readGaps - refGaps) {
			for(int to = 0; to < 4; to++) {
				int m = matches(to, rfm);
				if(m == 1) {
					// The assigned subject nucleotide matches the reference;
					// no penalty
					assert_lt(rfm, 16);
					tab[0][col].best[to].gaps_ = 0;
					tab[0][col].best[to].ns_ = 0;
					tab[0][col].best[to].score_ = 0;
					tab[0][col].mask[to].diag = 0xf;
				} else if(m == 0 && pen.snp <= penceil) {
					// The assigned subject nucleotide does not match the
					// reference nucleotide, so we add a SNP penalty
					tab[0][col].best[to].gaps_ = 0;
					tab[0][col].best[to].ns_ = 0;
					tab[0][col].best[to].score_ = -pen.snp;
					tab[0][col].mask[to].diag = 0xf;
				} else if(m == -1) {
					// The assigned subject nucleotide does not match the
					// reference nucleotide, so we add a SNP penalty
					tab[0][col].best[to].gaps_ = 0;
					tab[0][col].best[to].ns_ = 1;
					tab[0][col].best[to].score_ = -pen.n(30);
					tab[0][col].mask[to].diag = 0xf;
				} else {
					// Leave mask[to] cleared
				}
			}
		}
		// Calculate horizontals if barrier allows
		if(pa.gapBar <= 1 && col > 0) {
			tab[0][col].updateHoriz(tab[0][col-1], rfm, pen, nceil, penceil);
			res.swcups++;
		}
		assert(!tab[0][col].finalized);
		if(tab[0][col].finalize(penceil)) validInRow = true;
	}
	res.swrows++;
	if(!validInRow) {
		res.swskiprows += (rdf - rdi);
		assert(res.empty());
		return -1;
	}

	//
	// Calculate all subsequent rows
	//

	// Do rest of table
	for(int row = 1; row <= (int)(rdf-rdi); row++) {
		res.swrows++;
		tab.expand(); // add another row
		bool onlyDiagInto = (row+1 <= pa.gapBar || (int)(rdf-rdi)-row <= pa.gapBar);
		tab.back().resize(whi-wlo+1); // add enough space for columns
		assert_range(1, (int)qu.length(), row);
		assert_range(1, (int)rd.length(), row);
		int c = rd[row-1];   // read character in this row
		int q = QUAL(row-1); // quality for the read character; should be Phred
		assert_geq(q, 0);
		const int mmpen = ((c > 3) ? pen.n(30) : pen.mm(q));
		validInRow = false;		
		//
		// Handle col == wlo case before (and the col == whi case
		// after) the upcoming inner loop to reduce the number of
		// guards needed inside it.
		//
		assert_leq(wlo, whi);
		int col = wlo;
		TCell& cur = tab[row][0];
		cur.clear();
		if(!tab[row-1][0].empty) {
			const int fullcol = col + row;
			cur.updateDiag(
				tab[row-1][0],     // cell diagonally above and to the left
				rf[rfi + fullcol], // ref mask associated with destination cell
				c,                 // color being traversed
				mmpen,             // penalty to incur for color miscall
				pen,               // Penalties
				nceil,             // max # Ns allowed
				penceil);          // max penalty allowed
		}
		if(!onlyDiagInto && col < whi && !tab[row-1][1].empty) {
			cur.updateVert(
				tab[row-1][1],     // cell above
				c,                 // color being traversed
				mmpen,             // penalty to incur for color miscall
				pen,               // Penalties
				nceil,             // max # Ns allowed
				penceil);          // max penalty allowed
		}
		res.swcups++;
		// 'cur' is now initialized
		assert(!cur.finalized);
		if(cur.finalize(penceil)) validInRow = true;

		// Iterate from leftmost to rightmost inner diagonals
		for(col = wlo+1; col < whi; col++) {
			const int fullcol = col + row;
			int r = rf[rfi + fullcol];
			TCell& cur = tab[row][col-wlo];
			cur.clear();
			TCell& dg = tab[row-1][col-wlo];
			// The mismatch penalty is a function of the read character
			// in this row & the reference character in this column
			// (specifically: whether they match and whether either is
			// an N) as well as the quality value of the read
			// character.
			cur.updateDiag(
				dg,            // cell diagonally above and to the left
				r,             // ref mask associated with destination column
				c,             // nucleotide in destination row
				mmpen,         // penalty to incur for color miscall
				pen,           // Penalties
				nceil,         // max # Ns allowed
				penceil);      // max penalty allowed
			if(!onlyDiagInto) {
				TCell& up = tab[row-1][col-wlo+1];
				cur.updateVert(
					up,        // cell above
					c,         // nucleotide in destination row
					mmpen,     // penalty to incur for color miscall
					pen,       // Penalties
					nceil,     // max # Ns allowed
					penceil);  // max penalty allowed
				// Can do horizontal
				TCell& lf = tab[row][col-wlo-1];
				cur.updateHoriz(
					lf,        // cell to the left
					r,         // ref mask associated with destination column
					pen,       // Penalties
					nceil,     // max # Ns allowed
					penceil);  // max penalty allowed
			}
			res.swcups++;
			// 'cur' is now initialized
			assert(!cur.finalized);
			if(cur.finalize(penceil)) validInRow = true;
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
			const int fullcol = col + row;
			const int r = rf[rfi + fullcol];
			TCell& dg = tab[row-1][col-wlo];
			if(!dg.empty) {
				cur.updateDiag(
					dg,        // cell diagonally above and to the left
					r,         // ref mask associated with destination column
					c,         // nucleotide in destination row
					mmpen,     // penalty to incur for color miscall
					pen,       // Penalties
					nceil,     // max # Ns allowed
					penceil);  // max penalty allowed
			}
			TCell& lf = tab[row][col-wlo-1];
			if(!onlyDiagInto && !lf.empty) {
				cur.updateHoriz(
					lf,        // cell to the left
					r,         // ref mask associated with destination column
					pen,       // Penalties
					nceil,     // max # Ns allowed
					penceil);  // max penalty allowed
			}
			res.swcups++;
			// 'cur' is now initialized
			assert(!cur.finalized);
			if(cur.finalize(penceil)) validInRow = true;
		}
		if(!validInRow) {
			assert_geq((int)(rdf-rdi), row);
			res.swskiprows += (rdf - rdi - row);
			assert(res.empty());
			return -1;
		}
	}
	assert_eq(tab.size(), rd.length()+1);
	// Go hunting for cell to backtrace from; i.e. best score in the
	// bottom row and in the last readGaps*2+1 columns.
	AlignmentScore bestScore = AlignmentScore::INVALID();
	int btCol = -1; // column to backtrace from
	int btC = -1;
	int lastRow = rdf-rdi;
	for(int col = wlo; col <= whi; col++) {
		// Can we backtrace from this cell?  Depends on gaps.
		int fromEnd = whi - col;
		// greater than or equal to???
		if(fromEnd >= refGaps - readGaps && col >= readGaps - refGaps) {
			if(!tab[lastRow][col].empty) {
				assert(tab[lastRow][col].finalized);
				if(tab[lastRow][col].updateBest(bestScore, btC, penceil)) {
					assert(!VALID_AL_SCORE(bestScore) || abs(bestScore.score()) <= penceil);
					btCol = col;
				}
			}
		}
	}
	assert_range(wlo, whi, btCol);
	assert_range(0, 3, btC);
	assert_leq(abs(bestScore.score()), penceil);
	int off = backtrackColors(
		rd,        // read sequence
		qu,        // read qualities
		rdi,       // offset of first character within 'rd' to consider
		rdf,       // offset of last char (exclusive) in 'rd' to consider
		rf,        // reference sequence, as masks
		rfi,       // offset of first character within 'rf' to consider
		rff,       // offset of last char (exclusive) in 'rf' to consider
		bestScore, // best score of any cell that corresponded to a totally-aligned read
		readGaps,  // maximum number of read gaps allowed
		refGaps,   // maximum number of reference gaps allowed
		pa,        // parameters governing Smith-Waterman problem
		pen,       // penalties for various edits
		nup,       // upstream decoded nucleotide
		ndn,       // downstream decoded nucleotide
		res,       // destination for best alignment
		btCol+lastRow, // col
		btC,       // randomly selected char from best cell
		rnd);      // pseudo-random generator
	assert_geq(off, 0);

#ifndef NDEBUG
	for(int i = 0; i < (int)res.alres.ced().size(); i++) {
		assert_lt(res.alres.ced()[i].pos, rdf-rdi);
	}
#endif
	return off;
}

