/*
 *  aligner_sw.cpp
 */

#include "aligner_sw.h"
#include "aligner_sw_col.h"
#include "aligner_result.h"
#include "search_globals.h"
#include "scoring.h"
#include "mask.h"

/**
 * We have set all three of the the cell's intermediate scores; set the 'empty'
 * and 'finalized' fields appropriately.
 */
inline bool SwNucCell::finalize(TAlScore floorsc) {
	ASSERT_ONLY(finalized = true);
	assert(empty);
	assert(repOk());
	// Profiling shows cache misses on following line
	bool aboveFloor = oallBest.valid() && oallBest.score() >= floorsc;
	empty = true;
	if(!mask.empty() && aboveFloor) {
		assert(VALID_AL_SCORE(oallBest));
		empty = false;
	} else if(aboveFloor) {
		assert(VALID_AL_SCORE(oallBest));
		terminal = true;
	}
	return !empty;
}

/**
 * Check that cell is internally consistent
 */
bool SwNucCell::repOk() const {
	assert(oallBest >= rdgapBest || !rdgapBest.valid());
	assert(oallBest >= rfgapBest || !rfgapBest.valid());
	return true;
}

/**
 * Given a read, an alignment orientation, a range of characters in a referece
 * sequence, and a bit-encoded version of the reference, set up and execute the
 * corresponding dynamic programming problem.
 *
 * The caller has already narrowed down the relevant portion of the reference
 * using, e.g., the location of a seed hit, or the range of possible fragment
 * lengths if we're searching for the opposite mate in a pair.
 */
void SwAligner::init(
	const Read& rd,        // read to align
	size_t rdi,            // off of first char in 'rd' to consider
	size_t rdf,            // off of last char (exclusive) in 'rd' to consider
	bool fw,               // whether to align forward or revcomp read
	bool color,            // colorspace?
	uint32_t refidx,       // reference aligned against
	int64_t rfi,           // first reference base to SW align against
	int64_t rff,           // last ref base (exclusive) to SW align against
	const BitPairReference& refs, // Reference strings
	size_t reflen,         // length of reference sequence
	size_t width,          // # bands to do (width of parallelogram)
	EList<bool>* st,       // mask indicating which columns we can start in
	EList<bool>* en,       // mask indicating which columns we can end in
	const Scoring& pen,    // penalty scheme
	TAlScore minsc,        // minimum score for a valid alignment
	TAlScore floorsc,      // local-alignment score floor
	int nceil)             // max # Ns allowed in reference portion of aln
{
	assert_gt(rff, rfi);
	assert_gt(rdf, rdi);
	assert(color == rd.color);
	// Figure the number of Ns we're going to add to either side
	size_t leftNs  =
		(rfi >= 0               ? 0 : (size_t)std::abs(rfi));
	size_t rightNs =
		(rff <= (int64_t)reflen ? 0 : (size_t)std::abs(rff - (int64_t)reflen));
	// rflen = full length of the reference substring to consider, including
	// overhang off the boundaries of the reference sequence
	const size_t rflen = (size_t)(rff - rfi);
	// rflenInner = length of just the portion that doesn't overhang ref ends
	const size_t rflenInner = rflen - (leftNs + rightNs);
#ifndef NDEBUG
	bool haveRfbuf2 = false;
	EList<char> rfbuf2(rflen);
	// This is really slow, so only do it some of the time
	if((rand() % 10) == 0) {
		int64_t rfii = rfi;
		for(size_t i = 0; i < rflen; i++) {
			if(rfii < 0 || (size_t)rfii >= reflen) {
				rfbuf2.push_back(4);
			} else {
				rfbuf2.push_back(refs.getBase(refidx, (uint32_t)rfii));
			}
			rfii++;
		}
		haveRfbuf2 = true;
	}
#endif
	// rfbuf_ = uint32_t list large enough to accommodate both the reference
	// sequence and any Ns we might add to either side.
	rfwbuf_.resize((rflen + 16) / 4);
	int offset = refs.getStretch(
		rfwbuf_.ptr(),               // buffer to store words in
		refidx,                      // which reference
		(rfi < 0) ? 0 : (size_t)rfi, // starting offset (can't be < 0)
		rflenInner);                 // length to grab (exclude overhang)
	assert_leq(offset, 16);
	rf_ = (char*)rfwbuf_.ptr() + offset;
	// Shift ref chars away from 0 so we can stick Ns at the beginning
	if(leftNs > 0) {
		// Slide everyone down
		for(size_t i = rflenInner; i > 0; i--) {
			rf_[i+leftNs-1] = rf_[i-1];
		}
		// Add Ns
		for(size_t i = 0; i < leftNs; i++) {
			rf_[i] = 4;
		}
	}
	if(rightNs > 0) {
		// Add Ns to the end
		for(size_t i = 0; i < rightNs; i++) {
			rf_[i + leftNs + rflenInner] = 4;
		}
	}
	// Count Ns and convert reference characters into A/C/G/T masks.  Ambiguous
	// nucleotides (IUPAC codes) have more than one mask bit set.
	for(size_t i = 0; i < rflen; i++) {
		assert(!haveRfbuf2 || rf_[i] == rfbuf2[i]);
		// Make it into a mask
		assert_range(0, 4, (int)rf_[i]);
		// rf_[i] gets mask version of refence char, with N=16
		rf_[i] = (1 << rf_[i]);
	}
	RandomSource rnd(rd.seed);
	const BTDnaString& rdseq  = fw ? rd.patFw : rd.patRc;
	const BTString&    rdqual = fw ? rd.qual  : rd.qualRev;
	init(
		rdseq,       // read sequence
		rdqual,      // read qualities
		rdi,         // offset of first char in read to consdier
		rdf,         // offset of last char (exclusive) in read to consdier
		fw,          // true iff read sequence is original fw read
		color,       // true iff read is colorspace
		refidx,      // id of reference aligned against
		rfi,         // offset of upstream ref char aligned against
		rf_,         // reference sequence, wrapped up in BTString object
		0,           // use the whole thing
		rflen,       // ditto
		width,       // # bands to do (width of parallelogram)
		st,          // mask indicating which columns we can start in
		en,          // mask indicating which columns we can end in
		pen,         // scoring scheme
		minsc,       // minimum score for valid alignment
		floorsc);    // local-alignment score floor
	filter(nceil);
}

/**
 * Align read 'rd' to reference using read & reference information given
 * last time init() was called.  If the read is colorspace, the decoding is
 * determined simultaneously with alignment.  Uses dynamic programming.
 */
bool SwAligner::align(RandomSource& rnd) {
	assert(inited());
	assert_eq(STATE_INITED, state_);
	bool possible;
	nfills_++;
	state_ = STATE_ALIGNED;
	if(color_) {
		possible = alignColors(rnd);
	} else {
		possible = alignNucleotides(rnd);
	}
	if(VALID_SCORE(solbest_)) {
		nsucc_++;
	} else {
		nfail_++;
	}
	assert(repOk());
	cural_ = 0;
	return possible;
}

/**
 * Select a path for backtracking from the oall table version of this cell.
 * If there is a tie among eligible paths, break it randomly.  Return value
 * is a flag indicating the backtrack type (see enum defining SW_BT_*
 * above).
 */
int SwNucCellMask::randOverallBacktrack(
	RandomSource& rand,
	bool& branch,
	bool clear)
{
	ASSERT_ONLY(int num = numOverallPossible());
	int i = ((oall_diag != 0) << 0) |
			((oall_rfop != 0) << 1) |
			((oall_rfex != 0) << 2) |
			((oall_rdop != 0) << 3) |
			((oall_rdex != 0) << 4);
	int ret = randFromMask(rand, i) + SW_BT_OALL_DIAG;
	// If we're choosing from among >1 possibilities, inform caller so that
	// caller can add a frame to the backtrack stack.
	branch = alts5[i] > 1;
	assert_range((int)SW_BT_OALL_DIAG, (int)SW_BT_OALL_READ_EXTEND, ret);
	// Clear the bit associated with the path chosen
	if(clear) {
		switch(ret) {
			case SW_BT_OALL_DIAG:        oall_diag = 0; break;
			case SW_BT_OALL_REF_OPEN:    oall_rfop = 0; break;
			case SW_BT_OALL_REF_EXTEND:  oall_rfex = 0; break;
			case SW_BT_OALL_READ_OPEN:   oall_rdop = 0; break;
			case SW_BT_OALL_READ_EXTEND: oall_rdex = 0; break;
			default: throw 1; break;
		}
	}
	assert_eq(num - (clear ? 1 : 0), numOverallPossible());
	return ret;
}

/**
 * Select a path for backtracking from the rdgap table version of this cell.
 * If there is a tie among eligible paths, break it randomly.  Return value
 * is a flag indicating the backtrack type (see enum defining SW_BT_*
 * above).
 */
int SwNucCellMask::randReadGapBacktrack(
	RandomSource& rand,
	bool& branch,
	bool clear)
{
	ASSERT_ONLY(int num = numReadGapPossible());
	int i = ((rdgap_op != 0) << 0) |
			((rdgap_ex != 0) << 1);
	int ret = randFromMask(rand, i) + SW_BT_RDGAP_OPEN;
	// If we're choosing from among >1 possibilities, inform caller so that
	// caller can add a frame to the backtrack stack.
	branch = alts5[i] > 1;
	assert(ret == SW_BT_RDGAP_OPEN || ret == SW_BT_RDGAP_EXTEND);
	// Clear the bit associated with the path chosen
	if(clear) {
		switch(ret) {
			case SW_BT_RDGAP_OPEN:   rdgap_op = 0; break;
			case SW_BT_RDGAP_EXTEND: rdgap_ex = 0; break;
			default: throw 1; break;
		}
	}
	assert_eq(num - (clear ? 1 : 0), numReadGapPossible());
	return ret;
}

/**
 * Select a path for backtracking from the rfgap table version of this cell.
 * If there is a tie among eligible paths, break it randomly.  Return value
 * is a flag indicating the backtrack type (see enum defining SW_BT_*
 * above).
 */
int SwNucCellMask::randRefGapBacktrack(
	RandomSource& rand,
	bool& branch,
	bool clear)
{
	ASSERT_ONLY(int num = numRefGapPossible());
	int i = ((rfgap_op != 0) << 0) |
			((rfgap_ex != 0) << 1);
	int ret = randFromMask(rand, i) + SW_BT_RFGAP_OPEN;
	// If we're choosing from among >1 possibilities, inform caller so that
	// caller can add a frame to the backtrack stack.
	branch = alts5[i] > 1;
	assert(ret == SW_BT_RFGAP_OPEN || ret == SW_BT_RFGAP_EXTEND);
	// Clear the bit associated with the path chosen
	if(clear) {
		switch(ret) {
			case SW_BT_RFGAP_OPEN:   rfgap_op = 0; break;
			case SW_BT_RFGAP_EXTEND: rfgap_ex = 0; break;
			default: throw 1; break;
		}
	}
	assert_eq(num - (clear ? 1 : 0), numRefGapPossible());
	return ret;
}

/**
 * Given the dynamic programming table and a cell, trace backwards from the
 * cell and install the edits and score/penalty in the appropriate fields
 * of res.  The RandomSource is used to break ties among equally good ways
 * of tracing back.
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
bool SwAligner::backtrackNucleotides(
	TAlScore       escore, // score we expect to get over backtrack
	SwResult&      res,    // out: store results (edits and scores) here
	size_t&        off,    // out: store results (edits and scores) here
	size_t         row,    // start in this parallelogram row
	size_t         col,    // start in this parallelogram column
	RandomSource&  rand)   // pseudo-random generator
{
	typedef SwNucCell TCell;
	ELList<TCell>& tab = ntab_;
	assert_lt(row, rd_->length());
	btnstack_.clear();
	btcells_.clear();
	size_t tabcol = col - row;
	AlnScore score; score.score_ = 0;
	score.gaps_ = score.ns_ = 0;
	size_t origCol = col;
	size_t gaps = 0, readGaps = 0, refGaps = 0;
	res.alres.reset();
	EList<Edit>& ned = res.alres.ned();
	assert(ned.empty());
	assert_gt(dpRows(), row);
	size_t trimEnd = dpRows() - row - 1; 
	size_t trimBeg = 0;
	assert(!sc_->monotone || escore <= 0);
	assert(!sc_->monotone || score.score() >= escore);
	int ct = SW_BT_CELL_OALL; // cell type
	while((int)row >= 0) {
		nbts_++;
		assert_lt(row, tab.size());
		assert_leq(col, origCol);
		assert_geq(col, row);
		tabcol = col - row;
		assert_geq(tabcol, 0);
		TCell& curc = tab[row][tabcol];
		bool empty = curc.mask.numPossible(ct) == 0;
		assert_eq(gaps, Edit::numGaps(ned));
		assert_leq(gaps, rdgap_ + rfgap_);
		// Cell was involved in a previously-reported alignment?
		if(!curc.canMoveThrough(ct)) {
			if(!btnstack_.empty()) {
				// Remove all the cells from list back to and including the
				// cell where the branch occurred
				btcells_.resize(btnstack_.back().celsz);
				// Pop record off the top of the stack
				ned.resize(btnstack_.back().nedsz);
				//aed.resize(btnstack_.back().aedsz);
				row      = btnstack_.back().row;
				col      = btnstack_.back().col;
				gaps     = btnstack_.back().gaps;
				readGaps = btnstack_.back().readGaps;
				refGaps  = btnstack_.back().refGaps;
				score    = btnstack_.back().score;
				ct       = btnstack_.back().ct;
				btnstack_.pop_back();
				assert(!sc_->monotone || score.score() >= escore);
				continue;
			} else {
				// No branch points to revisit; just give up
				return false;
			}
		}
		assert(!curc.reportedThru_);
		assert(!sc_->monotone || score.score() >= minsc_);
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
		int cur = curc.mask.randBacktrack(ct, rand, branch, true);
		if(branch) {
			assert_gt(curc.mask.numPossible(ct), 0);
			// Add a frame to the backtrack stack
			btnstack_.expand();
			btnstack_.back().init(
				ned.size(),
				0,               // aed.size()
				btcells_.size(),
				row,
				col,
				gaps,
				readGaps,
				refGaps,
				score,
				ct);
		}
		btcells_.expand();
		btcells_.back().first = row;
		btcells_.back().second = tabcol;
		SizeTPair p = make_pair(row, rfi_ + tabcol);
		switch(cur) {
			// Move up and to the left.  If the reference nucleotide in the
			// source row mismatches the read nucleotide, penalize
			// it and add a nucleotide mismatch.
			case SW_BT_OALL_DIAG: {
				assert_gt(row, 0); assert_gt(col, 0);
				// Check for color mismatch
				int readC = (*rd_)[row];
				int refNmask = (int)rf_[rfi_+col];
				assert_gt(refNmask, 0);
				int m = matchesEx(readC, refNmask);
				ASSERT_ONLY(int lastct = ct);
				ct = SW_BT_CELL_OALL;
				if(m != 1) {
					Edit e(
						(int)row,
						mask2dna[refNmask],
						"ACGTN"[readC],
						EDIT_TYPE_MM);
					assert(e.repOk());
					assert(ned.empty() || ned.back().pos >= row);
					ned.push_back(e);
					int pen = QUAL2(row, col);
					assert_geq(tab[row-1][col-row].best(ct).score(),
					           tab[row][col-row].best(lastct).score());
					assert_eq(pen, tab[row-1][col-row].best(ct).score() -
					               tab[row][col-row].best(lastct).score());
					score.score_ -= pen;
					assert(!sc_->monotone || score.score() >= escore);
				} else {
					// Reward a match
					int64_t bonus = sc_->match(30);
					assert_geq(tab[row][col-row].best(lastct).score(),
					           tab[row-1][col-row].best(ct).score());
					assert_eq(bonus, tab[row][col-row].best(lastct).score() -
					                 tab[row-1][col-row].best(ct).score());
					score.score_ += bonus;
					assert(!sc_->monotone || score.score() >= escore);
				}
				if(m == -1) {
					score.ns_++;
				}
				row--; col--;
				assert_lt(row, tab.size());
				assert(VALID_AL_SCORE(score));
				break;
			}
			// Move up.  Add an edit encoding the ref gap.
			case SW_BT_OALL_REF_OPEN:
			case SW_BT_RFGAP_OPEN:
			{
				assert_gt(row, 0);
				Edit e(
					(int)row,
					'-',
					"ACGTN"[(int)(*rd_)[row]],
					EDIT_TYPE_REF_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				row--;
				ASSERT_ONLY(int lastct = ct);
				ct = SW_BT_CELL_OALL;
				int pen = sc_->refGapOpen();
				assert_geq(tab[row][col-row].best(ct).score(),
				           tab[row+1][col-(row+1)].best(lastct).score());
				assert_eq(pen, tab[row][col-row].best(ct).score() -
				               tab[row+1][col-(row+1)].best(lastct).score());
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; refGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				assert_lt(row, tab.size());
				assert_lt(col - row, tab[row].size());
				break;
			}
			// Move up.  Add an edit encoding the ref gap.
			case SW_BT_OALL_REF_EXTEND:
			case SW_BT_RFGAP_EXTEND:
			{
				assert_gt(row, 1);
				Edit e(
					(int)row,
					'-',
					"ACGTN"[(int)(*rd_)[row]],
					EDIT_TYPE_REF_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				row--;
				ASSERT_ONLY(int lastct = ct);
				ct = SW_BT_CELL_RFGAP;
				int pen = sc_->refGapExtend();
				assert_geq(tab[row][col-row].best(ct).score(),
				           tab[row+1][col-(row+1)].best(lastct).score());
				assert_eq(pen, tab[row][col-row].best(ct).score() -
				               tab[row+1][col-(row+1)].best(lastct).score());
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; refGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				assert_lt(row, tab.size());
				tabcol = col - row;
				assert_lt(col - row, tab[row].size()-1);
				break;
			}
			case SW_BT_OALL_READ_OPEN:
			case SW_BT_RDGAP_OPEN:
			{
				assert_gt(col, 0);
				Edit e(
					(int)row+1,
					mask2dna[(int)rf_[rfi_+col]],
					'-',
					EDIT_TYPE_READ_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				col--;
				ASSERT_ONLY(int lastct = ct);
				ct = SW_BT_CELL_OALL;
				int pen = sc_->readGapOpen();
				assert_geq(tab[row][col-row].best(ct).score(),
				           tab[row][col-row+1].best(lastct).score());
				assert_eq(pen, tab[row][col-row].best(ct).score() -
				               tab[row][col-row+1].best(lastct).score());
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; readGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				assert_lt(row, tab.size());
				tabcol = col - row;
				assert_lt(col - row, tab[row].size());
				break;
			}
			case SW_BT_OALL_READ_EXTEND:
			case SW_BT_RDGAP_EXTEND:
			{
				assert_gt(col, 1);
				Edit e(
					(int)row+1,
					mask2dna[(int)rf_[rfi_+col]],
					'-',
					EDIT_TYPE_READ_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				col--;
				ASSERT_ONLY(int lastct = ct);
				ct = SW_BT_CELL_RDGAP;
				int pen = sc_->readGapExtend();
				assert_geq(tab[row][col-row].best(ct).score(),
				           tab[row][col-row+1].best(lastct).score());
				assert_eq(pen, tab[row][col-row].best(ct).score() -
				               tab[row][col-row+1].best(lastct).score());
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; readGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				assert_lt(row, tab.size());
				tabcol = col - row;
				assert_lt(col - row, tab[row].size());
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
	// Set 'reported' flag on each cell
	for(size_t i = 0; i < btcells_.size(); i++) {
		size_t rw = btcells_[i].first;
		size_t cl = btcells_[i].second;
		assert(!tab[rw][cl].reportedThru_);
		tab[rw][cl].setReportedThrough();
	}
	int readC = (*rd_)[rdi_+row];      // get last char in read
	int refNmask = (int)rf_[rfi_+col]; // get last ref char ref involved in aln
	assert_gt(refNmask, 0);
	int m = matchesEx(readC, refNmask);
	if(m != 1) {
		Edit e((int)row, mask2dna[refNmask], "ACGTN"[readC], EDIT_TYPE_MM);
		assert(e.repOk());
		assert(ned.empty() || ned.back().pos >= row);
		ned.push_back(e);
		score.score_ -= QUAL2(row, col);
		assert_geq(score.score(), minsc_);
	} else {
		score.score_ += sc_->match(30);
	}
	if(m == -1) {
		score.ns_++;
	}
	res.reverse();
	assert(Edit::repOk(ned, (*rd_)));
#ifndef NDEBUG
	size_t gapsCheck = 0;
	for(size_t i = 0; i < ned.size(); i++) {
		if(ned[i].isGap()) gapsCheck++;
	}
	assert_eq(gaps, gapsCheck);
	BTDnaString refstr;
	for(size_t i = col; i <= origCol; i++) {
		refstr.append(firsts5[(int)rf_[rfi_+i]]);
	}
	BTDnaString editstr;
	Edit::toRef((*rd_), ned, editstr, trimBeg, trimEnd);
	if(refstr != editstr) {
		cerr << "Decoded nucleotides and edits don't match reference:" << endl;
		cerr << "           score: " << score.score()
		     << " (" << gaps << " gaps)" << endl;
		cerr << "           edits: ";
		Edit::print(cerr, ned);
		cerr << endl;
		cerr << "    decoded nucs: " << (*rd_) << endl;
		cerr << "     edited nucs: " << editstr << endl;
		cerr << "  reference nucs: " << refstr << endl;
		assert(0);
	}
#endif
	assert_eq(score.score(), escore);
	assert_leq(gaps, rdgap_ + rfgap_);
	off = (col - row);
	assert_lt(col + rfi_, rff_);
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
		true,                     // alignment trim soft?
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
 * Update a cell's oallBest, rdgapBest, and mask with respect to its neighbor
 * on the left.
 */
inline void SwAligner::updateNucHoriz(
	const SwNucCell& lc,
	SwNucCell&       dstc,
	int              rfm)
{
	assert(lc.finalized);
	if(lc.empty) return;
	IF_COUNT_N(bool ninvolved = rfm > 15);
	{
		AlnScore leftOallBest  = lc.oallBest;
		AlnScore leftRdgapBest = lc.rdgapBest;
		AlnScore& myOallBest   = dstc.oallBest;
		AlnScore& myRdgapBest  = dstc.rdgapBest;
		SwNucCellMask& myMask  = dstc.mask;
		assert(!myRdgapBest.valid());
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
		assert(dstc.repOk());
	}
}

/**
 * Update a SwCell's best[] and mask[] arrays with respect to its neighbor on
 * the left.
 */
inline void SwAligner::updateNucVert(
	const SwNucCell& uc,
	SwNucCell&       dstc,
	int              rdc)
{
	assert(uc.finalized);
	if(uc.empty) return;
	IF_COUNT_N(bool ninvolved = rdc > 3);
	{
		{
			AlnScore  upOallBest  = uc.oallBest;
			AlnScore  upRfgapBest = uc.rfgapBest;
			AlnScore& myOallBest  = dstc.oallBest;
			AlnScore& myRfgapBest = dstc.rfgapBest;
			SwNucCellMask& myMask = dstc.mask;
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
					myMask.oall_rfex = 1;
				}
				if(ex >= myRfgapBest) {
					if(ex > myRfgapBest) {
						myMask.clearRefGapMask();
						myRfgapBest = ex;
					}
					myMask.rfgap_ex = 1;
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
					myMask.oall_rfop = 1;
				}
				if(ex >= myRfgapBest) {
					if(ex > myRfgapBest) {
						myMask.clearRefGapMask();
						myRfgapBest = ex;
					}
					myMask.rfgap_op = 1;
				}
			}
			assert(dstc.repOk());
		}
	}
}

/**
 * Update a SwCell's best[] and mask[] arrays with respect to its
 * neighbor up and to the left.  SNPs are charged
 */
inline void SwAligner::updateNucDiag(
	const SwNucCell& dc,
	SwNucCell&       dstc,
	int              rdc,
	int              rfm,
	int              qpen,
	bool&            improved)
{
	assert(dc.finalized);
	if(dc.empty) return;
	bool match = matches(rdc, rfm);
	TAlScore add;
	if(match) {
		add = sc_->match(30);
		improved = true;
	} else {
		add = -qpen;
	}
	{
		{
			AlnScore& myOallBest  = dstc.oallBest;
			SwNucCellMask& myMask = dstc.mask;
			assert(!sc_->monotone || dc.oallBest.score() <= 0);
			if(!VALID_AL_SCORE(dc.oallBest)) {
				return;
			}
			COUNT_N(dc.oallBest);
			AlnScore ex = dc.oallBest + (int)add;
			if(ex.score_ >= floorsc_ && ex >= myOallBest) {
				if(ex > myOallBest) {
					myOallBest = ex;
					myMask.clear();
				}
				myMask.oall_diag = 1;
			}
			assert(dstc.repOk());
		}
	}
}

/**
 * Set elements of en_ to false if an ungapped alignment extending
 * diagonally back from the corresponding cell in the last row would
 * overlap too many Ns (more than nlim).
 */
size_t SwAligner::nfilter(size_t nlim) {
	size_t nrow = dpRows();
	size_t ns = 0;
	size_t filtered = 0;
	assert_gt(nrow, 0);
	assert_eq(rff_ - rfi_, width_ + nrow - 1);
	assert(en_ != NULL);
	assert_geq(rff_ - rfi_, nrow);
	assert_eq(rff_ - rfi_ - (nrow-1), en_->size());
	// For each reference character involved
	for(size_t i = rfi_; i < rff_; i++) {
		int rfc = rf_[i];
		assert_range(0, 16, rfc);
		size_t ii = i - rfi_;
		if(rfc == 16) {
			ns++;
		}
		if(ii >= (nrow-1)) {
			size_t col = ii - (nrow-1);
			if(ii >= nrow) {
				int prev_rfc = rf_[i - nrow];
				assert_range(0, 16, prev_rfc);
				if(prev_rfc == 16) {
					assert_gt(ns, 0);
					ns--;
				}
			}
			if(ns > nlim) {
				if((*en_)[col]) {
					(*en_)[col] = false;
					filtered++;
				}
			}
		}
	}
	return filtered;
}

// Assumes:
// 1. 'row' = 0-based row of DP table
// 2. 'col' = 0-based col of DP table
// 3. 'minsc_' = minimum AlignmentScore required for a cell to have a sol
// 4. 'solrowlo_' = if row >= this, cells in this row can have solutions
// 5. 'solrows_' is a pairs encoding the range of rows s.t. all rows with
//    solutions are in the range
// 6. 'solcols_' is a list of pairs, where each pair cor
// 7. 'solbest_'
// 8. 'solrowbest_'
#define UPDATE_SOLS(cur, row, col, improved) { \
	assert_lt(col, width_); \
	if(cur.oallBest.score() >= minsc_ && \
	   (en_ == NULL || (*en_)[col])) \
	{ \
		const bool local = !sc_->monotone; \
		/* For local alignment, a cell is only a solution candidate if */ \
		/* the score was improved when we moved into the cell. */ \
		if(!local || improved) { \
			/* Score is acceptable */ \
			if((int64_t)row >= solrowlo_) { \
				/* Both score and row are acceptable */ \
				if(row < solrows_.first)        solrows_.first = row; \
				if(row > solrows_.second)       solrows_.second = row; \
				if(col < solcols_[row].first) { \
					solcols_[row].first  = col; \
					assert(en_ == NULL || solcols_[row].first < en_->size()); \
				} \
				if(col > solcols_[row].second) { \
					solcols_[row].second = col; \
					assert(en_ == NULL || solcols_[row].second < en_->size()); \
				} \
				if(cur.oallBest.score() > solbest_) \
					solbest_ = cur.oallBest.score(); \
				assert(!sc_->monotone || solbest_ <= 0); \
				if(cur.oallBest.score() > solrowbest_[row]) \
					solrowbest_[row] = cur.oallBest.score(); \
				assert(!sc_->monotone || solrowbest_[row] <= 0); \
				nsols_++; \
				soldone_ = false; \
				assert(repOk()); \
			} \
		} \
	} \
}

#define FINALIZE_CELL(r, c, improved) { \
	assert_lt(c, width_); \
	assert(!tab[r][c].finalized); \
	ncups_++; \
	ASSERT_ONLY(tab[r][c].finalized = true); \
	assert(tab[r][c].empty); \
	assert(tab[r][c].repOk()); \
	bool aboveFloor = \
		tab[r][c].oallBest.valid() && \
		tab[r][c].oallBest.score() >= loBound; \
	if(!tab[r][c].mask.empty() && aboveFloor) { \
		assert(VALID_AL_SCORE(tab[r][c].oallBest)); \
		tab[r][c].empty = false; \
		UPDATE_SOLS(tab[r][c], r, c, improved); \
		assert(!tab[r][c].empty); \
		validInRow = true; \
	} else if(aboveFloor) { \
		assert(VALID_AL_SCORE(tab[r][c].oallBest)); \
		tab[r][c].terminal = true; \
	} \
}

/**
 * Align the nucleotide read in *rd_ against the reference string in rf_ using
 * a banded dynamic programming algorithm.  Depending on whether we're rounding
 * negative scores up to 0, and depending on whether we're giving a positive
 * score to matches, we might be doing local alignment, global alignment, or
 * something in between.
 *
 * ELList<SwNucCell> ntab_ is filled in with the scores.
 *
 * If an alignment is found, its offset relative to rdi is returned.
 * E.g. if an alignment is found that occurs starting at rdi, 0 is
 * returned.  If no alignment is found, -1 is returned.
 */
bool SwAligner::alignNucleotides(RandomSource& rnd) {
	typedef SwNucCell TCell;
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
	ELList<TCell>& tab = ntab_;
	tab.resize(1); // add first row to row list
	solrows_ = SwAligner::EXTREMES;        // at first, no rows have sols
	solcols_.resize(rdf_ - rdi_);
	solcols_.fill(SwAligner::EXTREMES);    // at first, no cols have sols
	solrowbest_.resize(rdf_ - rdi_);
	solrowbest_.fill(std::numeric_limits<TAlScore>::min()); // no row has sol
	const size_t wlo = 0;
	const size_t whi = (int)(width_ - 1);
	assert_lt(whi, rff_-rfi_);
	tab[0].resize(whi-wlo+1); // add columns to first row
	TAlScore loBound = sc_->monotone ? minsc_ : floorsc_;
	bool validInRow = !sc_->monotone;
	// Calculate starting values for the rest of the columns in the
	// first row.
	for(size_t col = 0; col <= whi; col++) {
		TCell& curc = tab[0][col];
		curc.clear(); // clear the cell; masks and scores
		assert(curc.repOk());
		int rdc = (*rd_)[rdi_+0];
		int rfm = rf_[rfi_+col];
		// Can we start from here?
		bool canStart = (st_ == NULL || (*st_)[col]);
		bool improved = false;
		if(canStart) {
			curc.oallBest.invalidate();
			curc.rdgapBest.invalidate();
			curc.rfgapBest.invalidate();
			int m = matchesEx(rdc, rfm);
			if(m == 1) {
				// The assigned subject nucleotide matches the reference;
				// no penalty
				assert_lt(rdc, 4);
				assert_lt(rfm, 16);
				curc.oallBest.score_ = sc_->match(30);
				curc.mask.oall_diag = 1;
				improved = true;
			} else if(m == 0 && QUAL2(0, col) <= -minsc_) {
				// Reference char mismatches
				curc.oallBest.score_ = -QUAL2(0, col);
				curc.mask.oall_diag = 1;
			} else if(m == -1) {
				int npen = sc_->n((int)(*qu_)[rdi_] - 33);
				if(npen <= -minsc_) {
					curc.oallBest.score_ = -npen;
					curc.oallBest.ns_ = 1;
					curc.mask.oall_diag = 1;
				}
			}
		}
		// Calculate horizontals if barrier allows
		if(sc_->gapbar < 1 && col > 0 && !tab[0][col-1].empty) {
#ifndef INLINE_CUPS
			updateNucHoriz(
				tab[0][col-1],
				curc,
				rfm);
#else
			assert(tab[0][col-1].finalized);
			IF_COUNT_N(bool ninvolved = rfm > 15);
			AlnScore leftOallBest  = tab[0][col-1].oallBest;
			AlnScore leftRdgapBest = tab[0][col-1].rdgapBest;
			AlnScore& myOallBest   = curc.oallBest;
			AlnScore& myRdgapBest  = curc.rdgapBest;
			SwNucCellMask& myMask  = curc.mask;
			assert(!myRdgapBest.valid());
			assert(!sc_->monotone || leftOallBest.score()  <= 0);
			assert(!sc_->monotone || leftRdgapBest.score() <= 0);
			assert_leq(leftRdgapBest, leftOallBest);
			if(VALID_AL_SCORE(leftOallBest)) {
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
				assert(curc.repOk());
			}
#endif
		}
		ncups_++;
		assert(tab[0][col].empty);
		assert(tab[0][col].repOk());
		assert(!tab[0][col].finalized);
		ASSERT_ONLY(tab[0][col].finalized = true);
		bool aboveFloor =
			tab[0][col].oallBest.valid() &&
			tab[0][col].oallBest.score() >= loBound;
		if(!tab[0][col].mask.empty() && aboveFloor) {
			assert(VALID_AL_SCORE(tab[0][col].oallBest));
			tab[0][col].empty = false;
			UPDATE_SOLS(tab[0][col], 0, col, improved);
			assert(!tab[0][col].empty);
			validInRow = true;
		} else if(aboveFloor) {
			assert(VALID_AL_SCORE(tab[0][col].oallBest));
			tab[0][col].terminal = true;
		}
	}
	nrowups_++;
	if(!validInRow) {
		nrowskips_ += (rdf_ - rdi_ - 1);
		return !soldone_;
	}

	//
	// Calculate all subsequent rows
	//

	// Do rest of table
	for(size_t row = 1; row < rdf_-rdi_; row++) {
		nrowups_++;
		tab.expand(); // add another row
		bool onlyDiagInto =
			(row+1 <= (size_t)sc_->gapbar ||
			 (int)(rdf_-rdi_)-row <= (size_t)sc_->gapbar);
		tab.back().resize(whi-wlo+1); // add enough space for columns
		assert_gt(row, 0);
		assert_lt(row, qu_->length());
		assert_lt(row, rd_->length());
		int c = (*rd_)[row];   // read character in this row
		validInRow = !sc_->monotone;
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
			// FIRST COLUMN OF A MIDDLE ROW - DIAGONAL UPDATE
			const size_t fc = col + row;
#ifndef INLINE_CUPS
			updateNucDiag(
				tab[row-1][0],     // cell diagonally above and to the left
				cur,               // destination cell
				c,                 // color being traversed
				rf_[rfi_ + fc],    // ref mask at destination cell
				QUAL2(row, fc),    // amt to penalize mismatch
				improved);         // =true if there was path of improvement
#else
			bool match = matches(c, rf_[rfi_ + fc]);
			TAlScore add;
			if(match) {
				add = sc_->match(30);
				improved = true;
			} else {
				add = -QUAL2(row, fc);
			}
			const TCell& dc = tab[row-1][0];
			TCell& dstc = cur;
			AlnScore& myOallBest  = dstc.oallBest;
			SwNucCellMask& myMask = dstc.mask;
			assert(!sc_->monotone || dc.oallBest.score() <= 0);
			if(VALID_AL_SCORE(dc.oallBest)) {
				COUNT_N(dc.oallBest);
				AlnScore ex = dc.oallBest + (int)add;
				if(ex.score_ >= floorsc_ && ex >= myOallBest) {
					if(ex > myOallBest) {
						myOallBest = ex;
						myMask.clear();
					}
					myMask.oall_diag = 1;
				}
				assert(dstc.repOk());
			}
#endif
		}
		if(!onlyDiagInto && col < whi && !tab[row-1][1].empty) {
			// FIRST COLUMN OF A MIDDLE ROW - VERTICAL UPDATE
#ifndef INLINE_CUPS
			updateNucVert(
				tab[row-1][1],     // cell diagonally above and to the left
				cur,               // destination cell
				c);                // color being traversed
#else
			IF_COUNT_N(bool ninvolved = c > 3);
			AlnScore  upOallBest  = tab[row-1][1].oallBest;
			AlnScore  upRfgapBest = tab[row-1][1].rfgapBest;
			AlnScore& myOallBest  = cur.oallBest;
			AlnScore& myRfgapBest = cur.rfgapBest;
			SwNucCellMask& myMask = cur.mask;
			assert(!sc_->monotone || upOallBest.score()  <= 0);
			assert(!sc_->monotone || upRfgapBest.score() <= 0);
			assert_leq(upRfgapBest, upOallBest);
			if(VALID_AL_SCORE(upOallBest)) {
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
						myMask.oall_rfex = 1;
					}
					if(ex >= myRfgapBest) {
						if(ex > myRfgapBest) {
							myMask.clearRefGapMask();
							myRfgapBest = ex;
						}
						myMask.rfgap_ex = 1;
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
						myMask.oall_rfop = 1;
					}
					if(ex >= myRfgapBest) {
						if(ex > myRfgapBest) {
							myMask.clearRefGapMask();
							myRfgapBest = ex;
						}
						myMask.rfgap_op = 1;
					}
				}
				assert(cur.repOk());
			}
#endif
		} // end if(skip_left == 0)
		// 'cur' is now initialized
		ncups_++;
		assert(tab[row][col].empty);
		assert(tab[row][col].repOk());
		assert(!tab[row][col].finalized);
		ASSERT_ONLY(tab[row][col].finalized = true);
		bool aboveFloor =
			tab[row][col].oallBest.valid() &&
			tab[row][col].oallBest.score() >= loBound;
		if(!tab[row][col].mask.empty() && aboveFloor) {
			assert(VALID_AL_SCORE(tab[row][col].oallBest));
			tab[row][col].empty = false;
			UPDATE_SOLS(tab[row][col], row, col, improved);
			assert(!tab[row][col].empty);
			validInRow = true;
		} else if(aboveFloor) {
			assert(VALID_AL_SCORE(tab[row][col].oallBest));
			tab[row][col].terminal = true;
		}
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
			improved = false;
#ifndef INLINE_CUPS
			const int mmpen = QUAL2(row, fullcol);
			updateNucDiag(
				dg,            // cell diagonally above and to the left
				cur,           // destination cell
				c,             // nucleotide in destination row
				r,             // ref mask associated with destination column
				mmpen,         // penalty to incur for color miscall
				improved);     // =true if there was path of improvement
#else
			bool match = matches(c, r);
			TAlScore add;
			if(match) {
				add = sc_->match(30);
				improved = true;
			} else {
				add = -QUAL2(row, fullcol);
			}
			AlnScore& myOallBest  = cur.oallBest;
			SwNucCellMask& myMask = cur.mask;
			assert(!sc_->monotone || dg.oallBest.score() <= 0);
			if(VALID_AL_SCORE(dg.oallBest)) {
				COUNT_N(dg.oallBest);
				AlnScore ex = dg.oallBest + (int)add;
				if(ex.score_ >= floorsc_ && ex >= myOallBest) {
					if(ex > myOallBest) {
						myOallBest = ex;
						myMask.clear();
					}
					myMask.oall_diag = 1;
				}
				assert(cur.repOk());
			}
#endif
			TCell& up = tab[row-1][col-wlo+1];
			if(!onlyDiagInto && !up.empty) {
#ifndef INLINE_CUPS
				updateNucVert(
					up,        // cell above
					cur,       // destination cell
					c);        // nucleotide in destination row
#else
				IF_COUNT_N(bool ninvolved = c > 3);
				AlnScore  upOallBest  = up.oallBest;
				AlnScore  upRfgapBest = up.rfgapBest;
				AlnScore& myOallBest  = cur.oallBest;
				AlnScore& myRfgapBest = cur.rfgapBest;
				SwNucCellMask& myMask = cur.mask;
				assert(!sc_->monotone || upOallBest.score()  <= 0);
				assert(!sc_->monotone || upRfgapBest.score() <= 0);
				assert_leq(upRfgapBest, upOallBest);
				if(VALID_AL_SCORE(upOallBest)) {
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
							myMask.oall_rfex = 1;
						}
						if(ex >= myRfgapBest) {
							if(ex > myRfgapBest) {
								myMask.clearRefGapMask();
								myRfgapBest = ex;
							}
							myMask.rfgap_ex = 1;
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
							myMask.oall_rfop = 1;
						}
						if(ex >= myRfgapBest) {
							if(ex > myRfgapBest) {
								myMask.clearRefGapMask();
								myRfgapBest = ex;
							}
							myMask.rfgap_op = 1;
						}
					}
					assert(cur.repOk());
				}
#endif
			}
			// Can do horizontal
			TCell& lf = tab[row][col-wlo-1];
			if(!onlyDiagInto && !lf.empty) {
				assert(lf.finalized);
#ifndef INLINE_CUPS
				updateNucHoriz(
					lf,        // cell to the left
					cur,       // destination cell
					r);        // ref mask associated with destination column
#else
				IF_COUNT_N(bool ninvolved = r > 15);
				AlnScore leftOallBest  = lf.oallBest;
				AlnScore leftRdgapBest = lf.rdgapBest;
				AlnScore& myOallBest   = cur.oallBest;
				AlnScore& myRdgapBest  = cur.rdgapBest;
				SwNucCellMask& myMask  = cur.mask;
				assert(!myRdgapBest.valid());
				assert(!sc_->monotone || leftOallBest.score()  <= 0);
				assert(!sc_->monotone || leftRdgapBest.score() <= 0);
				assert_leq(leftRdgapBest, leftOallBest);
				if(VALID_AL_SCORE(leftOallBest)) {
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
					assert(cur.repOk());
				}
#endif
			} // end loop over inner diagonals
			// 'cur' is now initialized
			ncups_++;
			assert(tab[row][col].empty);
			assert(tab[row][col].repOk());
			assert(!tab[row][col].finalized);
			ASSERT_ONLY(tab[row][col].finalized = true);
			bool aboveFloor =
				tab[row][col].oallBest.valid() &&
				tab[row][col].oallBest.score() >= loBound;
			if(!tab[row][col].mask.empty() && aboveFloor) {
				assert(VALID_AL_SCORE(tab[row][col].oallBest));
				tab[row][col].empty = false;
				UPDATE_SOLS(tab[row][col], row, col, improved);
				assert(!tab[row][col].empty);
				validInRow = true;
			} else if(aboveFloor) {
				assert(VALID_AL_SCORE(tab[row][col].oallBest));
				tab[row][col].terminal = true;
			}
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
			const int mmpenf = QUAL2(row, fullcol);
			TCell& dg = tab[row-1][col-wlo];
			improved = false;
			if(!dg.empty) {
				assert(dg.finalized);
#ifndef INLINE_CUPS
				updateNucDiag(
					dg,        // cell diagonally above and to the left
					cur,       // destination cell
					c,         // nucleotide in destination row
					r,         // ref mask associated with destination column
					mmpenf,    // penalty to incur for color miscall
					improved); // =true if there was path of improvement
#else
				bool match = matches(c, r);
				TAlScore add;
				if(match) {
					add = sc_->match(30);
					improved = true;
				} else {
					add = -mmpenf;
				}
				AlnScore& myOallBest  = cur.oallBest;
				SwNucCellMask& myMask = cur.mask;
				assert(!sc_->monotone || dg.oallBest.score() <= 0);
				if(VALID_AL_SCORE(dg.oallBest)) {
					COUNT_N(dg.oallBest);
					AlnScore ex = dg.oallBest + (int)add;
					if(ex.score_ >= floorsc_ && ex >= myOallBest) {
						if(ex > myOallBest) {
							myOallBest = ex;
							myMask.clear();
						}
						myMask.oall_diag = 1;
					}
					assert(cur.repOk());
				}
#endif
			}
			TCell& lf = tab[row][col-wlo-1];
			if(!onlyDiagInto && !lf.empty) {
				assert(lf.finalized);
#ifndef INLINE_CUPS
				updateNucHoriz(
					lf,        // cell to the left
					cur,       // destination cell
					r);        // ref mask associated with destination column
#else
				IF_COUNT_N(bool ninvolved = r > 15);
				AlnScore leftOallBest  = lf.oallBest;
				AlnScore leftRdgapBest = lf.rdgapBest;
				AlnScore& myOallBest   = cur.oallBest;
				AlnScore& myRdgapBest  = cur.rdgapBest;
				SwNucCellMask& myMask  = cur.mask;
				assert(!myRdgapBest.valid());
				assert(!sc_->monotone || leftOallBest.score()  <= 0);
				assert(!sc_->monotone || leftRdgapBest.score() <= 0);
				assert_leq(leftRdgapBest, leftOallBest);
				if(VALID_AL_SCORE(leftOallBest)) {
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
					assert(cur.repOk());
				}
#endif
			}
			// 'cur' is now initialized
			ncups_++;
			assert(tab[row][col].empty);
			assert(tab[row][col].repOk());
			assert(!tab[row][col].finalized);
			ASSERT_ONLY(tab[row][col].finalized = true);
			bool aboveFloor =
				tab[row][col].oallBest.valid() &&
				tab[row][col].oallBest.score() >= loBound;
			if(!tab[row][col].mask.empty() && aboveFloor) {
				assert(VALID_AL_SCORE(tab[row][col].oallBest));
				tab[row][col].empty = false;
				UPDATE_SOLS(tab[row][col], row, col, improved);
				assert(!tab[row][col].empty);
				validInRow = true;
			} else if(aboveFloor) {
				assert(VALID_AL_SCORE(tab[row][col].oallBest));
				tab[row][col].terminal = true;
			}
		}
		if(!validInRow) {
			assert_geq(rdf_-rdi_, row+1);
			nrowskips_ += (rdf_ - rdi_ - row - 1);
			return !soldone_;
		}
	}
	return !soldone_;
}

/**
 * Populate the given SwResult with information about the "next best"
 * alignment if there is one.  If there isn't one, false is returned.  Note
 * that false might be returned even though a call to done() would have
 * returned false.
 *
 * Which alignment is "next best" depends on the 'rowfirst_' setting.  If
 * it's true, then alignments in rows further down in the DP matrix get
 * higher priority than rows higher up.  Within a row, alignments are
 * ordered according to their score.  If 'rowfirst_' is false, then
 * priority is given first according to score, then according to row.
 */
bool SwAligner::nextAlignment(
	SwResult& res,
	RandomSource& rnd)
{
	assert(inited());
	assert_eq(STATE_ALIGNED, state_);
	assert(repOk());
	if(soldone_) {
		return false;
	}
	// Must be rows with solutions
	assert_geq(solrows_.second, solrows_.first);
	assert_geq((int64_t)solrows_.first, solrowlo_);
	assert_geq((int64_t)solrows_.second, solrowlo_);
	assert_gt(solrowbest_[solrows_.first ],
	          std::numeric_limits<TAlScore>::min()); // should be valid
	assert_gt(solrowbest_[solrows_.second],
	          std::numeric_limits<TAlScore>::min()); // should be valid
	// Look for the next cell to backtrace from
	size_t off = std::numeric_limits<size_t>::max();
	if(sc_->rowFirst) {
		// Iterate over rows that potentially have solution cells.  Start
		// at the bottom and move up.  Within each row, take solutions
		// from best score to worst score.
		bool done = false;
		while(!done) {
			// From low row to high row
			for(size_t ii = solrows_.second+1; ii >= solrows_.first+1; ii--) {
				size_t i = ii-1;
#ifndef NDEBUG
				// All lower rows should be exhausted
				for(size_t j = i+1; j < solrowbest_.size(); j++) {
					assert(!VALID_SCORE(solrowbest_[j]));
					assert_gt(solcols_[j].first, solcols_[j].second);
				}
#endif
				if(VALID_SCORE(solrowbest_[i])) {
					// This row could have a solution cell
					assert(VALID_SCORE(solbest_));
					assert_geq(solrowbest_[i], minsc_);
					assert_leq(solrowbest_[i], solbest_);
					if(color_) {
						done = nextAlignmentTryColorRow(
							res,             // put alignment result here
							off,             // ref offset of alignment
							i,               // row to trace back from
							solrowbest_[i],  // score to select
							rnd);            // pseudo-random source
					} else {
						done = nextAlignmentTryNucRow(
							res,             // put alignment result here
							off,             // ref offset of alignment
							i,               // row to trace back from
							solrowbest_[i],  // score to select
							rnd);            // pseudo-random source
					}
					if(done) {
						// Found the alignment we were looking for
						break;
					}
				}
				if(solrows_.first == std::numeric_limits<size_t>::max()) {
					break;
				}
			}
			if(!done) {
				// Finished with *all* distinct alignments
				soldone_ = true;
				assert(repOk());
				// Reset res because we have have added e.g. edits to it in the
				// process of discovering there were no alignments
				res.reset();
				return false;
			} else {
				// Got a solution from one of the rows
				assert(repOk());
			}
		}
	} else {
		// Iterate first over scores from best score.  Then, within each score
		// stratum, take solutions from latest row to earliest row.
		bool done = false;
		while(!done) {
			// From low row to high row
			while(!soldone_) {
				for(size_t ii = solrows_.second+1; ii >= solrows_.first+1; ii--) {
					size_t i = ii-1;
					if(VALID_SCORE(solrowbest_[i]) && solrowbest_[i] == solbest_) {
						// This row could have a solution cell
						assert(VALID_SCORE(solbest_));
						assert(VALID_SCORE(solrowbest_[i]));
						TAlScore oldbest = solbest_;
						if(color_) {
							done = nextAlignmentTryColorRow(
								res,             // put alignment result here
								off,             // ref offset of alignment
								i,               // row to trace back from
								solbest_,        // score to select
								rnd);            // pseudo-random source
						} else {
							done = nextAlignmentTryNucRow(
								res,             // put alignment result here
								off,             // ref offset of alignment
								i,               // row to trace back from
								solbest_,        // score to select
								rnd);            // pseudo-random source
						}
						if(done) {
							break; // Found the alignment
						}
						if(oldbest != solbest_) {
							// The best overall score changed.  When this
							// happens, we go back to the deepest row and start
							// our row scan over.
							ii = solrows_.second+1+1;
						}
					}
					if(solrows_.first == std::numeric_limits<size_t>::max()) {
						break;
					}
				}
				if(done) {
					break; // Found the alignment
				}
			}
			if(!done) {
				// Finished with *all* distinct alignments
				assert(!VALID_SCORE(solbest_));
				soldone_ = true;
				assert(repOk());
				// Reset res because we have have added e.g. edits to it in the
				// process of discovering there were no alignments
				res.reset();
				return false;
			} else {
				assert(repOk());
				break; // Found the alignment
			}
		}
	}
	// Alignment offset should have been set in call to nextAlignmentTry...
	assert_neq(std::numeric_limits<size_t>::max(), off);
	assert(!res.alres.empty());
	assert(res.repOk());
	if(!fw_) {
		// All edits are currently w/r/t upstream end; if read aligned
		// to Crick strand, we need to invert them so that they're
		// w/r/t the read's 5' end instead.
		res.alres.invertEdits();
	}
	if(color_) {
		assert_range(0, 3, res.nup);
		assert_range(0, 3, res.ndn);
		res.alres.setNucs(fw_, res.nup, res.ndn);
	}
	cural_++;
	return true;
}

/**
 * If there are empty rows on either end, trim them off.  If we're out of rows
 * with potential solutions, reset everything.
 */
void SwAligner::nextAlignmentUpdateBest() {
	assert_geq(solrows_.second, solrows_.first);
	assert_geq((int64_t)solrows_.first, solrowlo_);
	assert_geq((int64_t)solrows_.second, solrowlo_);
	ASSERT_ONLY(TAlScore oldbest = solbest_);
	INVALIDATE_SCORE(solbest_);
	size_t first = 0;
	size_t last = 0;
	bool setFirst = false;
	bool setLast = false;
	for(size_t i = solrows_.first; i <= solrows_.second; i++) {
		if(VALID_SCORE(solrowbest_[i])) {
			if(!setFirst) {
				first = i;
				setFirst = true;
			}
			setLast = true;
			last = i;
			if(solrowbest_[i] > solbest_) {
				solbest_ = solrowbest_[i];
				assert(!sc_->monotone || solbest_ <= 0);
			}
		}
	}
	if(setFirst) {
		// There are still some rows with potential solutions
		solrows_.first = first;
		solrows_.second = last;
		assert(!soldone_);
		assert(VALID_SCORE(solbest_));
	} else {
		// No more rows with potential solutions
		solrows_ = EXTREMES;
		soldone_ = true;
		assert(!VALID_SCORE(solbest_));
	}
	assert(solbest_ <= oldbest || !VALID_SCORE(oldbest));
}

/**
 * Given a range of columns in a row, update solrowbest_, solcols_, solbest_.
 * Return true iff there areone or more solution cells in the row.
 */
bool SwAligner::nextAlignmentUpdateColorRow(size_t row) {
	assert(repOk());
	// For each column that used to be viable
	INVALIDATE_SCORE(solrowbest_[row]);
	if(solcols_[row].second < solcols_[row].first) {
		// The columns were never refined
		solcols_[row].first = 0;
		solcols_[row].second = width_;
	}
	// Last column with a solution cell
	size_t last = std::numeric_limits<size_t>::max();
	bool first = true;
	// For each column that might have a solution cell...
	assert(en_ == NULL || solcols_[row].second <= en_->size());
	for(size_t i = solcols_[row].first; i <= solcols_[row].second; i++) {
		if(en_ != NULL && !(*en_)[i]) {
			// Skip b/c this diagonal cannot have a solution cell
			continue;
		}
		// Is it a solution cell?
		AlnScore best;
		if(ctab_[row][i].bestSolutionGeq(minsc_, best)) {
			// Yes
			if(first) {
				// If it's the first solution cell, set the left column
				// boundary to this column.
				solcols_[row].first = i;
				first = false;
			}
			if(best.score() > solrowbest_[row]) {
				// If it's the best in the row so far, update the row-best
				solrowbest_[row] = best.score();
				assert(!sc_->monotone || solrowbest_[row] <= 0);
				if(best.score() > solbest_) {
					// If it's the best overall, update the best overall
					solbest_ = best.score();
				}
			} else {
				// If it's not better than row-best, can't be better than
				// overall-best
				assert_leq(best.score(), solbest_);
			}
			last = i;
		}
	}
	// Did we find any solution cells?
	if(last < std::numeric_limits<size_t>::max()) {
		// Yes, found at least 1 solution cell.  Update the rightmost
		// column with a solution cell.
		solcols_[row].second = last;
		return true; // success
	} else {
		// No, no solution cells found in this row
		assert(!VALID_SCORE(solrowbest_[row]));
		solcols_[row] = SwAligner::EXTREMES;
		return false; // no solutions
	}
}

/**
 * Given a range of columns in a row, update solrowbest_, solcols_, solbest_.
 * Return true iff there areone or more solution cells in the row.
 */
bool SwAligner::nextAlignmentUpdateNucRow(size_t row) {
	assert(repOk());
	// For each column that used to be viable
	INVALIDATE_SCORE(solrowbest_[row]);
	if(solcols_[row].second < solcols_[row].first) {
		// The columns were never refined
		solcols_[row].first = 0;
		solcols_[row].second = width_;
	}
	// Last column with a solution cell
	size_t last = std::numeric_limits<size_t>::max();
	bool first = true;
	// For each column that might have a solution cell...
	assert(en_ == NULL || solcols_[row].second <= en_->size());
	for(size_t i = solcols_[row].first; i <= solcols_[row].second; i++) {
		if(en_ != NULL && !(*en_)[i]) {
			// Skip b/c this diagonal cannot have a solution cell
			continue;
		}
		// Is it a solution cell?
		AlnScore best;
		if(ntab_[row][i].bestSolutionGeq(minsc_, best)) {
			// Yes
			if(first) {
				// If it's the first solution cell, set the left column
				// boundary to this column.
				solcols_[row].first = i;
				first = false;
			}
			if(best.score() > solrowbest_[row]) {
				// If it's the best in the row so far, update the row-best
				solrowbest_[row] = best.score();
				assert(!sc_->monotone || solrowbest_[row] <= 0);
				if(best.score() > solbest_) {
					// If it's the best overall, update the best overall
					solbest_ = best.score();
				}
			} else {
				// If it's not better than row-best, can't be better than
				// overall-best
				assert_leq(best.score(), solbest_);
			}
			last = i;
		}
	}
	// Did we find any solution cells?
	if(last < std::numeric_limits<size_t>::max()) {
		// Yes, found at least 1 solution cell.  Update the rightmost
		// column with a solution cell.
		solcols_[row].second = last;
		return true; // success
	} else {
		// No, no solution cells found in this row
		assert(!VALID_SCORE(solrowbest_[row]));
		solcols_[row] = SwAligner::EXTREMES;
		return false; // no solutions
	}
}

/**
 * Try to report an alignment with the given score from the given row.
 * After each attempt, update information about the row.  Returns false
 * iff the row was exhausted.
 */
bool SwAligner::nextAlignmentTryColorRow(
	SwResult& res,       // install result here
	size_t& off,         // ref offset of alignment
	size_t row,          // row to try a solution from
	TAlScore sc,         // potential solutions must have this score
	RandomSource& rnd)   // pseudo-random generator for backtrack
{
	assert(repOk());
	assert_lt(row, solcols_.size());
	assert_lt(row, ctab_.size());
	assert_geq(solcols_[row].second, solcols_[row].first);
	assert(VALID_SCORE(solrowbest_[row]));
	// Try each column that might contain a solution cell.  Iterate over
	// rectangle columns.
	bool triedOne = false;
	for(size_t i = solcols_[row].first; i <= solcols_[row].second; i++) {
		assert_lt(i, ctab_[row].size());
		if(en_ != NULL && !(*en_)[i]) {
			// This can't be at either end of the range of valid row cells
			// TODO: had to comment these out - why?
			//assert_neq(i, solcols_[row].first);
			//assert_neq(i, solcols_[row].second);
			// Skip b/c this diagonal cannot have a solution cell
			continue;
		}
		// Does this cell have the desired score and is it legal to end in
		// this diagonal?
		int nuc = -1;
		if(ctab_[row][i].nextSolutionEq(sc, nuc)) {
			assert_range(0, 3, nuc);
			triedOne = true;
			bool ret = backtrackColors(
				sc,    // score we expect to get over backtrack
				res,   // store results (edits and scores) here
				off,   // store result offset here
				row,   // start in this parallelogram row
				i+row, // start in this parallelogram column
				nuc,   // final decoded nucleotide
				rnd);  // pseudo-random generator
			nextAlignmentUpdateColorRow(row);
			nextAlignmentUpdateBest();
			// If the row is exhausted, the caller will be able to see by
			// checking whether solrowbest_[row] is valid, for example
			if(ret) {
				return true; // found a result
			}
			if(!VALID_SCORE(solrowbest_[row])) {
				return false; // definitely won't find a result
			}
		}
	}
	if(!triedOne) {
		nextAlignmentUpdateColorRow(row);
		nextAlignmentUpdateBest();
	}
	return false; // no result
}

/**
 * Try to report an alignment with the given score from the given row.
 * After each attempt, update information about the row.
 */
bool SwAligner::nextAlignmentTryNucRow(
	SwResult& res,       // install result here
	size_t& off,         // ref offset of alignment
	size_t row,          // row to try a solution from
	TAlScore sc,         // potential solutions must have this score
	RandomSource& rnd)   // pseudo-random generator for backtrack
{
	assert(repOk());
	assert_lt(row, solcols_.size());
	assert_lt(row, ntab_.size());
	assert_geq(solcols_[row].second, solcols_[row].first);
	assert(VALID_SCORE(solrowbest_[row]));
	// Try each column that might contain a solution cell.  Iterate over
	// rectangle columns.
	bool triedOne = false;
	for(size_t i = solcols_[row].first; i <= solcols_[row].second; i++) {
		assert_lt(i, ntab_[row].size());
		if(en_ != NULL && !(*en_)[i]) {
			// This can't be at either end of the range of valid row cells
			// TODO: had to comment these out - why?
			//assert_neq(i, solcols_[row].first);
			//assert_neq(i, solcols_[row].second);
			// Skip b/c this diagonal cannot have a solution cell
			continue;
		}
		// Does this cell have the desired score and is it legal to end in
		// this diagonal?
		if(ntab_[row][i].nextSolutionEq(sc)) {
			triedOne = true;
			bool ret = backtrackNucleotides(
				sc,    // score we expect to get over backtrack
				res,   // store results (edits and scores) here
				off,   // store result offset here
				row,   // start in this parallelogram row
				i+row, // start in this parallelogram column
				rnd);  // pseudo-random generator
			nextAlignmentUpdateNucRow(row);
			nextAlignmentUpdateBest();
			// If the row is exhausted, the caller will be able to see by
			// checking whether solrowbest_[row] is valid, for example
			if(ret) {
				return true; // found a result
			}
			if(!VALID_SCORE(solrowbest_[row])) {
				return false; // definitely won't find a result
			}
		}
	}
	if(!triedOne) {
		nextAlignmentUpdateNucRow(row);
		nextAlignmentUpdateBest();
	}
	return false;
}

#ifdef MAIN_ALIGNER_SW

#include <sstream>
#include <utility>
#include <getopt.h>
#include "scoring.h"
#include "aligner_seed_policy.h"
#include "color.h"

int  gGapBarrier;
int  gSnpPhred;
bool gColor;
bool gColorExEnds;
static int bonusMatchType;   // how to reward matches
static int bonusMatch;       // constant if match bonus is a constant
static int penMmcType;       // how to penalize mismatches
static int penMmc;           // constant if mm pelanty is a constant
static int penSnp;           // penalty for nucleotide mm in decoded color aln
static int penNType;         // how to penalize Ns in the read
static int penN;             // constant if N pelanty is a constant
static bool nPairCat;        // true -> concatenate mates before N filter
static int penRdExConst;     // constant coeff for cost of gap in read
static int penRfExConst;     // constant coeff for cost of gap in ref
static int penRdExLinear;    // linear coeff for cost of gap in read
static int penRfExLinear;    // linear coeff for cost of gap in ref
static float costMinConst;   // constant coeff for min score w/r/t read len
static float costMinLinear;  // linear coeff for min score w/r/t read len
static float costFloorConst; // constant coeff for score floor w/r/t read len
static float costFloorLinear;// linear coeff for score floor w/r/t read len
static float nCeilConst;     // constant coeff for N ceiling w/r/t read len
static float nCeilLinear;    // linear coeff for N ceiling w/r/t read len
static bool  nCatPair;       // concat mates before applying N filter?
static int multiseedMms;     // mismatches permitted in a multiseed seed
static int multiseedLen;     // length of multiseed seeds
static int multiseedPeriod;  // space between multiseed seeds
static int multiseedIvalType;
static float multiseedIvalA;
static float multiseedIvalB;
static float posmin;
static float posfrac;
static float rowmult;

enum {
	ARG_TESTS = 256
};

static const char *short_opts = "s:m:r:d:i:";
static struct option long_opts[] = {
	{(char*)"snppen",       required_argument, 0, 's'},
	{(char*)"misspen",      required_argument, 0, 'm'},
	{(char*)"seed",         required_argument, 0, 'r'},
	{(char*)"color",        no_argument,       0, 'C'},
	{(char*)"align-policy", no_argument,       0, 'A'},
	{(char*)"test",         no_argument,       0, ARG_TESTS},
};

static void printUsage(ostream& os) {
	os << "Usage: aligner_sw <read-seq> <ref-nuc-seq> [options]*" << endl;
	os << "Options:" << endl;
	os << "  -C/--color          reads are colorspace" << endl;
	os << "  -s/--snppen <int>   penalty incurred by SNP; used for decoding"
	   << endl;
	os << "  -m/--misspen <int>  quality to use for read chars" << endl;
	os << "  -r/-seed <int>      seed for pseudo-random generator" << endl;
}

/**
 * Parse a T from a string 's'
 */
template<typename T>
T parse(const char *s) {
	T tmp;
	stringstream ss(s);
	ss >> tmp;
	return tmp;
}

static EList<bool> stbuf, enbuf;
static BTDnaString btread;
static BTString btqual;
static BTString btref;
static BTString btref2;

/**
 * Helper function for running a case consisting of a read (sequence
 * and quality), a reference string, and an offset that anchors the 0th
 * character of the read to a reference position.
 */
static void doTestCase(
	SwAligner&         al,
	const BTDnaString& read,
	const BTString&    qual,
	const BTString&    refin,
	int64_t            off,
	EList<bool>       *st,
	EList<bool>       *en,
	const Scoring&     sc,  
	TAlScore           minsc,
	TAlScore           floorsc,
	SwResult&          res,
	bool               color,
	bool               nsInclusive,
	bool               filterns,
	uint32_t           seed)
{
	RandomSource rnd(seed);
	btref2 = refin;
	assert_eq(read.length(), qual.length());
	size_t nrow = read.length();
	nrow += (color ? 1 : 0);
	int64_t rfi, rff;
	// Calculate the largest possible number of read and reference gaps given
	// 'minsc' and 'pens'
	int readGaps = sc.maxReadGaps(minsc, read.length());
	int refGaps = sc.maxRefGaps(minsc, read.length());
	assert_geq(readGaps, 0);
	assert_geq(refGaps, 0);
	int maxGaps = max(readGaps, refGaps);
	size_t nceil = sc.nCeil(read.length());
	size_t width = 1 + 2 * maxGaps;
	rfi = off;
	off = 0;
	// Pad the beginning of the reference with Ns if necessary
	if(rfi < maxGaps) {
		size_t beginpad = (size_t)(maxGaps - rfi);
		for(size_t i = 0; i < beginpad; i++) {
			btref2.insert('N', 0);
			off--;
		}
		rfi = 0;
	} else {
		rfi -= maxGaps;
	}
	assert_geq(rfi, 0);
	// Pad the end of the reference with Ns if necessary
	while(rfi + nrow + 2 * maxGaps > btref2.length()) {
		btref2.append('N');
	}
	rff = rfi + nrow + (2 * maxGaps);
	// Convert reference string to masks
	for(size_t i = 0; i < btref2.length(); i++) {
		if(toupper(btref2[i]) == 'N' && !nsInclusive) {
			btref2.set(16, i);
		} else {
			int num = 0;
			int alts[] = {4, 4, 4, 4};
			decodeNuc(toupper(btref2[i]), num, alts);
			assert_leq(num, 4);
			assert_gt(num, 0);
			btref2.set(0, i);
			for(int j = 0; j < num; j++) {
				btref2.set(btref2[i] | (1 << alts[j]), i);
			}
		}
	}
	bool fw = true;
	uint32_t refidx = 0;
	if(st == NULL) {
		stbuf.resize(width);
		stbuf.fill(true);
		st = &stbuf;
	}
	if(en == NULL) {
		enbuf.resize(width);
		enbuf.fill(true);
		en = &enbuf;
	}
	assert_geq(rfi, 0);
	assert_gt(rff, rfi);
	al.init(
		read,          // read sequence
		qual,          // read qualities
		0,             // offset of first character within 'read' to consider
		read.length(), // offset of last char (exclusive) in 'read' to consider
		fw,            // 'read' is forward version of read?
		color,         // whether read is nucleotide-space or colorspace
		refidx,        // id of reference aligned to
		off,           // offset of upstream ref char aligned against
		btref2.wbuf(), // reference sequence (masks)
		rfi,           // offset of first char in 'ref' to consider
		rff,           // offset of last char (exclusive) in 'ref' to consider
		width,         // # bands to do (width of parallelogram)
		st,            // mask indicating which columns we can start in
		en,            // mask indicating which columns we can end in
		sc,            // scoring scheme
		minsc,         // minimum score for valid alignment
		floorsc);      // local-alignment score floor
	if(filterns) {
		al.filter((int)nceil);
	}
	al.align(rnd);
}

/**
 * Another interface for running a case.
 */
static void doTestCase2(
	SwAligner&         al,
	const char        *read,
	const char        *qual,
	const char        *refin,
	int64_t            off,
	const Scoring&     sc,
	float              costMinConst,
	float              costMinLinear,
	SwResult&          res,
	bool               color,
	bool               nsInclusive = false,
	bool               filterns = false,
	uint32_t           seed = 0)
{
	btread.install(read, true, color);
	TAlScore minsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costMinConst,
		costMinLinear));
	TAlScore floorsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costFloorConst,
		costFloorLinear));
	btqual.install(qual);
	btref.install(refin);
	doTestCase(
		al,
		btread,
		btqual,
		btref,
		off,
		NULL,
		NULL,
		sc,  
		minsc,
		floorsc,
		res,
		color,
		nsInclusive,
		filterns,
		seed
	);
}

/**
 * Another interface for running a case.
 */
static void doTestCase3(
	SwAligner&         al,
	const char        *read,
	const char        *qual,
	const char        *refin,
	int64_t            off,
	Scoring&           sc,
	float              costMinConst,
	float              costMinLinear,
	float              nCeilConst,
	float              nCeilLinear,
	SwResult&          res,
	bool               color,
	bool               nsInclusive = false,
	bool               filterns = false,
	uint32_t           seed = 0)
{
	btread.install(read, true, color);
	// Calculate the penalty ceiling for the read
	TAlScore minsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costMinConst,
		costMinLinear));
	TAlScore floorsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costFloorConst,
		costFloorLinear));
	btqual.install(qual);
	btref.install(refin);
	sc.nCeilConst = nCeilConst;
	sc.nCeilLinear = nCeilLinear;
	doTestCase(
		al,
		btread,
		btqual,
		btref,
		off,
		NULL,
		NULL,
		sc,  
		minsc,
		floorsc,
		res,
		color,
		nsInclusive,
		filterns,
		seed
	);
}

/**
 * Do a set of unit tests.
 */
static void doTests() {
	bonusMatchType  = DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch      = DEFAULT_MATCH_BONUS;
	penMmcType      = DEFAULT_MM_PENALTY_TYPE;
	penMmc          = DEFAULT_MM_PENALTY;
	penSnp          = DEFAULT_SNP_PENALTY;
	penNType        = DEFAULT_N_PENALTY_TYPE;
	penN            = DEFAULT_N_PENALTY;
	nPairCat        = DEFAULT_N_CAT_PAIR;
	penRdExConst    = DEFAULT_READ_GAP_CONST;
	penRfExConst    = DEFAULT_REF_GAP_CONST;
	penRdExLinear   = DEFAULT_READ_GAP_LINEAR;
	penRfExLinear   = DEFAULT_REF_GAP_LINEAR;
	costMinConst    = DEFAULT_MIN_CONST;
	costMinLinear   = DEFAULT_MIN_LINEAR;
	costFloorConst  = DEFAULT_FLOOR_CONST;
	costFloorLinear = DEFAULT_FLOOR_LINEAR;
	nCeilConst      = 1.0f; // constant factor in N ceil w/r/t read len
	nCeilLinear     = 0.1f; // coeff of linear term in N ceil w/r/t read len
	multiseedMms    = DEFAULT_SEEDMMS;
	multiseedLen    = DEFAULT_SEEDLEN;
	multiseedPeriod = DEFAULT_SEEDPERIOD;
	// Set up penalities
	Scoring sc(
		bonusMatch,
		penMmcType,    // how to penalize mismatches
		penMmc,        // constant if mm pelanty is a constant
		penSnp,        // penalty for decoded SNP
		costMinConst,  // constant factor in N ceiling w/r/t read length
		costMinLinear, // coeff of linear term in N ceiling w/r/t read length
		costFloorConst,  // constant factor in N ceiling w/r/t read length
		costFloorLinear, // coeff of linear term in N ceiling w/r/t read length
		nCeilConst,    // constant factor in N ceiling w/r/t read length
		nCeilLinear,   // coeff of linear term in N ceiling w/r/t read length
		penNType,      // how to penalize Ns in the read
		penN,          // constant if N pelanty is a constant
		nPairCat,      // true -> concatenate mates before N filtering
		penRdExConst,  // constant coeff for cost of gap in read
		penRfExConst,  // constant coeff for cost of gap in ref
		penRdExLinear, // linear coeff for cost of gap in read
		penRfExLinear, // linear coeff for cost of gap in ref
		1,             // # rows at top/bot can only be entered diagonally
		-1,            // min row idx to backtrace from; -1 = no limit
		false          // sort results first by row then by score?
	);
	// Set up alternative penalities
	Scoring sc2(
		bonusMatch,
		COST_MODEL_QUAL, // how to penalize mismatches
		penMmc,          // constant if mm pelanty is a constant
		penSnp,          // penalty for decoded SNP
		costMinConst,  // constant factor in N ceiling w/r/t read length
		costMinLinear, // coeff of linear term in N ceiling w/r/t read length
		costFloorConst,  // constant factor in N ceiling w/r/t read length
		costFloorLinear, // coeff of linear term in N ceiling w/r/t read length
		1.0f,            // constant factor in N ceiling w/r/t read length
		1.0f,            // coeff of linear term in N ceiling w/r/t read length
		penNType,        // how to penalize Ns in the read
		penN,            // constant if N pelanty is a constant
		nPairCat,        // true -> concatenate mates before N filtering
		penRdExConst,    // constant coeff for cost of gap in read
		penRfExConst,    // constant coeff for cost of gap in ref
		penRdExLinear,   // linear coeff for cost of gap in read
		penRfExLinear,   // linear coeff for cost of gap in ref
		1,               // # rows at top/bot can only be entered diagonally
		-1,              // min row idx to backtrace from; -1 = no limit
		false            // sort results first by row then by score?
	);
	SwResult res;
	
	//
	// Basic nucleotide-space tests
	//
	cerr << "Running tests..." << endl;
	int tests = 1;
	bool nIncl = false;
	bool nfilter = false;

	SwAligner al;
	RandomSource rnd(73);
	for(int i = 0; i < 3; i++) {
		cerr << "  Test " << tests++ << " (nuc space, offset "
		     << (i*4) << ", exact)...";
		sc.rdGapConst = 40;
		sc.rfGapConst = 40;
		doTestCase2(
			al,
			"ACGTACGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty()); assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1mm allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTTCGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -30);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1mm allowed by minsc, check qual 1)...";
		doTestCase2(
			al,
			"ACGTTCGT",         // read
			"ABCDEFGH",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc2,                // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		size_t lo, hi;
		if(i == 0) {
			lo = 0; hi = 1;
		} else if(i == 1) {
			lo = 1; hi = 2;
		} else {
			lo = 2; hi = 3;
		}
		for(size_t j = lo; j < hi; j++) {
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert_eq(j*4, res.alres.refoff());
			assert_eq(8, res.alres.refExtent());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), -36);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
		}
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1mm allowed by minsc, check qual 2)...";
		doTestCase2(
			al,
			"ACGAACGT",         // read
			"ABCDEFGH",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc2,                // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -35);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(res.empty());
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1mm allowed by minsc, check qual )...";
		assert(res.empty());
		doTestCase2(
			al,
			"TCGTACGT",         // read
			"ABCDEFGH",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc2,                // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -32);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(res.empty());
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1mm at the beginning, allowed by minsc)...";
		doTestCase2(
			al,
			"CCGTACGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -30);
		assert_eq(res.alres.score().ns(), 0);
		assert_eq(1, res.alres.ned().size());
		assert_eq(0, res.alres.aed().size());
		assert_eq(0, res.alres.ced().size());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 n in read, allowed)...";
		doTestCase3(
			al,
			"ACGTNCGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			1.0f,               // allow 1 N
			0.0f,               // allow 1 N
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -1);
		assert_eq(res.alres.score().ns(), 1);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 2 n in read, allowed)...";
		doTestCase3(
			al,
			"ACGNNCGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			2.0f,               // const coeff for N ceiling
			0.0f,               // linear coeff for N ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -2);
		assert_eq(res.alres.score().ns(), 2);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 2 n in read, 1 at beginning, allowed)...";
		doTestCase2(
			al,
			"NCGTNCGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -2);
		assert_eq(res.alres.score().ns(), 2);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 n in ref, allowed)...";
		doTestCase2(
			al,
			"ACGTACGT",         // read
			"IIIIIIII",         // qual
			"ACGTNCGTACGTANGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -1);
		assert_eq(res.alres.score().ns(), 1);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1mm disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTTCGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-10.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		// Read gap with equal read and ref gap penalties
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", read gap allowed by minsc)...";
		assert(res.empty());
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 15;
		sc.rdGapLinear = 15;
		doTestCase2(
			al,
			"ACGTCGT",          // read
			"IIIIIII",          // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", read gap disallowed by minsc)...";
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 15;
		sc.rdGapLinear = 15;
		doTestCase2(
			al,
			"ACGTCGT",          // read
			"IIIIIII",          // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();

		cerr << "PASSED" << endl;
		// Ref gap with equal read and ref gap penalties
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", ref gap allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",        // read
			"IIIIIIIII",        // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", ref gap disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",        // read
			"IIIIIIIII",        // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(al.done());
		cerr << "PASSED" << endl;
		
		// Read gap with one read gap and zero ref gaps allowed
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 read gap, ref gaps disallowed by minsc)...";
		sc.rfGapConst = 35;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 20;
		sc.rdGapLinear = 10;
		doTestCase2(
			al,
			"ACGTCGT",          // read
			"IIIIIII",          // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -35);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", gaps disallowed by minsc)...";
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 10;
		sc.rdGapLinear = 10;
		doTestCase2(
			al,
			"ACGTCGT",          // read
			"IIIIIII",          // qual 
			"ACGTACGTACGTACGT", // ref 
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(res.empty());
		cerr << "PASSED" << endl;
		
		// Ref gap with one ref gap and zero read gaps allowed
		sc.rfGapConst = 25;
		sc.rdGapConst = 35;
		sc.rfGapLinear = 12;
		sc.rdGapLinear = 22;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 ref gap, read gaps disallowed by minsc)...";
		assert(res.empty());
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -37);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
			<< ", gaps disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		cerr << "PASSED" << endl;
		
		// Read gap with one read gap and two ref gaps allowed
		sc.rfGapConst = 20;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 10;
		sc.rdGapLinear = 15;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 read gap, 2 ref gaps allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", gaps disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		cerr << "PASSED" << endl;

		// Ref gap with one ref gap and two read gaps allowed
		sc.rfGapConst = 25;
		sc.rdGapConst = 10;
		sc.rfGapLinear = 15;
		sc.rdGapLinear = 10;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 ref gap, 2 read gaps allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", gaps disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(al.done());
		cerr << "PASSED" << endl;
		
		// Read gap with two read gaps and two ref gaps allowed
		sc.rfGapConst = 10;
		sc.rdGapConst = 10;
		sc.rfGapLinear = 10;
		sc.rdGapLinear = 10;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 2 ref gaps, 2 read gaps allowed by minsc)...";
		doTestCase3(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			1.0,                // const coeff for N ceiling
			0.0,                // linear coeff for N ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			true);              // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		if(!res.empty()) {
			al.printResultStacked(res, cerr); cerr << endl;
		}
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -20);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		if(!res.empty()) {
			al.printResultStacked(res, cerr); cerr << endl;
		}
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 ref gap, 1 read gap allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -20);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		
		// Ref gap with two ref gaps and zero read gaps allowed
		sc.rfGapConst = 15;
		sc.rdGapConst = 15;
		sc.rfGapLinear = 5;
		sc.rdGapLinear = 5;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 2 ref gaps, 2 read gaps allowed by minsc)...";
		// Careful: it might be possible for the read to align with overhang
		// instead of with a gap
		doTestCase3(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-35.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			1.0f,               // needed to avoid overhang alignments
			0.0f,               // needed to avoid overhang alignments
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			true);              // filter Ns
		if(i == 0) {
			lo = 0; hi = 2;
		} else if(i == 1) {
			lo = 1; hi = 3;
		} else {
			lo = 2; hi = 3;
		}
		for(size_t j = lo; j < hi; j++) {
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			al.printResultStacked(res, cerr); cerr << endl;
			assert(res.alres.refoff() == 0 ||
			       res.alres.refoff() == 4 ||
				   res.alres.refoff() == 8);
			assert_eq(8, res.alres.refExtent());
			assert_eq(res.alres.score().gaps(), 1);
			assert_eq(res.alres.score().score(), -20);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
		}
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 4;
		sc.rdGapLinear = 4;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 ref gap, 1 read gap allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -29);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", short read)...";
		doTestCase2(
			al,
			"A",
			"I",
			"AAAAAAAAAAAA",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;

		if(i == 0) {
			cerr << "  Test " << tests++
			     << " (nuc space, offset 0, short read & ref)...";
			doTestCase2(
				al,
				"A",
				"I",
				"A",
				0,                  // off
				sc,                 // scoring scheme
				-30.0f,             // const coeff for cost ceiling
				0.0f,               // linear coeff for cost ceiling
				res,                // result
				false,              // color
				nIncl,              // Ns inclusive (not mismatches)
				nfilter);           // filter Ns
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), 0);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
			cerr << "PASSED" << endl;
		}

		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", short read, many allowed gaps)...";
		doTestCase2(
			al,
			"A",
			"I",
			"AAAAAAAAAAAA",
			i*4,                // off
			sc,                 // scoring scheme
			-150.0f,            // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			false,              // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;

		if(i == 0) {
			cerr << "  Test " << tests++
			     << " (nuc space, offset 0, short read & ref, "
				 << "many allowed gaps)...";
			doTestCase2(
				al,
				"A",
				"I",
				"A",
				0,                  // off
				sc,                 // scoring scheme
				-150.0f,            // const coeff for cost ceiling
				0.0f,               // linear coeff for cost ceiling
				res,                // result
				false,              // color
				nIncl,              // Ns inclusive (not mismatches)
				nfilter);           // filter Ns
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), 0);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
			cerr << "PASSED" << endl;
		}
	}

	for(int i = 0; i < 3; i++) {
		if(i == 0) {
			cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
				 << ", exact)...";
			doTestCase2(
				al,
				"131",
				"III",
				"ACGT",
				i*4,                // off
				sc,                 // scoring scheme
				-30.0f,             // const coeff for cost ceiling
				0.0f,               // linear coeff for cost ceiling
				res,                // result
				true,               // color
				nIncl,              // Ns inclusive (not mismatches)
				nfilter);           // filter Ns
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), 0);
			assert(res.alres.ned().empty()); assert(res.alres.ced().empty());
			assert(res.alres.aed().empty());
			res.reset();
			cerr << "PASSED" << endl;
		}

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", exact)...";
		sc.rfGapConst = 40;
		sc.rdGapConst = 40;
		doTestCase2(
			al,
			"1313131",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert(res.alres.ned().empty()); assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", exact2)...";
		sc.rfGapConst = 40;
		sc.rdGapConst = 40;
		doTestCase2(
			al,
			"1313131",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", exact3)...";
		sc.rfGapConst = 50;
		sc.rdGapConst = 40;
		doTestCase2(
			al,
			"131313131",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-80.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", exact4)...";
		sc.rfGapConst = 40;
		sc.rdGapConst = 50;
		doTestCase2(
			al,
			"131313131",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-80.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty()); assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", 1 N allowed by minsc)...";
		sc.rfGapConst = 40;
		sc.rdGapConst = 40;
		sc.snp = 30;
		doTestCase3(
			al,
			"131.131",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			1.0f,               // allow 1 N
			0.0f,               // allow 1 N
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -1);
		assert_eq(1, res.alres.score().ns());
		assert_eq(0, res.alres.ned().size());
		assert_eq(1, res.alres.ced().size());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", 1 ref N and 1 read N)...";
		sc.rfGapConst = 40;
		sc.rdGapConst = 40;
		sc.snp = 30;
		doTestCase3(
			al,
			".313131",
			"IIIIIII",
			"ACGTACNTACGTACNT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			1.0f,               // allow 2 Ns
			0.2f,               // allow 2 Ns
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -2);
		assert_eq(res.alres.score().ns(), 2);
		assert_eq(2, res.alres.score().ns());
		assert_eq(1, res.alres.ned().size());
		assert_eq(1, res.alres.ced().size());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", 2 Nn allowed by minsc)...";
		sc.rfGapConst = 40;
		sc.rdGapConst = 40;
		sc.snp = 30;
		doTestCase2(
			al,
			".31313.",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -2);
		assert_eq(2, res.alres.score().ns());
		assert_eq(0, res.alres.ned().size());
		assert_eq(2, res.alres.ced().size());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		if(i == 0) {
			cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
			     << ", 2 Ns allowed by minsc)...";
			sc.rfGapConst = 40;
			sc.rdGapConst = 40;
			sc.snp = 30;
			doTestCase2(
				al,
				"..13131",
				"IIIIIII",
				"ANGTACGT",
				i*4,                // off
				sc,                 // scoring scheme
				-30.0f,             // const coeff for cost ceiling
				0.0f,               // linear coeff for cost ceiling
				res,                // result
				true,               // color
				nIncl,              // Ns inclusive (not mismatches)
				nfilter);           // filter Ns
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert_eq(res.alres.score().gaps(), 0);
			// TODO: should be -2
			assert_eq(res.alres.score().score(), -3);
			// 2 or 3?
			assert_eq(3, res.alres.score().ns());
			assert_eq(1, res.alres.ned().size());
			assert_eq(2, res.alres.ced().size());
			assert(res.alres.aed().empty());
			res.reset();
			cerr << "PASSED" << endl;
		}

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", 1 SNP allowed by minsc)...";
		sc.rfGapConst = 40;
		sc.rdGapConst = 40;
		sc.snp = 30;
		doTestCase2(
			al,
			"1310231",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -30);
		assert_eq(1, res.alres.ned().size());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", 1 SNP allowed by minsc)...";
		sc.rfGapConst = 40;
		sc.rdGapConst = 40;
		sc.snp = 30;
		doTestCase2(
			al,
			"0213131",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -30);
		assert_eq(1, res.alres.ned().size());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", 1 SNP allowed by minsc 2)...";
		sc.rfGapConst = 20;
		sc.rdGapConst = 20;
		sc.snp = 20;
		doTestCase2(
			al,
			"1310231",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-20.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -20);
		assert_eq(1, res.alres.ned().size());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", 1 MM allowed by minsc)...";
		sc.rfGapConst = 40;
		sc.rdGapConst = 40;
		sc.snp = 30;
		doTestCase2(
			al,
			"1310131",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -30);
		assert_eq(0, res.alres.ned().size());
		assert_eq(1, res.alres.ced().size());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", 1 MM allowed by minsc 2)...";
		sc.snp = 20;
		doTestCase2(
			al,
			"1310131",
			"5555555",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc2,                // scoring scheme
			-20.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -20);
		assert_eq(0, res.alres.ned().size());
		assert_eq(1, res.alres.ced().size());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", 1 read gap)...";
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 14;
		sc.rdGapLinear = 14;
		sc.snp = 30;
		doTestCase2(
			al,
			"131231",
			"IIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -39);
		assert_eq(1, res.alres.ned().size());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", 1 ref gap)...";
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 14;
		sc.rdGapLinear = 14;
		sc.snp = 30;
		doTestCase2(
			al,
			"13130131",
			"IIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -39);
		assert_eq(1, res.alres.ned().size());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4)
		     << ", 2 ref gaps)...";
		sc.rfGapConst  = 25;
		sc.rdGapConst  = 25;
		sc.rfGapLinear = 15;
		sc.rdGapLinear = 15;
		sc.snp = 30;
		doTestCase2(
			al,
			"131003131",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-60.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			true,               // color
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 2);
		assert_eq(res.alres.score().score(), -55);
		assert_eq(2, res.alres.ned().size());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;
	}

	// A test case where a valid alignment with a worse score should be
	// accepted over a valid alignment with a better score but too many
	// Ns
	cerr << "  Test " << tests++ << " (N ceiling 1)...";
	sc.mmcostType = COST_MODEL_CONSTANT;
	sc.mmcost = 10;
	sc.snp = 30;
	sc.nCeilConst  = 0.0f;
	sc.nCeilLinear = 0.0f;
	sc.rfGapConst  = 10;
	sc.rdGapLinear = 10;
	sc.rfGapConst  = 10;
	sc.rfGapLinear = 10;
	sc.setNPen(COST_MODEL_CONSTANT, 2);
	sc.gapbar = 1;
	// No Ns allowed, so this hit should be filtered
	doTestCase2(
		al,
		"ACGTACGT", // read seq
		"IIIIIIII", // read quals
		"NCGTACGT", // ref seq
		0,          // offset
		sc,         // scoring scheme
		-25.0f,     // const coeff for cost ceiling
		0.0f,       // linear coeff for cost ceiling
		res,        // result
		false,      // colorspace
		false,      // ns are in inclusive
		true,       // nfilter
		0);
	al.nextAlignment(res, rnd);
	assert(res.empty());
	cerr << "PASSED" << endl;
	res.reset();

	// No Ns allowed, so this hit should be filtered
	cerr << "  Test " << tests++ << " (N ceiling 2)...";
	doTestCase3(
		al,
		"ACGTACGT", // read seq
		"IIIIIIII", // read quals
		"NCGTACGT", // ref seq
		0,          // offset
		sc,         // scoring scheme
		-25.0f,     // const coeff for cost ceiling
		0.0f,       // linear coeff for cost ceiling
		1.0f,       // constant coefficient for # Ns allowed
		0.0f,       // linear coefficient for # Ns allowed
		res,        // result
		false,      // colorspace
		false,      // ns are in inclusive
		false,      // nfilter - NOTE: FILTER OFF
		0);
	al.nextAlignment(res, rnd);
	assert(!res.empty());
	assert_eq(0,  res.alres.score().gaps());
	assert_eq(-2, res.alres.score().score());
	assert_eq(1,  res.alres.score().ns());
	cerr << "PASSED" << endl;
	res.reset();

	// No Ns allowed, so this hit should be filtered
	cerr << "  Test " << tests++ << " (N ceiling 3)...";
	sc.nCeilConst = 1.0f;
	sc.nCeilLinear = 0.0f;
	doTestCase2(
		al,
		"ACGTACGT", // read seq
		"IIIIIIII", // read quals
		"NCGTACGT", // ref seq
		0,          // offset
		sc,         // scoring scheme
		-25.0f,     // const coeff for cost ceiling
		0.0f,       // linear coeff for cost ceiling
		res,        // result
		false,      // colorspace
		false,      // ns are in inclusive
		true,       // nfilter - NOTE: FILTER ON
		0);
	al.nextAlignment(res, rnd);
	assert(!res.empty());
	assert_eq(0,  res.alres.score().gaps());
	assert_eq(-2, res.alres.score().score());
	assert_eq(1,  res.alres.score().ns());
	cerr << "PASSED" << endl;
	res.reset();

	// No Ns allowed, so this hit should be filtered
	cerr << "  Test " << tests++ << " (redundant alignment elimination 1)...";
	sc.nCeilConst = 1.0f;
	sc.nCeilLinear = 0.0f;
	sc.rfGapConst  = 25;
	sc.rdGapLinear = 15;
	sc.rfGapConst  = 25;
	sc.rfGapLinear = 15;
	doTestCase2(
		al,
		//                   1         2         3         4
		//         01234567890123456789012345678901234567890123456
		          "AGGCTATGCCTCTGACGCGATATCGGCGCCCACTTCAGAGCTAACCG",
		          "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
		  "TTTTTTTTAGGCTATGCCTCTGACGCGATATCGGCGCCCACTTCAGAGCTAACCGTTTTTTT",
		// 01234567890123456789012345678901234567890123456789012345678901
		//           1         2         3         4         5         6
		8,          // offset
		sc,         // scoring scheme
		-25.0f,     // const coeff for cost ceiling
		-5.0f,      // linear coeff for cost ceiling
		res,        // result
		false,      // colorspace
		false,      // ns are in inclusive
		true,       // nfilter - NOTE: FILTER ON
		0);
	al.nextAlignment(res, rnd);
	assert(!res.empty());
	assert_eq(8, res.alres.refoff());
	assert_eq(47, res.alres.refExtent());
	assert_eq(0, res.alres.score().gaps());
	assert_eq(0, res.alres.score().score());
	assert_eq(0, res.alres.score().ns());
	res.reset();
	al.nextAlignment(res, rnd);
	assert(res.empty());
	assert(al.done());
	cerr << "PASSED" << endl;
	res.reset();
	
}

int main(int argc, char **argv) {
	int option_index = 0;
	int next_option;
	unsigned seed = 0;
	gGapBarrier = 1;
	gSnpPhred = 30;
	gColor = false;
	gColorExEnds = true;
	bool nsInclusive = false;
	bool nfilter = false;
	bonusMatchType  = DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch      = DEFAULT_MATCH_BONUS;
	penMmcType      = DEFAULT_MM_PENALTY_TYPE;
	penMmc          = DEFAULT_MM_PENALTY;
	penSnp          = DEFAULT_SNP_PENALTY;
	penNType        = DEFAULT_N_PENALTY_TYPE;
	penN            = DEFAULT_N_PENALTY;
	penRdExConst    = DEFAULT_READ_GAP_CONST;
	penRfExConst    = DEFAULT_REF_GAP_CONST;
	penRdExLinear   = DEFAULT_READ_GAP_LINEAR;
	penRfExLinear   = DEFAULT_REF_GAP_LINEAR;
	costMinConst    = DEFAULT_MIN_CONST;
	costMinLinear   = DEFAULT_MIN_LINEAR;
	costFloorConst  = DEFAULT_FLOOR_CONST;
	costFloorLinear = DEFAULT_FLOOR_LINEAR;
	nCeilConst      = 1.0f; // constant factor in N ceiling w/r/t read length
	nCeilLinear     = 1.0f; // coeff of linear term in N ceiling w/r/t read length
	nCatPair        = false;
	multiseedMms    = DEFAULT_SEEDMMS;
	multiseedLen    = DEFAULT_SEEDLEN;
	multiseedPeriod = DEFAULT_SEEDPERIOD;
	multiseedIvalType = DEFAULT_IVAL;
	multiseedIvalA    = DEFAULT_IVAL_A;
	multiseedIvalB    = DEFAULT_IVAL_B;
	posmin          = DEFAULT_POSMIN;
	posfrac         = DEFAULT_POSFRAC;
	rowmult         = DEFAULT_ROWMULT;
	do {
		next_option = getopt_long(argc, argv, short_opts, long_opts, &option_index);
		switch (next_option) {
			case 'C': gColor     = true; break;
			case 's': gSnpPhred  = parse<int>(optarg); break;
			case 'r': seed       = parse<unsigned>(optarg); break;
			case ARG_TESTS: {
				doTests();
				cout << "PASSED" << endl;
				return 0;
			}
			case 'A': {
				bool localAlign = false;
				bool noisyHpolymer = false;
				SeedAlignmentPolicy::parseString(
					optarg,
					localAlign,
					noisyHpolymer,
					bonusMatchType,
					bonusMatch,
					penMmcType,
					penMmc,
					penSnp,
					penNType,
					penN,
					penRdExConst,
					penRfExConst,
					penRdExLinear,
					penRfExLinear,
					costMinConst,
					costMinLinear,
					costFloorConst,
					costFloorLinear,
					nCeilConst,
					nCeilLinear,
					nCatPair,
					multiseedMms,
					multiseedLen,
					multiseedPeriod,
					multiseedIvalType,
					multiseedIvalA,
					multiseedIvalB,
					posmin,
					posfrac,
					rowmult);
				break;
			}
			case -1: break;
			default: {
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
			}
		}
	} while(next_option != -1);
	srand(seed);
	if(argc - optind < 4) {
		cerr << "Need at least 4 arguments" << endl;
		printUsage(cerr);
		exit(1);
	}
	BTDnaString read;
	BTString ref, qual;
	// Get read
	if(isalpha(argv[optind][0])) {
		read.installChars(argv[optind]);
	} else {
		read.installColors(argv[optind]);
	}
	// Get qualities
	qual.install(argv[optind+1]);
	assert_eq(read.length(), qual.length());
	// Get reference
	ref.install(argv[optind+2]);
	// Get reference offset
	size_t off = parse<size_t>(argv[optind+3]);
	// Set up penalities
	Scoring sc(
		false,         // local alignment?
		false,         // bad homopolymer?
		bonusMatchType,
		bonusMatch,
		penMmcType,    // how to penalize mismatches
		penMmc,        // constant if mm pelanty is a constant
		penSnp,        // penalty for nucleotide mismatch in decoded color alns
		costMinConst,
		costMinLinear,
		costFloorConst,
		costFloorLinear,
		nCeilConst,    // N ceiling constant coefficient
		nCeilLinear,   // N ceiling linear coefficient
		penNType,      // how to penalize Ns in the read
		penN,          // constant if N pelanty is a constant
		nCatPair,      // true -> concatenate mates before N filtering
		penRdExConst,  // constant cost of extending a gap in the read
		penRfExConst,  // constant cost of extending a gap in the reference
		penRdExLinear, // coeff of linear term for cost of gap extension in read
		penRfExLinear  // coeff of linear term for cost of gap extension in ref
	);
	// Calculate the penalty ceiling for the read
	TAlScore minsc = Scoring::linearFunc(
		read.length(),
		costMinConst,
		costMinLinear);
	TAlScore floorsc = Scoring::linearFunc(
		read.length(),
		costFloorConst,
		costFloorLinear);
	SwResult res;
	SwAligner al;
	doTestCase(
		al,
		read,
		qual,
		ref,
		off,
		NULL,
		NULL,
		sc,  
		minsc,
		floorsc,
		res,
		gColor,
		nsInclusive,
		nfilter,
		seed);
}
#endif /*MAIN_ALIGNER_SW*/
