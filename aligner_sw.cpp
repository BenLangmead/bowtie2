/*
 *  aligner_sw.cpp
 */

#include "aligner_sw.h"
#include "search_globals.h"
#include "penalty.h"
#include "mask.h"

/**
 * Initialize fields of SwParams by copying values of globals (see
 * search_globals.h).
 */
void SwParams::initFromGlobals() {
	gapBar = gGapBarrier;
	exEnds = gColorExEnds;
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
	const SwParams& pa,    // dynamic programming parameters
	const Penalties& pen,  // penalty scheme
	int penceil,           // max penalty we can incur for a valid alignment
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
		pa,          // dynamic programming parameters
		pen,         // penalties
		penceil);    // max penalty
	
	filter(nceil);
}

/**
 * Align read 'rd' to reference using read & reference information given
 * last time init() was called.  If the read is colorspace, the decoding is
 * determined simultaneously with alignment.  Uses dynamic programming.
 */
bool SwAligner::align(
	SwResult& res,
	RandomSource& rnd)
{
	assert(inited());
	int off = -1;
	if(color_) {
		off = alignColors(res, rnd);
	} else {
		off = alignNucleotides(res, rnd);
	}
	assert(off == -1 || !res.empty());
	assert(off != -1 ||  res.empty());
	if(off >= 0) {
		// Calculate the number of reference characters covered by the
		// alignment
		size_t extent = rdf_ - rdi_;
		const EList<Edit>& ned = res.alres.ned();
		for(size_t i = 0; i < ned.size(); i++) {
			if     (ned[i].isInsert()) extent++;
			else if(ned[i].isDelete()) extent--;
		}
		// Set the reference idx
		res.alres.setCoord(refidx_, off + rfi_ + refoff_, fw_, extent);
		assert(!res.alres.empty());
		assert(res.repOk());
		if(!fw_) {
			// All edits are currently w/r/t upstream end; if read aligned
			// to Crick strand, we need to invert them so that they're
			// w/r/t the read's 5' end instead.
			res.alres.invertEdits(rd_->length());
		}
		if(color_) {
			assert_range(0, 3, res.nup);
			assert_range(0, 3, res.ndn);
			res.alres.setNucs(fw_, res.nup, res.ndn);
		}
	}
	return off >= 0;
}

/**
 * Select a path for backtracking from this cell.  If there is a
 * tie among eligible paths, break it randomly.  Return value is
 * a pair where first = a flag indicating the backtrack type (see
 * enum defining SW_BT_* above), and second = a selection for
 * the read character for the next row up.  second should be
 * ignored if the backtrack type is a gap in the read.
 */
int SwNucCellMask::randBacktrack(RandomSource& rand) {
	ASSERT_ONLY(int num = numPossible());
	int i = ((diag != 0) << 0) |
			((rfop != 0) << 1) |
			((rfex != 0) << 2) |
			((rdop != 0) << 3) |
			((rdex != 0) << 4);
	int ret = randFromMask(rand, i);
	assert_range(0, 4, ret);
	switch(ret) {
		case SW_BT_DIAG:        diag = 0; break;
		case SW_BT_REF_OPEN:    rfop = 0; break;
		case SW_BT_REF_EXTEND:  rfex = 0; break;
		case SW_BT_READ_OPEN:   rdop = 0; break;
		case SW_BT_READ_EXTEND: rdex = 0; break;
		default: throw 1; break;
	}
	assert_eq(num-1, numPossible());
	return ret;
}

/**
 * Given the dynamic programming table, trace backwards from the given
 * cell and install the edits and the score/penalty in the appropriate
 * fields of res.
 *
 * The return value is the alignment's upstream-most reference
 * character's offset with respect to rfi.
 */
int SwAligner::backtrackNucleotides(
	AlignmentScore escore, // score we expect to get over backtrack
	SwResult&      res,    // store results (edits and scores) here
	int            col,    // start in this column (w/r/t the full matrix)
	RandomSource&  rand)   // pseudo-random generator
{
	ELList<SwNucCell>& tab = ntab_;
	int row = (int)rd_->length()-1;
	assert_eq(row, (int)tab.size()-1);
	AlignmentScore score; score.score_ = 0;
	score.gaps_ = score.ns_ = 0;
	bool refExtend = false, readExtend = false;
	ASSERT_ONLY(int origCol = col);
	ASSERT_ONLY(int gaps = 0);
	res.alres.reset();
	EList<Edit>& ned = res.alres.ned();
	while(row > 0) {
		res.swbts++;
		assert_range(0, (int)tab.size()-1, row);
		assert_range(0, origCol, col);
		int tabcol = col - row;
		assert_geq(tabcol, 0);
		ASSERT_ONLY(AlignmentScore scoreThisRound = score);
		ASSERT_ONLY(AlignmentScore bestThisRound = tab[row][tabcol].best);
		assert_gt(tab[row][tabcol].mask.numPossible(), 0);
		int cur = tab[row][tabcol].mask.randBacktrack(rand);
		switch(cur) {
			// Move up and to the left.  If the reference nucleotide in the
			// source row mismatches the read nucleotide, penalize
			// it and add a nucleotide mismatch.
			case SW_BT_DIAG: {
				refExtend = readExtend = false;
				assert_gt(row, 0); assert_gt(col, 0);
				// Check for color mismatch
				int readC = (*rd_)[row];
				int refNmask = (int)rf_[rfi_+col];
				int m = matchesEx(readC, refNmask);
				if(m != 1) {
					Edit e(row, mask2dna[refNmask], "ACGTN"[readC], EDIT_TYPE_MM);
					assert(e.repOk());
					ned.push_back(e);
					score.score_ -= QUAL2(row, col);
				}
				if(m == -1) {
					score.ns_++;
				}
				row--; col--;
				assert_range(0, (int)tab.size()-1, row);
				assert(VALID_AL_SCORE(score));
				ASSERT_ONLY(scoreThisRound = score - scoreThisRound);
				ASSERT_ONLY(bestThisRound -= tab[row][tabcol].best);
				// Make sure that both changed in the same way
				assert_eq(scoreThisRound, bestThisRound);
				break;
			}
			// Move up.  Add an edit encoding the ref gap.
			case SW_BT_REF_OPEN: {
				refExtend = true; readExtend = false;
				assert_gt(row, 0);
				Edit e(row, '-', "ACGTN"[(int)(*rd_)[row]], EDIT_TYPE_DEL);
				assert(e.repOk());
				ned.push_back(e);
				assert_geq(row, pa_->gapBar);
				assert_geq((int)(rdf_-rdi_-row-1), pa_->gapBar-1);
				row--;
				score.score_ -= pen_->refOpen;
				score.gaps_++;
#ifndef NDEBUG
				gaps++;
				assert_leq(score.gaps_, rdgap_ + rfgap_);
				assert_range(0, (int)tab.size()-1, row);
				tabcol = col - row;
				assert_range(0, (int)tab[row].size()-1, tabcol);
				scoreThisRound = score - scoreThisRound;
				bestThisRound -= tab[row][tabcol].best;
				assert_eq(scoreThisRound, bestThisRound);
#endif
				break;
			}
			// Move up.  Add an edit encoding the ref gap.
			case SW_BT_REF_EXTEND: {
				refExtend = true; readExtend = false;
				assert_gt(row, 1);
				Edit e(row, '-', "ACGTN"[(int)(*rd_)[row]], EDIT_TYPE_DEL);
				assert(e.repOk());
				ned.push_back(e);
				assert_geq(row, pa_->gapBar);
				assert_geq((int)(rdf_-rdi_-row-1), pa_->gapBar-1);
				row--;
				score.score_ -= pen_->refExConst;
				score.gaps_++;
#ifndef NDEBUG
				gaps++;
				assert_leq(score.gaps_, rdgap_ + rfgap_);
				assert_range(0, (int)tab.size()-1, row);
				tabcol = col - row;
				assert_range(0, (int)tab[row].size()-1, tabcol);
				scoreThisRound = score - scoreThisRound;
				bestThisRound -= tab[row][tabcol].best;
				assert_eq(scoreThisRound, bestThisRound);
#endif
				break;
			}
			case SW_BT_READ_OPEN: {
				refExtend = false; readExtend = true;
				assert_gt(col, 0);
				Edit e(row+1, mask2dna[(int)rf_[rfi_+col]], '-', EDIT_TYPE_INS);
				assert(e.repOk());
				ned.push_back(e);
				assert_geq(row, pa_->gapBar);
				assert_geq((int)(rdf_-rdi_-row-1), pa_->gapBar-1);
				col--;
				score.score_ -= pen_->readOpen;
				score.gaps_++;
#ifndef NDEBUG
				gaps++;
				assert_leq(score.gaps_, rdgap_ + rfgap_);
				assert_range(0, (int)tab.size()-1, row);
				tabcol = col - row;
				assert_range(0, (int)tab[row].size()-1, tabcol);
				scoreThisRound = score - scoreThisRound;
				bestThisRound -= tab[row][tabcol].best;
				assert_eq(scoreThisRound, bestThisRound);
#endif
				break;
			}
			case SW_BT_READ_EXTEND: {
				refExtend = false; readExtend = true;
				assert_gt(col, 1);
				Edit e(row+1, mask2dna[(int)rf_[rfi_+col]], '-', EDIT_TYPE_INS);
				assert(e.repOk());
				ned.push_back(e);
				assert_geq(row, pa_->gapBar);
				assert_geq((int)(rdf_-rdi_-row-1), pa_->gapBar-1);
				col--;
				score.score_ -= pen_->readExConst;
				score.gaps_++;
#ifndef NDEBUG
				gaps++;
				assert_leq(score.gaps_, rdgap_ + rfgap_);
				assert_range(0, (int)tab.size()-1, row);
				tabcol = col - row;
				assert_range(0, (int)tab[row].size()-1, tabcol);
				scoreThisRound = score - scoreThisRound;
				bestThisRound -= tab[row][tabcol].best;
				assert_eq(scoreThisRound, bestThisRound);
#endif
				break;
			}
			default: throw 1;
		}
	}
	assert_eq(0, row);
	int readC = (*rd_)[rdi_+row];         // get last char in read
	int refNmask = (int)rf_[rfi_+col]; // get last char in ref involved in alignment
	int m = matchesEx(readC, refNmask);
	if(m != 1) {
		Edit e(row, mask2dna[refNmask], "ACGTN"[readC], EDIT_TYPE_MM);
		assert(e.repOk());
		ned.push_back(e);
		score.score_ -= QUAL2(row, col);
	}
	if(m == -1) {
		score.ns_++;
	}
	res.reverse();
	assert(Edit::repOk(ned, (*rd_)));
#ifndef NDEBUG
	BTDnaString refstr;
	for(int i = col; i <= origCol; i++) {
		refstr.append(firsts5[(int)rf_[rfi_+i]]);
	}
	BTDnaString editstr;
	Edit::toRef((*rd_), ned, editstr);
	if(refstr != editstr) {
		cerr << "Decoded nucleotides and edits don't match reference:" << endl;
		cerr << "           score: " << score.score() << " (" << gaps << " gaps)" << endl;
		cerr << "           edits: ";
		Edit::print(cerr, ned);
		cerr << endl;
		cerr << "    decoded nucs: " << (*rd_) << endl;
		cerr << "     edited nucs: " << editstr << endl;
		cerr << "  reference nucs: " << refstr << endl;
		assert(0);
	}
#endif
	// done
	assert_eq(score.score(), escore.score());
	//assert_leq(score.gaps, escore.gaps);
	assert_leq(gaps, rdgap_ + rfgap_);
	// Dummy values for refid and fw
	res.alres.setScore(score);
	return col;
}

/**
 * Update a SwCell's best[] and mask[] arrays with respect to its
 * neighbor on the left.
 */
inline void SwAligner::updateNucHoriz(
	const SwNucCell& lc,
	SwNucCell&       dstc,
	int              rfm)
{
	assert(lc.finalized);
	if(lc.empty) return;
	bool ninvolved = rfm > 15;
	{
		AlignmentScore leftBest = lc.best;
		const SwNucCellMask& frMask = lc.mask;
		AlignmentScore& myBest = dstc.best;
		SwNucCellMask& myMask = dstc.mask;
		assert_leq(leftBest.score(), 0);
		if(ninvolved) leftBest.incNs(nceil_);
		if(!VALID_AL_SCORE(leftBest)) return;
		// *Don't* penalize for a nucleotide mismatch because we must
		// have already done that in a previous vertical or diagonal
		// step.
		if(frMask.readExtendPossible()) {
			// Read gap extension possible?
			AlignmentScore ex = leftBest;
			assert_leq(ex.score(), 0);
			assert(VALID_AL_SCORE(ex));
			ex.score_ -= pen_->readExConst;
			assert_leq(ex.score(), 0);
			if(-ex.score_ <= penceil_ && ex >= myBest) {
				if(ex > myBest) {
					myMask.clear();
					myBest = ex;
				}
				myMask.rdex = 1;
				assert(VALID_AL_SCORE(myBest));
			}
		}
		if(frMask.readOpenPossible()) {
			// Read gap open possible?
			AlignmentScore ex = leftBest;
			assert_leq(ex.score_, 0);
			assert(VALID_AL_SCORE(ex));
			ex.score_ -= pen_->readOpen;
			assert_leq(ex.score_, 0);
			if(-ex.score_ <= penceil_ && ex >= myBest) {
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
	{
		// Assuming that the read character in this row is 'to'...
		// No SNP penalty because destination read char aligns to a
		// gap in the reference.
		{
			AlignmentScore from = uc.best;
			assert(!uc.empty || !VALID_AL_SCORE(from));
			const SwNucCellMask& frMask = uc.mask;
			AlignmentScore& myBest = dstc.best;
			SwNucCellMask& myMask = dstc.mask;
			if(rdc > 3) {
				from.incNs(nceil_);
			}
			if(!VALID_AL_SCORE(from)) return;
			assert_leq(from.score_, 0);
			if(frMask.refExtendPossible()) {
				// Extend is possible
				from.score_ -= pen_->refExConst;
				assert_leq(from.score_, 0);
				if(-from.score_ <= penceil_ && from >= myBest) {
					if(from > myBest) {
						myBest = from;
						myMask.clear();
					}
					myMask.rfex = 1;
					assert(VALID_AL_SCORE(myBest));
				}
				// put it back
				from.score_ += pen_->refExConst;
			}
			if(frMask.refOpenPossible()){
				// Open is possible
				from.score_ -= pen_->refOpen;
				if(-from.score_ <= penceil_ && from >= myBest) {
					if(from > myBest) {
						myBest = from;
						myMask.clear();
					}
					myMask.rfop = 1;
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
inline void SwAligner::updateNucDiag(
	const SwNucCell& dc,
	SwNucCell&       dstc,
	int              rdc,
	int              rfm,
	int              qpen)
{
	assert(dc.finalized);
	if(dc.empty) return;
	bool ninvolved = (rdc > 3 || rfm > 15);
	int add = (matches(rdc, rfm) ? 0 : (ninvolved ? pen_->n(30) : qpen));
	AlignmentScore from = dc.best - add;
	{
		{
			AlignmentScore& myBest = dstc.best;
			SwNucCellMask& myMask = dstc.mask;
			if(ninvolved) from.incNs(nceil_);
			if(!VALID_AL_SCORE(from)) return;
			assert_leq(from.score_, 0);
			if(-from.score_ <= penceil_ && from >= myBest) {
				if(from > myBest) {
					myBest = from;
					myMask.clear();
					assert_eq(0, myMask.diag);
				}
				myMask.diag = 1;
				assert(VALID_AL_SCORE(myBest));
			}
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
	assert_eq(rff_ - rfi_, width_ + nrow - 1);
	// For each reference character involved
	for(size_t i = rfi_; i <= rff_; i++) {
		int rfc = rf_[i];
		assert_range(0, 16, rfc);
		if(i >= nrow) {
			size_t col = i - nrow;
			if(ns > nlim) {
				if((*en_)[col]) {
					(*en_)[col] = false;
					filtered++;
				}
			}
			int prev_rfc = rf_[i - nrow];
			assert_range(0, 16, prev_rfc);
			if(prev_rfc == 16) {
				assert_gt(ns, 0);
				ns--;
			}
		}
		if(rfc == 16) {
			ns++;
		}
	}
	return filtered;
}

/**
 * Align the nucleotide read 'read' as aligned against the reference
 * string 'rf' using a banded dynamic programming algorithm.
 *
 * If an alignment is found, its offset relative to rdi is returned.
 * E.g. if an alignment is found that occurs starting at rdi, 0 is
 * returned.  If no alignment is found, -1 is returned.
 */
int SwAligner::alignNucleotides(
	SwResult& res,
	RandomSource& rnd)
{
	typedef SwNucCell TCell;
	assert_leq(rdf_, rd_->length());
	assert_leq(rdf_, qu_->length());
	assert_lt(rfi_, rff_);
	assert_lt(rdi_, rdf_);
	assert_eq(rd_->length(), qu_->length());
	assert_geq(pa_->gapBar, 1);
	res.sws++;
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
	const int wlo = 0;
	const int whi = (int)(width_ - 1);
	assert_lt(whi, (int)(rff_-rfi_));
	tab[0].resize(whi-wlo+1); // add columns to first row
	bool validInRow = false;
	// Calculate starting values for the rest of the columns in the
	// first row.
	for(int col = 0; col <= whi; col++) {
		tab[0][col].clear(); // clear the cell; masks and scores
		int rdc = (*rd_)[rdi_+0];
		int rfm = rf_[rfi_+col];
		// Can we start from here?
		bool canStart = (st_ == NULL || (*st_)[col]);
		if(canStart) {
			tab[0][col].best.gaps_ = 0;
			tab[0][col].best.ns_ = 0;
			tab[0][col].best.score_ = 0;
			int m = matchesEx(rdc, rfm);
			if(m == 1) {
				// The assigned subject nucleotide matches the reference;
				// no penalty
				assert_lt(rdc, 4);
				assert_lt(rfm, 16);
				tab[0][col].mask.diag = 1;
			} else if(QUAL2(0, col) <= penceil_) {
				// Reference char mismatches
				int n = (rfm > 15 || rdc > 3) ? 1 : 0;
				tab[0][col].best.score_ -= QUAL2(0, col);
				tab[0][col].best.ns_ = n;
				tab[0][col].mask.diag = 1;
			} else {
				// Leave mask cleared
			}
		}
		// Calculate horizontals if barrier allows
		if(pa_->gapBar <= 1 && col > 0) {
			updateNucHoriz(
				tab[0][col-1],
				tab[0][col],
				rfm);
			res.swcups++;
		}
		assert(!tab[0][col].finalized);
		if(tab[0][col].finalize(penceil_)) validInRow = true;
	}
	res.swrows++;
	if(!validInRow) {
		res.swskiprows += (rdf_ - rdi_ - 1);
		assert(res.empty());
		return -1;
	}

	//
	// Calculate all subsequent rows
	//

	// Do rest of table
	for(int row = 1; row < (int)(rdf_-rdi_); row++) {
		res.swrows++;
		tab.expand(); // add another row
		bool onlyDiagInto = (row+1 <= pa_->gapBar || (int)(rdf_-rdi_)-row <= pa_->gapBar);
		tab.back().resize(whi-wlo+1); // add enough space for columns
		assert_range(1, (int)qu_->length()-1, row);
		assert_range(1, (int)rd_->length()-1, row);
		int c = (*rd_)[row];   // read character in this row
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
			const int fc = col + row;
			updateNucDiag(
				tab[row-1][0],     // cell diagonally above and to the left
				cur,               // destination cell
				c,                 // color being traversed
				rf_[rfi_ + fc], // ref mask at destination cell
				QUAL2(row, fc));   //
		}
		if(!onlyDiagInto && col < whi && !tab[row-1][1].empty) {
			updateNucVert(
				tab[row-1][1],     // cell diagonally above and to the left
				cur,               // destination cell
				c);                // color being traversed
		}
		res.swcups++;
		// 'cur' is now initialized
		assert(!cur.finalized);
		if(cur.finalize(penceil_)) validInRow = true;
		
		// Iterate from leftmost to rightmost inner diagonals
		for(col = wlo+1; col < whi; col++) {
			const int fullcol = col + row;
			int r = rf_[rfi_ + fullcol];
			TCell& cur = tab[row][col-wlo];
			cur.clear();
			TCell& dg = tab[row-1][col-wlo];
			// The mismatch penalty is a function of the read character
			// in this row & the reference character in this column
			// (specifically: whether they match and whether either is
			// an N) as well as the quality value of the read
			// character.
			const int mmpen = QUAL2(row, fullcol);
			updateNucDiag(
				dg,            // cell diagonally above and to the left
				cur,           // destination cell
				c,             // nucleotide in destination row
				r,             // ref mask associated with destination column
				mmpen);        // penalty to incur for color miscall
			if(!onlyDiagInto) {
				TCell& up = tab[row-1][col-wlo+1];
				updateNucVert(
					up,        // cell above
					cur,       // destination cell
					c);        // nucleotide in destination row
				// Can do horizontal
				TCell& lf = tab[row][col-wlo-1];
				updateNucHoriz(
					lf,        // cell to the left
					cur,       // destination cell
					r);        // ref mask associated with destination column
			} // end loop over inner diagonals
			res.swcups++;
			// 'cur' is now initialized
			assert(!cur.finalized);
			if(cur.finalize(penceil_)) validInRow = true;
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
			const int r = rf_[rfi_ + fullcol];
			const int mmpenf = QUAL2(row, fullcol);
			TCell& dg = tab[row-1][col-wlo];
			if(!dg.empty) {
				updateNucDiag(
					dg,        // cell diagonally above and to the left
					cur,       // destination cell
					c,         // nucleotide in destination row
					r,         // ref mask associated with destination column
					mmpenf);   // penalty to incur for color miscall
			}
			TCell& lf = tab[row][col-wlo-1];
			if(!onlyDiagInto && !lf.empty) {
				updateNucHoriz(
					lf,        // cell to the left
					cur,       // destination cell
					r);        // ref mask associated with destination column
			}
			res.swcups++;
			// 'cur' is now initialized
			assert(!cur.finalized);
			if(cur.finalize(penceil_)) validInRow = true;
		}
		if(!validInRow) {
			assert_geq((int)(rdf_-rdi_), row+1);
			res.swskiprows += (rdf_ - rdi_ - row - 1);
			assert(res.empty());
			return -1;
		}
	}
	assert_eq(tab.size(), rd_->length());
	// Go hunting for cell to backtrace from; i.e. best score in the
	// bottom row and in the last readGaps*2+1 columns.
	AlignmentScore bscore = AlignmentScore::INVALID();
	int btCol = -1; // column to backtrace from
	int lastRow = (int)(rdf_-rdi_-1);
	for(int col = wlo; col <= whi; col++) {
		// greater than or equal to???
		bool canEnd = (en_ == NULL || (*en_)[col-wlo]);
		if(canEnd && !tab[lastRow][col].empty) {
			assert(tab[lastRow][col].finalized);
			assert_leq(abs(tab[lastRow][col].best.score()), penceil_);
			if(tab[lastRow][col].updateBest(bscore, penceil_)) {
				assert(!VALID_AL_SCORE(bscore) || abs(bscore.score()) <= penceil_);
				btCol = col;
			}
		}
	}
	if(btCol == -1) {
		return -1;
	}
	//assert_range(wlo, whi, btCol);
	assert(VALID_AL_SCORE(bscore));
	assert_leq(abs(bscore.score()), penceil_);
	int off = backtrackNucleotides(
		bscore,        // score we expect to get over backtrack
		res,           // store results (edits and scores) here
		btCol+lastRow, // start in this column
		rnd);          // pseudo-random generator
	assert_geq(off, 0);
#ifndef NDEBUG
	EList<Edit>& ned = res.alres.ned();
	EList<Edit>& aed = res.alres.aed();
	assert(res.alres.ced().empty());
	for(int i = 0; i < (int)ned.size(); i++) {
		assert_lt(ned[i].pos, rd_->length());
	}
	for(int i = 0; i < (int)aed.size(); i++) {
		assert_lt(aed[i].pos, rd_->length());
	}
#endif
	res.alres.setColor(false);
	// Return offset of alignment with respect to rfi
	return off;
}

#ifdef MAIN_ALIGNER_SW

#include <sstream>
#include <utility>
#include <getopt.h>
#include "penalty.h"
#include "aligner_seed_policy.h"
#include "color.h"

int  gGapBarrier;
int  gSnpPhred;
bool gColor;
bool gColorExEnds;
static int penMmcType;       // how to penalize mismatches
static int penMmc;           // constant if mm pelanty is a constant
static int penSnp;           // penalty associated with nucleotide mismatch in decoded colorspace alignment
static int penNType;         // how to penalize Ns in the read
static int penN;             // constant if N pelanty is a constant
static bool nPairCat;        // true -> concatenate mates before N filter
static int penRdOpen;        // cost of opening a gap in the read
static int penRfOpen;        // cost of opening a gap in the reference
static int penRdExConst;     // constant cost of extending a gap in the read
static int penRfExConst;     // constant cost of extending a gap in the reference
static int penRdExLinear;    // coeff of linear term for cost of gap extension in read
static int penRfExLinear;    // coeff of linear term for cost of gap extension in ref
static float costCeilConst;  // constant factor in cost ceiling w/r/t read length
static float costCeilLinear; // coeff of linear term in cost ceiling w/r/t read length
static float nCeilConst;     // constant factor in N ceiling w/r/t read length
static float nCeilLinear;    // coeff of linear term in N ceiling w/r/t read length
static int multiseedMms;     // mismatches permitted in a multiseed seed
static int multiseedLen;     // length of multiseed seeds
static int multiseedPeriod;  // space between multiseed seeds
static int multiseedIvalType;
static float multiseedIvalA;
static float multiseedIvalB;

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
	os << "  -s/--snppen <int>   penalty incurred by SNP; used for decoding" << endl;
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

/**
 * Helper function for running a case consisting of a read (sequence
 * and quality), a reference string, and an offset that anchors the 0th
 * character of the read to a reference position.
 */
static void doTestCase(
	const BTDnaString& read,
	const BTString&    qual,
	const BTString&    refin,
	size_t             off,
	EList<bool>       *st,
	EList<bool>       *en,
	const SwParams&    pa,
	const Penalties&   pens,
	int                penceil,
	SwResult&          res,
	bool               color,
	bool               nsInclusive,
	bool               filterns,
	uint32_t           seed)
{
	SwAligner al;
	RandomSource rnd(seed);
	BTString ref = refin;
	assert_eq(read.length(), qual.length());
	size_t rfi = 0;
	size_t rff = ref.length();
	// Calculate the largest possible number of read and reference
	// gaps given 'penceil' and 'pens'
	int readGaps = pens.maxReadGaps(penceil);
	int refGaps = pens.maxRefGaps(penceil);
	assert_geq(readGaps, 0);
	assert_geq(refGaps, 0);
	int maxGaps = max(readGaps, refGaps);
	size_t nceil = pens.nCeil(read.length());
	size_t width = 1 + 2 * maxGaps;
	size_t refoff = off;
	// Pad the beginning of the reference with Ns if necessary
	while(refoff < maxGaps) {
		ref.insert('N', 0);
		refoff++;
	}
	rfi = refoff - maxGaps;
	// Pad the end of the reference with Ns if necessary
	int colori = color ? 1 : 0;
	while(refoff + read.length() + colori + maxGaps > ref.length()) {
		ref.append('N');
	}
	rff = refoff + read.length() + colori + maxGaps;
	// Convert reference string to masks
	for(size_t i = 0; i < ref.length(); i++) {
		if(toupper(ref[i]) == 'N' && !nsInclusive) {
			ref.set(16, i);
		} else {
			int num = 0;
			int alts[] = {4, 4, 4, 4};
			decodeNuc(toupper(ref[i]), num, alts);
			assert_leq(num, 4);
			assert_gt(num, 0);
			ref.set(0, i);
			for(int j = 0; j < num; j++) {
				ref.set(ref[i] | (1 << alts[j]), i);
			}
		}
	}
	bool fw = true;
	uint32_t refidx = 0;
	EList<bool> stbuf, enbuf;
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
	al.init(
		read,          // read sequence
		qual,          // read qualities
		0,             // offset of first character within 'read' to consider
		read.length(), // offset of last char (exclusive) in 'read' to consider
		fw,            // 'read' is forward version of read?
		color,         // whether read is nucleotide-space or colorspace
		refidx,        // id of reference aligned to
		off,           // offset of upstream ref char aligned against
		ref.wbuf(),    // reference sequence (masks)
		rfi,           // offset of first char in 'ref' to consider
		rff,           // offset of last char (exclusive) in 'ref' to consider
		width,         // # bands to do (width of parallelogram)
		st,            // mask indicating which columns we can start in
		en,            // mask indicating which columns we can end in
		pa,            // dynamic programming parameters
		pens,          // penalties
		penceil);      // max total penalty
	if(filterns) {
		al.filter((int)nceil);
	}
	al.align(
		res,
		rnd);
	if(!res.empty()) {
		cout << " Score: " << res.alres.score() << endl;
		if(color) {
			cout << "   Read colors: " << endl;
			cout << "     ";
			for(size_t i = 0; i < read.length(); i++) {
				printColor((int)read[i]);
			}
			cout << endl;
			cout << "   Color alignment (decoded): " << endl;
			Edit::printQAlign(cout, "     ", read, res.alres.ced());
		} else {
			cout << "   Nucleotide alignment: " << endl;
			Edit::printQAlign(cout, "     ", read, res.alres.ned());
		}
	} else {
		cout << "Empty result" << endl;
	}
}

/**
 * Another interface for running a case.
 */
static void doTestCase2(
	const char        *read,
	const char        *qual,
	const char        *refin,
	size_t             off,
	const SwParams&    pa,
	const Penalties&   pens,
	float              costCeilConst,
	float              costCeilLinear,
	SwResult&          res,
	bool               color,
	bool               nsInclusive = false,
	bool               filterns = false,
	uint32_t           seed = 0)
{
	BTDnaString btread(read, true, color);
	// Calculate the penalty ceiling for the read
	int penceil = Constraint::instantiate(
		btread.length(),
		costCeilConst,
		costCeilLinear);
	doTestCase(
		btread,
		BTString(qual),
		BTString(refin),
		off,
		NULL,
		NULL,
		pa,
		pens,
		penceil,
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
	penMmcType      = DEFAULT_MM_PENALTY_TYPE;
	penMmc          = DEFAULT_MM_PENALTY;
	penSnp          = DEFAULT_SNP_PENALTY;
	penNType        = DEFAULT_N_PENALTY_TYPE;
	penN            = DEFAULT_N_PENALTY;
	nPairCat        = DEFAULT_N_CAT_PAIR;
	penRdOpen       = DEFAULT_READ_OPEN;
	penRfOpen       = DEFAULT_REF_OPEN;
	penRdExConst    = DEFAULT_READ_EXTEND_CONST;
	penRfExConst    = DEFAULT_REF_EXTEND_CONST;
	penRdExLinear   = DEFAULT_READ_EXTEND_LINEAR;
	penRfExLinear   = DEFAULT_REF_EXTEND_LINEAR;
	costCeilConst   = DEFAULT_CEIL_CONST;
	costCeilLinear  = DEFAULT_CEIL_LINEAR;
	nCeilConst      = 1.0f; // constant factor in N ceiling w/r/t read length
	nCeilLinear     = 1.0f; // coeff of linear term in N ceiling w/r/t read length
	multiseedMms    = DEFAULT_SEEDMMS;
	multiseedLen    = DEFAULT_SEEDLEN;
	multiseedPeriod = DEFAULT_SEEDPERIOD;
	// Set up penalities
	Penalties pens(
		penMmcType,    // how to penalize mismatches
		penMmc,        // constant if mm pelanty is a constant
		penSnp,
		nCeilConst,    // constant factor in N ceiling w/r/t read length
		nCeilLinear,   // coeff of linear term in N ceiling w/r/t read length
		penNType,      // how to penalize Ns in the read
		penN,          // constant if N pelanty is a constant
		nPairCat,      // true -> concatenate mates before N filtering
		penRdOpen,     // cost of opening a gap in the read
		penRfOpen,     // cost of opening a gap in the reference
		penRdExConst,  // constant cost of extending a gap in the read
		penRfExConst,  // constant cost of extending a gap in the reference
		penRdExLinear, // coefficient of linear term for cost of gap extension in read
		penRfExLinear  // coefficient of linear term for cost of gap extension in ref
	);
	// Set up alternative penalities
	Penalties pens2(
		COST_MODEL_QUAL, // how to penalize mismatches
		penMmc,        // constant if mm pelanty is a constant
		penSnp,
		1.0f,
		1.0f,
		penNType,      // how to penalize Ns in the read
		penN,          // constant if N pelanty is a constant
		nPairCat,      // true -> concatenate mates before N filtering
		40,            // cost of opening a gap in the read
		40,            // cost of opening a gap in the reference
		penRdExConst,  // constant cost of extending a gap in the read
		penRfExConst,  // constant cost of extending a gap in the reference
		penRdExLinear, // coefficient of linear term for cost of gap extension in read
		penRfExLinear  // coefficient of linear term for cost of gap extension in ref
	);
	SwParams pa;
	SwResult res;
	
	//
	// Basic nucleotide-space tests
	//
	cerr << "Running tests..." << endl;
	int tests = 1;
	bool nIncl = false;
	bool nfilter = false;

	for(int i = 0; i < 3; i++) {
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", exact)...";
		pens.refOpen = 40;
		pens.readOpen = 40;
		doTestCase2(
			"ACGTACGT",
			"IIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,
			pa,
			pens,
			30.0f,
			0.0f,
			res,
			false,
			nIncl,
			nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty()); assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1mm allowed by penceil)...";
		doTestCase2("ACGTTCGT", "IIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -30);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1mm allowed by penceil, check qual 1)...";
		doTestCase2("ACGTTCGT", "ABCDEFGH", "ACGTACGTACGTACGT", i*4, pa, pens2, 40.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -36);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1mm allowed by penceil, check qual 2)...";
		doTestCase2("ACGAACGT", "ABCDEFGH", "ACGTACGTACGTACGT", i*4, pa, pens2, 40.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -35);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1mm allowed by penceil, check qual )...";
		doTestCase2("TCGTACGT", "ABCDEFGH", "ACGTACGTACGTACGT", i*4, pa, pens2, 40.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -32);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1mm at the beginning, allowed by penceil)...";
		doTestCase2("CCGTACGT", "IIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -30);
		assert_eq(res.alres.score().ns(), 0);
		assert_eq(1, res.alres.ned().size());
		assert_eq(0, res.alres.aed().size());
		assert_eq(0, res.alres.ced().size());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1 n in read, allowed)...";
		doTestCase2("ACGTNCGT", "IIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -1);
		assert_eq(res.alres.score().ns(), 1);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 2 n in read, allowed)...";
		doTestCase2("ACGNNCGT", "IIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -2);
		assert_eq(res.alres.score().ns(), 2);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 2 n in read, 1 at beginning, allowed)...";
		doTestCase2("NCGTNCGT", "IIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -2);
		assert_eq(res.alres.score().ns(), 2);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1 n in ref, allowed)...";
		doTestCase2("ACGTACGT", "IIIIIIII", "ACGTNCGTACGTANGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -1);
		assert_eq(res.alres.score().ns(), 1);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1mm disallowed by penceil)...";
		doTestCase2("ACGTTCGT", "IIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 10.0f, 0.0f, res, false, nIncl, nfilter);
		assert(res.empty());
		res.reset();
		// Read gap with equal read and ref gap penalties
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", read gap allowed by penceil)...";
		doTestCase2("ACGTCGT", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 40.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", read gap disallowed by penceil)...";
		doTestCase2("ACGTCGT", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		// Ref gap with equal read and ref gap penalties
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", ref gap allowed by penceil)...";
		doTestCase2("ACGTAACGT", "IIIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 40.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", ref gap disallowed by penceil)...";
		doTestCase2("ACGTAACGT", "IIIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		// Read gap with one read gap and zero ref gaps allowed
		pens.refOpen = 50;
		pens.readOpen = 40;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1 read gap, ref gaps disallowed by penceil)...";
		doTestCase2("ACGTCGT", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 40.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", gaps disallowed by penceil)...";
		doTestCase2("ACGTCGT", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		// Ref gap with one ref gap and zero read gaps allowed
		pens.refOpen = 40;
		pens.readOpen = 50;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1 ref gap, read gaps disallowed by penceil)...";
		doTestCase2("ACGTAACGT", "IIIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 40.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", gaps disallowed by penceil)...";
		doTestCase2("ACGTAACGT", "IIIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		// Read gap with one read gap and two ref gaps allowed
		pens.refOpen = 20;
		pens.readOpen = 40;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1 read gap, 2 ref gaps allowed by penceil)...";
		doTestCase2("ACGTCGT", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 40.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", gaps disallowed by penceil)...";
		doTestCase2("ACGTCGT", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(res.empty());
		res.reset();
		// Ref gap with one ref gap and two read gaps allowed
		pens.refOpen = 40;
		pens.readOpen = 20;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1 ref gap, 2 read gaps allowed by penceil)...";
		doTestCase2("ACGTAACGT", "IIIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 40.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", gaps disallowed by penceil)...";
		doTestCase2("ACGTAACGT", "IIIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		// Read gap with two read gaps and two ref gaps allowed
		pens.refOpen = 20;
		pens.readOpen = 20;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 2 ref gaps, 2 read gaps allowed by penceil)...";
		doTestCase2("ACGTCGT", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 40.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -20);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1 ref gap, 1 read gap allowed by penceil)...";
		doTestCase2("ACGTCGT", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -20);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		// Ref gap with two ref gaps and zero read gaps allowed
		pens.refOpen = 20;
		pens.readOpen = 20;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 2 ref gaps, 2 read gaps allowed by penceil)...";
		doTestCase2("ACGTAACGT", "IIIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 40.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -20);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", 1 ref gap, 1 read gap allowed by penceil)...";
		doTestCase2("ACGTAACGT", "IIIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -20);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", short read)...";
		doTestCase2("A", "I", "AAAAAAAAAAAA", i*4, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;

		if(i == 0) {
			cerr << "  Test " << tests++ << " (nuc space, offset 0, short read & ref)...";
			doTestCase2("A", "I", "A", 0, pa, pens, 30.0f, 0.0f, res, false, nIncl, nfilter);
			assert(!res.empty());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), 0);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
			cerr << "PASSED" << endl;
		}

		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", short read, many allowed gaps)...";
		doTestCase2("A", "I", "AAAAAAAAAAAA", i*4, pa, pens, 150.0f, 0.0f, res, false, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;

		if(i == 0) {
			cerr << "  Test " << tests++ << " (nuc space, offset 0, short read & ref, many allowed gaps)...";
			doTestCase2("A", "I", "A", 0, pa, pens, 150.0f, 0.0f, res, false, nIncl, nfilter);
			assert(!res.empty());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), 0);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
			cerr << "PASSED" << endl;
		}
	}

	for(int i = 0; i < 3; i++) {
		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", exact)...";
		pens.refOpen = 40;
		pens.readOpen = 40;
		pa.exEnds = false;
		doTestCase2("1313131", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, true, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert(res.alres.ned().empty()); assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", exact2)...";
		pens.refOpen = 40;
		pens.readOpen = 40;
		pa.exEnds = false;
		doTestCase2("1313131", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 40.0f, 0.0f, res, true, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", exact3)...";
		pens.refOpen = 50;
		pens.readOpen = 40;
		pa.exEnds = false;
		doTestCase2("131313131", "IIIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 80.0f, 0.0f, res, true, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0); assert_eq(res.alres.score().score(), 0);
		assert(res.alres.ned().empty()); assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", exact4)...";
		pens.refOpen = 40;
		pens.readOpen = 50;
		pa.exEnds = false;
		doTestCase2("131313131", "IIIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 80.0f, 0.0f, res, true, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty()); assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", 1 N allowed by penceil)...";
		pens.refOpen = 40;
		pens.readOpen = 40;
		pens.snp = 30;
		pa.exEnds = false;
		doTestCase2("131.131", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, true, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -1);
		assert_eq(1, res.alres.score().ns());
		assert_eq(0, res.alres.ned().size());
		assert_eq(1, res.alres.ced().size());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", 1 ref N and 1 read N)...";
		pens.refOpen = 40;
		pens.readOpen = 40;
		pens.snp = 30;
		pa.exEnds = false;
		doTestCase2(".313131", "IIIIIII", "ACGTACNTACGTACNT", i*4, pa, pens, 30.0f, 0.0f, res, true, nIncl, nfilter);
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

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", 2 Nn allowed by penceil)...";
		pens.refOpen = 40;
		pens.readOpen = 40;
		pens.snp = 30;
		pa.exEnds = false;
		doTestCase2(".31313.", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, true, nIncl, nfilter);
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
			cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", 2 Nn allowed by penceil)...";
			pens.refOpen = 40;
			pens.readOpen = 40;
			pens.snp = 30;
			pa.exEnds = false;
			doTestCase2("..13131", "IIIIIII", "ANGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, true, nIncl, nfilter);
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

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", 1 SNP allowed by penceil)...";
		pens.refOpen = 40;
		pens.readOpen = 40;
		pens.snp = 30;
		pa.exEnds = false;
		doTestCase2("1310231", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, true, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0); assert_eq(res.alres.score().score(), -30);
		assert_eq(1, res.alres.ned().size());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", 1 SNP allowed by penceil)...";
		pens.refOpen = 40;
		pens.readOpen = 40;
		pens.snp = 30;
		pa.exEnds = false;
		doTestCase2("0213131", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, true, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -30);
		assert_eq(1, res.alres.ned().size());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", 1 SNP allowed by penceil 2)...";
		pens.refOpen = 20;
		pens.readOpen = 20;
		pens.snp = 20;
		pa.exEnds = false;
		doTestCase2("1310231", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 20.0f, 0.0f, res, true, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0); assert_eq(res.alres.score().score(), -20);
		assert_eq(1, res.alres.ned().size());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", 1 MM allowed by penceil)...";
		pens.refOpen = 40;
		pens.readOpen = 40;
		pens.snp = 30;
		pa.exEnds = false;
		doTestCase2("1310131", "IIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 30.0f, 0.0f, res, true, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0); assert_eq(res.alres.score().score(), -30);
		assert_eq(0, res.alres.ned().size());
		assert_eq(1, res.alres.ced().size());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", 1 MM allowed by penceil 2)...";
		pens.snp = 20;
		pa.exEnds = false;
		doTestCase2("1310131", "5555555", "ACGTACGTACGTACGT", i*4, pa, pens2, 20.0f, 0.0f, res, true, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0); assert_eq(res.alres.score().score(), -20);
		assert_eq(0, res.alres.ned().size());
		assert_eq(1, res.alres.ced().size());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", 1 read gap)...";
		pens.refOpen = 40;
		pens.readOpen = 40;
		pens.snp = 30;
		pa.exEnds = false;
		doTestCase2("131231", "IIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 40.0f, 0.0f, res, true, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1); assert_eq(res.alres.score().score(), -40);
		assert_eq(1, res.alres.ned().size());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", 1 ref gap)...";
		pens.refOpen = 40;
		pens.readOpen = 40;
		pens.snp = 30;
		pa.exEnds = false;
		doTestCase2("13130131", "IIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 40.0f, 0.0f, res, true, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 1); assert_eq(res.alres.score().score(), -40);
		assert_eq(1, res.alres.ned().size());
		assert(res.alres.ced().empty());
		assert(res.alres.aed().empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (color space, offset " << (i*4) << ", 2 ref gaps)...";
		pens.refOpen = 40;
		pens.readOpen = 40;
		pens.refExConst = 15;
		pens.refExLinear = 0;
		pens.readOpen = 40;
		pens.readExConst = 15;
		pens.readExLinear = 0;
		pens.snp = 30;
		pa.exEnds = false;
		doTestCase2("131003131", "IIIIIIIII", "ACGTACGTACGTACGT", i*4, pa, pens, 80.0f, 0.0f, res, true, nIncl, nfilter);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 2); assert_eq(res.alres.score().score(), -55);
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
	pens.mmcostType = COST_MODEL_CONSTANT;
	pens.mmcost = 10;
	pens.snp = 30;
	pens.nCeilConst = 0.0f;
	pens.nCeilLinear = 0.0f;
	pens.readOpen = 10;
	pens.refOpen = 10;
	pens.readExConst = 10;
	pens.readExLinear = 10;
	pens.refExConst = 10;
	pens.refExLinear = 10;
	pens.setNPen(COST_MODEL_CONSTANT, 2);
	pa.exEnds = false;
	pa.gapBar = 1;
	// No Ns allowed, so this hit should be filtered
	doTestCase2(
		"ACGTACGT", // read seq
		"IIIIIIII", // read quals
		"NCGTACGT", // ref seq
		0,          // offset
		pa,         // alignment params
		pens,       // penalties
		25.0f,      //
		0.0f,
		res,
		false,  // colorspace
		false, // ns are in inclusive
		true, // nfilter
		0);
	assert(res.empty());
	cerr << "PASSED" << endl;
	res.reset();

	// No Ns allowed, so this hit should be filtered
	cerr << "  Test " << tests++ << " (N ceiling 2)...";
	pens.nCeilConst = 1.0f;
	pens.nCeilLinear = 0.0f;
	doTestCase2(
		"ACGTACGT", // read seq
		"IIIIIIII", // read quals
		"NCGTACGT", // ref seq
		0,          // offset
		pa,         // alignment params
		pens,       // penalties
		25.0f,      //
		0.0f,
		res,
		false,  // colorspace
		false, // ns are in inclusive
		false, // nfilter - NOTE: FILTER OFF
		0);
	assert(!res.empty());
	assert_eq(0,  res.alres.score().gaps());
	assert_eq(-2, res.alres.score().score());
	assert_eq(1,  res.alres.score().ns());
	cerr << "PASSED" << endl;
	res.reset();

	// No Ns allowed, so this hit should be filtered
	cerr << "  Test " << tests++ << " (N ceiling 3)...";
	pens.nCeilConst = 1.0f;
	pens.nCeilLinear = 0.0f;
	doTestCase2(
		"ACGTACGT", // read seq
		"IIIIIIII", // read quals
		"NCGTACGT", // ref seq
		0,          // offset
		pa,         // alignment params
		pens,       // penalties
		25.0f,      //
		0.0f,
		res,
		false,  // colorspace
		false,  // ns are in inclusive
		true,   // nfilter - NOTE: FILTER ON
		0);
	assert(!res.empty());
	assert_eq(0,  res.alres.score().gaps());
	assert_eq(-2, res.alres.score().score());
	assert_eq(1,  res.alres.score().ns());
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
	penMmcType      = DEFAULT_MM_PENALTY_TYPE;
	penMmc          = DEFAULT_MM_PENALTY;
	penSnp          = DEFAULT_SNP_PENALTY;
	penNType        = DEFAULT_N_PENALTY_TYPE;
	penN            = DEFAULT_N_PENALTY;
	penRdOpen       = DEFAULT_READ_OPEN;
	penRfOpen       = DEFAULT_REF_OPEN;
	penRdExConst    = DEFAULT_READ_EXTEND_CONST;
	penRfExConst    = DEFAULT_REF_EXTEND_CONST;
	penRdExLinear   = DEFAULT_READ_EXTEND_LINEAR;
	penRfExLinear   = DEFAULT_REF_EXTEND_LINEAR;
	costCeilConst   = DEFAULT_CEIL_CONST;
	costCeilLinear  = DEFAULT_CEIL_LINEAR;
	nCeilConst      = 1.0f; // constant factor in N ceiling w/r/t read length
	nCeilLinear     = 1.0f; // coeff of linear term in N ceiling w/r/t read length
	multiseedMms    = DEFAULT_SEEDMMS;
	multiseedLen    = DEFAULT_SEEDLEN;
	multiseedPeriod = DEFAULT_SEEDPERIOD;
	multiseedIvalType = DEFAULT_IVAL;
	multiseedIvalA    = DEFAULT_IVAL_A;
	multiseedIvalB    = DEFAULT_IVAL_B;
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
				bool noisyHpolymer = false;
				SeedAlignmentPolicy::parseString(
					optarg,
					noisyHpolymer,
					penMmcType,
					penMmc,
					penSnp,
					penNType,
					penN,
					penRdOpen,
					penRfOpen,
					penRdExConst,
					penRfExConst,
					penRdExLinear,
					penRfExLinear,
					costCeilConst,
					costCeilLinear,
					nCeilConst,
					nCeilLinear,
					multiseedMms,
					multiseedLen,
					multiseedPeriod,
					multiseedIvalType,
					multiseedIvalA,
					multiseedIvalB);
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
	Penalties pens(
		penMmcType,    // how to penalize mismatches
		penMmc,        // constant if mm pelanty is a constant
		penSnp,        // penalty for nucleotide mismatch in decoded colorspace als
		nCeilConst,    // N ceiling constant coefficient
		nCeilLinear,   // N ceiling linear coefficient
		penNType,      // how to penalize Ns in the read
		penN,          // constant if N pelanty is a constant
		nPairCat,      // true -> concatenate mates before N filtering
		penRdOpen,     // cost of opening a gap in the read
		penRfOpen,     // cost of opening a gap in the reference
		penRdExConst,  // constant cost of extending a gap in the read
		penRfExConst,  // constant cost of extending a gap in the reference
		penRdExLinear, // coefficient of linear term for cost of gap extension in read
		penRfExLinear  // coefficient of linear term for cost of gap extension in ref
	);
	// Calculate the penalty ceiling for the read
	int penceil = Constraint::instantiate(
		read.length(),
		costCeilConst,
		costCeilLinear);
	SwParams pa;
	SwResult res;
	doTestCase(
		read,
		qual,
		ref,
		off,
		NULL,
		NULL,
		pa,
		pens,
		penceil,
		res,
		gColor,
		nsInclusive,
		nfilter,
		seed);
}
#endif /*MAIN_ALIGNER_SW*/
