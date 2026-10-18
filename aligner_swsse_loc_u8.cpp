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

/**
 * aligner_sw_sse.cpp
 *
 * Versions of key alignment functions that use vector instructions to
 * accelerate dynamic programming.  Based chiefly on the striped Smith-Waterman
 * paper and implementation by Michael Farrar.  See:
 *
 * Farrar M. Striped Smith-Waterman speeds database searches six times over
 * other SIMD implementations. Bioinformatics. 2007 Jan 15;23(2):156-61.
 * http://sites.google.com/site/farrarmichael/smith-waterman
 *
 * While the paper describes an implementation of Smith-Waterman, we extend it
 * do end-to-end read alignment as well as local alignment.  The change
 * required for this is minor: we simply let vmax be the maximum element in the
 * score domain rather than the minimum.
 *
 * The vectorized dynamic programming implementation lacks some features that
 * make it hard to adapt to solving the entire dynamic-programming alignment
 * problem.  For instance:
 *
 * - It doesn't respect gap barriers on either end of the read
 * - It just gives a maximum; not enough information to backtrace without
 *   redoing some alignment
 * - It's a little difficult to handle st_ and en_, especially st_.
 * - The query profile mechanism makes handling of ambiguous reference bases a
 *   little tricky (16 cols in query profile lookup table instead of 5)
 *
 * Given the drawbacks, it is tempting to use SSE dynamic programming as a
 * filter rather than as an aligner per se.  Here are a few ideas for how it
 * can be extended to handle more of the alignment problem:
 *
 * - Save calculated scores to a big array as we go.  We return to this array
 *   to find and backtrace from good solutions.
 */

#include <limits>
#include "aligner_sw.h"

static const size_t NWORDS_PER_REG  = NBYTES_PER_REG;
//static const size_t NBITS_PER_WORD  = 8;
static const size_t NBYTES_PER_WORD = 1;

// In local mode, we start low (0) and go high (255).  Factoring in a query
// profile involves unsigned saturating addition.  All query profile elements
// should be expressed as a positive number; this is done by adding -min
// where min is the smallest (negative) score in the query profile.

typedef uint8_t TCScore;

/**
 * Build query profile look up tables for the read.  The query profile look
 * up table is organized as a 1D array indexed by [i][j] where i is the
 * reference character in the current DP column (0=A, 1=C, etc), and j is
 * the segment of the query we're currently working on.
 */
void SwAligner::buildQueryProfileLocalSseU8(bool fw) {
	bool& done = fw ? sseU8fwBuilt_ : sseU8rcBuilt_;
	if(done) {
		return;
	}
	done = true;
	const BTDnaString* rd = fw ? rdfw_ : rdrc_;
	const BTString* qu = fw ? qufw_ : qurc_;
	const size_t len = rd->length();
	const size_t seglen = (len + (NWORDS_PER_REG-1)) / NWORDS_PER_REG;
	// How many SSERegI's are needed
	size_t nsses =
		64 +                    // slack bytes, for alignment?
		(seglen * ALPHA_SIZE)   // query profile data
		* 2;                    // & gap barrier data
	assert_gt(nsses, 0);
	SSEData& d = fw ? sseU8fw_ : sseU8rc_;
	d.profbuf_.resizeNoCopy(nsses);
	assert(!d.profbuf_.empty());
	d.maxPen_      = d.maxBonus_ = 0;
	d.lastIter_    = d.lastWord_ = 0;
	d.qprofStride_ = d.gbarStride_ = 2;
	d.bias_ = 0;
	// The scorer's qualProf table holds sc_->score(readc, 1 << refc,
	// readq - 33) for the five reference characters at once, as bytes, with
	// the smallest score and the largest penalty and bonus among them.
	const char *rdbuf = rd->buf();
	const char *qubuf = qu->buf();
	// Calculate bias: minus the most negative score, or 0
	int minsc = 0;
	for(size_t j = 0; j < len; j++) {
		const int readc = rdbuf[j];
		const int readq = qubuf[j];
		assert_range(0, 4, readc);
		assert_range(0, 255, readq - 33);
		minsc = std::min(minsc, sc_->qualProf[readc][readq - 33].locMin);
	}
	d.bias_ = -minsc;
	assert_geq(d.bias_, 0);
	const int bias = d.bias_;
	if(len > 0) {
		// Remember which 128-bit word and which smaller word has the final
		// row: read position j is in segment j % seglen, word j / seglen
		d.lastIter_ = (len - 1) % seglen;
		d.lastWord_ = (len - 1) / seglen;
	}
	// For each reference character A, C, G, T, N there are seglen segments,
	// each a query-profile register followed by a gap-barrier register.
	// Word k of segment i holds read position j = k * seglen + i.
	SSERegI *pvec = d.profbuf_.ptr();
	uint8_t *prof = reinterpret_cast<uint8_t*>(pvec);
	const size_t segBytes = 2 * NBYTES_PER_REG;  // one segment's two registers
	const size_t refBytes = seglen * segBytes;   // one reference character's segments
	// The gap barrier depends on the read position alone: 0xff where j < len
	// and j or len - j - 1 is below gapbar, 0 elsewhere.  Zero it, then mark
	// the barrier positions at both ends.
	SSERegI vzero = sse_setzero_siall();
	for(size_t s = 0; s < seglen * ALPHA_SIZE; s++) {
		sse_store_siall(pvec + 2 * s + 1, vzero);
	}
	{
		const size_t gapbar = (size_t)sc_->gapbar;
		const size_t nbar = std::min(len, gapbar);
		// j < gapbar for j in [0, nbar); len - j - 1 < gapbar for j in
		// [len - nbar, len)
		const size_t ranges[2][2] = {
			{ 0, nbar },
			{ std::max(nbar, len - nbar), len } };
		for(size_t r = 0; r < 2; r++) {
			if(ranges[r][0] >= ranges[r][1]) {
				continue;
			}
			size_t i = ranges[r][0] % seglen, k = ranges[r][0] / seglen;
			for(size_t j = ranges[r][0]; j < ranges[r][1]; j++) {
				uint8_t *g = prof + i * segBytes + NBYTES_PER_REG + k;
				for(size_t refc = 0; refc < ALPHA_SIZE; refc++) {
					g[refc * refBytes] = 0xff;
				}
				if(++i == seglen) {
					i = 0;
					k++;
				}
			}
		}
	}
	// Query profile: (uint8_t)(score + bias), adding the bias to each of the
	// five bytes modulo 256 without carries between them; the bias alone (a
	// score of 0) past the read's end.
	const uint64_t hi = 0x8080808080808080ULL;
	const uint64_t biasw = (uint64_t)(uint8_t)bias * 0x0101010101010101ULL;
	size_t maxPen = 0, maxBonus = 0;
	for(size_t i = 0; i < seglen; i++) {
		uint8_t *p = prof + i * segBytes;
		size_t j = i, k = 0;
		for(; j < len; j += seglen, k++) {
			const int readc = rdbuf[j];
			const int readq = qubuf[j];
			const Scoring::QualProf& e = sc_->qualProf[readc][readq - 33];
			const uint64_t w = (((e.loc & ~hi) + (biasw & ~hi)) ^ ((e.loc ^ biasw) & hi));
			assert_eq((uint8_t)w,         (uint8_t)(sc_->score(readc, 1, readq - 33) + bias));
			assert_eq((uint8_t)(w >> 8),  (uint8_t)(sc_->score(readc, 2, readq - 33) + bias));
			assert_eq((uint8_t)(w >> 16), (uint8_t)(sc_->score(readc, 4, readq - 33) + bias));
			assert_eq((uint8_t)(w >> 24), (uint8_t)(sc_->score(readc, 8, readq - 33) + bias));
			assert_eq((uint8_t)(w >> 32), (uint8_t)(sc_->score(readc, 16, readq - 33) + bias));
			maxPen = std::max(maxPen, (size_t)e.locPen);
			maxBonus = std::max(maxBonus, (size_t)e.locBonus);
			p[k] = (uint8_t)w;
			p[k + refBytes] = (uint8_t)(w >> 8);
			p[k + 2 * refBytes] = (uint8_t)(w >> 16);
			p[k + 3 * refBytes] = (uint8_t)(w >> 24);
			p[k + 4 * refBytes] = (uint8_t)(w >> 32);
		}
		for(; k < NWORDS_PER_REG; k++) {
			for(size_t refc = 0; refc < ALPHA_SIZE; refc++) {
				p[k + refc * refBytes] = (uint8_t)bias;
			}
		}
	}
	d.maxPen_ = maxPen;
	d.maxBonus_ = maxBonus;
}

#ifndef NDEBUG
/**
 * Return true iff the cell has sane E/F/H values w/r/t its predecessors.
 */
static bool cellOkLocalU8(
	SSEData& d,
	size_t row,
	size_t col,
	int refc,
	int readc,
	int readq,
	const Scoring& sc)     // scoring scheme
{
	TCScore floorsc = 0;
	TCScore ceilsc = 255 - d.bias_ - 1;
	TAlScore offsetsc = 0;
	TAlScore sc_h_cur = (TAlScore)d.mat_.helt(row, col);
	TAlScore sc_e_cur = (TAlScore)d.mat_.eelt(row, col);
	TAlScore sc_f_cur = (TAlScore)d.mat_.felt(row, col);
	if(sc_h_cur > floorsc) {
		sc_h_cur += offsetsc;
	}
	if(sc_e_cur > floorsc) {
		sc_e_cur += offsetsc;
	}
	if(sc_f_cur > floorsc) {
		sc_f_cur += offsetsc;
	}
	bool gapsAllowed = true;
	size_t rowFromEnd = d.mat_.nrow() - row - 1;
	if(row < (size_t)sc.gapbar || rowFromEnd < (size_t)sc.gapbar) {
		gapsAllowed = false;
	}
	bool e_left_trans = false, h_left_trans = false;
	bool f_up_trans   = false, h_up_trans = false;
	bool h_diag_trans = false;
	if(gapsAllowed) {
		TAlScore sc_h_left = floorsc;
		TAlScore sc_e_left = floorsc;
		TAlScore sc_h_up   = floorsc;
		TAlScore sc_f_up   = floorsc;
		if(col > 0 && sc_e_cur > floorsc && sc_e_cur <= ceilsc) {
			sc_h_left = d.mat_.helt(row, col-1) + offsetsc;
			sc_e_left = d.mat_.eelt(row, col-1) + offsetsc;
			e_left_trans = (sc_e_left > floorsc && sc_e_cur == sc_e_left - sc.readGapExtend());
			h_left_trans = (sc_h_left > floorsc && sc_e_cur == sc_h_left - sc.readGapOpen());
			assert(e_left_trans || h_left_trans);
		}
		if(row > 0 && sc_f_cur > floorsc && sc_f_cur <= ceilsc) {
			sc_h_up = d.mat_.helt(row-1, col) + offsetsc;
			sc_f_up = d.mat_.felt(row-1, col) + offsetsc;
			f_up_trans = (sc_f_up > floorsc && sc_f_cur == sc_f_up - sc.refGapExtend());
			h_up_trans = (sc_h_up > floorsc && sc_f_cur == sc_h_up - sc.refGapOpen());
			assert(f_up_trans || h_up_trans);
		}
	} else {
		assert_geq(floorsc, sc_e_cur);
		assert_geq(floorsc, sc_f_cur);
	}
	if(col > 0 && row > 0 && sc_h_cur > floorsc && sc_h_cur <= ceilsc) {
		TAlScore sc_h_upleft = d.mat_.helt(row-1, col-1) + offsetsc;
		TAlScore sc_diag = sc.score(readc, (int)refc, readq - 33);
		h_diag_trans = sc_h_cur == sc_h_upleft + sc_diag;
	}
	assert(
		sc_h_cur <= floorsc ||
		e_left_trans ||
		h_left_trans ||
		f_up_trans   ||
		h_up_trans   ||
		h_diag_trans ||
		sc_h_cur > ceilsc ||
		row == 0 ||
		col == 0);
	return true;
}
#endif /*ndef NDEBUG*/

#ifdef NDEBUG

#define assert_all_eq0(x)
#define assert_all_gt(x, y)
#define assert_all_gt_lo(x)
#define assert_all_lt(x, y)
#define assert_all_lt_hi(x)

#else

#define assert_all_eq0(x) { \
	SSERegI z = sse_setzero_siall(); \
	SSERegI tmp = sse_setzero_siall(); \
	z = sse_xor_siall(z, z); \
	tmp = sse_cmpeq_epi16(x, z); \
	assert_eq(SSE_MASK_ALL, sse_movemask_epi8(tmp)); \
}

#define assert_all_gt(x, y) { \
	SSERegI tmp = sse_cmpgt_epu8(x, y); \
	assert_eq(SSE_MASK_ALL, sse_movemask_epi8(tmp)); \
}

#define assert_all_gt_lo(x) { \
	SSERegI z = sse_setzero_siall(); \
	SSERegI tmp = sse_setzero_siall(); \
	z = sse_xor_siall(z, z); \
	tmp = sse_cmpgt_epu8(x, z); \
	assert_eq(SSE_MASK_ALL, sse_movemask_epi8(tmp)); \
}

#define assert_all_lt(x, y) { \
	SSERegI z = sse_setzero_siall(); \
	z = sse_xor_siall(z, z); \
	SSERegI tmp = sse_subs_epu8(y, x); \
	tmp = sse_cmpeq_epi16(tmp, z); \
	assert_eq(0x0000, sse_movemask_epi8(tmp)); \
}

#define assert_all_lt_hi(x) { \
	SSERegI z = sse_setzero_siall(); \
	SSERegI tmp = sse_setzero_siall(); \
	z = sse_cmpeq_epu8(z, z); \
	z = sse_srli_epu8(z, 1); \
	tmp = sse_cmplt_epu8(x, z); \
	assert_eq(SSE_MASK_ALL, sse_movemask_epi8(tmp)); \
}
#endif

/**
 * Aligns by filling a dynamic programming matrix with the SSE-accelerated,
 * banded DP approach of Farrar.  As it goes, it determines which cells we
 * might backtrace from and tallies the best (highest-scoring) N backtrace
 * candidate cells per diagonal.  Also returns the alignment score of the best
 * alignment in the matrix.
 *
 * This routine does *not* maintain a matrix holding the entire matrix worth of
 * scores, nor does it maintain any other dense O(mn) data structure, as this
 * would quickly exhaust memory for queries longer than about 10,000 kb.
 * Instead, in the fill stage it maintains two columns worth of scores at a
 * time (current/previous, or right/left) - these take O(m) space.  When
 * finished with the current column, it determines which cells from the
 * previous column, if any, are candidates we might backtrace from to find a
 * full alignment.  A candidate cell has a score that rises above the threshold
 * and isn't improved upon by a match in the next column.  The best N
 * candidates per diagonal are stored in a O(m + n) data structure.
 */
TAlScore SwAligner::alignGatherLoc8(int& flag, bool debug) {
	assert_leq(rdf_, rd_->length());
	assert_leq(rdf_, qu_->length());
	assert_lt(rfi_, rff_);
	assert_lt(rdi_, rdf_);
	assert_eq(rd_->length(), qu_->length());
	assert_geq(sc_->gapbar, 1);
	assert_gt(minsc_, 0);
	assert(repOk());
#ifndef NDEBUG
	for(size_t i = (size_t)rfi_; i < (size_t)rff_; i++) {
		assert_range(0, 16, (int)rf_[i]);
	}
#endif

	SSEData& d = fw_ ? sseU8fw_ : sseU8rc_;
	SSEMetrics& met = extend_ ? sseU8ExtendMet_ : sseU8MateMet_;
	if(!debug) met.dp++;
	buildQueryProfileLocalSseU8(fw_);
	assert(!d.profbuf_.empty());
	assert_gt(d.bias_, 0);
	assert_lt(d.bias_, 127);

	assert_gt(d.maxBonus_, 0);
	size_t iter =
		(dpRows() + (NWORDS_PER_REG-1)) / NWORDS_PER_REG; // iter = segLen

	// Now set up the score vectors.  We just need two columns worth, which
	// we'll call "left" and "right".
	d.vecbuf_.resize(ROWSTRIDE_2COL * iter * 2);
	d.vecbuf_.zero();
	SSERegI *vbuf_l = d.vecbuf_.ptr();
	SSERegI *vbuf_r = d.vecbuf_.ptr() + (ROWSTRIDE_2COL * iter);

	// This is the data structure that holds candidate cells per diagonal.
	const size_t ndiags = rff_ - rfi_ + dpRows() - 1;
	if(!debug) {
		btdiag_.init(ndiags, 2);
	}

	// Data structure that holds checkpointed anti-diagonals
	TAlScore perfectScore = sc_->perfectScore(dpRows());
	bool checkpoint = true;
	bool cpdebug = false;
#ifndef NDEBUG
	cpdebug = dpRows() < 1000;
#endif
	cper_.init(
		dpRows(),      // # rows
		rff_ - rfi_,   // # columns
		cperPerPow2_,  // checkpoint every 1 << perpow2 diags (& next)
		perfectScore,  // perfect score (for sanity checks)
		true,          // matrix cells have 8-bit scores?
		cperTri_,      // triangular mini-fills?
		true,          // alignment is local?
		cpdebug);      // save all cells for debugging?

	// Many thanks to Michael Farrar for releasing his striped Smith-Waterman
	// implementation:
	//
	//  http://sites.google.com/site/farrarmichael/smith-waterman
	//
	// Much of the implmentation below is adapted from Michael's code.

	// Set all elts to reference gap open penalty
	SSERegI rfgapo   = sse_setzero_siall();
	SSERegI rfgape   = sse_setzero_siall();
	SSERegI rdgapo   = sse_setzero_siall();
	SSERegI rdgape   = sse_setzero_siall();
	SSERegI vlo      = sse_setzero_siall();
	SSERegI vhi      = sse_setzero_siall();
	SSERegI vmax     = sse_setzero_siall();
	SSERegI vcolmax  = sse_setzero_siall();
	SSERegI ve       = sse_setzero_siall();
	SSERegI vf       = sse_setzero_siall();
	SSERegI vh       = sse_setzero_siall();
	SSERegI vhd      = sse_setzero_siall();
	SSERegI vhdtmp   = sse_setzero_siall();
	SSERegI vtmp     = sse_setzero_siall();
	SSERegI vzero    = sse_setzero_siall();
	SSERegI vbias    = sse_setzero_siall();
	SSERegI vbiasm1  = sse_setzero_siall();

	assert_gt(sc_->refGapOpen(), 0);
	assert_leq(sc_->refGapOpen(), MAX_U8);
	sse_fill_u8(sc_->refGapOpen(), rfgapo);

	// Set all elts to reference gap extension penalty
	assert_gt(sc_->refGapExtend(), 0);
	assert_leq(sc_->refGapExtend(), MAX_U8);
	assert_leq(sc_->refGapExtend(), sc_->refGapOpen());
	sse_fill_u8(sc_->refGapExtend(), rfgape);

	// Set all elts to read gap open penalty
	assert_gt(sc_->readGapOpen(), 0);
	assert_leq(sc_->readGapOpen(), MAX_U8);
	sse_fill_u8(sc_->readGapOpen(), rdgapo);

	// Set all elts to read gap extension penalty
	assert_gt(sc_->readGapExtend(), 0);
	assert_leq(sc_->readGapExtend(), MAX_U8);
	assert_leq(sc_->readGapExtend(), sc_->readGapOpen());
	sse_fill_u8(sc_->readGapExtend(), rdgape);

	// Set all elts to minimum score threshold.  Actually, to 1 less than the
	// threshold so we can use gt instead of geq.
	sse_fill_u8((d.bias_ - 1), vbiasm1);

	sse_fill_u8_opt(0xff, vhi);
	sse_fill_u8_opt(0, vlo);
	vmax = vlo;

	// Make a vector of bias offsets
	sse_fill_u8(d.bias_, vbias);

	// Points to a long vector of SSERegI where each element is a block of
	// contiguous cells in the E, F or H matrix.  If the index % 3 == 0, then
	// the block of cells is from the E matrix.  If index % 3 == 1, they're
	// from the F matrix.  If index % 3 == 2, then they're from the H matrix.
	// Blocks of cells are organized in the same interleaved manner as they are
	// calculated by the Farrar algorithm.
	const SSERegI *pvScore; // points into the query profile

	const size_t colstride = ROWSTRIDE_2COL * iter;

	// Initialize the H and E vectors in the first matrix column
	SSERegI *pvELeft = vbuf_l + 0; SSERegI *pvERight = vbuf_r + 0;
	/* SSERegI *pvFLeft = vbuf_l + 1; */ SSERegI *pvFRight = vbuf_r + 1;
	SSERegI *pvHLeft = vbuf_l + 2; SSERegI *pvHRight = vbuf_r + 2;

	for(size_t i = 0; i < iter; i++) {
		// start low in local mode
		sse_store_siall(pvERight, vlo); pvERight += ROWSTRIDE_2COL;
		sse_store_siall(pvHRight, vlo); pvHRight += ROWSTRIDE_2COL;
	}

	assert_gt(sc_->gapbar, 0);
	size_t nfixup = 0;
	TAlScore matchsc = sc_->match(30);
	TAlScore leftmax = MIN_I64;

	// Fill in the table as usual but instead of using the same gap-penalty
	// vector for each iteration of the inner loop, load words out of a
	// pre-calculated gap vector parallel to the query profile.  The pre-
	// calculated gap vectors enforce the gap barrier constraint by making it
	// infinitely costly to introduce a gap in barrier rows.
	//
	// AND use a separate loop to fill in the first row of the table, enforcing
	// the st_ constraints in the process.  This is awkward because it
	// separates the processing of the first row from the others and might make
	// it difficult to use the first-row results in the next row, but it might
	// be the simplest and least disruptive way to deal with the st_ constraint.

	size_t off = MAX_SIZE_T, lastoff;
	bool bailed = false;
	for(size_t i = (size_t)rfi_; i < (size_t)rff_; i++) {
		// Swap left and right; vbuf_l is the vector on the left, which we
		// generally load from, and vbuf_r is the vector on the right, which we
		// generally store to.
		swap(vbuf_l, vbuf_r);
		pvELeft = vbuf_l + 0; pvERight = vbuf_r + 0;
		/* pvFLeft = vbuf_l + 1; */ pvFRight = vbuf_r + 1;
		pvHLeft = vbuf_l + 2; pvHRight = vbuf_r + 2;

		// Fetch this column's reference mask
		const int refm = (int)rf_[i];

		// Fetch the appropriate query profile
		lastoff = off;
		off = (size_t)firsts5[refm] * iter * 2;
		pvScore = d.profbuf_.ptr() + off; // even elts = query profile, odd = gap barrier

		// Load H vector from the final row of the previous column.
		// ??? perhaps we should calculate the next iter's F instead of the
		// current iter's?  The way we currently do it, seems like it will
		// almost always require at least one fixup loop iter (to recalculate
		// this topmost F).
		vh = sse_load_siall(pvHLeft + colstride - ROWSTRIDE_2COL);

		// Set all cells to low value
		vf = sse_xor_siall(vf, vf);
		// vf now contains the vertical contribution

		// Store cells in F, calculated previously
		// No need to veto ref gap extensions, they're all 0x00s
		sse_store_siall(pvFRight, vf);
		pvFRight += ROWSTRIDE_2COL;

		// Shift down so that topmost (least sig) cell gets 0
		vh = sse_slli_siall(vh, NBYTES_PER_WORD);

		// We pull out one loop iteration to make it easier to veto values in the top row

		// Load cells from E, calculated previously
		ve = sse_load_siall(pvELeft);
		vhd = sse_load_siall(pvHLeft);
		assert_all_lt(ve, vhi);
		pvELeft += ROWSTRIDE_2COL;
		// ve now contains the horizontal contribution

		// Factor in query profile (matches and mismatches)
		vh = sse_adds_epu8(vh, pvScore[0]);
		vh = sse_subs_epu8(vh, vbias);
		// vh now contains the diagonal contribution

		vhdtmp = vhd;
		vhd = sse_subs_epu8(vhd, rdgapo);
		vhd = sse_subs_epu8(vhd, pvScore[1]); // veto some read gap opens
		ve = sse_subs_epu8(ve, rdgape);
		ve = sse_max_epu8(ve, vhd);

		vh = sse_max_epu8(vh, ve);
		vf = vh;

		// Update highest score so far
		vcolmax = vh;

		// Save the new vH values
		sse_store_siall(pvHRight, vh);

		vh = vhdtmp;
		assert_all_lt(ve, vhi);
		pvHRight += ROWSTRIDE_2COL;
		pvHLeft += ROWSTRIDE_2COL;

		// Save E values
		sse_store_siall(pvERight, ve);
		pvERight += ROWSTRIDE_2COL;

		// Update vf value
		vf = sse_subs_epu8(vf, rfgapo);
		assert_all_lt(vf, vhi);

		pvScore += 2; // move on to next query profile

		// For each character in the reference text:
		size_t j;
		for(j = 1; j < iter; j++) {
			// Load cells from E, calculated previously
			ve = sse_load_siall(pvELeft);
			vhd = sse_load_siall(pvHLeft);
			assert_all_lt(ve, vhi);
			pvELeft += ROWSTRIDE_2COL;

			// Store cells in F, calculated previously
			vf = sse_subs_epu8(vf, pvScore[1]); // veto some ref gap extensions
			sse_store_siall(pvFRight, vf);
			pvFRight += ROWSTRIDE_2COL;

			// Factor in query profile (matches and mismatches)
			vh = sse_adds_epu8(vh, pvScore[0]);
			vh = sse_subs_epu8(vh, vbias);

			// Update H, factoring in E and F
			vh = sse_max_epu8(vh, vf);

			vhdtmp = vhd;
			vhd = sse_subs_epu8(vhd, rdgapo);
			vhd = sse_subs_epu8(vhd, pvScore[1]); // veto some read gap opens
			ve = sse_subs_epu8(ve, rdgape);
			ve = sse_max_epu8(ve, vhd);

			vh = sse_max_epu8(vh, ve);
			vtmp = vh;

			// Update highest score encountered this far
			vcolmax = sse_max_epu8(vcolmax, vh);

			// Save the new vH values
			sse_store_siall(pvHRight, vh);

			vh = vhdtmp;

			assert_all_lt(ve, vhi);
			pvHRight += ROWSTRIDE_2COL;
			pvHLeft += ROWSTRIDE_2COL;

			// Save E values
			sse_store_siall(pvERight, ve);
			pvERight += ROWSTRIDE_2COL;

			// Update vf value
			vtmp = sse_subs_epu8(vtmp, rfgapo);
			vf = sse_subs_epu8(vf, rfgape);
			assert_all_lt(vf, vhi);
			vf = sse_max_epu8(vf, vtmp);

			pvScore += 2; // move on to next query profile / gap veto
		}
		// pvHStore, pvELoad, pvEStore have all rolled over to the next column
		pvFRight -= colstride; // reset to start of column
		vtmp = sse_load_siall(pvFRight);

		pvHRight -= colstride; // reset to start of column
		vh = sse_load_siall(pvHRight);

		pvScore = d.profbuf_.ptr() + off + 1; // reset veto vector

		// vf from last row gets shifted down by one to overlay the first row
		// rfgape has already been subtracted from it.
		vf = sse_slli_siall(vf, NBYTES_PER_WORD);

		vf = sse_subs_epu8(vf, *pvScore); // veto some ref gap extensions
		vf = sse_max_epu8(vtmp, vf);
		// TODO: We're testing whether F changed.  Can't we just assume that F
		// did change and instead check whether H changed?  Might save us from
		// entering the fixup loop.
		vtmp = sse_subs_epu8(vf, vtmp);
		vtmp = sse_cmpeq_epi8(vtmp, vzero);
		int cmp = sse_movemask_epi8(vtmp);

		// If any element of vtmp is greater than H - gap-open...
		j = 0;
		while(cmp != SSE_MASK_ALL) {
			// Store this vf
			sse_store_siall(pvFRight, vf);
			pvFRight += ROWSTRIDE_2COL;

			// Update vh w/r/t new vf
			vh = sse_max_epu8(vh, vf);

			// Save vH values
			sse_store_siall(pvHRight, vh);
			pvHRight += ROWSTRIDE_2COL;

			// Update highest score encountered so far.
			vcolmax = sse_max_epu8(vcolmax, vh);

			pvScore += 2;

			assert_lt(j, iter);
			if(++j == iter) {
				pvFRight -= colstride;
				vtmp = sse_load_siall(pvFRight);   // load next vf ASAP
				pvHRight -= colstride;
				vh = sse_load_siall(pvHRight);     // load next vh ASAP
				pvScore = d.profbuf_.ptr() + off + 1;
				j = 0;
				vf = sse_slli_siall(vf, NBYTES_PER_WORD);
			} else {
				vtmp = sse_load_siall(pvFRight);   // load next vf ASAP
				vh = sse_load_siall(pvHRight);     // load next vh ASAP
			}

			// Update F with another gap extension
			vf = sse_subs_epu8(vf, rfgape);
			vf = sse_subs_epu8(vf, *pvScore); // veto some ref gap extensions
			vf = sse_max_epu8(vtmp, vf);
			vtmp = sse_subs_epu8(vf, vtmp);
			vtmp = sse_cmpeq_epi8(vtmp, vzero);
			cmp = sse_movemask_epi8(vtmp);
			nfixup++;
		}

		// Now we'd like to know exactly which cells in the left column are
		// candidates we might backtrace from.  First question is: did *any*
		// elements in the column exceed the minimum score threshold?
		if(!debug && leftmax >= minsc_) {
			// Yes.  Next question is: which cells are candidates?  We have to
			// allow matches in the right column to override matches above and
			// to the left in the left column.
			assert_gt(i - rfi_, 0);
			pvHLeft  = vbuf_l + 2;
			assert_lt(lastoff, MAX_SIZE_T);
			pvScore = d.profbuf_.ptr() + lastoff; // even elts = query profile, odd = gap barrier
			for(size_t k = 0; k < iter; k++) {
				vh = sse_load_siall(pvHLeft);
				vtmp = sse_cmpgt_epi8(pvScore[0], vbiasm1);
				int cmp = sse_movemask_epi8(vtmp);
				if(cmp != 0) {
					// At least one candidate in this mask.  Now iterate
					// through vm/vh to evaluate individual cells.
					for(size_t m = 0; m < NWORDS_PER_REG; m++) {
						size_t row = k + m * iter;
						if(row >= dpRows()) {
							break;
						}
						if(((TCScore *)&vtmp)[m] > 0 && ((TCScore *)&vh)[m] >= minsc_) {
							TCScore sc = ((TCScore *)&vh)[m];
							assert_geq(sc, minsc_);
							// Add to data structure holding all candidates
							size_t col = i - rfi_ - 1; // -1 b/c prev col
							size_t frombot = dpRows() - row - 1;
							DpBtCandidate cand(row, col, sc);
							btdiag_.add(frombot + col, cand);
						}
					}
				}
				pvHLeft += ROWSTRIDE_2COL;
				pvScore += 2;
			}
		}

		// Save some elements to checkpoints
		if(checkpoint) {

			SSERegI *pvE = vbuf_r + 0;
			SSERegI *pvF = vbuf_r + 1;
			SSERegI *pvH = vbuf_r + 2;
			size_t coli = i - rfi_;
			if(coli < cper_.locol_) cper_.locol_ = coli;
			if(coli > cper_.hicol_) cper_.hicol_ = coli;
			if(cperTri_) {
				// Checkpoint for triangular mini-fills
				size_t rc_mod = coli & cper_.lomask_;
				assert_lt(rc_mod, cper_.per_);
				int64_t row = -rc_mod-1;
				int64_t row_mod = row;
				int64_t row_div = 0;
				size_t idx = coli >> cper_.perpow2_;
				size_t idxrow = idx * cper_.nrow_;
				assert_eq(4, ROWSTRIDE_2COL);
				bool done = false;
				while(true) {
					row += (cper_.per_ - 2);
					row_mod += (cper_.per_ - 2);
					for(size_t j = 0; j < 2; j++) {
						row++;
						row_mod++;
						if(row >= 0 && (size_t)row < cper_.nrow_) {
							// Update row divided by iter_ and mod iter_
							while(row_mod >= (int64_t)iter) {
								row_mod -= (int64_t)iter;
								row_div++;
							}
							size_t delt = idxrow + row;
							size_t vecoff = (row_mod << (ROWSTRIDE_LOG2+BYTES_LOG2_PER_REG)) + row_div;
							assert_lt(row_div, NBYTES_PER_REG);
							int16_t h_sc = ((uint8_t*)pvH)[vecoff];
							int16_t e_sc = ((uint8_t*)pvE)[vecoff];
							int16_t f_sc = ((uint8_t*)pvF)[vecoff];
							assert_leq(h_sc, cper_.perf_);
							assert_leq(e_sc, cper_.perf_);
							assert_leq(f_sc, cper_.perf_);
							CpQuad *qdiags = ((j == 0) ? cper_.qdiag1s_.ptr() : cper_.qdiag2s_.ptr());
							qdiags[delt].sc[0] = h_sc;
							qdiags[delt].sc[1] = e_sc;
							qdiags[delt].sc[2] = f_sc;
						} // if(row >= 0 && row < nrow_)
						else if(row >= 0 && (size_t)row >= cper_.nrow_) {
							done = true;
							break;
						}
					} // for(size_t j = 0; j < 2; j++)
					if(done) {
						break;
					}
					idx++;
					idxrow += cper_.nrow_;
				} // while(true)
			} else {
				// Checkpoint for square mini-fills

				// If this is the first column, take this opportunity to
				// pre-calculate the coordinates of the elements we're going to
				// checkpoint.
				if(coli == 0) {
					size_t cpi    = cper_.per_-1;
					size_t cpimod = cper_.per_-1;
					size_t cpidiv = 0;
					cper_.commitMap_.clear();
					while(cpi < cper_.nrow_) {
						while(cpimod >= iter) {
							cpimod -= iter;
							cpidiv++;
						}
						size_t vecoff = (cpimod << (ROWSTRIDE_LOG2+BYTES_LOG2_PER_REG)) + cpidiv;
						cper_.commitMap_.push_back(vecoff);
						cpi += cper_.per_;
						cpimod += cper_.per_;
					}
				}
				// Save all the rows
				size_t rowoff = 0;
				size_t sz = cper_.commitMap_.size();
				for(size_t i = 0; i < sz; i++, rowoff += cper_.ncol_) {
					size_t vecoff = cper_.commitMap_[i];
					int16_t h_sc = ((uint8_t*)pvH)[vecoff];
					//int16_t e_sc = ((uint8_t*)pvE)[vecoff];
					int16_t f_sc = ((uint8_t*)pvF)[vecoff];
					assert_leq(h_sc, cper_.perf_);
					//assert_leq(e_sc, cper_.perf_);
					assert_leq(f_sc, cper_.perf_);
					CpQuad& dst = cper_.qrows_[rowoff + coli];
					dst.sc[0] = h_sc;
					//dst.sc[1] = e_sc;
					dst.sc[2] = f_sc;
				}
				// Is this a column we'd like to checkpoint?
				if((coli & cper_.lomask_) == cper_.lomask_) {
					// Save the column using memcpys
					assert_gt(coli, 0);
					size_t wordspercol = cper_.niter_ * ROWSTRIDE_2COL;
					size_t coloff = (coli >> cper_.perpow2_) * wordspercol;
					SSERegI *dst = cper_.qcols_.ptr() + coloff;
					memcpy(dst, vbuf_r, sizeof(SSERegI) * wordspercol);
				}
			}
			if(cper_.debug_) {
				// Save the column using memcpys
				size_t wordspercol = cper_.niter_ * ROWSTRIDE_2COL;
				size_t coloff = coli * wordspercol;
				SSERegI *dst = cper_.qcolsD_.ptr() + coloff;
				memcpy(dst, vbuf_r, sizeof(SSERegI) * wordspercol);
			}
		}

		// Store column maximum vector in first element of tmp
		vmax = sse_max_epu8(vmax, vcolmax);

		{
			// Get single largest score in this column
			int score = 0;
			sse_max_score_u8(vcolmax, score);

			// Could we have saturated?
			if(score + d.bias_ >= 255) {
				flag = -2; // yes
				if(!debug) met.dpsat++;
				return MIN_I64;
			}

			if(score < minsc_) {
				size_t ncolleft = rff_ - i - 1;
				if(score + (TAlScore)ncolleft * matchsc < minsc_) {
					// Bail!  There can't possibly be a valid alignment that
					// passes through this column.
					bailed = true;
					break;
				}
			}

			leftmax = score;
		}
	}

	lastoff = off;

	// Now we'd like to know exactly which cells in the *rightmost* column are
	// candidates we might backtrace from.  Did *any* elements exceed the
	// minimum score threshold?
	if(!debug && !bailed && leftmax >= minsc_) {
		// Yes.  Next question is: which cells are candidates?  We have to
		// allow matches in the right column to override matches above and
		// to the left in the left column.
		pvHLeft  = vbuf_r + 2;
		assert_lt(lastoff, MAX_SIZE_T);
		pvScore = d.profbuf_.ptr() + lastoff; // even elts = query profile, odd = gap barrier
		for(size_t k = 0; k < iter; k++) {
			vh = sse_load_siall(pvHLeft);
			vtmp = sse_cmpgt_epi8(pvScore[0], vbiasm1);
			int cmp = sse_movemask_epi8(vtmp);
			if(cmp != 0) {
				// At least one candidate in this mask.  Now iterate
				// through vm/vh to evaluate individual cells.
				for(size_t m = 0; m < NWORDS_PER_REG; m++) {
					size_t row = k + m * iter;
					if(row >= dpRows()) {
						break;
					}
					if(((TCScore *)&vtmp)[m] > 0 && ((TCScore *)&vh)[m] >= minsc_) {
						TCScore sc = ((TCScore *)&vh)[m];
						assert_geq(sc, minsc_);
						// Add to data structure holding all candidates
						size_t col = rff_ - rfi_ - 1; // -1 b/c prev col
						size_t frombot = dpRows() - row - 1;
						DpBtCandidate cand(row, col, sc);
						btdiag_.add(frombot + col, cand);
					}
				}
			}
			pvHLeft += ROWSTRIDE_2COL;
			pvScore += 2;
		}
	}

	// Update metrics
	if(!debug) {
		size_t ninner = (rff_ - rfi_) * iter;
		met.col   += (rff_ - rfi_);             // DP columns
		met.cell  += (ninner * NWORDS_PER_REG); // DP cells
		met.inner += ninner;                    // DP inner loop iters
		met.fixup += nfixup;                    // DP fixup loop iters
	}

	// Find largest score in vmax
	int score = 0;
	sse_max_score_u8(vmax, score);
	flag = 0;

	// Could we have saturated?
	if(score + d.bias_ >= 255) {
		flag = -2; // yes
		if(!debug) met.dpsat++;
		return MIN_I64;
	}

	// Did we find a solution?
	if(score == MIN_U8 || score < minsc_) {
		flag = -1; // no
		if(!debug) met.dpfail++;
		return (TAlScore)score;
	}

	// Now take all the backtrace candidates in the btdaig_ structure and
	// dump them into the btncand_ array.  They'll be sorted later.
	if(!debug) {
		assert(!btdiag_.empty());
		btdiag_.dump(btncand_);
		assert(!btncand_.empty());
	}

	// Return largest score
	if(!debug) met.dpsucc++;
	return (TAlScore)score;
}

#if NBYTES_PER_REG == 32
/**
 * The column loop of alignNucleotidesLocalSseU8 for a column of ITER = 2 to 4
 * 256-bit vectors (33 to 128 rows).  It computes the general loop's saturating
 * u8 recurrence (see the comments there) and stores the same E, F, H and
 * column-maximum vectors in the same places, and stops at the same column.
 * The general loop's speed is set by the chain of dependent instructions from
 * one column to the next: H from the last vector row through the diagonal
 * shift, the first pass's F through every row, and Farrar's lazy-F loop, which
 * carries F one lane further per pass and whose exit is a data-dependent
 * branch.  Here the chain is a prefix scan with in-half shifts and a few
 * instructions, and every value passed between columns stays in a register.
 * This is alignEE8V256Short (aligner_swsse_ee_u8.cpp) with local mode's
 * diagonal: H0[j] = max(diag + profile - bias, E[j]), saturating at both ends,
 * and 0 shifted into row 0.
 *
 * Row j of lane L is read row L * ITER + j.  For each column:
 *  - within a lane, F flows from row j to row j + 1 at the cost of rfgapo
 *    from H0 or of the barrier and rfgape from F (rfgapo >= rfgape); HL and
 *    FL are H and F with these carries only;
 *  - the carry c into each lane from the lanes above is the prefix scan of
 *    the lanes' carry-outs, done within each 128-bit half and then across
 *    the halves once;
 *  - each row takes c less the cost from the top of its lane (vpen): F =
 *    max(FL, that), H = max(HL, that), and the next column's E = max(E -
 *    rdgape, H - rdgapo - barrier).
 * The next column's carry-out, the scan's input, is max_j(H0'[j] - cpen[j]).
 * Every operation is monotone, so it is also the max of terms in HL and E,
 * which are ready while the scan runs, and of terms in c and in c shifted
 * down a lane, each taking its least cost (zc, wc, kb).  The chain from one
 * column to the next is then the scan and one subtraction.
 *
 * sat(sat(x - a) - b) = sat(x - sat(a + b)), so costs may be summed before
 * they are subtracted.  The diagonal adds the profile, a gain where the read
 * matches: a term in c through it is sat(sat(sat(c - vpen) + P) - bias -
 * cpen).  While sat(c - vpen) + P < 255 that is sat(c - (vpen + cpen + bias -
 * P)), and the caller makes vpen + cpen + bias - P >= 0 (the largest bonus is
 * at most rfgapo) and computes cpen + bias - P without saturating (rfgapo +
 * (ITER-1) * rfgape + bias <= 255).  Where sat(c - vpen) + P >= 255, the next
 * column's H0 reaches 255 - bias, so that column saturates in both loops,
 * which return at it with the same flag.  Its own cells may then differ; the
 * caller does not use the matrix of a saturated fill.
 *
 * Returns -2 if a column saturated, else 0 with the maximum of H over the
 * filled columns in vmax.
 */
template<size_t ITER>
static int alignLoc8V256Short(
	const __m256i *prof,    // query profile: per reference character, ITER (score, barrier) pairs
	__m256i *mat,           // matrix: column 0, row 0, E (SSEMatrix layout)
	size_t colstride,       // vectors per matrix column
	const char *rf,         // reference
	size_t rfi,             // first column
	size_t rff,             // last column (excl.)
	const Scoring& sc,      // gap penalties
	int bias,               // added to every profile score
	TAlScore minsc,         // minimum score of a valid alignment
	TAlScore matchsc,       // score of a match, for the bail-out test
	size_t& colstop,        // out: columns filled (as colstop_)
	size_t& lastsolcol,     // out: last column reaching minsc (as lastsolcol_)
	__m256i& vmax)          // out: maximum of H over the filled columns
{
	const size_t E = SSEMatrix::E, F = SSEMatrix::F, H = SSEMatrix::H, TMP = SSEMatrix::TMP;
	const __m256i zero = _mm256_setzero_si256();
	const __m256i vrfgapo = _mm256_set1_epi8((char)sc.refGapOpen());
	const __m256i vrfgape = _mm256_set1_epi8((char)sc.refGapExtend());
	const __m256i vrdgapo = _mm256_set1_epi8((char)sc.readGapOpen());
	const __m256i vrdgape = _mm256_set1_epi8((char)sc.readGapExtend());
	const __m256i vbias = _mm256_set1_epi8((char)bias);
	// Gap barrier rows (the same for every reference character), and the costs
	// derived from them, per vector row j:
	//  bpe:   F carried from row j to row j + 1 (barrier, then rfgape);
	//  cpen:  H0 at row j to the lane's carry-out (rfgapo, then bpe below);
	//  vpen:  the carry into the lane to row j (as the general loop's vpen);
	//  epen:  H to the next column's E (rdgapo, then barrier);
	//  below: 0xff where a row below row j in the lane is in the barrier.
	__m256i bar[ITER], bpe[ITER], cpen[ITER], vpen[ITER], epen[ITER], below[ITER];
	__m256i vbar = zero;
	for(size_t j = 0; j < ITER; j++) {
		bar[j] = _mm256_load_si256(prof + 2 * j + 1);
		bpe[j] = _mm256_adds_epu8(bar[j], vrfgape);
		epen[j] = _mm256_adds_epu8(vrdgapo, bar[j]);
		vpen[j] = j == 0 ? bar[0] : _mm256_adds_epu8(_mm256_adds_epu8(vpen[j-1], vrfgape), bar[j]);
		vbar = _mm256_or_si256(vbar, bar[j]);
	}
	cpen[ITER-1] = vrfgapo;
	below[ITER-1] = zero;
	for(size_t j = ITER-1; j > 0; j--) {
		cpen[j-1] = _mm256_adds_epu8(cpen[j], bpe[j]);
		below[j-1] = _mm256_or_si256(below[j], bar[j]);
	}
	// The scan's costs: passing lane L whole is pass(L); qN(L) is the cost of
	// passing the N lanes that end at L, within a half.  chi(L) is the cost of
	// passing lanes 16 to L - 1, for the carry from lane 15 into lane L > 16.
	size_t passPen = ITER * (size_t)sc.refGapExtend();
	if(passPen > 0xff) {
		passPen = 0xff;
	}
	const __m256i q1 = _mm256_or_si256(_mm256_set1_epi8((char)passPen), vbar);
	const __m256i q2 = _mm256_adds_epu8(q1, _mm256_slli_si256(q1, 1));
	const __m256i q3 = _mm256_adds_epu8(q2, _mm256_slli_si256(q1, 2));
	const __m256i q4 = _mm256_adds_epu8(q2, _mm256_slli_si256(q2, 2));
	const __m256i q8 = _mm256_adds_epu8(q4, _mm256_slli_si256(q4, 4));
	const __m256i q12 = _mm256_adds_epu8(q8, _mm256_slli_si256(q4, 8));
	const __m256i q16 = _mm256_adds_epu8(q8, _mm256_slli_si256(q8, 8));
	const __m256i chi = _mm256_slli_si256(q16, 1);
	// Selects byte 15 of the upper half (lane 15, once moved there) into
	// every byte of the upper half, and zero into the lower half
	const __m256i vidx15 = _mm256_set_epi64x(
		0x0f0f0f0f0f0f0f0fLL, 0x0f0f0f0f0f0f0f0fLL,
		(long long)0x8080808080808080ULL, (long long)0x8080808080808080ULL);
	// Costs to the next column's carry-out, max_j(H0'[j] - cpen[j]), where
	// H0'[j] = max(diag[j] + profile'[j] - bias, E'[j]), E'[j] = max(E[j] -
	// rdgape, H[j] - epen[j]), H[j] = max(HL[j], c - vpen[j]) and diag[j] =
	// H[j - 1] or, for j = 0, H[ITER - 1] shifted down a lane (0 into lane 0):
	//  rc[j]: from E[j] (rdgape, cpen); ec[j]: from H[j] (epen, cpen);
	//  bc[j]: from diag[j] after the profile is added (bias, cpen);
	//  zc[r]: from c, the least over its paths, for reference character r;
	//  wc[r]: from c shifted down a lane; kb[r]: from lane 15 of the scan,
	//  carried into the upper half, for both (c and c shifted are taken from
	//  the scan's values directly).
	// Through the diagonal, a path's cost is vpen + cpen + bias - profile.
	// cpen + bias - profile is computed from the scalar cpen of a lane without
	// a barrier below (the caller keeps it at most 255), as 0xff where there
	// is one.
	__m256i rc[ITER], ec[ITER], bc[ITER], zc[5], wc[5], kb[5];
	// chi for c shifted down a lane: lane 16 takes lane 14's value, with no
	// carry from lane 15
	const __m256i chish = _mm256_or_si256(_mm256_slli_si256(chi, 1),
		_mm256_set_epi64x(0, 0xff, 0, 0));
	const __m256i vpenSh = _mm256_alignr_epi8(vpen[ITER-1],
		_mm256_permute2x128_si256(vpen[ITER-1], vpen[ITER-1], 0x08), 15);
	__m256i cbias[ITER];
	for(size_t j = 0; j < ITER; j++) {
		rc[j] = _mm256_adds_epu8(vrdgape, cpen[j]);
		ec[j] = _mm256_adds_epu8(epen[j], cpen[j]);
		bc[j] = _mm256_adds_epu8(vbias, cpen[j]);
		cbias[j] = _mm256_set1_epi8((char)(sc.refGapOpen() +
			(int)(ITER - 1 - j) * sc.refGapExtend() + bias));
	}
	for(size_t r = 0; r < 5; r++) {
		const __m256i *pr = prof + r * ITER * 2;
		__m256i z = _mm256_adds_epu8(vpen[0], ec[0]);
		for(size_t j = 1; j < ITER; j++) {
			const __m256i cps = _mm256_or_si256(
				_mm256_subs_epu8(cbias[j], _mm256_load_si256(pr + 2 * j)), below[j]);
			z = _mm256_min_epu8(z, _mm256_adds_epu8(vpen[j-1], cps));
			z = _mm256_min_epu8(z, _mm256_adds_epu8(vpen[j], ec[j]));
		}
		zc[r] = z;
		const __m256i cps0 = _mm256_or_si256(
			_mm256_subs_epu8(cbias[0], _mm256_load_si256(pr)), below[0]);
		wc[r] = _mm256_adds_epu8(vpenSh, cps0);
		kb[r] = _mm256_min_epu8(_mm256_adds_epu8(chi, zc[r]), _mm256_adds_epu8(chish, wc[r]));
	}
	// lpen[j]: H0 at row j to F of the lane's last row (its barrier included)
	__m256i lpen[ITER-1];
	lpen[ITER-2] = vrfgapo;
	for(size_t j = ITER-2; j > 0; j--) {
		lpen[j-1] = _mm256_adds_epu8(lpen[j], bpe[j]);
	}
	for(size_t j = 0; j + 1 < ITER; j++) {
		lpen[j] = _mm256_adds_epu8(lpen[j], bar[ITER-1]);
	}
	// H0 and E of the column, and its carry-out: for the first, from the
	// staged initial column (zero), so every diagonal is 0
	__m256i h0[ITER], ecur[ITER], co = zero;
	{
		const __m256i *pvS = prof + (size_t)firsts5[(int)rf[rfi]] * ITER * 2;
		for(size_t j = 0; j < ITER; j++) {
			h0[j] = _mm256_subs_epu8(_mm256_load_si256(pvS + 2 * j), vbias);
			ecur[j] = zero;
			co = _mm256_max_epu8(co, _mm256_subs_epu8(h0[j], cpen[j]));
		}
	}
	__m256i vm = zero;
	size_t ls = 0;
	colstop = rff - rfi;
	__m256i *pvCol = mat;
	for(size_t i = rfi; i < rff; i++) {
		// The next column's reference character.  (After the last column,
		// the values computed with it are not used.)
		const size_t rn = (size_t)firsts5[(int)rf[i + 1 < rff ? i + 1 : i]];
		const __m256i *pvS = prof + rn * ITER * 2;
		// The scan within each half, four lanes per step: over the 4 lanes
		// ending at each lane, then over 16
		__m256i y = co;
		y = _mm256_max_epu8(
			_mm256_max_epu8(y, _mm256_subs_epu8(_mm256_slli_si256(y, 1), q1)),
			_mm256_max_epu8(_mm256_subs_epu8(_mm256_slli_si256(y, 2), q2),
			                _mm256_subs_epu8(_mm256_slli_si256(y, 3), q3)));
		y = _mm256_max_epu8(
			_mm256_max_epu8(y, _mm256_subs_epu8(_mm256_slli_si256(y, 4), q4)),
			_mm256_max_epu8(_mm256_subs_epu8(_mm256_slli_si256(y, 8), q8),
			                _mm256_subs_epu8(_mm256_slli_si256(y, 12), q12)));
		// The carry into each lane is the scan's value one lane up (y15), and
		// into the upper half also lane 15's (yc), carried down
		const __m256i py = _mm256_permute2x128_si256(y, y, 0x08);
		const __m256i y15 = _mm256_alignr_epi8(y, py, 15);
		const __m256i y14 = _mm256_alignr_epi8(y, py, 14);
		const __m256i yc = _mm256_shuffle_epi8(py, vidx15);
		const __m256i c = _mm256_max_epu8(y15, _mm256_subs_epu8(yc, chi));
		// The next column's carry-out: its terms in c and in c shifted down a
		// lane (from y14 and yc)
		const __m256i cc = _mm256_max_epu8(
			_mm256_subs_epu8(y14, wc[rn]), _mm256_subs_epu8(yc, kb[rn]));
		// c shifted down a lane (0 into lane 0), for the next column's
		// diagonal
		const __m256i csh = _mm256_max_epu8(y14, _mm256_subs_epu8(yc, chish));
		// F and H with the carries from within the lane only
		__m256i fl[ITER], hl[ITER];
		fl[0] = zero;
		hl[0] = h0[0];
		{
			__m256i vf = _mm256_subs_epu8(h0[0], vrfgapo);
			for(size_t j = 1; j < ITER; j++) {
				// (into the last row of ITER > 2, vf has the barrier, from lpen)
				fl[j] = j + 1 == ITER && ITER > 2 ? vf : _mm256_subs_epu8(vf, bar[j]);
				hl[j] = _mm256_max_epu8(h0[j], fl[j]);
				if(j + 2 < ITER) {
					vf = _mm256_max_epu8(
						_mm256_subs_epu8(vf, bpe[j]), _mm256_subs_epu8(h0[j], vrfgapo));
				} else if(j + 2 == ITER) {
					// into the last row, as a maximum over the rows above
					__m256i fz = _mm256_subs_epu8(h0[0], lpen[0]);
					for(size_t k = 1; k < ITER-1; k++) {
						fz = _mm256_max_epu8(fz, _mm256_subs_epu8(h0[k], lpen[k]));
					}
					vf = fz;
				}
			}
		}
		// ... and its other terms, the latest (from the last row of HL) last.
		// hlsh: HL's last row shifted down a lane (0 into lane 0)
		__m256i hlsh;
		{
			__m256i a = _mm256_max_epu8(
				_mm256_subs_epu8(ecur[0], rc[0]), _mm256_subs_epu8(ecur[1], rc[1]));
			for(size_t j = 2; j < ITER; j++) {
				a = _mm256_max_epu8(a, _mm256_subs_epu8(ecur[j], rc[j]));
			}
			for(size_t j = 0; j + 1 < ITER; j++) {
				a = _mm256_max_epu8(a, _mm256_max_epu8(
					_mm256_subs_epu8(hl[j], ec[j]),
					_mm256_subs_epu8(_mm256_adds_epu8(hl[j], _mm256_load_si256(pvS + 2 * (j+1))), bc[j+1])));
			}
			hlsh = _mm256_alignr_epi8(hl[ITER-1],
				_mm256_permute2x128_si256(hl[ITER-1], hl[ITER-1], 0x08), 15);
			a = _mm256_max_epu8(a, _mm256_max_epu8(
				_mm256_subs_epu8(hl[ITER-1], ec[ITER-1]),
				_mm256_subs_epu8(_mm256_adds_epu8(hlsh, _mm256_load_si256(pvS)), bc[0])));
			co = _mm256_max_epu8(_mm256_max_epu8(a, _mm256_subs_epu8(y15, zc[rn])), cc);
		}
		// F, H and the next column's E with the carry from the lanes above
		__m256i hcur[ITER];
		__m256i colmax = zero;
		__m256i *pv = pvCol;
		for(size_t j = 0; j < ITER; j++) {
			const __m256i s = _mm256_subs_epu8(c, vpen[j]);
			const __m256i f = j == 0 ? s : _mm256_max_epu8(fl[j], s);
			const __m256i h = _mm256_max_epu8(hl[j], s);
			const __m256i e = _mm256_max_epu8(
				_mm256_subs_epu8(ecur[j], vrdgape), _mm256_subs_epu8(h, epen[j]));
			_mm256_store_si256(pv + F, f);
			_mm256_store_si256(pv + H, h);
			_mm256_store_si256(pv + colstride + E, e);
			pv += ROWSTRIDE;
			colmax = _mm256_max_epu8(colmax, h);
			hcur[j] = h;
			ecur[j] = e;
		}
		// The column's maximum, its saturation and bail-out tests, as the
		// general loop's
		_mm256_store_si256(pvCol + TMP, colmax);
		vm = _mm256_max_epu8(vm, colmax);
		{
			// The largest byte of colmax (as sse_max_score_u8): 255 less the
			// smallest byte of its complement, taken by phminposuw over the
			// minima of byte pairs
			int score;
			{
				const __m128i m = _mm_max_epu8(_mm256_castsi256_si128(colmax),
					_mm256_extracti128_si256(colmax, 1));
				const __m128i n = _mm_xor_si128(m, _mm_set1_epi8((char)0xff));
				const __m128i p = _mm_min_epu8(n, _mm_srli_epi16(n, 8));
				score = 255 - (_mm_cvtsi128_si32(_mm_minpos_epu16(p)) & 0xffff);
			}
			if(score + bias >= 255) {
				lastsolcol = ls;
				vmax = vm;
				return -2;
			}
			if(score < minsc) {
				size_t ncolleft = rff - i - 1;
				if(score + (TAlScore)ncolleft * matchsc < minsc) {
					colstop = (i+1) - rfi;
					break;
				}
			} else {
				ls = i - rfi;
			}
		}
		pvCol += colstride;
		// H0 of the next column.  Its diagonal is row j - 1 of this column,
		// or for j = 0 the last vector row shifted down a lane, with 0 in
		// lane 0.  That is max(HL, c - vpen) shifted, taken from hlsh and csh
		// so that the shift is not after c.
		{
			const __m256i d0 = _mm256_max_epu8(hlsh, _mm256_subs_epu8(csh, vpenSh));
			h0[0] = _mm256_max_epu8(_mm256_subs_epu8(
				_mm256_adds_epu8(d0, _mm256_load_si256(pvS)), vbias), ecur[0]);
		}
		for(size_t j = 1; j < ITER; j++) {
			h0[j] = _mm256_max_epu8(_mm256_subs_epu8(
				_mm256_adds_epu8(hcur[j-1], _mm256_load_si256(pvS + 2 * j)), vbias), ecur[j]);
		}
	}
	lastsolcol = ls;
	vmax = vm;
	return 0;
}
#endif

/**
 * Solve the current alignment problem using SSE instructions that operate on 16
 * unsigned 8-bit values packed into a single 128-bit register.
 */
TAlScore SwAligner::alignNucleotidesLocalSseU8(int& flag, bool debug) {
	assert_leq(rdf_, rd_->length());
	assert_leq(rdf_, qu_->length());
	assert_lt(rfi_, rff_);
	assert_lt(rdi_, rdf_);
	assert_eq(rd_->length(), qu_->length());
	assert_geq(sc_->gapbar, 1);
	assert(repOk());
#ifndef NDEBUG
	for(size_t i = (size_t)rfi_; i < (size_t)rff_; i++) {
		assert_range(0, 16, (int)rf_[i]);
	}
#endif

	SSEData& d = fw_ ? sseU8fw_ : sseU8rc_;
	SSEMetrics& met = extend_ ? sseU8ExtendMet_ : sseU8MateMet_;
	if(!debug) met.dp++;
	buildQueryProfileLocalSseU8(fw_);
	assert(!d.profbuf_.empty());
	assert_geq(d.bias_, 0);

	assert_gt(d.maxBonus_, 0);
	size_t iter =
		(dpRows() + (NWORDS_PER_REG-1)) / NWORDS_PER_REG; // iter = segLen

	// Many thanks to Michael Farrar for releasing his striped Smith-Waterman
	// implementation:
	//
	//  http://sites.google.com/site/farrarmichael/smith-waterman
	//
	// Much of the implmentation below is adapted from Michael's code.

	// Set all elts to reference gap open penalty
	SSERegI rfgapo   = sse_setzero_siall();
	SSERegI rfgape   = sse_setzero_siall();
	SSERegI rdgapo   = sse_setzero_siall();
	SSERegI rdgape   = sse_setzero_siall();
	SSERegI vlo      = sse_setzero_siall();
	SSERegI vhi      = sse_setzero_siall();
	SSERegI vmax     = sse_setzero_siall();
	SSERegI vcolmax  = sse_setzero_siall();
	SSERegI ve       = sse_setzero_siall();
	SSERegI vf       = sse_setzero_siall();
	SSERegI vh       = sse_setzero_siall();
	SSERegI vtmp     = sse_setzero_siall();
	SSERegI vzero    = sse_setzero_siall();
	SSERegI vbias    = sse_setzero_siall();

	assert_gt(sc_->refGapOpen(), 0);
	assert_leq(sc_->refGapOpen(), MAX_U8);
	sse_fill_u8(sc_->refGapOpen(), rfgapo);

	// Set all elts to reference gap extension penalty
	assert_gt(sc_->refGapExtend(), 0);
	assert_leq(sc_->refGapExtend(), MAX_U8);
	assert_leq(sc_->refGapExtend(), sc_->refGapOpen());
	sse_fill_u8(sc_->refGapExtend(), rfgape);

	// Set all elts to read gap open penalty
	assert_gt(sc_->readGapOpen(), 0);
	assert_leq(sc_->readGapOpen(), MAX_U8);
	sse_fill_u8(sc_->readGapOpen(), rdgapo);

	// Set all elts to read gap extension penalty
	assert_gt(sc_->readGapExtend(), 0);
	assert_leq(sc_->readGapExtend(), MAX_U8);
	assert_leq(sc_->readGapExtend(), sc_->readGapOpen());
	sse_fill_u8(sc_->readGapExtend(), rdgape);

	sse_fill_u8_opt(0xff, vhi);
	sse_fill_u8_opt(0, vlo);
	vmax = vlo;

	// Make a vector of bias offsets
	sse_fill_u8(d.bias_, vbias);

	// Points to a long vector of SSERegI where each element is a block of
	// contiguous cells in the E, F or H matrix.  If the index % 3 == 0, then
	// the block of cells is from the E matrix.  If index % 3 == 1, they're
	// from the F matrix.  If index % 3 == 2, then they're from the H matrix.
	// Blocks of cells are organized in the same interleaved manner as they are
	// calculated by the Farrar algorithm.
	const SSERegI *pvScore; // points into the query profile

	d.mat_.init(dpRows(), rff_ - rfi_, NWORDS_PER_REG);
	const size_t colstride = d.mat_.colstride();
	//const size_t rowstride = d.mat_.rowstride();
	assert_eq(ROWSTRIDE, colstride / iter);

	// Initialize the H and E vectors in the first matrix column
	SSERegI *pvHTmp = d.mat_.tmpvec(0, 0);
	SSERegI *pvETmp = d.mat_.evec(0, 0);

	for(size_t i = 0; i < iter; i++) {
		sse_store_siall(pvETmp, vlo);
		sse_store_siall(pvHTmp, vlo); // start low in local mode
		pvETmp += ROWSTRIDE;
		pvHTmp += ROWSTRIDE;
	}
	// These are swapped just before the innermost loop
	SSERegI *pvHStore = d.mat_.hvec(0, 0);
	SSERegI *pvHLoad  = d.mat_.tmpvec(0, 0);
	SSERegI *pvELoad  = d.mat_.evec(0, 0);
	SSERegI *pvEStore = d.mat_.evecUnsafe(0, 1);
	SSERegI *pvFStore = d.mat_.fvec(0, 0);
	SSERegI *pvFTmp   = NULL;

	assert_gt(sc_->gapbar, 0);
	size_t nfixup = 0;
	TAlScore matchsc = sc_->match(30);

	// Fill in the table as usual but instead of using the same gap-penalty
	// vector for each iteration of the inner loop, load words out of a
	// pre-calculated gap vector parallel to the query profile.  The pre-
	// calculated gap vectors enforce the gap barrier constraint by making it
	// infinitely costly to introduce a gap in barrier rows.
	//
	// AND use a separate loop to fill in the first row of the table, enforcing
	// the st_ constraints in the process.  This is awkward because it
	// separates the processing of the first row from the others and might make
	// it difficult to use the first-row results in the next row, but it might
	// be the simplest and least disruptive way to deal with the st_ constraint.

	colstop_ = rff_ - rfi_;
	lastsolcol_ = 0;
#if NBYTES_PER_REG == 32
	// 33 to 128 rows.  The loop below spends most of each column on the chain
	// of dependent instructions from one column to the next and on the lazy-F
	// loop; alignLoc8V256Short computes the same cells with a shorter chain.
	// It needs every path cost through the diagonal to be non-negative (the
	// largest match bonus at most rfgapo) and the costs it derives from rfgapo
	// and rfgape not to saturate.
	if(iter >= 2 && iter <= 4 &&
	   d.maxBonus_ <= (size_t)sc_->refGapOpen() &&
	   sc_->refGapExtend() <= sc_->refGapOpen() &&
	   sc_->refGapOpen() + (int)(iter - 1) * sc_->refGapExtend() + d.bias_ <= 255)
	{
		size_t colstop = 0, lastsolcol = 0;
		int ret;
		if(iter == 2) {
			ret = alignLoc8V256Short<2>(d.profbuf_.ptr(), d.mat_.ptr(), colstride,
				rf_, (size_t)rfi_, (size_t)rff_, *sc_, d.bias_, minsc_, matchsc,
				colstop, lastsolcol, vmax);
		} else if(iter == 3) {
			ret = alignLoc8V256Short<3>(d.profbuf_.ptr(), d.mat_.ptr(), colstride,
				rf_, (size_t)rfi_, (size_t)rff_, *sc_, d.bias_, minsc_, matchsc,
				colstop, lastsolcol, vmax);
		} else {
			ret = alignLoc8V256Short<4>(d.profbuf_.ptr(), d.mat_.ptr(), colstride,
				rf_, (size_t)rfi_, (size_t)rff_, *sc_, d.bias_, minsc_, matchsc,
				colstop, lastsolcol, vmax);
		}
		lastsolcol_ = lastsolcol;
		if(ret == -2) {
			flag = -2; // saturated
			if(!debug) met.dpsat++;
			return MIN_I64;
		}
		colstop_ = colstop;
		// The scan does the work of one fix-up pass per column
		nfixup += colstop_ * iter;
#ifndef NDEBUG
		for(size_t i = (size_t)rfi_; i < (size_t)rfi_ + colstop_; i++) {
			if((rand() & 15) == 0) {
				for(size_t k = 0; k < dpRows(); k++) {
					assert(cellOkLocalU8(
						d,
						k,                   // row
						i - rfi_,            // col
						(int)rf_[i],         // reference mask
						(int)(*rd_)[rdi_+k], // read char
						(int)(*qu_)[rdi_+k], // read quality
						*sc_));              // scoring scheme
				}
			}
		}
#endif
	} else
#endif
	for(size_t i = (size_t)rfi_; i < (size_t)rff_; i++) {
		assert(pvFStore == d.mat_.fvec(0, i - rfi_));
		assert(pvHStore == d.mat_.hvec(0, i - rfi_));

		// Fetch this column's reference mask
		const int refm = (int)rf_[i];

		// Fetch the appropriate query profile
		size_t off = (size_t)firsts5[refm] * iter * 2;
		pvScore = d.profbuf_.ptr() + off; // even elts = query profile, odd = gap barrier

		// Load H vector from the final row of the previous column
		vh = sse_load_siall(pvHLoad + colstride - ROWSTRIDE);

		// Set all cells to low value
		vf = sse_xor_siall(vf, vf);

		// Store cells in F, calculated previously
		// No need to veto ref gap extensions, they're all 0x00s
		sse_store_siall(pvFStore, vf);
		pvFStore += ROWSTRIDE;

		// Shift down so that topmost (least sig) cell gets 0
		vh = sse_slli_siall(vh, NBYTES_PER_WORD);

		// We pull out one loop iteration to make it easier to veto values in the top row

		// Load cells from E, calculated previously
		ve = sse_load_siall(pvELoad);
		assert_all_lt(ve, vhi);
		pvELoad += ROWSTRIDE;

		// Factor in query profile (matches and mismatches)
		vh = sse_adds_epu8(vh, pvScore[0]);
		vh = sse_subs_epu8(vh, vbias);

		// Update H, factoring in E and F
		vh = sse_max_epu8(vh, ve);
		vh = sse_max_epu8(vh, vf);

		// Update highest score so far
		vcolmax = sse_xor_siall(vcolmax, vcolmax);
		vcolmax = sse_max_epu8(vcolmax, vh);

		// Save the new vH values
		sse_store_siall(pvHStore, vh);
		pvHStore += ROWSTRIDE;

		// Update vE value
		vf = vh;
		vh = sse_subs_epu8(vh, rdgapo);
		vh = sse_subs_epu8(vh, pvScore[1]); // veto some read gap opens
		ve = sse_subs_epu8(ve, rdgape);
		ve = sse_max_epu8(ve, vh);
		assert_all_lt(ve, vhi);

		// Load the next h value
		vh = sse_load_siall(pvHLoad);
		pvHLoad += ROWSTRIDE;

		// Save E values
		sse_store_siall(pvEStore, ve);
		pvEStore += ROWSTRIDE;

		// Update vf value
		vf = sse_subs_epu8(vf, rfgapo);
		assert_all_lt(vf, vhi);

		pvScore += 2; // move on to next query profile

		// For each character in the reference text:
		size_t j;
		for(j = 1; j < iter; j++) {
			// Load cells from E, calculated previously
			ve = sse_load_siall(pvELoad);
			assert_all_lt(ve, vhi);
			pvELoad += ROWSTRIDE;

			// Store cells in F, calculated previously
			vf = sse_subs_epu8(vf, pvScore[1]); // veto some ref gap extensions
			sse_store_siall(pvFStore, vf);
			pvFStore += ROWSTRIDE;

			// Factor in query profile (matches and mismatches)
			vh = sse_adds_epu8(vh, pvScore[0]);
			vh = sse_subs_epu8(vh, vbias);

			// Update H, factoring in E and F
			vh = sse_max_epu8(vh, ve);
			vh = sse_max_epu8(vh, vf);

			// Update highest score encountered this far
			vcolmax = sse_max_epu8(vcolmax, vh);

			// Save the new vH values
			sse_store_siall(pvHStore, vh);
			pvHStore += ROWSTRIDE;

			// Update vE value
			vtmp = vh;
			vh = sse_subs_epu8(vh, rdgapo);
			vh = sse_subs_epu8(vh, pvScore[1]); // veto some read gap opens
			ve = sse_subs_epu8(ve, rdgape);
			ve = sse_max_epu8(ve, vh);
			assert_all_lt(ve, vhi);

			// Load the next h value
			vh = sse_load_siall(pvHLoad);
			pvHLoad += ROWSTRIDE;

			// Save E values
			sse_store_siall(pvEStore, ve);
			pvEStore += ROWSTRIDE;

			// Update vf value
			vtmp = sse_subs_epu8(vtmp, rfgapo);
			vf = sse_subs_epu8(vf, rfgape);
			assert_all_lt(vf, vhi);
			vf = sse_max_epu8(vf, vtmp);

			pvScore += 2; // move on to next query profile / gap veto
		}
		// pvHStore, pvELoad, pvEStore have all rolled over to the next column
		pvFTmp = pvFStore;
		pvFStore -= colstride; // reset to start of column
		vtmp = sse_load_siall(pvFStore);

		pvHStore -= colstride; // reset to start of column
		vh = sse_load_siall(pvHStore);

		pvEStore -= colstride; // reset to start of column
		ve = sse_load_siall(pvEStore);

		pvHLoad = pvHStore;    // new pvHLoad = pvHStore
		pvScore = d.profbuf_.ptr() + off + 1; // reset veto vector

		// vf from last row gets shifted down by one to overlay the first row
		// rfgape has already been subtracted from it.
		vf = sse_slli_siall(vf, NBYTES_PER_WORD);

		vf = sse_subs_epu8(vf, *pvScore); // veto some ref gap extensions
		vf = sse_max_epu8(vtmp, vf);
		vtmp = sse_subs_epu8(vf, vtmp);
		vtmp = sse_cmpeq_epi8(vtmp, vzero);
		int cmp = sse_movemask_epi8(vtmp);

		// If any element of vtmp is greater than H - gap-open...
		j = 0;
		while(cmp != SSE_MASK_ALL) {
			// Store this vf
			sse_store_siall(pvFStore, vf);
			pvFStore += ROWSTRIDE;

			// Update vh w/r/t new vf
			vh = sse_max_epu8(vh, vf);

			// Save vH values
			sse_store_siall(pvHStore, vh);
			pvHStore += ROWSTRIDE;

			// Update highest score encountered this far
			vcolmax = sse_max_epu8(vcolmax, vh);

			// Update E in case it can be improved using our new vh
			vh = sse_subs_epu8(vh, rdgapo);
			vh = sse_subs_epu8(vh, *pvScore); // veto some read gap opens
			ve = sse_max_epu8(ve, vh);
			sse_store_siall(pvEStore, ve);
			pvEStore += ROWSTRIDE;
			pvScore += 2;

			assert_lt(j, iter);
			if(++j == iter) {
				pvFStore -= colstride;
				vtmp = sse_load_siall(pvFStore);   // load next vf ASAP
				pvHStore -= colstride;
				vh = sse_load_siall(pvHStore);     // load next vh ASAP
				pvEStore -= colstride;
				ve = sse_load_siall(pvEStore);     // load next ve ASAP
				pvScore = d.profbuf_.ptr() + off + 1;
				j = 0;
				vf = sse_slli_siall(vf, NBYTES_PER_WORD);
			} else {
				vtmp = sse_load_siall(pvFStore);   // load next vf ASAP
				vh = sse_load_siall(pvHStore);     // load next vh ASAP
				ve = sse_load_siall(pvEStore);     // load next vh ASAP
			}

			// Update F with another gap extension
			vf = sse_subs_epu8(vf, rfgape);
			vf = sse_subs_epu8(vf, *pvScore); // veto some ref gap extensions
			vf = sse_max_epu8(vtmp, vf);
			vtmp = sse_subs_epu8(vf, vtmp);
			vtmp = sse_cmpeq_epi8(vtmp, vzero);
			cmp = sse_movemask_epi8(vtmp);
			nfixup++;
		}

#ifndef NDEBUG
		if((rand() & 15) == 0) {
			// This is a work-intensive sanity check; each time we finish filling
			// a column, we check that each H, E, and F is sensible.
			for(size_t k = 0; k < dpRows(); k++) {
				assert(cellOkLocalU8(
					d,
					k,                   // row
					i - rfi_,            // col
					refm,                // reference mask
					(int)(*rd_)[rdi_+k], // read char
					(int)(*qu_)[rdi_+k], // read quality
					*sc_));              // scoring scheme
			}
		}
#endif

		// Store column maximum vector in first element of tmp
		vmax = sse_max_epu8(vmax, vcolmax);

		sse_store_siall(d.mat_.tmpvec(0, i - rfi_), vcolmax);

		{
			// Get single largest score in this column
			int score = 0;
			sse_max_score_u8(vcolmax, score);

			// Could we have saturated?
			if(score + d.bias_ >= 255) {
				flag = -2; // yes
				if(!debug) met.dpsat++;
				return MIN_I64;
			}

			if(score < minsc_) {
				size_t ncolleft = rff_ - i - 1;
				if(score + (TAlScore)ncolleft * matchsc < minsc_) {
					// Bail!  We're guaranteed not to see a valid alignment in
					// the rest of the matrix
					colstop_ = (i+1) - rfi_;
					break;
				}
			} else {
				lastsolcol_ = i - rfi_;
			}
		}

		// pvELoad and pvHLoad are already where they need to be

		// Adjust the load and store vectors here.
		pvHStore = pvHLoad + colstride;
		pvEStore = pvELoad + colstride;
		pvFStore = pvFTmp;
	}

	// Update metrics
	if(!debug) {
		size_t ninner = (rff_ - rfi_) * iter;
		met.col   += (rff_ - rfi_);             // DP columns
		met.cell  += (ninner * NWORDS_PER_REG); // DP cells
		met.inner += ninner;                    // DP inner loop iters
		met.fixup += nfixup;                    // DP fixup loop iters
	}

	// Find largest score in vmax
	int score = 0;
	sse_max_score_u8(vmax, score);
	flag = 0;

	// Could we have saturated?
	if(score + d.bias_ >= 255) {
		flag = -2; // yes
		if(!debug) met.dpsat++;
		return MIN_I64;
	}

	// Did we find a solution?
	if(score == MIN_U8 || score < minsc_) {
		flag = -1; // no
		if(!debug) met.dpfail++;
		return (TAlScore)score;
	}

	// Return largest score
	if(!debug) met.dpsucc++;
	return (TAlScore)score;
}

/**
 * Given a filled-in DP table, populate the btncand_ list with candidate cells
 * that might be at the ends of valid alignments.  No need to do this unless
 * the maximum score returned by the align*() func is >= the minimum.
 *
 * We needn't consider cells that have no chance of reaching any of the core
 * diagonals.  These are the cells that are more than 'maxgaps' cells away from
 * a core diagonal.
 *
 * We need to be careful to consider that the rectangle might be truncated on
 * one or both ends.
 *
 * The seed extend case looks like this:
 *
 *      |Rectangle|   0: seed diagonal
 *      **OO0oo----   o: "RHS gap" diagonals
 *      -**OO0oo---   O: "LHS gap" diagonals
 *      --**OO0oo--   *: "LHS extra" diagonals
 *      ---**OO0oo-   -: cells that can't possibly be involved in a valid
 *      ----**OO0oo      alignment that overlaps one of the core diagonals
 *
 * The anchor-to-left case looks like this:
 *
 *   |Anchor|  | ---- Rectangle ---- |
 *   o---------OO0000000000000oo------  0: mate diagonal (also core diags!)
 *   -o---------OO0000000000000oo-----  o: "RHS gap" diagonals
 *   --o---------OO0000000000000oo----  O: "LHS gap" diagonals
 *   ---oo--------OO0000000000000oo---  *: "LHS extra" diagonals
 *   -----o--------OO0000000000000oo--  -: cells that can't possibly be
 *   ------o--------OO0000000000000oo-     involved in a valid alignment that
 *   -------o--------OO0000000000000oo     overlaps one of the core diagonals
 *                     XXXXXXXXXXXXX
 *                     | RHS Range |
 *                     ^           ^
 *                     rl          rr
 *
 * The anchor-to-right case looks like this:
 *
 *    ll          lr
 *    v           v
 *    | LHS Range |
 *    XXXXXXXXXXXXX          |Anchor|
 *  OO0000000000000oo--------o--------  0: mate diagonal (also core diags!)
 *  -OO0000000000000oo--------o-------  o: "RHS gap" diagonals
 *  --OO0000000000000oo--------o------  O: "LHS gap" diagonals
 *  ---OO0000000000000oo--------oo----  *: "LHS extra" diagonals
 *  ----OO0000000000000oo---------o---  -: cells that can't possibly be
 *  -----OO0000000000000oo---------o--     involved in a valid alignment that
 *  ------OO0000000000000oo---------o-     overlaps one of the core diagonals
 *  | ---- Rectangle ---- |
 */
bool SwAligner::gatherCellsNucleotidesLocalSseU8(TAlScore best) {
	// What's the minimum number of rows that can possibly be spanned by an
	// alignment that meets the minimum score requirement?
	assert(sse8succ_);
	size_t bonus = (size_t)sc_->match(30);
	const size_t ncol = lastsolcol_ + 1;
	const size_t nrow = dpRows();
	assert_gt(nrow, 0);
	btncand_.clear();
	btncanddone_.clear();
	SSEData& d = fw_ ? sseU8fw_ : sseU8rc_;
	SSEMetrics& met = extend_ ? sseU8ExtendMet_ : sseU8MateMet_;
	assert(!d.profbuf_.empty());
	//const size_t rowstride = d.mat_.rowstride();
	//const size_t colstride = d.mat_.colstride();
	size_t iter = (dpRows() + (NWORDS_PER_REG - 1)) / NWORDS_PER_REG;
	assert_gt(iter, 0);
	assert_geq(minsc_, 0);
	assert_gt(bonus, 0);
	size_t minrow = (size_t)(((minsc_ + bonus - 1) / bonus) - 1);
	// A candidate is a cell in a row from minrow on whose H is at least
	// minsc_, whose read character matches the reference character, and
	// whose successor on the diagonal (down and to the right) does not match.
	// No H byte can reach a minsc_ above 255, and no row can be a candidate
	// when minrow is past the last one.  Otherwise the tests are done on a
	// whole vector of H at a time, and only the lanes that pass all of them
	// are visited, in the order of a scan by column, then vector row, then
	// lane.
	if(minsc_ <= 0xff && minrow < nrow) {
		SSERegI vminsc;
		SSERegI vzero = sse_setzero_siall();
		sse_fill_u8((int)minsc_, vminsc);
		// For each vector row i, two vectors with one byte per lane k, i.e.
		// per row i + k*iter: the row's read character as a mask (1 << c),
		// and the next row's.  Rows before minrow, and padding rows past the
		// read, get 0 in both, which matches no reference character.  The
		// last row has no successor and gets 0 in the second.
		d.vecbuf_.resizeNoCopy(2 * iter);
		uint8_t *rdmask = reinterpret_cast<uint8_t*>(d.vecbuf_.ptr());
		for(size_t i = 0; i < iter; i++) {
			uint8_t *rdc = rdmask + (2 * i) * NWORDS_PER_REG; // this row
			uint8_t *rdn = rdc + NWORDS_PER_REG;              // next row
			for(size_t k = 0; k < NWORDS_PER_REG; k++) {
				const size_t rdoff = i + k * iter;
				rdc[k] = rdn[k] = 0;
				if(rdoff >= minrow && rdoff < nrow) {
					assert_lt((int)(*rd_)[rdoff], 8);
					rdc[k] = (uint8_t)(1 << (*rd_)[rdoff]);
					if(rdoff < dpRows()-1) {
						assert_lt((int)(*rd_)[rdoff+1], 8);
						rdn[k] = (uint8_t)(1 << (*rd_)[rdoff+1]);
					}
				}
			}
		}
		for(size_t j = 0; j < ncol; j++) {
			// First, check if there is a cell in this column with a score
			// above the score threshold: the column maxima are in the first
			// TMP vector
			SSERegI vcolmax = sse_load_siall(d.mat_.tmpvec(0, j));
			int colpass = sse_movemask_epi8(
				sse_cmpeq_epi8(sse_max_epu8(vcolmax, vminsc), vcolmax));
#ifndef NDEBUG
			{
				int score = 0;
				sse_max_score_u8(vcolmax, score);
				// Start in upper vector row and move down
				TAlScore max = 0;
				SSERegI *pvH = d.mat_.hvec(0, j);
				for(size_t i = 0; i < iter; i++) {
					for(size_t k = 0; k < NWORDS_PER_REG; k++) {
						TAlScore sc = (TAlScore)((TCScore*)pvH)[k];
						if(sc > max) {
							max = sc;
						}
					}
					pvH += ROWSTRIDE;
				}
				assert_eq(max, score);
				assert_eq(colpass != 0, (TAlScore)score >= minsc_);
			}
#endif
			if(colpass == 0) {
				// Scores in column aren't good enough
				continue;
			}
			// This column's reference character and the next one's (the
			// reference has one character past the rectangle), as masks
			assert_range(0, 16, (int)rf_[j + rfi_]);
			assert_range(0, 16, (int)rf_[j + rfi_ + 1]);
			SSERegI vref, vrefSucc;
			sse_fill_u8((uint8_t)rf_[j + rfi_], vref);
			sse_fill_u8((uint8_t)rf_[j + rfi_ + 1], vrefSucc);
			// Get pointer to first cell in column to examine:
			SSERegI *pvH = d.mat_.hvec(0, j);
			const SSERegI *pvRd = d.vecbuf_.ptr();
			// Start in upper vector row and move down
			for(size_t i = 0; i < iter; i++) {
				SSERegI vh = sse_load_siall(pvH);
				// Lanes whose H is at least minsc_ ...
				uint32_t lanes = (uint32_t)sse_movemask_epi8(
					sse_cmpeq_epi8(sse_max_epu8(vh, vminsc), vh));
				// ... whose read character matches the reference character ...
				lanes &= ~(uint32_t)sse_movemask_epi8(
					sse_cmpeq_epi8(sse_and_siall(vref, pvRd[0]), vzero));
				// ... and whose diagonal successor does not match
				lanes &= (uint32_t)sse_movemask_epi8(
					sse_cmpeq_epi8(sse_and_siall(vrefSucc, pvRd[1]), vzero));
				while(lanes != 0) {
					const size_t k = (size_t)__builtin_ctz(lanes);
					lanes &= lanes - 1;
					const size_t rdoff = i + k * iter;
					assert_geq(rdoff, minrow);
					assert_lt(rdoff, nrow);
					TAlScore sc = (TAlScore)((TCScore*)pvH)[k];
					assert_leq(sc, best);
					assert_geq(sc, minsc_);
					// Yes, this is legit
					met.gathsol++;
					btncand_.expand();
					btncand_.back().init(rdoff, j, sc);
				}
				pvH += ROWSTRIDE;
				pvRd += 2;
			}
		}
	}
#ifndef NDEBUG
	{
		// The same candidates, in the same order, by a scalar scan of every
		// cell
		size_t n = 0;
		for(size_t j = 0; j < ncol; j++) {
			SSERegI *pvH = d.mat_.hvec(0, j);
			for(size_t i = 0; i < iter; i++) {
				for(size_t k = 0; k < NWORDS_PER_REG; k++) {
					size_t rdoff = i + k * iter;
					if(rdoff >= nrow) {
						break;
					}
					TAlScore sc = (TAlScore)((TCScore*)pvH)[k];
					if(rdoff < minrow || sc < minsc_) {
						continue;
					}
					bool match = ((rf_[j + rfi_] & (1 << (*rd_)[rdoff])) != 0);
					bool matchSucc = rdoff < dpRows()-1 &&
						((rf_[j + rfi_ + 1] & (1 << (*rd_)[rdoff+1])) != 0);
					if(match && !matchSucc) {
						assert_lt(n, btncand_.size());
						assert_eq(rdoff, btncand_[n].row);
						assert_eq(j, btncand_[n].col);
						assert_eq(sc, btncand_[n].score);
						n++;
					}
				}
				pvH += ROWSTRIDE;
			}
		}
		assert_eq(n, btncand_.size());
	}
#endif
	if(!btncand_.empty()) {
		d.mat_.initMasks();
	}
	return !btncand_.empty();
}

/*
 * The backtrace addresses the striped matrix by one byte pointer per cell:
 * btcur points at the current cell's E byte, element row / nvecrow_ of the
 * E vector in vector row btrowvec = row % nvecrow_ of its column.  The F and
 * H bytes are one and two vectors further on (btF, btH), the cell to the
 * left is one column (btLeft bytes) back, and the cell above is one vector
 * row back (btUp bytes), or, from vector row 0, the previous element of the
 * last vector row (btUpWrap bytes on).  These are the moves the parent's
 * MOVE_VEC_PTR_UP and MOVE_VEC_PTR_LEFT made on each of its four pointers.
 */
#define BT_SET_ROW_COL(row, col) { \
	btrowvec = (row) % btnvecrow; \
	btcur = btmatbuf + ((col) * btcolstride + btrowvec * ROWSTRIDE) * NBYTES_PER_REG + (row) / btnvecrow; \
}

/**
 * Given the dynamic programming table and a cell, trace backwards from the
 * cell and install the edits and score/penalty in the appropriate fields
 * of SwResult res, which contains an AlnRes.  The RandomSource is used to
 * break ties among equally good ways of tracing back.
 *
 * Upon entering a cell, we check if the read/ref coordinates of the cell
 * correspond to a cell we traversed constructing a previous alignment.  If so,
 * we backtrack to the last decision point, mask out the path that led to the
 * previously observed cell, and continue along a different path; or, if there
 * are no more paths to try, we give up.
 *
 * An alignment found is subject to a filtering step designed to remove
 * alignments that could spuriously trump a better alignment falling partially
 * outside the rectangle.
 *
 *          1
 *      67890123456   0: seed diagonal
 *      **OO0oo----   o: right-hand "gap" diagonals: band of 'maxgap' diags
 *      -**OO0oo---   O: left-hand "gap" diagonals: band of 'maxgap' diags
 *      --**OO0oo--   *: "extra" diagonals: additional band of 'maxgap' diags
 *      ---**OO0oo-   +: cells not in any of the above
 *      ----**OO0oo
 *            |-|
 *   Gotta touch one of these diags
 *
 * Basically, the filtering step removes alignments that do not at some point
 * touch a cell labeled '0' or 'O' in the diagram above.
 *
 */
bool SwAligner::backtraceNucleotidesLocalSseU8(
	TAlScore       escore, // in: expected score
	SwResult&      res,    // out: store results (edits and scores) here
	size_t&        off,    // out: store diagonal projection of origin
	size_t&        nbts,   // out: # backtracks
	size_t         row,    // start in this row
	size_t         col,    // start in this column
	RandomSource&  rnd)    // random gen, to choose among equal paths
{
	assert_lt(row, dpRows());
	assert_lt(col, (size_t)(rff_ - rfi_));
	SSEData& d = fw_ ? sseU8fw_ : sseU8rc_;
	SSEMetrics& met = extend_ ? sseU8ExtendMet_ : sseU8MateMet_;
	met.bt++;
	// Loop invariants, kept in locals (the build has -fno-strict-aliasing,
	// so the compiler would otherwise reload them after every store): the
	// vector matrix's geometry in bytes (BT_SET_ROW_COL), the scoring
	// parameters, the read, quality and reference, and the mask, reset and
	// dirty-interval arrays (resetRow changes a row's contents, not these)
	const size_t btnvecrow = d.mat_.nvecrow_;
	const size_t btcolstride = d.mat_.colstride_;
	uint8_t * const btmatbuf = (uint8_t*)d.mat_.matbuf_.ptr();
	const size_t btLeft = btcolstride * NBYTES_PER_REG;
	const size_t btUp = ROWSTRIDE * NBYTES_PER_REG;
	const size_t btUpWrap = (btcolstride - ROWSTRIDE) * NBYTES_PER_REG - 1;
	const size_t btF = SSEMatrix::F * NBYTES_PER_REG;
	const size_t btH = SSEMatrix::H * NBYTES_PER_REG;
	assert_eq(0, SSEMatrix::E);
	EList<uint16_t> * const btmasks = d.mat_.masks_.ptr();
	const bool * const btreset = d.mat_.reset_.ptr();
	size_t * const btdirtylo = d.mat_.dirtylo_.ptr();
	size_t * const btdirtyhi = d.mat_.dirtyhi_.ptr();
	const size_t btnrow = d.mat_.nrow();
	const size_t gapbar = (size_t)sc_->gapbar;
	const int rdgapo = sc_->readGapOpen();
	const int rdgape = sc_->readGapExtend();
	const int rfgapo = sc_->refGapOpen();
	const int rfgape = sc_->refGapExtend();
	const int64_t bonus = sc_->match(30);
	const Scoring * const sc = sc_;
	const char * const rdbuf = rd_->buf();
	const char * const qubuf = qu_->buf();
	const char * const rfbuf = rf_ + rfi_;
	const size_t rdoff = rdi_;
	size_t ncell = 0; // cells visited, added to met.btcell and nbts on return
	// Branch points the parent code pushed on its backtrack stack (btnstack_).
	// A frame's cell is marked reported through before the frame is pushed,
	// so a popped frame re-enters a reported cell, which pops the next frame:
	// once a walk reaches a cell it cannot move through, it fails after one
	// more cell visit per frame, and no frame ever yields a second path.
	// Only their number is kept.
	size_t nframes = 0;
	assert(!d.profbuf_.empty());
	assert_lt(row, rd_->length());
	btcells_.clear();  // empty the cells-so-far list
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
	size_t ct = SSEMatrix::H; // cell type
	// The current cell in the SSE vector matrix (BT_SET_ROW_COL)
	size_t btrowvec;
	uint8_t *btcur;
	BT_SET_ROW_COL(row, col);
	while((int)row >= 0) {
		ncell++;
		int readc = rdbuf[rdoff + row];
		int refm  = (int)rfbuf[col];
		int readq = qubuf[row];
		assert_eq(readc, (*rd_)[rdi_ + row]);
		assert_eq(refm, (int)rf_[rfi_ + col]);
		assert_eq(readq, (*qu_)[row]);
		assert_leq(col, origCol);
		// The cell above; used only when row > 0
		const bool btwrap = (btrowvec == 0);
		uint8_t * const btup = btwrap ? btcur + (ptrdiff_t)btUpWrap : btcur - btUp; // (size_t)-1 when nvecrow_ is 1
		assert_eq(btcur[0], ((TCScore*)(d.mat_.matbuf_.ptr() + col * btcolstride + btrowvec * ROWSTRIDE))[row / btnvecrow]);
		assert(row == 0 || btup[btH] == d.mat_.helt(row - 1, col));
		// Get score in this cell
		bool empty = false;
		bool reportedThru, canMoveThru, branch = false;
		int cur = SSEMatrix::H;
		if(!btreset[row]) {
			d.mat_.resetRow(row);
		}
		// This cell's mask word (see SSEMatrix).  Bit 0, reported through,
		// is the only bit any code reads.  This backtrace writes a word only
		// when it enters a cell whose bit 0 is clear, and then writes 1, so
		// a word is 0 or 1.  The parent code also stored the remaining H, E
		// and F masks in bits 1-12 and read them back only when bit 0 was
		// clear, where its words were 0 too.
		uint16_t& mk = btmasks[row][col];
		assert_eq(&mk, &d.mat_.masks_[row][col]);
		reportedThru = (mk & (1 << 0)) != 0;
		assert(reportedThru || mk == 0);
		canMoveThru = true;
		if(reportedThru) {
			canMoveThru = false;
		} else {
			empty = false;
			if(row > 0) {
				assert_gt(row, 0);
				size_t rowFromEnd = btnrow - row - 1;
				bool gapsAllowed = true;
				if(row < gapbar ||
				   rowFromEnd < gapbar)
				{
					gapsAllowed = false;
				}
				const int floorsc = 0;
				const int offsetsc = 0;
				// Move to beginning of column/row
				if(ct == SSEMatrix::E) { // AKA rdgap
					assert_gt(col, 0);
					const uint8_t *btleft = btcur - btLeft;
					TAlScore sc_cur = ((TCScore)btcur[0]) + offsetsc;
					assert(gapsAllowed);
					// Currently in the E matrix; incoming transition must come from the
					// left.  It's either a gap open from the H matrix or a gap extend from
					// the E matrix.
					// TODO: save and restore origMask as well as mask
					int origMask = 0, mask = 0;
					// Get H score of cell to the left
					TAlScore sc_h_left = ((TCScore)btleft[btH]) + offsetsc;
					if(sc_h_left > 0 && sc_h_left - rdgapo == sc_cur) {
						mask |= (1 << 0);
					}
					// Get E score of cell to the left
					TAlScore sc_e_left = ((TCScore)btleft[0]) + offsetsc;
					if(sc_e_left > 0 && sc_e_left - rdgape == sc_cur) {
						mask |= (1 << 1);
					}
					origMask = mask;
					assert(origMask > 0 || sc_cur <= sc_->match());
					if(mask == 3) {
						// Pick H -> E cell
						cur = SW_BT_OALL_READ_OPEN;
						branch = true;
					} else if(mask == 2) {
						// I chose the E cell
						cur = SW_BT_RDGAP_EXTEND;
					} else if(mask == 1) {
						// I chose the H cell
						cur = SW_BT_OALL_READ_OPEN;
					} else {
						empty = true;
						// It's empty, so the only question left is whether we should be
						// allowed in terimnate in this cell.  If it's got a valid score
						// then we *shouldn't* be allowed to terminate here because that
						// means it's part of a larger alignment that was already reported.
						canMoveThru = (origMask == 0);
					}
					assert(!empty || !canMoveThru);
				} else if(ct == SSEMatrix::F) { // AKA rfgap
					assert_gt(row, 0);
					assert(gapsAllowed);
					TAlScore sc_h_up = ((TCScore)btup[btH]) + offsetsc;
					TAlScore sc_f_up = ((TCScore)btup[btF]) + offsetsc;
					TAlScore sc_cur  = ((TCScore)btcur[btF]) + offsetsc;
					// Currently in the F matrix; incoming transition must come from above.
					// It's either a gap open from the H matrix or a gap extend from the F
					// matrix.
					// TODO: save and restore origMask as well as mask
					int origMask = 0, mask = 0;
					// Get H score of cell above
					if(sc_h_up > floorsc && sc_h_up - rfgapo == sc_cur) {
						mask |= (1 << 0);
					}
					// Get F score of cell above
					if(sc_f_up > floorsc && sc_f_up - rfgape == sc_cur) {
						mask |= (1 << 1);
					}
					origMask = mask;
					assert(origMask > 0 || sc_cur <= sc_->match());
					if(mask == 3) {
						// I chose the H cell
						cur = SW_BT_OALL_REF_OPEN;
						branch = true;
					} else if(mask == 2) {
						// I chose the F cell
						cur = SW_BT_RFGAP_EXTEND;
					} else if(mask == 1) {
						// I chose the H cell
						cur = SW_BT_OALL_REF_OPEN;
					} else {
						empty = true;
						// It's empty, so the only question left is whether we should be
						// allowed in terimnate in this cell.  If it's got a valid score
						// then we *shouldn't* be allowed to terminate here because that
						// means it's part of a larger alignment that was already reported.
						canMoveThru = (origMask == 0);
					}
					assert(!empty || !canMoveThru);
				} else {
					assert_eq(SSEMatrix::H, ct);
					TAlScore sc_cur      = ((TCScore)btcur[btH]) + offsetsc;
					TAlScore sc_f_up     = ((TCScore)btup[btF])  + offsetsc;
					TAlScore sc_h_up     = ((TCScore)btup[btH])  + offsetsc;
					TAlScore sc_h_left   = col > 0 ? (((TCScore)(btcur - btLeft)[btH]) + offsetsc) : floorsc;
					TAlScore sc_e_left   = col > 0 ? (((TCScore)(btcur - btLeft)[0])   + offsetsc) : floorsc;
					TAlScore sc_h_upleft = col > 0 ? (((TCScore)(btup  - btLeft)[btH]) + offsetsc) : floorsc;
					TAlScore sc_diag     = sc->score(readc, refm, readq - 33);
					// TODO: save and restore origMask as well as mask
					int origMask = 0, mask = 0;
					if(gapsAllowed) {
						if(sc_h_up     > floorsc && sc_cur == sc_h_up   - rfgapo) {
							mask |= (1 << 0);
						}
						if(sc_h_left   > floorsc && sc_cur == sc_h_left - rdgapo) {
							mask |= (1 << 1);
						}
						if(sc_f_up     > floorsc && sc_cur == sc_f_up   - rfgape) {
							mask |= (1 << 2);
						}
						if(sc_e_left   > floorsc && sc_cur == sc_e_left - rdgape) {
							mask |= (1 << 3);
						}
					}
					if(sc_h_upleft > floorsc && sc_cur == sc_h_upleft + sc_diag) {
						mask |= (1 << 4);
					}
					origMask = mask;
					assert(origMask > 0 || sc_cur <= sc_->match());
					assert(gapsAllowed || mask == (1 << 4) || mask == 0);
					if(mask != 0) {
						// The move and the mask left behind, from the table
						// (see btHSelect)
						int sel = btHSelect[mask];
						cur = sel & 7;
						int rest = sel >> 3;
						assert_eq(rest, mask & ~btHSelectBit(mask));
						assert_eq(rest != 0, alts5[mask] > 1);
						assert_eq(cur, btHSelectBit(mask) == 16 ? (int)SW_BT_OALL_DIAG :
						               btHSelectBit(mask) ==  1 ? (int)SW_BT_OALL_REF_OPEN :
						               btHSelectBit(mask) ==  2 ? (int)SW_BT_OALL_READ_OPEN :
						               btHSelectBit(mask) ==  4 ? (int)SW_BT_RFGAP_EXTEND :
						                                          (int)SW_BT_RDGAP_EXTEND);
						if(rest != 0) {
							branch = true;
						}
					} else {
						empty = true;
						// It's empty, so the only question left is whether we should be
						// allowed in terimnate in this cell.  If it's got a valid score
						// then we *shouldn't* be allowed to terminate here because that
						// means it's part of a larger alignment that was already reported.
						canMoveThru = (origMask == 0);
					}
				}
				assert(!empty || !canMoveThru || ct == SSEMatrix::H);
			}
		}
		//cerr << "reportedThrough rejected (" << row << ", " << col << ")" << endl;
		if(!reportedThru) {
			// Mark the cell reported through, and record it in the row's
			// dirty interval (SSEMatrix::markDirty).  A cell already
			// reported through was marked in this fill, after its row was
			// reset, so its word is 1 and its column is in the interval.
			mk = (1 << 0);
			btdirtylo[row] = (col < btdirtylo[row]) ? col : btdirtylo[row];
			btdirtyhi[row] = (col >= btdirtyhi[row]) ? col + 1 : btdirtyhi[row];
		}
		assert(d.mat_.reportedThrough(row, col));
		assert_eq(gaps, Edit::numGaps(ned));
		assert_leq(gaps, rdgap_ + rfgap_);
		// Cell was involved in a previously-reported alignment?
		if(!canMoveThru) {
			// Give up.  The parent code first popped each frame of its
			// backtrack stack and re-entered the frame's cell, which it had
			// marked reported through before pushing the frame, so each pop
			// counted one more cell and led to the next pop.
			ncell += nframes;
			met.btcell += ncell;
			nbts += ncell;
			res.reset();
			met.btfail++; // DP backtraces failed
			return false;
		}
		assert(!reportedThru);
		assert(!sc_->monotone || score.score() >= minsc_);
		if(empty || row == 0) {
			assert_eq(SSEMatrix::H, ct);
			btcells_.expand();
			btcells_.back().first = row;
			btcells_.back().second = col;
			// This cell is at the end of a legitimate alignment
			trimBeg = row;
			assert_eq(btcells_.size(), dpRows() - trimBeg - trimEnd + readGaps);
			break;
		}
		if(branch) {
			// The parent code added a frame to its backtrack stack here
			nframes++;
		}
		btcells_.expand();
		btcells_.back().first = row;
		btcells_.back().second = col;
		switch(cur) {
			// Move up and to the left.  If the reference nucleotide in the
			// source row mismatches the read nucleotide, penalize
			// it and add a nucleotide mismatch.
			case SW_BT_OALL_DIAG: {
				assert_gt(row, 0); assert_gt(col, 0);
				int readC = rdbuf[row];
				int refNmask = (int)rfbuf[col];
				assert_eq(readC, (*rd_)[row]);
				assert_gt(refNmask, 0);
				int m = matchesEx(readC, refNmask);
				ct = SSEMatrix::H;
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
					score.score_ -= pen;
					assert(!sc_->monotone || score.score() >= escore);
				} else {
					// Reward a match
					score.score_ += bonus;
					assert(!sc_->monotone || score.score() >= escore);
				}
				if(m == -1) {
					score.ns_++;
				}
				row--; col--;
				btcur = btup - btLeft;
				btrowvec = btwrap ? btnvecrow - 1 : btrowvec - 1;
				assert(VALID_AL_SCORE(score));
				break;
			}
			// Move up.  Add an edit encoding the ref gap.
			case SW_BT_OALL_REF_OPEN:
			{
				assert_gt(row, 0);
				Edit e(
					(int)row,
					'-',
					"ACGTN"[(int)rdbuf[row]],
					EDIT_TYPE_REF_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				row--;
				ct = SSEMatrix::H;
				int pen = rfgapo;
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; refGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				btcur = btup;
				btrowvec = btwrap ? btnvecrow - 1 : btrowvec - 1;
				break;
			}
			// Move up.  Add an edit encoding the ref gap.
			case SW_BT_RFGAP_EXTEND:
			{
				assert_gt(row, 1);
				Edit e(
					(int)row,
					'-',
					"ACGTN"[(int)rdbuf[row]],
					EDIT_TYPE_REF_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				row--;
				ct = SSEMatrix::F;
				int pen = rfgape;
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; refGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				btcur = btup;
				btrowvec = btwrap ? btnvecrow - 1 : btrowvec - 1;
				break;
			}
			case SW_BT_OALL_READ_OPEN:
			{
				assert_gt(col, 0);
				Edit e(
					(int)row+1,
					mask2dna[(int)rfbuf[col]],
					'-',
					EDIT_TYPE_READ_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				col--;
				ct = SSEMatrix::H;
				int pen = rdgapo;
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; readGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				btcur -= btLeft;
				break;
			}
			case SW_BT_RDGAP_EXTEND:
			{
				assert_gt(col, 1);
				Edit e(
					(int)row+1,
					mask2dna[(int)rfbuf[col]],
					'-',
					EDIT_TYPE_READ_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				col--;
				ct = SSEMatrix::E;
				int pen = rdgape;
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; readGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				btcur -= btLeft;
				break;
			}
			default: throw 1;
		}
	} // while((int)row > 0)
	met.btcell += ncell;
	nbts += ncell;
	assert_geq(col, 0);
	assert_eq(SSEMatrix::H, ct);
	// The number of cells in the backtracs should equal the number of read
	// bases after trimming plus the number of gaps
	assert_eq(btcells_.size(), dpRows() - trimBeg - trimEnd + readGaps);
	// Check whether we went through a core diagonal and set 'reported' flag on
	// each cell
	bool overlappedCoreDiag = false;
	for(size_t i = 0; i < btcells_.size(); i++) {
		size_t rw = btcells_[i].first;
		size_t cl = btcells_[i].second;
		// Calculate the diagonal within the *trimmed* rectangle, i.e. the
		// rectangle we dealt with in align, gather and backtrack.
		int64_t diagi = cl - rw;
		// Now adjust to the diagonal within the *untrimmed* rectangle by
		// adding on the amount trimmed from the left.
		diagi += rect_->triml;
		if(diagi >= 0) {
			size_t diag = (size_t)diagi;
			if(diag >= rect_->corel && diag <= rect_->corer) {
				overlappedCoreDiag = true;
				break;
			}
		}
#ifndef NDEBUG
		//assert(!d.mat_.reportedThrough(rw, cl));
		//d.mat_.setReportedThrough(rw, cl);
		assert(d.mat_.reportedThrough(rw, cl));
#endif
	}
	if(!overlappedCoreDiag) {
		// Must overlap a core diagonal.  Otherwise, we run the risk of
		// reporting an alignment that overlaps (and trumps) a higher-scoring
		// alignment that lies partially outside the dynamic programming
		// rectangle.
		res.reset();
		met.corerej++;
		return false;
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
	if(score.ns_ > nceil_) {
		// Alignment has too many Ns in it!
		res.reset();
		met.nrej++;
		return false;
	}
	res.reverse();
	assert(Edit::repOk(ned, (*rd_)));
	assert_eq(score.score(), escore);
	assert_leq(gaps, rdgap_ + rfgap_);
	off = col;
	assert_lt(col + (size_t)rfi_, (size_t)rff_);
	score.gaps_ = gaps;
	score.edits_ = (int)ned.size();
	score.basesAligned_ = (int)(rdf_ - rdi_ - trimBeg - trimEnd - score.edits_);
	res.alres.setScore(score);
	res.alres.setShape(
		refidx_,                  // ref id
		off + rfi_ + rect_->refl, // 0-based ref offset
		reflen_,                  // reference length
		fw_,                      // aligned to Watson?
		rdf_ - rdi_,              // read length
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
	assert(Edit::repOk(ned, (*rd_), true, trimBeg, trimEnd));
	assert(res.repOk());
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
	Edit::toRef((*rd_), ned, editstr, true, trimBeg, trimEnd);
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
	met.btsucc++; // DP backtraces succeeded
	return true;
}
