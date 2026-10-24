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

// In end-to-end mode, we start high (255) and go low (0).  Factoring in
// a query profile involves unsigned saturating subtraction, so all the
// query profile elements should be expressed as a positive penalty rather
// than a negative score.

typedef uint8_t TCScore;

/**
 * Build query profile look up tables for the read.  The query profile look
 * up table is organized as a 1D array indexed by [i][j] where i is the
 * reference character in the current DP column (0=A, 1=C, etc), and j is
 * the segment of the query we're currently working on.
 */
void SwAligner::buildQueryProfileEnd2EndSseU8(bool fw) {
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
	d.bias_ = 0; // no bias needed for end-to-end alignment; just use subtraction
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
	// Query profile: sc_->score(readc, 1 << refc, readq - 33), negated to a
	// positive penalty that fits in an unsigned byte; 0 past the read's end.
	// The scorer's qualProf table holds the five bytes for every reference
	// character at once, and the largest penalty among them.
	const char *rdbuf = rd->buf();
	const char *qubuf = qu->buf();
	size_t maxPen = 0;
	for(size_t i = 0; i < seglen; i++) {
		uint8_t *p = prof + i * segBytes;
		size_t j = i, k = 0;
		for(; j < len; j += seglen, k++) {
			const int readc = rdbuf[j];
			const int readq = qubuf[j];
			assert_range(0, 4, readc);
			assert_range(0, 255, readq - 33);
			const Scoring::QualProf& e = sc_->qualProf[readc][readq - 33];
			const uint64_t w = e.ee;
			assert_eq((uint8_t)w,         (uint8_t)-sc_->score(readc, 1, readq - 33));
			assert_eq((uint8_t)(w >> 8),  (uint8_t)-sc_->score(readc, 2, readq - 33));
			assert_eq((uint8_t)(w >> 16), (uint8_t)-sc_->score(readc, 4, readq - 33));
			assert_eq((uint8_t)(w >> 24), (uint8_t)-sc_->score(readc, 8, readq - 33));
			assert_eq((uint8_t)(w >> 32), (uint8_t)-sc_->score(readc, 16, readq - 33));
			maxPen = std::max(maxPen, (size_t)e.eeMax);
			p[k] = (uint8_t)w;
			p[k + refBytes] = (uint8_t)(w >> 8);
			p[k + 2 * refBytes] = (uint8_t)(w >> 16);
			p[k + 3 * refBytes] = (uint8_t)(w >> 24);
			p[k + 4 * refBytes] = (uint8_t)(w >> 32);
		}
		for(; k < NWORDS_PER_REG; k++) {
			for(size_t refc = 0; refc < ALPHA_SIZE; refc++) {
				p[k + refc * refBytes] = 0;
			}
		}
	}
	d.maxPen_ = maxPen;
}

#ifndef NDEBUG
/**
 * Return true iff the cell has sane E/F/H values w/r/t its predecessors.
 */
static bool cellOkEnd2EndU8(
	SSEData& d,
	size_t row,
	size_t col,
	int refc,
	int readc,
	int readq,
	const Scoring& sc)     // scoring scheme
{
	TCScore floorsc = 0;
	TAlScore ceilsc = MAX_I64;
	TAlScore offsetsc = -0xff;
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
			// Check that we couldn't have got a better E score
			assert_geq(sc_e_cur, sc_e_left - sc.readGapExtend());
			assert_geq(sc_e_cur, sc_h_left - sc.readGapOpen());
		}
		if(row > 0 && sc_f_cur > floorsc && sc_f_cur <= ceilsc) {
			sc_h_up = d.mat_.helt(row-1, col) + offsetsc;
			sc_f_up = d.mat_.felt(row-1, col) + offsetsc;
			f_up_trans = (sc_f_up > floorsc && sc_f_cur == sc_f_up - sc.refGapExtend());
			h_up_trans = (sc_h_up > floorsc && sc_f_cur == sc_h_up - sc.refGapOpen());
			assert(f_up_trans || h_up_trans);
			// Check that we couldn't have got a better F score
			assert_geq(sc_f_cur, sc_f_up - sc.refGapExtend());
			assert_geq(sc_f_cur, sc_h_up - sc.refGapOpen());
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
TAlScore SwAligner::alignGatherEE8(int& flag, bool debug) {
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
	buildQueryProfileEnd2EndSseU8(fw_);
	assert(!d.profbuf_.empty());

	assert_eq(0, d.maxBonus_);
	size_t iter =
		(dpRows() + (NWORDS_PER_REG-1)) / NWORDS_PER_REG; // iter = segLen
	
	// Now set up the score vectors.  We just need two columns worth, which
	// we'll call "left" and "right".
	d.vecbuf_.resize(4 * 2 * iter);
	d.vecbuf_.zero();
	SSERegI *vbuf_l = d.vecbuf_.ptr();
	SSERegI *vbuf_r = d.vecbuf_.ptr() + (4 * iter);

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
		false,         // alignment is local?
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
	SSERegI ve       = sse_setzero_siall();
	SSERegI vf       = sse_setzero_siall();
	SSERegI vh       = sse_setzero_siall();
	SSERegI vhd      = sse_setzero_siall();
	SSERegI vhdtmp   = sse_setzero_siall();
	SSERegI vtmp     = sse_setzero_siall();
	SSERegI vzero    = sse_setzero_siall();
	SSERegI vhilsw   = sse_setzero_siall();

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

#ifndef NDEBUG
	SSERegI vhi      = sse_setzero_siall();
	sse_fill_u8_opt(0xff, vhi)
#endif
	sse_fill_u8_opt(0, vlo);

	// vhilsw: topmost (least sig) word set to 0xff, all other words=0
	sse_set_low_u8(0xff, vhilsw);
	
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
	
	// Maximum score in final row
	bool found = false;
	TCScore lrmax = MIN_U8;
	
	for(size_t i = 0; i < iter; i++) {
		sse_store_siall(pvERight, vlo); pvERight += ROWSTRIDE_2COL;
		// Could initialize Hs to high or low.  If high, cells in the lower
		// triangle will have somewhat more legitiate scores, but still won't
		// be exhaustively scored.
		sse_store_siall(pvHRight, vlo); pvHRight += ROWSTRIDE_2COL;
	}
	
	assert_gt(sc_->gapbar, 0);
	size_t nfixup = 0;

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
	
	for(size_t i = (size_t)rfi_; i < (size_t)rff_; i++) {
		// Swap left and right; vbuf_l is the vector on the left, which we
		// generally load from, and vbuf_r is the vector on the right, which we
		// generally store to.
		swap(vbuf_l, vbuf_r);
		pvELeft = vbuf_l + 0; pvERight = vbuf_r + 0;
		/* pvFLeft = vbuf_l + 1; */ pvFRight = vbuf_r + 1;
		pvHLeft = vbuf_l + 2; pvHRight = vbuf_r + 2;
		
		// Fetch the appropriate query profile.  Note that elements of rf_ must
		// be numbers, not masks.
		const int refc = (int)rf_[i];
		
		// Fetch the appropriate query profile
		size_t off = (size_t)firsts5[refc] * iter * 2;
		pvScore = d.profbuf_.ptr() + off; // even elts = query profile, odd = gap barrier
		
		// Set all cells to low value
		vf = sse_xor_siall(vf, vf);

		// Load H vector from the final row of the previous column
		vh = sse_load_siall(pvHLeft + colstride - ROWSTRIDE_2COL);
		// Shift 2 bytes down so that topmost (least sig) cell gets 0
		vh = sse_slli_siall(vh, NBYTES_PER_WORD);
		// Fill topmost (least sig) cell with high value
		vh = sse_or_siall(vh, vhilsw);
		
		// For each character in the reference text:
		size_t j;
		for(j = 0; j < iter; j++) {
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
			vh = sse_subs_epu8(vh, pvScore[0]);
			
			// Update H, factoring in E and F
			vh = sse_max_epu8(vh, vf);
			
			// Update vE value
			vhdtmp = vhd;
			vhd = sse_subs_epu8(vhd, rdgapo);
			vhd = sse_subs_epu8(vhd, pvScore[1]); // veto some read gap opens
			ve = sse_subs_epu8(ve, rdgape);
			ve = sse_max_epu8(ve, vhd);
			vh = sse_max_epu8(vh, ve);
			
			// Save the new vH values
			sse_store_siall(pvHRight, vh);
			pvHRight += ROWSTRIDE_2COL;
			vtmp = vh;
			assert_all_lt(ve, vhi);
			
			// Load the next h value
			vh = vhdtmp;
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
		
		// Check in the last row for the maximum so far
		SSERegI *vtmp = vbuf_r + 2 /* H */ + (d.lastIter_ * ROWSTRIDE_2COL);
		// Note: we may not want to extract from the final row
		TCScore lr = ((TCScore*)(vtmp))[d.lastWord_];
		found = true;
		if(lr > lrmax) {
			lrmax = lr;
		}
		
		// Now we'd like to know whether the bottommost element of the right
		// column is a candidate we might backtrace from.  First question is:
		// did it exceed the minimum score threshold?
		TAlScore score = (TAlScore)(lr - 0xff);
		if(lr == MIN_U8) {
			score = MIN_I64;
		}
		if(!debug && score >= minsc_) {
			DpBtCandidate cand(dpRows() - 1, i - rfi_, score);
			btdiag_.add(i - rfi_, cand);
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
							if(h_sc == 0) h_sc = MIN_I16;
							else h_sc -= 0xff;
							if(e_sc == 0) e_sc = MIN_I16;
							else e_sc -= 0xff;
							if(f_sc == 0) f_sc = MIN_I16;
							else f_sc -= 0xff;
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
					} // end of loop over anti-diags
					if(done) {
						break;
					}
					idx++;
					idxrow += cper_.nrow_;
				}
			} else {
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
					if(h_sc == 0) h_sc = MIN_I16;
					else h_sc -= 0xff;
					//if(e_sc == 0) e_sc = MIN_I16;
					//else e_sc -= 0xff;
					if(f_sc == 0) f_sc = MIN_I16;
					else f_sc -= 0xff;
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
	}
	
	// Update metrics
	if(!debug) {
		size_t ninner = (rff_ - rfi_) * iter;
		met.col   += (rff_ - rfi_);             // DP columns
		met.cell  += (ninner * NWORDS_PER_REG); // DP cells
		met.inner += ninner;                    // DP inner loop iters
		met.fixup += nfixup;                    // DP fixup loop iters
	}

	flag = 0;

	// Did we find a solution?
	TAlScore score = MIN_I64;
	if(!found) {
		flag = -1; // no
		if(!debug) met.dpfail++;
		return MIN_I64;
	} else {
		score = (TAlScore)(lrmax - 0xff);
		if(score < minsc_) {
			flag = -1; // no
			if(!debug) met.dpfail++;
			return score;
		}
	}
	
	// Could we have saturated?
	if(lrmax == MIN_U8) {
		flag = -2; // yes
		if(!debug) met.dpsat++;
		return MIN_I64;
	}

	// Now take all the backtrace candidates in the btdaig_ structure and
	// dump them into the btncand_ array.  They'll be sorted later.
	if(!debug) {
		btdiag_.dump(btncand_);
		assert(!btncand_.empty());
	}
	
	// Return largest score
	if(!debug) met.dpsucc++;
	return score;
}

#if NBYTES_PER_REG == 32
/**
 * The fill of alignNucleotidesEnd2EndSseU8 for a column of ITER = 2 to 4
 * 256-bit vectors (33 to 128 rows).  It computes the general loop's saturating
 * u8 recurrence (see the comments there) and stores the same E, F and H
 * vectors in the same places.  The general loop's speed is set by the chain
 * of dependent instructions from one column to the next: the first pass's F
 * through every row, a scan with a shift across the 128-bit halves in every
 * step, and the fix-up pass.  Here that chain is the scan, with in-half
 * shifts, and a few instructions, and every value passed from one column to
 * the next stays in a register.
 *
 * Row j of lane L is read row L * ITER + j.  For each column:
 *  - H0[j] = max(diag - profile, E[j]), H before the reference gaps (F);
 *  - within a lane, F flows from row j to row j + 1 at the cost of rfgapo
 *    from H0 (opening from H = max(H0, F) costs no less than extending F, as
 *    rfgapo >= rfgape) or of the barrier and rfgape from F; HL and FL are H
 *    and F with these carries only;
 *  - the carry c into each lane from the lanes above is the general loop's
 *    prefix scan of the lanes' carry-outs, done within each 128-bit half (two
 *    steps of four lanes, with shifts that do not cross the halves) and then
 *    across the halves once;
 *  - each row takes c less the cost from the top of its lane (vpen): F =
 *    max(FL, that), H = max(HL, that), and the next column's E = max(E -
 *    rdgape, H - rdgapo - barrier).
 * The next column's carry-out, the scan's input, is max_j(H0[j] - cpen[j]).
 * Subtraction distributes over max, so it is also the max of terms in HL
 * and E, which are ready while the scan runs, and of terms in c and in c
 * shifted down a lane, each with its smallest cost (zc, wc, kb).  The chain
 * from one column to the next is then the scan and one subtraction.  All
 * arithmetic is saturating u8, and sat(sat(x - a) - b) = sat(x - sat(a + b)),
 * so costs may be summed before they are subtracted.  Returns the maximum of
 * H over the final row.
 */
template<size_t ITER>
static uint8_t alignEE8V256Short(
	const __m256i *prof,    // query profile: per reference character, ITER (score, barrier) pairs
	__m256i *mat,           // matrix: column 0, row 0, E (SSEMatrix layout)
	size_t colstride,       // vectors per matrix column
	const char *rf,         // reference
	size_t rfi,             // first column
	size_t rff,             // last column (excl.)
	const Scoring& sc,      // gap penalties
	size_t lastWordIdx)     // final row: byte offset from a column's first H vector
{
	const size_t E = SSEMatrix::E, F = SSEMatrix::F, H = SSEMatrix::H;
	const __m256i zero = _mm256_setzero_si256();
	const __m256i vrfgapo = _mm256_set1_epi8((char)sc.refGapOpen());
	const __m256i vrfgape = _mm256_set1_epi8((char)sc.refGapExtend());
	const __m256i vrdgapo = _mm256_set1_epi8((char)sc.readGapOpen());
	const __m256i vrdgape = _mm256_set1_epi8((char)sc.readGapExtend());
	// Gap barrier rows (the same for every reference character), and the costs
	// derived from them, per vector row j:
	//  bpe:  F carried from row j to row j + 1 (barrier, then rfgape);
	//  cpen: H0 at row j to the lane's carry-out (rfgapo, then bpe below);
	//  vpen: the carry into the lane to row j (as the general loop's vpen);
	//  epen: H to the next column's E (rdgapo, then barrier).
	__m256i bar[ITER], bpe[ITER], cpen[ITER], vpen[ITER], epen[ITER];
	__m256i vbar = zero;
	for(size_t j = 0; j < ITER; j++) {
		bar[j] = _mm256_load_si256(prof + 2 * j + 1);
		bpe[j] = _mm256_adds_epu8(bar[j], vrfgape);
		epen[j] = _mm256_adds_epu8(vrdgapo, bar[j]);
		vpen[j] = j == 0 ? bar[0] : _mm256_adds_epu8(_mm256_adds_epu8(vpen[j-1], vrfgape), bar[j]);
		vbar = _mm256_or_si256(vbar, bar[j]);
	}
	cpen[ITER-1] = vrfgapo;
	for(size_t j = ITER-1; j > 0; j--) {
		cpen[j-1] = _mm256_adds_epu8(cpen[j], bpe[j]);
	}
	// The scan's costs: passing lane L whole is pass(L) (as the general loop's
	// vpass); qN(L) is the cost of passing the N lanes that end at L, within
	// a half (lanes outside the half are never read).  chi(L) is the cost of
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
	// Byte 15 of the lower half is 0xff: shifted into lane 0, it is row 0's
	// diagonal
	const __m256i vtop = _mm256_set_epi64x(0, 0, (long long)0xff00000000000000ULL, 0);
	// Costs to the next column's carry-out, max_j(H0'[j] - cpen[j]), where
	// H0'[j] = max(diag[j] - profile'[j], E'[j]), E'[j] = max(E[j] - rdgape,
	// H[j] - epen[j]), H[j] = max(HL[j], c - vpen[j]) and diag[j] = H[j - 1]
	// or, for j = 0, H[ITER - 1] shifted down a lane:
	//  rc[j]: from E[j] (rdgape, cpen); ec[j]: from H[j] (epen, cpen);
	//  pc[r][j]: from diag[j], for reference character r (profile, cpen);
	//  zc[r]: from c, the least over its paths; wc[r]: from c shifted down a
	//  lane; kb[r]: from lane 15 of the scan, carried into the upper half,
	//  for both (c and c shifted are taken from the scan's values directly).
	__m256i rc[ITER], ec[ITER], pc[5][ITER], zc[5], wc[5], kb[5];
	// chi for c shifted down a lane: lane 16 takes lane 14's value, with no
	// carry from lane 15
	const __m256i chish = _mm256_or_si256(_mm256_slli_si256(chi, 1),
		_mm256_set_epi64x(0, 0xff, 0, 0));
	const __m256i vpenSh = _mm256_alignr_epi8(vpen[ITER-1],
		_mm256_permute2x128_si256(vpen[ITER-1], vpen[ITER-1], 0x08), 15);
	for(size_t j = 0; j < ITER; j++) {
		rc[j] = _mm256_adds_epu8(vrdgape, cpen[j]);
		ec[j] = _mm256_adds_epu8(epen[j], cpen[j]);
	}
	for(size_t r = 0; r < 5; r++) {
		__m256i z = _mm256_adds_epu8(vpen[0], ec[0]);
		for(size_t j = 0; j < ITER; j++) {
			pc[r][j] = _mm256_adds_epu8(_mm256_load_si256(prof + (r * ITER + j) * 2), cpen[j]);
			if(j > 0) {
				z = _mm256_min_epu8(z, _mm256_adds_epu8(vpen[j-1], pc[r][j]));
				z = _mm256_min_epu8(z, _mm256_adds_epu8(vpen[j], ec[j]));
			}
		}
		zc[r] = z;
		wc[r] = _mm256_adds_epu8(vpenSh, pc[r][0]);
		kb[r] = _mm256_min_epu8(_mm256_adds_epu8(chi, zc[r]), _mm256_adds_epu8(chish, wc[r]));
	}
	// lpen[j]: H0 at row j to the carry into the lane's last row
	__m256i lpen[ITER-1];
	lpen[ITER-2] = vrfgapo;
	for(size_t j = ITER-2; j > 0; j--) {
		lpen[j-1] = _mm256_adds_epu8(lpen[j], bpe[j]);
	}
	// H0 and E of the column, and its carry-out: for the first, from the
	// staged initial column (zero), whose diagonal is 0xff in row 0 only
	__m256i h0[ITER], ecur[ITER], co;
	{
		const __m256i *pvS = prof + (size_t)firsts5[(int)rf[rfi]] * ITER * 2;
		const __m256i p = _mm256_permute2x128_si256(zero, vtop, 0x02);
		const __m256i d0 = _mm256_alignr_epi8(zero, p, 15);
		h0[0] = _mm256_subs_epu8(d0, _mm256_load_si256(pvS));
		co = _mm256_subs_epu8(h0[0], cpen[0]);
		for(size_t j = 1; j < ITER; j++) {
			h0[j] = zero;
		}
		for(size_t j = 0; j < ITER; j++) {
			ecur[j] = zero;
		}
	}
	uint8_t lrmax = 0;
	__m256i *pvCol = mat;
	for(size_t i = rfi; i < rff; i++) {
		// The next column's reference character.  (After the last column,
		// the values computed with it are not used.)
		const size_t rn = (size_t)firsts5[(int)rf[i + 1 < rff ? i + 1 : i]];
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
		// F and H with the carries from within the lane only
		__m256i fl[ITER], hl[ITER];
		fl[0] = zero;
		hl[0] = h0[0];
		{
			__m256i vf = _mm256_subs_epu8(h0[0], vrfgapo);
			for(size_t j = 1; j < ITER; j++) {
				fl[j] = _mm256_subs_epu8(vf, bar[j]);
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
		// ... and its other terms, the latest (from the last row of HL) last
		{
			__m256i a = _mm256_max_epu8(
				_mm256_subs_epu8(ecur[0], rc[0]), _mm256_subs_epu8(ecur[1], rc[1]));
			for(size_t j = 2; j < ITER; j++) {
				a = _mm256_max_epu8(a, _mm256_subs_epu8(ecur[j], rc[j]));
			}
			for(size_t j = 0; j + 1 < ITER; j++) {
				a = _mm256_max_epu8(a, _mm256_max_epu8(
					_mm256_subs_epu8(hl[j], ec[j]), _mm256_subs_epu8(hl[j], pc[rn][j+1])));
			}
			const __m256i hlsh = _mm256_alignr_epi8(hl[ITER-1],
				_mm256_permute2x128_si256(hl[ITER-1], vtop, 0x02), 15);
			a = _mm256_max_epu8(a, _mm256_max_epu8(
				_mm256_subs_epu8(hl[ITER-1], ec[ITER-1]), _mm256_subs_epu8(hlsh, pc[rn][0])));
			co = _mm256_max_epu8(_mm256_max_epu8(a, _mm256_subs_epu8(y15, zc[rn])), cc);
		}
		// F, H and the next column's E with the carry from the lanes above
		__m256i hcur[ITER];
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
			hcur[j] = h;
			ecur[j] = e;
		}
		const uint8_t lr = reinterpret_cast<const uint8_t*>(pvCol + H)[lastWordIdx];
		if(lr > lrmax) {
			lrmax = lr;
		}
		pvCol += colstride;
		// H0 of the next column.  Its diagonal is row j - 1 of this column,
		// or for j = 0 the last vector row shifted down a lane, with 0xff in
		// lane 0.
		const __m256i *pvS = prof + rn * ITER * 2;
		{
			const __m256i p = _mm256_permute2x128_si256(hcur[ITER-1], vtop, 0x02);
			const __m256i d0 = _mm256_alignr_epi8(hcur[ITER-1], p, 15);
			h0[0] = _mm256_max_epu8(_mm256_subs_epu8(d0, _mm256_load_si256(pvS)), ecur[0]);
		}
		for(size_t j = 1; j < ITER; j++) {
			h0[j] = _mm256_max_epu8(
				_mm256_subs_epu8(hcur[j-1], _mm256_load_si256(pvS + 2 * j)), ecur[j]);
		}
	}
	return lrmax;
}

/**
 * How many columns before its first exact column alignEE8V256ShortSplit
 * starts its second chain.  A fill started at column s computes every cell
 * as the maximum over the paths that begin in row 0 at column s or later; the
 * true fill also has paths from earlier columns.  Such a path reaches column
 * s + W only after advancing W + 1 columns or more while descending at most
 * 32 * ITER - 1 rows (padding lanes included), so it takes at least
 * W + 2 - 32 * ITER read-gap steps, each costing rdgape or more.  With
 * W = 32 * ITER - 2 + ceil(255 / rdgape) its score has saturated to 0, which
 * no cell is below, so from column s + W on the two fills are equal.  (With
 * free diagonals and rdgapo = rdgape the bound is reached: one column fewer
 * is not exact.)
 */
template<size_t ITER>
static inline size_t alignEE8V256ShortWarmup(const Scoring& sc) {
	const size_t rdgape = (size_t)sc.readGapExtend();
	assert_gt(rdgape, 0);
	return 32 * ITER - 2 + (255 + rdgape - 1) / rdgape;
}

/**
 * alignEE8V256Short for a fill of many columns (e.g. a mate-rescue window),
 * run as two chains of columns at once: chain A over the first n columns
 * from rfi, and chain B over the last n columns to rff, started as a fresh
 * fill W = alignEE8V256ShortWarmup() columns before the first column A does
 * not cover, so it is exact from there on (see alignEE8V256ShortWarmup).
 * Where the two overlap, A writes each column at a later step than B (A
 * started earlier), so the matrix ends with A's values there; B's values in
 * the warm-up are at most the true ones, so they do not raise the final-row
 * maximum.
 *
 * Each step computes alignEE8V256Short's column, with one difference: the
 * next column's carry-out, the scan's input, is taken directly as
 * max_j(H0'[j] - cpen[j]) from the next column's H0, rather than as the
 * maximum of early terms (in HL and E) and late terms (in the scan's
 * values).  That puts H, E and H0 on the chain from one column to the next,
 * but needs about 25 fewer instructions per column; with two chains in
 * flight the fill is bound by issue rather than by one chain's latency, and
 * the alignEE8V256Short form of two chains was no faster than one chain.
 * Requires rff - rfi > W + 1.
 */
template<size_t ITER>
static uint8_t alignEE8V256ShortSplit(
	const __m256i *prof,    // query profile: per reference character, ITER (score, barrier) pairs
	__m256i *mat,           // matrix: column 0, row 0, E (SSEMatrix layout)
	size_t colstride,       // vectors per matrix column
	const char *rf,         // reference
	size_t rfi,             // first column
	size_t rff,             // last column (excl.)
	const Scoring& sc,      // gap penalties
	size_t lastWordIdx,     // final row: byte offset from a column's first H vector
	size_t W)               // alignEE8V256ShortWarmup<ITER>(sc)
{
	assert_gt(rff - rfi, W + 1);
	const size_t E = SSEMatrix::E, F = SSEMatrix::F, H = SSEMatrix::H;
	const __m256i zero = _mm256_setzero_si256();
	const __m256i vrfgapo = _mm256_set1_epi8((char)sc.refGapOpen());
	const __m256i vrfgape = _mm256_set1_epi8((char)sc.refGapExtend());
	const __m256i vrdgapo = _mm256_set1_epi8((char)sc.readGapOpen());
	const __m256i vrdgape = _mm256_set1_epi8((char)sc.readGapExtend());
	// The costs of alignEE8V256Short, per vector row j (bpe, cpen, vpen, epen
	// as there; lpen: H0 at row j to the carry into the lane's last row)
	__m256i bar[ITER], bpe[ITER], cpen[ITER], vpen[ITER], epen[ITER];
	__m256i vbar = zero;
	for(size_t j = 0; j < ITER; j++) {
		bar[j] = _mm256_load_si256(prof + 2 * j + 1);
		bpe[j] = _mm256_adds_epu8(bar[j], vrfgape);
		epen[j] = _mm256_adds_epu8(vrdgapo, bar[j]);
		vpen[j] = j == 0 ? bar[0] : _mm256_adds_epu8(_mm256_adds_epu8(vpen[j-1], vrfgape), bar[j]);
		vbar = _mm256_or_si256(vbar, bar[j]);
	}
	cpen[ITER-1] = vrfgapo;
	for(size_t j = ITER-1; j > 0; j--) {
		cpen[j-1] = _mm256_adds_epu8(cpen[j], bpe[j]);
	}
	__m256i lpen[ITER-1];
	lpen[ITER-2] = vrfgapo;
	for(size_t j = ITER-2; j > 0; j--) {
		lpen[j-1] = _mm256_adds_epu8(lpen[j], bpe[j]);
	}
	// The scan's costs (q1-q12, chi) and constants, as in alignEE8V256Short
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
	const __m256i vidx15 = _mm256_set_epi64x(
		0x0f0f0f0f0f0f0f0fLL, 0x0f0f0f0f0f0f0f0fLL,
		(long long)0x8080808080808080ULL, (long long)0x8080808080808080ULL);
	const __m256i vtop = _mm256_set_epi64x(0, 0, (long long)0xff00000000000000ULL, 0);
	const size_t n = (rff - rfi + W + 1) / 2; // columns per chain
	const size_t sb = rff - n;                // B's first column
	// One chain's state between columns: H0 and E of its next column, and
	// that column's carry-out
	__m256i h0a[ITER], ea[ITER], coa, h0b[ITER], eb[ITER], cob;
	// A chain's first column, from a staged initial column (zero), whose
	// diagonal is 0xff in row 0 only
	auto start = [&](__m256i *h0, __m256i *ecur, __m256i& co, size_t col) __attribute__((always_inline)) {
		const __m256i *pvS = prof + (size_t)firsts5[(int)rf[col]] * ITER * 2;
		const __m256i p = _mm256_permute2x128_si256(zero, vtop, 0x02);
		const __m256i d0 = _mm256_alignr_epi8(zero, p, 15);
		h0[0] = _mm256_subs_epu8(d0, _mm256_load_si256(pvS));
		co = _mm256_subs_epu8(h0[0], cpen[0]);
		for(size_t j = 1; j < ITER; j++) {
			h0[j] = zero;
		}
		for(size_t j = 0; j < ITER; j++) {
			ecur[j] = zero;
		}
	};
	// One column of a chain: stores its F, H and the next column's E at
	// pvCol, advances the state to the next column (reference character rn)
	// and returns the column's final-row H
	auto step = [&](__m256i *h0, __m256i *ecur, __m256i& co, size_t rn,
		__m256i *pvCol) __attribute__((always_inline)) -> uint8_t
	{
		// The scan within each half, then the carry into each lane (as
		// alignEE8V256Short)
		__m256i y = co;
		y = _mm256_max_epu8(
			_mm256_max_epu8(y, _mm256_subs_epu8(_mm256_slli_si256(y, 1), q1)),
			_mm256_max_epu8(_mm256_subs_epu8(_mm256_slli_si256(y, 2), q2),
			                _mm256_subs_epu8(_mm256_slli_si256(y, 3), q3)));
		y = _mm256_max_epu8(
			_mm256_max_epu8(y, _mm256_subs_epu8(_mm256_slli_si256(y, 4), q4)),
			_mm256_max_epu8(_mm256_subs_epu8(_mm256_slli_si256(y, 8), q8),
			                _mm256_subs_epu8(_mm256_slli_si256(y, 12), q12)));
		const __m256i py = _mm256_permute2x128_si256(y, y, 0x08);
		const __m256i y15 = _mm256_alignr_epi8(y, py, 15);
		const __m256i yc = _mm256_shuffle_epi8(py, vidx15);
		const __m256i c = _mm256_max_epu8(y15, _mm256_subs_epu8(yc, chi));
		// F and H with the carries from within the lane only
		__m256i fl[ITER], hl[ITER];
		fl[0] = zero;
		hl[0] = h0[0];
		{
			__m256i vf = _mm256_subs_epu8(h0[0], vrfgapo);
			for(size_t j = 1; j < ITER; j++) {
				fl[j] = _mm256_subs_epu8(vf, bar[j]);
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
		// F, H and the next column's E with the carry from the lanes above
		__m256i hcur[ITER];
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
			hcur[j] = h;
			ecur[j] = e;
		}
		const uint8_t lr = reinterpret_cast<const uint8_t*>(pvCol + H)[lastWordIdx];
		// H0 of the next column, and its carry-out
		const __m256i *pvS = prof + rn * ITER * 2;
		{
			const __m256i p = _mm256_permute2x128_si256(hcur[ITER-1], vtop, 0x02);
			const __m256i d0 = _mm256_alignr_epi8(hcur[ITER-1], p, 15);
			h0[0] = _mm256_max_epu8(_mm256_subs_epu8(d0, _mm256_load_si256(pvS)), ecur[0]);
		}
		for(size_t j = 1; j < ITER; j++) {
			h0[j] = _mm256_max_epu8(
				_mm256_subs_epu8(hcur[j-1], _mm256_load_si256(pvS + 2 * j)), ecur[j]);
		}
		__m256i m = _mm256_max_epu8(
			_mm256_subs_epu8(h0[0], cpen[0]), _mm256_subs_epu8(h0[1], cpen[1]));
		for(size_t j = 2; j + 1 < ITER; j += 2) {
			m = _mm256_max_epu8(m, _mm256_max_epu8(
				_mm256_subs_epu8(h0[j], cpen[j]), _mm256_subs_epu8(h0[j+1], cpen[j+1])));
		}
		if(ITER % 2 == 1) {
			m = _mm256_max_epu8(m, _mm256_subs_epu8(h0[ITER-1], cpen[ITER-1]));
		}
		co = m;
		return lr;
	};
	start(h0a, ea, coa, rfi);
	start(h0b, eb, cob, sb);
	uint8_t lrmax = 0;
	__m256i *pa = mat, *pb = mat + colstride * (sb - rfi);
	for(size_t k = 0; k < n; k++) {
		// The next columns' reference characters.  A's next column is at
		// most rfi + n < rff; after B's last column, the values computed with
		// its character are not used.
		const size_t ra = (size_t)firsts5[(int)rf[rfi + k + 1]];
		const size_t ib = sb + k + 1 < rff ? sb + k + 1 : sb + k;
		const size_t rb = (size_t)firsts5[(int)rf[ib]];
		const uint8_t la = step(h0a, ea, coa, ra, pa);
		const uint8_t lb = step(h0b, eb, cob, rb, pb);
		const uint8_t l = la > lb ? la : lb;
		lrmax = l > lrmax ? l : lrmax;
		pa += colstride;
		pb += colstride;
	}
	return lrmax;
}

/**
 * The fill of alignNucleotidesEnd2EndSseU8 for a column of one 256-bit vector
 * (up to 32 rows), whose lanes are the rows in order.  The column is held as
 * two 128-bit halves (rows 0-15 and 16-31), so that a shift across them is a
 * single in-lane instruction.  It computes the general loop's saturating u8
 * recurrence and stores the same E, F and H vectors in the same places:
 *  - H0 = max(diag - profile, E), with F = 0 from the first pass, as every
 *    row is the first of its lane;
 *  - F is the general loop's prefix scan of the carries H0 - rfgapo, run
 *    within each half, four rows per step, and then carried from row 15 into
 *    the upper half; rfgapo is subtracted after the scan, with the barrier;
 *  - H = max(H0, F), and the next column's E = max(E - rdgape, H - rdgapo -
 *    barrier).
 * The speed is set by the chain of dependent instructions from one column to
 * the next.  The next column's H0, the scan's input, is max(diag' -
 * profile', E') with diag' = H shifted down a row; as subtraction
 * distributes over max, it is also the max of terms in this column's H0 and
 * E, which are ready while the scan runs, and of terms in the scan's values
 * shifted down one and two rows, each with its whole cost.  The chain is
 * then the scan, one shift and one subtraction.  sat(sat(x - a) - b) =
 * sat(x - sat(a + b)), so costs may be summed before they are subtracted.
 * Returns the maximum of H over the final row.
 */
static uint8_t alignEE8V256Row32(
	const __m128i *pvProf,  // query profile: per reference character, the score and barrier halves
	__m128i *pvCol,         // matrix: column 0, row 0, E (SSEMatrix layout, in halves)
	size_t colstride,       // 256-bit vectors per matrix column
	const char *rf,         // reference
	size_t rfi,             // first column
	size_t rff,             // last column (excl.)
	const Scoring& sc,      // gap penalties
	size_t lastWordIdx)     // final row
{
	const size_t E = SSEMatrix::E, F = SSEMatrix::F, H = SSEMatrix::H;
	// Gap barrier rows (the same for every reference character)
	const __m128i barlo = _mm_load_si128(pvProf + 2);
	const __m128i barhi = _mm_load_si128(pvProf + 3);
	// Cost of a carry passing a row: rfgape, or all of it in the barrier.
	// passN: the cost of passing the N rows that end at each row, within a
	// half.
	const __m128i vrfgape = _mm_set1_epi8((char)sc.refGapExtend());
	const __m128i pass1lo = _mm_or_si128(vrfgape, barlo);
	const __m128i pass1hi = _mm_or_si128(vrfgape, barhi);
	const __m128i pass2lo = _mm_adds_epu8(pass1lo, _mm_slli_si128(pass1lo, 1));
	const __m128i pass2hi = _mm_adds_epu8(pass1hi, _mm_slli_si128(pass1hi, 1));
	const __m128i pass3lo = _mm_adds_epu8(pass2lo, _mm_slli_si128(pass1lo, 2));
	const __m128i pass3hi = _mm_adds_epu8(pass2hi, _mm_slli_si128(pass1hi, 2));
	const __m128i pass4lo = _mm_adds_epu8(pass2lo, _mm_slli_si128(pass2lo, 2));
	const __m128i pass4hi = _mm_adds_epu8(pass2hi, _mm_slli_si128(pass2hi, 2));
	const __m128i pass8lo = _mm_adds_epu8(pass4lo, _mm_slli_si128(pass4lo, 4));
	const __m128i pass8hi = _mm_adds_epu8(pass4hi, _mm_slli_si128(pass4hi, 4));
	const __m128i pass12lo = _mm_adds_epu8(pass8lo, _mm_slli_si128(pass4lo, 8));
	const __m128i pass12hi = _mm_adds_epu8(pass8hi, _mm_slli_si128(pass4hi, 8));
	const __m128i pass16hi = _mm_adds_epu8(pass8hi, _mm_slli_si128(pass8hi, 8));
	// rfgapo or rdgapo, plus the gap barrier
	const __m128i vrfgapo = _mm_set1_epi8((char)sc.refGapOpen());
	const __m128i vrdgapo = _mm_set1_epi8((char)sc.readGapOpen());
	const __m128i fpenlo = _mm_adds_epu8(vrfgapo, barlo);
	const __m128i fpenhi = _mm_adds_epu8(vrfgapo, barhi);
	const __m128i epenlo = _mm_adds_epu8(vrdgapo, barlo);
	const __m128i epenhi = _mm_adds_epu8(vrdgapo, barhi);
	// The carry from row 15 to row 16 + k, k >= 1, loses the cost of
	// passing rows 16 to 15 + k, then rfgapo and row 16 + k's barrier; row
	// 16 itself takes row 15's value through the shift.  carryhiSh: the same
	// for F shifted down a row.
	const __m128i carryhi = _mm_adds_epu8(
		_mm_or_si128(_mm_slli_si128(pass16hi, 1), _mm_cvtsi32_si128(0xff)), fpenhi);
	const __m128i carryhiSh = _mm_or_si128(_mm_slli_si128(carryhi, 1), _mm_cvtsi32_si128(0xff));
	const __m128i vrdgape = _mm_set1_epi8((char)sc.readGapExtend());
	// Byte 15 is 0xff: shifted into row 0, it is row 0's diagonal
	const __m128i vtop = _mm_slli_si128(_mm_cvtsi32_si128(0xff), 15);
	const __m128i vrow15 = _mm_set1_epi8(15);
	// Costs to the next column's H0 from the scan's values:
	//  k1: one row down (F, then rdgapo and the barrier into E');
	//  k2[r]: two rows down (F one row up, then the diagonal into H0 for
	//  reference character r);
	//  kbc[r]: from row 15 carried into the upper half, for both.
	const __m128i k1lo = _mm_adds_epu8(fpenlo, epenlo);
	const __m128i k1hi = _mm_adds_epu8(fpenhi, epenhi);
	const __m128i fpenShlo = _mm_slli_si128(fpenlo, 1);
	const __m128i fpenShhi = _mm_alignr_epi8(fpenhi, fpenlo, 15);
	__m128i k2lo[5], k2hi[5], kbc[5];
	for(size_t r = 0; r < 5; r++) {
		k2lo[r] = _mm_adds_epu8(fpenShlo, _mm_load_si128(pvProf + 4 * r));
		k2hi[r] = _mm_adds_epu8(fpenShhi, _mm_load_si128(pvProf + 4 * r + 1));
		kbc[r] = _mm_min_epu8(_mm_adds_epu8(carryhi, epenhi),
			_mm_adds_epu8(carryhiSh, _mm_load_si128(pvProf + 4 * r + 1)));
	}
	// H0 and E of the column: for the first, from the staged initial column
	// (zero)
	__m128i h0lo, h0hi;
	{
		const size_t r = (size_t)firsts5[(int)rf[rfi]];
		h0lo = _mm_subs_epu8(_mm_cvtsi32_si128(0xff), _mm_load_si128(pvProf + 4 * r));
		h0hi = _mm_setzero_si128();
	}
	__m128i elo = _mm_setzero_si128(), ehi = _mm_setzero_si128();
	// Maximum of H over the columns, for the final row's maximum
	__m128i hmaxlo = _mm_setzero_si128(), hmaxhi = _mm_setzero_si128();
	for(size_t i = rfi; i < rff; i++) {
		// The next column's reference character.  (After the last column,
		// the values computed with it are not used.)
		const size_t rn = (size_t)firsts5[(int)rf[i + 1 < rff ? i + 1 : i]];
		// The scan within each half, four rows per step: over the 4 rows
		// ending at each row, then over 16
		__m128i ylo = h0lo, yhi = h0hi;
		ylo = _mm_max_epu8(
			_mm_max_epu8(ylo, _mm_subs_epu8(_mm_slli_si128(ylo, 1), pass1lo)),
			_mm_max_epu8(_mm_subs_epu8(_mm_slli_si128(ylo, 2), pass2lo),
			             _mm_subs_epu8(_mm_slli_si128(ylo, 3), pass3lo)));
		yhi = _mm_max_epu8(
			_mm_max_epu8(yhi, _mm_subs_epu8(_mm_slli_si128(yhi, 1), pass1hi)),
			_mm_max_epu8(_mm_subs_epu8(_mm_slli_si128(yhi, 2), pass2hi),
			             _mm_subs_epu8(_mm_slli_si128(yhi, 3), pass3hi)));
		ylo = _mm_max_epu8(
			_mm_max_epu8(ylo, _mm_subs_epu8(_mm_slli_si128(ylo, 4), pass4lo)),
			_mm_max_epu8(_mm_subs_epu8(_mm_slli_si128(ylo, 8), pass8lo),
			             _mm_subs_epu8(_mm_slli_si128(ylo, 12), pass12lo)));
		yhi = _mm_max_epu8(
			_mm_max_epu8(yhi, _mm_subs_epu8(_mm_slli_si128(yhi, 4), pass4hi)),
			_mm_max_epu8(_mm_subs_epu8(_mm_slli_si128(yhi, 8), pass8hi),
			             _mm_subs_epu8(_mm_slli_si128(yhi, 12), pass12hi)));
		// The scan's values one row down (the carry into each row), two rows
		// down, and row 15's, for the upper half
		const __m128i y1lo = _mm_slli_si128(ylo, 1);
		const __m128i y2lo = _mm_slli_si128(ylo, 2);
		const __m128i y1hi = _mm_alignr_epi8(yhi, ylo, 15);
		const __m128i y2hi = _mm_alignr_epi8(yhi, ylo, 14);
		const __m128i y15 = _mm_shuffle_epi8(ylo, vrow15);
		// The next column's H0: its terms in this column's H0 and E ...
		const __m128i glo = _mm_max_epu8(
			_mm_max_epu8(
				_mm_subs_epu8(_mm_alignr_epi8(h0lo, vtop, 15), _mm_load_si128(pvProf + 4 * rn)),
				_mm_subs_epu8(elo, vrdgape)),
			_mm_subs_epu8(h0lo, epenlo));
		const __m128i ghi = _mm_max_epu8(
			_mm_max_epu8(
				_mm_subs_epu8(_mm_alignr_epi8(h0hi, h0lo, 15), _mm_load_si128(pvProf + 4 * rn + 1)),
				_mm_subs_epu8(ehi, vrdgape)),
			_mm_subs_epu8(h0hi, epenhi));
		// ... and in F, from the scan
		const __m128i nlo = _mm_max_epu8(
			_mm_max_epu8(glo, _mm_subs_epu8(y1lo, k1lo)), _mm_subs_epu8(y2lo, k2lo[rn]));
		const __m128i nhi = _mm_max_epu8(
			_mm_max_epu8(ghi, _mm_subs_epu8(y1hi, k1hi)),
			_mm_max_epu8(_mm_subs_epu8(y2hi, k2hi[rn]), _mm_subs_epu8(y15, kbc[rn])));
		// This column's F, H and the next column's E
		const __m128i flo = _mm_subs_epu8(y1lo, fpenlo);
		const __m128i fhi = _mm_max_epu8(
			_mm_subs_epu8(y1hi, fpenhi), _mm_subs_epu8(y15, carryhi));
		const __m128i hlo = _mm_max_epu8(h0lo, flo);
		const __m128i hhi = _mm_max_epu8(h0hi, fhi);
		elo = _mm_max_epu8(_mm_subs_epu8(elo, vrdgape), _mm_subs_epu8(hlo, epenlo));
		ehi = _mm_max_epu8(_mm_subs_epu8(ehi, vrdgape), _mm_subs_epu8(hhi, epenhi));
		_mm_store_si128(pvCol + 2 * F, flo);
		_mm_store_si128(pvCol + 2 * F + 1, fhi);
		_mm_store_si128(pvCol + 2 * H, hlo);
		_mm_store_si128(pvCol + 2 * H + 1, hhi);
		pvCol += 2 * colstride;
		_mm_store_si128(pvCol + 2 * E, elo);
		_mm_store_si128(pvCol + 2 * E + 1, ehi);
		hmaxlo = _mm_max_epu8(hmaxlo, hlo);
		hmaxhi = _mm_max_epu8(hmaxhi, hhi);
		h0lo = nlo;
		h0hi = nhi;
	}
	uint8_t hmax[32];
	_mm_storeu_si128(reinterpret_cast<__m128i*>(hmax), hmaxlo);
	_mm_storeu_si128(reinterpret_cast<__m128i*>(hmax + 16), hmaxhi);
	return hmax[lastWordIdx];
}

/**
 * How many columns before its first exact column alignEE8V256Row32Split
 * starts its second chain.  A fill started at column s computes every cell
 * as the maximum over the paths that begin in row 0 at column s or later; the
 * true fill also has paths from earlier columns.  Such a path reaches column
 * s + W only after advancing W + 1 columns or more while descending at most
 * 31 rows, so it takes at least W - 30 read-gap steps, each costing rdgape or
 * more.  With W - 30 >= 255 / rdgape its score has saturated to 0, which no
 * cell is below, so from column s + W on the two fills are equal.
 */
static inline size_t alignEE8V256Row32Warmup(const Scoring& sc) {
	const size_t rdgape = (size_t)sc.readGapExtend();
	assert_gt(rdgape, 0);
	return 32 + (255 + rdgape - 1) / rdgape;
}

/**
 * alignEE8V256Row32 for a fill of many columns.  alignEE8V256Row32 is bound
 * by the chain of dependent instructions from one column to the next, and
 * uses only one 128-bit lane of the vector units.  Here the fill runs as two
 * chains, one in each 128-bit lane of 256-bit registers, doing the same
 * column arithmetic: chain A in lane 0 over the first n columns from rfi,
 * and chain B in lane 1 over the last n columns to rff.  B starts as a fresh
 * fill W = alignEE8V256Row32Warmup() columns before the first column A does
 * not cover, so it is exact from there on (see alignEE8V256Row32Warmup).
 * Where the two overlap, A writes each column at a later step than B (A
 * started earlier), so the matrix ends with A's values there; B's values in
 * the warm-up are at most the true ones, so they do not raise the final-row
 * maximum.  Register "lo" holds rows 0-15 of both chains and "hi" rows
 * 16-31, so every shift of alignEE8V256Row32's column stays within a lane.
 * The per-column costs that depend on the reference character are looked up
 * per pair (A's character, B's).  In each iteration the next column's scan
 * comes before this column's F, H, E and stores, so that it is older than
 * them when both are ready and goes first.  Requires rff - rfi > W + 1.
 */
static uint8_t alignEE8V256Row32Split(
	const __m128i *pvProf,  // query profile: per reference character, the score and barrier halves
	__m128i *pvCol,         // matrix: column 0, row 0, E (SSEMatrix layout, in halves)
	size_t colstride,       // 256-bit vectors per matrix column
	const char *rf,         // reference
	size_t rfi,             // first column
	size_t rff,             // last column (excl.)
	const Scoring& sc,      // gap penalties
	size_t lastWordIdx,     // final row
	size_t W)               // alignEE8V256Row32Warmup(sc)
{
	assert_gt(rff - rfi, W + 1);
	const size_t E = SSEMatrix::E, F = SSEMatrix::F, H = SSEMatrix::H;
	// The costs of alignEE8V256Row32, for one chain (128 bits)
	const __m128i barlo = _mm_load_si128(pvProf + 2);
	const __m128i barhi = _mm_load_si128(pvProf + 3);
	const __m128i vrfgape = _mm_set1_epi8((char)sc.refGapExtend());
	const __m128i pass1lo = _mm_or_si128(vrfgape, barlo);
	const __m128i pass1hi = _mm_or_si128(vrfgape, barhi);
	const __m128i pass2lo = _mm_adds_epu8(pass1lo, _mm_slli_si128(pass1lo, 1));
	const __m128i pass2hi = _mm_adds_epu8(pass1hi, _mm_slli_si128(pass1hi, 1));
	const __m128i pass3lo = _mm_adds_epu8(pass2lo, _mm_slli_si128(pass1lo, 2));
	const __m128i pass3hi = _mm_adds_epu8(pass2hi, _mm_slli_si128(pass1hi, 2));
	const __m128i pass4lo = _mm_adds_epu8(pass2lo, _mm_slli_si128(pass2lo, 2));
	const __m128i pass4hi = _mm_adds_epu8(pass2hi, _mm_slli_si128(pass2hi, 2));
	const __m128i pass8lo = _mm_adds_epu8(pass4lo, _mm_slli_si128(pass4lo, 4));
	const __m128i pass8hi = _mm_adds_epu8(pass4hi, _mm_slli_si128(pass4hi, 4));
	const __m128i pass12lo = _mm_adds_epu8(pass8lo, _mm_slli_si128(pass4lo, 8));
	const __m128i pass12hi = _mm_adds_epu8(pass8hi, _mm_slli_si128(pass4hi, 8));
	const __m128i pass16hi = _mm_adds_epu8(pass8hi, _mm_slli_si128(pass8hi, 8));
	const __m128i vrfgapo = _mm_set1_epi8((char)sc.refGapOpen());
	const __m128i vrdgapo = _mm_set1_epi8((char)sc.readGapOpen());
	const __m128i fpenlo = _mm_adds_epu8(vrfgapo, barlo);
	const __m128i fpenhi = _mm_adds_epu8(vrfgapo, barhi);
	const __m128i epenlo = _mm_adds_epu8(vrdgapo, barlo);
	const __m128i epenhi = _mm_adds_epu8(vrdgapo, barhi);
	const __m128i carryhi = _mm_adds_epu8(
		_mm_or_si128(_mm_slli_si128(pass16hi, 1), _mm_cvtsi32_si128(0xff)), fpenhi);
	const __m128i carryhiSh = _mm_or_si128(_mm_slli_si128(carryhi, 1), _mm_cvtsi32_si128(0xff));
	const __m128i k1lo = _mm_adds_epu8(fpenlo, epenlo);
	const __m128i k1hi = _mm_adds_epu8(fpenhi, epenhi);
	const __m128i fpenShlo = _mm_slli_si128(fpenlo, 1);
	const __m128i fpenShhi = _mm_alignr_epi8(fpenhi, fpenlo, 15);
	__m128i k2lo[5], k2hi[5], kbc[5];
	for(size_t r = 0; r < 5; r++) {
		k2lo[r] = _mm_adds_epu8(fpenShlo, _mm_load_si128(pvProf + 4 * r));
		k2hi[r] = _mm_adds_epu8(fpenShhi, _mm_load_si128(pvProf + 4 * r + 1));
		kbc[r] = _mm_min_epu8(_mm_adds_epu8(carryhi, epenhi),
			_mm_adds_epu8(carryhiSh, _mm_load_si128(pvProf + 4 * r + 1)));
	}
	// The same costs in both lanes
	const __m256i wpass1lo = _mm256_broadcastsi128_si256(pass1lo);
	const __m256i wpass1hi = _mm256_broadcastsi128_si256(pass1hi);
	const __m256i wpass2lo = _mm256_broadcastsi128_si256(pass2lo);
	const __m256i wpass2hi = _mm256_broadcastsi128_si256(pass2hi);
	const __m256i wpass3lo = _mm256_broadcastsi128_si256(pass3lo);
	const __m256i wpass3hi = _mm256_broadcastsi128_si256(pass3hi);
	const __m256i wpass4lo = _mm256_broadcastsi128_si256(pass4lo);
	const __m256i wpass4hi = _mm256_broadcastsi128_si256(pass4hi);
	const __m256i wpass8lo = _mm256_broadcastsi128_si256(pass8lo);
	const __m256i wpass8hi = _mm256_broadcastsi128_si256(pass8hi);
	const __m256i wpass12lo = _mm256_broadcastsi128_si256(pass12lo);
	const __m256i wpass12hi = _mm256_broadcastsi128_si256(pass12hi);
	const __m256i wfpenlo = _mm256_broadcastsi128_si256(fpenlo);
	const __m256i wfpenhi = _mm256_broadcastsi128_si256(fpenhi);
	const __m256i wepenlo = _mm256_broadcastsi128_si256(epenlo);
	const __m256i wepenhi = _mm256_broadcastsi128_si256(epenhi);
	const __m256i wcarryhi = _mm256_broadcastsi128_si256(carryhi);
	const __m256i wk1lo = _mm256_broadcastsi128_si256(k1lo);
	const __m256i wk1hi = _mm256_broadcastsi128_si256(k1hi);
	const __m256i wrdgape = _mm256_set1_epi8((char)sc.readGapExtend());
	// Byte 15 of each lane is 0xff: shifted into row 0, it is row 0's diagonal
	const __m256i wtop = _mm256_broadcastsi128_si256(_mm_slli_si128(_mm_cvtsi32_si128(0xff), 15));
	const __m256i wrow15 = _mm256_set1_epi8(15);
	// Per pair of reference characters (A's, B's): the profile's two halves,
	// k2lo, k2hi and kbc, A's in lane 0 and B's in lane 1
	__m256i pair[25][5];
	for(size_t a = 0; a < 5; a++) {
		for(size_t b = 0; b < 5; b++) {
			const __m128i va[5] = { _mm_load_si128(pvProf + 4 * a), _mm_load_si128(pvProf + 4 * a + 1),
				k2lo[a], k2hi[a], kbc[a] };
			const __m128i vb[5] = { _mm_load_si128(pvProf + 4 * b), _mm_load_si128(pvProf + 4 * b + 1),
				k2lo[b], k2hi[b], kbc[b] };
			for(size_t m = 0; m < 5; m++) {
				pair[a * 5 + b][m] = _mm256_inserti128_si256(_mm256_castsi128_si256(va[m]), vb[m], 1);
			}
		}
	}
	const size_t n = (rff - rfi + W + 1) / 2; // columns per chain
	const size_t sb = rff - n;                // B's first column
	// alignEE8V256Row32's scan (F within each half), for both chains
	auto scan = [&](__m256i ylo, __m256i yhi, __m256i& outlo, __m256i& outhi) {
		ylo = _mm256_max_epu8(
			_mm256_max_epu8(ylo, _mm256_subs_epu8(_mm256_slli_si256(ylo, 1), wpass1lo)),
			_mm256_max_epu8(_mm256_subs_epu8(_mm256_slli_si256(ylo, 2), wpass2lo),
			                _mm256_subs_epu8(_mm256_slli_si256(ylo, 3), wpass3lo)));
		yhi = _mm256_max_epu8(
			_mm256_max_epu8(yhi, _mm256_subs_epu8(_mm256_slli_si256(yhi, 1), wpass1hi)),
			_mm256_max_epu8(_mm256_subs_epu8(_mm256_slli_si256(yhi, 2), wpass2hi),
			                _mm256_subs_epu8(_mm256_slli_si256(yhi, 3), wpass3hi)));
		outlo = _mm256_max_epu8(
			_mm256_max_epu8(ylo, _mm256_subs_epu8(_mm256_slli_si256(ylo, 4), wpass4lo)),
			_mm256_max_epu8(_mm256_subs_epu8(_mm256_slli_si256(ylo, 8), wpass8lo),
			                _mm256_subs_epu8(_mm256_slli_si256(ylo, 12), wpass12lo)));
		outhi = _mm256_max_epu8(
			_mm256_max_epu8(yhi, _mm256_subs_epu8(_mm256_slli_si256(yhi, 4), wpass4hi)),
			_mm256_max_epu8(_mm256_subs_epu8(_mm256_slli_si256(yhi, 8), wpass8hi),
			                _mm256_subs_epu8(_mm256_slli_si256(yhi, 12), wpass12hi)));
	};
	// H0 and E of each chain's first column, from the staged initial column
	__m256i h0lo, h0hi;
	{
		const size_t ra = (size_t)firsts5[(int)rf[rfi]];
		const size_t rb = (size_t)firsts5[(int)rf[sb]];
		const __m256i first = _mm256_broadcastsi128_si256(_mm_cvtsi32_si128(0xff));
		h0lo = _mm256_subs_epu8(first, pair[ra * 5 + rb][0]);
		h0hi = _mm256_setzero_si256();
	}
	__m256i elo = _mm256_setzero_si256(), ehi = _mm256_setzero_si256();
	__m256i hmaxlo = _mm256_setzero_si256(), hmaxhi = _mm256_setzero_si256();
	__m128i *pa = pvCol, *pb = pvCol + 2 * colstride * (sb - rfi);
	__m256i ylo, yhi;
	scan(h0lo, h0hi, ylo, yhi);
	for(size_t k = 0; k < n; k++) {
		// The next columns' reference characters.  A's next column is at
		// most rfi + n < rff; after B's last column, the values computed with
		// its character are not used.
		const size_t ib = sb + k + 1 < rff ? sb + k + 1 : sb + k;
		const __m256i *t = pair[(size_t)firsts5[(int)rf[rfi + k + 1]] * 5 +
		                        (size_t)firsts5[(int)rf[ib]]];
		const __m256i y1lo = _mm256_slli_si256(ylo, 1);
		const __m256i y2lo = _mm256_slli_si256(ylo, 2);
		const __m256i y1hi = _mm256_alignr_epi8(yhi, ylo, 15);
		const __m256i y2hi = _mm256_alignr_epi8(yhi, ylo, 14);
		const __m256i y15 = _mm256_shuffle_epi8(ylo, wrow15);
		const __m256i glo = _mm256_max_epu8(
			_mm256_max_epu8(
				_mm256_subs_epu8(_mm256_alignr_epi8(h0lo, wtop, 15), t[0]),
				_mm256_subs_epu8(elo, wrdgape)),
			_mm256_subs_epu8(h0lo, wepenlo));
		const __m256i ghi = _mm256_max_epu8(
			_mm256_max_epu8(
				_mm256_subs_epu8(_mm256_alignr_epi8(h0hi, h0lo, 15), t[1]),
				_mm256_subs_epu8(ehi, wrdgape)),
			_mm256_subs_epu8(h0hi, wepenhi));
		const __m256i nlo = _mm256_max_epu8(
			_mm256_max_epu8(glo, _mm256_subs_epu8(y1lo, wk1lo)), _mm256_subs_epu8(y2lo, t[2]));
		const __m256i nhi = _mm256_max_epu8(
			_mm256_max_epu8(ghi, _mm256_subs_epu8(y1hi, wk1hi)),
			_mm256_max_epu8(_mm256_subs_epu8(y2hi, t[3]), _mm256_subs_epu8(y15, t[4])));
		// The next column's scan, ahead of this column's F, H and E
		__m256i ynlo, ynhi;
		scan(nlo, nhi, ynlo, ynhi);
		const __m256i flo = _mm256_subs_epu8(y1lo, wfpenlo);
		const __m256i fhi = _mm256_max_epu8(
			_mm256_subs_epu8(y1hi, wfpenhi), _mm256_subs_epu8(y15, wcarryhi));
		const __m256i hlo = _mm256_max_epu8(h0lo, flo);
		const __m256i hhi = _mm256_max_epu8(h0hi, fhi);
		elo = _mm256_max_epu8(_mm256_subs_epu8(elo, wrdgape), _mm256_subs_epu8(hlo, wepenlo));
		ehi = _mm256_max_epu8(_mm256_subs_epu8(ehi, wrdgape), _mm256_subs_epu8(hhi, wepenhi));
		// Lane 0 to A's column, lane 1 to B's
		_mm256_store_si256(reinterpret_cast<__m256i*>(pa + 2 * F), _mm256_permute2x128_si256(flo, fhi, 0x20));
		_mm256_store_si256(reinterpret_cast<__m256i*>(pb + 2 * F), _mm256_permute2x128_si256(flo, fhi, 0x31));
		_mm256_store_si256(reinterpret_cast<__m256i*>(pa + 2 * H), _mm256_permute2x128_si256(hlo, hhi, 0x20));
		_mm256_store_si256(reinterpret_cast<__m256i*>(pb + 2 * H), _mm256_permute2x128_si256(hlo, hhi, 0x31));
		pa += 2 * colstride;
		pb += 2 * colstride;
		_mm256_store_si256(reinterpret_cast<__m256i*>(pa + 2 * E), _mm256_permute2x128_si256(elo, ehi, 0x20));
		_mm256_store_si256(reinterpret_cast<__m256i*>(pb + 2 * E), _mm256_permute2x128_si256(elo, ehi, 0x31));
		hmaxlo = _mm256_max_epu8(hmaxlo, hlo);
		hmaxhi = _mm256_max_epu8(hmaxhi, hhi);
		h0lo = nlo;
		h0hi = nhi;
		ylo = ynlo;
		yhi = ynhi;
	}
	uint8_t hmax[32];
	_mm_storeu_si128(reinterpret_cast<__m128i*>(hmax), _mm_max_epu8(
		_mm256_castsi256_si128(hmaxlo), _mm256_extracti128_si256(hmaxlo, 1)));
	_mm_storeu_si128(reinterpret_cast<__m128i*>(hmax + 16), _mm_max_epu8(
		_mm256_castsi256_si128(hmaxhi), _mm256_extracti128_si256(hmaxhi, 1)));
	return hmax[lastWordIdx];
}
#endif

/**
 * Solve the current alignment problem using SSE instructions that operate on 16
 * unsigned 8-bit values packed into a single 128-bit register.
 */
TAlScore SwAligner::alignNucleotidesEnd2EndSseU8(int& flag, bool debug) {
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
	buildQueryProfileEnd2EndSseU8(fw_);
	assert(!d.profbuf_.empty());

	assert_eq(0, d.maxBonus_);
	const size_t iter =
		(dpRows() + (NWORDS_PER_REG-1)) / NWORDS_PER_REG; // iter = segLen

        const size_t lastWordIdx = NWORDS_PER_REG*(d.lastIter_*ROWSTRIDE)+d.lastWord_;

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
	SSERegI ve       = sse_setzero_siall();
	SSERegI vf       = sse_setzero_siall();
	SSERegI vh       = sse_setzero_siall();
#if 0
	SSERegI vhd      = sse_setzero_siall();
	SSERegI vhdtmp   = sse_setzero_siall();
#endif
	SSERegI vtmp     = sse_setzero_siall();
	SSERegI vhilsw   = sse_setzero_siall();

	SSERegI vs0     = sse_setzero_siall();
	SSERegI vs1     = sse_setzero_siall();

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
	
	// vhilsw: topmost (least sig) word set to 0xff, all other words=0
	sse_set_low_u8(0xff, vhilsw);	
	
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
	
	// Maximum score in final row
	bool found = false;
	TCScore lrmax = MIN_U8;
	
	for(size_t i = 0; i < iter; i++) {
		sse_store_siall(pvETmp, vlo);
		sse_store_siall(pvHTmp, vlo); // start high in end-to-end mode
		pvETmp += ROWSTRIDE;
		pvHTmp += ROWSTRIDE;
	}
	// These are swapped just before the innermost loop
	SSERegI *pvHStore = d.mat_.hvec(0, 0);
	SSERegI *pvHLoad  = d.mat_.tmpvec(0, 0);
	SSERegI *pvELoad  = d.mat_.evec(0, 0);
	SSERegI *pvEStore = d.mat_.evecUnsafe(0, 1);
	SSERegI *pvFStore = d.mat_.fvec(0, 0);

	assert_gt(sc_->gapbar, 0);
	size_t nfixup = 0;

	// F (reference gap) values flow down a column, and in the striped layout
	// the first pass over a column only carries them down within each lane.
	// Rather than Farrar's lazy-F loop, which carries them one lane further
	// per pass until nothing changes, the carries are propagated across all
	// lanes at once with a prefix scan (max-plus, saturating).  A carry passes
	// through a lane at the cost of iter gap extensions, unless a row in the
	// lane is inside the gap barrier, which sets F to 0 and so blocks it.
	// vpass1 holds that cost for the lane below each lane; vpass2, 4, ... are
	// its sums over 2, 4, ... consecutive lanes, for the scan's steps.  The
	// barrier vectors are the same for every reference character.
	SSERegI vpass1, vpass2, vpass4, vpass8;
#if NBYTES_PER_REG == 32
	SSERegI vpass16;
#endif
	{
		const SSERegI *pvBar = d.profbuf_.ptr() + 1;
		SSERegI vbar = sse_setzero_siall();
		for(size_t j = 0; j < iter; j++) {
			vbar = sse_or_siall(vbar, pvBar[2 * j]);
		}
		size_t passPen = iter * (size_t)sc_->refGapExtend();
		if(passPen > 0xff) {
			passPen = 0xff;
		}
		SSERegI vpass;
		sse_fill_u8((int)passPen, vpass);
		vpass  = sse_or_siall(vpass, vbar);
		vpass1 = sse_slli_siall(vpass, NBYTES_PER_WORD);
		vpass2 = sse_adds_epu8(vpass1, sse_slli_siall(vpass1, 1 * NBYTES_PER_WORD));
		vpass4 = sse_adds_epu8(vpass2, sse_slli_siall(vpass2, 2 * NBYTES_PER_WORD));
		vpass8 = sse_adds_epu8(vpass4, sse_slli_siall(vpass4, 4 * NBYTES_PER_WORD));
#if NBYTES_PER_REG == 32
		vpass16 = sse_adds_epu8(vpass8, sse_slli_siall(vpass8, 8 * NBYTES_PER_WORD));
#endif
	}

	// The final row of H vectors of the previous column, kept in a register
	// from one column to the next; before the first column, the staged
	// initial column (vlo).
	SSERegI vhlast = vlo;

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

	colstop_ = rff_ - 1;
	lastsolcol_ = 0;

#if NBYTES_PER_REG == 32
	if(iter <= 4) {
		// Up to 128 rows.  The loop below spends most of each column on the
		// chain of dependent instructions from one column to the next, much
		// of it byte shifts across the two 128-bit halves (vperm2i128 then
		// vpalignr); these fills compute the same cells with a shorter chain.
		if(iter == 1) {
			assert_eq(0, d.lastIter_);
			// A fill of many columns (e.g. a mate-rescue window) runs as two
			// chains at once, the second started this many columns early
			const size_t warmup = alignEE8V256Row32Warmup(*sc_);
			if((size_t)(rff_ - rfi_) >= warmup + 32) {
				lrmax = alignEE8V256Row32Split(
					reinterpret_cast<const __m128i*>(d.profbuf_.ptr()),
					reinterpret_cast<__m128i*>(d.mat_.evec(0, 0)), colstride, rf_,
					(size_t)rfi_, (size_t)rff_, *sc_, lastWordIdx, warmup);
			} else {
				lrmax = alignEE8V256Row32(
					reinterpret_cast<const __m128i*>(d.profbuf_.ptr()),
					reinterpret_cast<__m128i*>(d.mat_.evec(0, 0)), colstride, rf_,
					(size_t)rfi_, (size_t)rff_, *sc_, lastWordIdx);
			}
		} else {
			// A fill of many columns (e.g. a mate-rescue window) runs as two
			// chains at once, the second started W columns early.  Below about
			// 3 W columns (3.4 W at ITER = 4) the warm-up costs more than the
			// second chain gains.
			const size_t ncol = (size_t)(rff_ - rfi_);
			if(iter == 2) {
				const size_t W = alignEE8V256ShortWarmup<2>(*sc_);
				lrmax = ncol >= 4 * W ?
					alignEE8V256ShortSplit<2>(d.profbuf_.ptr(), d.mat_.ptr(), colstride,
						rf_, (size_t)rfi_, (size_t)rff_, *sc_, lastWordIdx, W) :
					alignEE8V256Short<2>(d.profbuf_.ptr(), d.mat_.ptr(), colstride,
						rf_, (size_t)rfi_, (size_t)rff_, *sc_, lastWordIdx);
			} else if(iter == 3) {
				const size_t W = alignEE8V256ShortWarmup<3>(*sc_);
				lrmax = ncol >= 4 * W ?
					alignEE8V256ShortSplit<3>(d.profbuf_.ptr(), d.mat_.ptr(), colstride,
						rf_, (size_t)rfi_, (size_t)rff_, *sc_, lastWordIdx, W) :
					alignEE8V256Short<3>(d.profbuf_.ptr(), d.mat_.ptr(), colstride,
						rf_, (size_t)rfi_, (size_t)rff_, *sc_, lastWordIdx);
			} else {
				const size_t W = alignEE8V256ShortWarmup<4>(*sc_);
				lrmax = ncol >= 4 * W ?
					alignEE8V256ShortSplit<4>(d.profbuf_.ptr(), d.mat_.ptr(), colstride,
						rf_, (size_t)rfi_, (size_t)rff_, *sc_, lastWordIdx, W) :
					alignEE8V256Short<4>(d.profbuf_.ptr(), d.mat_.ptr(), colstride,
						rf_, (size_t)rfi_, (size_t)rff_, *sc_, lastWordIdx);
			}
		}
		nfixup += (rff_ - rfi_) * iter;
		found = true;
#ifndef NDEBUG
		for(size_t i = (size_t)rfi_; i < (size_t)rff_; i++) {
			if((rand() & 15) == 0) {
				for(size_t k = 0; k < dpRows(); k++) {
					assert(cellOkEnd2EndU8(
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
		
		// Fetch the appropriate query profile.  Note that elements of rf_ must
		// be numbers, not masks.
		const int refc = (int)rf_[i];
		size_t off = (size_t)firsts5[refc] * iter * 2;
		pvScore = d.profbuf_.ptr() + off; // even elts = query profile, odd = gap barrier
		
		// Set all cells to low value
		vf = sse_xor_siall(vf, vf);

		// H vector from the final row of the previous column
		// Shift 2 bytes down so that topmost (least sig) cell gets 0
		vh = sse_slli_siall(vhlast, NBYTES_PER_WORD);
		// Fill topmost (least sig) cell with high value
		vh = sse_or_siall(vh, vhilsw);
		
		// For each character in the reference text:
		size_t j;
		for(j = 0; j < iter; j++) {
			vs0 = sse_load_siall(pvScore);
                        pvScore++;
			vs1 = sse_load_siall(pvScore);
                        pvScore++;
			// Load cells from E, calculated previously
			ve = sse_load_siall(pvELoad);
#if 0
			vhd = sse_load_siall(pvHLoad);
#endif
			assert_all_lt(ve, vhi);
			pvELoad += ROWSTRIDE;
			
			// Store cells in F, calculated previously
			vf = sse_subs_epu8(vf, vs1); // veto some ref gap extensions
			sse_store_siall(pvFStore, vf);
			pvFStore += ROWSTRIDE;
			
			// Factor in query profile (matches and mismatches)
			vh = sse_subs_epu8(vh, vs0);
			
			// Update H, factoring in E and F
			vh = sse_max_epu8(vh, ve);
			vh = sse_max_epu8(vh, vf);
			
			// Save the new vH values
			sse_store_siall(pvHStore, vh);
			pvHStore += ROWSTRIDE;
			
			// Update vE value
			vtmp = vh;
#if 0
			vhdtmp = vhd;
			vhd = sse_subs_epu8(vhd, rdgapo);
			vhd = sse_subs_epu8(vhd, vs1); // veto some read gap opens
			ve = sse_subs_epu8(ve, rdgape);
			ve = sse_max_epu8(ve, vhd);
#else
			vh = sse_subs_epu8(vh, rdgapo);
			vh = sse_subs_epu8(vh, vs1); // veto some read gap opens
			ve = sse_subs_epu8(ve, rdgape);
			ve = sse_max_epu8(ve, vh);
#endif
			assert_all_lt(ve, vhi);
			
			// Load the next h value
#if 0
			vh = vhdtmp;
#else
			vh = sse_load_siall(pvHLoad);
#endif
			pvHLoad += ROWSTRIDE;
			
			// Save E values
			sse_store_siall(pvEStore, ve);
			pvEStore += ROWSTRIDE;
			
			// Update vf value
			vtmp = sse_subs_epu8(vtmp, rfgapo);
			vf = sse_subs_epu8(vf, rfgape);
			assert_all_lt(vf, vhi);
			vf = sse_max_epu8(vf, vtmp);
		}
		// vf holds, in each lane, the F value the first pass carries out of
		// the lane's last row: max(F - rfgape, H - rfgapo).  Lane L's first
		// row receives the carry of lane L-1, or of any lane above it less the
		// cost of passing through the lanes in between.  vf becomes that
		// carry for every lane: a prefix scan in log2(lanes) steps.
		vf = sse_slli_siall(vf, NBYTES_PER_WORD);
		vtmp = sse_subs_epu8(sse_slli_siall(vf, 1 * NBYTES_PER_WORD), vpass1);
		vf = sse_max_epu8(vf, vtmp);
		vtmp = sse_subs_epu8(sse_slli_siall(vf, 2 * NBYTES_PER_WORD), vpass2);
		vf = sse_max_epu8(vf, vtmp);
		vtmp = sse_subs_epu8(sse_slli_siall(vf, 4 * NBYTES_PER_WORD), vpass4);
		vf = sse_max_epu8(vf, vtmp);
		vtmp = sse_subs_epu8(sse_slli_siall(vf, 8 * NBYTES_PER_WORD), vpass8);
		vf = sse_max_epu8(vf, vtmp);
#if NBYTES_PER_REG == 32
		vtmp = sse_subs_epu8(sse_slli_siall(vf, 16 * NBYTES_PER_WORD), vpass16);
		vf = sse_max_epu8(vf, vtmp);
#endif

		// One more pass down the column extends each lane's carry through its
		// rows and folds it into F, then H, then the next column's E: the
		// values the lazy-F loop converges to, in exactly iter steps.  The
		// carry reaching row j has lost j gap extensions, or everything if a
		// row up to j is inside the gap barrier; vpen accumulates that
		// penalty (saturating) apart from the carry, so that each row's value
		// is one subtraction from the carry rather than a chain of them.
		{
			const SSERegI *pvBar = pvScore - (iter*2) + 1; // gap barrier, row 0
			SSERegI *pvF = pvFStore - colstride;            // this column
			SSERegI *pvE = pvELoad;                         // next column
			// The first column's pvHLoad started in the staging vectors, so
			// take this column's H from pvHStore instead.
			pvHLoad = pvHStore - colstride;
			SSERegI *pvH = pvHLoad;                         // this column
			SSERegI vpen = sse_setzero_siall();
			for(j = 0; j < iter; j++) {
				vs1 = sse_load_siall(pvBar);
				pvBar += 2;
				vpen = sse_adds_epu8(vpen, vs1); // veto some ref gap extensions
				vtmp = sse_load_siall(pvF);
				vtmp = sse_max_epu8(vtmp, sse_subs_epu8(vf, vpen));
				sse_store_siall(pvF, vtmp);
				pvF += ROWSTRIDE;
				vh = sse_load_siall(pvH);
				vh = sse_max_epu8(vh, vtmp);
				sse_store_siall(pvH, vh);
				pvH += ROWSTRIDE;
				vhlast = vh;
				// Update E in case it can be improved using our new vh
				ve = sse_load_siall(pvE);
				vh = sse_subs_epu8(vh, rdgapo);
				vh = sse_subs_epu8(vh, vs1); // veto some read gap opens
				ve = sse_max_epu8(ve, vh);
				sse_store_siall(pvE, ve);
				pvE += ROWSTRIDE;
				vpen = sse_adds_epu8(vpen, rfgape);
			}
			nfixup += iter;
		}
		
#ifndef NDEBUG
		if(true && (rand() & 15) == 0) {
			// This is a work-intensive sanity check; each time we finish filling
			// a column, we check that each H, E, and F is sensible.
			for(size_t k = 0; k < dpRows(); k++) {
				assert(cellOkEnd2EndU8(
					d,
					k,                   // row
					i - rfi_,            // col
					refc,                // reference mask
					(int)(*rd_)[rdi_+k], // read char
					(int)(*qu_)[rdi_+k], // read quality
					*sc_));              // scoring scheme
			}
		}
#endif
		
		// Note: we may not want to extract from the final row
		TCScore lr = ((TCScore*)(pvHLoad))[lastWordIdx];
		found = true;
		if(lr > lrmax) {
			lrmax = lr;
		}

		// pvHLoad (this column's H), pvELoad (next column's E), pvHStore and
		// pvFStore (next column) and pvEStore (the column after) are already
		// where the next column needs them.
	}
	
	// Update metrics
	if(!debug) {
		size_t ninner = (rff_ - rfi_) * iter;
		met.col   += (rff_ - rfi_);             // DP columns
		met.cell  += (ninner * NWORDS_PER_REG); // DP cells
		met.inner += ninner;                    // DP inner loop iters
		met.fixup += nfixup;                    // DP fixup loop iters
	}
	
	flag = 0;
	
	// Did we find a solution?
	TAlScore score = MIN_I64;
	if(!found) {
		flag = -1; // no
		if(!debug) met.dpfail++;
		return MIN_I64;
	} else {
		score = (TAlScore)(lrmax - 0xff);
		if(score < minsc_) {
			flag = -1; // no
			if(!debug) met.dpfail++;
			return score;
		}
	}
	
	// Could we have saturated?
	if(lrmax == MIN_U8) {
		flag = -2; // yes
		if(!debug) met.dpsat++;
		return MIN_I64;
	}
	
	// Return largest score
	if(!debug) met.dpsucc++;
	return score;
}

/**
 * Given a filled-in DP table, populate the btncand_ list with candidate cells
 * that might be at the ends of valid alignments.  No need to do this unless
 * the maximum score returned by the align*() func is >= the minimum.
 *
 * Only cells that are exhaustively scored are candidates.  Those are the
 * cells inside the shape made of o's in this:
 *
 *  |-maxgaps-|
 *  *********************************    -
 *   ********************************    |
 *    *******************************    |
 *     ******************************    |
 *      *****************************    |
 *       **************************** read len
 *        ***************************    |
 *         **************************    |
 *          *************************    |
 *           ************************    |
 *            ***********oooooooooooo    -
 *            |-maxgaps-|
 *  |-readlen-|
 *  |-------skip--------|
 *
 * And it's possible for the shape to be truncated on the left and right sides.
 *
 * 
 */
bool SwAligner::gatherCellsNucleotidesEnd2EndSseU8(TAlScore best) {
	// What's the minimum number of rows that can possibly be spanned by an
	// alignment that meets the minimum score requirement?
	assert(sse8succ_);
	const size_t ncol = rff_ - rfi_;
	const size_t nrow = dpRows();
	assert_gt(nrow, 0);
	btncand_.clear();
	btncanddone_.clear();
	SSEData& d = fw_ ? sseU8fw_ : sseU8rc_;
	SSEMetrics& met = extend_ ? sseU8ExtendMet_ : sseU8MateMet_;
	assert(!d.profbuf_.empty());
	const size_t colstride = d.mat_.colstride();
	ASSERT_ONLY(bool sawbest = false);
	SSERegI *pvH = d.mat_.hvec(d.lastIter_, 0);
	for(size_t j = 0; j < ncol; j++) {
		TAlScore sc = (TAlScore)(((TCScore*)pvH)[d.lastWord_] - 0xff);
		assert_leq(sc, best);
		ASSERT_ONLY(sawbest = (sawbest || sc == best));
		if(sc >= minsc_) {
			// Yes, this is legit
			met.gathsol++;
			btncand_.expand();
			btncand_.back().init(nrow-1, j, sc);
		}
		pvH += colstride;
	}
	assert(sawbest);
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
bool SwAligner::backtraceNucleotidesEnd2EndSseU8(
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
	assert(!d.profbuf_.empty());
	assert_lt(row, rd_->length());
	btnstack_.clear(); // empty the backtrack stack
	btcells_.clear();  // empty the cells-so-far list
	AlnScore score;
	score.reset();
	score.basesAligned_ = std::numeric_limits<int>::min();
	score.edits_ = std::numeric_limits<int>::min();
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
		// This cell's mask word (see SSEMatrix): bit 0 reported through,
		// bit 1 H mask set and bits 2-6 the H mask, bit 7 E mask set and
		// bits 8-9 the E mask, bit 10 F mask set and bits 11-12 the F mask
		uint16_t& mk = btmasks[row][col];
		assert_eq(&mk, &d.mat_.masks_[row][col]);
		reportedThru = (mk & (1 << 0)) != 0;
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
				const TAlScore floorsc = MIN_I64;
				const int offsetsc = -0xff;
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
					if(sc_h_left > floorsc && sc_h_left - rdgapo == sc_cur) {
						mask |= (1 << 0);
					}
					// Get E score of cell to the left
					TAlScore sc_e_left = ((TCScore)btleft[0]) + offsetsc;
					if(sc_e_left > floorsc && sc_e_left - rdgape == sc_cur) {
						mask |= (1 << 1);
					}
					origMask = mask;
					assert(origMask > 0 || sc_cur <= sc_->match());
					if((mk & (1 << 7)) != 0) {
						mask = (mk >> 8) & 3;
					}
					if(mask == 3) {
						// Pick H -> E cell
						cur = SW_BT_OALL_READ_OPEN;
						{ mk &= ~(7 << 7); mk |= (1 << 7 | (2) << 8); } // might choose E later
						branch = true;
					} else if(mask == 2) {
						// I chose the E cell
						cur = SW_BT_RDGAP_EXTEND;
						{ mk &= ~(7 << 7); mk |= (1 << 7 | (0) << 8); } // done
					} else if(mask == 1) {
						// I chose the H cell
						cur = SW_BT_OALL_READ_OPEN;
						{ mk &= ~(7 << 7); mk |= (1 << 7 | (0) << 8); } // done
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
					if((mk & (1 << 10)) != 0) {
						mask = (mk >> 11) & 3;
					}
					if(mask == 3) {
						// I chose the H cell
						cur = SW_BT_OALL_REF_OPEN;
						{ mk &= ~(7 << 10); mk |= (1 << 10 | (2) << 11); } // might choose E later
						branch = true;
					} else if(mask == 2) {
						// I chose the F cell
						cur = SW_BT_RFGAP_EXTEND;
						{ mk &= ~(7 << 10); mk |= (1 << 10 | (0) << 11); } // done
					} else if(mask == 1) {
						// I chose the H cell
						cur = SW_BT_OALL_REF_OPEN;
						{ mk &= ~(7 << 10); mk |= (1 << 10 | (0) << 11); } // done
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
					if((mk & (1 << 1)) != 0) {
						mask = (mk >> 2) & 31;
					}
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
						{ mk &= ~(31 << 1); mk |= (1 << 1 | rest << 2); }
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
		// Mark the cell reported through, and record it in the row's dirty
		// interval (SSEMatrix::markDirty)
		mk |= (1 << 0);
		btdirtylo[row] = (col < btdirtylo[row]) ? col : btdirtylo[row];
		btdirtyhi[row] = (col >= btdirtyhi[row]) ? col + 1 : btdirtyhi[row];
		assert_eq(gaps, Edit::numGaps(ned));
		assert_leq(gaps, rdgap_ + rfgap_);
		// Cell was involved in a previously-reported alignment?
		if(!canMoveThru) {
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
				BT_SET_ROW_COL(row, col);
				continue;
			} else {
				// No branch points to revisit; just give up
				met.btcell += ncell;
				nbts += ncell;
				res.reset();
				met.btfail++; // DP backtraces failed
				return false;
			}
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
			assert_eq(0, trimBeg);
			assert_eq(btcells_.size(), dpRows() - trimBeg - trimEnd + readGaps);
			break;
		}
		if(branch) {
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
				(int)ct);
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
	assert_eq(0, trimBeg);
	assert_eq(0, trimEnd);
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
		reflen_,                  // length of entire reference
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
