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

static const size_t NBYTES_PER_REG  = 16;
static const size_t NWORDS_PER_REG  = 8;
static const size_t NBITS_PER_WORD  = 16;
static const size_t NBYTES_PER_WORD = 2;

// In 16-bit local mode, we have the option of using signed saturated
// arithmetic.  Because we have signed arithmetic, there's no need to
// add/subtract bias when building an applying the query profile.  The lowest
// value we can use is 0x8000, greatest is 0x7fff.

typedef int16_t TCScore;

/**
 * Build query profile look up tables for the read.  The query profile look
 * up table is organized as a 1D array indexed by [i][j] where i is the
 * reference character in the current DP column (0=A, 1=C, etc), and j is
 * the segment of the query we're currently working on.
 */
void SwAligner::buildQueryProfileLocalSseI16(bool fw) {
	bool& done = fw ? sseI16fwBuilt_ : sseI16rcBuilt_;
	if(done) {
		return;
	}
	done = true;
	const BTDnaString* rd = fw ? rdfw_ : rdrc_;
	const BTString* qu = fw ? qufw_ : qurc_;
	const size_t len = rd->length();
	const size_t seglen = (len + (NWORDS_PER_REG-1)) / NWORDS_PER_REG;
	// How many __m128i's are needed
	size_t n128s =
		64 +                    // slack bytes, for alignment?
		(seglen * ALPHA_SIZE)   // query profile data
		* 2;                    // & gap barrier data
	assert_gt(n128s, 0);
	SSEData& d = fw ? sseI16fw_ : sseI16rc_;
	d.profbuf_.resizeNoCopy(n128s);
	assert(!d.profbuf_.empty());
	d.maxPen_      = d.maxBonus_ = 0;
	d.lastIter_    = d.lastWord_ = 0;
	d.qprofStride_ = d.gbarStride_ = 2;
	d.bias_ = 0; // no bias when words are signed
	// For each reference character A, C, G, T, N ...
	for(size_t refc = 0; refc < ALPHA_SIZE; refc++) {
		// For each segment ...
		for(size_t i = 0; i < seglen; i++) {
			size_t j = i;
			int16_t *qprofWords =
				reinterpret_cast<int16_t*>(d.profbuf_.ptr() + (refc * seglen * 2) + (i * 2));
			int16_t *gbarWords =
				reinterpret_cast<int16_t*>(d.profbuf_.ptr() + (refc * seglen * 2) + (i * 2) + 1);
			// For each sub-word (byte) ...
			for(size_t k = 0; k < NWORDS_PER_REG; k++) {
				int sc = 0;
				*gbarWords = 0;
				if(j < len) {
					int readc = (*rd)[j];
					int readq = (*qu)[j];
					sc = sc_->score(readc, (int)(1 << refc), readq - 33);
					size_t j_from_end = len - j - 1;
					if(j < (size_t)sc_->gapbar ||
					   j_from_end < (size_t)sc_->gapbar)
					{
						// Inside the gap barrier
						*gbarWords = 0x8000; // add this twice
					}
				}
				if(refc == 0 && j == len-1) {
					// Remember which 128-bit word and which smaller word has
					// the final row
					d.lastIter_ = i;
					d.lastWord_ = k;
				}
				if(sc < 0) {
					if((size_t)(-sc) > d.maxPen_) {
						d.maxPen_ = (size_t)(-sc);
					}
				} else {
					if((size_t)sc > d.maxBonus_) {
						d.maxBonus_ = (size_t)sc;
					}
				}
				*qprofWords = (int16_t)sc;
				gbarWords++;
				qprofWords++;
				j += seglen; // update offset into query
			}
		}
	}
}

#ifndef NDEBUG
/**
 * Return true iff the cell has sane E/F/H values w/r/t its predecessors.
 */
static bool cellOkLocalI16(
	SSEData& d,
	size_t row,
	size_t col,
	int refc,
	int readc,
	int readq,
	const Scoring& sc)     // scoring scheme
{
	TCScore floorsc = MIN_I16;
	TCScore ceilsc = MIN_I16-1;
	TAlScore offsetsc = 0x8000;
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
	__m128i z = _mm_setzero_si128(); \
	__m128i tmp = _mm_setzero_si128(); \
	z = _mm_xor_si128(z, z); \
	tmp = _mm_cmpeq_epi16(x, z); \
	assert_eq(0xffff, _mm_movemask_epi8(tmp)); \
}

#define assert_all_gt(x, y) { \
	__m128i tmp = _mm_cmpgt_epi16(x, y); \
	assert_eq(0xffff, _mm_movemask_epi8(tmp)); \
}

#define assert_all_gt_lo(x) { \
	__m128i z = _mm_setzero_si128(); \
	__m128i tmp = _mm_setzero_si128(); \
	z = _mm_xor_si128(z, z); \
	tmp = _mm_cmpgt_epi16(x, z); \
	assert_eq(0xffff, _mm_movemask_epi8(tmp)); \
}

#define assert_all_lt(x, y) { \
	__m128i tmp = _mm_cmplt_epi16(x, y); \
	assert_eq(0xffff, _mm_movemask_epi8(tmp)); \
}

#define assert_all_leq(x, y) { \
	__m128i tmp = _mm_cmpgt_epi16(x, y); \
	assert_eq(0x0000, _mm_movemask_epi8(tmp)); \
}

#define assert_all_lt_hi(x) { \
	__m128i z = _mm_setzero_si128(); \
	__m128i tmp = _mm_setzero_si128(); \
	z = _mm_cmpeq_epi16(z, z); \
	z = _mm_srli_epi16(z, 1); \
	tmp = _mm_cmplt_epi16(x, z); \
	assert_eq(0xffff, _mm_movemask_epi8(tmp)); \
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
TAlScore SwAligner::alignGatherLoc16(int& flag, bool debug) {
	assert_leq(rdf_, rd_->length());
	assert_leq(rdf_, qu_->length());
	assert_lt(rfi_, rff_);
	assert_lt(rdi_, rdf_);
	assert_eq(rd_->length(), qu_->length());
	assert_geq(sc_->gapbar, 1);
	assert_gt(minsc_, 0);
	assert_leq(minsc_, MAX_I16);
	assert(repOk());
#ifndef NDEBUG
	for(size_t i = rfi_; i < rff_; i++) {
		assert_range(0, 16, (int)rf_[i]);
	}
#endif

	SSEData& d = fw_ ? sseI16fw_ : sseI16rc_;
	SSEMetrics& met = extend_ ? sseI16ExtendMet_ : sseI16MateMet_;
	if(!debug) met.dp++;
	buildQueryProfileLocalSseI16(fw_);
	assert(!d.profbuf_.empty());

	assert_gt(d.maxBonus_, 0);
	size_t iter =
		(dpRows() + (NWORDS_PER_REG-1)) / NWORDS_PER_REG; // iter = segLen
	
	// Now set up the score vectors.  We just need two columns worth, which
	// we'll call "left" and "right".
	d.vecbuf_.resize(ROWSTRIDE_2COL * iter * 2);
	d.vecbuf_.zero();
	__m128i *vbuf_l = d.vecbuf_.ptr();
	__m128i *vbuf_r = d.vecbuf_.ptr() + (ROWSTRIDE_2COL * iter);
	
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
		false,         // matrix cells have 8-bit scores?
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
	__m128i rfgapo   = _mm_setzero_si128();
	__m128i rfgape   = _mm_setzero_si128();
	__m128i rdgapo   = _mm_setzero_si128();
	__m128i rdgape   = _mm_setzero_si128();
	__m128i vlo      = _mm_setzero_si128();
	__m128i vhi      = _mm_setzero_si128();
	__m128i vlolsw   = _mm_setzero_si128();
	__m128i vmax     = _mm_setzero_si128();
	__m128i vcolmax  = _mm_setzero_si128();
	__m128i vmaxtmp  = _mm_setzero_si128();
	__m128i ve       = _mm_setzero_si128();
	__m128i vf       = _mm_setzero_si128();
	__m128i vh       = _mm_setzero_si128();
	__m128i vhd      = _mm_setzero_si128();
	__m128i vhdtmp   = _mm_setzero_si128();
	__m128i vtmp     = _mm_setzero_si128();
	__m128i vzero    = _mm_setzero_si128();
	__m128i vminsc   = _mm_setzero_si128();

	assert_gt(sc_->refGapOpen(), 0);
	assert_leq(sc_->refGapOpen(), MAX_I16);
	rfgapo = _mm_insert_epi16(rfgapo, sc_->refGapOpen(), 0);
	rfgapo = _mm_shufflelo_epi16(rfgapo, 0);
	rfgapo = _mm_shuffle_epi32(rfgapo, 0);
	
	// Set all elts to reference gap extension penalty
	assert_gt(sc_->refGapExtend(), 0);
	assert_leq(sc_->refGapExtend(), MAX_I16);
	assert_leq(sc_->refGapExtend(), sc_->refGapOpen());
	rfgape = _mm_insert_epi16(rfgape, sc_->refGapExtend(), 0);
	rfgape = _mm_shufflelo_epi16(rfgape, 0);
	rfgape = _mm_shuffle_epi32(rfgape, 0);

	// Set all elts to read gap open penalty
	assert_gt(sc_->readGapOpen(), 0);
	assert_leq(sc_->readGapOpen(), MAX_I16);
	rdgapo = _mm_insert_epi16(rdgapo, sc_->readGapOpen(), 0);
	rdgapo = _mm_shufflelo_epi16(rdgapo, 0);
	rdgapo = _mm_shuffle_epi32(rdgapo, 0);
	
	// Set all elts to read gap extension penalty
	assert_gt(sc_->readGapExtend(), 0);
	assert_leq(sc_->readGapExtend(), MAX_I16);
	assert_leq(sc_->readGapExtend(), sc_->readGapOpen());
	rdgape = _mm_insert_epi16(rdgape, sc_->readGapExtend(), 0);
	rdgape = _mm_shufflelo_epi16(rdgape, 0);
	rdgape = _mm_shuffle_epi32(rdgape, 0);
	
	// Set all elts to minimum score threshold.  Actually, to 1 less than the
	// threshold so we can use gt instead of geq.
	vminsc = _mm_insert_epi16(vminsc, (int)minsc_-1, 0);
	vminsc = _mm_shufflelo_epi16(vminsc, 0);
	vminsc = _mm_shuffle_epi32(vminsc, 0);

	// Set all elts to 0x8000 (min value for signed 16-bit)
	vlo = _mm_cmpeq_epi16(vlo, vlo);             // all elts = 0xffff
	vlo = _mm_slli_epi16(vlo, NBITS_PER_WORD-1); // all elts = 0x8000
	
	// Set all elts to 0x7fff (max value for signed 16-bit)
	vhi = _mm_cmpeq_epi16(vhi, vhi);             // all elts = 0xffff
	vhi = _mm_srli_epi16(vhi, 1);                // all elts = 0x7fff
	
	// Set all elts to 0x8000 (min value for signed 16-bit)
	vmax = vlo;
	
	// vlolsw: topmost (least sig) word set to 0x8000, all other words=0
	vlolsw = _mm_shuffle_epi32(vlo, 0);
	vlolsw = _mm_srli_si128(vlolsw, NBYTES_PER_REG - NBYTES_PER_WORD);
	
	// Points to a long vector of __m128i where each element is a block of
	// contiguous cells in the E, F or H matrix.  If the index % 3 == 0, then
	// the block of cells is from the E matrix.  If index % 3 == 1, they're
	// from the F matrix.  If index % 3 == 2, then they're from the H matrix.
	// Blocks of cells are organized in the same interleaved manner as they are
	// calculated by the Farrar algorithm.
	const __m128i *pvScore; // points into the query profile

	const size_t colstride = ROWSTRIDE_2COL * iter;
	
	// Initialize the H and E vectors in the first matrix column
	__m128i *pvELeft = vbuf_l + 0; __m128i *pvERight = vbuf_r + 0;
	//__m128i *pvFLeft = vbuf_l + 1;
	__m128i *pvFRight = vbuf_r + 1;
	__m128i *pvHLeft = vbuf_l + 2; __m128i *pvHRight = vbuf_r + 2;
	
	for(size_t i = 0; i < iter; i++) {
		// start low in local mode
		_mm_store_si128(pvERight, vlo); pvERight += ROWSTRIDE_2COL;
		_mm_store_si128(pvHRight, vlo); pvHRight += ROWSTRIDE_2COL;
		// Note: right and left are going to be swapped as soon as we enter
		// the outer loop below
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
	for(size_t i = rfi_; i < rff_; i++) {
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
		vh = _mm_load_si128(pvHLeft + colstride - ROWSTRIDE_2COL);
		
		// Set all F cells to low value
		vf = _mm_cmpeq_epi16(vf, vf);
		vf = _mm_slli_epi16(vf, NBITS_PER_WORD-1);
		vf = _mm_or_si128(vf, vlolsw);
		// vf now contains the vertical contribution

		// Store cells in F, calculated previously
		// No need to veto ref gap extensions, they're all 0x8000s
		_mm_store_si128(pvFRight, vf);
		pvFRight += ROWSTRIDE_2COL;
		
		// Shift down so that topmost (least sig) cell gets 0
		vh = _mm_slli_si128(vh, NBYTES_PER_WORD);
		// Fill topmost (least sig) cell with low value
		vh = _mm_or_si128(vh, vlolsw);
		
		// We pull out one loop iteration to make it easier to veto values in the top row
		
		// Load cells from E, calculated previously
		ve = _mm_load_si128(pvELeft);
		vhd = _mm_load_si128(pvHLeft);
		assert_all_lt(ve, vhi);
		pvELeft += ROWSTRIDE_2COL;
		// ve now contains the horizontal contribution
		
		// Factor in query profile (matches and mismatches)
		vh = _mm_adds_epi16(vh, pvScore[0]);
		// vh now contains the diagonal contribution
		
		// Update vE value
		vhdtmp = vhd;
		vhd = _mm_subs_epi16(vhd, rdgapo);
		vhd = _mm_adds_epi16(vhd, pvScore[1]); // veto some read gap opens
		vhd = _mm_adds_epi16(vhd, pvScore[1]); // veto some read gap opens
		ve = _mm_subs_epi16(ve, rdgape);
		ve = _mm_max_epi16(ve, vhd);

		// Update H, factoring in E and F
		vh = _mm_max_epi16(vh, ve);
		// F won't change anything!

		vf = vh;

		// Update highest score so far
		vcolmax = vh;
		
		// Save the new vH values
		_mm_store_si128(pvHRight, vh);

		assert_all_lt(ve, vhi);

		vh = vhdtmp;

		assert_all_lt(ve, vhi);
		pvHRight += ROWSTRIDE_2COL;
		pvHLeft += ROWSTRIDE_2COL;
		
		// Save E values
		_mm_store_si128(pvERight, ve);
		pvERight += ROWSTRIDE_2COL;
		
		// Update vf value
		vf = _mm_subs_epi16(vf, rfgapo);
		assert_all_lt(vf, vhi);
		
		pvScore += 2; // move on to next query profile

		// For each character in the reference text:
		size_t j;
		for(j = 1; j < iter; j++) {
			// Load cells from E, calculated previously
			ve = _mm_load_si128(pvELeft);
			vhd = _mm_load_si128(pvHLeft);
			assert_all_lt(ve, vhi);
			pvELeft += ROWSTRIDE_2COL;
			
			// Store cells in F, calculated previously
			vf = _mm_adds_epi16(vf, pvScore[1]); // veto some ref gap extensions
			vf = _mm_adds_epi16(vf, pvScore[1]); // veto some ref gap extensions
			_mm_store_si128(pvFRight, vf);
			pvFRight += ROWSTRIDE_2COL;
			
			// Factor in query profile (matches and mismatches)
			vh = _mm_adds_epi16(vh, pvScore[0]);
			vh = _mm_max_epi16(vh, vf);
			
			// Update vE value
			vhdtmp = vhd;
			vhd = _mm_subs_epi16(vhd, rdgapo);
			vhd = _mm_adds_epi16(vhd, pvScore[1]); // veto some read gap opens
			vhd = _mm_adds_epi16(vhd, pvScore[1]); // veto some read gap opens
			ve = _mm_subs_epi16(ve, rdgape);
			ve = _mm_max_epi16(ve, vhd);
			
			vh = _mm_max_epi16(vh, ve);
			vtmp = vh;
			
			// Update highest score encountered this far
			vcolmax = _mm_max_epi16(vcolmax, vh);
			
			// Save the new vH values
			_mm_store_si128(pvHRight, vh);

			vh = vhdtmp;

			assert_all_lt(ve, vhi);
			pvHRight += ROWSTRIDE_2COL;
			pvHLeft += ROWSTRIDE_2COL;
			
			// Save E values
			_mm_store_si128(pvERight, ve);
			pvERight += ROWSTRIDE_2COL;
			
			// Update vf value
			vtmp = _mm_subs_epi16(vtmp, rfgapo);
			vf = _mm_subs_epi16(vf, rfgape);
			assert_all_lt(vf, vhi);
			vf = _mm_max_epi16(vf, vtmp);
			
			pvScore += 2; // move on to next query profile / gap veto
		}
		// pvHStore, pvELoad, pvEStore have all rolled over to the next column
		pvFRight -= colstride; // reset to start of column
		vtmp = _mm_load_si128(pvFRight);
		
		pvHRight -= colstride; // reset to start of column
		vh = _mm_load_si128(pvHRight);
		
		pvScore = d.profbuf_.ptr() + off + 1; // reset veto vector
		
		// vf from last row gets shifted down by one to overlay the first row
		// rfgape has already been subtracted from it.
		vf = _mm_slli_si128(vf, NBYTES_PER_WORD);
		vf = _mm_or_si128(vf, vlolsw);
		
		vf = _mm_adds_epi16(vf, *pvScore); // veto some ref gap extensions
		vf = _mm_adds_epi16(vf, *pvScore); // veto some ref gap extensions
		vf = _mm_max_epi16(vtmp, vf);
		vtmp = _mm_cmpgt_epi16(vf, vtmp);
		int cmp = _mm_movemask_epi8(vtmp);
		
		// If any element of vtmp is greater than H - gap-open...
		j = 0;
		while(cmp != 0x0000) {
			// Store this vf
			_mm_store_si128(pvFRight, vf);
			pvFRight += ROWSTRIDE_2COL;
			
			// Update vh w/r/t new vf
			vh = _mm_max_epi16(vh, vf);
			
			// Save vH values
			_mm_store_si128(pvHRight, vh);
			pvHRight += ROWSTRIDE_2COL;
			
			// Update highest score encountered so far.
			vcolmax = _mm_max_epi16(vcolmax, vh);

			pvScore += 2;
			
			assert_lt(j, iter);
			if(++j == iter) {
				pvFRight -= colstride;
				vtmp = _mm_load_si128(pvFRight);   // load next vf ASAP
				pvHRight -= colstride;
				vh = _mm_load_si128(pvHRight);     // load next vh ASAP
				pvScore = d.profbuf_.ptr() + off + 1;
				j = 0;
				vf = _mm_slli_si128(vf, NBYTES_PER_WORD);
				vf = _mm_or_si128(vf, vlolsw);
			} else {
				vtmp = _mm_load_si128(pvFRight);   // load next vf ASAP
				vh = _mm_load_si128(pvHRight);     // load next vh ASAP
			}
			
			// Update F with another gap extension
			vf = _mm_subs_epi16(vf, rfgape);
			vf = _mm_adds_epi16(vf, *pvScore); // veto some ref gap extensions
			vf = _mm_adds_epi16(vf, *pvScore); // veto some ref gap extensions
			vf = _mm_max_epi16(vtmp, vf);
			vtmp = _mm_cmpgt_epi16(vf, vtmp);
			cmp = _mm_movemask_epi8(vtmp);
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
				vh = _mm_load_si128(pvHLeft);
				vtmp = _mm_cmpgt_epi16(pvScore[0], vzero);
				int cmp = _mm_movemask_epi8(vtmp);
				if(cmp != 0) {
					// At least one candidate in this mask.  Now iterate
					// through vm/vh to evaluate individual cells.
					for(size_t m = 0; m < NWORDS_PER_REG; m++) {
						size_t row = k + m * iter;
						if(row >= dpRows()) {
							break;
						}
						TAlScore sc = (TAlScore)(((TCScore *)&vh)[m] + 0x8000);
						if(sc >= minsc_) {
							if(((TCScore *)&vtmp)[m] != 0) {
								// Add to data structure holding all candidates
								size_t col = i - rfi_ - 1; // -1 b/c prev col
								size_t frombot = dpRows() - row - 1;
								DpBtCandidate cand(row, col, sc);
								btdiag_.add(frombot + col, cand);
							}
						}
					}
				}
				pvHLeft += ROWSTRIDE_2COL;
				pvScore += 2;
			}
		}

		// Save some elements to checkpoints
		if(checkpoint) {
			
			__m128i *pvE = vbuf_r + 0;
			__m128i *pvF = vbuf_r + 1;
			__m128i *pvH = vbuf_r + 2;
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
							size_t vecoff = (row_mod << 5) + row_div;
							assert_lt(row_div, 8);
							int16_t h_sc = ((int16_t*)pvH)[vecoff];
							int16_t e_sc = ((int16_t*)pvE)[vecoff];
							int16_t f_sc = ((int16_t*)pvF)[vecoff];
							h_sc += 0x8000; assert_geq(h_sc, 0);
							e_sc += 0x8000; assert_geq(e_sc, 0);
							f_sc += 0x8000; assert_geq(f_sc, 0);
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
						size_t vecoff = (cpimod << 5) + cpidiv;
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
					int16_t h_sc = ((int16_t*)pvH)[vecoff];
					//int16_t e_sc = ((int16_t*)pvE)[vecoff];
					int16_t f_sc = ((int16_t*)pvF)[vecoff];
					h_sc += 0x8000; assert_geq(h_sc, 0);
					//e_sc += 0x8000; assert_geq(e_sc, 0);
					f_sc += 0x8000; assert_geq(f_sc, 0);
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
					__m128i *dst = cper_.qcols_.ptr() + coloff;
					memcpy(dst, vbuf_r, sizeof(__m128i) * wordspercol);
				}
			}
			if(cper_.debug_) {
				// Save the column using memcpys
				size_t wordspercol = cper_.niter_ * ROWSTRIDE_2COL;
				size_t coloff = coli * wordspercol;
				__m128i *dst = cper_.qcolsD_.ptr() + coloff;
				memcpy(dst, vbuf_r, sizeof(__m128i) * wordspercol);
			}
		}

		vmax = _mm_max_epi16(vmax, vcolmax);
		{
			// Get single largest score in this column
			vmaxtmp = vcolmax;
			vtmp = _mm_srli_si128(vmaxtmp, 8);
			vmaxtmp = _mm_max_epi16(vmaxtmp, vtmp);
			vtmp = _mm_srli_si128(vmaxtmp, 4);
			vmaxtmp = _mm_max_epi16(vmaxtmp, vtmp);
			vtmp = _mm_srli_si128(vmaxtmp, 2);
			vmaxtmp = _mm_max_epi16(vmaxtmp, vtmp);
			int16_t ret = _mm_extract_epi16(vmaxtmp, 0);
			TAlScore score = (TAlScore)(ret + 0x8000);
			if(ret == MIN_I16) {
				score = MIN_I64;
			}
			
			if(score < minsc_) {
				size_t ncolleft = rff_ - i - 1;
				if(max<TAlScore>(score, 0) + (TAlScore)ncolleft * matchsc < minsc_) {
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
			vh = _mm_load_si128(pvHLeft);
			vtmp = _mm_cmpgt_epi16(pvScore[0], vzero);
			int cmp = _mm_movemask_epi8(vtmp);
			if(cmp != 0) {
				// At least one candidate in this mask.  Now iterate
				// through vm/vh to evaluate individual cells.
				for(size_t m = 0; m < NWORDS_PER_REG; m++) {
					size_t row = k + m * iter;
					if(row >= dpRows()) {
						break;
					}
					TAlScore sc = (TAlScore)(((TCScore *)&vh)[m] + 0x8000);
					if(sc >= minsc_) {
						if(((TCScore *)&vtmp)[m] != 0) {
							// Add to data structure holding all candidates
							size_t col = rff_ - rfi_ - 1; // -1 b/c prev col
							size_t frombot = dpRows() - row - 1;
							DpBtCandidate cand(row, col, sc);
							btdiag_.add(frombot + col, cand);
						}
					}
				}
			}
			pvHLeft += ROWSTRIDE_2COL;
			pvScore += 2;
		}
	}

	// Find largest score in vmax
	vtmp = _mm_srli_si128(vmax, 8);
	vmax = _mm_max_epi16(vmax, vtmp);
	vtmp = _mm_srli_si128(vmax, 4);
	vmax = _mm_max_epi16(vmax, vtmp);
	vtmp = _mm_srli_si128(vmax, 2);
	vmax = _mm_max_epi16(vmax, vtmp);
	int16_t ret = _mm_extract_epi16(vmax, 0);

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
	if(ret == MIN_I16) {
		flag = -1; // no
		if(!debug) met.dpfail++;
		return MIN_I64;
	} else {
		score = (TAlScore)(ret + 0x8000);
		if(score < minsc_) {
			flag = -1; // no
			if(!debug) met.dpfail++;
			return score;
		}
	}
	
	// Could we have saturated?
	if(ret == MAX_I16) {
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

/**
 * Solve the current alignment problem using SSE instructions that operate on 8
 * signed 16-bit values packed into a single 128-bit register.
 */
TAlScore SwAligner::alignNucleotidesLocalSseI16(int& flag, bool debug) {
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

	SSEData& d = fw_ ? sseI16fw_ : sseI16rc_;
	SSEMetrics& met = extend_ ? sseI16ExtendMet_ : sseI16MateMet_;
	if(!debug) met.dp++;
	buildQueryProfileLocalSseI16(fw_);
	assert(!d.profbuf_.empty());

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
	__m128i rfgapo   = _mm_setzero_si128();
	__m128i rfgape   = _mm_setzero_si128();
	__m128i rdgapo   = _mm_setzero_si128();
	__m128i rdgape   = _mm_setzero_si128();
	__m128i vlo      = _mm_setzero_si128();
	__m128i vhi      = _mm_setzero_si128();
	__m128i vlolsw   = _mm_setzero_si128();
	__m128i vmax     = _mm_setzero_si128();
	__m128i vcolmax  = _mm_setzero_si128();
	__m128i vmaxtmp  = _mm_setzero_si128();
	__m128i ve       = _mm_setzero_si128();
	__m128i vf       = _mm_setzero_si128();
	__m128i vh       = _mm_setzero_si128();
	__m128i vtmp     = _mm_setzero_si128();

	assert_gt(sc_->refGapOpen(), 0);
	assert_leq(sc_->refGapOpen(), MAX_I16);
	rfgapo = _mm_insert_epi16(rfgapo, sc_->refGapOpen(), 0);
	rfgapo = _mm_shufflelo_epi16(rfgapo, 0);
	rfgapo = _mm_shuffle_epi32(rfgapo, 0);
	
	// Set all elts to reference gap extension penalty
	assert_gt(sc_->refGapExtend(), 0);
	assert_leq(sc_->refGapExtend(), MAX_I16);
	assert_leq(sc_->refGapExtend(), sc_->refGapOpen());
	rfgape = _mm_insert_epi16(rfgape, sc_->refGapExtend(), 0);
	rfgape = _mm_shufflelo_epi16(rfgape, 0);
	rfgape = _mm_shuffle_epi32(rfgape, 0);

	// Set all elts to read gap open penalty
	assert_gt(sc_->readGapOpen(), 0);
	assert_leq(sc_->readGapOpen(), MAX_I16);
	rdgapo = _mm_insert_epi16(rdgapo, sc_->readGapOpen(), 0);
	rdgapo = _mm_shufflelo_epi16(rdgapo, 0);
	rdgapo = _mm_shuffle_epi32(rdgapo, 0);
	
	// Set all elts to read gap extension penalty
	assert_gt(sc_->readGapExtend(), 0);
	assert_leq(sc_->readGapExtend(), MAX_I16);
	assert_leq(sc_->readGapExtend(), sc_->readGapOpen());
	rdgape = _mm_insert_epi16(rdgape, sc_->readGapExtend(), 0);
	rdgape = _mm_shufflelo_epi16(rdgape, 0);
	rdgape = _mm_shuffle_epi32(rdgape, 0);

	// Set all elts to 0x8000 (min value for signed 16-bit)
	vlo = _mm_cmpeq_epi16(vlo, vlo);             // all elts = 0xffff
	vlo = _mm_slli_epi16(vlo, NBITS_PER_WORD-1); // all elts = 0x8000
	
	// Set all elts to 0x7fff (max value for signed 16-bit)
	vhi = _mm_cmpeq_epi16(vhi, vhi);             // all elts = 0xffff
	vhi = _mm_srli_epi16(vhi, 1);                // all elts = 0x7fff
	
	// Set all elts to 0x8000 (min value for signed 16-bit)
	vmax = vlo;
	
	// vlolsw: topmost (least sig) word set to 0x8000, all other words=0
	vlolsw = _mm_shuffle_epi32(vlo, 0);
	vlolsw = _mm_srli_si128(vlolsw, NBYTES_PER_REG - NBYTES_PER_WORD);
	
	// Points to a long vector of __m128i where each element is a block of
	// contiguous cells in the E, F or H matrix.  If the index % 3 == 0, then
	// the block of cells is from the E matrix.  If index % 3 == 1, they're
	// from the F matrix.  If index % 3 == 2, then they're from the H matrix.
	// Blocks of cells are organized in the same interleaved manner as they are
	// calculated by the Farrar algorithm.
	const __m128i *pvScore; // points into the query profile

	d.mat_.init(dpRows(), rff_ - rfi_, NWORDS_PER_REG);
	const size_t colstride = d.mat_.colstride();
	//const size_t rowstride = d.mat_.rowstride();
	assert_eq(ROWSTRIDE, colstride / iter);
	
	// Initialize the H and E vectors in the first matrix column
	__m128i *pvHTmp = d.mat_.tmpvec(0, 0);
	__m128i *pvETmp = d.mat_.evec(0, 0);
	
	for(size_t i = 0; i < iter; i++) {
		_mm_store_si128(pvETmp, vlo);
		_mm_store_si128(pvHTmp, vlo); // start low in local mode
		pvETmp += ROWSTRIDE;
		pvHTmp += ROWSTRIDE;
	}
	// These are swapped just before the innermost loop
	__m128i *pvHStore = d.mat_.hvec(0, 0);
	__m128i *pvHLoad  = d.mat_.tmpvec(0, 0);
	__m128i *pvELoad  = d.mat_.evec(0, 0);
	__m128i *pvEStore = d.mat_.evecUnsafe(0, 1);
	__m128i *pvFStore = d.mat_.fvec(0, 0);
	__m128i *pvFTmp   = NULL;
	
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
	for(size_t i = rfi_; i < rff_; i++) {
		assert(pvFStore == d.mat_.fvec(0, i - rfi_));
		assert(pvHStore == d.mat_.hvec(0, i - rfi_));
		
		// Fetch this column's reference mask
		const int refm = (int)rf_[i];
		
		// Fetch the appropriate query profile
		size_t off = (size_t)firsts5[refm] * iter * 2;
		pvScore = d.profbuf_.ptr() + off; // even elts = query profile, odd = gap barrier
		
		// Load H vector from the final row of the previous column
		vh = _mm_load_si128(pvHLoad + colstride - ROWSTRIDE);
		
		// Set all F cells to low value
		vf = _mm_cmpeq_epi16(vf, vf);
		vf = _mm_slli_epi16(vf, NBITS_PER_WORD-1);
		vf = _mm_or_si128(vf, vlolsw);
		// vf now contains the vertical contribution

		// Store cells in F, calculated previously
		// No need to veto ref gap extensions, they're all 0x8000s
		_mm_store_si128(pvFStore, vf);
		pvFStore += ROWSTRIDE;
		
		// Shift down so that topmost (least sig) cell gets 0
		vh = _mm_slli_si128(vh, NBYTES_PER_WORD);
		// Fill topmost (least sig) cell with low value
		vh = _mm_or_si128(vh, vlolsw);
		
		// We pull out one loop iteration to make it easier to veto values in the top row
		
		// Load cells from E, calculated previously
		ve = _mm_load_si128(pvELoad);
		assert_all_lt(ve, vhi);
		pvELoad += ROWSTRIDE;
		// ve now contains the horizontal contribution
		
		// Factor in query profile (matches and mismatches)
		vh = _mm_adds_epi16(vh, pvScore[0]);
		// vh now contains the diagonal contribution
		
		// Update H, factoring in E and F
		vtmp = _mm_max_epi16(vh, ve);
		// F won't change anything!
		
		vh = vtmp;
		
		// Update highest score so far
		vcolmax = vlo;
		vcolmax = _mm_max_epi16(vcolmax, vh);
		
		// Save the new vH values
		_mm_store_si128(pvHStore, vh);
		pvHStore += ROWSTRIDE;
		
		// Update vE value
		vf = vh;
		vh = _mm_subs_epi16(vh, rdgapo);
		vh = _mm_adds_epi16(vh, pvScore[1]); // veto some read gap opens
		vh = _mm_adds_epi16(vh, pvScore[1]); // veto some read gap opens
		ve = _mm_subs_epi16(ve, rdgape);
		ve = _mm_max_epi16(ve, vh);
		assert_all_lt(ve, vhi);
		
		// Load the next h value
		vh = _mm_load_si128(pvHLoad);
		pvHLoad += ROWSTRIDE;
		
		// Save E values
		_mm_store_si128(pvEStore, ve);
		pvEStore += ROWSTRIDE;
		
		// Update vf value
		vf = _mm_subs_epi16(vf, rfgapo);
		assert_all_lt(vf, vhi);
		
		pvScore += 2; // move on to next query profile

		// For each character in the reference text:
		size_t j;
		for(j = 1; j < iter; j++) {
			// Load cells from E, calculated previously
			ve = _mm_load_si128(pvELoad);
			assert_all_lt(ve, vhi);
			pvELoad += ROWSTRIDE;
			
			// Store cells in F, calculated previously
			vf = _mm_adds_epi16(vf, pvScore[1]); // veto some ref gap extensions
			vf = _mm_adds_epi16(vf, pvScore[1]); // veto some ref gap extensions
			_mm_store_si128(pvFStore, vf);
			pvFStore += ROWSTRIDE;
			
			// Factor in query profile (matches and mismatches)
			vh = _mm_adds_epi16(vh, pvScore[0]);
			
			// Update H, factoring in E and F
			vh = _mm_max_epi16(vh, ve);
			vh = _mm_max_epi16(vh, vf);
			
			// Update highest score encountered this far
			vcolmax = _mm_max_epi16(vcolmax, vh);
			
			// Save the new vH values
			_mm_store_si128(pvHStore, vh);
			pvHStore += ROWSTRIDE;
			
			// Update vE value
			vtmp = vh;
			vh = _mm_subs_epi16(vh, rdgapo);
			vh = _mm_adds_epi16(vh, pvScore[1]); // veto some read gap opens
			vh = _mm_adds_epi16(vh, pvScore[1]); // veto some read gap opens
			ve = _mm_subs_epi16(ve, rdgape);
			ve = _mm_max_epi16(ve, vh);
			assert_all_lt(ve, vhi);
			
			// Load the next h value
			vh = _mm_load_si128(pvHLoad);
			pvHLoad += ROWSTRIDE;
			
			// Save E values
			_mm_store_si128(pvEStore, ve);
			pvEStore += ROWSTRIDE;
			
			// Update vf value
			vtmp = _mm_subs_epi16(vtmp, rfgapo);
			vf = _mm_subs_epi16(vf, rfgape);
			assert_all_lt(vf, vhi);
			vf = _mm_max_epi16(vf, vtmp);
			
			pvScore += 2; // move on to next query profile / gap veto
		}
		// pvHStore, pvELoad, pvEStore have all rolled over to the next column
		pvFTmp = pvFStore;
		pvFStore -= colstride; // reset to start of column
		vtmp = _mm_load_si128(pvFStore);
		
		pvHStore -= colstride; // reset to start of column
		vh = _mm_load_si128(pvHStore);
		
		pvEStore -= colstride; // reset to start of column
		ve = _mm_load_si128(pvEStore);
		
		pvHLoad = pvHStore;    // new pvHLoad = pvHStore
		pvScore = d.profbuf_.ptr() + off + 1; // reset veto vector
		
		// vf from last row gets shifted down by one to overlay the first row
		// rfgape has already been subtracted from it.
		vf = _mm_slli_si128(vf, NBYTES_PER_WORD);
		vf = _mm_or_si128(vf, vlolsw);
		
		vf = _mm_adds_epi16(vf, *pvScore); // veto some ref gap extensions
		vf = _mm_adds_epi16(vf, *pvScore); // veto some ref gap extensions
		vf = _mm_max_epi16(vtmp, vf);
		vtmp = _mm_cmpgt_epi16(vf, vtmp);
		int cmp = _mm_movemask_epi8(vtmp);
		
		// If any element of vtmp is greater than H - gap-open...
		j = 0;
		while(cmp != 0x0000) {
			// Store this vf
			_mm_store_si128(pvFStore, vf);
			pvFStore += ROWSTRIDE;
			
			// Update vh w/r/t new vf
			vh = _mm_max_epi16(vh, vf);
			
			// Save vH values
			_mm_store_si128(pvHStore, vh);
			pvHStore += ROWSTRIDE;
			
			// Update highest score encountered this far
			vcolmax = _mm_max_epi16(vcolmax, vh);
			
			// Update E in case it can be improved using our new vh
			vh = _mm_subs_epi16(vh, rdgapo);
			vh = _mm_adds_epi16(vh, *pvScore); // veto some read gap opens
			vh = _mm_adds_epi16(vh, *pvScore); // veto some read gap opens
			ve = _mm_max_epi16(ve, vh);
			_mm_store_si128(pvEStore, ve);
			pvEStore += ROWSTRIDE;
			pvScore += 2;
			
			assert_lt(j, iter);
			if(++j == iter) {
				pvFStore -= colstride;
				vtmp = _mm_load_si128(pvFStore);   // load next vf ASAP
				pvHStore -= colstride;
				vh = _mm_load_si128(pvHStore);     // load next vh ASAP
				pvEStore -= colstride;
				ve = _mm_load_si128(pvEStore);     // load next ve ASAP
				pvScore = d.profbuf_.ptr() + off + 1;
				j = 0;
				vf = _mm_slli_si128(vf, NBYTES_PER_WORD);
				vf = _mm_or_si128(vf, vlolsw);
			} else {
				vtmp = _mm_load_si128(pvFStore);   // load next vf ASAP
				vh = _mm_load_si128(pvHStore);     // load next vh ASAP
				ve = _mm_load_si128(pvEStore);     // load next vh ASAP
			}
			
			// Update F with another gap extension
			vf = _mm_subs_epi16(vf, rfgape);
			vf = _mm_adds_epi16(vf, *pvScore); // veto some ref gap extensions
			vf = _mm_adds_epi16(vf, *pvScore); // veto some ref gap extensions
			vf = _mm_max_epi16(vtmp, vf);
			vtmp = _mm_cmpgt_epi16(vf, vtmp);
			cmp = _mm_movemask_epi8(vtmp);
			nfixup++;
		}
		
#ifndef NDEBUG
		if((rand() & 15) == 0) {
			// This is a work-intensive sanity check; each time we finish filling
			// a column, we check that each H, E, and F is sensible.
			for(size_t k = 0; k < dpRows(); k++) {
				assert(cellOkLocalI16(
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
		vmax = _mm_max_epi16(vmax, vcolmax);
		_mm_store_si128(d.mat_.tmpvec(0, i - rfi_), vcolmax);

		{
			// Get single largest score in this column
			vmaxtmp = vcolmax;
			vtmp = _mm_srli_si128(vmaxtmp, 8);
			vmaxtmp = _mm_max_epi16(vmaxtmp, vtmp);
			vtmp = _mm_srli_si128(vmaxtmp, 4);
			vmaxtmp = _mm_max_epi16(vmaxtmp, vtmp);
			vtmp = _mm_srli_si128(vmaxtmp, 2);
			vmaxtmp = _mm_max_epi16(vmaxtmp, vtmp);
			int16_t ret = _mm_extract_epi16(vmaxtmp, 0);
			TAlScore score = (TAlScore)(ret + 0x8000);
			
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

	// Find largest score in vmax
	vtmp = _mm_srli_si128(vmax, 8);
	vmax = _mm_max_epi16(vmax, vtmp);
	vtmp = _mm_srli_si128(vmax, 4);
	vmax = _mm_max_epi16(vmax, vtmp);
	vtmp = _mm_srli_si128(vmax, 2);
	vmax = _mm_max_epi16(vmax, vtmp);
	int16_t ret = _mm_extract_epi16(vmax, 0);

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
	if(ret == MIN_I16) {
		flag = -1; // no
		if(!debug) met.dpfail++;
		return MIN_I64;
	} else {
		score = (TAlScore)(ret + 0x8000);
		if(score < minsc_) {
			flag = -1; // no
			if(!debug) met.dpfail++;
			return score;
		}
	}
	
	// Could we have saturated?
	if(ret == MAX_I16) {
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
bool SwAligner::gatherCellsNucleotidesLocalSseI16(TAlScore best) {
	// What's the minimum number of rows that can possibly be spanned by an
	// alignment that meets the minimum score requirement?
	assert(sse16succ_);
	size_t bonus = (size_t)sc_->match(30);
	const size_t ncol = lastsolcol_ + 1;
	const size_t nrow = dpRows();
	assert_gt(nrow, 0);
	btncand_.clear();
	btncanddone_.clear();
	SSEData& d = fw_ ? sseI16fw_ : sseI16rc_;
	SSEMetrics& met = extend_ ? sseI16ExtendMet_ : sseI16MateMet_;
	assert(!d.profbuf_.empty());
	//const size_t rowstride = d.mat_.rowstride();
	//const size_t colstride = d.mat_.colstride();
	size_t iter = (dpRows() + (NWORDS_PER_REG - 1)) / NWORDS_PER_REG;
	assert_gt(iter, 0);
	assert_geq(minsc_, 0);
	assert_gt(bonus, 0);
	size_t minrow = (size_t)(((minsc_ + bonus - 1) / bonus) - 1);
	for(size_t j = 0; j < ncol; j++) {
		// Establish the range of rows where a backtrace from the cell in this
		// row/col is close enough to one of the core diagonals that it could
		// conceivably count
		size_t nrow_lo = MIN_SIZE_T;
		size_t nrow_hi = nrow;
		// First, check if there is a cell in this column with a score
		// above the score threshold
		__m128i vmax = *d.mat_.tmpvec(0, j);
		__m128i vtmp = _mm_srli_si128(vmax, 8);
		vmax = _mm_max_epi16(vmax, vtmp);
		vtmp = _mm_srli_si128(vmax, 4);
		vmax = _mm_max_epi16(vmax, vtmp);
		vtmp = _mm_srli_si128(vmax, 2);
		vmax = _mm_max_epi16(vmax, vtmp);
		TAlScore score = (TAlScore)((int16_t)_mm_extract_epi16(vmax, 0) + 0x8000);
		assert_geq(score, 0);
#ifndef NDEBUG
		{
			// Start in upper vector row and move down
			TAlScore max = 0;
			vmax = *d.mat_.tmpvec(0, j);
			__m128i *pvH = d.mat_.hvec(0, j);
			for(size_t i = 0; i < iter; i++) {
				for(size_t k = 0; k < NWORDS_PER_REG; k++) {
					TAlScore sc = (TAlScore)(((TCScore*)pvH)[k] + 0x8000);
					TAlScore scm = (TAlScore)(((TCScore*)&vmax)[k] + 0x8000);
					assert_leq(sc, scm);
					if(sc > max) {
						max = sc;
					}
				}
				pvH += ROWSTRIDE;
			}
			assert_eq(max, score);
		}
#endif
		if(score < minsc_) {
			// Scores in column aren't good enough
			continue;
		}
		// Get pointer to first cell in column to examine:
		__m128i *pvHorig = d.mat_.hvec(0, j);
		__m128i *pvH     = pvHorig;
		// Get pointer to the vector in the following column that corresponds
		// to the cells diagonally down and to the right from the cells in pvH
		__m128i *pvHSucc = (j < ncol-1) ? d.mat_.hvec(0, j+1) : NULL;
		// Start in upper vector row and move down
		for(size_t i = 0; i < iter; i++) {
			if(pvHSucc != NULL) {
				pvHSucc += ROWSTRIDE;
				if(i == iter-1) {
					pvHSucc = d.mat_.hvec(0, j+1);
				}
			}
			// Which elements of this vector are exhaustively scored?
			size_t rdoff = i;
			for(size_t k = 0; k < NWORDS_PER_REG; k++) {
				// Is this row, col one that we can potential backtrace from?
				// I.e. are we close enough to a core diagonal?
				if(rdoff >= nrow_lo && rdoff < nrow_hi) {
					// This cell has been exhaustively scored
					if(rdoff >= minrow) {
						// ... and it could potentially score high enough
						TAlScore sc = (TAlScore)(((TCScore*)pvH)[k] + 0x8000);
						assert_leq(sc, best);
						if(sc >= minsc_) {
							// This is a potential solution
							bool matchSucc = false;
							int readc = (*rd_)[rdoff];
							int refc = rf_[j + rfi_];
							bool match = ((refc & (1 << readc)) != 0);
							if(rdoff < dpRows()-1) {
								int readcSucc = (*rd_)[rdoff+1];
								int refcSucc = rf_[j + rfi_ + 1];
								assert_range(0, 16, refcSucc);
								matchSucc = ((refcSucc & (1 << readcSucc)) != 0);
							}
							if(match && !matchSucc) {
								// Yes, this is legit
								met.gathsol++;
								btncand_.expand();
								btncand_.back().init(rdoff, j, sc);
							}
						}
					}
				} else {
					// Already saw every element in the vector that's been
					// exhaustively scored
					break;
				}
				rdoff += iter;
			}
			pvH += ROWSTRIDE;
		}
	}
	if(!btncand_.empty()) {
		d.mat_.initMasks();
	}
	return !btncand_.empty();
}

#define MOVE_VEC_PTR_UP(vec, rowvec, rowelt) { \
	if(rowvec == 0) { \
		rowvec += d.mat_.nvecrow_; \
		vec += d.mat_.colstride_; \
		rowelt--; \
	} \
	rowvec--; \
	vec -= ROWSTRIDE; \
}

#define MOVE_VEC_PTR_LEFT(vec, rowvec, rowelt) { vec -= d.mat_.colstride_; }

#define MOVE_VEC_PTR_UPLEFT(vec, rowvec, rowelt) { \
 	MOVE_VEC_PTR_UP(vec, rowvec, rowelt); \
 	MOVE_VEC_PTR_LEFT(vec, rowvec, rowelt); \
}

#define MOVE_ALL_LEFT() { \
	MOVE_VEC_PTR_LEFT(cur_vec, rowvec, rowelt); \
	MOVE_VEC_PTR_LEFT(left_vec, left_rowvec, left_rowelt); \
	MOVE_VEC_PTR_LEFT(up_vec, up_rowvec, up_rowelt); \
	MOVE_VEC_PTR_LEFT(upleft_vec, upleft_rowvec, upleft_rowelt); \
}

#define MOVE_ALL_UP() { \
	MOVE_VEC_PTR_UP(cur_vec, rowvec, rowelt); \
	MOVE_VEC_PTR_UP(left_vec, left_rowvec, left_rowelt); \
	MOVE_VEC_PTR_UP(up_vec, up_rowvec, up_rowelt); \
	MOVE_VEC_PTR_UP(upleft_vec, upleft_rowvec, upleft_rowelt); \
}

#define MOVE_ALL_UPLEFT() { \
	MOVE_VEC_PTR_UPLEFT(cur_vec, rowvec, rowelt); \
	MOVE_VEC_PTR_UPLEFT(left_vec, left_rowvec, left_rowelt); \
	MOVE_VEC_PTR_UPLEFT(up_vec, up_rowvec, up_rowelt); \
	MOVE_VEC_PTR_UPLEFT(upleft_vec, upleft_rowvec, upleft_rowelt); \
}

#define NEW_ROW_COL(row, col) { \
	rowelt = row / d.mat_.nvecrow_; \
	rowvec = row % d.mat_.nvecrow_; \
	eltvec = (col * d.mat_.colstride_) + (rowvec * ROWSTRIDE); \
	cur_vec = d.mat_.matbuf_.ptr() + eltvec; \
	left_vec = cur_vec; \
	left_rowelt = rowelt; \
	left_rowvec = rowvec; \
	MOVE_VEC_PTR_LEFT(left_vec, left_rowvec, left_rowelt); \
	up_vec = cur_vec; \
	up_rowelt = rowelt; \
	up_rowvec = rowvec; \
	MOVE_VEC_PTR_UP(up_vec, up_rowvec, up_rowelt); \
	upleft_vec = up_vec; \
	upleft_rowelt = up_rowelt; \
	upleft_rowvec = up_rowvec; \
	MOVE_VEC_PTR_LEFT(upleft_vec, upleft_rowvec, upleft_rowelt); \
}

/**
 * Given the dynamic programming table and a cell, trace backwards from the
 * cell and install the edits and score/penalty in the appropriate fields of
 * res.  The RandomSource is used to break ties among equally good ways of
 * tracing back.
 *
 * Whenever we enter a cell, we check if its read/ref coordinates correspond to
 * a cell we traversed constructing a previous alignment.  If so, we backtrack
 * to the last decision point, mask out the path that led to the previously
 * observed cell, and continue along a different path.  If there are no more
 * paths to try, we stop.
 *
 * If an alignment is found, 'off' is set to the alignment's upstream-most
 * reference character's offset and true is returned.  Otherwise, false is
 * returned.
 *
 * In local alignment mode, this method is liable to be slow, especially for
 * long reads.  This is chiefly because if there is one valid solution
 * (especially if it is pretty high scoring), then many, many paths shooting
 * off that solution's path will also have valid solutions.
 */
bool SwAligner::backtraceNucleotidesLocalSseI16(
	TAlScore       escore, // in: expected score
	SwResult&      res,    // out: store results (edits and scores) here
	size_t&        off,    // out: store diagonal projection of origin
	size_t&        nbts,   // out: # backtracks
	size_t         row,    // start in this row
	size_t         col,    // start in this column
	RandomSource&  rnd)    // random gen, to choose among equal paths
{
	assert_lt(row, dpRows());
	assert_lt(col, rff_-rfi_);
	SSEData& d = fw_ ? sseI16fw_ : sseI16rc_;
	SSEMetrics& met = extend_ ? sseI16ExtendMet_ : sseI16MateMet_;
	met.bt++;
	assert(!d.profbuf_.empty());
	assert_lt(row, rd_->length());
	btnstack_.clear(); // empty the backtrack stack
	btcells_.clear();  // empty the cells-so-far list
	AlnScore score;
	score.score_ = score.gaps_ = score.ns_ = 0;
	size_t origCol = col;
	size_t gaps = 0, readGaps = 0, refGaps = 0;
	res.alres.reset();
	EList<Edit>& ned = res.alres.ned();
	assert(ned.empty());
	assert_gt(dpRows(), row);
	size_t trimEnd = dpRows() - row - 1; 
	size_t trimBeg = 0;
	size_t ct = SSEMatrix::H; // cell type
	// Row and col in terms of where they fall in the SSE vector matrix
	size_t rowelt, rowvec, eltvec;
	size_t left_rowelt, up_rowelt, upleft_rowelt;
	size_t left_rowvec, up_rowvec, upleft_rowvec;
	__m128i *cur_vec, *left_vec, *up_vec, *upleft_vec;
	const size_t gbar = sc_->gapbar;
	NEW_ROW_COL(row, col);
	// If 'backEliminate' is true, then every time we visit a cell, we remove
	// edges into the cell.  We do this to avoid some of the thrashing around
	// that occurs when there are lots of valid candidates in the same DP
	// problem.
	//const bool backEliminate = true;
	while((int)row >= 0) {
		// TODO: As soon as we enter a cell, set it as being reported through,
		// *and* mark all cells that point into this cell as being reported
		// through.  This will save us from having to consider quite so many
		// candidates.
		
		met.btcell++;
		nbts++;
		int readc = (*rd_)[rdi_ + row];
		int refm  = (int)rf_[rfi_ + col];
		int readq = (*qu_)[row];
		assert_leq(col, origCol);
		// Get score in this cell
		bool empty, reportedThru, canMoveThru, branch = false;
		int cur = SSEMatrix::H;
		if(!d.mat_.reset_[row]) {
			d.mat_.resetRow(row);
		}
		reportedThru = d.mat_.reportedThrough(row, col);
		canMoveThru = true;
		if(reportedThru) {
			canMoveThru = false;
		} else {
			empty = false;
			if(row > 0) {
				size_t rowFromEnd = d.mat_.nrow() - row - 1;
				bool gapsAllowed = !(row < gbar || rowFromEnd < gbar);
				const int floorsc = 0;
				const int offsetsc = 0x8000;
				// Move to beginning of column/row
				if(ct == SSEMatrix::E) { // AKA rdgap
					assert_gt(col, 0);
					TAlScore sc_cur = ((TCScore*)(cur_vec + SSEMatrix::E))[rowelt] + offsetsc;
					assert(gapsAllowed);
					// Currently in the E matrix; incoming transition must come from the
					// left.  It's either a gap open from the H matrix or a gap extend from
					// the E matrix.
					// TODO: save and restore origMask as well as mask
					int origMask = 0, mask = 0;
					// Get H score of cell to the left
					TAlScore sc_h_left = ((TCScore*)(left_vec + SSEMatrix::H))[left_rowelt] + offsetsc;
					if(sc_h_left > floorsc && sc_h_left - sc_->readGapOpen() == sc_cur) {
						mask |= (1 << 0); // horiz H -> E move possible
					}
					// Get E score of cell to the left
					TAlScore sc_e_left = ((TCScore*)(left_vec + SSEMatrix::E))[left_rowelt] + offsetsc;
					if(sc_e_left > floorsc && sc_e_left - sc_->readGapExtend() == sc_cur) {
						mask |= (1 << 1); // horiz E -> E move possible
					}
					origMask = mask;
					assert(origMask > 0 || sc_cur <= sc_->match());
					if(d.mat_.isEMaskSet(row, col)) {
						mask = (d.mat_.masks_[row][col] >> 8) & 3;
					}
					if(mask == 3) {
						// Horiz H -> E or horiz E -> E moves possible
#if 1
						// Pick H -> E cell
						cur = SW_BT_OALL_READ_OPEN;
						d.mat_.eMaskSet(row, col, 2); // might choose E later
#else
						if(rnd.nextU2()) {
							// Pick H -> E cell
							cur = SW_BT_OALL_READ_OPEN;
							d.mat_.eMaskSet(row, col, 2); // might choose E later
						} else {
							// Pick E -> E cell
							cur = SW_BT_RDGAP_EXTEND;
							d.mat_.eMaskSet(row, col, 1); // might choose H later
						}
#endif
						branch = true;
					} else if(mask == 2) {
						// Only horiz E -> E move possible, pick it
						cur = SW_BT_RDGAP_EXTEND;
						d.mat_.eMaskSet(row, col, 0); // done
					} else if(mask == 1) {
						// I chose the H cell
						cur = SW_BT_OALL_READ_OPEN;
						d.mat_.eMaskSet(row, col, 0); // done
					} else {
						empty = true;
						// It's empty, so the only question left is whether we should be
						// allowed in terimnate in this cell.  If it's got a valid score
						// then we *shouldn't* be allowed to terminate here because that
						// means it's part of a larger alignment that was already reported.
						canMoveThru = (origMask == 0);
					}
					if(!branch) {
						// Is this where we can eliminate some incoming paths as well?
					}
					assert(!empty || !canMoveThru);
				} else if(ct == SSEMatrix::F) { // AKA rfgap
					assert_gt(row, 0);
					assert(gapsAllowed);
					TAlScore sc_h_up = ((TCScore*)(up_vec  + SSEMatrix::H))[up_rowelt] + offsetsc;
					TAlScore sc_f_up = ((TCScore*)(up_vec  + SSEMatrix::F))[up_rowelt] + offsetsc;
					TAlScore sc_cur  = ((TCScore*)(cur_vec + SSEMatrix::F))[rowelt] + offsetsc;
					// Currently in the F matrix; incoming transition must come from above.
					// It's either a gap open from the H matrix or a gap extend from the F
					// matrix.
					// TODO: save and restore origMask as well as mask
					int origMask = 0, mask = 0;
					// Get H score of cell above
					if(sc_h_up > floorsc && sc_h_up - sc_->refGapOpen() == sc_cur) {
						mask |= (1 << 0);
					}
					// Get F score of cell above
					if(sc_f_up > floorsc && sc_f_up - sc_->refGapExtend() == sc_cur) {
						mask |= (1 << 1);
					}
					origMask = mask;
					assert(origMask > 0 || sc_cur <= sc_->match());
					if(d.mat_.isFMaskSet(row, col)) {
						mask = (d.mat_.masks_[row][col] >> 11) & 3;
					}
					if(mask == 3) {
#if 1
						// I chose the H cell
						cur = SW_BT_OALL_REF_OPEN;
						d.mat_.fMaskSet(row, col, 2); // might choose E later
#else
						if(rnd.nextU2()) {
							// I chose the H cell
							cur = SW_BT_OALL_REF_OPEN;
							d.mat_.fMaskSet(row, col, 2); // might choose E later
						} else {
							// I chose the F cell
							cur = SW_BT_RFGAP_EXTEND;
							d.mat_.fMaskSet(row, col, 1); // might choose E later
						}
#endif
						branch = true;
					} else if(mask == 2) {
						// I chose the F cell
						cur = SW_BT_RFGAP_EXTEND;
						d.mat_.fMaskSet(row, col, 0); // done
					} else if(mask == 1) {
						// I chose the H cell
						cur = SW_BT_OALL_REF_OPEN;
						d.mat_.fMaskSet(row, col, 0); // done
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
					TAlScore sc_cur      = ((TCScore*)(cur_vec + SSEMatrix::H))[rowelt]    + offsetsc;
					TAlScore sc_f_up     = ((TCScore*)(up_vec  + SSEMatrix::F))[up_rowelt] + offsetsc;
					TAlScore sc_h_up     = ((TCScore*)(up_vec  + SSEMatrix::H))[up_rowelt] + offsetsc;
					TAlScore sc_h_left   = col > 0 ? (((TCScore*)(left_vec   + SSEMatrix::H))[left_rowelt]   + offsetsc) : floorsc;
					TAlScore sc_e_left   = col > 0 ? (((TCScore*)(left_vec   + SSEMatrix::E))[left_rowelt]   + offsetsc) : floorsc;
					TAlScore sc_h_upleft = col > 0 ? (((TCScore*)(upleft_vec + SSEMatrix::H))[upleft_rowelt] + offsetsc) : floorsc;
					TAlScore sc_diag     = sc_->score(readc, refm, readq - 33);
					// TODO: save and restore origMask as well as mask
					int origMask = 0, mask = 0;
					if(gapsAllowed) {
						if(sc_h_up     > floorsc && sc_cur == sc_h_up   - sc_->refGapOpen()) {
							mask |= (1 << 0);
						}
						if(sc_h_left   > floorsc && sc_cur == sc_h_left - sc_->readGapOpen()) {
							mask |= (1 << 1);
						}
						if(sc_f_up     > floorsc && sc_cur == sc_f_up   - sc_->refGapExtend()) {
							mask |= (1 << 2);
						}
						if(sc_e_left   > floorsc && sc_cur == sc_e_left - sc_->readGapExtend()) {
							mask |= (1 << 3);
						}
					}
					if(sc_h_upleft > floorsc && sc_cur == sc_h_upleft + sc_diag) {
						mask |= (1 << 4); // diagonal is 
					}
					origMask = mask;
					assert(origMask > 0 || sc_cur <= sc_->match());
					if(d.mat_.isHMaskSet(row, col)) {
						mask = (d.mat_.masks_[row][col] >> 2) & 31;
					}
					assert(gapsAllowed || mask == (1 << 4) || mask == 0);
					int opts = alts5[mask];
					int select = -1;
					if(opts == 1) {
						select = firsts5[mask];
						assert_geq(mask, 0);
						d.mat_.hMaskSet(row, col, 0);
					} else if(opts > 1) {
#if 1
						if(       (mask & 16) != 0) {
							select = 4; // H diag
						} else if((mask & 1) != 0) {
							select = 0; // H up
						} else if((mask & 4) != 0) {
							select = 2; // F up
						} else if((mask & 2) != 0) {
							select = 1; // H left
						} else if((mask & 8) != 0) {
							select = 3; // E left
						}
#else
						select = randFromMask(rnd, mask);
#endif
						assert_geq(mask, 0);
						mask &= ~(1 << select);
						assert(gapsAllowed || mask == (1 << 4) || mask == 0);
						d.mat_.hMaskSet(row, col, mask);
						branch = true;
					} else { /* No way to backtrack! */ }
					if(select != -1) {
						if(select == 4) {
							cur = SW_BT_OALL_DIAG;
						} else if(select == 0) {
							cur = SW_BT_OALL_REF_OPEN;
						} else if(select == 1) {
							cur = SW_BT_OALL_READ_OPEN;
						} else if(select == 2) {
							cur = SW_BT_RFGAP_EXTEND;
						} else {
							assert_eq(3, select)
							cur = SW_BT_RDGAP_EXTEND;
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
			} // if(row > 0)
		} // else clause of if(reportedThru)
		if(!reportedThru) {
			d.mat_.setReportedThrough(row, col);
		}
		assert(d.mat_.reportedThrough(row, col));
		//if(backEliminate && row < d.mat_.nrow()-1) {
		//	// Possibly pick off neighbors below and to the right if the
		//	// neighbor's only way of backtracking is through this cell.
		//}
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
				NEW_ROW_COL(row, col);
				continue;
			} else {
				// No branch points to revisit; just give up
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
				int readC = (*rd_)[row];
				int refNmask = (int)rf_[rfi_+col];
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
					int64_t bonus = sc_->match(30);
					score.score_ += bonus;
					assert(!sc_->monotone || score.score() >= escore);
				}
				if(m == -1) {
					score.ns_++;
				}
				row--; col--;
				MOVE_ALL_UPLEFT();
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
					"ACGTN"[(int)(*rd_)[row]],
					EDIT_TYPE_REF_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				row--;
				ct = SSEMatrix::H;
				int pen = sc_->refGapOpen();
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; refGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				MOVE_ALL_UP();
				break;
			}
			// Move up.  Add an edit encoding the ref gap.
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
				ct = SSEMatrix::F;
				int pen = sc_->refGapExtend();
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; refGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				MOVE_ALL_UP();
				break;
			}
			case SW_BT_OALL_READ_OPEN:
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
				ct = SSEMatrix::H;
				int pen = sc_->readGapOpen();
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; readGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				MOVE_ALL_LEFT();
				break;
			}
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
				ct = SSEMatrix::E;
				int pen = sc_->readGapExtend();
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; readGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				MOVE_ALL_LEFT();
				break;
			}
			default: throw 1;
		}
	} // while((int)row > 0)
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
	assert_lt(col + rfi_, rff_);
	score.gaps_ = gaps;
	res.alres.setScore(score);
	res.alres.setShape(
		refidx_,                  // ref id
		off + rfi_ + rect_->refl, // 0-based ref offset
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
