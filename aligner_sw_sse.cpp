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

#ifndef NO_SSE

#include "aligner_sw.h"

/**
 * Build query profile look up tables for the read.  The query profile look
 * up table is organized as a 1D array indexed by [i][j] where i is the
 * reference character in the current DP column (0=A, 1=C, etc), and j is
 * the segment of the query we're currently working on.
 */
void SwAligner::buildQueryProfileSseU8(bool fw) {
	const BTDnaString* rd = fw ? rdfw_ : rdrc_;
	const BTString* qu = fw ? qufw_ : qurc_;
	const size_t len = rd->length();
	const size_t seglen = (len + 15) / 16;
	// How many __m128i's are needed for the SSEU8
	size_t n128s =
		64 +                    // slack bytes, for alignment?
		(seglen * ALPHA_SIZE) + // query profile data
		(seglen * 3);           // pvHStore, pvHLoad, pvE
	assert_gt(n128s, 0);
	SSEData& d = fw ? sseU8fw_ : sseU8rc_;
	d.buf_.resize(n128s * sizeof(__m128i));
	// Get a 16-byte aligned pointer toward the beginning of the buffer.
	size_t aligned = ((size_t)d.buf_.ptr() + 15) & ~(0x0f);
	// Set up pointers into the buffer for fw query
	d.qprof_    = reinterpret_cast<__m128i*>(aligned);
	d.pvHStore_ = d.qprof_    + (seglen * ALPHA_SIZE);
	d.pvHLoad_  = d.pvHStore_ + seglen;
	d.pvE_      = d.pvHLoad_  + seglen;
	d.maxPen_   = d.maxBonus_ = 0;
	d.lastIter_ = d.lastWord_ = 0;
	uint8_t *qprofBytes = reinterpret_cast<uint8_t*>(d.qprof_);
	// For each reference character A, C, G, T, N ...
	for(size_t refc = 0; refc < ALPHA_SIZE; refc++) {
		size_t h = 0;
		// For each segment ...
		for(size_t i = 0; i < seglen; i++) {
			size_t j = i;
			// For each sub-word (byte) ...
			for(size_t k = 0; k < 16; k++) {
				int sc = 0;
				if(j < len) {
					int readc = (*rd)[j];
					int readq = (*qu)[j];
					sc = sc_->score(readc, (int)(1 << refc), readq);
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
				qprofBytes[refc * (seglen*16) + h] = (uint8_t)sc;
				h++;
				j += seglen;
			}
		}
	}
}

/**
 * Build query profile look up tables for the read.  The query profile look
 * up table is organized as a 1D array indexed by [i][j] where i is the
 * reference character in the current DP column (0=A, 1=C, etc), and j is
 * the segment of the query we're currently working on.
 */
void SwAligner::buildQueryProfileSseI16(bool fw) {
	const BTDnaString* rd = fw ? rdfw_ : rdrc_;
	const BTString* qu = fw ? qufw_ : qurc_;
	const size_t len = rd->length();
	const size_t seglen = (len + 7) / 8;
	// How many __m128i's are needed for the SSEI16
	size_t n128s =
		64 +                    // slack bytes, for alignment?
		(seglen * ALPHA_SIZE) + // query profile data
		(seglen * 3);           // pvHStore, pvHLoad, pvE
	assert_gt(n128s, 0);
	SSEData& d = fw ? sseI16fw_ : sseI16rc_;
	d.buf_.resize(n128s * sizeof(__m128i));
	// Get a 16-byte aligned pointer toward the beginning of the buffer.
	size_t aligned = ((size_t)d.buf_.ptr() + 15) & ~(0x0f);
	// Set up pointers into the buffer for fw query
	d.qprof_    = reinterpret_cast<__m128i*>(aligned);
	d.pvHStore_ = d.qprof_    + (seglen * ALPHA_SIZE);
	d.pvHLoad_  = d.pvHStore_ + seglen;
	d.pvE_      = d.pvHLoad_  + seglen;
	d.maxPen_   = d.maxBonus_ = 0;
	d.lastIter_ = d.lastWord_ = 0;
	int16_t *qprofBytes = reinterpret_cast<int16_t*>(d.qprof_);
	// For each reference character A, C, G, T, N ...
	for(size_t refc = 0; refc < ALPHA_SIZE; refc++) {
		size_t h = 0;
		// For each segment ...
		for(size_t i = 0; i < seglen; i++) {
			size_t j = i;
			for(size_t k = 0; k < 8; k++) {
				int sc = 0;
				if(j < len) {
					int readc = (*rd)[j];
					int readq = (*qu)[j];
					sc = sc_->score(readc, (int)(1 << refc), readq - 33);
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
				qprofBytes[refc * (seglen*8) + h] = (int16_t)sc;
				h++;
				j += seglen;
			}
		}
	}
}

/**
 *
 */
uint8_t SwAligner::alignNucleotidesSseU8() {
	return 0;
}

/**
 * Solve the current alignment problem using SSE instructions that operate on 8
 * signed 16-bit values packed into a single 128-bit register.
 */
int16_t SwAligner::alignNucleotidesSseI16() {
	const size_t NWORDS_PER_REG  = 8;
	const size_t NBYTES_PER_REG  = 16;
	const size_t NBITS_PER_WORD  = 16;
	const size_t NBYTES_PER_WORD = 2;

	assert_leq(rdf_, rd_->length());
	assert_leq(rdf_, qu_->length());
	assert_lt(rfi_, rff_);
	assert_lt(rdi_, rdf_);
	assert_eq(rd_->length(), qu_->length());
	assert_geq(sc_->gapbar, 1);
	assert(en_ == NULL || en_->size() == width_);
	assert(st_ == NULL || st_->size() == width_);
	assert(repOk());
#ifndef NDEBUG
	for(size_t i = rfi_; i < rff_; i++) {
		assert_range(0, 16, (int)rf_[i]);
	}
#endif

	SSEData& d = fw_ ? sseI16fw_ : sseI16rc_;
	assert(!d.buf_.empty());
	assert(d.qprof_ != NULL);
	assert(d.pvHStore_ != NULL);
	assert(d.pvHLoad_ != NULL);
	assert(d.pvE_ != NULL);

	bool monotone = (d.maxBonus_ == 0);
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
	__m128i vhilsw   = _mm_setzero_si128();
	__m128i vfilllsw = _mm_setzero_si128();
	__m128i vmax     = _mm_setzero_si128();
	__m128i ve       = _mm_setzero_si128();
	__m128i vf       = _mm_setzero_si128();
	__m128i vh       = _mm_setzero_si128();
	__m128i vtmp     = _mm_setzero_si128();

	assert_gt(sc_->refGapOpen(), 0);
	rfgapo = _mm_insert_epi16(rfgapo, sc_->refGapOpen(), 0);
	rfgapo = _mm_shufflelo_epi16(rfgapo, 0);
	rfgapo = _mm_shuffle_epi32(rfgapo, 0);
	
	// Set all elts to reference gap extension penalty
	assert_gt(sc_->refGapExtend(), 0);
	assert_leq(sc_->refGapExtend(), sc_->refGapOpen());
	rfgape = _mm_insert_epi16(rfgape, sc_->refGapExtend(), 0);
	rfgape = _mm_shufflelo_epi16(rfgape, 0);
	rfgape = _mm_shuffle_epi32(rfgape, 0);

	// Set all elts to read gap open penalty
	assert_gt(sc_->readGapOpen(), 0);
	rdgapo = _mm_insert_epi16(rdgapo, sc_->readGapOpen(), 0);
	rdgapo = _mm_shufflelo_epi16(rdgapo, 0);
	rdgapo = _mm_shuffle_epi32(rdgapo, 0);
	
	// Set all elts to read gap extension penalty
	assert_gt(sc_->readGapExtend(), 0);
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

	// vhilsw: topmost (least sig) word set to 0x7fff, all other words=0
	vhilsw = _mm_shuffle_epi32(vhi, 0);
	vhilsw = _mm_srli_si128(vhilsw, NBYTES_PER_REG - NBYTES_PER_WORD);
	
	if(monotone) {
		vfilllsw = vhilsw;
	} else {
		vfilllsw = vlolsw;
	}
	
	__m128i *pvScore;                // points into the query profile
	__m128i *pvtmp;                  // for swapping pvHStore/pvHLoad
	__m128i *pvE      = d.pvE_;      // E vectors across outer loop iters
	__m128i *pvHStore = d.pvHStore_; // H vectors across outer loop iters
	__m128i *pvHLoad  = d.pvHLoad_;  // H vectors across outer loop iters
	__m128i *qprof    = d.qprof_;    // query profile
	
	// Maximum score in final row
	int lrmax = std::numeric_limits<short>::min();
	
	// Initialize the vectors that carry H and E values over outer-loop iters
	if(monotone) {
		for(size_t i = 0; i < iter; i++) {
			_mm_store_si128(pvE + i, vlo);
			_mm_store_si128(pvHStore + i, vhi); // start high
		}
	} else {
		for(size_t i = 0; i < iter; i++) {
			_mm_store_si128(pvE + i, vlo);
			_mm_store_si128(pvHStore + i, vlo); // start low
		}
	}
	
	for(size_t i = rfi_; i < rff_; i++) {
		// Calculate distance from RHS of paralellogram; helps us decide
		// whether a solution can end in this column 
		size_t fromend = rff_ - i - 1;
		
		// Fetch the appropriate query profile.  Note that elements of rf_ must
		// be numbers, not masks.
		pvScore = qprof + ((size_t)firsts5[(int)rf_[i]] * iter);
		
		// Set all elts to 0x8000 (min value for signed 16-bit)
		vf = _mm_cmpeq_epi16(vf, vf);
		vf = _mm_slli_epi16(vf, NBITS_PER_WORD-1);

		// Load the next h value
		vh = _mm_load_si128(pvHStore + iter - 1);
		// Shift 2 bytes down so that topmost (least sig) word gets 0
		vh = _mm_slli_si128(vh, NBYTES_PER_WORD);
		// Now set topmost word to vlolsw
		if(i - rfi_ >= width_ || (st_ != NULL && !(*st_)[i - rfi_])) {
			// Fill with large negative value
			vh = _mm_or_si128(vh, vlolsw);
		} else {
			// Fill with normal initial value
			vh = _mm_or_si128(vh, vfilllsw);
		}

		// Swap load and store arrays
		pvtmp = pvHLoad;
		pvHLoad = pvHStore;
		pvHStore = pvtmp;
		
		size_t j;
		for(j = 0; j < iter; j++) {
			// Load value ve from previous col
			ve = _mm_load_si128(pvE + j);

			// Factor in score from next query profile elt
			vh = _mm_adds_epi16(vh, *pvScore++);

			// Update highest score encountered this far
			vmax = _mm_max_epi16(vmax, vh);

			// Update H, factoring in E and F
			vh = _mm_max_epi16(vh, ve);
			vh = _mm_max_epi16(vh, vf);

			// Save the new vH values
			_mm_store_si128 (pvHStore + j, vh);

			// Slight modification from the swsse2 code here to allow separate
			// read gap and reference gap penalties

			// Update vE value
			vtmp = vh;
			vh = _mm_subs_epi16(vh, rdgapo);
			ve = _mm_subs_epi16(ve, rdgape);
			ve = _mm_max_epi16(ve, vh);

			// Update vf value
			vtmp = _mm_subs_epi16(vtmp, rfgapo);
			vf = _mm_subs_epi16 (vf, rfgape);
			vf = _mm_max_epi16 (vf, vtmp);

			// Save E values
			_mm_store_si128(pvE + j, ve);

			// Load the next h value
			vh = _mm_load_si128(pvHLoad + j);
		}
		
		// Reset pointers to the start of the saved data
		j = 0;
		vh = _mm_load_si128(pvHStore + j);
		
		// Now I have the topmost vh and the bottommost vf

		// vf from last row gets shifted down by one to overlay the first row
		vf = _mm_slli_si128(vf, NBYTES_PER_WORD);
		vf = _mm_or_si128(vf, vlolsw);
		
		// Charge a gap open to vh - this is one way we might update vf
		vtmp = _mm_subs_epi16(vh, rfgapo);
		vtmp = _mm_cmpgt_epi16(vf, vtmp);
		int cmp = _mm_movemask_epi8(vtmp);
		
		// If any element of vtmp is greater than H - gap-open...
		while (cmp != 0x0000) {
		
			// Grab saved E vector
			// What's going on here?
			ve = _mm_load_si128(pvE + j);
			vh = _mm_max_epi16(vh, vf);

			// Save vH values
			_mm_store_si128 (pvHStore + j, vh);

			// Update E in case it can be improved using our new vh
			vtmp = _mm_subs_epi16(vh, rdgapo);
			ve = _mm_max_epi16(ve, vtmp);
			_mm_store_si128(pvE + j, ve);

			/* update vf value */
			vf = _mm_subs_epi16(vf, rfgape);

			j++;
			if(j >= iter) {
				j = 0;
				vf = _mm_slli_si128(vf, NBYTES_PER_WORD);
				// ? - not in pseudocode
				vf = _mm_or_si128(vf, vlolsw);
			}

			vh = _mm_load_si128(pvHStore + j);
			vtmp = _mm_subs_epi16(vh, rfgapo);
			vtmp = _mm_cmpgt_epi16(vf, vtmp);
			cmp  = _mm_movemask_epi8(vtmp);
		}
		
		if(monotone && fromend < width_) {
			if(en_ == NULL || (*en_)[width_ - fromend - 1]) {
				// Update lrmax
				vh = _mm_load_si128(pvHStore + d.lastIter_);
				// Note: we may not want to extract from the final row
				int16_t lr = ((int16_t*)(&vh))[d.lastWord_];
				if(lr > lrmax) {
					lrmax = lr;
				}
			}
		}
	}

	if(monotone) {
		if(lrmax == std::numeric_limits<short>::min()) {
			return (int16_t)0x8000;
		}
		return (int16_t)lrmax - 0x7fff;
	} else {
		// Find largest score in vmax
		vtmp = _mm_srli_si128(vmax, 8);
		vmax = _mm_max_epi16(vmax, vtmp);
		vtmp = _mm_srli_si128(vmax, 4);
		vmax = _mm_max_epi16(vmax, vtmp);
		vtmp = _mm_srli_si128(vmax, 2);
		vmax = _mm_max_epi16(vmax, vtmp);

		// Return largest score
		return (int16_t)_mm_extract_epi16(vmax, 0) + 0x8000;
	}
}

#endif
