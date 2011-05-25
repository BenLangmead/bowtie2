/**
 * aligner_sw_sse.cpp
 *
 * Versions of key alignment functions that use vector instructions to
 * accelerate dynamic programming.  Based chiefly on the striped Smith-Waterman
 * paper and implementation by Michael Farrar.  See:
 *
 * Farrar M. Striped Smith-Waterman speeds database searches six times over
 * other SIMD implementations. Bioinformatics. 2007 Jan 15;23(2):156-61.
 *
 * http://sites.google.com/site/farrarmichael/smith-waterman
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
	uint8_t *qprofBytes = reinterpret_cast<uint8_t*>(d.qprof_);
	// For each reference character A, C, G, T, N ...
	for(size_t refc = 0; refc < ALPHA_SIZE; refc++) {
		size_t h = 0;
		// For each segment ...
		for(size_t i = 0; i < seglen; i++) {
			size_t j = i;
			for(size_t k = 0; k < 16; k++) {
				int sc = 0;
				if(j < len) {
					// Set query profile for fw query
					int readc = (*rd)[j];
					if(readc > 3 || (int)refc > 3) {
						sc = sc_->n(30);
					} else if(readc != (int)refc) {
						sc = sc_->mm(30);
					} else {
						// leave it at 0
					}
				}
				qprofBytes[refc * seglen + h] = (uint8_t)sc;
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
	uint8_t *qprofBytes = reinterpret_cast<uint8_t*>(d.qprof_);
	// For each reference character A, C, G, T, N ...
	for(size_t refc = 0; refc < ALPHA_SIZE; refc++) {
		size_t h = 0;
		// For each segment ...
		for(size_t i = 0; i < seglen; i++) {
			size_t j = i;
			for(size_t k = 0; k < 8; k++) {
				int sc = 0;
				if(j < len) {
					// Set query profile for fw query
					int readc = (*rd)[j];
					if(readc > 3 || (int)refc > 3) {
						sc = sc_->n(30);
					} else if(readc != (int)refc) {
						sc = sc_->mm(30);
					} else {
						// leave it at 0
					}
				}
				qprofBytes[refc * seglen + h] = (uint8_t)sc;
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
 *
 */
int16_t SwAligner::alignNucleotidesSseI16() {
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

	SSEData& d = fw_ ? sseI16fw_ : sseI16rc_;
	assert(!d.buf_.empty());
	assert(d.qprof_ != NULL);
	assert(d.pvHStore_ != NULL);
	assert(d.pvHLoad_ != NULL);
	assert(d.pvE_ != NULL);

	size_t iter = (dpRows() + 7) / 8; // iter = segLen

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
	__m128i vmax     = _mm_setzero_si128();
	__m128i vminimus = _mm_setzero_si128();
	__m128i vmin     = _mm_setzero_si128();
	__m128i ve       = _mm_setzero_si128();
	__m128i vf       = _mm_setzero_si128();
	__m128i vh       = _mm_setzero_si128();
	__m128i vtmp     = _mm_setzero_si128();

	assert_geq(sc_->rfGapConst, 0);
	rfgapo = _mm_insert_epi16(rfgapo, sc_->rfGapConst, 0);
	rfgapo = _mm_shufflelo_epi16(rfgapo, 0);
	rfgapo = _mm_shuffle_epi32(rfgapo, 0);
	
	// Set all elts to reference gap extension penalty
	assert_gt(sc_->rfGapLinear, 0);
	rfgape = _mm_insert_epi16(rfgape, sc_->rfGapLinear, 0);
	rfgape = _mm_shufflelo_epi16(rfgape, 0);
	rfgape = _mm_shuffle_epi32(rfgape, 0);

	// Set all elts to read gap open penalty
	assert_geq(sc_->rdGapConst, 0);
	rdgapo = _mm_insert_epi16(rdgapo, sc_->rdGapConst, 0);
	rdgapo = _mm_shufflelo_epi16(rdgapo, 0);
	rdgapo = _mm_shuffle_epi32(rdgapo, 0);
	
	// Set all elts to read gap extension penalty
	assert_gt(sc_->rdGapLinear, 0);
	rdgape = _mm_insert_epi16(rdgape, sc_->rdGapLinear, 0);
	rdgape = _mm_shufflelo_epi16(rdgape, 0);
	rdgape = _mm_shuffle_epi32(rdgape, 0);
	
	// Set all elts to 0x8000 (min value for signed 16-bit)
	vmax = _mm_cmpeq_epi16(vmax, vmax);    // all elts = 0xffff
	vmax = _mm_slli_epi16(vmax, 15);       // all elts = 0x8000
	
	// Set all elts to 0x8000 (min value for signed 16-bit)
	vminimus = _mm_shuffle_epi32(vmax, 0); // all elts = 0x8000

	// Set all elts to 0x0002
	vmin = _mm_shuffle_epi32(vmax, 0);     // all elts = 0x8000
	vmin = _mm_srli_si128(vmin, 14);       // all elts = 0x0002
	
	__m128i *pvScore;      // points into the query profile
	__m128i *pvtmp;        // for swapping pvHStore/pvHLoad
	__m128i *pvE      = d.pvE_;      // E vectors across outer loop iters
	__m128i *pvHStore = d.pvHStore_; // H vectors across outer loop iters
	__m128i *pvHLoad  = d.pvHLoad_;  // H vectors across outer loop iters
	__m128i *qprof    = d.qprof_;    // query profile
	
	// Initialize the vectors that carry H and E values over outer-loop iters
	for(size_t i = 0; i < iter; i++) {
		_mm_store_si128(pvE + i, vmax);
		_mm_store_si128(pvHStore + i, vmax);
	}
	
	for(size_t i = rfi_; i < rff_; i++) {
		// Fetch the appropriate query profile
		pvScore = qprof + ((size_t)rf_[i] * iter);
		
		// Set all elts to 0x8000 (min value for signed 16-bit)
		vf = _mm_cmpeq_epi16(vf, vf);
		vf = _mm_slli_epi16(vf, 15);

		// Load the next h value
		vh = _mm_load_si128(pvHStore + iter - 1);
		vh = _mm_slli_si128(vh, 2);
		vh = _mm_or_si128(vh, vmin);

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
			vh = _mm_max_epi16 (vh, ve);
			vh = _mm_max_epi16 (vh, vf);

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
			vf = _mm_max_epi16 (vf, vh);

			// Save E values
			_mm_store_si128(pvE + j, ve);

			// Load the next h value
			vh = _mm_load_si128(pvHLoad + j);
		}
		
		// Reset pointers to the start of the saved data
		j = 0;
		vh = _mm_load_si128(pvHStore + j);

		/*  the computed vf value is for the given column.  since */
		/*  we are at the end, we need to shift the vf value over */
		/*  to the next column. */
		vf = _mm_slli_si128(vf, 2);
		vf = _mm_or_si128(vf, vmin);
		vtmp = _mm_subs_epi16(vh, rfgapo);
		vtmp = _mm_cmpgt_epi16(vf, vtmp);
		int cmp = _mm_movemask_epi8(vtmp);
		
		// If any element of vtmp is greater than H - gap-open...
		while (cmp != 0x0000) {
		
			// Grab saved E vector
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
				
				vf = _mm_slli_si128(vf, 2);
				// ? - not in pseudocode
				vf = _mm_or_si128(vf, vmin);
			}

			vh = _mm_load_si128(pvHStore + j);
			vtmp = _mm_subs_epi16(vh, rfgapo);
			vtmp = _mm_cmpgt_epi16(vf, vtmp);
			cmp  = _mm_movemask_epi8(vtmp);
		}
	}

	// Find largest score in vmax
	vtmp = _mm_srli_si128(vmax, 8);
	vmax = _mm_max_epi16(vmax, vtmp);
	vtmp = _mm_srli_si128(vmax, 4);
	vmax = _mm_max_epi16(vmax, vtmp);
	vtmp = _mm_srli_si128(vmax, 2);
	vmax = _mm_max_epi16(vmax, vtmp);

	// store in temporary variable
	int16_t score = (int16_t)_mm_extract_epi16(vmax, 0);

	/* return largest score */
	int SHORT_BIAS = 0;
	return score + SHORT_BIAS;
}

#endif
