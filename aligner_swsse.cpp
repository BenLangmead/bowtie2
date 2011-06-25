//
//  aligner_swsse_mat.cpp
//  bowtie2
//
//  Created by Benjamin Langmead on 6/2/11.
//  Copyright 2011 Johns Hopkins University. All rights reserved.
//

#ifndef NO_SSE

#include "aligner_sw_common.h"
#include "aligner_swsse.h"

/**
 * Given a number of rows (nrow), a number of columns (ncol), and the
 * number of words to fit inside a single __m128i vector, initialize the
 * matrix buffer to accomodate the needed configuration of vectors.
 */
void SSEMatrix::init(
	size_t nrow,
	size_t ncol,
	size_t wperv)
{
	nrow_ = nrow;
	ncol_ = ncol;
	wperv_ = wperv;
	nvecPerCol_ = (nrow + (wperv-1)) / wperv;
	// The +1 is so that we don't have to special-case the final column;
	// instead, we just write off the end of the useful part of the table
	// with pvEStore.
	buf_.resize((ncol+1) * nvecPerCell_ * nvecPerCol_ + 16);
	//bzero(buf_.ptr(), sizeof(__m128i) * ((ncol+1) * nvecPerCell_ * nvecPerCol_ + 16));
	// Get a 16-byte aligned pointer toward the beginning of the buffer.
	size_t aligned = ((size_t)buf_.ptr() + 15) & ~(0x0f);
	// Set up pointers into the buffer for fw query
	bufal_ = reinterpret_cast<__m128i*>(aligned);
	assert(wperv_ == 8 || wperv_ == 16);
	vecshift_ = (wperv_ == 8) ? 3 : 4;
	nvecrow_ = (nrow + (wperv_-1)) >> vecshift_;
	nveccol_ = ncol;
	colstride_ = nvecPerCol_ * nvecPerCell_;
	rowstride_ = nvecPerCell_;
	inited_ = true;
}

/**
 * Initialize the matrix of masks and backtracking flags.
 */
void SSEMatrix::initMasks() {
	assert_gt(nrow_, 0);
	assert_gt(ncol_, 0);
	masks_.resize(nrow_ * ncol_);
	bzero(masks_.ptr(), sizeof(uint16_t) * nrow_ * ncol_);
}

/**
 * Given a row, col and matrix (i.e. E, F or H), return the corresponding
 * element.
 */
int SSEMatrix::eltSlow(size_t row, size_t col, size_t mat) const {
	assert_lt(row, nrow_);
	assert_lt(col, ncol_);
	assert_leq(mat, 3);
	// Move to beginning of column/row
	size_t rowelt = row / nvecrow_;
	size_t rowvec = row % nvecrow_;
	size_t eltvec = (col * colstride_) + (rowvec * rowstride_) + mat;
	if(wperv_ == 16) {
		return (int)((uint8_t*)&bufal_[eltvec])[rowelt];
	} else {
		assert_eq(8, wperv_);
		return (int)((int16_t*)&bufal_[eltvec])[rowelt];
	}
}

#endif
