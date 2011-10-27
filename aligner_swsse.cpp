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

#include <string.h>
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
	matbuf_.resizeNoCopy((ncol+1) * nvecPerCell_ * nvecPerCol_);
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
	masks_.resize(nrow_);
	reset_.resizeNoCopy(nrow_);
	reset_.fill(false);
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
		return (int)((uint8_t*)(matbuf_.ptr() + eltvec))[rowelt];
	} else {
		assert_eq(8, wperv_);
		return (int)((int16_t*)(matbuf_.ptr() + eltvec))[rowelt];
	}
}
