//
//  aligner_swsse_mat.cpp
//  bowtie2
//
//  Created by Benjamin Langmead on 6/2/11.
//  Copyright 2011 Johns Hopkins University. All rights reserved.
//

#include "aligner_swsse.h"

/**
 * Given a row, col and matrix (i.e. E, F or H), return the corresponding
 * element.
 */
int SSEMatrix::eltSlow(size_t row, size_t col, size_t mat) const {
	assert(inited_);
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
