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

#include <string.h>
#include "aligner_sw_common.h"
#include "aligner_swsse.h"

/**
 * For each nonzero 5-bit H-cell backtrace mask, the move the u8 backtraces
 * take from the cell (an SW_BT_* value, bits 0-2) and the mask they leave
 * in the cell (bits 3-7).  Mask bits: 0 H above (ref gap open), 1 H left
 * (read gap open), 2 F above (ref gap extend), 3 E left (read gap extend), 4
 * H diagonal.  A mask with one bit set takes that move and leaves 0 (alts5,
 * firsts5); one with several takes the first of diagonal, H above, F above,
 * H left and E left, and leaves the others.  btHSelectBit is the same choice
 * made the parent's way, for the asserts.
 */
const uint8_t btHSelect[32] = {
	0,                                 // 0: no move (not used)
	SW_BT_OALL_REF_OPEN  | ( 0 << 3), //  1
	SW_BT_OALL_READ_OPEN | ( 0 << 3), //  2
	SW_BT_OALL_REF_OPEN  | ( 2 << 3), //  3
	SW_BT_RFGAP_EXTEND   | ( 0 << 3), //  4
	SW_BT_OALL_REF_OPEN  | ( 4 << 3), //  5
	SW_BT_RFGAP_EXTEND   | ( 2 << 3), //  6
	SW_BT_OALL_REF_OPEN  | ( 6 << 3), //  7
	SW_BT_RDGAP_EXTEND   | ( 0 << 3), //  8
	SW_BT_OALL_REF_OPEN  | ( 8 << 3), //  9
	SW_BT_OALL_READ_OPEN | ( 8 << 3), // 10
	SW_BT_OALL_REF_OPEN  | (10 << 3), // 11
	SW_BT_RFGAP_EXTEND   | ( 8 << 3), // 12
	SW_BT_OALL_REF_OPEN  | (12 << 3), // 13
	SW_BT_RFGAP_EXTEND   | (10 << 3), // 14
	SW_BT_OALL_REF_OPEN  | (14 << 3), // 15
	SW_BT_OALL_DIAG      | ( 0 << 3), // 16
	SW_BT_OALL_DIAG      | ( 1 << 3), // 17
	SW_BT_OALL_DIAG      | ( 2 << 3), // 18
	SW_BT_OALL_DIAG      | ( 3 << 3), // 19
	SW_BT_OALL_DIAG      | ( 4 << 3), // 20
	SW_BT_OALL_DIAG      | ( 5 << 3), // 21
	SW_BT_OALL_DIAG      | ( 6 << 3), // 22
	SW_BT_OALL_DIAG      | ( 7 << 3), // 23
	SW_BT_OALL_DIAG      | ( 8 << 3), // 24
	SW_BT_OALL_DIAG      | ( 9 << 3), // 25
	SW_BT_OALL_DIAG      | (10 << 3), // 26
	SW_BT_OALL_DIAG      | (11 << 3), // 27
	SW_BT_OALL_DIAG      | (12 << 3), // 28
	SW_BT_OALL_DIAG      | (13 << 3), // 29
	SW_BT_OALL_DIAG      | (14 << 3), // 30
	SW_BT_OALL_DIAG      | (15 << 3)  // 31
};

/**
 * Given a number of rows (nrow), a number of columns (ncol), and the
 * number of words to fit inside a single SSERegI vector, initialize the
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
	try {
		matbuf_.resizeNoCopy((ncol+1) * nvecPerCell_ * nvecPerCol_);
	} catch(exception& e) {
		cerr << "Tried to allocate DP matrix with " << (ncol+1)
		     << " columns, " << nvecPerCol_
			 << " vectors per column, and and " << nvecPerCell_
			 << " vectors per cell" << endl;
		throw e;
	}
	{
		// compute vecshift_=log2(wperv)
		size_t wp = wperv;
		assert(wp == (NBYTES_PER_REG/2) || wp == NBYTES_PER_REG);
	        vecshift_ = 0;
        	while( wp>>=1 ) vecshift_++;
	}
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
	// The dirty intervals only grow, so every row masks_ has ever exposed,
	// and so every row that may hold a buffer, keeps its record.  A row
	// without one yet gets [0, SIZE_MAX): resetRow then clears all of it.
	size_t ndirty = dirtylo_.size();
	if(ndirty < nrow_) {
		dirtylo_.resize(nrow_);
		dirtyhi_.resize(nrow_);
		for(size_t i = ndirty; i < nrow_; i++) {
			dirtylo_[i] = 0;
			dirtyhi_[i] = std::numeric_limits<size_t>::max();
		}
	}
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
	if(wperv_ == NBYTES_PER_REG) {
		return (int)((uint8_t*)(matbuf_.ptr() + eltvec))[rowelt];
	} else {
		assert_eq((NBYTES_PER_REG/2), wperv_);
		return (int)((int16_t*)(matbuf_.ptr() + eltvec))[rowelt];
	}
}
