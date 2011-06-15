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
#include "mask.h"

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

/**
 * Return true iff the H mask has been set with a previous call to hMaskSet.
 */
inline bool SSEMatrix::isHMaskSet(
	size_t row,          // current row
	size_t col) const    // current column
{
	return ((masks_[row * ncol_ + col] & (1 << 1)) != 0);
}

/**
 * Set the given cell's H mask.  This is the mask of remaining legal ways to
 * backtrack from the H cell at this coordinate.  It's 5 bits long and has
 * offset=2 into the 16-bit field.
 */
inline void SSEMatrix::hMaskSet(
	size_t row,          // current row
	size_t col,          // current column
	int mask)
{
	assert_lt(mask, 32);
	masks_[row * ncol_ + col] &= ~(31 << 1);
	masks_[row * ncol_ + col] |= (1 << 1 | mask << 2);
}

/**
 * Return true iff the E mask has been set with a previous call to eMaskSet.
 */
inline bool SSEMatrix::isEMaskSet(
	size_t row,          // current row
	size_t col) const    // current column
{
	return ((masks_[row * ncol_ + col] & (1 << 7)) != 0);
}

/**
 * Set the given cell's E mask.  This is the mask of remaining legal ways to
 * backtrack from the E cell at this coordinate.  It's 2 bits long and has
 * offset=8 into the 16-bit field.
 */
inline void SSEMatrix::eMaskSet(
	size_t row,          // current row
	size_t col,          // current column
	int mask)
{
	assert_lt(mask, 4);
	masks_[row * ncol_ + col] &= ~(7 << 7);
	masks_[row * ncol_ + col] |=  (1 << 7 | mask << 8);
}

/**
 * Return true iff the F mask has been set with a previous call to fMaskSet.
 */
inline bool SSEMatrix::isFMaskSet(
	size_t row,          // current row
	size_t col) const    // current column
{
	return ((masks_[row * ncol_ + col] & (1 << 10)) != 0);
}

/**
 * Set the given cell's F mask.  This is the mask of remaining legal ways to
 * backtrack from the F cell at this coordinate.  It's 2 bits long and has
 * offset=11 into the 16-bit field.
 */
inline void SSEMatrix::fMaskSet(
	size_t row,          // current row
	size_t col,          // current column
	int mask)
{
	assert_lt(mask, 4);
	masks_[row * ncol_ + col] &= ~(7 << 10);
	masks_[row * ncol_ + col] |=  (1 << 10 | mask << 11);
}

// The various ways that one might backtrack from a later cell (either oall,
// rdgap or rfgap) to an earlier cell
// enum {
//	SW_BT_OALL_DIAG,         // from oall cell to oall cell
//	SW_BT_OALL_REF_OPEN,     // from oall cell to oall cell
//	SW_BT_OALL_READ_OPEN,    // from oall cell to oall cell
//	SW_BT_RDGAP_EXTEND,      // from rdgap cell to rdgap cell
//	SW_BT_RFGAP_EXTEND       // from rfgap cell to rfgap cell
// };

/**
 * Analyze a cell in the SSE-filled dynamic programming matrix.  Determine &
 * memorize ways that we can backtrack from the cell.  If there is at least one
 * way to backtrack, select one at random and return the selection.
 */
void SSEMatrix::analyzeCell(
	size_t row,          // current row
	size_t col,          // current column
	size_t ct,           // current cell type: E/F/H
	int refc,
	int readc,
	int readq,
	const Scoring& sc,   // scoring scheme
	TAlScore offsetsc,   // add to matrix elements to get actual score
	TAlScore floorsc,    // local-alignment score floor
	RandomSource& rand,  // rand gen for choosing among equal options
	bool& empty,         // out: =true iff no way to backtrace
	int& cur,            // out: =type of transition
	bool& branch,        // out: =true iff we chose among >1 options
	bool& canMoveThru,   // out: =true iff ...
	bool& reportedThru)  // out: =true iff ...
{
	TAlScore sc_cur = helt(row, col) + offsetsc;
	reportedThru = reportedThrough(row, col);
	canMoveThru = true;
	if(reportedThru) {
		canMoveThru = false;
		return;
	}
	empty = false;
	if(row == 0) {
		return;
	}
	assert_gt(row, 0);
	size_t rowFromEnd = nrow_ - row - 1;
	bool gapsAllowed = true;
	if(row < (size_t)sc.gapbar || rowFromEnd < (size_t)sc.gapbar) {
		gapsAllowed = false;
	}
	if(ct == E) { // AKA rdgap
		assert(gapsAllowed);
		// Currently in the E matrix; incoming transition must come from the
		// left.  It's either a gap open from the H matrix or a gap extend from
		// the E matrix.
		assert_gt(col, 0);
		// TODO: save and restore origMask as well as mask
		int origMask = 0, mask = 0;
		// Get H score of cell to the left
		TAlScore sc_h_left = helt(row, col-1) + offsetsc;
		if(sc_h_left > floorsc && sc_h_left - sc.readGapOpen() == sc_cur) {
			mask |= (1 << 0);
		}
		// Get E score of cell to the left
		TAlScore sc_e_left = eelt(row, col-1) + offsetsc;
		if(sc_e_left > floorsc && sc_e_left - sc.readGapExtend() == sc_cur) {
			mask |= (1 << 1);
		}
		origMask = mask;
		if(isEMaskSet(row, col)) {
			mask = (masks_[row * ncol_ + col] >> 8) & 3;
		}
		if(mask == 3) {
			if(rand.nextU2()) {
				// I chose the H cell
				cur = SW_BT_OALL_READ_OPEN;
				eMaskSet(row, col, 2); // might choose E later
			} else {
				// I chose the E cell
				cur = SW_BT_RDGAP_EXTEND;
				eMaskSet(row, col, 1); // might choose H later
			}
			branch = true;
		} else if(mask == 2) {
			// I chose the E cell
			cur = SW_BT_RDGAP_EXTEND;
			eMaskSet(row, col, 0); // done
		} else if(mask == 1) {
			// I chose the H cell
			cur = SW_BT_OALL_READ_OPEN;
			eMaskSet(row, col, 0); // done
		} else {
			empty = true;
			// It's empty, so the only question left is whether we should be
			// allowed in terimnate in this cell.  If it's got a valid score
			// then we *shouldn't* be allowed to terminate here because that
			// means it's part of a larger alignment that was already reported.
			canMoveThru = (origMask == 0);
		}
	} else if(ct == F) { // AKA rfgap
		assert(gapsAllowed);
		// Currently in the F matrix; incoming transition must come from above.
		// It's either a gap open from the H matrix or a gap extend from the F
		// matrix.
		assert_gt(row, 0);
		// TODO: save and restore origMask as well as mask
		int origMask = 0, mask = 0;
		// Get H score of cell above
		TAlScore sc_h_up = helt(row-1, col) + offsetsc;
		if(sc_h_up > floorsc && sc_h_up - sc.refGapOpen() == sc_cur) {
			mask |= (1 << 0);
		}
		// Get F score of cell above
		TAlScore sc_f_up = felt(row-1, col) + offsetsc;
		if(sc_f_up > floorsc && sc_f_up - sc.refGapExtend() == sc_cur) {
			mask |= (1 << 1);
		}
		origMask = mask;
		if(isFMaskSet(row, col)) {
			mask = (masks_[row * ncol_ + col] >> 11) & 3;
		}
		if(mask == 3) {
			if(rand.nextU2()) {
				// I chose the H cell
				cur = SW_BT_OALL_REF_OPEN;
				fMaskSet(row, col, 2); // might choose E later
			} else {
				// I chose the F cell
				cur = SW_BT_RFGAP_EXTEND;
				fMaskSet(row, col, 1); // might choose E later
			}
			branch = true;
		} else if(mask == 2) {
			// I chose the F cell
			cur = SW_BT_RFGAP_EXTEND;
			fMaskSet(row, col, 0); // done
		} else if(mask == 1) {
			// I chose the H cell
			cur = SW_BT_OALL_REF_OPEN;
			fMaskSet(row, col, 0); // done
		} else {
			empty = true;
			// It's empty, so the only question left is whether we should be
			// allowed in terimnate in this cell.  If it's got a valid score
			// then we *shouldn't* be allowed to terminate here because that
			// means it's part of a larger alignment that was already reported.
			canMoveThru = (origMask == 0);
		}
	} else {
		assert_eq(H, ct);
		// TODO: save and restore origMask as well as mask
		int origMask = 0, mask = 0;
		TAlScore sc_f_up     = felt(row-1, col) + offsetsc;
		TAlScore sc_h_up     = helt(row-1, col) + offsetsc;
		TAlScore sc_h_left   = col > 0 ? (helt(row, col-1) + offsetsc) : floorsc;
		TAlScore sc_e_left   = col > 0 ? (eelt(row, col-1) + offsetsc) : floorsc;
		TAlScore sc_h_upleft = col > 0 ? (helt(row-1, col-1) + offsetsc) : floorsc;
		TAlScore sc_diag     = sc.score(readc, (int)refc, readq - 33);
		if(gapsAllowed) {
			if(sc_h_up     > floorsc && sc_cur == sc_h_up   - sc.refGapOpen()) {
				mask |= (1 << 0);
			}
			if(sc_h_left   > floorsc && sc_cur == sc_h_left - sc.readGapOpen()) {
				mask |= (1 << 1);
			}
			if(sc_f_up     > floorsc && sc_cur == sc_f_up   - sc.refGapExtend()) {
				mask |= (1 << 2);
			}
			if(sc_e_left   > floorsc && sc_cur == sc_e_left - sc.readGapExtend()) {
				mask |= (1 << 3);
			}
		}
		if(sc_h_upleft > floorsc && sc_cur == sc_h_upleft + sc_diag) {
			mask |= (1 << 4);
		}
		origMask = mask;
		if(isHMaskSet(row, col)) {
			mask = (masks_[row * ncol_ + col] >> 2) & 31;
		}
		assert(gapsAllowed || mask == (1 << 4) || mask == 0);
		int opts = alts5[mask];
		int select = -1;
		if(opts == 1) {
			select = firsts5[mask];
			assert_geq(mask, 0);
			hMaskSet(row, col, 0);
		} else if(opts > 1) {
			select = randFromMask(rand, mask);
			assert_geq(mask, 0);
			mask &= ~(1 << select);
			assert(gapsAllowed || mask == (1 << 4) || mask == 0);
			hMaskSet(row, col, mask);
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
	
#if 0
	if(empty) {
		cout << "EMPTY" << endl;
	} else {
		if(cur == SW_BT_OALL_DIAG) {
			cout << "H -> diag -> H from [" << row << ", " << col << "]" << endl;
		}
		if(cur == SW_BT_OALL_REF_OPEN) {
			cout << "H -> ref open -> H from [" << row << ", " << col << "]" << endl;
		}
		if(cur == SW_BT_OALL_READ_OPEN) {
			cout << "H -> read open -> H from [" << row << ", " << col << "]" << endl;
		}
		if(cur == SW_BT_RFGAP_EXTEND) {
			cout << "F -> ref extend -> H from [" << row << ", " << col << "]" << endl;
		}
		if(cur == SW_BT_RDGAP_EXTEND) {
			cout << "E -> read extend -> H from [" << row << ", " << col << "]" << endl;
		}
		if(branch) {
			cout << "  BRANCHED" << endl;
		}
	}
#endif
}

#endif
