/*
 *  aligner_sw_col.h
 */

#ifndef ALIGNER_SW_COL_H_
#define ALIGNER_SW_COL_H_

#include <stdint.h>
#include "ds.h"
#include "threading.h"
#include "aligner_seed.h"
#include "reference.h"
#include "random_source.h"
#include "mem_ids.h"
#include "aligner_result.h"

/**
 * A bitmask encoding which backtracking paths out of a particular cell
 * correspond to optimal subpaths.  This struct is tailored to the
 * colorspace case where each transition down in the dynamic
 * programming matrix is also associated with a nucleotide assignment
 * to the source (upper) row.
 */
struct SwColorCellMask {

	/**
	 * Set all flags to 0, indicating there is no way to backtrack from
	 * this cell to an optimal answer.
	 */
	void clear() {
		*((uint16_t*)this) = 0;
	}

	/**
	 * Return true iff there are no backward paths recorded in this
	 * mask.
	 */
	inline bool empty() const {
		return *((uint16_t*)this) == 0;
	}

	/**
	 * Return true iff it's possible to extend a gap in the reference
	 * in the cell below this one.
	 */
	inline bool refExtendPossible() const {
		return rfop || rfex;
	}

	/**
	 * Return true iff it's possible to open a gap in the reference
	 * in the cell below this one (false implies that only extension
	 * is possible).
	 */
	inline bool refOpenPossible() const {
		return diag || rfop || rfex;
	}

	/**
	 * Return true iff it's possible to extend a gap in the read
	 * in the cell to the right of this one.
	 */
	inline bool readExtendPossible() const {
		return rdop || rdex;
	}

	/**
	 * Return true iff it's possible to open a gap in the read in the
	 * cell to the right of this one (false implies that only extension
	 * is possible).
	 */
	inline bool readOpenPossible() const {
		return diag || rdop || rdex;
	}
	
	/**
	 * Return true iff there is >0 possible way to backtrack from this
	 * cell.
	 */
	inline int numPossible() const {
		int num = 0;
		num += mask2popcnt[diag];
		num += mask2popcnt[rfop];
		num += mask2popcnt[rfex];
		num += rdop;
		num += rdex;
		return num;
	}

	/**
	 * Select a path for backtracking from this cell.  If there is a
	 * tie among eligible paths, break it randomly.  Return value is
	 * a pair where first = a flag indicating the backtrack type (see
	 * enum defining SW_BT_* above), and second = a selection for
	 * the read character for the next row up.  second should be
	 * ignored if the backtrack type is a gap in the read.
	 */
	std::pair<int, int> randBacktrack(RandomSource& rand);

	uint16_t diag     : 4;
	uint16_t rfop     : 4;
	uint16_t rfex     : 4;
	uint16_t rdop     : 1;
	uint16_t rdex     : 1;
	uint16_t reserved : 2;
};

/**
 * Encapsulates all information needed to encode the optimal subproblem
 * at a cell in a colorspace SW matrix.
 */
struct SwColorCell {

	inline void updateHoriz(
		const SwColorCell& lc,
		int refMask,
		const Penalties& pen,
		size_t nceil,
		int penceil);

	inline void updateDiag(
		const SwColorCell& uc,
		int refMask,
		int prevColor, // color b/t this row, one above
		int prevQual,  // quality of color
		const Penalties& pens,
		size_t nceil,
		int penceil);   // penalty ceiling

	inline void updateVert(
		const SwColorCell& uc,
		int prevColor, // color b/t this row, one above
		int prevQual,  // quality of color
		const Penalties& pen,
		size_t nceil,
		int penceil);

	/**
	 * Clear this cell so that it's ready for updates.
	 */
	void clear() {
		// Initially, best scores are all invalid
		best[0] = best[1] = best[2] = best[3] = AlignmentScore::INVALID();
		// Initially, there's no way to backtrack from this cell
		mask[0].clear();
		mask[1].clear();
		mask[2].clear();
		mask[3].clear();
		empty = true;
		ASSERT_ONLY(finalized = false);
	}
	
	/**
	 * Return true if any bests are valid.
	 */
	bool valid() const {
		for(int i = 0; i < 4; i++) {
			if(VALID_AL_SCORE(best[i])) {
				assert(!mask[i].empty());
				return true;
			} else {
				assert(mask[i].empty());
			}
		}
		return false;
	}
	
	/**
	 * Caller supplies a current-best and current-second best score and
	 * we update them according to the incoming scores for this cell.
	 * If there is a tie, take the one with more gaps; this is just so
	 * that we can sanity-check the backtrack by rejecting if it has
	 * more gaps.  Likewise for Ns.
	 */
	bool updateBest(AlignmentScore& bestSc, int& c, int penceil) const {
		assert(finalized);
		bool ret = false;
		for(int i = 0; i < 4; i++) {
			assert_leq(abs(best[i].score()), penceil);
			if(best[i] > bestSc) {
				bestSc = best[i];
				ret = true;
				c = i;
			}
		}
		return ret;
	}
	
	/**
	 * We finished updating the cell; set empty and finalized
	 * appropriately.
	 */
	bool finalize(int penceil) {
		ASSERT_ONLY(finalized = true);
		for(int i = 0; i < 4; i++) {
			if(!mask[i].empty()) {
				assert(VALID_AL_SCORE(best[i]));
				assert_leq(abs(best[i].score()), penceil);
				empty = false;
#ifdef NDEBUG
				break;
#endif
			}
		}
		return !empty;
	}

	// Best incoming score for each 'to' character
	AlignmentScore best[4];
	// Mask for tied-for-best incoming paths for each 'to' character
	SwColorCellMask mask[4];
	
	bool empty;
	ASSERT_ONLY(bool finalized);
};

#endif /*ndef ALIGNER_SW_COL_H_*/
