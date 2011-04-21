/**
 * range.h
 */

#ifndef RANGE_H_
#define RANGE_H_

#include <stdint.h>
#include "ds.h"
#include "sstring.h"

/**
 * A Burrows-Wheeler range, along with information about the set of
 * alignments it represents.
 */
class Range {

public:

	Range() :
		top_(0xffffffff), bot_(0), cost_(0), stratum_(0), rlen_(0),
		qlen_(0), fw_(true), mate1_(true), edits_(128), ebwt_(NULL)
	{ }

	void init(uint32_t tp, uint32_t bt, uint16_t cst,
	          uint32_t strtm, uint32_t qln, bool f,
	          const EList<Edit>& edts,
	          const Ebwt *eb)
	{
		top_ = tp;
		bot_ = bt;
		cost_ = cst;
		stratum_ = strtm;
		qlen_ = qln;
		fw_ = f;
		mate1_ = true;
		edits_ = edts; // copy
		ebwt_ = eb;
		calcRlen();
		assert(Range::finalized(edits_));
		assert(repOk());
		assert(valid());
	}

	bool valid() const { return top_ < 0xffffffff; }
	void invalidate()  { top_ = 0xffffffff; }

	/**
	 * Clear relevant state and set to be invalid.
	 */
	void clear() {
		edits_.clear();
		invalidate();
	}

	/**
	 * Check if the Range is finalized.  Assert if not.
	 */
	static bool finalized(const EList<Edit>& edits) {
		for(size_t i = 1; i < edits.size(); i++) {
			assert_geq(edits[i].pos, edits[i-1].pos);
			if(edits[i].pos == edits[i-1].pos) {
				assert_eq(EDIT_TYPE_READ_GAP, edits[i-1].type);
			}
		}
		return true;
	}

	/**
	 * Set whether this range is mate #1 (or unpaired).
	 */
	void setMate1(bool mate1) { mate1_ = mate1; }

#ifndef NDEBUG
	/**
	 * Check that Range is internally consistent.
	 */
	bool repOk() const {
		assert_leq(stratum_, (int)edits_.size());
		assert_gt(bot_, top_);
		assert_neq(0xffff, cost_);
		assert_lt(stratum_, 4);
		assert(Range::finalized(edits_));
		// Non-zero cost implies at least one edit
		if(cost_ != 0 || stratum_ != 0) {
			assert_gt(edits_.size(), 0)
		}
		// Reference char can't be same as read char
		for(int i = 0; i < (int)edits_.size(); i++) {
			const Edit& e = edits_[i];
			assert_in(e.chr, "ACGT-");
			assert_in(e.qchr, "ACGTN-");
			assert_neq(e.chr, e.qchr);
		}
		return true;
	}
#endif

	const EList<Edit>& edits() const { return edits_; }
	uint16_t cost() const { return cost_; }
	uint32_t top() const { return top_; }
	uint32_t bot() const { return bot_; }
	int stratum() const { return stratum_; }
	uint32_t rlen() const { return rlen_; }
	uint32_t qlen() const { return qlen_; }
	bool fw() const { return fw_; }
	bool mate1() const { return mate1_; }
	const Ebwt* ebwt() const { return ebwt_; }

protected:

	uint32_t top_;     // top of range
	uint32_t bot_;     // bottom of range
	uint16_t cost_;    // cost
	int stratum_;      // stratum
	uint32_t rlen_;    // length of reference side of alignment
	uint32_t qlen_;    // length of query side of alignment
	bool fw_;          // the forward orientation of read aligned?
	bool mate1_;       // read aligned is #1 mate/single?
	EList<Edit> edits_; // edits
	const Ebwt *ebwt_;

	/**
	 * Calculate the length of the reference side of the alignment
	 * given the length of the query side and the list of all edits.
	 */
	void calcRlen() {
		assert_gt(qlen_, 0);
		assert(Range::finalized(edits_));
		rlen_ = qlen_;
		for(size_t i = 0; i < edits_.size(); i++) {
			if(edits_[i].isReadGap()) rlen_++;
			else if(edits_[i].isRefGap()) rlen_--;
		}
		assert_gt(rlen_, 0);
	}
};

#endif /* RANGE_H_ */
