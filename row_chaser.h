/*
 * row_chaser.h
 */

#ifndef ROW_CHASER_H_
#define ROW_CHASER_H_

#include <iostream>
#include "ebwt.h"
#include "aligner_metrics.h"

/**
 * A class that statefully converts a row index to a reference
 * location.  There is a large memory-latency penalty usually
 * associated with calling the Ebwt object's mapLF method, which this
 * object does repeatedly in order to resolve the reference offset.
 * The "statefulness" in how the computation is organized here allows
 * some or all of that penalty to be hidden using prefetching.
 */
class RowChaser {

	typedef std::pair<uint32_t,uint32_t> U32Pair;
	typedef Ebwt TEbwt;

public:
	RowChaser(AlignerMetrics *metrics = NULL) :
		done(false),
		prepped_(false),
		ebwt_(NULL),
		qlen_(0),
		eh_(NULL),
		row_(0xffffffff),
		jumps_(0),
		sideloc_(),
		off_(0xffffffff),
		tlen_(0),
		metrics_(metrics)
	{ }

	/**
	 * Convert a row to a joined reference offset.  This has to be
	 * converted to understand where it is w/r/t the reference hit and
	 * offset within it.
	 */
	static uint32_t toFlatRefOff(const TEbwt* ebwt, uint32_t qlen, uint32_t row) {
		RowChaser rc;
		rc.setRow(row, qlen, ebwt);
		while(!rc.done) {
			rc.advance();
		}
		return rc.flatOff();
	}

	/**
	 * Convert a row to a reference offset.
	 */
	static U32Pair toRefOff(const TEbwt* ebwt, uint32_t qlen, uint32_t row) {
		RowChaser rc;
		rc.setRow(row, qlen, ebwt);
		while(!rc.done) {
			rc.advance();
		}
		return rc.off();
	}

	/**
	 * Set the next row for us to "chase" (i.e. map to a reference
	 * location using the BWT step-left operation).
	 */
	void setRow(uint32_t row, uint32_t qlen, const TEbwt* ebwt) {
		assert_neq(0xffffffff, row);
		assert_gt(qlen, 0);
		assert(ebwt != NULL);
		ebwt_ = ebwt;
		eh_ = &ebwt->_eh;
		row_ = row;
		qlen_ = qlen;
		ASSERT_ONLY(sideloc_.invalidate());
		if(row_ == ebwt_->_zOff) {
			// We arrived at the extreme left-hand end of the reference
			off_ = 0;
			done = true;
			return;
		} else if((row_ & eh_->_offMask) == row_) {
			// We arrived at a marked row
			off_ = ebwt_->offs()[row_ >> eh_->_offRate];
			done = true;
			return;
		}
		done = false;
		jumps_ = 0;
		off_ = 0xffffffff;
		prepped_ = false;
		prep();
	}

	/**
	 * Advance the step-left process by one step.  Check if we're done.
	 */
	void advance() {
		// Advance by 1
		assert(!done);
		while(!done) {
			assert(prepped_);
			prepped_ = false;
			if(metrics_ != NULL) metrics_->curBwtOps_++;
			uint32_t newrow = ebwt_->mapLF(sideloc_);
			ASSERT_ONLY(sideloc_.invalidate());
			jumps_++;
			assert_neq(newrow, row_);
			// Update row_ field
			row_ = newrow;
			if(row_ == ebwt_->_zOff) {
				// We arrived at the extreme left-hand end of the reference
				off_ = jumps_;
				done = true;
			} else if((row_ & eh_->_offMask) == row_) {
				// We arrived at a marked row
				off_ = ebwt_->offs()[row_ >> eh_->_offRate] + jumps_;
				done = true;
			}
			prep();
		}
	}

	/**
	 * Prepare for the next call to advance() by prefetching the
	 * appropriate portions of the index.  The caller should make sure
	 * that the
	 */
	void prep() {
		if(!done) {
			assert(!prepped_);
			assert(!sideloc_.valid());
			assert_leq(row_, eh_->_len);
			sideloc_.initFromRow(row_, *eh_, (const uint8_t*)ebwt_->ebwt());
			assert(sideloc_.valid());
		}
		prepped_ = true;
	}

	/**
	 * Get the calculated offset.  This has to be converted with a call
	 * to Ebwt::joinedToTextOff() to understand where it is w/r/t the
	 * reference hit and offset within it.
	 */
	uint32_t flatOff() const {
		return off_;
	}

	/**
	 * Get the calculated offset.
	 */
	U32Pair off() {
		uint32_t off = flatOff();
		assert_neq(0xffffffff, off);
		uint32_t tidx;
		uint32_t textoff = 0xffffffff;
		ebwt_->joinedToTextOff(qlen_, off, tidx, textoff, tlen_);
		// Note: tidx may be 0xffffffff, if alignment overlaps a
		// reference boundary
		return make_pair(tidx, textoff);
	}

	uint32_t tlen() const {
		return tlen_;
	}

	bool done;               /// true = chase is done & answer is in off_
	bool prepped_; /// true = prefetch is issued and it's OK to call advance()

protected:

	const TEbwt* ebwt_;      /// index to resolve row in
	uint32_t qlen_;          /// length of read; needed to convert to ref. coordinates
	const EbwtParams* eh_;   /// eh field from index
	uint32_t row_;           /// current row
	uint32_t jumps_;         /// # steps so far
	SideLocus sideloc_;      /// current side locus
	uint32_t off_;           /// calculated offset (0xffffffff if not done)
	uint32_t tlen_;          /// hit text length
	AlignerMetrics *metrics_;
};

#endif /* ROW_CHASER_H_ */
