/*
 * aligner_23mm.h
 */

#ifndef ALIGNER_23MM_H_
#define ALIGNER_23MM_H_

#include <utility>
#include "aligner.h"
#include "hit.h"
#include "row_chaser.h"
#include "range_chaser.h"
#include "ref_aligner.h"
#include "bowtie1_range_source.h"
#include "ds.h"

/**
 * Concrete factory for constructing unpaired 2- or 3-mismatch aligners.
 */
class Unpaired23mmAlignerV1Factory : public AlignerFactory {
	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
	typedef EList<TRangeSrcDr*> TRangeSrcDrPtrVec;
public:
	Unpaired23mmAlignerV1Factory(
			Ebwt& ebwtFw,
			Ebwt* ebwtBw,
			bool two,
			HitSink& sink,
			const HitSinkPerThreadFactory& sinkPtFactory,
			RangeCache *cacheFw,
			RangeCache *cacheBw,
			uint32_t cacheLimit,
			ChunkPool *pool,
			BitPairReference* refs,
			EList<String<Dna5> >& os,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
			ebwtBw_(ebwtBw),
			two_(two),
			sink_(sink),
			sinkPtFactory_(sinkPtFactory),
			cacheFw_(cacheFw),
			cacheBw_(cacheBw),
			cacheLimit_(cacheLimit),
			pool_(pool),
			refs_(refs),
			os_(os),
			seed_(seed)
	{
		assert(ebwtFw.isInMemory());
		assert(ebwtBw != NULL);
		assert(ebwtBw->isInMemory());
	}

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {

		HitSinkPerThread* sinkPt = sinkPtFactory_.create();
		EbwtSearchParams* params =
			new EbwtSearchParams(*sinkPt, os_);

		const bool seeded = false;

		EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
			 ebwtBw_, true, 0xffffffff, true, 0, seeded);
		EbwtRangeSource *rFw_Fw = new EbwtRangeSource(
			&ebwtFw_, true, 0xffffffff, false, 0, seeded);
		EbwtRangeSource *rFw_BwHalf = new EbwtRangeSource(
			 ebwtBw_, true, 0xffffffff, false, 2,  seeded);
		EbwtRangeSource *rFw_FwHalf = NULL;
		if(!two_) {
			rFw_FwHalf = new EbwtRangeSource(
				&ebwtFw_, true, 0xffffffff, false, 3,  seeded);
		}

		// Driver wrapper for rFw_Bw
		EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
			*params, rFw_Bw, true, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, true, pool_, NULL);
		// Driver wrapper for rFw_Fw
		EbwtRangeSourceDriver * drFw_Fw = new EbwtRangeSourceDriver(
			*params, rFw_Fw, true, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, true, pool_, NULL);
		// Driver wrapper for rFw_Fw
		EbwtRangeSourceDriver * drFw_BwHalf = new EbwtRangeSourceDriver(
			*params, rFw_BwHalf, true, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_BEGINNING,    // nothing's unrevisitable
			PIN_TO_HI_HALF_EDGE,
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, true, pool_, NULL);
		// Driver wrapper for rFw_Fw
		EbwtRangeSourceDriver * drFw_FwHalf = NULL;
		if(!two_) {
			drFw_FwHalf = new EbwtRangeSourceDriver(
				*params, rFw_FwHalf, true, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, true, pool_, NULL);
		}

		TRangeSrcDrPtrVec *drVec = new TRangeSrcDrPtrVec();
		if(!gNofw) {
			drVec->push_back(drFw_Bw);
			drVec->push_back(drFw_Fw);
			drVec->push_back(drFw_BwHalf);
			if(!two_) {
				drVec->push_back(drFw_FwHalf);
			}
		}

		EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
			&ebwtFw_, false, 0xffffffff, true, 0, seeded);
		EbwtRangeSource *rRc_Bw = new EbwtRangeSource(
			 ebwtBw_, false, 0xffffffff, false, 0, seeded);
		EbwtRangeSource *rRc_FwHalf = new EbwtRangeSource(
			&ebwtFw_, false, 0xffffffff, false, 2,  seeded);
		EbwtRangeSource *rRc_BwHalf = NULL;
		if(!two_) {
			rRc_BwHalf = new EbwtRangeSource(
				 ebwtBw_, false, 0xffffffff, false, 3,  seeded);
		}

		// Driver wrapper for rRc_Fw
		EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
			*params, rRc_Fw, false, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, true, pool_, NULL);
		// Driver wrapper for rRc_Bw
		EbwtRangeSourceDriver * drRc_Bw = new EbwtRangeSourceDriver(
			*params, rRc_Bw, false, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, true, pool_, NULL);
		// Driver wrapper for rRc_Fw
		EbwtRangeSourceDriver * drRc_FwHalf = new EbwtRangeSourceDriver(
			*params, rRc_FwHalf, false, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_BEGINNING,    // nothing's unrevisitable
			PIN_TO_HI_HALF_EDGE,
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, true, pool_, NULL);
		EbwtRangeSourceDriver * drRc_BwHalf = NULL;
		if(!two_) {
			drRc_BwHalf = new EbwtRangeSourceDriver(
				*params, rRc_BwHalf, false, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, true, pool_, NULL);
		}
		if(!gNorc) {
			drVec->push_back(drRc_Fw);
			drVec->push_back(drRc_Bw);
			drVec->push_back(drRc_FwHalf);
			if(!two_) {
				drVec->push_back(drRc_BwHalf);
			}
		}
		TCostAwareRangeSrcDr* dr = new TCostAwareRangeSrcDr(drVec, false);
		delete drVec;

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(cacheLimit_, cacheFw_, cacheBw_);

		ColorspaceDecoder *dec = new ColorspaceDecoder();

		return new UnpairedAlignerV2<EbwtRangeSource>(
			params, dr, rchase,
			sink_, sinkPtFactory_, sinkPt, os_, refs_,
			INT_MAX, pool_, dec, NULL, NULL);
	}

private:
	Ebwt& ebwtFw_;
	Ebwt* ebwtBw_;
	bool two_;
	HitSink& sink_;
	const HitSinkPerThreadFactory& sinkPtFactory_;
	RangeCache *cacheFw_;
	RangeCache *cacheBw_;
	const uint32_t cacheLimit_;
	ChunkPool *pool_;
	BitPairReference* refs_;
	EList<String<Dna5> >& os_;
	uint32_t seed_;
};

/**
 * Concrete factory for constructing paired 2- or 3-mismatch aligners.
 */
class Paired23mmAlignerV1Factory : public AlignerFactory {
	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef ListRangeSourceDriver<EbwtRangeSource> TListRangeSrcDr;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
	typedef EList<TRangeSrcDr*> TRangeSrcDrPtrVec;
public:
	Paired23mmAlignerV1Factory(
			Ebwt& ebwtFw,
			Ebwt* ebwtBw,
			bool two,
			HitSink& sink,
			const HitSinkPerThreadFactory& sinkPtFactory,
			bool dontReconcile,
			uint32_t symCeil,
			uint32_t mixedThresh,
			uint32_t mixedAttemptLim,
			RangeCache *cacheFw,
			RangeCache *cacheBw,
			uint32_t cacheLimit,
			ChunkPool *pool,
			BitPairReference* refs,
			EList<String<Dna5> >& os,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
			ebwtBw_(ebwtBw),
			two_(two),
			sink_(sink),
			sinkPtFactory_(sinkPtFactory),
			dontReconcile_(dontReconcile),
			symCeil_(symCeil),
			mixedThresh_(mixedThresh),
			mixedAttemptLim_(mixedAttemptLim),
			cacheFw_(cacheFw),
			cacheBw_(cacheBw),
			cacheLimit_(cacheLimit),
			pool_(pool),
			refs_(refs), os_(os),
			seed_(seed)
	{
		assert(ebwtBw != NULL);
		assert(ebwtFw.isInMemory());
		assert(ebwtBw->isInMemory());
	}

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {
		HitSinkPerThread* sinkPt = sinkPtFactory_.createMult(2);
		HitSinkPerThread* sinkPtSe1 = NULL, * sinkPtSe2 = NULL;
		EbwtSearchParams* params = new EbwtSearchParams(*sinkPt, os_);
		EbwtSearchParams* paramsSe1 = NULL, * paramsSe2 = NULL;
		if(gReportSe) {
			sinkPtSe1 = sinkPtFactory_.create();
			sinkPtSe2 = sinkPtFactory_.create();
			paramsSe1 =
				new EbwtSearchParams(*sinkPtSe1, os_);
			paramsSe2 =
				new EbwtSearchParams(*sinkPtSe2, os_);
		}

		const bool seeded = false;

		bool do1Fw = true;
		bool do1Rc = true;
		bool do2Fw = true;
		bool do2Rc = true;
		if(gNofw) {
			if(gMate1fw) do1Fw = false;
			else        do1Rc = false;
			if(gMate2fw) do2Fw = false;
			else        do2Rc = false;
		}
		if(gNorc) {
			if(gMate1fw) do1Rc = false;
			else        do1Fw = false;
			if(gMate2fw) do2Rc = false;
			else        do2Fw = false;
		}

		TRangeSrcDrPtrVec *dr1FwVec = new TRangeSrcDrPtrVec();

		if(do1Fw) {
			EbwtRangeSource *r1Fw_Bw = new EbwtRangeSource(
				 ebwtBw_, true, 0xffffffff, true,  0, seeded);
			EbwtRangeSource *r1Fw_Fw = new EbwtRangeSource(
				&ebwtFw_, true, 0xffffffff, false, 0, seeded);
			EbwtRangeSource *r1Fw_BwHalf = new EbwtRangeSource(
				 ebwtBw_, true, 0xffffffff, false, 2, seeded);
			EbwtRangeSource *r1Fw_FwHalf = two_ ? NULL : new EbwtRangeSource(
				&ebwtFw_, true, 0xffffffff, false, 3, seeded);

			// Driver wrapper for rFw_Bw
			EbwtRangeSourceDriver * dr1Fw_Bw = new EbwtRangeSourceDriver(
				*params, r1Fw_Bw, true, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, true, pool_, NULL);
			// Driver wrapper for rFw_Fw
			EbwtRangeSourceDriver * dr1Fw_Fw = new EbwtRangeSourceDriver(
				*params, r1Fw_Fw, true, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, true, pool_, NULL);
			// Driver wrapper for rFw_Fw
			EbwtRangeSourceDriver * dr1Fw_BwHalf = new EbwtRangeSourceDriver(
				*params, r1Fw_BwHalf, true, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, true, pool_, NULL);
			dr1FwVec->push_back(dr1Fw_Bw);
			dr1FwVec->push_back(dr1Fw_Fw);
			dr1FwVec->push_back(dr1Fw_BwHalf);
			if(!two_) {
				// Driver wrapper for rFw_Fw
				EbwtRangeSourceDriver * dr1Fw_FwHalf = two_ ? NULL : new EbwtRangeSourceDriver(
					*params, r1Fw_FwHalf, true, false, sink_, sinkPt,
					0,          // seedLen (0 = whole read is seed)
					false,       // nudgeLeft (true for Fw index, false for Bw)
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_BEGINNING,
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_LEN,
					os_, true, pool_, NULL);
				dr1FwVec->push_back(dr1Fw_FwHalf);
			}
		}

		TRangeSrcDrPtrVec *dr1RcVec;
		dr1RcVec = dr1FwVec;

		if(do1Rc) {
			EbwtRangeSource *r1Rc_Fw = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, true,  0, seeded);
			EbwtRangeSource *r1Rc_Bw = new EbwtRangeSource(
				 ebwtBw_, false, 0xffffffff, false, 0, seeded);
			EbwtRangeSource *r1Rc_FwHalf = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, false, 2, seeded);
			EbwtRangeSource *r1Rc_BwHalf = two_ ? NULL : new EbwtRangeSource(
				 ebwtBw_, false, 0xffffffff, false, 3, seeded);

			// Driver wrapper for rRc_Fw
			EbwtRangeSourceDriver * dr1Rc_Fw = new EbwtRangeSourceDriver(
				*params, r1Rc_Fw, false, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, true, pool_, NULL);
			// Driver wrapper for rRc_Bw
			EbwtRangeSourceDriver * dr1Rc_Bw = new EbwtRangeSourceDriver(
				*params, r1Rc_Bw, false, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, true, pool_, NULL);
			// Driver wrapper for rRc_Fw
			EbwtRangeSourceDriver * dr1Rc_FwHalf = new EbwtRangeSourceDriver(
				*params, r1Rc_FwHalf, false, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_BEGINNING,
				PIN_TO_HI_HALF_EDGE,
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, true, pool_, NULL);
			dr1RcVec->push_back(dr1Rc_Fw);
			dr1RcVec->push_back(dr1Rc_Bw);
			dr1RcVec->push_back(dr1Rc_FwHalf);
			if(!two_) {
				// Driver wrapper for rRc_Bw
				EbwtRangeSourceDriver * dr1Rc_BwHalf = two_ ? NULL : new EbwtRangeSourceDriver(
					*params, r1Rc_BwHalf, false, false, sink_, sinkPt,
					0,          // seedLen (0 = whole read is seed)
					false,       // nudgeLeft (true for Fw index, false for Bw)
					PIN_TO_BEGINNING,
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_LEN,
					os_, true, pool_, NULL);
				dr1RcVec->push_back(dr1Rc_BwHalf);
			}
		}

		TRangeSrcDrPtrVec *dr2FwVec;
		dr2FwVec = dr1FwVec;

		if(do2Fw) {
			EbwtRangeSource *r2Fw_Bw = new EbwtRangeSource(
				 ebwtBw_, true, 0xffffffff, true,  0, seeded);
			EbwtRangeSource *r2Fw_Fw = new EbwtRangeSource(
				&ebwtFw_, true, 0xffffffff, false, 0, seeded);
			EbwtRangeSource *r2Fw_BwHalf = new EbwtRangeSource(
				 ebwtBw_, true, 0xffffffff, false, 2, seeded);
			EbwtRangeSource *r2Fw_FwHalf = two_ ? NULL : new EbwtRangeSource(
				&ebwtFw_, true, 0xffffffff, false, 3, seeded);

			// Driver wrapper for rFw_Bw
			EbwtRangeSourceDriver * dr2Fw_Bw = new EbwtRangeSourceDriver(
				*params, r2Fw_Bw, true, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, false, pool_, NULL);
			// Driver wrapper for rFw_Fw
			EbwtRangeSourceDriver * dr2Fw_Fw = new EbwtRangeSourceDriver(
				*params, r2Fw_Fw, true, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, false, pool_, NULL);
			// Driver wrapper for rFw_Fw
			EbwtRangeSourceDriver * dr2Fw_BwHalf = new EbwtRangeSourceDriver(
				*params, r2Fw_BwHalf, true, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, false, pool_, NULL);
			dr2FwVec->push_back(dr2Fw_Bw);
			dr2FwVec->push_back(dr2Fw_Fw);
			dr2FwVec->push_back(dr2Fw_BwHalf);
			if(!two_) {
				EbwtRangeSourceDriver * dr2Fw_FwHalf = two_ ? NULL : new EbwtRangeSourceDriver(
					*params, r2Fw_FwHalf, true, false, sink_, sinkPt,
					0,          // seedLen (0 = whole read is seed)
					false,       // nudgeLeft (true for Fw index, false for Bw)
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_BEGINNING,
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_LEN,
					os_, false, pool_, NULL);
				dr2FwVec->push_back(dr2Fw_FwHalf);
			}
		}

		TRangeSrcDrPtrVec *dr2RcVec;
		dr2RcVec = dr1FwVec;

		if(do2Rc) {
			EbwtRangeSource *r2Rc_Fw = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, true,  0, seeded);
			EbwtRangeSource *r2Rc_Bw = new EbwtRangeSource(
				 ebwtBw_, false, 0xffffffff, false, 0, seeded);
			EbwtRangeSource *r2Rc_FwHalf = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, false, 2,  seeded);
			EbwtRangeSource *r2Rc_BwHalf = two_ ? NULL : new EbwtRangeSource(
				 ebwtBw_, false, 0xffffffff, false, 3,  seeded);

			// Driver wrapper for rRc_Fw
			EbwtRangeSourceDriver * dr2Rc_Fw = new EbwtRangeSourceDriver(
				*params, r2Rc_Fw, false, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, false, pool_, NULL);
			// Driver wrapper for rRc_Bw
			EbwtRangeSourceDriver * dr2Rc_Bw = new EbwtRangeSourceDriver(
				*params, r2Rc_Bw, false, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, false, pool_, NULL);
			// Driver wrapper for rRc_Fw
			EbwtRangeSourceDriver * dr2Rc_FwHalf = new EbwtRangeSourceDriver(
				*params, r2Rc_FwHalf, false, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, false, pool_, NULL);
			dr2RcVec->push_back(dr2Rc_Fw);
			dr2RcVec->push_back(dr2Rc_Bw);
			dr2RcVec->push_back(dr2Rc_FwHalf);
			if(!two_) {
				EbwtRangeSourceDriver * dr2Rc_BwHalf = two_ ? NULL : new EbwtRangeSourceDriver(
					*params, r2Rc_BwHalf, false, false, sink_, sinkPt,
					0,          // seedLen (0 = whole read is seed)
					false,       // nudgeLeft (true for Fw index, false for Bw)
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_BEGINNING,
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_LEN,
					os_, false, pool_, NULL);
				dr2RcVec->push_back(dr2Rc_BwHalf);
			}
		}

		RefAligner<String<Dna5> >* refAligner;
		if(two_) {
			refAligner = new TwoMMRefAligner<String<Dna5> >();
		} else {
			refAligner = new ThreeMMRefAligner<String<Dna5> >();
		}

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(cacheLimit_, cacheFw_, cacheBw_);

		ColorspaceDecoder *dec = new ColorspaceDecoder();

		PairedBWAlignerV2<EbwtRangeSource>* al = new PairedBWAlignerV2<EbwtRangeSource>(
			params, paramsSe1, paramsSe2,
			new TCostAwareRangeSrcDr(dr1FwVec, true),
			refAligner, rchase,
			sink_, sinkPtFactory_,
			sinkPt, sinkPtSe1, sinkPtSe2,
			mixedAttemptLim_, refs_,
			INT_MAX, pool_, dec, NULL);
		delete dr1FwVec;
		return al;
	}

private:
	Ebwt& ebwtFw_;
	Ebwt* ebwtBw_;
	bool two_;
	HitSink& sink_;
	const HitSinkPerThreadFactory& sinkPtFactory_;
	const bool dontReconcile_;
	const uint32_t symCeil_;
	const uint32_t mixedThresh_;
	const uint32_t mixedAttemptLim_;
	RangeCache *cacheFw_;
	RangeCache *cacheBw_;
	const uint32_t cacheLimit_;
	ChunkPool *pool_;
	BitPairReference* refs_;
	EList<String<Dna5> >& os_;
	const uint32_t seed_;
};

#endif /* ALIGNER_23MM_H_ */
