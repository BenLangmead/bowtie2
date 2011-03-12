/*
 * aligner_1mm.h
 */

#ifndef ALIGNER_1MM_H_
#define ALIGNER_1MM_H_

#include <utility>
#include "aligner.h"
#include "hit.h"
#include "row_chaser.h"
#include "range_chaser.h"
#include "ref_aligner.h"
#include "bowtie1_range_source.h"
#include "ds.h"

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class Unpaired1mmAlignerV1Factory : public AlignerFactory {
	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
	typedef EList<TRangeSrcDr*> TRangeSrcDrPtrVec;
public:
	Unpaired1mmAlignerV1Factory(
			Ebwt& ebwtFw,
			Ebwt* ebwtBw,
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

		const int halfAndHalf = 0;
		const bool seeded = false;

		EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
			 ebwtBw_, true, 0xffffffff, true,  halfAndHalf, seeded);
		EbwtRangeSource *rFw_Fw = new EbwtRangeSource(
			&ebwtFw_, true, 0xffffffff, false, halfAndHalf, seeded);

		EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
			*params, rFw_Bw, true, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, true, pool_, NULL);
		//
		EbwtRangeSourceDriver * drFw_Fw = new EbwtRangeSourceDriver(
			*params, rFw_Fw, true, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, true, pool_, NULL);
		TRangeSrcDrPtrVec *drVec = new TRangeSrcDrPtrVec();
		if(!gNofw) {
			drVec->push_back(drFw_Bw);
			drVec->push_back(drFw_Fw);
		}

		EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
			&ebwtFw_, false, 0xffffffff, true,  halfAndHalf, seeded);
		EbwtRangeSource *rRc_Bw = new EbwtRangeSource(
			 ebwtBw_, false, 0xffffffff, false, halfAndHalf, seeded);

		EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
			*params, rRc_Fw, false, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, true, pool_, NULL);
		//
		EbwtRangeSourceDriver * drRc_Bw = new EbwtRangeSourceDriver(
			*params, rRc_Bw, false, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, true, pool_, NULL);
		if(!gNorc) {
			drVec->push_back(drRc_Fw);
			drVec->push_back(drRc_Bw);
		}
		TCostAwareRangeSrcDr* dr = new TCostAwareRangeSrcDr(drVec, false);
		delete drVec;

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(cacheLimit_, cacheFw_, cacheBw_);

		ColorspaceDecoder *dec = new ColorspaceDecoder();

		// Set up the aligner
		return new UnpairedAlignerV2<EbwtRangeSource>(
			params, dr, rchase,
			sink_, sinkPtFactory_, sinkPt, os_, refs_,
			INT_MAX, pool_, dec, NULL, NULL);
	}

private:
	Ebwt& ebwtFw_;
	Ebwt* ebwtBw_;
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
 * Concrete factory class for constructing unpaired exact aligners.
 */
class Paired1mmAlignerV1Factory : public AlignerFactory {
	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
	typedef EList<TRangeSrcDr*> TRangeSrcDrPtrVec;
public:
	Paired1mmAlignerV1Factory(
			Ebwt& ebwtFw,
			Ebwt* ebwtBw,
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
			refs_(refs),
			os_(os),
			seed_(seed)
	{
		assert(ebwtBw != NULL);
		assert(ebwtFw.isInMemory());
		assert(ebwtBw->isInMemory());
	}

	/**EbwtSearchParams
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {
		HitSinkPerThread* sinkPt = sinkPtFactory_.createMult(2);
		HitSinkPerThread* sinkPtSe1 = NULL, * sinkPtSe2 = NULL;
		EbwtSearchParams* params =
			new EbwtSearchParams(*sinkPt, os_);
		EbwtSearchParams* paramsSe1 = NULL, * paramsSe2 = NULL;
		if(gReportSe) {
			sinkPtSe1 = sinkPtFactory_.create();
			sinkPtSe2 = sinkPtFactory_.create();
			paramsSe1 =
				new EbwtSearchParams(*sinkPtSe1, os_);
			paramsSe2 =
				new EbwtSearchParams(*sinkPtSe2, os_);
		}

		const int halfAndHalf = 0;
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

		TRangeSrcDrPtrVec *dr1FwVec;
		dr1FwVec = new TRangeSrcDrPtrVec();
		if(do1Fw) {
			EbwtRangeSource *r1Fw_Bw = new EbwtRangeSource(
				 ebwtBw_, true, 0xffffffff, true,  halfAndHalf, seeded);
			EbwtRangeSource *r1Fw_Fw = new EbwtRangeSource(
				&ebwtFw_, true, 0xffffffff, false, halfAndHalf, seeded);

			EbwtRangeSourceDriver * dr1Fw_Bw = new EbwtRangeSourceDriver(
				*params, r1Fw_Bw, true, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_LEN, // allow 1 mismatch in rest of read
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, true, pool_, NULL);
			EbwtRangeSourceDriver * dr1Fw_Fw = new EbwtRangeSourceDriver(
				*params, r1Fw_Fw, true, false, sink_, sinkPt,
				0,          // seedLen
				false,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right-hand half alignment is unrevisitable
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, true, pool_, NULL);

			dr1FwVec->push_back(dr1Fw_Bw);
			dr1FwVec->push_back(dr1Fw_Fw);
		}

		TRangeSrcDrPtrVec *dr1RcVec;
		dr1RcVec = dr1FwVec;
		if(do1Rc) {
			EbwtRangeSource *r1Rc_Fw = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, true,  halfAndHalf, seeded);
			EbwtRangeSource *r1Rc_Bw = new EbwtRangeSource(
				 ebwtBw_, false, 0xffffffff, false, halfAndHalf, seeded);

			EbwtRangeSourceDriver * dr1Rc_Fw = new EbwtRangeSourceDriver(
				*params, r1Rc_Fw, false, false, sink_, sinkPt,
				0,          // seedLen
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right-hand half alignment is unrevisitable
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, true, pool_, NULL);
			EbwtRangeSourceDriver * dr1Rc_Bw = new EbwtRangeSourceDriver(
				*params, r1Rc_Bw, false, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_LEN, // allow 1 mismatch in rest of read
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, true, pool_, NULL);
			dr1RcVec->push_back(dr1Rc_Fw);
			dr1RcVec->push_back(dr1Rc_Bw);
		}

		TRangeSrcDrPtrVec *dr2FwVec;
		dr2FwVec = dr1FwVec;
		if(do2Fw) {
			EbwtRangeSource *r2Fw_Bw = new EbwtRangeSource(
				 ebwtBw_, true, 0xffffffff, true,  halfAndHalf, seeded);
			EbwtRangeSource *r2Fw_Fw = new EbwtRangeSource(
				&ebwtFw_, true, 0xffffffff, false, halfAndHalf, seeded);

			EbwtRangeSourceDriver * dr2Fw_Bw = new EbwtRangeSourceDriver(
				*params, r2Fw_Bw, true, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_LEN, // allow 1 mismatch in rest of read
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, false, pool_, NULL);
			EbwtRangeSourceDriver * dr2Fw_Fw = new EbwtRangeSourceDriver(
				*params, r2Fw_Fw, true, false, sink_, sinkPt,
				0,          // seedLen
				false,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right-hand half alignment is unrevisitable
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, false, pool_, NULL);
			dr2FwVec->push_back(dr2Fw_Bw);
			dr2FwVec->push_back(dr2Fw_Fw);
		}

		TRangeSrcDrPtrVec *dr2RcVec;
		dr2RcVec = dr1FwVec;
		if(do2Rc) {
			EbwtRangeSource *r2Rc_Fw = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, true,  halfAndHalf, seeded);
			EbwtRangeSource *r2Rc_Bw = new EbwtRangeSource(
				 ebwtBw_, false, 0xffffffff, false, halfAndHalf, seeded);

			EbwtRangeSourceDriver * dr2Rc_Fw = new EbwtRangeSourceDriver(
				*params, r2Rc_Fw, false, false, sink_, sinkPt,
				0,          // seedLen
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right-hand half alignment is unrevisitable
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, false, pool_, NULL);
			EbwtRangeSourceDriver * dr2Rc_Bw = new EbwtRangeSourceDriver(
				*params, r2Rc_Bw, false, false, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_LEN, // allow 1 mismatch in rest of read
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, false, pool_, NULL);
			dr2RcVec->push_back(dr2Rc_Fw);
			dr2RcVec->push_back(dr2Rc_Bw);
		}

		RefAligner<String<Dna5> >* refAligner =
			new OneMMRefAligner<String<Dna5> >();

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

#endif /* ALIGNER_1MM_H_ */
