/*
 * aligner_0mm.h
 */

#ifndef ALIGNER_0MM_H_
#define ALIGNER_0MM_H_

#include <utility>
#include "aligner.h"
#include "hit.h"
#include "row_chaser.h"
#include "range_chaser.h"
#include "bowtie1_range_source.h"
#include "ds.h"

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class UnpairedExactAlignerV1Factory : public AlignerFactory {

	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef EList<TRangeSrcDr*> TRangeSrcDrPtrVec;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;

public:
	UnpairedExactAlignerV1Factory(
			Ebwt& ebwtFw,
			Ebwt* ebwtBw,
			HitSink& sink,
			const HitSinkPerThreadFactory& sinkPtFactory,
			RangeCache* cacheFw,
			RangeCache* cacheBw,
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
	}

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {
		HitSinkPerThread* sinkPt = sinkPtFactory_.create();
		EbwtSearchParams* params =
			new EbwtSearchParams(*sinkPt, os_, true, true);

		const int halfAndHalf = 0;
		const bool seeded = false;

		EbwtRangeSource *rFw = new EbwtRangeSource(
			&ebwtFw_, true,  0xffffffff, true, halfAndHalf, seeded);
		EbwtRangeSource *rRc = new EbwtRangeSource(
			&ebwtFw_, false, 0xffffffff, true, halfAndHalf, seeded);

		EbwtRangeSourceDriver * driverFw = new EbwtRangeSourceDriver(
			*params, rFw, true, false, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_LEN, // whole alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, true, pool_, NULL);
		EbwtRangeSourceDriver * driverRc = new EbwtRangeSourceDriver(
			*params, rRc, false, false, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_LEN, // whole alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, true, pool_, NULL);
		TRangeSrcDrPtrVec *drVec = new TRangeSrcDrPtrVec();
		if(!gNofw) drVec->push_back(driverFw);
		if(!gNorc) drVec->push_back(driverRc);
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
class PairedExactAlignerV1Factory : public AlignerFactory {
	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
	typedef EList<TRangeSrcDr*> TRangeSrcDrPtrVec;
public:
	PairedExactAlignerV1Factory(
			Ebwt& ebwtFw,
			Ebwt* ebwtBw,
			HitSink& sink,
			const HitSinkPerThreadFactory& sinkPtFactory,
			bool dontReconcile,
			uint32_t symCeil,
			uint32_t mixedThresh,
			uint32_t mixedAttemptLim,
			RangeCache* cacheFw,
			RangeCache* cacheBw,
			uint32_t cacheLimit,
			ChunkPool *pool,
			BitPairReference* refs,
			EList<String<Dna5> >& os,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
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
		assert(ebwtFw.isInMemory());
	}

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {
		HitSinkPerThread* sinkPt = sinkPtFactory_.createMult(2);
		HitSinkPerThread* sinkPtSe1 = NULL, * sinkPtSe2 = NULL;
		EbwtSearchParams* params =
			new EbwtSearchParams(*sinkPt, os_, true, true);
		EbwtSearchParams* paramsSe1 = NULL, * paramsSe2 = NULL;
		if(gReportSe) {
			sinkPtSe1 = sinkPtFactory_.create();
			sinkPtSe2 = sinkPtFactory_.create();
			paramsSe1 =
				new EbwtSearchParams(*sinkPtSe1, os_, true, true);
			paramsSe2 =
				new EbwtSearchParams(*sinkPtSe2, os_, true, true);
		}

		const int halfAndHalf = 0;
		const bool seeded = false;

		bool do1Fw = true;
		bool do1Rc = true;
		bool do2Fw = true;
		bool do2Rc = true;
		if(gNofw) {
			if(gMate1fw) do1Fw = false;
			else         do1Rc = false;
			if(gMate2fw) do2Fw = false;
			else         do2Rc = false;
		}
		if(gNorc) {
			if(gMate1fw) do1Rc = false;
			else         do1Fw = false;
			if(gMate2fw) do2Rc = false;
			else         do2Fw = false;
		}

		EbwtRangeSource *r1Fw = NULL;
		EbwtRangeSource *r1Rc = NULL;
		TRangeSrcDr * driver1Fw = NULL;
		TRangeSrcDr * driver1Rc = NULL;
		EbwtRangeSource *r2Fw = NULL;
		EbwtRangeSource *r2Rc = NULL;
		TRangeSrcDr * driver2Fw = NULL;
		TRangeSrcDr * driver2Rc = NULL;
		if(do1Fw) {
			r1Fw = new EbwtRangeSource(
				&ebwtFw_, true,  0xffffffff, true, halfAndHalf, seeded);
			driver1Fw = new EbwtRangeSourceDriver(
				*params, r1Fw, true, false, sink_, sinkPt,
				0,          // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_LEN, // whole alignment is unrevisitable
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, true, pool_, NULL);
		}
		if(do2Fw) {
			r2Fw = new EbwtRangeSource(
				&ebwtFw_, true,  0xffffffff, true, halfAndHalf, seeded);
			driver2Fw = new EbwtRangeSourceDriver(
				*params, r2Fw, true, false, sink_, sinkPt,
				0,          // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_LEN, // whole alignment is unrevisitable
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, false, pool_, NULL);
		}
		if(do1Rc) {
			r1Rc = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, true, halfAndHalf, seeded);
			driver1Rc = new EbwtRangeSourceDriver(
				*params, r1Rc, false, false, sink_, sinkPt,
				0,          // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_LEN, // whole alignment is unrevisitable
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, true, pool_, NULL);
		}
		if(do2Rc) {
			r2Rc = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, true, halfAndHalf, seeded);
			driver2Rc = new EbwtRangeSourceDriver(
				*params, r2Rc, false, false, sink_, sinkPt,
				0,          // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_LEN, // whole alignment is unrevisitable
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, false, pool_, NULL);
		}

		RefAligner<String<Dna5> >* refAligner
			= new ExactRefAligner<String<Dna5> >();

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(cacheLimit_, cacheFw_, cacheBw_);

		ColorspaceDecoder *dec = new ColorspaceDecoder();

		TRangeSrcDrPtrVec *drVec = new TRangeSrcDrPtrVec();
		if(driver1Fw != NULL) drVec->push_back(driver1Fw);
		if(driver1Rc != NULL) drVec->push_back(driver1Rc);
		if(driver2Fw != NULL) drVec->push_back(driver2Fw);
		if(driver2Rc != NULL) drVec->push_back(driver2Rc);
		PairedBWAlignerV2<EbwtRangeSource>* al = new PairedBWAlignerV2<EbwtRangeSource>(
			params, paramsSe1, paramsSe2,
			new TCostAwareRangeSrcDr(drVec, true),
			refAligner,
			rchase, sink_, sinkPtFactory_, sinkPt,
			sinkPtSe1, sinkPtSe2,
			mixedAttemptLim_, refs_,
			INT_MAX, pool_, dec, NULL);
		delete drVec;
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

#endif /* ALIGNER_0MM_H_ */
