/*
 * aligner_seed_mm.h
 */

#ifndef ALIGNER_SEED_MM_H_
#define ALIGNER_SEED_MM_H_

#include <utility>
#include "aligner.h"
#include "hit.h"
#include "row_chaser.h"
#include "range_chaser.h"
#include "aligner_metrics.h"
#include "bowtie1_range_source.h"
#include "ds.h"

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class UnpairedSeedAlignerFactory : public AlignerFactory {

	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef EList<TRangeSrcDr*> TRangeSrcDrPtrVec;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;

public:
	UnpairedSeedAlignerFactory(
			Ebwt& ebwtFw,
			Ebwt* ebwtBw,
			uint32_t seedMms,
			uint32_t seedLen,
			int qualCutoff,
			int maxBts,
			HitSink& sink,
			const HitSinkPerThreadFactory& sinkPtFactory,
			RangeCache* cacheFw,
			RangeCache* cacheBw,
			uint32_t cacheLimit,
			ChunkPool *pool,
			BitPairReference* refs,
			EList<String<Dna5> >& os,
			uint32_t seed,
			AlignerMetrics *metrics) :
			ebwtFw_(ebwtFw),
			ebwtBw_(ebwtBw),
			seedMms_(seedMms),
			seedLen_(seedLen),
			qualCutoff_(qualCutoff),
			maxBts_(maxBts),
			sink_(sink),
			sinkPtFactory_(sinkPtFactory),
			cacheFw_(cacheFw),
			cacheBw_(cacheBw),
			cacheLimit_(cacheLimit),
			pool_(pool),
			refs_(refs),
			os_(os),
			metrics_(metrics)
	{
		assert(ebwtFw.isInMemory());
	}

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {
		HitSinkPerThread* sinkPt = sinkPtFactory_.create();
		EbwtSearchParams* params =
			new EbwtSearchParams(*sinkPt, os_);
		int *btCnt = new int[1];
		*btCnt = maxBts_;

		TRangeSrcDrPtrVec *drVec = new TRangeSrcDrPtrVec();
		if(seedMms_ == 0) {
			const int halfAndHalf = 0;
			bool mate1 = true;
			EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
				 ebwtBw_, true,  qualCutoff_, true,
				 halfAndHalf, false, metrics_);
			EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
				&ebwtFw_, false, qualCutoff_, true,
				halfAndHalf, false, metrics_);
			EbwtRangeSourceDriver * driverFw = new EbwtRangeSourceDriver(
				*params, rFw_Bw, true, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
				PIN_TO_SEED_EDGE, // "
				PIN_TO_SEED_EDGE, // "
				PIN_TO_SEED_EDGE, // "
				os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtRangeSourceDriver * driverRc = new EbwtRangeSourceDriver(
				*params, rRc_Fw, false, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
				PIN_TO_SEED_EDGE, // "
				PIN_TO_SEED_EDGE, // "
				PIN_TO_SEED_EDGE, // "
				os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			if(!gNofw) drVec->push_back(driverFw);
			if(!gNorc) drVec->push_back(driverRc);

		} else if(seedMms_ == 1) {
			const int halfAndHalf = 0;
			bool fw = true;
			bool mate1 = true;

			EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, true,
				 halfAndHalf, false, metrics_);
			EbwtRangeSourceFactory *rFw_BwSeed = new EbwtRangeSourceFactory(
				 ebwtBw_, fw,  qualCutoff_, true,
				 halfAndHalf, false, metrics_);
			EbwtRangeSource *rFw_FwSeedGen = new EbwtRangeSource(
				&ebwtFw_, fw,  qualCutoff_, false,
				 halfAndHalf, true,  metrics_);

			EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
				*params, rFw_Bw, fw, false,
				sink_, sinkPt, seedLen_,
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtRangeSourceDriverFactory * drFw_BwSeed = new EbwtRangeSourceDriverFactory(
				*params, rFw_BwSeed, fw, false,
				sink_, sinkPt, seedLen_,
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtRangeSourceDriver * drFw_FwSeedGen = new EbwtRangeSourceDriver(
				*params, rFw_FwSeedGen, fw, true,
				sink_, sinkPt, seedLen_,
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtSeededRangeSourceDriver * drFw_Seed = new EbwtSeededRangeSourceDriver(
				drFw_BwSeed, drFw_FwSeedGen, fw, seedLen_, mate1);

			fw = false;

			EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
				&ebwtFw_, fw, qualCutoff_, true,
				halfAndHalf, false, metrics_);
			EbwtRangeSourceFactory *rRc_FwSeed = new EbwtRangeSourceFactory(
				&ebwtFw_, fw, qualCutoff_, true,
				halfAndHalf, false, metrics_);
			EbwtRangeSource *rRc_BwSeedGen = new EbwtRangeSource(
				 ebwtBw_, fw, qualCutoff_, false,
				halfAndHalf, true,  metrics_);

			EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
				*params, rRc_Fw, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtRangeSourceDriverFactory * drRc_FwSeed = new EbwtRangeSourceDriverFactory(
				*params, rRc_FwSeed, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtRangeSourceDriver * drRc_BwSeedGen = new EbwtRangeSourceDriver(
				*params, rRc_BwSeedGen, fw, true, sink_, sinkPt,
				seedLen_,   // seedLen
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtSeededRangeSourceDriver * drRc_Seed = new EbwtSeededRangeSourceDriver(
				drRc_FwSeed, drRc_BwSeedGen, fw, seedLen_, true);

			if(!gNofw) {
				drVec->push_back(drFw_Bw);
				drVec->push_back(drFw_Seed);
			}
			if(!gNorc) {
				drVec->push_back(drRc_Fw);
				drVec->push_back(drRc_Seed);
			}
		} else if(seedMms_ == 2) {

			bool fw = true;
			bool mate1 = true;

			EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, true,  0, false, metrics_);
			EbwtRangeSourceFactory *rFw_BwSeed = new EbwtRangeSourceFactory(
				 ebwtBw_, fw,  qualCutoff_, true,  0, false, metrics_);
			EbwtRangeSource *rFw_FwSeedGen = new EbwtRangeSource(
				&ebwtFw_, fw,  qualCutoff_, false, 0, true, metrics_);
			EbwtRangeSource *rFw_BwHalf = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, false, 2,  false, metrics_);

			EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
				*params, rFw_Bw, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);
			EbwtRangeSourceDriverFactory * drFw_BwSeed = new EbwtRangeSourceDriverFactory(
				*params, rFw_BwSeed, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);
			EbwtRangeSourceDriver * drFw_FwSeedGen = new EbwtRangeSourceDriver(
				*params, rFw_FwSeedGen, fw, true, sink_, sinkPt,
				seedLen_,   // seedLen
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);
			EbwtSeededRangeSourceDriver * drFw_Seed = new EbwtSeededRangeSourceDriver(
				drFw_BwSeed, drFw_FwSeedGen, fw, seedLen_, mate1);
			EbwtRangeSourceDriver * drFw_BwHalf = new EbwtRangeSourceDriver(
				*params, rFw_BwHalf, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);

			fw = false;

			EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
				&ebwtFw_, fw, qualCutoff_, true,  0, false, metrics_);
			EbwtRangeSourceFactory *rRc_FwSeed = new EbwtRangeSourceFactory(
				&ebwtFw_, fw, qualCutoff_, true,  0, false, metrics_);
			EbwtRangeSource *rRc_BwSeedGen = new EbwtRangeSource(
				 ebwtBw_, fw, qualCutoff_, false, 0, true, metrics_);
			EbwtRangeSource *rRc_FwHalf = new EbwtRangeSource(
				&ebwtFw_, fw, qualCutoff_, false, 2,  false, metrics_);

			EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
				*params, rRc_Fw, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE, // no mismatches in hi half
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,    // up to 2 in lo half
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);
			EbwtRangeSourceDriverFactory * drRc_FwSeed = new EbwtRangeSourceDriverFactory(
				*params, rRc_FwSeed, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);
			EbwtRangeSourceDriver * drRc_BwSeedGen = new EbwtRangeSourceDriver(
				*params, rRc_BwSeedGen, fw, true, sink_, sinkPt,
				seedLen_,   // seedLen
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE, // no mismatches in lo half
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,    // up to 2 in hi half
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);
			EbwtSeededRangeSourceDriver * drRc_Seed = new EbwtSeededRangeSourceDriver(
				drRc_FwSeed, drRc_BwSeedGen, fw, seedLen_, true);
			EbwtRangeSourceDriver * drRc_FwHalf = new EbwtRangeSourceDriver(
				*params, rRc_FwHalf, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);

			if(!gNofw) {
				drVec->push_back(drFw_Bw);
				drVec->push_back(drFw_Seed);
				drVec->push_back(drFw_BwHalf);
			}
			if(!gNorc) {
				drVec->push_back(drRc_Fw);
				drVec->push_back(drRc_Seed);
				drVec->push_back(drRc_FwHalf);
			}
		} else if(seedMms_ > 2) {

			bool fw = true;
			bool mate1 = true;

			EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, true,
				 0, false, metrics_);

			// Partial and full aligners for alignments with 0
			// mismatches in the lo-half and up to 3 mismatches in the
			// hi-half
			EbwtRangeSourceFactory *rFw_BwSeed03 = new EbwtRangeSourceFactory(
				 ebwtBw_, fw,  qualCutoff_, true,
				 0, false, metrics_);
			EbwtRangeSource *rFw_FwSeedGen03 = new EbwtRangeSource(
				&ebwtFw_, fw,  qualCutoff_, false,
				 0, true,  metrics_);

			// Partial and full aligners for alignments with 1
			// mismatch in the lo-half and up to 2 mismatches in the
			// hi-half
			EbwtRangeSourceFactory *rFw_BwSeed12 = new EbwtRangeSourceFactory(
				 ebwtBw_, fw,  qualCutoff_, true,
				 0, false, metrics_);
			// Note: the following is half-and-half (unlike the 03 version)
			EbwtRangeSource *rFw_FwSeedGen12 = new EbwtRangeSource(
				&ebwtFw_, fw,  qualCutoff_, false,
				 3,  true,  metrics_);

			EbwtRangeSource *rFw_BwHalf12 = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, false,
				 2,  false, metrics_);

			EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
				*params, rFw_Bw, fw, false,
				sink_, sinkPt, seedLen_,
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE, // 0 mismatches in hi-half
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,    // up to 3 mismatches in lo-half
				os_, mate1, pool_, btCnt);

			EbwtRangeSourceDriverFactory * drFw_BwSeed03 = new EbwtRangeSourceDriverFactory(
				*params, rFw_BwSeed03, fw, false,
				sink_, sinkPt, seedLen_,
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);
			EbwtRangeSourceDriver * drFw_FwSeedGen03 = new EbwtRangeSourceDriver(
				*params, rFw_FwSeedGen03, fw, true,
				sink_, sinkPt, seedLen_,
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);
			EbwtSeededRangeSourceDriver * drFw_Seed03 = new EbwtSeededRangeSourceDriver(
				drFw_BwSeed03, drFw_FwSeedGen03, fw, seedLen_, mate1);

			EbwtRangeSourceDriverFactory * drFw_BwSeed12 = new EbwtRangeSourceDriverFactory(
				*params, rFw_BwSeed12, fw, false,
				sink_, sinkPt, seedLen_,
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);
			EbwtRangeSourceDriver * drFw_FwSeedGen12 = new EbwtRangeSourceDriver(
				*params, rFw_FwSeedGen12, fw, true,
				sink_, sinkPt, seedLen_,
				false,      // nudgeLeft
				PIN_TO_BEGINNING,
				PIN_TO_HI_HALF_EDGE, // 1-mismatch in lo-half
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,    // 1 or 2 mismatches in hi-half
				os_, mate1, pool_, btCnt);
			EbwtSeededRangeSourceDriver * drFw_Seed12 = new EbwtSeededRangeSourceDriver(
				drFw_BwSeed12, drFw_FwSeedGen12, fw, seedLen_, mate1);

			EbwtRangeSourceDriver * drFw_BwHalf12 = new EbwtRangeSourceDriver(
				*params, rFw_BwHalf12, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);

			fw = false;

			EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
				&ebwtFw_, fw,  qualCutoff_, true,
				0, false, metrics_);

			// Partial and full aligners for alignments with 0
			// mismatches in the lo-half and up to 3 mismatches in the
			// hi-half
			EbwtRangeSourceFactory *rRc_FwSeed03 = new EbwtRangeSourceFactory(
				&ebwtFw_, fw,  qualCutoff_, true,
				0, false, metrics_);
			EbwtRangeSource *rRc_BwSeedGen03 = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, false,
				0, true,   metrics_);

			// Partial and full aligners for alignments with 1
			// mismatch in the lo-half and up to 2 mismatches in the
			// hi-half
			EbwtRangeSourceFactory *rRc_FwSeed12 = new EbwtRangeSourceFactory(
				&ebwtFw_, fw,  qualCutoff_, true,
				0, false, metrics_);
			EbwtRangeSource *rRc_BwSeedGen12 = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, false,
				3,  true, metrics_);

			EbwtRangeSource *rRc_FwHalf12 = new EbwtRangeSource(
				&ebwtFw_, fw,  qualCutoff_, false,
				2,  false, metrics_);

			EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
				*params, rRc_Fw, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);

			EbwtRangeSourceDriverFactory * drRc_FwSeed03 = new EbwtRangeSourceDriverFactory(
				*params, rRc_FwSeed03, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);
			EbwtRangeSourceDriver * drRc_BwSeedGen03 = new EbwtRangeSourceDriver(
				*params, rRc_BwSeedGen03, fw, true, sink_, sinkPt,
				seedLen_,   // seedLen
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);
			EbwtSeededRangeSourceDriver * drRc_Seed03 = new EbwtSeededRangeSourceDriver(
				drRc_FwSeed03, drRc_BwSeedGen03, fw, seedLen_, mate1);

			EbwtRangeSourceDriverFactory * drRc_FwSeed12 = new EbwtRangeSourceDriverFactory(
				*params, rRc_FwSeed12, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);
			EbwtRangeSourceDriver * drRc_BwSeedGen12 = new EbwtRangeSourceDriver(
				*params, rRc_BwSeedGen12, fw, true, sink_, sinkPt,
				seedLen_,   // seedLen
				false,      // nudgeLeft
				PIN_TO_BEGINNING,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);
			EbwtSeededRangeSourceDriver * drRc_Seed12 = new EbwtSeededRangeSourceDriver(
				drRc_FwSeed12, drRc_BwSeedGen12, fw, seedLen_, mate1);

			EbwtRangeSourceDriver * drRc_FwHalf12 = new EbwtRangeSourceDriver(
				*params, rRc_FwHalf12, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				os_, mate1, pool_, btCnt);

			if(!gNofw) {
				drVec->push_back(drFw_Bw);
				drVec->push_back(drFw_Seed03);
				drVec->push_back(drFw_Seed12);
				drVec->push_back(drFw_BwHalf12);
			}
			if(!gNorc) {
				drVec->push_back(drRc_Fw);
				drVec->push_back(drRc_Seed03);
				drVec->push_back(drRc_Seed12);
				drVec->push_back(drRc_FwHalf12);
			}
		} else {
			cerr << "Unsupported --stateful mode: " << seedMms_ << endl;
		}
		TCostAwareRangeSrcDr* dr = new TCostAwareRangeSrcDr(drVec, false);
		delete drVec;

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(cacheLimit_, cacheFw_, cacheBw_, metrics_);

		ColorspaceDecoder *dec = new ColorspaceDecoder();

		return new UnpairedAlignerV2<EbwtRangeSource>(
			params, dr, rchase,
			sink_, sinkPtFactory_, sinkPt, os_, refs_,
			maxBts_, pool_, dec, btCnt, metrics_);
	}

private:
	Ebwt& ebwtFw_;
	Ebwt* ebwtBw_;
	const uint32_t seedMms_;
	const uint32_t seedLen_;
	const int qualCutoff_;
	const int maxBts_;
	HitSink& sink_;
	const HitSinkPerThreadFactory& sinkPtFactory_;
	RangeCache *cacheFw_;
	RangeCache *cacheBw_;
	const uint32_t cacheLimit_;
	ChunkPool *pool_;
	BitPairReference* refs_;
	EList<String<Dna5> >& os_;
	AlignerMetrics *metrics_;
};

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class PairedSeedAlignerFactory : public AlignerFactory {
	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef EList<TRangeSrcDr*> TRangeSrcDrPtrVec;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
public:
	PairedSeedAlignerFactory(
			Ebwt& ebwtFw,
			Ebwt* ebwtBw,
			uint32_t seedMms,
			uint32_t seedLen,
			int qualCutoff,
			int maxBts,
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
			ebwtBw_(ebwtBw),
			seedMms_(seedMms),
			seedLen_(seedLen),
			qualCutoff_(qualCutoff),
			maxBts_(maxBts),
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
			os_(os)
	{
		assert(ebwtFw.isInMemory());
		assert(ebwtBw->isInMemory());
	}

	/**
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
		RefAligner<String<Dna5> >* refAligner = NULL;
		int *btCnt = new int[1];
		*btCnt = maxBts_;
		if(seedMms_ == 0) {
			refAligner = new Seed0RefAligner<String<Dna5> >(seedLen_, qualCutoff_);
		} else if(seedMms_ == 1) {
			refAligner = new Seed1RefAligner<String<Dna5> >(seedLen_, qualCutoff_);
		} else if(seedMms_ == 2) {
			refAligner = new Seed2RefAligner<String<Dna5> >(seedLen_, qualCutoff_);
		} else {
			refAligner = new Seed3RefAligner<String<Dna5> >(seedLen_, qualCutoff_);
		}
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
		TRangeSrcDrPtrVec *dr1RcVec;
		TRangeSrcDrPtrVec *dr2FwVec;
		TRangeSrcDrPtrVec *dr2RcVec;
		dr1RcVec = dr1FwVec;
		dr2FwVec = dr1FwVec;
		dr2RcVec = dr1FwVec;
		if(seedMms_ == 0) {
			const int halfAndHalf = 0;
			if(do1Fw) {
				bool mate1 = true;
				bool fw = true;
				EbwtRangeSource *r1Fw_Bw = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, true,
					 halfAndHalf, false);
				EbwtRangeSourceDriver *dr1Fw_Bw = new EbwtRangeSourceDriver(
					*params, r1Fw_Bw, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft (not applicable)
					PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				dr1FwVec->push_back(dr1Fw_Bw);
			}
			if(do2Fw) {
				bool mate1 = false;
				bool fw = true;
				EbwtRangeSource *r2Fw_Bw = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, true, halfAndHalf,
					 false);
				EbwtRangeSourceDriver *dr2Fw_Bw = new EbwtRangeSourceDriver(
					*params, r2Fw_Bw, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft (not applicable)
					PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				dr2FwVec->push_back(dr2Fw_Bw);
			}
			if(do1Rc) {
				bool mate1 = true;
				bool fw = false;
				EbwtRangeSource *r1Rc_Fw = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, true,
					halfAndHalf, false);
				EbwtRangeSourceDriver *dr1Rc_Fw = new EbwtRangeSourceDriver(
					*params, r1Rc_Fw, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft (not applicable)
					PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				dr1RcVec->push_back(dr1Rc_Fw);
			}
			if(do2Rc) {
				bool mate1 = false;
				bool fw = false;
				EbwtRangeSource *r2Rc_Fw = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, true,
					halfAndHalf, false);
				EbwtRangeSourceDriver *dr2Rc_Fw = new EbwtRangeSourceDriver(
					*params, r2Rc_Fw, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft (not applicable)
					PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				dr2RcVec->push_back(dr2Rc_Fw);
			}
		} else if(seedMms_ == 1) {
			const int halfAndHalf = 0;
			if(do1Fw) {
				bool mate1 = true;
				bool fw = true;
				EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, true,
					halfAndHalf, false);
				EbwtRangeSourceFactory *rFw_BwSeed = new EbwtRangeSourceFactory(
					 ebwtBw_, fw,  qualCutoff_, true,
					halfAndHalf, false);
				EbwtRangeSource *rFw_FwSeedGen = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, false,
					halfAndHalf, true);
				EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
					*params, rFw_Bw, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriverFactory * drFw_BwSeed = new EbwtRangeSourceDriverFactory(
					*params, rFw_BwSeed, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriver * drFw_FwSeedGen = new EbwtRangeSourceDriver(
					*params, rFw_FwSeedGen, fw, true, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtSeededRangeSourceDriver * drFw_Seed = new EbwtSeededRangeSourceDriver(
					drFw_BwSeed, drFw_FwSeedGen, fw, seedLen_, mate1);
				dr1FwVec->push_back(drFw_Bw);
				dr1FwVec->push_back(drFw_Seed);
			}
			if(do2Fw) {
				bool mate1 = false;
				bool fw = true;
				EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, true,
					halfAndHalf, false);
				EbwtRangeSourceFactory *rFw_BwSeed = new EbwtRangeSourceFactory(
					 ebwtBw_, fw,  qualCutoff_, true,
					halfAndHalf, false);
				EbwtRangeSource *rFw_FwSeedGen = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, false,
					halfAndHalf, true);
				EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
					*params, rFw_Bw, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriverFactory * drFw_BwSeed = new EbwtRangeSourceDriverFactory(
					*params, rFw_BwSeed, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriver * drFw_FwSeedGen = new EbwtRangeSourceDriver(
					*params, rFw_FwSeedGen, fw, true, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtSeededRangeSourceDriver * drFw_Seed = new EbwtSeededRangeSourceDriver(
					drFw_BwSeed, drFw_FwSeedGen, fw, seedLen_, mate1);
				dr2FwVec->push_back(drFw_Bw);
				dr2FwVec->push_back(drFw_Seed);
			}
			if(do1Rc) {
				bool mate1 = true;
				bool fw = false;
				EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, true,
					halfAndHalf, false);
				EbwtRangeSourceFactory *rRc_FwSeed = new EbwtRangeSourceFactory(
					&ebwtFw_, fw,  qualCutoff_, true,
					halfAndHalf, false);
				EbwtRangeSource *rRc_BwSeedGen = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, false,
					halfAndHalf, true);
				EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
					*params, rRc_Fw, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriverFactory * drRc_FwSeed = new EbwtRangeSourceDriverFactory(
					*params, rRc_FwSeed, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriver * drRc_BwSeedGen = new EbwtRangeSourceDriver(
					*params, rRc_BwSeedGen, fw, true, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtSeededRangeSourceDriver * drRc_Seed = new EbwtSeededRangeSourceDriver(
					drRc_FwSeed, drRc_BwSeedGen, fw, seedLen_, mate1);
				dr1RcVec->push_back(drRc_Fw);
				dr1RcVec->push_back(drRc_Seed);
			}
			if(do2Rc) {
				bool mate1 = false;
				bool fw = false;
				EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, true,
					halfAndHalf, false);
				EbwtRangeSourceFactory *rRc_FwSeed = new EbwtRangeSourceFactory(
					&ebwtFw_, fw,  qualCutoff_, true,
					halfAndHalf, false);
				EbwtRangeSource *rRc_BwSeedGen = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, false,
					halfAndHalf, true);
				EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
					*params, rRc_Fw, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriverFactory * drRc_FwSeed = new EbwtRangeSourceDriverFactory(
					*params, rRc_FwSeed, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriver * drRc_BwSeedGen = new EbwtRangeSourceDriver(
					*params, rRc_BwSeedGen, fw, true, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtSeededRangeSourceDriver * drRc_Seed = new EbwtSeededRangeSourceDriver(
					drRc_FwSeed, drRc_BwSeedGen, fw, seedLen_, mate1);
				dr2RcVec->push_back(drRc_Fw);
				dr2RcVec->push_back(drRc_Seed);
			}
		} else if(seedMms_ > 1) {
			bool two = seedMms_ == 2;
			if(do1Fw) {
				bool mate1 = true;
				bool fw = true;
				EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, true,  0, false);
				EbwtRangeSourceFactory *rFw_BwSeed = new EbwtRangeSourceFactory(
					 ebwtBw_, fw,  qualCutoff_, true,  0, false);
				EbwtRangeSource *rFw_FwSeedGen = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, false, 0, true);
				EbwtRangeSource *rFw_BwHalf = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, false, 2,  false);
				EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
					*params, rFw_Bw, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				EbwtRangeSourceDriverFactory * drFw_BwSeed = new EbwtRangeSourceDriverFactory(
					*params, rFw_BwSeed, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				EbwtRangeSourceDriver * drFw_FwSeedGen = new EbwtRangeSourceDriver(
					*params, rFw_FwSeedGen, fw, true, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				EbwtSeededRangeSourceDriver * drFw_Seed = new EbwtSeededRangeSourceDriver(
					drFw_BwSeed, drFw_FwSeedGen, fw, seedLen_, mate1);
				EbwtRangeSourceDriverFactory * drFw_BwSeed12 = NULL;
				if(!two) {
					EbwtRangeSourceFactory *rFw_BwSeed12 = new EbwtRangeSourceFactory(
						 ebwtBw_, fw,  qualCutoff_, true,  0, false);
					drFw_BwSeed12 = new EbwtRangeSourceDriverFactory(
						*params, rFw_BwSeed12, fw, false, sink_, sinkPt,
						seedLen_,   // seedLen
						true,       // nudgeLeft
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						os_, mate1, pool_, btCnt);
				}
				EbwtRangeSourceDriver * drFw_FwSeedGen12 = NULL;
				if(!two) {
					EbwtRangeSource *rFw_FwSeedGen12 = new EbwtRangeSource(
						&ebwtFw_, fw,  qualCutoff_, false, 3,  true);
					drFw_FwSeedGen12 = new EbwtRangeSourceDriver(
						*params, rFw_FwSeedGen12, fw, true, sink_, sinkPt,
						seedLen_,   // seedLen
						false,      // nudgeLeft
						PIN_TO_BEGINNING,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_SEED_EDGE,
						os_, mate1, pool_, btCnt);
				}
				EbwtSeededRangeSourceDriver * drFw_Seed12 = NULL;
				if(!two) {
					drFw_Seed12 = new EbwtSeededRangeSourceDriver(
						drFw_BwSeed12, drFw_FwSeedGen12, fw, seedLen_, mate1);
				}
				EbwtRangeSourceDriver * drFw_BwHalf = new EbwtRangeSourceDriver(
					*params, rFw_BwHalf, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				dr1FwVec->push_back(drFw_Bw);
				dr1FwVec->push_back(drFw_Seed);
				if(drFw_Seed12 != NULL) {
					dr1FwVec->push_back(drFw_Seed12);
				}
				dr1FwVec->push_back(drFw_BwHalf);
			}
			if(do2Fw) {
				bool mate1 = false;
				bool fw = true;
				EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, true,  0, false);
				EbwtRangeSourceFactory *rFw_BwSeed = new EbwtRangeSourceFactory(
					 ebwtBw_, fw,  qualCutoff_, true,  0, false);
				EbwtRangeSource *rFw_FwSeedGen = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, false, 0, true);
				EbwtRangeSource *rFw_BwHalf = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, false, 2,  false);

				EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
					*params, rFw_Bw, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				EbwtRangeSourceDriverFactory * drFw_BwSeed = new EbwtRangeSourceDriverFactory(
					*params, rFw_BwSeed, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				EbwtRangeSourceDriver * drFw_FwSeedGen = new EbwtRangeSourceDriver(
					*params, rFw_FwSeedGen, fw, true, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				EbwtSeededRangeSourceDriver * drFw_Seed = new EbwtSeededRangeSourceDriver(
					drFw_BwSeed, drFw_FwSeedGen, fw, seedLen_, mate1);
				EbwtRangeSourceDriverFactory * drFw_BwSeed12 = NULL;
				if(!two) {
					EbwtRangeSourceFactory *rFw_BwSeed12 = new EbwtRangeSourceFactory(
						 ebwtBw_, fw,  qualCutoff_, true,  0, false);
					drFw_BwSeed12 = new EbwtRangeSourceDriverFactory(
						*params, rFw_BwSeed12, fw, false, sink_, sinkPt,
						seedLen_,   // seedLen
						true,       // nudgeLeft
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						os_, mate1, pool_, btCnt);
				}
				EbwtRangeSourceDriver * drFw_FwSeedGen12 = NULL;
				if(!two) {
					EbwtRangeSource *rFw_FwSeedGen12 = new EbwtRangeSource(
						&ebwtFw_, fw,  qualCutoff_, false, 3,  true);
					drFw_FwSeedGen12 = new EbwtRangeSourceDriver(
						*params, rFw_FwSeedGen12, fw, true, sink_, sinkPt,
						seedLen_,   // seedLen
						false,      // nudgeLeft
						PIN_TO_BEGINNING,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_SEED_EDGE,
						os_, mate1, pool_, btCnt);
				}
				EbwtSeededRangeSourceDriver * drFw_Seed12 = NULL;
				if(!two) {
					drFw_Seed12 = new EbwtSeededRangeSourceDriver(
						drFw_BwSeed12, drFw_FwSeedGen12, fw, seedLen_, mate1);
				}
				EbwtRangeSourceDriver * drFw_BwHalf = new EbwtRangeSourceDriver(
					*params, rFw_BwHalf, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				dr2FwVec->push_back(drFw_Bw);
				dr2FwVec->push_back(drFw_Seed);
				if(drFw_Seed12 != NULL) {
					dr2FwVec->push_back(drFw_Seed12);
				}
				dr2FwVec->push_back(drFw_BwHalf);
			}
			if(do1Rc) {
				bool mate1 = true;
				bool fw = false;
				EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, true,  0, false);
				EbwtRangeSourceFactory *rRc_FwSeed = new EbwtRangeSourceFactory(
					&ebwtFw_, fw,  qualCutoff_, true,  0, false);
				EbwtRangeSource *rRc_BwSeedGen = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, false, 0, true);
				EbwtRangeSource *rRc_FwHalf = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, false, 2,  false);

				EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
					*params, rRc_Fw, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				EbwtRangeSourceDriverFactory * drRc_FwSeed = new EbwtRangeSourceDriverFactory(
					*params, rRc_FwSeed, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				EbwtRangeSourceDriver * drRc_BwSeedGen = new EbwtRangeSourceDriver(
					*params, rRc_BwSeedGen, fw, true, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				EbwtSeededRangeSourceDriver * drRc_Seed = new EbwtSeededRangeSourceDriver(
					drRc_FwSeed, drRc_BwSeedGen, fw, seedLen_, mate1);
				EbwtRangeSourceDriverFactory * drRc_FwSeed12 = NULL;
				if(!two) {
					EbwtRangeSourceFactory *rRc_FwSeed12 = new EbwtRangeSourceFactory(
						&ebwtFw_, fw,  qualCutoff_, true,  0, false);
					drRc_FwSeed12 = new EbwtRangeSourceDriverFactory(
						*params, rRc_FwSeed12, fw, false, sink_, sinkPt,
						seedLen_,   // seedLen
						true,       // nudgeLeft
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						os_, mate1, pool_, btCnt);
				}
				EbwtRangeSourceDriver * drRc_BwSeedGen12 = NULL;
				if(!two) {
					EbwtRangeSource *rRc_BwSeedGen12 = new EbwtRangeSource(
						 ebwtBw_, fw,  qualCutoff_, false, 3,  true);
					drRc_BwSeedGen12 = new EbwtRangeSourceDriver(
						*params, rRc_BwSeedGen12, fw, true, sink_, sinkPt,
						seedLen_,   // seedLen
						false,      // nudgeLeft
						PIN_TO_BEGINNING,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_SEED_EDGE,
						os_, mate1, pool_, btCnt);
				}
				EbwtSeededRangeSourceDriver * drRc_Seed12 = NULL;
				if(!two) {
					drRc_Seed12 = new EbwtSeededRangeSourceDriver(
						drRc_FwSeed12, drRc_BwSeedGen12, fw, seedLen_, mate1);
				}
				EbwtRangeSourceDriver * drRc_FwHalf = new EbwtRangeSourceDriver(
					*params, rRc_FwHalf, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				dr1RcVec->push_back(drRc_Fw);
				dr1RcVec->push_back(drRc_Seed);
				if(drRc_Seed12 != NULL) {
					dr1RcVec->push_back(drRc_Seed12);
				}
				dr1RcVec->push_back(drRc_FwHalf);
			}
			if(do2Rc) {
				bool mate1 = false;
				bool fw = false;
				EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, true,  0, false);
				EbwtRangeSourceFactory *rRc_FwSeed = new EbwtRangeSourceFactory(
					&ebwtFw_, fw,  qualCutoff_, true,  0, false);
				EbwtRangeSource *rRc_BwSeedGen = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, false, 0, true);
				EbwtRangeSource *rRc_FwHalf = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, false, 2,  false);

				EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
					*params, rRc_Fw, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				EbwtRangeSourceDriverFactory * drRc_FwSeed = new EbwtRangeSourceDriverFactory(
					*params, rRc_FwSeed, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				EbwtRangeSourceDriver * drRc_BwSeedGen = new EbwtRangeSourceDriver(
					*params, rRc_BwSeedGen, fw, true, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				EbwtSeededRangeSourceDriver * drRc_Seed = new EbwtSeededRangeSourceDriver(
					drRc_FwSeed, drRc_BwSeedGen, fw, seedLen_, mate1);
				EbwtRangeSourceDriverFactory * drRc_FwSeed12 = NULL;
				if(!two) {
					EbwtRangeSourceFactory *rRc_FwSeed12 = new EbwtRangeSourceFactory(
						&ebwtFw_, fw,  qualCutoff_, true,  0, false);
					drRc_FwSeed12 = new EbwtRangeSourceDriverFactory(
						*params, rRc_FwSeed12, fw, false, sink_, sinkPt,
						seedLen_,   // seedLen
						true,       // nudgeLeft
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						os_, mate1, pool_, btCnt);
				}
				EbwtRangeSourceDriver * drRc_BwSeedGen12 = NULL;
				if(!two) {
					EbwtRangeSource *rRc_BwSeedGen12 = new EbwtRangeSource(
						 ebwtBw_, fw,  qualCutoff_, false, 3,  true);
					drRc_BwSeedGen12 = new EbwtRangeSourceDriver(
						*params, rRc_BwSeedGen12, fw, true, sink_, sinkPt,
						seedLen_,   // seedLen
						false,      // nudgeLeft
						PIN_TO_BEGINNING,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_SEED_EDGE,
						os_, mate1, pool_, btCnt);
				}
				EbwtSeededRangeSourceDriver * drRc_Seed12 = NULL;
				if(!two) {
					drRc_Seed12 = new EbwtSeededRangeSourceDriver(
						drRc_FwSeed12, drRc_BwSeedGen12, fw, seedLen_, mate1);
				}
				EbwtRangeSourceDriver * drRc_FwHalf = new EbwtRangeSourceDriver(
					*params, rRc_FwHalf, fw, false, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, mate1, pool_, btCnt);
				dr2RcVec->push_back(drRc_Fw);
				dr2RcVec->push_back(drRc_Seed);
				if(drRc_Seed12 != NULL) {
					dr2RcVec->push_back(drRc_Seed12);
				}
				dr2RcVec->push_back(drRc_FwHalf);
			}
		} else {
			cerr << "Unsupported --stateful mode: " << seedMms_ << endl;
		}
		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(cacheLimit_, cacheFw_, cacheBw_);

		ColorspaceDecoder *dec = new ColorspaceDecoder();

		// We dumped all the drivers into dr1FwVec
		PairedBWAlignerV2<EbwtRangeSource>* al = new PairedBWAlignerV2<EbwtRangeSource>(
			params, paramsSe1, paramsSe2,
			new TCostAwareRangeSrcDr(dr1FwVec, true),
			refAligner, rchase, sink_, sinkPtFactory_, sinkPt,
			sinkPtSe1, sinkPtSe2,
			mixedAttemptLim_, refs_, maxBts_, pool_, dec, btCnt);
		delete dr1FwVec;
		return al;
	}

private:
	Ebwt& ebwtFw_;
	Ebwt* ebwtBw_;
	const uint32_t seedMms_;
	const uint32_t seedLen_;
	const int qualCutoff_;
	const int maxBts_;
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
};

#endif /* ALIGNER_SEED_MM_H_ */
