/**
 * aligner.h
 *
 * A generic class providing a stateful way to find alignments.
 */

#ifndef ALIGNER_H_
#define ALIGNER_H_

#include <iostream>
#include <set>
#include <stdint.h>
#include "assert_helpers.h"
#include "ebwt.h"
#include "pat.h"
#include "range.h"
#include "range_source.h"
#include "range_chaser.h"
#include "ref_aligner.h"
#include "reference.h"
#include "aligner_metrics.h"
#include "search_globals.h"
#include "sstring.h"
#include "ds.h"
#include "read.h"

/**
 * State machine for carrying out an alignment, which usually consists
 * of a series of phases that conduct different alignments using
 * different backtracking constraints.
 *
 * Each Aligner should have a dedicated PatternSourcePerThread.
 */
class Aligner {
public:
	Aligner(bool _done) :
		done(_done), patsrc_(NULL), bufa_(NULL), bufb_(NULL) { }

	virtual ~Aligner() { }
	/// Advance the range search by one memory op
	virtual bool advance() = 0;

	/// Prepare Aligner for the next read
	virtual void setQuery(PatternSourcePerThread *patsrc) {
		assert(patsrc != NULL);
		patsrc_ = patsrc;
		bufa_ = &patsrc->bufa();
		assert(bufa_ != NULL);
		bufb_ = &patsrc->bufb();
		alen_ = (uint32_t)bufa_->length();
		blen_ = (uint32_t)((bufb_ != NULL) ? bufb_->length() : 0);
		rand_.init(bufa_->seed);
	}

	/**
	 * Set to true if all searching w/r/t the current query is
	 * finished or if there is no current query.
	 */
	bool done;

protected:

	// Current read pair
	PatternSourcePerThread* patsrc_;
	Read* bufa_;
	uint32_t alen_;
	Read* bufb_;
	uint32_t blen_;
	RandomSource rand_;
};

/**
 * Abstract parent factory class for constructing aligners of all kinds.
 */
class AlignerFactory {
public:
	virtual ~AlignerFactory() { }
	virtual Aligner* create() const = 0;

	/**
	 * Allocate a vector of n Aligners; use destroy(std::vector...) to
	 * free the memory.
	 */
	virtual EList<Aligner*>* create(uint32_t n) const {
		EList<Aligner*>* v = new EList<Aligner*>;
		for(uint32_t i = 0; i < n; i++) {
			v->push_back(create());
			assert(v->back() != NULL);
		}
		return v;
	}

	/// Free memory associated with the aligner
	virtual void destroy(Aligner* al) const {
		assert(al != NULL);
		// Free the Aligner
		delete al;
	}

	/// Free memory associated with an aligner list
	virtual void destroy(EList<Aligner*>* als) const {
		assert(als != NULL);
		// Free all of the Aligners
		for(size_t i = 0; i < als->size(); i++) {
			if((*als)[i] != NULL) {
				delete (*als)[i];
				(*als)[i] = NULL;
			}
		}
		// Free the vector
		delete als;
	}
};

/**
 * Coordinates multiple aligners of the same type (i.e. either all
 * single-end or all paired-end).
 */
class MultiAligner {
public:
	MultiAligner(
			uint32_t n,
			uint32_t qUpto,
			const AlignerFactory& alignFact,
			const PatternSourcePerThreadFactory& patsrcFact) :
			n_(n), qUpto_(qUpto),
			alignFact_(alignFact), patsrcFact_(patsrcFact),
			aligners_(NULL), patsrcs_(NULL)
	{
		aligners_ = alignFact_.create(n_);
		assert(aligners_ != NULL);
		patsrcs_ = patsrcFact_.create(n_);
		assert(patsrcs_ != NULL);
	}

	/// Free memory associated with the aligners and their pattern sources.
	virtual ~MultiAligner() {
		alignFact_.destroy(aligners_);
		patsrcFact_.destroy(patsrcs_);
	}

	/**
	 * Advance an array of aligners in parallel, using prefetches to
	 * try to hide all the latency.
	 */
	void run() {
		bool done = false;
		while(!done) {
			done = true;
			for(uint32_t i = 0; i < n_; i++) {
				if(!(*aligners_)[i]->done) {
					// Advance an aligner already in progress
					done = false;
					(*aligners_)[i]->advance();
				} else {
					// Get a new read and initialize an aligner with it
					(*patsrcs_)[i]->nextReadPair(true);
					if(!(*patsrcs_)[i]->empty() && (*patsrcs_)[i]->patid() < qUpto_) {
						(*aligners_)[i]->setQuery((*patsrcs_)[i]);
						assert(!(*aligners_)[i]->done);
						done = false;
					} else {
						// No more reads; if done == true, it remains
						// true
					}
				}
			}
		}
	}

protected:
	uint32_t n_;     /// Number of aligners
	uint32_t qUpto_; /// Number of reads to align before stopping
	const AlignerFactory&                  alignFact_;
	const PatternSourcePerThreadFactory&   patsrcFact_;
	EList<Aligner *>*                aligners_;
	EList<PatternSourcePerThread *>* patsrcs_;
};

/**
 * Coordinates multiple single-end and paired-end aligners, routing
 * reads to one or the other type as appropriate.
 */
class MixedMultiAligner {
public:
	MixedMultiAligner(
			uint32_t n,
			uint32_t qUpto,
			const AlignerFactory& alignSEFact,
			const AlignerFactory& alignPEFact,
			const PatternSourcePerThreadFactory& patsrcFact) :
			n_(n), qUpto_(qUpto),
			alignSEFact_(alignSEFact),
			alignPEFact_(alignPEFact),
			patsrcFact_(patsrcFact),
			alignersSE_(NULL),
			alignersPE_(NULL),
			seOrPe_(NULL),
			patsrcs_(NULL)
	{
		// Instantiate all single-end aligners
		alignersSE_ = alignSEFact_.create(n_);
		assert(alignersSE_ != NULL);
		// Instantiate all paired-end aligners
		alignersPE_ = alignPEFact_.create(n_);
		assert(alignersPE_ != NULL);
		// Allocate array of boolean flags indicating whether each of
		// the slots is currently using the single-end or paired-end
		// aligner
		seOrPe_ = new bool[n_];
		for(uint32_t i = 0; i < n_; i++) {
			seOrPe_[i] = true;
		}
		// Instantiate all read sources
		patsrcs_ = patsrcFact_.create(n_);
		assert(patsrcs_ != NULL);
	}

	/// Free memory associated with the aligners and their pattern sources.
	virtual ~MixedMultiAligner() {
		alignSEFact_.destroy(alignersSE_);
		alignPEFact_.destroy(alignersPE_);
		patsrcFact_.destroy(patsrcs_);
		delete[] seOrPe_;
	}

	/**
	 * Advance an array of aligners in parallel, using prefetches to
	 * try to hide all the latency.
	 */
	void run(bool verbose = false) {
		bool done = false;
		bool first = true;
		if(n_ == 1) {
			Aligner *al = seOrPe_[0] ? (*alignersSE_)[0] : (*alignersPE_)[0];
			PatternSourcePerThread *ps = (*patsrcs_)[0];
			while(!done) {
				done = true;
				if(!first && !al->done) {
					// Advance an aligner already in progress; this is
					// the common case
					done = false;
					al->advance();
				} else {
					// Get a new read
					ps->nextReadPair(true);
					if(ps->patid() < qUpto_ && !ps->empty()) {
						assert_eq(ps->bufa().color, gColor);
						if(ps->paired()) {
							// Read currently in buffer is paired-end
							assert_eq(ps->bufb().color, gColor);
							(*alignersPE_)[0]->setQuery(ps);
							al = (*alignersPE_)[0];
							seOrPe_[0] = false; // false -> paired
						} else {
							// Read currently in buffer is single-end
							(*alignersSE_)[0]->setQuery(ps);
							al = (*alignersSE_)[0];
							seOrPe_[0] = true; // true = unpaired
						}
						done = false;
					} else {
						// No more reads; if done == true, it remains
						// true
					}
				}
				first = false;
			}
		} else {
			while(!done) {
				done = true;
				for(uint32_t i = 0; i < n_; i++) {
					Aligner *al = seOrPe_[i] ? (*alignersSE_)[i] :
											   (*alignersPE_)[i];
					if(!first && !al->done) {
						// Advance an aligner already in progress; this is
						// the common case
						done = false;
						al->advance();
					} else {
						// Feed a new read to a vacant aligner
						PatternSourcePerThread *ps = (*patsrcs_)[i];
						// Get a new read
						ps->nextReadPair(true);
						if(ps->patid() < qUpto_ && !ps->empty()) {
							if(ps->paired()) {
								// Read currently in buffer is paired-end
								(*alignersPE_)[i]->setQuery(ps);
								seOrPe_[i] = false; // false -> paired
							} else {
								// Read currently in buffer is single-end
								(*alignersSE_)[i]->setQuery(ps);
								seOrPe_[i] = true; // true = unpaired
							}
							done = false;
						} else {
							// No more reads; if done == true, it remains
							// true
						}
					}
				}
				first = false;
			}
		}
	}

protected:
	uint32_t n_;     /// Number of aligners
	uint32_t qUpto_; /// Number of reads to align before stopping
	const AlignerFactory&                  alignSEFact_;
	const AlignerFactory&                  alignPEFact_;
	const PatternSourcePerThreadFactory&   patsrcFact_;
	EList<Aligner *>*                alignersSE_;
	EList<Aligner *>*                alignersPE_;
	bool *                                 seOrPe_;
	EList<PatternSourcePerThread *>* patsrcs_;
};

/**
 * An aligner for finding exact matches of unpaired reads.  Always
 * tries the forward-oriented version of the read before the reverse-
 * oriented read.
 */
template<typename TRangeSource>
class UnpairedAlignerV2 : public Aligner {
	typedef RangeSourceDriver<TRangeSource> TDriver;
public:
	UnpairedAlignerV2(
		EbwtSearchParams* params,
		TDriver* driver,
		RangeChaser* rchase,
		HitSink& sink,
		const HitSinkPerThreadFactory& sinkPtFactory,
		HitSinkPerThread* sinkPt,
		EList<SString<char> >& os,
		const BitPairReference* refs,
		int maxBts,
		ChunkPool *pool,
		ColorspaceDecoder * dec,
		int *btCnt,
		AlignerMetrics *metrics) :
		Aligner(true),
		refs_(refs),
		doneFirst_(true),
		firstIsFw_(true),
		chase_(false),
		sinkPtFactory_(sinkPtFactory),
		sinkPt_(sinkPt),
		params_(params),
		rchase_(rchase),
		driver_(driver),
		maxBts_(maxBts),
		pool_(pool),
		btCnt_(btCnt),
		metrics_(metrics),
		dec_(dec)
	{
		assert(pool_   != NULL);
		assert(sinkPt_ != NULL);
		assert(params_ != NULL);
		assert(driver_ != NULL);
		assert(dec_    != NULL);
	}

	virtual ~UnpairedAlignerV2() {
		delete driver_;  driver_  = NULL;
		delete params_;  params_  = NULL;
		delete rchase_;  rchase_  = NULL;
		delete[] btCnt_; btCnt_   = NULL;
		delete dec_;     dec_     = NULL;
		sinkPtFactory_.destroy(sinkPt_); sinkPt_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc) {
		Aligner::setQuery(patsrc); // set fields & random seed
		if(metrics_ != NULL) {
			metrics_->nextRead(patsrc->bufa().patFw);
		}
		pool_->reset(&patsrc->bufa().name, patsrc->patid());
		if(patsrc->bufa().length() < 4) {
			if(!gQuiet) {
				cerr << "Warning: Skipping read " << patsrc->bufa().name
				     << " because it is less than 4 characters long" << endl;
			}
			this->done = true;
			sinkPt_->finishRead(*patsrc_, true, true);
			return;
		}
		driver_->setQuery(patsrc, NULL);
		this->done = driver_->done;
		doneFirst_ = false;
		// Reset #-backtrack countdown
		if(btCnt_ != NULL) *btCnt_ = maxBts_;
		if(sinkPt_->setHits(*patsrc->bufa().hitset)) {
			this->done = true;
			sinkPt_->finishRead(*patsrc_, true, true);
		}
		// Grab a bit from the pseudo-random seed to determine whether
		// to start with forward or reverse complement
		firstIsFw_ = ((patsrc->bufa().seed & 0x10) == 0);
		chase_ = false;
	}

	/**
	 * Helper for reporting an alignment.
	 */
	inline bool report(const Range& ra,
	                   uint32_t first,
	                   uint32_t second,
	                   uint32_t tlen)
	{
		bool ebwtFw = ra.ebwt()->fw();
		assert_eq(bufa_->color, gColor);
		Hit hit;
		return params_->reportHit(
				ra.fw() ? (ebwtFw? bufa_->patFw    : bufa_->patFwRev) :
				          (ebwtFw? bufa_->patRc    : bufa_->patRcRev),
				ra.fw() ? (ebwtFw? &bufa_->qual    : &bufa_->qualRev) :
				          (ebwtFw? &bufa_->qualRev : &bufa_->qual),
				&bufa_->name,
				refs_,
				ra.ebwt()->rmap(),
				dec_,
				rand_,
				ra.fw(),
				ebwtFw,
				ra.edits(),
				make_pair(first, second), // position
				make_pair(ra.top(), ra.bot()),// arrows
				tlen,                     // textlen
				alen_,                    // qlen
				ra.stratum(),             // alignment stratum
				ra.cost(),                // cost, including qual penalty
				ra.bot() - ra.top() - 1,  // # other hits
				patsrc_->patid(),         // pattern id
				bufa_->seed,              // pseudo-random seed
				0,                        // mate (0 = unpaired)
				hit
				ASSERT_ONLY(, tmp_destU32_));
	}

	/**
	 * Advance the aligner.  Return true iff we're
	 * done with this read.
	 */
	virtual bool advance() {
		assert(!this->done);
		if(chase_) {
			assert(!gRangeMode);
			assert(driver_->foundRange);
			assert(!sinkPt_->irrelevantCost(driver_->range().cost()));
			if(!rchase_->foundOff() && !rchase_->done) {
				rchase_->advance();
				return false;
			}
			if(rchase_->foundOff()) {
				this->done = report(driver_->range(), rchase_->off().first,
				                    rchase_->off().second, rchase_->tlen());
				rchase_->reset();
			} else {
				assert(rchase_->done);
				// Forget this range; keep looking for ranges
				chase_ = false;
				driver_->foundRange = false;
				this->done = driver_->done;
			}
		}
		// Still advancing a
		if(!this->done && !chase_) {
			assert(!driver_->done || driver_->foundRange);
			if(driver_->foundRange) {
				const Range& ra = driver_->range();
				assert(!sinkPt_->irrelevantCost(ra.cost()));
				assert(ra.repOk());
				if(gRangeMode) {
					this->done = report(ra, ra.top(), ra.bot(), 0);
					driver_->foundRange = false;
				} else {
					rchase_->setTopBot(ra.top(), ra.bot(), ra.rlen(), rand_, ra.ebwt());
					if(rchase_->foundOff()) {
						this->done = report(
								ra, rchase_->off().first,
								rchase_->off().second, rchase_->tlen());
						rchase_->reset();
					}
					if(!rchase_->done && !sinkPt_->irrelevantCost(ra.cost())) {
						// Keep chasing this range
						chase_ = true;
					} else {
						driver_->foundRange = false;
					}
				}
			} else {
				this->done = sinkPt_->irrelevantCost(driver_->minCost);
				if(!this->done) {
					driver_->advance(ADV_COST_CHANGES);
				} else {
					// No longer necessarily true with chain input
					//assert(!sinkPt_->spanStrata());
				}
			}
			if(driver_->done && !driver_->foundRange && !chase_) {
				this->done = true;
			}
		}
		if(this->done) {
			sinkPt_->finishRead(*patsrc_, true, true);
		}
		return this->done;
	}

protected:

	// Reference sequences (needed for colorspace decoding)
	const BitPairReference* refs_;

	// Progress state
	bool doneFirst_;
	bool firstIsFw_;
	bool chase_;

	// Temporary HitSink; to be deleted
	const HitSinkPerThreadFactory& sinkPtFactory_;
	HitSinkPerThread* sinkPt_;

	// State for alignment
	EbwtSearchParams* params_;

	// State for getting alignments from ranges statefully
	RangeChaser* rchase_;

	// Range-finding state
	TDriver* driver_;

	const int maxBts_;
	ChunkPool *pool_;
	int *btCnt_;
	AlignerMetrics *metrics_;
	ColorspaceDecoder * dec_;
	
	ASSERT_ONLY(SStringExpandable<uint32_t> tmp_destU32_);
};

/**
 * Helper struct that holds a Range together with the coordinates where it al
 */
struct RangeWithCoords {
	Range r;
	U32Pair h;
};

/**
 * An aligner for finding paired alignments while operating entirely
 * within the Burrows-Wheeler domain.
 */
template<typename TRangeSource>
class PairedBWAlignerV2 : public Aligner {

	typedef std::pair<uint32_t,uint32_t> U32Pair;
	typedef EList<U32Pair> U32PairVec;
	typedef EList<Range> TRangeVec;
	typedef RangeSourceDriver<TRangeSource> TDriver;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:
	PairedBWAlignerV2(
		EbwtSearchParams* params,
		EbwtSearchParams* paramsSe1,
		EbwtSearchParams* paramsSe2,
		TDriver* driver,
		RefAligner<SString<char> >* refAligner,
		RangeChaser* rchase,
		HitSink& sink,
		const HitSinkPerThreadFactory& sinkPtFactory,
		HitSinkPerThread* sinkPt,
		HitSinkPerThread* sinkPtSe1,
		HitSinkPerThread* sinkPtSe2,
		uint32_t mixedAttemptLim,
		const BitPairReference* refs,
		int maxBts,
		ChunkPool *pool,
		ColorspaceDecoder * dec,
		int *btCnt) :
		Aligner(true),
		refs_(refs),
		patsrc_(NULL),
		qlen1_(0), qlen2_(0),
		chase_(false),
		donePe_(false),
		doneSe1_(false),
		doneSe2_(false),
		refAligner_(refAligner),
		sinkPtFactory_(sinkPtFactory),
		sinkPt_(sinkPt),
		sinkPtSe1_(sinkPtSe1),
		sinkPtSe2_(sinkPtSe2),
		params_(params),
		paramsSe1_(paramsSe1),
		paramsSe2_(paramsSe2),
		mixedAttemptLim_(mixedAttemptLim),
		mixedAttempts_(0),
		rchase_(rchase),
		driver_(driver),
		pool_(pool),
		maxBts_(maxBts),
		btCnt_(btCnt),
		dec_(dec)
	{
		assert(sinkPt_ != NULL);
		assert(params_ != NULL);
		assert(driver_ != NULL);
		assert(dec_    != NULL);
	}

	virtual ~PairedBWAlignerV2() {
		delete driver_; driver_ = NULL;
		delete params_; params_ = NULL;
		if(paramsSe1_ != NULL) {
			delete paramsSe1_; paramsSe1_ = NULL;
			delete paramsSe2_; paramsSe2_ = NULL;
		}
		delete rchase_; rchase_ = NULL;
		delete[] btCnt_; btCnt_ = NULL;
		delete refAligner_; refAligner_ = NULL;
		delete dec_; dec_ = NULL;
		sinkPtFactory_.destroy(sinkPt_); sinkPt_ = NULL;
		if(sinkPtSe1_ != NULL) {
			sinkPtFactory_.destroy(sinkPtSe1_); sinkPtSe1_ = NULL;
			sinkPtFactory_.destroy(sinkPtSe2_); sinkPtSe2_ = NULL;
		}
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc) {
		assert(!patsrc->bufa().empty());
		Aligner::setQuery(patsrc); // set fields & random seed
		assert(!patsrc->bufb().empty());
		// Give all of the drivers pointers to the relevant read info
		patsrc_ = patsrc;
		pool_->reset(&patsrc->bufa().name, patsrc->patid());
		if(patsrc->bufa().length() < 4 || patsrc->bufb().length() < 4) {
			if(!gQuiet) {
				cerr << "Warning: Skipping pair " << patsrc->bufa().name
				     << " because a mate is less than 4 characters long" << endl;
			}
			this->done = true;
			sinkPt_->finishRead(*patsrc_, true, true);
			return;
		}
		driver_->setQuery(patsrc, NULL);
		qlen1_ = (uint32_t)patsrc_->bufa().length();
		qlen2_ = (uint32_t)patsrc_->bufb().length();
		if(btCnt_ != NULL) (*btCnt_) = maxBts_;
		mixedAttempts_ = 0;
		// Neither orientation is done
		this->done = false;
		// No ranges are being chased yet
		chase_ = false;
		donePe_ = doneSe1_ = doneSe2_ = false;
		pairs_fw_.clear();
		pairs_rc_.clear();
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 *
	 * A call to this function does one of many things:
	 * 1. Advance a RangeSourceDriver and check if it found a new range
	 * 2. Advance a RowChaseDriver and check if it found a reference
	 *    offset for a an alignment in a range
	 */
	virtual bool advance() {
		assert(!this->done);
		if(chase_) {
			assert(!gRangeMode); // chasing ranges
			if(!rchase_->foundOff() && !rchase_->done) {
				rchase_->advance();
				return false;
			}
			assert(rchase_->foundOff() || rchase_->done);
			if(rchase_->foundOff()) {
				const Range& r = driver_->range();
				assert(r.repOk());
				resolveOutstanding(
					rchase_->off(),
					r.ebwt()->plen()[rchase_->off().first], r);
				rchase_->reset();
			} else {
				assert(rchase_->done);
				// Forget this range; keep looking for ranges
				chase_ = false;
				this->done = driver_->done;
			}
		}

		if(!this->done && !chase_) {
			// Search for more ranges for whichever mate currently has
			// fewer candidate alignments
			if(!driver_->done) {
				if(!this->done) {
					//
					// Check whether any of the PE/SE possibilities
					// have become impossible due to the minCost
					//
					if(!donePe_) {
						assert(!this->done);
						donePe_ = sinkPt_->irrelevantCost(driver_->minCost);
						if(donePe_ && (!sinkPt_->empty() || sinkPtSe1_ == NULL)) {
							// Paired-end alignment(s) were found, no
							// more will be found, and no unpaired
							// alignments are requested, so stop
							this->done = true;
						}
						if(donePe_ && sinkPtSe1_ != NULL) {
							// Note: removeMate affects minCost
							if(doneSe1_) driver_->removeMate(1);
							if(doneSe2_) driver_->removeMate(2);
						}
					}
					if(!this->done && sinkPtSe1_ != NULL) {
						if(!doneSe1_) {
							doneSe1_ = sinkPtSe1_->irrelevantCost(driver_->minCost);
							if(doneSe1_ && donePe_) driver_->removeMate(1);
						}
						if(!doneSe2_) {
							doneSe2_ = sinkPtSe2_->irrelevantCost(driver_->minCost);
							if(doneSe2_ && donePe_) driver_->removeMate(2);
						}
						// Do Se1 again, because removing Se2 may have
						// nudged minCost over the threshold
						if(!doneSe1_) {
							doneSe1_ = sinkPtSe1_->irrelevantCost(driver_->minCost);
							if(doneSe1_ && donePe_) driver_->removeMate(1);
						}
						if(doneSe1_ && doneSe2_) assert(donePe_);
						this->done = donePe_ && doneSe1_ && doneSe2_;
					}

					if(!this->done) {
						if(sinkPtSe1_ != NULL) {
							assert(doneSe1_ || !sinkPtSe1_->irrelevantCost(driver_->minCost));
							assert(doneSe2_ || !sinkPtSe2_->irrelevantCost(driver_->minCost));
						}
						assert(donePe_ || !sinkPt_->irrelevantCost(driver_->minCost));
						driver_->advance(ADV_COST_CHANGES);
					}
				}
				if(driver_->foundRange) {
					// Use Burrows-Wheeler for this pair (as usual)
					chase_ = true;
					driver_->foundRange = false;
					const Range& r = driver_->range();
					assert(r.repOk());
					rchase_->setTopBot(r.top(), r.bot(), r.rlen(), rand_, r.ebwt());
				}
			} else {
				this->done = true;
			}
		}

		if(this->done) {
			bool reportedPe = (sinkPt_->finishRead(*patsrc_, true, true) > 0);
			if(sinkPtSe1_ != NULL) {
				sinkPtSe1_->finishRead(*patsrc_, !reportedPe, false);
				sinkPtSe2_->finishRead(*patsrc_, !reportedPe, false);
			}
		}
		return this->done;
	}

protected:

	/**
	 * Helper for reporting a pair of alignments.  As of now, we report
	 * a paired alignment by reporting two consecutive alignments, one
	 * for each mate.
	 */
	bool report(const Range& rL, // range for upstream mate
	            const Range& rR, // range for downstream mate
	            uint32_t first,  // ref idx
	            uint32_t upstreamOff, // offset for upstream mate
	            uint32_t dnstreamOff, // offset for downstream mate
	            uint32_t tlen, // length of ref
	            bool pairFw,   // whether the pair is being mapped to fw strand
	            bool ebwtFwL,
	            bool ebwtFwR,
	            const ReferenceMap *rmap)
	{
		assert_lt(upstreamOff, dnstreamOff);
		uint32_t spreadL = rL.bot() - rL.top();
		uint32_t spreadR = rR.bot() - rR.top();
		uint32_t oms = min(spreadL, spreadR) - 1;
		Read* bufL = pairFw ? bufa_ : bufb_;
		Read* bufR = pairFw ? bufb_ : bufa_;
		uint32_t lenL = pairFw ? alen_ : blen_;
		uint32_t lenR = pairFw ? blen_ : alen_;
		assert(!params_->sink().exceededOverThresh());
		assert_eq(bufL->color, gColor);
		// Print upstream mate first
		Hit hitL, hitR;
		assert_eq(bufR->color, gColor);
		return params_->reportHitPair(
				rL.fw() ? (ebwtFwL?  bufL->patFw  :  bufL->patFwRev) :
				          (ebwtFwL?  bufL->patRc  :  bufL->patRcRev),
				rL.fw() ? (ebwtFwL? &bufL->qual    : &bufL->qualRev) :
				          (ebwtFwL? &bufL->qualRev : &bufL->qual),
				&bufL->name,                  // read name
				rL.fw(),                      // fw?
				ebwtFwL,                      // index was fw?
				rL.edits(),                   // edits
				make_pair(first, upstreamOff),// position
				make_pair(rL.top(), rL.bot()),// arrows
				tlen,                         // textlen
				lenL,                         // qlen
				rL.stratum(),                 // alignment stratum
				rL.cost(),                    // cost, including quality penalty
				oms,                          // # other hits
				bufL->patid,
				bufL->seed,
				pairFw ? 1 : 2,
				rR.fw() ? (ebwtFwR?  bufR->patFw  :  bufR->patFwRev) :
				          (ebwtFwR?  bufR->patRc  :  bufR->patRcRev),
				rR.fw() ? (ebwtFwR? &bufR->qual    : &bufR->qualRev) :
				          (ebwtFwR? &bufR->qualRev : &bufR->qual),
				&bufR->name,                  // read name
				rR.fw(),                      // fw?
				ebwtFwR,                      // index was fw?
				rR.edits(),                   // edits
				make_pair(first, dnstreamOff),// position
				make_pair(rR.top(), rR.bot()),// arrows
				tlen,                         // textlen
				lenR,                         // qlen
				rR.stratum(),                 // alignment stratum
				rR.cost(),                    // cost, including quality penalty
				oms,                          // # other hits
				bufR->patid,
				bufR->seed,
				pairFw ? 2 : 1,
				refs_,                        // BitPairReference
				rmap,                         // ReferenceMap
				dec_,                         // ColorspaceDecoder
				rand_,                        // Pseudo-random generator
				hitL,
				hitR
				ASSERT_ONLY(, tmp_destU32_));
	}

	/**
	 * Helper for reporting a pair of alignments.  As of now, we report
	 * a paired alignment by reporting two consecutive alignments, one
	 * for each mate.
	 */
	void reportSe(const Range& r, U32Pair h, uint32_t tlen) {
		EbwtSearchParams*params = (r.mate1() ? paramsSe1_ : paramsSe2_);
		assert(!(r.mate1() ? doneSe1_ : doneSe2_));
		Read* buf = r.mate1() ? bufa_ : bufb_;
		bool ebwtFw = r.ebwt()->fw();
		uint32_t len = r.mate1() ? alen_ : blen_;
		assert_eq(buf->color, gColor);
		// Print upstream mate first
		Hit hit;
		if(params->reportHit(
			r.fw() ? (ebwtFw?  buf->patFw   :  buf->patFwRev) :
			         (ebwtFw?  buf->patRc   :  buf->patRcRev),
			r.fw() ? (ebwtFw? &buf->qual    : &buf->qualRev) :
			         (ebwtFw? &buf->qualRev : &buf->qual),
			&buf->name,
			refs_,
			r.ebwt()->rmap(),
			dec_,
			rand_,
			r.fw(),
			ebwtFw,
			r.edits(),               // edits
			h,                       // position
			make_pair(r.top(), r.bot()), // arrows
			tlen,                    // textlen
			len,                     // qlen
			r.stratum(),             // alignment stratum
			r.cost(),                // cost, including quality penalty
			r.bot() - r.top() - 1,   // # other hits
			buf->patid,
			buf->seed,
			0,
			hit
			ASSERT_ONLY(, tmp_destU32_)))
		{
			if(r.mate1()) doneSe1_ = true;
			else          doneSe2_ = true;
			if(donePe_)   driver_->removeMate(r.mate1() ? 1 : 2);
		}
	}

	void resolveOutstanding(const U32Pair& off,
	                        const uint32_t tlen,
	                        const Range& range)
	{
		assert(!this->done);
		if(!donePe_) {
			bool ret = resolveOutstandingInRef(off, tlen, range);
			if(++mixedAttempts_ > mixedAttemptLim_ || ret) {
				// Give up on this pair
				donePe_ = true;
				if(sinkPtSe1_ != NULL) {
					if(doneSe1_) driver_->removeMate(1);
					if(doneSe2_) driver_->removeMate(2);
				}
			}
			this->done = (donePe_ && (!sinkPt_->empty() || sinkPtSe1_ == NULL || (doneSe1_ && doneSe2_)));
		}
		if(!this->done && sinkPtSe1_ != NULL) {
			bool doneSe = (range.mate1() ? doneSe1_ : doneSe2_);
			if(!doneSe) {
				// Hold onto this single-end alignment in case we don't
				// find any paired alignments
				reportSe(range, off, tlen);
			}
			this->done = doneSe1_ && doneSe2_ && donePe_;
		}
	}

	/**
	 * Given a vector of reference positions where one of the two mates
	 * (the "anchor" mate) has aligned, look directly at the reference
	 * sequence for instances where the other mate (the "outstanding"
	 * mate) aligns such that mating constraint is satisfied.
	 *
	 * This function picks up to 'pick' anchors at random from the
	 * 'offs' array.  It returns the number that it actually picked.
	 */
	bool resolveOutstandingInRef(const U32Pair& off,
	                             const uint32_t tlen,
	                             const Range& range)
	{
		assert(!donePe_);
		assert(refs_->loaded());
		assert_lt(off.first, refs_->numRefs());
		// pairFw = true if the anchor indicates that the pair will
		// align in its forward orientation (i.e. with mate1 to the
		// left of mate2)
		bool pairFw = (range.mate1())? (range.fw() == gMate1fw) : (range.fw() == gMate2fw);
		// matchRight = true, if the opposite mate will be to the right
		// of the anchor mate
		bool matchRight = (pairFw ? range.mate1() : !range.mate1());
		// fw = orientation of the opposite mate
		bool fw = range.mate1() ? gMate2fw : gMate1fw; // whether outstanding mate is fw/rc
		if(!pairFw) fw = !fw;
		// 'seq' = sequence for opposite mate
		const BTDnaString& seq  =
			fw ? (range.mate1() ? patsrc_->bufb().patFw   :
		                          patsrc_->bufa().patFw)  :
		         (range.mate1() ? patsrc_->bufb().patRc   :
		                          patsrc_->bufa().patRc);
		// 'qual' = qualities for opposite mate
		const BTString& qual =
			fw ? (range.mate1() ? patsrc_->bufb().qual  :
			                      patsrc_->bufa().qual) :
			     (range.mate1() ? patsrc_->bufb().qualRev :
			                      patsrc_->bufa().qualRev);
		uint32_t qlen = (uint32_t)seq.length();  // length of outstanding mate
		uint32_t alen = (uint32_t)(range.mate1() ? patsrc_->bufa().length() :
		                                           patsrc_->bufb().length());
		int minins = gMinInsert;
		int maxins = gMaxInsert;
		if(gMate1fw) {
			minins = max<int>(0, minins - patsrc_->bufa().trimmed5);
			maxins = max<int>(0, maxins - patsrc_->bufa().trimmed5);
		} else {
			minins = max<int>(0, minins - patsrc_->bufa().trimmed3);
			maxins = max<int>(0, maxins - patsrc_->bufa().trimmed3);
		}
		if(gMate2fw) {
			minins = max<int>(0, minins - patsrc_->bufb().trimmed3);
			maxins = max<int>(0, maxins - patsrc_->bufb().trimmed3);
		} else {
			minins = max<int>(0, minins - patsrc_->bufb().trimmed5);
			maxins = max<int>(0, maxins - patsrc_->bufb().trimmed5);
		}
		assert_geq(minins, 0);
		assert_geq(maxins, 0);
		// Don't even try if either of the mates is longer than the
		// maximum insert size.
		if((uint32_t)maxins <= max(qlen, alen)) {
			return false;
		}
		const uint32_t tidx = off.first;  // text id where anchor mate hit
		const uint32_t toff = off.second; // offset where anchor mate hit
		// Set begin/end to the range of reference positions where
		// outstanding mate may align while fulfilling insert-length
		// constraints.
		uint32_t begin, end;
		assert_geq(maxins, minins);
		uint32_t insDiff = maxins - minins;
		if(matchRight) {
			end = toff + maxins;
			begin = toff + 1;
			if(qlen < alen) begin += alen-qlen;
			if(end > insDiff + qlen) {
				begin = max<uint32_t>(begin, end - insDiff - qlen);
			}
			end = min<uint32_t>(refs_->approxLen(tidx), end);
			begin = min<uint32_t>(refs_->approxLen(tidx), begin);
		} else {
			if(toff + alen < (uint32_t)maxins) {
				begin = 0;
			} else {
				begin = toff + alen - maxins;
			}
			uint32_t mi = min<uint32_t>(alen, qlen);
			end = toff + mi - 1;
			end = min<uint32_t>(end, toff + alen - minins + qlen);
			if(toff + alen + qlen < (uint32_t)(minins)) end = 0;
		}
		// Check if there's not enough space in the range to fit an
		// alignment for the outstanding mate.
		if(end - begin < qlen) return false;
		results_.clear();
		refAligner_->find(1, tidx, refs_, seq, qual, begin, end, results_,
		                  pairFw ? &pairs_fw_ : &pairs_rc_, toff, fw);
		for(size_t i = 0; i < results_.size(); i++) {
			RefAlignerHit& r = results_[i];
			assert(Range::finalized(r.edits));
			if(!fw) {
				// Reverse pos's in the r.edits so so that they're
				// w/r/t the 5' end
				Edit::invertPoss(r.edits, qlen);
				assert(Range::finalized(r.edits));
			}
			Range rn;
			rn.init(range.top(), range.bot(), r.cost, r.stratum, qlen,
			        fw, r.edits, range.ebwt());
			bool ebwtLFw = matchRight ? range.ebwt()->fw() : true;
			bool ebwtRFw = matchRight ? true : range.ebwt()->fw();
			if(report(
				matchRight ? range : rn, // range for upstream mate
				matchRight ? rn : range, // range for downstream mate
				tidx,                    // ref idx
				matchRight ? toff : r.off, // upstream offset
				matchRight ? r.off : toff, // downstream offset
				tlen,       // length of ref
				pairFw,     // whether the pair is being mapped to fw strand
				ebwtLFw,
				ebwtRFw,
				range.ebwt()->rmap())) return true;
		}
		return false;
	}

	const BitPairReference* refs_;

	PatternSourcePerThread *patsrc_;
	uint32_t qlen1_, qlen2_;
	bool chase_;

	// true -> we're no longer shooting for paired-end alignments;
	// just collecting single-end ones
	bool donePe_, doneSe1_, doneSe2_;

	// For searching for outstanding mates
	RefAligner<SString<char> >* refAligner_;

	// Temporary HitSink; to be deleted
	const HitSinkPerThreadFactory& sinkPtFactory_;
	HitSinkPerThread* sinkPt_;
	HitSinkPerThread* sinkPtSe1_, * sinkPtSe2_;

	// State for alignment
	EbwtSearchParams* params_;
	// for single-end:
	EbwtSearchParams* paramsSe1_, * paramsSe2_;

	const uint32_t mixedAttemptLim_;
	uint32_t mixedAttempts_;

	// State for getting alignments from ranges statefully
	RangeChaser* rchase_;

	// Range-finding state for first mate
	TDriver* driver_;

	// Pool for distributing chunks of best-first path descriptor memory
	ChunkPool *pool_;

	int maxBts_; // maximum allowed # backtracks
	int *btCnt_; // current backtrack count

	/// For keeping track of paired alignments that have already been
	/// found for the forward and reverse-comp pair orientations
	TSetPairs pairs_fw_, pairs_rc_;

	EList<RefAlignerHit> results_;

	ColorspaceDecoder * dec_;

	ASSERT_ONLY(SStringExpandable<uint32_t> tmp_destU32_);
};

#endif /* ALIGNER_H_ */
