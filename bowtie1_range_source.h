#ifndef EBWT_SEARCH_BACKTRACK_H_
#define EBWT_SEARCH_BACKTRACK_H_

#include <stdexcept>
#include "pat.h"
#include "qual.h"
#include "range.h"
#include "range_source.h"
#include "aligner_metrics.h"
#include "search_globals.h"
#include "sstring.h"

/**
 * Class that coordinates quality- and quantity-aware backtracking over
 * some range of a read sequence.
 *
 * The creator can configure the BacktrackManager to treat different
 * stretches of the read differently.
 */
class EbwtRangeSource : public RangeSource {
	typedef std::pair<int, int> TIntPair;
public:
	EbwtRangeSource(
			const Ebwt* ebwt,
			bool         fw,
			uint32_t     qualLim,
			bool         reportExacts,
			int          halfAndHalf,
			bool         partial,
			AlignerMetrics *metrics = NULL) :
		RangeSource(),
		qry_(NULL),
		qlen_(0), rlen_(0),
		qual_(NULL),
		name_(NULL),
		altQry_(NULL),
		altQual_(NULL),
		alts_(0),
		ref_(NULL),
		ebwt_(ebwt),
		fw_(fw),
		qualLim_(qualLim),
		reportExacts_(reportExacts),
		halfAndHalf_(halfAndHalf),
		partial_(partial),
		qdepth5_(0), rdepth5_(0),
		qdepth3_(0), rdepth3_(0),
		skippingThisRead_(false),
		metrics_(metrics),
		loGapDepth_(0),
		hiGapDepth_(0)
	{
		curEbwt_ = ebwt_;
		memset(qdeps_, 0, sizeof(uint16_t)*4);
		memset(rdeps_, 0, sizeof(uint16_t)*4);
	}

	/**
	 * Set a new query read.
	 */
	virtual void setQuery(Read& r, Range *seedRange) {
		const bool ebwtFw = ebwt_->fw();
		if(ebwtFw) {
			qry_     = fw_ ? &r.patFw   : &r.patRc;
			qual_    = fw_ ? &r.qual    : &r.qualRev;
			altQry_  = fw_ ? r.altPatFw : r.altPatRc;
			altQual_ = fw_ ? r.altQual  : r.altQualRev;
		} else {
			qry_     = fw_ ? &r.patFwRev   : &r.patRcRev;
			qual_    = fw_ ? &r.qualRev    : &r.qual;
			altQry_  = fw_ ? r.altPatFwRev : r.altPatRcRev;
			altQual_ = fw_ ? r.altQualRev  : r.altQual;
		}
		patid_ = r.patid;
		name_ = &r.name;
		alts_ = r.alts;
		if(seedRange != NULL) seedRange_ = *seedRange;
		else                  seedRange_.invalidate();
		cmap_.clear();
		qryDelDeps_.clear();
		qlen_ = qry_->length();
		refBuf_.clear();
		// Apply edits from the partial alignment to the query pattern
		if(seedRange_.valid()) {
			assert(seedRange_.repOk());
			const EList<Edit>& es = seedRange_.edits();
			const size_t esz = es.size();
			assert_gt(esz, 0);
			size_t ecur = 0;
			for(int i = 0; i < (int)qlen_; i++) {
				bool mm = false, ins = false, del = false;
				int chrs = 0;
				cmap_.q2r.append((int)refBuf_.length());
				while(ecur < esz && (int)es[ecur].pos == i) {
					if(es[ecur].isMismatch() || es[ecur].isReadGap()) {
						assert(!del);
						refBuf_.appendChar(es[ecur].chr);
						cmap_.r2q.append(i);
						chrs++;
						if(es[ecur].isMismatch()) mm = true;
						else {
							assert(!mm);
							ins = true;
						}
					} else {
						assert(es[ecur].isRefGap());
						del = true;
					}
					ecur++;
				}
				if(!mm && !del) {
					refBuf_.append((*qry_)[qlen_ - i - 1]);
					cmap_.r2q.append(i);
					chrs++;
				}
				qryDelDeps_.append(del);
			}
			assert_eq(refBuf_.length(), cmap_.r2q.length());
			assert_eq(qlen_, cmap_.q2r.length());
			refBuf_.reverse();
			ref_ = &refBuf_;
		} else {
			// No partial alignment, so ref is the same as qry
			ref_ = qry_;
			assert_gt(ref_->length(), 0);
			for(size_t i = 0; i < qlen_; i++) {
				cmap_.q2r.append((int)i);
				cmap_.r2q.append((int)i);
				qryDelDeps_.append(false);
			}
		}
		// Make sure every qual is a valid qual ASCII character (>= 33)
		for(size_t i = 0; i < qual_->length(); i++) {
			assert_geq((*qual_)[i], 33);
			for(int j = 0; j < alts_; j++) {
				assert_geq(altQual_[j][i], 33);
			}
		}
		rlen_ = ref_->length();
		if(rlen_ != qlen_) {
			assert(seedRange_.valid());
			for(int i = 0; i < 4; i++) {
				rdeps_[i] = qdeps_[i] + (rlen_ - qlen_);
			}
			rdepth3_ = qdepth3_ + (uint32_t)(rlen_ - qlen_);
			rdepth5_ = qdepth5_ + (uint32_t)(rlen_ - qlen_);
		}
		assert_geq(qual_->length(), qlen_);
		this->done = false;
		this->foundRange = false;
		skippingThisRead_ = false;
		rand_.init(r.seed);
		if(gGaps) {
			// alignment must be gapless as long as qdepth < loGapDepth_
			loGapDepth_ = gGapBarrier;
			loGapDepth_ = max<int>(loGapDepth_, qdeps_[0]);
			loGapDepth_ = min<int>((int)qlen_, 0);
			// alignment must be gapless as long as qdepth >= loGapDepth_
			hiGapDepth_ = (int)(qry_->length() - gGapBarrier + 1);
			hiGapDepth_ = max<int>((int)hiGapDepth_, 0);
			hiGapDepth_ = min<int>((int)hiGapDepth_, (int)qlen_);
			// Calculate the largest possible number of characters that the
			// reference side of the alignment can have, given
			size_t pen = min<size_t>(gInsOpen, gInsExtend);
			pen = min<size_t>(pen, gDelOpen);
			pen = min<size_t>(pen, gDelExtend);
			maxGaps_ = 0;
			if(qualLim_ > 0 && pen <= qualLim_) {
				maxGaps_ = (int)(qualLim_ / pen);
			}
		} else {
			loGapDepth_ = hiGapDepth_ = maxGaps_ = 0;
		}
		maxEdits_ = (int)(qlen_ + maxGaps_);
	}

	/**
	 * Set backtracking constraints.  Called by
	 * EbwtRangeSourceDriver::initRangeSource(), which is usually
	 * called right after setQuery().
	 */
	void setOffs(uint32_t depth5,   // depth of far edge of hi-half
	             uint32_t depth3,   // depth of far edge of lo-half
	             uint32_t unrevOff, // depth above which we cannot backtrack
	             uint32_t revOff1,  // depth above which we may backtrack just once
	             uint32_t revOff2,  // depth above which we may backtrack just twice
	             uint32_t revOff3)  // depth above which we may backtrack just three times
	{
		qdepth5_ = rdepth5_ = depth5;
		qdepth3_ = rdepth3_ = depth3;
		assert_geq(qdepth3_, qdepth5_);
		qdeps_[0] = rdeps_[0] = unrevOff;
		qdeps_[1] = rdeps_[1] = revOff1;
		qdeps_[2] = rdeps_[2] = revOff2;
		qdeps_[3] = rdeps_[3] = revOff3;
		if(rdeps_[3] == qlen_) {
			// The whole read is governed by the seed constraints, so
			// we can calculate a tighter maxEdits_
			maxEdits_ = 0;
			if(qdeps_[0] != qdeps_[1]) maxEdits_ = 1;
			if(qdeps_[1] != qdeps_[2]) maxEdits_ = 2;
			if(qdeps_[2] != qdeps_[3]) maxEdits_ = 3;
			maxGaps_ = min<int>(maxGaps_, maxEdits_);
		}
		if(seedRange_.valid()) {
#ifndef NDEBUG
			for(int i = 0; i < 3; i++) {
				assert_eq(qdeps_[i], qdeps_[i+1]);
			}
			assert_eq(qdepth3_, qdeps_[3]);
#endif
			if(rlen_ != qlen_) {
				for(int i = 0; i < 4; i++) {
					rdeps_[i] += (rlen_ - qlen_);
				}
				rdepth3_ = qdepth3_ + (uint32_t)(rlen_ - qlen_);
				rdepth5_ = qdepth5_ + (uint32_t)(rlen_ - qlen_);
			}
		}
	}

	/**
	 * Return true iff this RangeSource is allowed to report exact
	 * alignments (exact = no edits).
	 */
	bool reportExacts() const {
		return reportExacts_;
	}

	/// Return the current range
	virtual Range& range() {
		return curRange_;
	}

	/**
	 * Set qlen_ according to parameter, except don't let it fall below
	 * the length of the query.
	 */
	void setQlen(uint32_t qlen) {
		assert(qry_ != NULL);
		assert(!seedRange_.valid());
		qlen_ = rlen_ = qlen;
	}

	/**
	 * Install a "root branch" (bad metaphor) into the PathManager so
	 * that the next advance() begins
	 * a new search.  Note that contMan is empty upon return if there
	 * are no valid continuations to begin with.  Also note that
	 * calling initConts() may result in finding a range (i.e., if we
	 * immediately jump to a valid range using the ftab).
	 */
	virtual void
	initBranch(PathManager& pm) {
		assert(curEbwt_ != NULL);
		assert_gt(qry_->length(), 0);
		assert_leq(qlen_, qry_->length());
		assert_leq(rlen_, ref_->length());
		assert_geq(qual_->length(), qry_->length());
		const Ebwt& ebwt = *ebwt_;
		int ftabChars = ebwt._eh._ftabChars;
		this->foundRange = false;
		int nsInSeed = 0; int nsInFtab = 0;
		ASSERT_ONLY(allTops_.clear());
		if(skippingThisRead_) {
			this->done = true;
			return;
		}
		if(qlen_ < 4) {
			uint32_t maxmms = 0;
			if(qdeps_[0] != qdeps_[1]) maxmms = 1;
			if(qdeps_[1] != qdeps_[2]) maxmms = 2;
			if(qdeps_[2] != qdeps_[3]) maxmms = 3;
			if(qlen_ <= maxmms) {
				if(!gQuiet) {
					cerr << "Warning: Read (" << (*name_) << ") is less than "
					     << (maxmms+1) << " characters long; skipping..." << endl;
				}
				this->done = true;
				skippingThisRead_ = true;
				return;
			}
		}
		if(!tallyNs(nsInSeed, nsInFtab)) {
			// No alignments are possible because of the distribution
			// of Ns in the read in combination with the backtracking
			// constraints.
			assert(!seedRange_.valid()); // shouldn't be any Ns in seed
			return;
		}
		// icost = total cost penalty (major bits = stratum, minor bits =
		// quality penalty) incurred so far by partial alignment
		uint16_t icost = (seedRange_.valid()) ? seedRange_.cost() : 0;
		// iham = total quality penalty incurred so far by partial alignment
		uint16_t iham = (seedRange_.valid()) ? (seedRange_.cost() & ~0xc000): 0;
		assert_leq(iham, qualLim_);
		// m = depth beyond which ftab must not extend or else we might
		// miss some legitimate paths
		uint32_t m = min<uint32_t>(rdeps_[0], (uint32_t)qlen_);
		// Let skipInvalidExact = true if using the ftab would be a
		// waste because it would jump directly to an alignment we
		// couldn't use.
		bool ftabSkipsToEnd = (rlen_ == (uint32_t)ftabChars);
		bool skipInvalidExact = (!reportExacts_ && ftabSkipsToEnd);

		// If it's OK to use the ftab...
		if(nsInFtab == 0 && m >= (uint32_t)ftabChars && !skipInvalidExact) {
			// Use the ftab to jump 'ftabChars' chars into the read
			// from the right
			uint32_t ftabOff = calcFtabOff();
			uint32_t top = ebwt.ftabHi(ftabOff);
			uint32_t bot = ebwt.ftabLo(ftabOff+1);
			if(qlen_ == (uint32_t)ftabChars && bot > top) {
				// We found a range with 0 mismatches immediately.  Set
				// fields to indicate we found a range.
				assert(reportExacts_);
				edits_.clear();
				addPartialEdits();
				curRange_.init(top, bot, icost, (icost >> 14), (uint32_t)qlen_,
				               fw_, edits_, ebwt_);
				// Lump in the edits from the partial alignment
				assert(curRange_.repOk());
				this->foundRange  = true;
				//this->done = true;
				return;
			} else if (bot > top) {
				// We have a range to extend
				assert_leq(top, ebwt._eh._len);
				assert_leq(bot, ebwt._eh._len);
				// Create the root of the branch tree
				Branch *b = pm.bpool.alloc();
				if(b == NULL) {
					assert(pm.empty());
					return;
				}
				if(!b->init(
					pm.rpool, visit_, pm.bpool.lastId(), NULL, 0,
					Edit(), 0, (uint32_t)rlen_, rdeps_, 0, ftabChars, icost,
					iham, top, bot, ebwt._eh, ebwt.ebwt()))
				{
					// Negative result from b->init() indicates we ran
					// out of best-first chunk memory
					assert(pm.empty());
					return;
				}
				assert(!b->curtailed_);
				assert(!b->exhausted_);
				assert_gt(b->deps_[3], 0);
				visit_.reset(loGapDepth_, hiGapDepth_, maxGaps_);
				pm.enqueue(b); // insert into priority queue
				assert(!pm.empty());
			} else {
				// The arrows are already closed within the
				// unrevisitable region; give up
			}
		} else {
			// We can't use the ftab, so we start from the rightmost
			// position and use _fchr
			Branch *b = pm.bpool.alloc();
			if(b == NULL) {
				assert(pm.empty());
				return;
			}
			if(!b->init(pm.rpool, visit_, pm.bpool.lastId(), NULL, 0,
			            Edit(), 0, (uint32_t)rlen_, rdeps_, 0, 0, icost, iham, 0,
			            ebwt.fchr()[4], ebwt._eh, ebwt.ebwt()))
			{
				// Negative result from b->init() indicates we ran
				// out of best-first chunk memory
				assert(pm.empty());
				return;
			}
			assert(!b->curtailed_);
			assert(!b->exhausted_);
			assert_gt(b->deps_[3], 0);
			visit_.reset(loGapDepth_, hiGapDepth_, maxGaps_);
			pm.enqueue(b); // insert into priority queue
			assert(!pm.empty());
		}
		return;
	}

	/**
	 * Advance along the lowest-cost branch managed by the given
	 * PathManager.  Keep advancing until condition 'until' is
	 * satisfied.  Typically, the stopping condition 'until' is
	 * set to stop whenever pm's minCost changes.
	 */
	virtual void
	advanceBranch(int until, uint16_t minCost, PathManager& pm) {
		assert(curEbwt_ != NULL);
		this->foundRange = false;
		assert_gt(ref_->length(), 0);
		assert_gt(qry_->length(), 0);
		assert_leq(qlen_, qry_->length());
		assert_geq(qual_->length(), qry_->length());
		assert(!pm.empty());
		do {
			// Get the highest-priority branch according to the priority
			// queue in 'pm'
			Branch* br = pm.front(true);
			// Shouldn't be curtailed or exhausted
			assert(!br->exhausted_);
			assert(!br->curtailed_);
			assert_gt(br->deps_[3], 0);
			assert_leq(br->ham_, qualLim_);
			if(gVerbose) {
				br->print((*ref_), minCost, cout, (halfAndHalf_>0),
				          partial_ ? (uint32_t)qlen_ : 0, fw_, ebwt_->fw());
				br->printEdits((uint32_t)qlen_, ebwt_->fw());
			}
			assert(br->repOk((uint32_t)rlen_));
			assert_geq(qualLim_, br->ham_);

			ASSERT_ONLY(int stratum = br->cost_ >> 14); // shift the stratum over
			assert_lt(stratum, 4);
			int rdepth = br->tipDepth(); // depth into reference string
			int qdepth = (rdepth < (int)rlen_ ? cmap_.r2q[rdepth] : (uint32_t)qlen_);
			if(br->edit_.initialized()) {
				assert_gt(br->numEdits_, 0);
				assert(br->edit_.isReadGap() || (int)br->edit_.pos < rdepth);
			}
			assert(seedRange_.valid() || rdepth == qdepth);
			const Ebwt& ebwt = *ebwt_;
			if(halfAndHalf_ > 0) {
				assert(!seedRange_.valid());
				assert_gt(qdepth3_, qdepth5_);
				assert_gt(rdepth3_, rdepth5_);
			}
			bool reportedPartial = false;
			bool invalidExact = false;
			bool empty = false;
			bool hit = false;
			uint16_t cost = br->cost_;
			int rcur = 0, qcur = 0;
			bool hasEdits = false;
			bool noInsert = false;
			bool delExtend = false;
			// How far are we from the near end of the read?
			int endoff = qdepth;
			if(qry_->length() > qlen_) {
				// The near end may actually be partway in
				endoff += qry_->length() - qlen_;
			}
			if(halfAndHalf_ && !hhCheckTop(br, rdepth, 0)) {
				if(rdepth < (int)rlen_) {
					br->backPos()->eliminated_ = true;
				}
				assert(!seedRange_.valid());
				// Stop extending this branch because it violates a half-
				// and-half constraint
				if(metrics_ != NULL) metrics_->curBacktracks_++;
				pm.curtail(br, rdepth3_);
				goto bail;
			}
			rcur = (int)(rlen_ - rdepth - 1); // current offset into ref_
			qcur = (int)(qlen_ - qdepth - 1); // current offset into qry_
			assert_leq(rdepth, (int)ref_->length());
			assert_leq(rcur, (int)ref_->length());
			assert_leq(qdepth, (int)qry_->length());
			assert_leq(qcur, (int)qry_->length());
			assert_geq(rcur, -1); assert_geq(qcur, -1);
			assert(seedRange_.valid() || rcur == qcur);
			if(visit_.isJoinCell(qdepth, qdeps_[0])) {
				if(visit_.isVisited(qdepth, qdeps_[0], br->gapOff_, br->top_)) {
					assert(pm.triedGap());
					// The branch couldn't be extended further
					if(metrics_ != NULL) metrics_->curBacktracks_++;
					if(br->len_ > 0) br->retract();
					pm.curtail(br, rdepth3_ /* seed len, for stratum calc */);
					goto bail;
				} else {
					visit_.setVisited(qdepth, qdeps_[0], br->gapOff_, br->top_);
				}
			}
			// There is a center insertion position covered by more
			// than one search direction.  Set noInsert = true if we're
			// searching in the reverse index to make sure it's not
			// covered twice.
			noInsert = (rdepth == (int)rdepth5_ && !ebwt.fw());
			if(rdepth < (int)rlen_) {
				// Determine whether ranges at this location are candidates
				// for backtracking
				int c = (*ref_)[rcur]; // get char at this position
				assert_range(0, 4, c);
				// If any uncalled base's penalty is still under
				// the ceiling, then this position is an alternative
				uint8_t q[4] = {'!', '!', '!', '!'};
				uint8_t bestq = 0;
				bool revisitable = (rdepth >= br->deps_[0]);
				if(revisitable) {
					bestq = penaltiesAt(qcur, q, alts_, *qual_, altQry_, altQual_);
				}

				// The current query position is a legit alternative if it a) is
				// not in the unrevisitable region, and b) its selection would
				// not necessarily cause the quality ceiling (if one exists) to
				// be exceeded
				bool curIsAlternative =
					revisitable && (br->ham_ + bestq <= qualLim_);
				uint32_t otop = br->top_;

				// If c is 'N', then it's a mismatch
				if(c == 4 && rdepth > 0) {
					// Force the 'else if(curIsAlternative)' or 'else'
					// branches below
					br->top_ = br->bot_ = 1;
				} else if(c == 4) {
					// We'll take the 'if(br->top == 0 && br->bot == 0)'
					// branch below
					br->top_ = br->bot_ = 0;
				}

				// Get the state for this position
				Pos *rs = br->backPos();
				assert(rs != NULL);
				delExtend = qryDelDeps_[qdepth3_-1] && qdepth == (int)qdepth3_;

				if(br->justDeleted()) {
					ASSERT_ONLY(int r =)
					br->installRanges(
						c,                    // next character from the ref
						visit_,               // record of visited paths
						endoff,               // offset from near end of read
						(uint32_t)qry_->length(), // full length of read
						qdepth,               // depth into query region
						qdeps_[0],            // globally unrevisitable
					    qualLim_ - br->ham_,  // remaining quality budget
						q,                    // penalty for this mm
						noInsert,             // avoid redundant insertion
						delExtend);           // true -> deletion here extends
					assert(r < 9 || c == 4);
					// Update top and bot
					if(c < 4) {
						br->top_ = rs->tops[c];
						br->bot_ = rs->bots[c];
					} else {
						br->top_ = br->bot_ = 1;
					}
					assert(rs->repOk());
				} else {
					// Not just after a delete, in which case we have to
					// calculate tops and bots
					assert_eq(0, rs->tops[0]); assert_eq(0, rs->bots[0]);
					assert_eq(0, rs->tops[1]); assert_eq(0, rs->bots[1]);
					assert_eq(0, rs->tops[2]); assert_eq(0, rs->bots[2]);
					assert_eq(0, rs->tops[3]); assert_eq(0, rs->bots[3]);
					// Calculate the ranges for this position
					if(br->tipDepth() == 0 && !br->edit_.initialized()) {
						assert_eq(0, br->itop_);
						assert_eq(ebwt.fchr()[4], br->ibot_);
						// Calculate first quartet of ranges using the _fchr[]
						// array
									  rs->tops[0] = ebwt.fchr()[0];
						rs->bots[0] = rs->tops[1] = ebwt.fchr()[1];
						rs->bots[1] = rs->tops[2] = ebwt.fchr()[2];
						rs->bots[2] = rs->tops[3] = ebwt.fchr()[3];
						rs->bots[3]               = ebwt.fchr()[4];
						ASSERT_ONLY(int r =)
						br->installRanges(
							c,               // next character from the read
							visit_,          // record of visited paths
							endoff,          // offset from near end of read
							(uint32_t)qry_->length(),  // full length of read
							qdepth,          // depth into query region
							qdeps_[0],       // globally unrevisitable
						    qualLim_ - br->ham_, // remaining quality budget
							q,               // penalty for this mm
							noInsert,        // avoid redundant insertion
							delExtend);      // true -> deletion here extends
						assert(r < 9 || c == 4);
						// Update top and bot
						if(c < 4) {
							br->top_ = rs->tops[c];
							br->bot_ = rs->bots[c];
						}
						assert(rs->repOk());
					} else if(curIsAlternative && (br->bot_ > br->top_ || c == 4)) {
						// Calculate next quartet of ranges.  We hope that the
						// appropriate cache lines are prefetched.
						assert(br->ltop_.valid());
									  rs->tops[0] =
						rs->bots[0] = rs->tops[1] =
						rs->bots[1] = rs->tops[2] =
						rs->bots[2] = rs->tops[3] =
						rs->bots[3]               = 0;
						if(metrics_ != NULL) metrics_->curBwtOps_++;
						if(br->lbot_.valid()) {
							// The range delimited by ltop_/lbot_ has size >1
							ebwt.mapLFEx(br->ltop_, br->lbot_, rs->tops, rs->bots);
						} else {
							// The range delimited by ltop_/lbot_ has size 1
							int cc = ebwt.mapLF1(otop, br->ltop_);
							br->top_ = otop;
							assert(cc == -1 || (cc >= 0 && cc < 4));
							if(cc >= 0) {
								assert_lt(cc, 4);
								rs->tops[cc] = br->top_;
								rs->bots[cc] = (br->top_ + 1);
							}
						}
						ASSERT_ONLY(int r =)
						br->installRanges(
							c,               // next character from the read
							visit_,          // record of visited paths
							endoff,          // offset from near end of read
							(uint32_t)qry_->length(),  // full length of read
							qdepth,          // depth into query region
							qdeps_[0],       // global unrevisitable
							qualLim_ - br->ham_, // remaining quality budget
							q,               // penalty for this mm
							noInsert,        // avoid redundant insertion
							delExtend);      // true -> deletion here extends
						assert(r < 9 || c == 4);
						// Update top and bot
						if(c < 4) {
							br->top_ = rs->tops[c];
							br->bot_ = rs->bots[c];
						} else {
							br->top_ = br->bot_ = 1;
						}
						assert(rs->repOk());
					} else if(br->bot_ > br->top_) {
						// This read position is not a legitimate backtracking
						// alternative.  No need to do the bookkeeping for the
						// entire quartet, just do c.  We hope that the
						// appropriate cache lines are prefetched before now;
						// otherwise, we're about to take an expensive cache
						// miss.
						assert(br->ltop_.valid());
						if(revisitable) {
							rs->eliminated_ = true;
							assert(br->eliminated(br->len_));
						}
						if(c < 4) {
							rs->flags.qchr = c;
							if(br->top_ + 1 == br->bot_) {
								if(metrics_ != NULL) metrics_->curBwtOps_++;
								rs->tops[c] = rs->bots[c] =
									br->bot_ = br->top_ =
										ebwt.mapLF1(br->top_, br->ltop_, c);
								if(br->bot_ != 0xffffffff) {
									br->bot_ = rs->bots[c] = br->bot_ + 1;
								}
							} else {
								if(metrics_ != NULL) metrics_->curBwtOps_++;
								rs->tops[c] = br->top_ = ebwt.mapLF(br->ltop_, c);
								assert(br->lbot_.valid());
								if(metrics_ != NULL) metrics_->curBwtOps_++;
								rs->bots[c] = br->bot_ = ebwt.mapLF(br->lbot_, c);
							}
						}
						if(revisitable) assert(rs->repOk());
					} else {
						// No outgoing paths
						rs->eliminated_ = true;
						if(revisitable) {
							assert(br->eliminated(br->len_));
							assert(rs->repOk());
						}
					}
				}
				assert_geq(br->bot_, br->top_);
				// br->top_ and br->bot_ now contain the next top and bot
			} else {
				// We've just finished processing every character in
				// the read.
				assert_eq((int)rlen_, rdepth);
				rcur = 0; // was -1
				assert_geq(br->bot_, br->top_);
			}
			assert_geq(br->bot_, br->top_);
			empty = (br->top_ == br->bot_);
			hit = (rcur == 0 && !empty);

			// Check whether we've obtained an exact alignment when
			// we've been instructed not to report exact alignments
			hasEdits = br->edit_.initialized();
			invalidExact = (hit && !hasEdits && !reportExacts_);
			assert_leq(br->ham_, qualLim_);

			// Set this to true if the only way to make legal progress
			// is via one or more additional backtracks.
			if(halfAndHalf_ && !hhCheck(br, rdepth, empty)) {
				// This alignment doesn't satisfy the half-and-half
				// requirements; reject it
				if(metrics_ != NULL) metrics_->curBacktracks_++;
				pm.curtail(br, rdepth3_);
				goto bail;
			}

			if(hit &&            // there is a range to report
			   !invalidExact &&  // not disqualified by no-exact-hits setting
			   !reportedPartial) // not an already-reported partial alignment
			{
				if(gVerbose) {
					if(partial_) cout << " Partial alignment:" << endl;
					else         cout << " Final alignment:"   << endl;
					br->len_++;
					br->print((*ref_), minCost, cout, halfAndHalf_ > 0,
					          partial_ ? (uint32_t)qlen_ : 0, fw_, ebwt_->fw());
					br->printEdits((uint32_t)qlen_, ebwt_->fw());
					br->len_--;
					cout << endl;
				}
				assert_gt(br->bot_, br->top_);
				assert_leq(br->ham_, qualLim_);
				int strat = br->cost_ >> 14;
				assert_leq((int)(br->cost_ & ~0xc000), qualLim_);
				if(metrics_ != NULL) metrics_->setReadHasRange();
				edits_.clear();
				Branch *cur = br;
				Edit e = cur->edit_;
				while(e.initialized()) {
					assert_neq(e.qchr, e.chr);
					if(rlen_ != qlen_) {
						// If we extended a partial alignment that had
						// net gaps != 0, transform back into read
						// coordinates
						assert(seedRange_.valid());
						ASSERT_ONLY(int rc = (*ref_)[rlen_ - e.pos - 1]);
						e.pos += (qlen_ - rlen_);
						ASSERT_ONLY(int qc = (*qry_)[qlen_ - e.pos - 1]);
						assert_eq(rc, qc);
					}
					if(e.qchr != '-') {
						assert_eq((char)e.qchr, qry_->toChar(qlen_ - e.pos - 1));
					}
					assert_lt(e.pos, qlen_);
					edits_.push_back(e);
					cur = cur->parent_;
					assert(cur != NULL);
					e = cur->edit_;
				}
				assert(cur->parent_ == NULL);
				// We added the edits in reverse order, so now we
				// re-reverse them.
				edits_.reverse();
				assert(Range::finalized(edits_));
				if(ebwt_->fw() == fw_) {
					// Re-re-reverse the edits while also changing
					// their .pos field to be w/r/t the 5' end.
					Edit::invertPoss(edits_, qlen_);
					assert(Range::finalized(edits_));
				}
				// Make e.pos be w/r/t 5' end
				addPartialEdits();
				curRange_.init(br->top_, br->bot_, br->cost_, strat,
				               (uint32_t)qlen_, fw_, edits_, ebwt_);
				// edits_ is now invalid, having been swapped out
				this->foundRange = true;
#ifndef NDEBUG
				if(gAllowRedundant <= 1) {
					int64_t top2 = (int64_t)br->top_;
					top2++; // ensure it's not 0
					if(ebwt_->fw()) top2 = -top2;
					int len = (int)(rlen_ + br->gapOff_);
					top2 ^= ((int64_t)len << 32);
					assert(allTops_.find(top2) == allTops_.end());
					allTops_.insert(top2);
				}
#endif
				assert(curRange_.repOk());
				// Must curtail because we've consumed the whole pattern
				if(metrics_ != NULL) metrics_->curBacktracks_++;
				pm.curtail(br, rdepth3_ /* seed len, for stratum calc */);
			} else if(empty || rcur == 0) {
				// The branch couldn't be extended further
				if(metrics_ != NULL) metrics_->curBacktracks_++;
				pm.curtail(br, rdepth3_ /* seed len, for stratum calc */);
			} else {
				// Extend the branch by one position; no change to its cost
				// so there's no need to reconsider where it lies in the
				// priority queue
				assert_neq(0, rcur);
				br->extend();
			}
		bail:
			// Make sure the front element of the priority queue is
			// extendable (i.e. not curtailed) and then prep it.
			if(!pm.splitAndPrep(visit_, cmap_, rand_, (uint32_t)rlen_, qdeps_[0],
			                    qualLim_, rdepth3_, ebwt_->_eh, ebwt_->ebwt()))
			{
				pm.reset(0);
				assert(pm.empty());
			}
			if(pm.empty()) {
				// No more branches
				break;
			}
			assert(!pm.front(false)->curtailed_);
			assert(!pm.front(false)->exhausted_);
			if(until == ADV_COST_CHANGES && pm.front(false)->cost_ != cost) break;
			else if(until == ADV_STEP) break;
		} while(!this->foundRange);
		if(!pm.empty()) {
			assert(!pm.front(false)->curtailed_);
			assert(!pm.front(false)->exhausted_);
		}
	}

	/**
	 * Return true iff we're enforcing a half-and-half constraint
	 * (forced edits in both seed halves).
	 */
	int halfAndHalf() const {
		return halfAndHalf_;
	}

protected:

	/**
	 * Lump all the seed-alignment edits from the seedRange_ range
	 * found previously to the curRange_ range just found.
	 */
	void addPartialEdits() {
		// Lump in the edits from the partial alignment; these have
		// already been adjusted such that .pos is w/r/t 5' end.
		if(seedRange_.valid()) {
			if(edits_.empty()) {
				edits_ = seedRange_.edits();
				assert(Range::finalized(edits_));
			} else {
				edits_.insert(seedRange_.edits(), 0);
				assert(Range::finalized(edits_));
			}
		}
	}

	/**
	 * Return true iff we're OK to continue after considering which
	 * half-seed boundary we're passing through, together with the
	 * number of mismatches accumulated so far.  Return false if we
	 * should stop because a half-and-half constraint is violated.  If
	 * we're not currently passing a half-seed boundary, just return
	 * true.
	 */
	bool hhCheck(Branch *b, uint32_t refdep, bool empty) {
		assert(!seedRange_.valid());
		const uint32_t nedits = b->numEdits_;
		ASSERT_ONLY(uint32_t lim3 = (rdeps_[3] == rdeps_[2])? 2 : 3);
		ASSERT_ONLY(uint32_t lim5 = (rdeps_[1] == rdeps_[0])? 2 : 1);
		if(b->edit_.isReadGap()) {
			// tipDepth() wasn't advanced when this topmost insertion
			// was added, so we have to add 1 to refdep to avoid
			// thinking that we just crossed the border, when in fact
			// we added an insert just after crossing the border.
			refdep++;
		}
		if((refdep == (rdepth5_-1)) && !empty) {
			// We're crossing the boundary separating the hi-half
			// from the lo-half
			// We should induce a mismatch if we haven't mismatched
			// yet, so that we don't waste time pursuing a match
			// that was covered by a previous phase
			assert_leq(nedits, lim5);
			return nedits > 0;
		} else if((refdep == (rdepth3_-1)) && !empty) {
			// We're crossing the boundary separating the lo-half
			// from the non-seed portion of the read
			assert_leq(nedits, lim3);
			assert_gt(nedits, 0);
			// Count the mismatches in the lo and hi halves
			uint32_t loHalfMms = 0, hiHalfMms = 0;
			Branch *cur = b;
			Edit e = cur->edit_;
			while(e.initialized()) {
				if     (e.pos < rdepth5_) hiHalfMms++;
				else if(e.pos < rdepth3_) loHalfMms++;
				else assert(false);
				cur = cur->parent_;
				assert(cur != NULL);
				e = cur->edit_;
			}
			assert(cur->parent_ == NULL);
			assert_leq(loHalfMms + hiHalfMms, lim3);
			bool invalidHalfAndHalf = (loHalfMms == 0 || hiHalfMms == 0);
			return (nedits >= (uint32_t)halfAndHalf_ && !invalidHalfAndHalf);
		}
#ifndef NDEBUG
		if(refdep < rdepth5_-1) {
			assert_leq(nedits, lim5);
		}
		else if(refdep >= rdepth5_ && refdep < rdepth3_-1) {
			assert_gt(nedits, 0);
			assert_leq(nedits, lim3);
		}
#endif
		return true;
	}

	/**
	 * Return true iff the state of the backtracker as encoded by
	 * stackDepth, d and iham is compatible with the current half-and-
	 * half alignment mode.
	 */
	bool hhCheckTop(Branch* b, uint32_t refdep, uint32_t iham) {
		// Crossing from the hi-half into the lo-half
		assert(!seedRange_.valid());
		ASSERT_ONLY(uint32_t lim3 = (rdeps_[3] == rdeps_[2])? 2 : 3);
		ASSERT_ONLY(uint32_t lim5 = (rdeps_[1] == rdeps_[0])? 2 : 1);
		const uint32_t nedits = b->numEdits_;
		if(b->edit_.isReadGap()) {
			// tipDepth() wasn't advanced when this topmost insertion
			// was added, so we have to add 1 to refdep to avoid
			// thinking that we just crossed the border, when in fact
			// we added an insert just after crossing the border.
			refdep++;
		}
		if(refdep == rdepth5_) {
			if(nedits == 0) {
				return false;
			}
			assert_leq(nedits, lim5);
		} else if(refdep == rdepth3_) {
			assert_leq(nedits, lim3);
			if(nedits < (uint32_t)halfAndHalf_) {
				return false;
			}
		}
#ifndef NDEBUG
		else {
			// We didn't just cross a boundary, so do an in-between check
			if(refdep >= rdepth5_) {
				assert_geq(nedits, 1);
			} else if(refdep >= rdepth3_) {
				assert_geq(nedits, lim3);
			}
		}
#endif
		return true;
	}

	/**
	 * Tally how many Ns occur in the seed region and in the ftab-
	 * jumpable region of the read.  Check whether the mismatches
	 * induced by the Ns already violates the current policy.  Return
	 * false if the policy is already violated, true otherwise.
	 */
	bool tallyNs(int& nsInSeed, int& nsInFtab) {
		const Ebwt& ebwt = *ebwt_;
		int ftabChars = ebwt._eh._ftabChars;
		// Count Ns in the seed region of the read and short-circuit if
		// the configuration of Ns guarantees that there will be no
		// valid alignments given the backtracking constraints.
		for(size_t i = 0; i < rdeps_[3]; i++) {
			if((int)(*ref_)[rlen_-i-1] == 4) {
				nsInSeed++;
				if(nsInSeed == 1) {
					if(i < rdeps_[0]) {
						return false; // Exceeded edit budget on Ns alone
					}
				} else if(nsInSeed == 2) {
					if(i < rdeps_[1]) {
						return false; // Exceeded edit budget on Ns alone
					}
				} else if(nsInSeed == 3) {
					if(i < rdeps_[2]) {
						return false; // Exceeded edit budget on Ns alone
					}
				} else {
					assert_gt(nsInSeed, 3);
					return false;     // Exceeded edit budget on Ns alone
				}
			}
		}
		// Calculate the number of Ns there are in the region that
		// would get jumped over if the ftab were used.
		for(size_t i = 0; i < (size_t)ftabChars && i < rlen_; i++) {
			if((int)(*ref_)[rlen_-i-1] == 4) nsInFtab++;
		}
		assert(!seedRange_.valid() || nsInSeed == 0);
		return true;
	}

	/**
	 * Calculate the offset into the ftab for the rightmost 'ftabChars'
	 * characters of the current query. Rightmost char gets least
	 * significant bit-pair.
	 */
	uint32_t calcFtabOff() {
		const Ebwt& ebwt = *ebwt_;
		int ftabChars = ebwt._eh._ftabChars;
		uint32_t ftabOff = (*ref_)[rlen_ - ftabChars];
		assert_lt(ftabOff, 4);
		assert_lt(ftabOff, ebwt._eh._ftabLen - 1);
		for(int i = ftabChars - 1; i > 0; i--) {
			ftabOff <<= 2;
			assert_lt((uint32_t)(*ref_)[rlen_ - i], 4);
			ftabOff |= (uint32_t)(*ref_)[rlen_ - i];
			assert_lt(ftabOff, ebwt._eh._ftabLen - 1);
		}
		assert_lt(ftabOff, ebwt._eh._ftabLen - 1);
		return ftabOff;
	}

	BTDnaString *qry_;     // sequence
	size_t       qlen_;    // sequence length (might be < qry_->length())
	size_t       rlen_;    // sequence length (might be < qry_->length())
	BTDnaString  qryBuf_;  // sequence with partial-alignment edits applied
	BTString    *qual_;    // quality
	BTString    *name_;    // name
	TReadId      patid_;   // pattern id
	BTDnaString *altQry_;  // alternate sequences (up to 3)
	BTString    *altQual_; // alternate qualities (up to 3)
	int          alts_;    // max # alternatives
	BTDnaString *ref_;     // reference characters already aligned to
	BTDnaString  refBuf_;  // reference characters already aligned to
	SStringExpandable<bool> qryDelDeps_;
	CoordMap     cmap_;
	const Ebwt*  ebwt_;   // Ebwt to search in
	bool         fw_;
	uint16_t     qdeps_[4]; // unrevisitable chunks w/r/t read
	uint16_t     rdeps_[4]; // unrevisitable chunks w/r/t reference
	/// Reject alignments where sum of qualities at mismatched
	/// positions is greater than qualLim_
	uint32_t            qualLim_;
	/// Report exact alignments iff this is true
	bool                reportExacts_;
	/// Whether to use the _os array together with a naive matching
	/// algorithm to double-check reported alignments (or the lack
	/// thereof)
	int                 halfAndHalf_;
	/// Whether we're generating partial alignments for a longer
	/// alignment in the opposite index.
	bool                partial_;
	/// Depth of 5'-seed-half border
	uint32_t            qdepth5_; // w/r/t read
	uint32_t            rdepth5_; // w/r/t reference
	/// Depth of 3'-seed-half border
	uint32_t            qdepth3_; // w/r/t read
	uint32_t            rdepth3_; // w/r/t reference
	/// Source of pseudo-random numbers
	RandomSource        rand_;
	// Current range to expose to consumers
	Range               curRange_;
	// Range for the partial alignment we're extending (NULL if we
	// aren't extending a partial)
	Range               seedRange_;
	// Starts as false; set to true as soon as we know we want to skip
	// all further processing of this read
	bool                skippingThisRead_;
	// Object encapsulating metrics
	AlignerMetrics*     metrics_;
	// Temporary expandable list for edits
	EList<Edit> edits_;
	// Matrix of sets of BW ranges visited
	Visit visit_;
	// lo/hiGapDepth_ are relatively tight bounds on the range of
	// BWT-SW *rows* where gaps (i.e. horizontal and vertical
	// movements) are permitted.
	int loGapDepth_;
	int hiGapDepth_;
	// Correct but not-necessarily tight upper bound on possible # gaps
	int maxGaps_;
	// Correct but not-necessarily tight upper bound on possible # edits
	int maxEdits_;
#ifndef NDEBUG
	std::set<int64_t>   allTops_;
#endif
};

/**
 * Concrete factory for EbwtRangeSource objects.
 */
class EbwtRangeSourceFactory {

public:
	EbwtRangeSourceFactory(
			const Ebwt* ebwt,
			bool         fw,
			uint32_t     qualThresh,
			bool         reportExacts,
			bool         halfAndHalf,
			bool         seeded,
			AlignerMetrics *metrics = NULL) :
			ebwt_(ebwt),
			fw_(fw),
			qualThresh_(qualThresh),
			reportExacts_(reportExacts),
			halfAndHalf_(halfAndHalf),
			seeded_(seeded),
			metrics_(metrics) { }

	/**
	 * Return new EbwtRangeSource with predefined params.s
	 */
	EbwtRangeSource *create() {
		return new EbwtRangeSource(ebwt_, fw_, qualThresh_,
		                           reportExacts_,
		                           halfAndHalf_, seeded_,
		                           metrics_);
	}

protected:
	const Ebwt*  ebwt_;
	bool         fw_;
	uint32_t     qualThresh_;
	bool         reportExacts_;
	bool         halfAndHalf_;
	bool         seeded_;
	AlignerMetrics *metrics_;
};

/**
 * What boundary within the alignment to "pin" a particular
 * backtracking constraint to.
 */
enum SearchConstraintExtent {
	PIN_TO_BEGINNING = 1, // depth 0; i.e., constraint is inactive
	PIN_TO_LEN,           // constraint applies to while alignment
	PIN_TO_HI_HALF_EDGE,  // constraint applies to hi-half of seed region
	PIN_TO_SEED_EDGE      // constraint applies to entire seed region
};

/**
 * Concrete RangeSourceDriver that deals properly with
 * GreedyDFSRangeSource by calling setOffs() with the appropriate
 * parameters when initializing it;
 */
class EbwtRangeSourceDriver :
	public SingleRangeSourceDriver<EbwtRangeSource>
{
public:
	EbwtRangeSourceDriver(
			EbwtSearchParams& params,
			EbwtRangeSource* rs,
			bool fw,
			bool seed,
			HitSink& sink,
			HitSinkPerThread* sinkPt,
			uint32_t seedLen,
			bool nudgeLeft,
			SearchConstraintExtent rev0Off,
			SearchConstraintExtent rev1Off,
			SearchConstraintExtent rev2Off,
			SearchConstraintExtent rev3Off,
			EList<SString<char> >& os,
			bool mate1,
			ChunkPool* pool,
			int *btCnt) :
			SingleRangeSourceDriver<EbwtRangeSource>(
				params, rs, fw, sink, sinkPt, os, mate1, 0, pool, btCnt),
			seed_(seed),
			rs_(rs), seedLen_(seedLen),
			nudgeLeft_(nudgeLeft),
			rev0Off_(rev0Off), rev1Off_(rev1Off),
			rev2Off_(rev2Off), rev3Off_(rev3Off)
	{
		assert(!seed_ || seedLen > 0);
	}

	virtual ~EbwtRangeSourceDriver() { }

	bool seed() const { return seed_; }

	bool ebwtFw() const { return rs_->curEbwt()->fw(); }

	/**
	 * Called every time setQuery() is called in the parent class,
	 * after setQuery() has been called on the RangeSource but before
	 * initConts() has been called.
	 */
	virtual void initRangeSource(const BTString& qual,
	                             int alts,
	                             const BTString* altQuals)
	{
		// If seedLen_ is huge, then it will always cover the whole
		// alignment
		assert_eq(len_, qual.length());
		uint32_t s = (seedLen_ > 0 ? min<uint32_t>((uint32_t)seedLen_, (uint32_t)len_) : (uint32_t)len_);
		uint32_t sLeft  = s >> 1;
		uint32_t sRight = s >> 1;
		// If seed has odd length, then nudge appropriate half up by 1
		if((s & 1) != 0) { if(nudgeLeft_) sLeft++; else sRight++; }
		uint32_t rev0Off = cextToDepth(rev0Off_, sRight, s, (uint32_t)len_);
		uint32_t rev1Off = cextToDepth(rev1Off_, sRight, s, (uint32_t)len_);
		uint32_t rev2Off = cextToDepth(rev2Off_, sRight, s, (uint32_t)len_);
		uint32_t rev3Off = cextToDepth(rev3Off_, sRight, s, (uint32_t)len_);
		// Truncate the pattern if necessary
		uint32_t qlen = (uint32_t)qual.length();
		if(seed_) {
			if(len_ > s) {
				rs_->setQlen(s);
				qlen = s;
			}
			assert(!rs_->reportExacts());
		}
		// If there are any Ns in the unrevisitable region, then this
		// driver is guaranteed to yield no fruit.
		uint16_t minCost = 0;
		if(rs_->reportExacts()) {
			// Keep minCost at 0
		} else if (!rs_->halfAndHalf() && rev0Off < s) {
			// Exacts not allowed, so there must be at least 1 mismatch
			// outside of the unrevisitable area
			minCost = 1 << 14;
			uint8_t lowQual = 0xff;
			for(uint32_t d = rev0Off; d < s; d++) {
				uint8_t lowAtPos;
				lowAtPos = loPenaltyAt(qlen-d-1, alts, qual, altQuals);
				if(lowAtPos < lowQual) lowQual = lowAtPos;
			}
			assert_lt(lowQual, 0xff);
			if(gGaps) {
				lowQual = min<uint32_t>(lowQual, gInsOpen);
				lowQual = min<uint32_t>(lowQual, gDelOpen);
			}
			minCost += lowQual;
		} else if(rs_->halfAndHalf() && sRight > 0 && sRight < (s-1)) {
			// Half-and-half constraints are active, so there must be
			// at least 1 mismatch in both halves of the seed
			assert(rs_->halfAndHalf());
			minCost = (seed_ ? 3 : 2) << 14;
			assert(rs_->halfAndHalf() == 2 || rs_->halfAndHalf() == 3);
			uint8_t lowQual1 = 0xff;
			for(uint32_t d = 0; d < sRight; d++) {
				uint8_t lowAtPos;
				lowAtPos = loPenaltyAt(qlen-d-1, alts, qual, altQuals);
				if(lowAtPos < lowQual1) lowQual1 = lowAtPos;
			}
			assert_lt(lowQual1, 0xff);
			if(gGaps) {
				lowQual1 = min<uint32_t>(lowQual1, gInsOpen);
				lowQual1 = min<uint32_t>(lowQual1, gDelOpen);
			}
			minCost += lowQual1;
			uint8_t lowQual2_1 = 0xff;
			uint8_t lowQual2_2 = 0xff;
			for(uint32_t d = sRight; d < s; d++) {
				uint8_t lowAtPos;
				lowAtPos = loPenaltyAt(qlen-d-1, alts, qual, altQuals);
				if(lowAtPos < lowQual2_1) {
					if(lowQual2_1 != 0xff) {
						lowQual2_2 = lowQual2_1;
					}
					lowQual2_1 = lowAtPos;
				} else if(lowAtPos < lowQual2_2) {
					lowQual2_2 = lowAtPos;
				}
			}
			assert_lt(lowQual2_1, 0xff);
			if(gGaps) {
				lowQual2_1 = min<uint32_t>(lowQual2_1, gInsOpen);
				lowQual2_1 = min<uint32_t>(lowQual2_1, gDelOpen);
				lowQual2_1 = min<uint32_t>(lowQual2_1, gDelExtend);
			}
			minCost += lowQual2_1;
			if(rs_->halfAndHalf() > 2 && lowQual2_2 != 0xff) {
				if(gGaps) {
					lowQual2_2 = min<uint32_t>(lowQual2_2, gInsOpen);
					lowQual2_2 = min<uint32_t>(lowQual2_2, gDelOpen);
					lowQual2_2 = min<uint32_t>(lowQual2_2, gInsExtend);
					lowQual2_2 = min<uint32_t>(lowQual2_2, gDelExtend);
				}
				minCost += lowQual2_2;
			}
		}
		if(gVerbose) cout << "initRangeSource minCost: " << minCost << endl;
		this->minCostAdjustment_ = minCost;
		rs_->setOffs(
			sRight,   // depth of far edge of hi-half (only matters where half-and-half is possible)
			s,        // depth of far edge of lo-half (only matters where half-and-half is possible)
			rev0Off,  // depth above which we cannot backtrack
			rev1Off,  // depth above which we may backtrack just once
			rev2Off,  // depth above which we may backtrack just twice
			rev3Off); // depth above which we may backtrack just three times
	}

protected:

	/**
	 * Convert a search constraint extent to an actual depth into the
	 * read.
	 */
	inline uint32_t cextToDepth(SearchConstraintExtent cext,
	                            uint32_t sRight,
	                            uint32_t s,
	                            uint32_t len)
	{
		if(cext == PIN_TO_SEED_EDGE)    return s;
		if(cext == PIN_TO_HI_HALF_EDGE) return sRight;
		if(cext == PIN_TO_BEGINNING)    return 0;
		if(cext == PIN_TO_LEN)          return len;
		cerr << "Bad SearchConstraintExtent: " << cext;
		throw 1;
	}

	bool seed_;
	EbwtRangeSource* rs_;
	uint32_t seedLen_;
	bool nudgeLeft_;
	SearchConstraintExtent rev0Off_;
	SearchConstraintExtent rev1Off_;
	SearchConstraintExtent rev2Off_;
	SearchConstraintExtent rev3Off_;
};

/**
 * Create appropriately-configured instances of EbwtRangeSourceDriver
 * on demand.
 */
class EbwtRangeSourceDriverFactory {
public:
	EbwtRangeSourceDriverFactory(
			EbwtSearchParams& params,
			EbwtRangeSourceFactory* rs,
			bool fw,
			bool seed,
			HitSink& sink,
			HitSinkPerThread* sinkPt,
			uint32_t seedLen,
			bool nudgeLeft,
			SearchConstraintExtent rev0Off,
			SearchConstraintExtent rev1Off,
			SearchConstraintExtent rev2Off,
			SearchConstraintExtent rev3Off,
			EList<SString<char> >& os,
			bool mate1,
			ChunkPool* pool,
			int *btCnt = NULL) :
			params_(params),
			rs_(rs),
			fw_(fw),
			seed_(seed),
			sink_(sink),
			sinkPt_(sinkPt),
			seedLen_(seedLen),
			nudgeLeft_(nudgeLeft),
			rev0Off_(rev0Off),
			rev1Off_(rev1Off),
			rev2Off_(rev2Off),
			rev3Off_(rev3Off),
			os_(os),
			mate1_(mate1),
			pool_(pool),
			btCnt_(btCnt)
	{ }

	~EbwtRangeSourceDriverFactory() {
		delete rs_; rs_ = NULL;
	}

	/**
	 * Return a newly-allocated EbwtRangeSourceDriver with the given
	 * parameters.
	 */
	EbwtRangeSourceDriver *create() const {
		return new EbwtRangeSourceDriver(
				params_, rs_->create(), fw_, seed_,
				sink_, sinkPt_, seedLen_, nudgeLeft_,
				rev0Off_, rev1Off_, rev2Off_, rev3Off_, os_,
				mate1_, pool_, btCnt_);
	}

protected:
	EbwtSearchParams& params_;
	EbwtRangeSourceFactory* rs_;
	bool fw_;
	bool seed_;
	HitSink& sink_;
	HitSinkPerThread* sinkPt_;
	uint32_t seedLen_;
	bool nudgeLeft_;
	SearchConstraintExtent rev0Off_;
	SearchConstraintExtent rev1Off_;
	SearchConstraintExtent rev2Off_;
	SearchConstraintExtent rev3Off_;
	EList<SString<char> >& os_;
	bool mate1_;
	ChunkPool* pool_;
	int *btCnt_;
};

/**
 * A RangeSourceDriver that manages two child EbwtRangeSourceDrivers,
 * one for searching for seed strings with mismatches in the hi-half,
 * and one for extending those seed strings toward the 3' end.
 */
class EbwtSeededRangeSourceDriver : public RangeSourceDriver<EbwtRangeSource> {
	typedef RangeSourceDriver<EbwtRangeSourceDriver>* TRangeSrcDrPtr;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
public:
	EbwtSeededRangeSourceDriver(
			EbwtRangeSourceDriverFactory* rsFact,
			EbwtRangeSourceDriver* rsSeed,
			bool fw,
			uint32_t seedLen,
			bool mate1) :
			RangeSourceDriver<EbwtRangeSource>(true, 0),
			rsFact_(rsFact), rsFull_(NULL, true),
			rsSeed_(rsSeed), patsrc_(NULL), seedLen_(seedLen), fw_(fw),
			mate1_(mate1), seedRange_(0)
	{
		assert(rsSeed_->seed());
	}

	virtual ~EbwtSeededRangeSourceDriver() {
		delete rsFact_; rsFact_ = NULL;
		delete rsSeed_; rsSeed_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQueryImpl(PatternSourcePerThread* patsrc, Range *partial) {
		this->done = false;
		rsSeed_->setQuery(patsrc, partial);
		this->minCostAdjustment_ = max(rsSeed_->minCostAdjustment_, rsSeed_->minCost);
		this->minCost = this->minCostAdjustment_;
		rsFull_.clearSources();
		rsFull_.setQuery(patsrc, partial);
		rsFull_.minCost = this->minCost;
		assert_gt(rsFull_.minCost, 0);
		patsrc_ = patsrc;
		// The minCostAdjustment comes from the seed range source
		// driver, based on Ns and quals in the hi-half
		this->foundRange = false;
		assert_eq(this->minCost, min<uint16_t>(rsSeed_->minCost, rsFull_.minCost));
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual void advance(int until) {
		assert(!this->foundRange);
		until = max<int>(until, ADV_COST_CHANGES);
		ASSERT_ONLY(uint16_t preCost = this->minCost);
		advanceImpl(until);
		if(this->foundRange) {
			assert_eq(range().cost(), preCost);
		}
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual void advanceImpl(int until) {
		assert(!this->done);
		assert(!this->foundRange);
		assert_gt(rsFull_.minCost, 0);
		// Advance the seed range source
		if(rsSeed_->done && rsFull_.done &&
		   !rsSeed_->foundRange && !rsFull_.foundRange)
		{
			this->done = true;
			return;
		}
		if(rsSeed_->done && !rsSeed_->foundRange) {
			rsSeed_->minCost = 0xffff;
			if(rsFull_.minCost > this->minCost) {
				this->minCost = rsFull_.minCost;
				// Cost changed, so return
				return;
			}
		}
		if(rsFull_.done && !rsFull_.foundRange) {
			rsFull_.minCost = 0xffff;
			if(rsSeed_->minCost > this->minCost) {
				this->minCost = rsSeed_->minCost;
				// Cost changed, so return
				return;
			}
		}
		assert(rsSeed_->minCost != 0xffff || rsFull_.minCost != 0xffff);
		// Extend a partial alignment
		ASSERT_ONLY(uint16_t oldMinCost = this->minCost);
		assert_eq(this->minCost, min<uint16_t>(rsSeed_->minCost, rsFull_.minCost));
		bool doFull = rsFull_.minCost <= rsSeed_->minCost;
		if(!doFull) {
			// Advance the partial-alignment generator
			assert_eq(rsSeed_->minCost, this->minCost);
			if(!rsSeed_->foundRange) {
				rsSeed_->advance(until);
			}
			if(rsSeed_->foundRange) {
				assert_eq(this->minCost, rsSeed_->range().cost());
				assert_eq(oldMinCost, rsSeed_->range().cost());
				seedRange_ = &rsSeed_->range();
				rsSeed_->foundRange = false;
				assert_geq(seedRange_->cost(), this->minCostAdjustment_);
				this->minCostAdjustment_ = seedRange_->cost();
				assert_gt(seedRange_->edits().size(), 0);
				// Keep the range for the hi-half partial alignment so
				// that the driver can (a) modify the pattern string
				// and (b) modify results from the RangeSource to
				// include these edits.
				EbwtRangeSourceDriver *partial = rsFact_->create();
				partial->minCost = seedRange_->cost();
				rsFull_.minCost = seedRange_->cost();
				rsFull_.addSource(partial, seedRange_);
				if(rsFull_.foundRange) {
					this->foundRange = true;
					rsFull_.foundRange = false;
					assert(rsFull_.range().repOk());
					assert_eq(range().cost(), oldMinCost);
				}
			}
			if(rsSeed_->minCost > this->minCost) {
				this->minCost = rsSeed_->minCost;
				if(!rsFull_.done) {
					this->minCost = min(this->minCost, rsFull_.minCost);
					assert_eq(this->minCost, min<uint16_t>(rsSeed_->minCost, rsFull_.minCost));
				}
			}
		} else {
			// Extend a full alignment
			assert(!rsFull_.done);
			assert(!rsFull_.foundRange);
			uint16_t oldFullCost = rsFull_.minCost;
			if(!rsFull_.foundRange) {
				rsFull_.advance(until);
			}
			// Found a minimum-cost range
			if(rsFull_.foundRange) {
				this->foundRange = true;
				rsFull_.foundRange = false;
				assert(rsFull_.range().repOk());
				assert_eq(range().cost(), oldMinCost);
			}
			assert_geq(rsFull_.minCost, oldFullCost);
			// Did the min cost change?
			if(rsFull_.minCost > oldFullCost) {
				// If a range was found, hold on to it and save it for
				// later.  Update the minCost.
				assert(!rsSeed_->done || rsSeed_->minCost == 0xffff);
				this->minCost = min(rsFull_.minCost, rsSeed_->minCost);
			}
		}
	}

	/**
	 * Return the range found.
	 */
	virtual Range& range() {
		Range& r = rsFull_.range();
		r.setMate1(mate1_);
		return r;
	}

	/**
	 * Return whether we're generating ranges for the first or the
	 * second mate.
	 */
	virtual bool mate1() const {
		return mate1_;
	}

	/**
	 * Return true iff current pattern is forward-oriented.
	 */
	virtual bool fw() const {
		return fw_;
	}

protected:

	EbwtRangeSourceDriverFactory* rsFact_;
	TCostAwareRangeSrcDr rsFull_;
	EbwtRangeSourceDriver* rsSeed_;
	PatternSourcePerThread* patsrc_;
	uint32_t seedLen_;
	bool fw_;
	bool mate1_;
	bool generating_;
	Range *seedRange_;
};

#endif /*EBWT_SEARCH_BACKTRACK_H_*/
