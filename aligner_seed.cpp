/*
 *  aligner_seed.cpp
 */

#include "aligner_cache.h"
#include "aligner_seed.h"
#include "aligner_counters.h"
#include "search_globals.h"
#include "ebwt.h"
#include <algorithm>

using namespace std;

/**
 * Construct a constraint with no edits of any kind allowed.
 */
Constraint Constraint::exact() {
	Constraint c;
	c.edits = c.mms = c.ins = c.dels = c.penalty = 0;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint.
 */
Constraint Constraint::penaltyBased(int pen) {
	Constraint c;
	c.penalty = pen;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint related to the length of the read.
 */
Constraint Constraint::penaltyFuncBased(float penConst, float penLinear) {
	Constraint c;
	c.penConst = penConst;
	c.penLinear = penLinear;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint.
 */
Constraint Constraint::mmBased(int mms) {
	Constraint c;
	c.mms = mms;
	c.edits = c.dels = c.ins = 0;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint.
 */
Constraint Constraint::editBased(int edits) {
	Constraint c;
	c.edits = edits;
	c.dels = c.ins = c.mms = 0;
	return c;
}

//
// Some static methods for constructing some standard SeedPolicies
//

/**
 * Given a read, depth and orientation, extract a seed data structure
 * from the read and fill in the steps & zones arrays.  The Seed
 * contains the sequence and quality values.
 */
bool
Seed::instantiate(
	const Read& read,
	const BTDnaString& seq, // seed read sequence
	const BTString& qual,   // seed quality sequence
	const Scoring& pens,
	int depth,
	int seedoffidx,
	int seedtypeidx,
	bool fw,
	InstantiatedSeed& is) const
{
	assert(overall != NULL);
	int seedlen = len;
	if((int)read.length() < seedlen) {
		// Shrink seed length to fit read if necessary
		seedlen = (int)read.length();
	}
	assert_gt(seedlen, 0);
	is.steps.resize(seedlen);
	is.zones.resize(seedlen);
	// Fill in 'steps' and 'zones'
	//
	// The 'steps' list indicates which read character should be
	// incorporated at each step of the search process.  Often we will
	// simply proceed from one end to the other, in which case the
	// 'steps' list is ascending or descending.  In some cases (e.g.
	// the 2mm case), we might want to switch directions at least once
	// during the search, in which case 'steps' will jump in the
	// middle.  When an element of the 'steps' list is negative, this
	// indicates that the next
	//
	// The 'zones' list indicates which zone constraint is active at
	// each step.  Each element of the 'zones' list is a pair; the
	// first pair element indicates the applicable zone when
	// considering either mismatch or delete (ref gap) events, while
	// the second pair element indicates the applicable zone when
	// considering insertion (read gap) events.  When either pair
	// element is a negative number, that indicates that we are about
	// to leave the zone for good, at which point we may need to
	// evaluate whether we have reached the zone's budget.
	//
	switch(type) {
		case SEED_TYPE_EXACT: {
			for(int k = 0; k < seedlen; k++) {
				is.steps[k] = -(seedlen - k);
				// Zone 0 all the way
				is.zones[k].first = is.zones[k].second = 0;
			}
			break;
		}
		case SEED_TYPE_LEFT_TO_RIGHT: {
			for(int k = 0; k < seedlen; k++) {
				is.steps[k] = k+1;
				// Zone 0 from 0 up to ceil(len/2), then 1
				is.zones[k].first = is.zones[k].second = ((k < (seedlen+1)/2) ? 0 : 1);
			}
			// Zone 1 ends at the RHS
			is.zones[seedlen-1].first = is.zones[seedlen-1].second = -1;
			break;
		}
		case SEED_TYPE_RIGHT_TO_LEFT: {
			for(int k = 0; k < seedlen; k++) {
				is.steps[k] = -(seedlen - k);
				// Zone 0 from 0 up to floor(len/2), then 1
				is.zones[k].first  = ((k < seedlen/2) ? 0 : 1);
				// Inserts: Zone 0 from 0 up to ceil(len/2)-1, then 1
				is.zones[k].second = ((k < (seedlen+1)/2+1) ? 0 : 1);
			}
			is.zones[seedlen-1].first = is.zones[seedlen-1].second = -1;
			break;
		}
		case SEED_TYPE_INSIDE_OUT: {
			// Zone 0 from ceil(N/4) up to N-floor(N/4)
			int step = 0;
			for(int k = (seedlen+3)/4; k < seedlen - (seedlen/4); k++) {
				is.zones[step].first = is.zones[step].second = 0;
				is.steps[step++] = k+1;
			}
			// Zone 1 from N-floor(N/4) up
			for(int k = seedlen - (seedlen/4); k < seedlen; k++) {
				is.zones[step].first = is.zones[step].second = 1;
				is.steps[step++] = k+1;
			}
			// No Zone 1 if seedlen is short (like 2)
			//assert_eq(1, is.zones[step-1].first);
			is.zones[step-1].first = is.zones[step-1].second = -1;
			// Zone 2 from ((seedlen+3)/4)-1 down to 0
			for(int k = ((seedlen+3)/4)-1; k >= 0; k--) {
				is.zones[step].first = is.zones[step].second = 2;
				is.steps[step++] = -(k+1);
			}
			assert_eq(2, is.zones[step-1].first);
			is.zones[step-1].first = is.zones[step-1].second = -2;
			assert_eq(seedlen, step);
			break;
		}
		default:
			throw 1;
	}
	// Instantiate constraints
	for(int i = 0; i < 3; i++) {
		is.cons[i] = zones[i];
		is.cons[i].instantiate(read.length());
	}
	is.overall = *overall;
	is.overall.instantiate(read.length());
	// Take a sweep through the seed sequence.  Consider where the Ns
	// occur and how zones are laid out.  Calculate the maximum number
	// of positions we can jump over initially (e.g. with the ftab) and
	// perhaps set this function's return value to false, indicating
	// that the arrangements of Ns prevents the seed from aligning.
	bool streak = true;
	is.maxjump = 0;
	bool ret = true;
	bool ltr = (is.steps[0] > 0); // true -> left-to-right
	for(size_t i = 0; i < is.steps.size(); i++) {
		assert_neq(0, is.steps[i]);
		int off = is.steps[i];
		off = abs(off)-1;
		Constraint& cons = is.cons[abs(is.zones[i].first)];
		int c = seq[off];  assert_range(0, 4, c);
		int q = qual[off];
		if(ltr != (is.steps[i] > 0) || // changed direction
		   is.zones[i].first != 0 ||   // changed zone
		   is.zones[i].second != 0)    // changed zone
		{
			streak = false;
		}
		if(c == 4) {
			// Induced mismatch
			if(cons.canN(q, pens)) {
				cons.chargeN(q, pens);
			} else {
				// Seed disqualified due to arrangement of Ns
				return false;
			}
		}
		if(streak) is.maxjump++;
	}
	is.seedoff = depth;
	is.seedoffidx = seedoffidx;
	is.fw = fw;
	is.s = *this;
	return ret;
}

/**
 * Return a set consisting of 1 seed encapsulating an exact matching
 * strategy.
 */
void
Seed::zeroMmSeeds(int ln, EList<Seed>& pols, Constraint& oall) {
	oall.init();
	// Seed policy 1: left-to-right search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_EXACT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::exact();
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
}

/**
 * Return a set of 2 seeds encapsulating a half-and-half 1mm strategy.
 */
void
Seed::oneMmSeeds(int ln, EList<Seed>& pols, Constraint& oall) {
	oall.init();
	// Seed policy 1: left-to-right search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
	// Seed policy 2: right-to-left search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[1].mmsCeil = 0;
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
}

/**
 * Return a set of 3 seeds encapsulating search roots for:
 *
 * 1. Starting from the left-hand side and searching toward the
 *    right-hand side allowing 2 mismatches in the right half.
 * 2. Starting from the right-hand side and searching toward the
 *    left-hand side allowing 2 mismatches in the left half.
 * 3. Starting (effectively) from the center and searching out toward
 *    both the left and right-hand sides, allowing one mismatch on
 *    either side.
 *
 * This is not exhaustive.  There are 2 mismatch cases mised; if you
 * imagine the seed as divided into four successive quarters A, B, C
 * and D, the cases we miss are when mismatches occur in A and C or B
 * and D.
 */
void
Seed::twoMmSeeds(int ln, EList<Seed>& pols, Constraint& oall) {
	oall.init();
	// Seed policy 1: left-to-right search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(2);
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
	// Seed policy 2: right-to-left search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(2);
	pols.back().zones[1].mmsCeil = 1; // Must have used at least 1 mismatch
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
	// Seed policy 3: inside-out search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_INSIDE_OUT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[1].mmsCeil = 0; // Must have used at least 1 mismatch
	pols.back().zones[2] = Constraint::mmBased(1);
	pols.back().zones[2].mmsCeil = 0; // Must have used at least 1 mismatch
	pols.back().overall = &oall;
}

/**
 * Types of actions that can be taken by the SeedAligner.
 */
enum {
	SA_ACTION_TYPE_RESET = 1,
	SA_ACTION_TYPE_SEARCH_SEED, // 2
	SA_ACTION_TYPE_FTAB,        // 3
	SA_ACTION_TYPE_FCHR,        // 4
	SA_ACTION_TYPE_MATCH,       // 5
	SA_ACTION_TYPE_EDIT         // 6
};

// Action for when a new search is begun
struct SAActionReset : public SAAction {
	SAActionReset() : SAAction() {
		type = SA_ACTION_TYPE_RESET;
	}
};

// Action for when a new search seed is begun
struct SAActionSearchSeed : public SAAction {
	SAActionSearchSeed(int of, int sd) : SAAction() {
		type = SA_ACTION_TYPE_SEARCH_SEED; seedoff = of; seed = sd;
	}
};

// Action for when an ftab jump is made
struct SAActionFtab : public SAAction {
	SAActionFtab(int of, int sd, int ps, int ln, bool lr) : SAAction() {
		type = SA_ACTION_TYPE_FTAB; seedoff = of; seed = sd; pos = ps; len = ln; ltr = lr;
	}
};

// Action for when an fchr jump is made
struct SAActionFchr : public SAAction {
	SAActionFchr(int of, int sd, int ps, bool lr) : SAAction() {
		type = SA_ACTION_TYPE_FCHR; seedoff = of; seed = sd; pos = ps; len = 1; ltr = lr;
	}
};

// Action for when we advance by matching the read character with a
// reference character
struct SAActionMatch : public SAAction {
	SAActionMatch(int of, int sd, int ps, bool lr, int dp) : SAAction() {
		type = SA_ACTION_TYPE_MATCH; seedoff = of; seed = sd; pos = ps; depth = dp; len = 1; ltr = lr;
	}
};

// Action for when we advance by editing the read character into a
// reference character
struct SAActionEdit : public SAAction {
	SAActionEdit(int of, int sd, int ps, bool lr, int dp, Edit e) : SAAction() {
		type = SA_ACTION_TYPE_EDIT; seedoff = of; seed = sd; pos = ps; depth = dp; len = 1; ltr = lr; edit = e;
	}
};

#define MIN(x, y) ((x < y) ? x : y)

#ifdef ALIGNER_SEED_TRACE
// Counters
static SACounters counters;
// History of actions performed by the seed aligner each time
// searchAllSeeds is called
static EList<SAAction> actions;
#define SA_RESET() { \
	if(counterSinks_ != NULL && !counterSinks_->empty()) { \
		for(size_t i = 0; i < counterSinks_->size(); i++) { \
			(*counterSinks_)[i]->reportCounters(*read_, counters); \
	}	} \
	if(actionSinks_ != NULL && !actionSinks_->empty()) { \
		for(size_t i = 0; i < actionSinks_->size(); i++) { \
			(*actionSinks_)[i]->reportActions(*read_, actions); \
	}	} \
	counters.reset(); \
	actions.clear(); \
	actions.push_back(SAActionReset()); \
}
#define SA_SEARCH_SEED(x, y) counters.seed++; actions.push_back(SAActionSearchSeed(x, y))
#define SA_FTAB(x, y, z) counters.ftab++; actions.push_back(SAActionFtab(s.seedoff, s.s.type, x, y, z))
#define SA_FCHR(x, y) counters.fchr++; actions.push_back(SAActionFchr(s.seedoff, s.s.type, x, y))
#define SA_MATCH(x, y, z) counters.match++; counters.matchd[MIN(z, 3)]++; actions.push_back(SAActionMatch(s.seedoff, s.s.type, x, y, z))
#define SA_EDIT(x, y, z, w) counters.match++; counters.editd[MIN(z, 3)]++; actions.push_back(SAActionEdit(s.seedoff, s.s.type, x, y, z, w))
#define SA_HIT() { \
	counters.hits++; \
	if(sinks_ != NULL && !sinks_->empty()) { \
		for(size_t i = 0; i < sinks_->size(); i++) { \
			(*sinks_)[i]->reportSeedHit(*read_, *seq_); } } }
#else
#define SA_RESET()
#define SA_SEARCH_SEED(x, y)
#define SA_FTAB(x, y, z)
#define SA_FCHR(x, y)
#define SA_MATCH(x, y, z)
#define SA_EDIT(x, y, z, w)
#define SA_HIT()
#endif /*ndef ALIGNER_SEED_TRACE*/

/**
 * Given a read and a few coordinates that describe a substring of the
 * read (or its reverse complement), fill in 'seq' and 'qual' objects
 * with the seed sequence and qualities.
 */
void
SeedAligner::instantiateSeq(
	const Read& read, // input read
	BTDnaString& seq, // output sequence
	BTString& qual,   // output qualities
	int len,          // seed length
	int depth,        // seed's 0-based offset from 5' end
	bool fw) const    // seed's orientation
{
	// Fill in 'seq' and 'qual'
	int seedlen = len;
	if((int)read.length() < seedlen) seedlen = (int)read.length();
	seq.resize(len);
	qual.resize(len);
	// If fw is false, we take characters starting at the 3' end of the
	// reverse complement of the read.
	for(int i = 0; i < len; i++) {
		seq.set(read.patFw.windowGetDna(i, fw, read.color, depth, len), i);
		qual.set(read.qual.windowGet(i, fw, depth, len), i);
	}
}

/**
 * We assume that all seeds are the same length.
 *
 * For each seed, tnstantiate the seed, retracting if necessary.
 */
pair<int, int> SeedAligner::instantiateSeeds(
	const EList<Seed>& seeds,  // search seeds
	int per,                   // interval between seeds
	const Read& read,          // read to align
	const Scoring& pens,       // scoring scheme
	float nCeilConst,          // ceil on # Ns w/r/t read len, const coeff
	float nCeilLinear,         // ceil on # Ns w/r/t read len, linear coeff
	bool nofw,                 // don't align forward read
	bool norc,                 // don't align revcomp read
	AlignmentCacheIface& cache,// holds some seed hits from previous reads
	SeedResults& sr,           // holds all the seed hits
	SeedSearchMetrics& met)    // metrics
{
	assert(!seeds.empty());
	assert_gt(read.length(), 0);
	// Check whether read has too many Ns
	offIdx2off_.clear();
	int len = seeds[0].len;
#ifndef NDEBUG
	for(size_t i = 1; i < seeds.size(); i++) {
		assert_eq(len, seeds[i].len);
	}
#endif
	// Calc # seeds within read interval
	int nseeds = 1;
	if((int)read.length() > len) {
		nseeds += ((int)read.length() - len) / per;
	}
	for(int i = 0; i < nseeds; i++) {
		offIdx2off_.push_back(per * i);
	}
	pair<int, int> ret;
	ret.first = 0;  // # seeds that require alignment
	ret.second = 0; // # seeds that hit in cache with non-empty results
	sr.reset(read, offIdx2off_, nseeds);
	assert(sr.repOk(&cache.current(), true)); // require that SeedResult be initialized
	// For each seed position
	for(int fwi = 0; fwi < 2; fwi++) {
		bool fw = (fwi == 0);
		if((fw && nofw) || (!fw && norc)) {
			// Skip this orientation b/c user specified --nofw or --norc
			continue;
		}
		// For each seed position
		for(int i = 0; i < nseeds; i++) {
			int depth = i * per;
			int seedlen = seeds[0].len;
			// Extract the seed sequence at this offset
			// If fw == true, we extract the characters from i*per to
			// i*(per-1) (exclusive).  If fw == false, 
			instantiateSeq(
				read,
				sr.seqs(fw)[i],
				sr.quals(fw)[i],
				std::min<int>((int)seedlen, (int)read.length()),
				depth,
				fw);
			QKey qk(sr.seqs(fw)[i]);
			QVal* qv;
			if(qk.cacheable() && (qv = cache.queryCopy(qk)) != NULL) {
				assert(qv->repOk(cache.current()));
				sr.add(*qv, cache.current(), i, fw, len);
				met.interhit++;
				if(!qv->empty()) ret.second++;
			} else {
				// For each search strategy
				EList<InstantiatedSeed>& iss = sr.instantiatedSeeds(fw, i);
				for(int j = 0; j < (int)seeds.size(); j++) {
					iss.expand();
					assert_eq(seedlen, seeds[j].len);
					InstantiatedSeed* is = &iss.back();
 					if(seeds[j].instantiate(
						read,
						sr.seqs(fw)[i],
						sr.quals(fw)[i],
						pens,
						depth,
						i,
						j,
						fw,
						*is))
					{
						// Can we fill this seed hit in from the cache?
						ret.first++;
					} else {
						// Seed may fail to instantiate if there are Ns
						// that prevent it from matching
						met.filteredseed++;
						iss.pop_back();
					}
				}
			}
		}
	}
	return ret;
}

/**
 * We assume that all seeds are the same length.
 *
 * For each seed:
 *
 * 1. Instantiate all seeds, retracting them if necessary.
 * 2. Calculate zone boundaries for each seed
 */
void SeedAligner::searchAllSeeds(
	const EList<Seed>& seeds,    // search seeds
	const Ebwt* ebwtFw,          // BWT index
	const Ebwt* ebwtBw,          // BWT' index
	const Read& read,            // read to align
	const Scoring& pens,         // scoring scheme
	AlignmentCacheIface& cache,  // local cache for seed alignments
	SeedResults& sr,             // holds all the seed hits
	SeedSearchMetrics& met,      // metrics
	EList<ReadCounterSink*>* readCounterSink,// if non-NULL, list of sinks to send per-read counter updates to
	EList<SeedHitSink*>*     sinks,          // if non-NULL, list of sinks to send hits to
	EList<SeedCounterSink*>* counterSinks,   // if non-NULL, list of sinks to send SACounters to
	EList<SeedActionSink*>*  actionSinks)    // if non-NULL, list of sinks to send SAActions to
{
	assert(!seeds.empty());
	assert(ebwtFw != NULL);
	assert(ebwtFw->isInMemory());
	assert(sr.repOk(&cache.current()));
	int len = seeds[0].len;
	ebwtFw_ = ebwtFw;
	ebwtBw_ = ebwtBw;
	sc_ = &pens;
	sinks_ = sinks;
	counterSinks_ = counterSinks;
	actionSinks_ = actionSinks;
	read_ = &read;
	ca_ = &cache;
	SA_RESET();
	bwops_ = bwedits_ = 0;
	uint64_t possearches = 0, seedsearches = 0, intrahits = 0, interhits = 0, ooms = 0;
	// For each instantiated seed
	for(int i = 0; i < (int)sr.numOffs(); i++) {
		for(int fwi = 0; fwi < 2; fwi++) {
			bool fw = (fwi == 0);
			assert(sr.repOk(&cache.current()));
			EList<InstantiatedSeed>& iss = sr.instantiatedSeeds(fw, i);
			if(iss.empty()) {
				// Cache hit in an across-read cache
				continue;
			}
			QVal qv;
			seq_ = &sr.seqs(fw)[i];
			qual_ = &sr.quals(fw)[i];
			// Tell the cache that we've started aligning, so the cache can
			// expect a series of on-the-fly updates
			int ret = cache.beginAlign(*seq_, *qual_, qv);
			ASSERT_ONLY(hits_.clear());
			if(ret == -1) {
				// Out of memory when we tried to add key to map
				ooms++;
				continue;
			}
			bool abort = false;
			if(ret == 0) {
				assert(cache.aligning());
				possearches++;
				for(size_t j = 0; j < iss.size(); j++) {
					// Set seq_ and qual_ appropriately, using the seed sequences
					// and qualities already installed in SeedResults
					assert_eq(fw, iss[j].fw);
					assert_eq(i, (int)iss[j].seedoffidx);
					s_ = &iss[j];
					SA_SEARCH_SEED(i, s_->seedtypeidx);
					// Do the search!
					if(!searchSeedBi()) {
						// Memory exhausted during search
						ooms++;
						abort = true;
						break;
					}
					seedsearches++;
					assert(cache.aligning());
				}
				if(!abort) {
					qv = cache.finishAlign();
				}
			} else {
				assert_eq(1, ret);
				assert(qv.valid());
				intrahits++;
			}
			assert(abort || !cache.aligning());
			if(qv.valid()) {
				sr.add(
					qv,    // range of ranges in cache
					cache.current(), // cache
					i,     // seed index (from 5' end)
					fw,    // whether seed is from forward read
					len);  // length of each seed
			}
		}
	}
	met.seedsearch += seedsearches;
	met.possearch += possearches;
	met.intrahit += intrahits;
	met.interhit += interhits;
	met.ooms += ooms;
	met.bwops += bwops_;
	met.bweds += bwedits_;
}

/**
 * Wrapper for initial invcation of searchSeed.
 */
bool
SeedAligner::searchSeedBi() {
	return searchSeedBi(
		0, 0,
		0, 0, 0, 0,
		SideLocus(), SideLocus(),
		s_->cons[0], s_->cons[1], s_->cons[2], s_->overall,
		NULL);
}

/**
 * Get tloc, bloc ready for the next step.  If the new range is under
 * the ceiling.
 */
inline void
SeedAligner::nextLocsBi(
	SideLocus& tloc,            // top locus
	SideLocus& bloc,            // bot locus
	uint32_t topf,              // top in BWT
	uint32_t botf,              // bot in BWT
	uint32_t topb,              // top in BWT'
	uint32_t botb,              // bot in BWT'
	int step                    // step to get ready for
#if 0
	, const SABWOffTrack* prevOt, // previous tracker
	SABWOffTrack& ot            // current tracker
#endif
	)
{
	assert_gt(botf, 0);
	assert(ebwtBw_ == NULL || botb > 0);
	assert_gt(step, 0); // next step can't be first one
	assert(ebwtBw_ == NULL || botf-topf == botb-topb);
	if(step == (int)s_->steps.size()) return; // no more steps!
	// Which direction are we going in next?
	if(s_->steps[step] > 0) {
		// Left to right; use BWT'
		if(botb - topb == 1) {
			// Already down to 1 row; just init top locus
			tloc.initFromRow(topb, ebwtBw_->eh(), ebwtBw_->ebwt());
			bloc.invalidate();
		} else {
			SideLocus::initFromTopBot(
				topb, botb, ebwtBw_->eh(), ebwtBw_->ebwt(), tloc, bloc);
			assert(bloc.valid());
		}
	} else {
		// Right to left; use BWT
		if(botf - topf == 1) {
			// Already down to 1 row; just init top locus
			tloc.initFromRow(topf, ebwtFw_->eh(), ebwtFw_->ebwt());
			bloc.invalidate();
		} else {
			SideLocus::initFromTopBot(
				topf, botf, ebwtFw_->eh(), ebwtFw_->ebwt(), tloc, bloc);
			assert(bloc.valid());
		}
	}
	// Check if we should update the tracker with this refinement
#if 0
	if(botf-topf <= BW_OFF_TRACK_CEIL) {
		if(ot.size() == 0 && prevOt != NULL && prevOt->size() > 0) {
			// Inherit state from the predecessor
			ot = *prevOt;
		}
		bool ltr = s_->steps[step-1] > 0;
		int adj = abs(s_->steps[step-1])-1;
		const Ebwt* ebwt = ltr ? ebwtBw_ : ebwtFw_;
		ot.update(
			ltr ? topb : topf,    // top
			ltr ? botb : botf,    // bot
			adj,                  // adj (to be subtracted from offset)
			ebwt->offs(),         // offs array
			ebwt->eh().offRate(), // offrate (sample = every 1 << offrate elts)
			NULL                  // dead
		);
		assert_gt(ot.size(), 0);
	}
#endif
	assert(botf - topf == 1 ||  bloc.valid());
	assert(botf - topf > 1  || !bloc.valid());
}

/**
 * Report a seed hit found by searchSeedBi() by adding it to the
 * SeedResults data structure.  Return false if the hit could not be
 * reported because of, e.g., cache exhaustion.
 */
bool
SeedAligner::reportHit(
	uint32_t topf,                     // top in BWT
	uint32_t botf,                     // bot in BWT
	uint32_t topb,                     // top in BWT'
	uint32_t botb,                     // bot in BWT'
	uint16_t len,                      // length of hit
	DoublyLinkedList<Edit> *prevEdit)  // previous edit
{
	// Add information about the seed hit to AlignmentCache.  This
	// information eventually makes its way back to the SeedResults
	// object when we call finishAlign(...).  finishAlign() will
	// populate an ACRangeRange appropriately.
	BTDnaString rf;
	edits_.clear();
	if(prevEdit != NULL) {
		prevEdit->toList(edits_);
		Edit::sort(edits_);
		assert(Edit::repOk(edits_, *seq_));
		Edit::toRef(*seq_, edits_, rf);
	} else {
		rf = *seq_;
	}
	// Sanity check: shouldn't add the same hit twice.  If this
	// happens, it may be because our zone Constraints are not set up
	// properly and erroneously return true from acceptable() when they
	// should return false in some cases.
	assert_eq(hits_.size(), ca_->curNumRanges());
	assert(hits_.insert(rf));
	if(!ca_->addOnTheFly(rf, topf, botf)) {
		return false;
	}
	assert_eq(hits_.size(), ca_->curNumRanges());
	SA_HIT();
#ifndef NDEBUG
	// Sanity check that the topf/botf and topb/botb ranges really
	// correspond to the reference sequence aligned to
	{
		BTDnaString rfr;
		uint32_t tpf, btf, tpb, btb;
		tpf = btf = tpb = btb = 0;
		assert(ebwtFw_->contains(rf, &tpf, &btf));
		if(ebwtBw_ != NULL) {
			rfr = rf;
			rfr.reverse();
			assert(ebwtBw_->contains(rfr, &tpb, &btb));
			assert_eq(tpf, topf);
			assert_eq(btf, botf);
			assert_eq(tpb, topb);
			assert_eq(btb, botb);
		}
	}
#endif
	return true;
}

/**
 * Given a seed, search.  Assumes zone 0 = no backtracking.
 *
 * Return a list of Seed hits.
 * 1. Edits
 * 2. Bidirectional BWT range(s) on either end
 */
bool
SeedAligner::searchSeedBi(
	int step,             // depth into steps_[] array
	int depth,            // recursion depth
	uint32_t topf,        // top in BWT
	uint32_t botf,        // bot in BWT
	uint32_t topb,        // top in BWT'
	uint32_t botb,        // bot in BWT'
	SideLocus tloc,       // locus for top (perhaps unititialized)
	SideLocus bloc,       // locus for bot (perhaps unititialized)
	Constraint c0,        // constraints to enforce in seed zone 0
	Constraint c1,        // constraints to enforce in seed zone 1
	Constraint c2,        // constraints to enforce in seed zone 2
	Constraint overall,   // overall constraints to enforce
	DoublyLinkedList<Edit> *prevEdit  // previous edit
#if 0
	, const SABWOffTrack* prevOt // prev off tracker (if tracking started)
#endif
	)
{
	assert(s_ != NULL);
	const InstantiatedSeed& s = *s_;
	assert_gt(s.steps.size(), 0);
	assert(ebwtBw_ == NULL || ebwtBw_->eh().ftabChars() == ebwtFw_->eh().ftabChars());
#ifndef NDEBUG
	for(int i = 0; i < 4; i++) {
		assert(ebwtBw_ == NULL || ebwtBw_->fchr()[i] == ebwtFw_->fchr()[i]);
	}
#endif
	if(step == (int)s.steps.size()) {
		// Finished aligning seed
		assert(c0.acceptable());
		assert(c1.acceptable());
		assert(c2.acceptable());
		if(!reportHit(topf, botf, topb, botb, seq_->length(), prevEdit)) {
			return false; // Memory exhausted
		}
		return true;
	}
#ifndef NDEBUG
	if(depth > 0) {
		assert(botf - topf == 1 ||  bloc.valid());
		assert(botf - topf > 1  || !bloc.valid());
	}
#endif
	int off;
	uint32_t tp[4], bp[4]; // dest BW ranges for "prime" index
	if(step == 0) {
		// Just starting
		assert(prevEdit == NULL);
		assert(!tloc.valid());
		assert(!bloc.valid());
		off = s.steps[0];
		bool ltr = off > 0;
		off = abs(off)-1;
		// Check whether/how far we can jump using ftab or fchr
		int ftabLen = ebwtFw_->eh().ftabChars();
		if(ftabLen > 1 && ftabLen <= s.maxjump) {
			if(!ltr) {
				assert_geq(off+1, ftabLen-1);
				off = off - ftabLen + 1;
			}
			ebwtFw_->ftabLoHi(*seq_, off, topf, botf);
			#ifdef NDEBUG
			if(botf - topf == 0) return true;
			#endif
			#ifdef NDEBUG
			if(ebwtBw_ != NULL) {
				topb = ebwtBw_->ftabHi(*seq_, off);
				botb = topb + (botf-topf);
			}
			#else
			if(ebwtBw_ != NULL) {
				ebwtBw_->ftabLoHi(*seq_, off, topb, botb);
				assert_eq(botf-topf, botb-topb);
			}
			if(botf - topf == 0) return true;
			#endif
			SA_FTAB(step, ftabLen, true);
			step += ftabLen;
		} else if(s.maxjump > 0) {
			// Use fchr
			int c = (*seq_)[off];
			assert_range(0, 3, c);
			topf = topb = ebwtFw_->fchr()[c];
			botf = botb = ebwtFw_->fchr()[c+1];
			if(botf - topf == 0) return true;
			SA_FCHR(step, true);
			step++;
		} else {
			assert_eq(0, s.maxjump);
			// TODO: handle this gracefully
			throw 1;
		}
		if(step == (int)s.steps.size()) {
			// Finished aligning seed
			assert(c0.acceptable());
			assert(c1.acceptable());
			assert(c2.acceptable());
			if(!reportHit(topf, botf, topb, botb, seq_->length(), prevEdit)) {
				return false; // Memory exhausted
			}
			return true;
		}
		nextLocsBi(tloc, bloc, topf, botf, topb, botb, step);
		assert(tloc.valid());
	} else assert(prevEdit != NULL);
	assert(tloc.valid());
	assert(botf - topf == 1 ||  bloc.valid());
	assert(botf - topf > 1  || !bloc.valid());
	assert_gt(step, 0);
	uint32_t t[4], b[4]; // dest BW ranges
	Constraint* zones[3] = { &c0, &c1, &c2 };
	ASSERT_ONLY(uint32_t lasttot = botf - topf);
	for(int i = step; i < (int)s.steps.size(); i++) {
		assert_gt(botf, topf);
		assert(botf - topf == 1 ||  bloc.valid());
		assert(botf - topf > 1  || !bloc.valid());
		assert(ebwtBw_ == NULL || botf-topf == botb-topb);
		assert(tloc.valid());
		off = s.steps[i];
		bool ltr = off > 0;
		const Ebwt* ebwt = ltr ? ebwtBw_ : ebwtFw_;
		assert(ebwt != NULL);
		if(ltr) {
			tp[0] = tp[1] = tp[2] = tp[3] = topf;
			bp[0] = bp[1] = bp[2] = bp[3] = botf;
		} else {
			tp[0] = tp[1] = tp[2] = tp[3] = topb;
			bp[0] = bp[1] = bp[2] = bp[3] = botb;
		}
		t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
		if(bloc.valid()) {
			// Range delimited by tloc/bloc has size >1.  If size == 1,
			// we use a simpler query (see if(!bloc.valid()) blocks below)
			bwops_++;
			ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
			ASSERT_ONLY(uint32_t tot = (b[0]-t[0])+(b[1]-t[1])+(b[2]-t[2])+(b[3]-t[3]));
			ASSERT_ONLY(uint32_t totp = (bp[0]-tp[0])+(bp[1]-tp[1])+(bp[2]-tp[2])+(bp[3]-tp[3]));
			assert_eq(tot, totp);
			assert_leq(tot, lasttot);
			ASSERT_ONLY(lasttot = tot);
		}
		uint32_t *tf = ltr ? tp : t, *tb = ltr ? t : tp;
		uint32_t *bf = ltr ? bp : b, *bb = ltr ? b : bp;
		off = abs(off)-1;
		//
		bool leaveZone = s.zones[i].first < 0;
		//bool leaveZoneIns = zones_[i].second < 0;
		Constraint& cons    = *zones[abs(s.zones[i].first)];
		Constraint& insCons = *zones[abs(s.zones[i].second)];
		int c = (*seq_)[off];  assert_range(0, 4, c);
		int q = (*qual_)[off];
		// Is it legal for us to advance on characters other than 'c'?
		if(!(cons.mustMatch() && !overall.mustMatch()) || c == 4) {
			// There may be legal edits
			bool bail = false;
			if(!bloc.valid()) {
				// Range delimited by tloc/bloc has size 1
				uint32_t ntop = ltr ? topb : topf;
				bwops_++;
				int cc = ebwt->mapLF1(ntop, tloc);
				assert_range(-1, 3, cc);
				if(cc < 0) bail = true;
				else { t[cc] = ntop; b[cc] = ntop+1; }
			}
			if(!bail) {
				if((cons.canMismatch(q, *sc_) && overall.canMismatch(q, *sc_)) || c == 4) {
					Constraint oldCons = cons, oldOvCons = overall;
					SideLocus oldTloc = tloc, oldBloc = bloc;
					if(c != 4) {
						cons.chargeMismatch(q, *sc_);
						overall.chargeMismatch(q, *sc_);
					}
					// Can leave the zone as-is
					if(!leaveZone || (cons.acceptable() && overall.acceptable())) {
						for(int j = 0; j < 4; j++) {
							if(j == c || b[j] == t[j]) continue;
							// Potential mismatch
							nextLocsBi(tloc, bloc, tf[j], bf[j], tb[j], bb[j], i+1);
							int loff = off;
							if(!ltr) loff = (int)(s.steps.size() - loff - 1);
							assert(prevEdit == NULL || prevEdit->next == NULL);
							Edit edit(off, j, c, EDIT_TYPE_MM, false);
							DoublyLinkedList<Edit> editl;
							editl.payload = edit;
							if(prevEdit != NULL) {
								prevEdit->next = &editl;
								editl.prev = prevEdit;
							}
							assert(editl.next == NULL);
							SA_EDIT(i, ltr, depth, edit);
							bwedits_++;
							if(!searchSeedBi(
								i+1,     // depth into steps_[] array
								depth+1, // recursion depth
								tf[j],   // top in BWT
								bf[j],   // bot in BWT
								tb[j],   // top in BWT'
								bb[j],   // bot in BWT'
								tloc,    // locus for top (perhaps unititialized)
								bloc,    // locus for bot (perhaps unititialized)
								c0,      // constraints to enforce in seed zone 0
								c1,      // constraints to enforce in seed zone 1
								c2,      // constraints to enforce in seed zone 2
								overall, // overall constraints to enforce
								&editl))  // latest edit
							{
								return false;
							}
							if(prevEdit != NULL) prevEdit->next = NULL;
						}
					} else {
						// Not enough edits to make this path
						// non-redundant with other seeds
					}
					cons = oldCons;
					overall = oldOvCons;
					tloc = oldTloc;
					bloc = oldBloc;
				}
				if(cons.canGap() && overall.canGap()) {
					throw 1; // TODO
					int delEx = 0;
					if(cons.canDelete(delEx, *sc_) && overall.canDelete(delEx, *sc_)) {
						// Try delete
					}
					int insEx = 0;
					if(insCons.canInsert(insEx, *sc_) && overall.canInsert(insEx, *sc_)) {
						// Try insert
					}
				}
			} // if(!bail)
		}
		if(c == 4) {
			return true; // couldn't handle the N
		}
		if(leaveZone && (!cons.acceptable() || !overall.acceptable())) {
			// Not enough edits to make this path non-redundant with
			// other seeds
			return true;
		}
		if(!bloc.valid()) {
			assert(ebwtBw_ == NULL || bp[c] == tp[c]+1);
			// Range delimited by tloc/bloc has size 1
			uint32_t top = ltr ? topb : topf;
			bwops_++;
			t[c] = ebwt->mapLF1(top, tloc, c);
			if(t[c] == 0xffffffff) {
				return true;
			}
			assert_geq(t[c], ebwt->fchr()[c]);
			assert_lt(t[c],  ebwt->fchr()[c+1]);
			b[c] = t[c]+1;
			assert_gt(b[c], 0);
		}
		assert(ebwtBw_ == NULL || bf[c]-tf[c] == bb[c]-tb[c]);
		assert_leq(bf[c]-tf[c], lasttot);
		ASSERT_ONLY(lasttot = bf[c]-tf[c]);
		if(b[c] == t[c]) {
			return true;
		}
		SA_MATCH(i, ltr, depth);
		topf = tf[c]; botf = bf[c];
		topb = tb[c]; botb = bb[c];
		if(i+1 == (int)s.steps.size()) {
			// Finished aligning seed
			assert(c0.acceptable());
			assert(c1.acceptable());
			assert(c2.acceptable());
			if(!reportHit(topf, botf, topb, botb, seq_->length(), prevEdit)) {
				return false; // Memory exhausted
			}
			return true;
		}
		nextLocsBi(tloc, bloc, tf[c], bf[c], tb[c], bb[c], i+1);
	}
	return true;
}

#ifdef ALIGNER_SEED_MAIN

#include <getopt.h>
#include <string>

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(const char *errmsg, const char *arg) {
	long l;
	char *endPtr = NULL;
	l = strtol(arg, &endPtr, 10);
	if (endPtr != NULL) {
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	throw 1;
	return -1;
}

enum {
	ARG_NOFW = 256,
	ARG_NORC,
	ARG_MM,
	ARG_SHMEM,
	ARG_TESTS,
	ARG_RANDOM_TESTS,
	ARG_SEED
};

static const char *short_opts = "vCt";
static struct option long_opts[] = {
	{(char*)"verbose",  no_argument,       0, 'v'},
	{(char*)"color",    no_argument,       0, 'C'},
	{(char*)"timing",   no_argument,       0, 't'},
	{(char*)"nofw",     no_argument,       0, ARG_NOFW},
	{(char*)"norc",     no_argument,       0, ARG_NORC},
	{(char*)"mm",       no_argument,       0, ARG_MM},
	{(char*)"shmem",    no_argument,       0, ARG_SHMEM},
	{(char*)"tests",    no_argument,       0, ARG_TESTS},
	{(char*)"random",   required_argument, 0, ARG_RANDOM_TESTS},
	{(char*)"seed",     required_argument, 0, ARG_SEED},
};

static void printUsage(ostream& os) {
	os << "Usage: ac [options]* <index> <patterns>" << endl;
	os << "Options:" << endl;
	os << "  --mm                memory-mapped mode" << endl;
	os << "  --shmem             shared memory mode" << endl;
	os << "  --nofw              don't align forward-oriented read" << endl;
	os << "  --norc              don't align reverse-complemented read" << endl;
	os << "  -t/--timing         show timing information" << endl;
	os << "  -C/--color          colorspace mode" << endl;
	os << "  -v/--verbose        talkative mode" << endl;
}

bool gNorc = false;
bool gNofw = false;
bool gColor = false;
int gVerbose = 0;
int gGapBarrier = 1;
bool gColorExEnds = true;
int gSnpPhred = 30;
bool gReportOverhangs = true;

extern void aligner_seed_tests();
extern void aligner_random_seed_tests(
	int num_tests,
	uint32_t qslo,
	uint32_t qshi,
	uint32_t maxrowslo,
	uint32_t maxrowshi,
	uint32_t maxposslo,
	uint32_t maxposshi,
	bool color,
	uint32_t seed);

/**
 * A way of feeding simply tests to the seed alignment infrastructure.
 */
int main(int argc, char **argv) {
	bool useMm = false;
	bool useShmem = false;
	bool mmSweep = false;
	bool noRefNames = false;
	bool sanity = false;
	bool timing = false;
	int option_index = 0;
	int seed = 777;
	int next_option;
	do {
		next_option = getopt_long(
			argc, argv, short_opts, long_opts, &option_index);
		switch (next_option) {
			case 'v':       gVerbose = true; break;
			case 'C':       gColor   = true; break;
			case 't':       timing   = true; break;
			case ARG_NOFW:  gNofw    = true; break;
			case ARG_NORC:  gNorc    = true; break;
			case ARG_MM:    useMm    = true; break;
			case ARG_SHMEM: useShmem = true; break;
			case ARG_SEED:  seed = parseInt("", optarg); break;
			case ARG_TESTS: {
				aligner_seed_tests();
				aligner_random_seed_tests(
					100,     // num references
					100,   // queries per reference lo
					400,   // queries per reference hi
					1,     // maxrows lo
					100,   // maxrows hi
					1,     // maxposs lo
					100,   // maxposs hi
					false, // true -> generate colorspace reference/reads
					18);   // pseudo-random seed
				return 0;
			}
			case ARG_RANDOM_TESTS: {
				seed = parseInt("", optarg);
				aligner_random_seed_tests(
					100,   // num references
					100,   // queries per reference lo
					400,   // queries per reference hi
					1,     // maxrows lo
					100,   // maxrows hi
					1,     // maxposs lo
					100,   // maxposs hi
					false, // true -> generate colorspace reference/reads
					seed); // pseudo-random seed
				return 0;
			}
			case -1: break;
			default: {
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
			}
		}
	} while(next_option != -1);
	char *reffn;
	if(optind >= argc) {
		cerr << "No reference; quitting..." << endl;
		return 1;
	}
	reffn = argv[optind++];
	if(optind >= argc) {
		cerr << "No reads; quitting..." << endl;
		return 1;
	}
	string ebwtBase(reffn);
	BitPairReference ref(
		ebwtBase,    // base path
		gColor,      // whether we expect it to be colorspace
		sanity,      // whether to sanity-check reference as it's loaded
		NULL,        // fasta files to sanity check reference against
		NULL,        // another way of specifying original sequences
		false,       // true -> infiles (2 args ago) contains raw seqs
		useMm,       // use memory mapping to load index?
		useShmem,    // use shared memory (not memory mapping)
		mmSweep,     // touch all the pages after memory-mapping the index
		gVerbose,    // verbose
		gVerbose);   // verbose but just for startup messages
	Timer *t = new Timer(cerr, "Time loading fw index: ", timing);
	Ebwt ebwtFw(
		ebwtBase,
		gColor,      // index is colorspace
		0,           // don't need entireReverse for fw index
		true,        // index is for the forward direction
		5,           // offrate (irrelevant)
		useMm,       // whether to use memory-mapped files
		useShmem,    // whether to use shared memory
		mmSweep,     // sweep memory-mapped files
		!noRefNames, // load names?
		false,       // load SA sample?
		true,        // load ftab?
		true,        // load rstarts?
		NULL,        // reference map, or NULL if none is needed
		gVerbose,    // whether to be talkative
		gVerbose,    // talkative during initialization
		false,       // handle memory exceptions, don't pass them up
		sanity);
	delete t;
	t = new Timer(cerr, "Time loading bw index: ", timing);
	Ebwt ebwtBw(
		ebwtBase + ".rev",
		gColor,      // index is colorspace
		1,           // need entireReverse
		false,       // index is for the backward direction
		5,           // offrate (irrelevant)
		useMm,       // whether to use memory-mapped files
		useShmem,    // whether to use shared memory
		mmSweep,     // sweep memory-mapped files
		!noRefNames, // load names?
		false,       // load SA sample?
		true,        // load ftab?
		false,       // load rstarts?
		NULL,        // reference map, or NULL if none is needed
		gVerbose,    // whether to be talkative
		gVerbose,    // talkative during initialization
		false,       // handle memory exceptions, don't pass them up
		sanity);
	delete t;
	for(int i = optind; i < argc; i++) {
	}
}
#endif
