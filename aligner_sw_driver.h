/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 *  aligner_sw_driver.h
 *
 * REDUNDANT SEED HITS
 *
 * We say that two seed hits are redundant if they trigger identical
 * seed-extend dynamic programming problems.  Put another way, they both lie on
 * the same diagonal of the overall read/reference dynamic programming matrix.
 * Detecting redundant seed hits is simple when the seed hits are ungapped.  We
 * do this after offset resolution but before the offset is converted to genome
 * coordinates (see uses of the seenDiags1_/seenDiags2_ fields for examples).
 *
 * REDUNDANT ALIGNMENTS
 *
 * In an unpaired context, we say that two alignments are redundant if they
 * share any cells in the global DP table.  Roughly speaking, this is like
 * saying that two alignments are redundant if any read character aligns to the
 * same reference character (same reference sequence, same strand, same offset)
 * in both alignments.
 *
 * In a paired-end context, we say that two paired-end alignments are redundant
 * if the mate #1s are redundant and the mate #2s are redundant.
 *
 * How do we enforce this?  In the unpaired context, this is relatively simple:
 * the cells from each alignment are checked against a set containing all cells
 * from all previous alignments.  Given a new alignment, for each cell in the
 * new alignment we check whether it is in the set.  If there is any overlap,
 * the new alignment is rejected as redundant.  Otherwise, the new alignment is
 * accepted and its cells are added to the set.
 *
 * Enforcement in a paired context is a little trickier.  Consider the
 * following approaches:
 *
 * 1. Skip anchors that are redundant with any previous anchor or opposite
 *    alignment.  This is sufficient to ensure no two concordant alignments
 *    found are redundant.
 *
 * 2. Same as scheme 1, but with a "transitive closure" scheme for finding all
 *    concordant pairs in the vicinity of an anchor.  Consider the AB/AC
 *    scenario from the previous paragraph.  If B is the anchor alignment, we
 *    will find AB but not AC.  But under this scheme, once we find AB we then
 *    let B be a new anchor and immediately look for its opposites.  Likewise,
 *    if we find any opposite, we make them anchors and continue searching.  We
 *    don't stop searching until every opposite is used as an anchor.
 *
 * 3. Skip anchors that are redundant with any previous anchor alignment (but
 *    allow anchors that are redundant with previous opposite alignments).
 *    This isn't sufficient to avoid redundant concordant alignments.  To avoid
 *    redundant concordants, we need an additional procedure that checks each
 *    new concordant alignment one-by-one against a list of previous concordant
 *    alignments to see if it is redundant.
 *
 * We take approach 1.
 */

#ifndef ALIGNER_SW_DRIVER_H_
#define ALIGNER_SW_DRIVER_H_

#include <utility>
#include "ds.h"
#include "aligner_seed.h"
#include "aligner_sw.h"
#include "aligner_cache.h"
#include "reference.h"
#include "group_walk.h"
#include "bt2_idx.h"
#include "mem_ids.h"
#include "aln_sink.h"
#include "pe.h"
#include "ival_list.h"
#include "simple_func.h"
#include "random_util.h"

struct SeedPos {

	SeedPos() : fw(false), offidx(0), rdoff(0), seedlen(0) { }

	SeedPos(
		bool fw_,
		uint32_t offidx_,
		uint32_t rdoff_,
		uint32_t seedlen_)
	{
		init(fw_, offidx_, rdoff_, seedlen_);
	}
	
	void init(
		bool fw_,
		uint32_t offidx_,
		uint32_t rdoff_,
		uint32_t seedlen_)
	{
		fw      = fw_;
		offidx  = offidx_;
		rdoff   = rdoff_;
		seedlen = seedlen_;
	}
	
	bool operator<(const SeedPos& o) const {
		if(offidx < o.offidx)   return true;
		if(offidx > o.offidx)   return false;
		if(rdoff < o.rdoff)     return true;
		if(rdoff > o.rdoff)     return false;
		if(seedlen < o.seedlen) return true;
		if(seedlen > o.seedlen) return false;
		if(fw && !o.fw)         return true;
		if(!fw && o.fw)         return false;
		return false;
	}
	
	bool operator>(const SeedPos& o) const {
		if(offidx < o.offidx)   return false;
		if(offidx > o.offidx)   return true;
		if(rdoff < o.rdoff)     return false;
		if(rdoff > o.rdoff)     return true;
		if(seedlen < o.seedlen) return false;
		if(seedlen > o.seedlen) return true;
		if(fw && !o.fw)         return false;
		if(!fw && o.fw)         return true;
		return false;
	}
	
	bool operator==(const SeedPos& o) const {
		return fw == o.fw && offidx == o.offidx &&
		       rdoff == o.rdoff && seedlen == o.seedlen;
	}

	bool fw;
	uint32_t offidx;
	uint32_t rdoff;
	uint32_t seedlen;
};

/**
 * An SATuple along with the associated seed position.
 */
struct SATupleAndPos {
	
	SATuple sat;    // result for this seed hit
	SeedPos pos;    // seed position that yielded the range this was taken from
	size_t  origSz; // size of range this was taken from
	size_t  nlex;   // # position we can extend seed hit to left w/o edit
	size_t  nrex;   // # position we can extend seed hit to right w/o edit
	
	bool operator<(const SATupleAndPos& o) const {
		if(sat < o.sat) return true;
		if(sat > o.sat) return false;
		return pos < o.pos;
	}

	bool operator==(const SATupleAndPos& o) const {
		return sat == o.sat && pos == o.pos;
	}
};

/**
 * Encapsulates the weighted random sampling scheme we want to use to pick
 * which seed hit range to sample a row from.
 */
class RowSampler {

public:

	RowSampler(int cat = 0) : elim_(cat), masses_(cat) { 
		mass_ = 0.0f;
	}
	
	/**
	 * Initialze sampler with respect to a range of elements in a list of
	 * SATupleAndPos's.
	 */
	void init(
		const EList<SATupleAndPos, 16>& salist,
		size_t sai,
		size_t saf,
		bool lensq, // whether to square the numerator, which = extended length
		bool szsq)  // whether to square denominator, which = 
	{
		assert_gt(saf, sai);
		elim_.resize(saf - sai);
		elim_.fill(false);
		// Initialize mass
		mass_ = 0.0f;
		masses_.resize(saf - sai);
		for(size_t i = sai; i < saf; i++) {
			size_t len = salist[i].nlex + salist[i].nrex + 1; // + salist[i].sat.key.len;
			double num = (double)len;
			if(lensq) {
				num *= num;
			}
			double denom = (double)salist[i].sat.size();
			if(szsq) {
				denom *= denom;
			}
			masses_[i - sai] = num / denom;
			mass_ += masses_[i - sai];
		}
	}
	
	/**
	 * Caller is indicating that the bin at index i is exhausted and we should
	 * exclude it from our sampling from now on.
	 */
	void finishedRange(size_t i) {
		assert_lt(i, masses_.size());
		elim_[i] = true;
		mass_ -= masses_[i];
	}
	
	/**
	 * Sample randomly from the mass.
	 */
	size_t next(RandomSource& rnd) {
		// Throw the dart
		double rd = rnd.nextFloat() * mass_;
		double mass_sofar = 0.0f;
		size_t sz = masses_.size();
		size_t last_unelim = std::numeric_limits<size_t>::max();
		for(size_t i = 0; i < sz; i++) {
			if(!elim_[i]) {
				last_unelim = i;
				mass_sofar += masses_[i];
				if(rd < mass_sofar) {
					// This is the one we hit
					return i;
				}
			}
		}
		assert_neq(std::numeric_limits<size_t>::max(), last_unelim);
		return last_unelim;
	}

protected:
	double        mass_;    // total probability mass to throw darts at
	EList<bool>   elim_;    // whether the range is eliminated
	EList<double> masses_;  // mass of each range
};

/**
 * Return values from extendSeeds and extendSeedsPaired.
 */
enum {
	// All end-to-end and seed hits were examined
	// The policy does not need us to look any further
	EXTEND_EXHAUSTED_CANDIDATES = 1,
	EXTEND_POLICY_FULFILLED,
	// We stopped because we reached a point where the only remaining
	// alignments of interest have perfect scores, but we already investigated
	// perfect alignments
	EXTEND_PERFECT_SCORE,
	// We stopped because we ran up against a limit on how much work we should
	// do for one set of seed ranges, e.g. the limit on number of consecutive
	// unproductive DP extensions
	EXTEND_EXCEEDED_SOFT_LIMIT,
	// We stopped because we ran up against a limit on how much work we should
	// do for overall before giving up on a mate
	EXTEND_EXCEEDED_HARD_LIMIT
};

/**
 * Data structure encapsulating a range that's been extended out in two
 * directions.
 */
struct ExtendRange {

	void init(size_t off_, size_t len_, size_t sz_) {
		off = off_; len = len_; sz = sz_;
	}

	size_t off; // offset of extended region
	size_t len; // length between extremes of extended region
	size_t sz;  // # of elements in SA range
};

class SwDriver {

	typedef PList<uint32_t, CACHE_PAGE_SZ> TSAList;

public:

	SwDriver(size_t bytes) :
		satups_(DP_CAT),
		gws_(DP_CAT),
		seenDiags1_(DP_CAT),
		seenDiags2_(DP_CAT),
		redAnchor_(DP_CAT),
		redMate1_(DP_CAT),
		redMate2_(DP_CAT),
		pool_(bytes, CACHE_PAGE_SZ, DP_CAT),
		salistEe_(DP_CAT),
		gwstate_(GW_CAT) { }

	/**
	 * Given a collection of SeedHits for a single read, extend seed alignments
	 * into full alignments.  Where possible, try to avoid redundant offset
	 * lookups and dynamic programming problems.  Optionally report alignments
	 * to a AlnSinkWrap object as they are discovered.
	 *
	 * If 'reportImmediately' is true, returns true iff a call to
	 * mhs->report() returned true (indicating that the reporting
	 * policy is satisfied and we can stop).  Otherwise, returns false.
	 */
	int extendSeeds(
		Read& rd,                    // read to align
		bool mate1,                  // true iff rd is mate #1
		SeedResults& sh,             // seed hits to extend into full alignments
		const Ebwt& ebwtFw,          // BWT
		const Ebwt* ebwtBw,          // BWT'
		const BitPairReference& ref, // Reference strings
		SwAligner& swa,              // dynamic programming aligner
		const Scoring& sc,           // scoring scheme
		int seedmms,                 // # mismatches allowed in seed
		int seedlen,                 // length of seed
		int seedival,                // interval between seeds
		TAlScore& minsc,             // minimum score for anchor
		int nceil,                   // maximum # Ns permitted in ref portion
		size_t maxhalf,              // maximum width on one side of DP table
		bool doUngapped,             // do ungapped alignment
		size_t maxIters,             // stop after this many seed-extend loop iters
		size_t maxUg,                // max # ungapped extends
		size_t maxDp,                // max # DPs
		size_t maxUgStreak,          // stop after streak of this many ungap fails
		size_t maxDpStreak,          // stop after streak of this many dp fails
		bool doExtend,               // do seed extension
		bool enable8,                // use 8-bit SSE where possible
		size_t cminlen,              // use checkpointer if read longer than this
		size_t cpow2,                // interval between diagonals to checkpoint
		bool doTri,                  // triangular mini-fills
		int tighten,                 // -M score tightening mode
		AlignmentCacheIface& ca,     // alignment cache for seed hits
		RandomSource& rnd,           // pseudo-random source
		WalkMetrics& wlm,            // group walk left metrics
		SwMetrics& swmSeed,          // DP metrics for seed-extend
		PerReadMetrics& prm,         // per-read metrics
		AlnSinkWrap* mhs,            // HitSink for multiseed-style aligner
		bool reportImmediately,      // whether to report hits immediately to mhs
		bool& exhaustive);

	/**
	 * Given a collection of SeedHits for a read pair, extend seed
	 * alignments into full alignments and then look for the opposite
	 * mate using dynamic programming.  Where possible, try to avoid
	 * redundant offset lookups.  Optionally report alignments to a
	 * AlnSinkWrap object as they are discovered.
	 *
	 * If 'reportImmediately' is true, returns true iff a call to
	 * mhs->report() returned true (indicating that the reporting
	 * policy is satisfied and we can stop).  Otherwise, returns false.
	 */
	int extendSeedsPaired(
		Read& rd,                    // mate to align as anchor
		Read& ord,                   // mate to align as opposite
		bool anchor1,                // true iff anchor mate is mate1
		bool oppFilt,                // true iff opposite mate was filtered out
		SeedResults& sh,             // seed hits for anchor
		const Ebwt& ebwtFw,          // BWT
		const Ebwt* ebwtBw,          // BWT'
		const BitPairReference& ref, // Reference strings
		SwAligner& swa,              // dyn programming aligner for anchor
		SwAligner& swao,             // dyn programming aligner for opposite
		const Scoring& sc,           // scoring scheme
		const PairedEndPolicy& pepol,// paired-end policy
		int seedmms,                 // # mismatches allowed in seed
		int seedlen,                 // length of seed
		int seedival,                // interval between seeds
		TAlScore& minsc,             // minimum score for anchor
		TAlScore& ominsc,            // minimum score for opposite
		int nceil,                   // max # Ns permitted in ref for anchor
		int onceil,                  // max # Ns permitted in ref for opposite
		bool nofw,                   // don't align forward read
		bool norc,                   // don't align revcomp read
		size_t maxhalf,              // maximum width on one side of DP table
		bool doUngapped,             // do ungapped alignment
		size_t maxIters,             // stop after this many seed-extend loop iters
		size_t maxUg,                // max # ungapped extends
		size_t maxDp,                // max # DPs
		size_t maxEeStreak,          // stop after streak of this many end-to-end fails
		size_t maxUgStreak,          // stop after streak of this many ungap fails
		size_t maxDpStreak,          // stop after streak of this many dp fails
		size_t maxMateStreak,        // stop seed range after N mate-find fails
		bool doExtend,               // do seed extension
		bool enable8,                // use 8-bit SSE where possible
		size_t cminlen,              // use checkpointer if read longer than this
		size_t cpow2,                // interval between diagonals to checkpoint
		bool doTri,                  // triangular mini-fills
		int tighten,                 // -M score tightening mode
		AlignmentCacheIface& cs,     // alignment cache for seed hits
		RandomSource& rnd,           // pseudo-random source
		WalkMetrics& wlm,            // group walk left metrics
		SwMetrics& swmSeed,          // DP metrics for seed-extend
		SwMetrics& swmMate,          // DP metrics for mate finidng
		PerReadMetrics& prm,         // per-read metrics for anchor
		AlnSinkWrap* msink,          // AlnSink wrapper for multiseed-style aligner
		bool swMateImmediately,      // whether to look for mate immediately
		bool reportImmediately,      // whether to report hits immediately to msink
		bool discord,                // look for discordant alignments?
		bool mixed,                  // look for unpaired as well as paired alns?
		bool& exhaustive);

	/**
	 * Prepare for a new read.
	 */
	void nextRead(bool paired, size_t mate1len, size_t mate2len) {
		redAnchor_.reset();
		seenDiags1_.reset();
		seenDiags2_.reset();
		seedExRangeFw_[0].clear(); // mate 1 fw
		seedExRangeFw_[1].clear(); // mate 2 fw
		seedExRangeRc_[0].clear(); // mate 1 rc
		seedExRangeRc_[1].clear(); // mate 2 rc
		size_t maxlen = mate1len;
		if(paired) {
			redMate1_.reset();
			redMate1_.init(mate1len);
			redMate2_.reset();
			redMate2_.init(mate2len);
			if(mate2len > maxlen) {
				maxlen = mate2len;
			}
		}
		redAnchor_.init(maxlen);
	}

protected:

	bool eeSaTups(
		const Read& rd,              // read
		SeedResults& sh,             // seed hits to extend into full alignments
		const Ebwt& ebwt,            // BWT
		const BitPairReference& ref, // Reference strings
		RandomSource& rnd,           // pseudo-random generator
		WalkMetrics& wlm,            // group walk left metrics
		SwMetrics& swmSeed,          // metrics for seed extensions
		size_t& nelt_out,            // out: # elements total
        size_t maxelts,              // max # elts to report
		bool all);                   // report all hits?

	void extend(
		const Read& rd,       // read
		const Ebwt& ebwtFw,   // Forward Bowtie index
		const Ebwt* ebwtBw,   // Backward Bowtie index
		uint32_t topf,        // top in fw index
		uint32_t botf,        // bot in fw index
		uint32_t topb,        // top in bw index
		uint32_t botb,        // bot in bw index
		bool fw,              // seed orientation
		size_t off,           // seed offset from 5' end
		size_t len,           // seed length
		PerReadMetrics& prm,  // per-read metrics
		size_t& nlex,         // # positions we can extend to left w/o edit
		size_t& nrex);        // # positions we can extend to right w/o edit

	void prioritizeSATups(
		const Read& rd,              // read
		SeedResults& sh,             // seed hits to extend into full alignments
		const Ebwt& ebwtFw,          // BWT
		const Ebwt* ebwtBw,          // BWT'
		const BitPairReference& ref, // Reference strings
		int seedmms,                 // # seed mismatches allowed
		size_t maxelt,               // max elts we'll consider
		bool doExtend,               // extend out seeds
		bool lensq,                  // square extended length
		bool szsq,                   // square SA range size
		size_t nsm,                  // if range as <= nsm elts, it's "small"
		AlignmentCacheIface& ca,     // alignment cache for seed hits
		RandomSource& rnd,           // pseudo-random generator
		WalkMetrics& wlm,            // group walk left metrics
		PerReadMetrics& prm,         // per-read metrics
		size_t& nelt_out,            // out: # elements total
		bool all);                   // report all hits?

	Random1toN               rand_;    // random number generators
	EList<Random1toN, 16>    rands_;   // random number generators
	EList<Random1toN, 16>    rands2_;  // random number generators
	EList<EEHit, 16>         eehits_;  // holds end-to-end hits
	EList<SATupleAndPos, 16> satpos_;  // holds SATuple, SeedPos pairs
	EList<SATupleAndPos, 16> satpos2_; // holds SATuple, SeedPos pairs
	EList<SATuple, 16>       satups_;  // holds SATuples to explore elements from
	EList<GroupWalk2S<TSlice, 16> > gws_;   // list of GroupWalks; no particular order
	EList<size_t>            mateStreaks_; // mate-find fail streaks
	RowSampler               rowsamp_;     // row sampler
	
	// Ranges that we've extended through when extending seed hits
	EList<ExtendRange> seedExRangeFw_[2];
	EList<ExtendRange> seedExRangeRc_[2];

	// Data structures encapsulating the diagonals that have already been used
	// to seed alignment for mate 1 and mate 2.
	EIvalMergeListBinned seenDiags1_;
	EIvalMergeListBinned seenDiags2_;

	// For weeding out redundant alignments
	RedundantAlns  redAnchor_;  // database of cells used for anchor alignments
	RedundantAlns  redMate1_;   // database of cells used for mate 1 alignments
	RedundantAlns  redMate2_;   // database of cells used for mate 2 alignments

	// For holding results for anchor (res_) and opposite (ores_) mates
	SwResult       resGap_;    // temp holder for alignment result
	SwResult       oresGap_;   // temp holder for alignment result, opp mate
	SwResult       resUngap_;  // temp holder for ungapped alignment result
	SwResult       oresUngap_; // temp holder for ungap. aln. opp mate
	SwResult       resEe_;     // temp holder for ungapped alignment result
	SwResult       oresEe_;    // temp holder for ungap. aln. opp mate
	
	Pool           pool_;      // memory pages for salistExact_
	TSAList        salistEe_;  // PList for offsets for end-to-end hits
	GroupWalkState gwstate_;   // some per-thread state shared by all GroupWalks
	
	// For AlnRes::matchesRef:
	ASSERT_ONLY(SStringExpandable<char>     raw_refbuf_);
	ASSERT_ONLY(SStringExpandable<uint32_t> raw_destU32_);
	ASSERT_ONLY(EList<bool>                 raw_matches_);
	ASSERT_ONLY(BTDnaString                 tmp_rf_);
	ASSERT_ONLY(BTDnaString                 tmp_rdseq_);
	ASSERT_ONLY(BTString                    tmp_qseq_);
};

#endif /*ALIGNER_SW_DRIVER_H_*/
