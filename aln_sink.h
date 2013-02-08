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

#ifndef ALN_SINK_H_
#define ALN_SINK_H_

#include <limits>
#include "read.h"
#include "unique.h"
#include "sam.h"
#include "ds.h"
#include "simple_func.h"
#include "outq.h"
#include <utility>

// Forward decl	
class SeedResults;

enum {
	OUTPUT_SAM = 1
};

/**
 * Metrics summarizing the work done by the reporter and summarizing
 * the number of reads that align, that fail to align, and that align
 * non-uniquely.
 */
struct ReportingMetrics {

	ReportingMetrics():mutex_m() {
	    reset();
	}

	void reset() {
		init(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	}

	void init(
		uint64_t nread_,
		uint64_t npaired_,
		uint64_t nunpaired_,
		uint64_t nconcord_uni_,
		uint64_t nconcord_uni1_,
		uint64_t nconcord_uni2_,
		uint64_t nconcord_rep_,
		uint64_t nconcord_0_,
		uint64_t ndiscord_,
		uint64_t nunp_0_uni_,
		uint64_t nunp_0_uni1_,
		uint64_t nunp_0_uni2_,
		uint64_t nunp_0_rep_,
		uint64_t nunp_0_0_,
		uint64_t nunp_rep_uni_,
		uint64_t nunp_rep_uni1_,
		uint64_t nunp_rep_uni2_,
		uint64_t nunp_rep_rep_,
		uint64_t nunp_rep_0_,
		uint64_t nunp_uni_,
		uint64_t nunp_uni1_,
		uint64_t nunp_uni2_,
		uint64_t nunp_rep_,
		uint64_t nunp_0_,
		uint64_t sum_best1_,
		uint64_t sum_best2_,
		uint64_t sum_best_)
	{
		nread         = nread_;
		
		npaired       = npaired_;
		nunpaired     = nunpaired_;
		
		nconcord_uni  = nconcord_uni_;
		nconcord_uni1 = nconcord_uni1_;
		nconcord_uni2 = nconcord_uni2_;
		nconcord_rep  = nconcord_rep_;
		nconcord_0    = nconcord_0_;
		
		ndiscord      = ndiscord_;
		
		nunp_0_uni    = nunp_0_uni_;
		nunp_0_uni1   = nunp_0_uni1_;
		nunp_0_uni2   = nunp_0_uni2_;
		nunp_0_rep    = nunp_0_rep_;
		nunp_0_0      = nunp_0_0_;

		nunp_rep_uni  = nunp_rep_uni_;
		nunp_rep_uni1 = nunp_rep_uni1_;
		nunp_rep_uni2 = nunp_rep_uni2_;
		nunp_rep_rep  = nunp_rep_rep_;
		nunp_rep_0    = nunp_rep_0_;

		nunp_uni      = nunp_uni_;
		nunp_uni1     = nunp_uni1_;
		nunp_uni2     = nunp_uni2_;
		nunp_rep      = nunp_rep_;
		nunp_0        = nunp_0_;

		sum_best1     = sum_best1_;
		sum_best2     = sum_best2_;
		sum_best      = sum_best_;
	}
	
	/**
	 * Merge (add) the counters in the given ReportingMetrics object
	 * into this object.  This is the only safe way to update a
	 * ReportingMetrics shared by multiple threads.
	 */
	void merge(const ReportingMetrics& met, bool getLock = false) {
        ThreadSafe ts(&mutex_m, getLock);
		nread         += met.nread;

		npaired       += met.npaired;
		nunpaired     += met.nunpaired;

		nconcord_uni  += met.nconcord_uni;
		nconcord_uni1 += met.nconcord_uni1;
		nconcord_uni2 += met.nconcord_uni2;
		nconcord_rep  += met.nconcord_rep;
		nconcord_0    += met.nconcord_0;

		ndiscord      += met.ndiscord;

		nunp_0_uni    += met.nunp_0_uni;
		nunp_0_uni1   += met.nunp_0_uni1;
		nunp_0_uni2   += met.nunp_0_uni2;
		nunp_0_rep    += met.nunp_0_rep;
		nunp_0_0      += met.nunp_0_0;

		nunp_rep_uni  += met.nunp_rep_uni;
		nunp_rep_uni1 += met.nunp_rep_uni1;
		nunp_rep_uni2 += met.nunp_rep_uni2;
		nunp_rep_rep  += met.nunp_rep_rep;
		nunp_rep_0    += met.nunp_rep_0;

		nunp_uni      += met.nunp_uni;
		nunp_uni1     += met.nunp_uni1;
		nunp_uni2     += met.nunp_uni2;
		nunp_rep      += met.nunp_rep;
		nunp_0        += met.nunp_0;

		sum_best1     += met.sum_best1;
		sum_best2     += met.sum_best2;
		sum_best      += met.sum_best;
	}

	uint64_t  nread;         // # reads
	uint64_t  npaired;       // # pairs
	uint64_t  nunpaired;     // # unpaired reads
	
	// Paired
	
	//  Concordant
	uint64_t  nconcord_uni;  // # pairs with unique concordant alns
	uint64_t  nconcord_uni1; // # pairs with exactly 1 concordant alns
	uint64_t  nconcord_uni2; // # pairs with >1 concordant aln, still unique
	uint64_t  nconcord_rep;  // # pairs with repetitive concordant alns
	uint64_t  nconcord_0;    // # pairs with 0 concordant alns
	//  Discordant
	uint64_t  ndiscord;      // # pairs with 1 discordant aln
	
	//  Unpaired from failed pairs
	uint64_t  nunp_0_uni;    // # unique from nconcord_0_ - ndiscord_
	uint64_t  nunp_0_uni1;   // # pairs with exactly 1 concordant alns
	uint64_t  nunp_0_uni2;   // # pairs with >1 concordant aln, still unique
	uint64_t  nunp_0_rep;    // # repetitive from 
	uint64_t  nunp_0_0;      // # with 0 alignments

	//  Unpaired from repetitive pairs
	uint64_t  nunp_rep_uni;  // # pairs with unique concordant alns
	uint64_t  nunp_rep_uni1; // # pairs with exactly 1 concordant alns
	uint64_t  nunp_rep_uni2; // # pairs with >1 concordant aln, still unique
	uint64_t  nunp_rep_rep;  // # pairs with repetitive concordant alns
	uint64_t  nunp_rep_0;    // # pairs with 0 concordant alns
	
	// Unpaired
	
	uint64_t  nunp_uni;      // # unique from nconcord_0_ - ndiscord_
	uint64_t  nunp_uni1;     // # pairs with exactly 1 concordant alns
	uint64_t  nunp_uni2;     // # pairs with >1 concordant aln, still unique
	uint64_t  nunp_rep;      // # repetitive from 
	uint64_t  nunp_0;        // # with 0 alignments

	
	uint64_t  sum_best1;     // Sum of all the best alignment scores
	uint64_t  sum_best2;     // Sum of all the second-best alignment scores
	uint64_t  sum_best;      // Sum of all the best and second-best

	MUTEX_T mutex_m;
};

// Type for expression numbers of hits
typedef int64_t THitInt;

/**
 * Parameters affecting reporting of alignments, specifically -k & -a,
 * -m & -M.
 */
struct ReportingParams {

	explicit ReportingParams(
		THitInt khits_,
		THitInt mhits_,
		THitInt pengap_,
		bool msample_,
		bool discord_,
		bool mixed_)
	{
		init(khits_, mhits_, pengap_, msample_, discord_, mixed_);
	}

	void init(
		THitInt khits_,
		THitInt mhits_,
		THitInt pengap_,
		bool msample_,
		bool discord_,
		bool mixed_)
	{
		khits   = khits_;     // -k (or high if -a)
		mhits   = ((mhits_ == 0) ? std::numeric_limits<THitInt>::max() : mhits_);
		pengap  = pengap_;
		msample = msample_;
		discord = discord_;
		mixed   = mixed_;
	}
	
#ifndef NDEBUG
	/**
	 * Check that reporting parameters are internally consistent.
	 */
	bool repOk() const {
		assert_geq(khits, 1);
		assert_geq(mhits, 1);
		return true;
	}
#endif
	
	/**
	 * Return true iff a -m or -M limit was set by the user.
	 */
	inline bool mhitsSet() const {
		return mhits < std::numeric_limits<THitInt>::max();
	}
	
	/**
	 * Return a multiplier that indicates how many alignments we might look for
	 * (max).  We can use this to boost parameters like ROWM and POSF
	 * appropriately.
	 */
	inline THitInt mult() const {
		if(mhitsSet()) {
			return mhits+1;
		}
		return khits;
	}

	/**
	 * Given ROWM, POSF thresholds, boost them according to mult().
	 */
	void boostThreshold(SimpleFunc& func) {
		THitInt mul = mult();
		assert_gt(mul, 0);
		if(mul == std::numeric_limits<THitInt>::max()) {
			func.setMin(std::numeric_limits<double>::max());
		} else if(mul > 1) {
			func.mult(mul);
		}
	}
	
	/**
	 * Return true iff we are reporting all hits.
	 */
	bool allHits() const {
		return khits == std::numeric_limits<THitInt>::max();
	}

	// Number of alignments to report
	THitInt khits;
	
	// Read is non-unique if mhits-1 next-best alignments are within
	// pengap of the best alignment
	THitInt mhits, pengap;
	
	// true if -M is specified, meaning that if the -M ceiling is
	// exceeded, we should report 'khits' alignments chosen at random
	// from those found
	bool msample;
	
	// true iff we should seek and report discordant paired-end alignments for
	// paired-end reads.
	bool discord;

	// true iff we should seek and report unpaired mate alignments when there
	// are paired-end alignments for a paired-end read, or if the number of
	// paired-end alignments exceeds the -m ceiling.
	bool mixed;
};

/**
 * A state machine keeping track of the number and type of alignments found so
 * far.  Its purpose is to inform the caller as to what stage the alignment is
 * in and what categories of alignment are still of interest.  This information
 * should allow the caller to short-circuit some alignment work.  Another
 * purpose is to tell the AlnSinkWrap how many and what type of alignment to
 * report.
 *
 * TODO: This class does not keep accurate information about what
 * short-circuiting took place.  If a read is identical to a previous read,
 * there should be a way to query this object to determine what work, if any,
 * has to be re-done for the new read.
 */
class ReportingState {

public:

	enum {
		NO_READ = 1,        // haven't got a read yet
		CONCORDANT_PAIRS,   // looking for concordant pairs
		DISCORDANT_PAIRS,   // looking for discordant pairs
		UNPAIRED,           // looking for unpaired
		DONE                // finished looking
	};

	// Flags for different ways we can finish out a category of potential
	// alignments.
	
	enum {
		EXIT_DID_NOT_EXIT = 1,        // haven't finished
		EXIT_DID_NOT_ENTER,           // never tried search	
		EXIT_SHORT_CIRCUIT_k,         // -k exceeded
		EXIT_SHORT_CIRCUIT_M,         // -M exceeded
		EXIT_SHORT_CIRCUIT_TRUMPED,   // made irrelevant
		EXIT_CONVERTED_TO_DISCORDANT, // unpair became discord
		EXIT_NO_ALIGNMENTS,           // none found
		EXIT_WITH_ALIGNMENTS          // some found
	};
	
	ReportingState(const ReportingParams& p) : p_(p) { reset(); }
	
	/**
	 * Set all state to uninitialized defaults.
	 */
	void reset() {
		state_ = ReportingState::NO_READ;
		paired_ = false;
		nconcord_ = 0;
		ndiscord_ = 0;
		nunpair1_ = 0;
		nunpair2_ = 0;
		doneConcord_ = false;
		doneDiscord_ = false;
		doneUnpair_  = false;
		doneUnpair1_ = false;
		doneUnpair2_ = false;
		exitConcord_ = ReportingState::EXIT_DID_NOT_ENTER;
		exitDiscord_ = ReportingState::EXIT_DID_NOT_ENTER;
		exitUnpair1_ = ReportingState::EXIT_DID_NOT_ENTER;
		exitUnpair2_ = ReportingState::EXIT_DID_NOT_ENTER;
		done_ = false;
	}
	
	/**
	 * Return true iff this ReportingState has been initialized with a call to
	 * nextRead() since the last time reset() was called.
	 */
	bool inited() const { return state_ != ReportingState::NO_READ; }

	/**
	 * Initialize state machine with a new read.  The state we start in depends
	 * on whether it's paired-end or unpaired.
	 */
	void nextRead(bool paired);

	/**
	 * Caller uses this member function to indicate that one additional
	 * concordant alignment has been found.
	 */
	bool foundConcordant();

	/**
	 * Caller uses this member function to indicate that one additional
	 * discordant alignment has been found.
	 */
	bool foundUnpaired(bool mate1);
	
	/**
	 * Called to indicate that the aligner has finished searching for
	 * alignments.  This gives us a chance to finalize our state.
	 *
	 * TODO: Keep track of short-circuiting information.
	 */
	void finish();
	
	/**
	 * Populate given counters with the number of various kinds of alignments
	 * to report for this read.  Concordant alignments are preferable to (and
	 * mutually exclusive with) discordant alignments, and paired-end
	 * alignments are preferable to unpaired alignments.
	 *
	 * The caller also needs some additional information for the case where a
	 * pair or unpaired read aligns repetitively.  If the read is paired-end
	 * and the paired-end has repetitive concordant alignments, that should be
	 * reported, and 'pairMax' is set to true to indicate this.  If the read is
	 * paired-end, does not have any conordant alignments, but does have
	 * repetitive alignments for one or both mates, then that should be
	 * reported, and 'unpair1Max' and 'unpair2Max' are set accordingly.
	 *
	 * Note that it's possible in the case of a paired-end read for the read to
	 * have repetitive concordant alignments, but for one mate to have a unique
	 * unpaired alignment.
	 */
	void getReport(
		uint64_t& nconcordAln, // # concordant alignments to report
		uint64_t& ndiscordAln, // # discordant alignments to report
		uint64_t& nunpair1Aln, // # unpaired alignments for mate #1 to report
		uint64_t& nunpair2Aln, // # unpaired alignments for mate #2 to report
		bool& pairMax,         // repetitive concordant alignments
		bool& unpair1Max,      // repetitive alignments for mate #1
		bool& unpair2Max)      // repetitive alignments for mate #2
		const;

	/**
	 * Return an integer representing the alignment state we're in.
	 */
	inline int state() const { return state_; }
	
	/**
	 * If false, there's no need to solve any more dynamic programming problems
	 * for finding opposite mates.
	 */
	inline bool doneConcordant() const { return doneConcord_; }
	
	/**
	 * If false, there's no need to seek any more discordant alignment.
	 */
	inline bool doneDiscordant() const { return doneDiscord_; }
	
	/**
	 * If false, there's no need to seek any more unpaired alignments for the
	 * specified mate.  Note: this doesn't necessarily mean we can stop looking
	 * for alignments for the mate, since this might be necessary for finding
	 * concordant and discordant alignments.
	 */
	inline bool doneUnpaired(bool mate1) const {
		return mate1 ? doneUnpair1_ : doneUnpair2_;
	}
	
	/**
	 * If false, no further consideration of the given mate is necessary.  It's
	 * not needed for *any* class of alignment: concordant, discordant or
	 * unpaired.
	 */
	inline bool doneWithMate(bool mate1) const {
		bool doneUnpair = mate1 ? doneUnpair1_ : doneUnpair2_;
		uint64_t nun = mate1 ? nunpair1_ : nunpair2_;
		if(!doneUnpair || !doneConcord_) {
			return false; // still needed for future concordant/unpaired alns
		}
		if(!doneDiscord_ && nun == 0) {
			return false; // still needed for future discordant alignments
		}
		return true; // done
	}

	/**
	 * Return true iff there's no need to seek any more unpaired alignments.
	 */
	inline bool doneUnpaired() const { return doneUnpair_; }
	
	/**
	 * Return true iff all alignment stages have been exited.
	 */
	inline bool done() const { return done_; }

	inline uint64_t numConcordant() const { return nconcord_; }
	inline uint64_t numDiscordant() const { return ndiscord_; }
	inline uint64_t numUnpaired1()  const { return nunpair1_; }
	inline uint64_t numUnpaired2()  const { return nunpair2_; }

	inline int exitConcordant() const { return exitConcord_; }
	inline int exitDiscordant() const { return exitDiscord_; }
	inline int exitUnpaired1()  const { return exitUnpair1_; }
	inline int exitUnpaired2()  const { return exitUnpair2_; }

#ifndef NDEBUG
	/**
	 * Check that ReportingState is internally consistent.
	 */
	bool repOk() const {
		assert(p_.discord || doneDiscord_);
		assert(p_.mixed   || !paired_ || doneUnpair_);
		assert(doneUnpair_ || !doneUnpair1_ || !doneUnpair2_);
		if(p_.mhitsSet()) {
			assert_leq(numConcordant(), (uint64_t)p_.mhits+1);
			assert_leq(numDiscordant(), (uint64_t)p_.mhits+1);
			assert(paired_ || numUnpaired1() <= (uint64_t)p_.mhits+1);
			assert(paired_ || numUnpaired2() <= (uint64_t)p_.mhits+1);
		}
		assert(done() || !doneWithMate(true) || !doneWithMate(false));
		return true;
	}
#endif

	/**
	 * Return ReportingParams object governing this ReportingState.
	 */
	const ReportingParams& params() const {
		return p_;
	}

protected:

	/**
	 * Update state to reflect situation after converting two unique unpaired
	 * alignments, one for mate 1 and one for mate 2, into a single discordant
	 * alignment.
	 */
	void convertUnpairedToDiscordant() {
		assert_eq(1, numUnpaired1());
		assert_eq(1, numUnpaired2());
		assert_eq(0, numDiscordant());
		exitUnpair1_ = exitUnpair2_ = ReportingState::EXIT_CONVERTED_TO_DISCORDANT;
		nunpair1_ = nunpair2_ = 0;
		ndiscord_ = 1;
		assert_eq(1, numDiscordant());
	}

	/**
	 * Given the number of alignments in a category, check whether we
	 * short-circuited out of the category.  Set the done and exit arguments to
	 * indicate whether and how we short-circuited.
	 */
	inline void areDone(
		uint64_t cnt,     // # alignments in category
		bool& done,       // out: whether we short-circuited out of category
		int& exit) const; // out: if done, how we short-circuited (-k? -m? etc)
	
	/**
	 * Update done_ field to reflect whether we're totally done now.
	 */
	inline void updateDone() {
		doneUnpair_ = doneUnpair1_ && doneUnpair2_;
		done_ = doneUnpair_ && doneDiscord_ && doneConcord_;
	}

	const ReportingParams& p_;  // reporting parameters
	int state_;          // state we're currently in
	bool paired_;        // true iff read we're currently handling is paired
	uint64_t nconcord_;  // # concordants found so far
	uint64_t ndiscord_;  // # discordants found so far
	uint64_t nunpair1_;  // # unpaired alignments found so far for mate 1
	uint64_t nunpair2_;  // # unpaired alignments found so far for mate 2
	bool doneConcord_;   // true iff we're no longner interested in concordants
	bool doneDiscord_;   // true iff we're no longner interested in discordants
	bool doneUnpair_;    // no longner interested in unpaired alns
	bool doneUnpair1_;   // no longner interested in unpaired alns for mate 1
	bool doneUnpair2_;   // no longner interested in unpaired alns for mate 2
	int exitConcord_;    // flag indicating how we exited concordant state
	int exitDiscord_;    // flag indicating how we exited discordant state
	int exitUnpair1_;    // flag indicating how we exited unpaired 1 state
	int exitUnpair2_;    // flag indicating how we exited unpaired 2 state
	bool done_;          // done with all alignments
};

/**
 * Global hit sink for hits from the MultiSeed aligner.  Encapsulates
 * all aspects of the MultiSeed aligner hitsink that are global to all
 * threads.  This includes aspects relating to:
 *
 * (a) synchronized access to the output stream
 * (b) the policy to be enforced by the per-thread wrapper
 *
 * TODO: Implement splitting up of alignments into separate files
 * according to genomic coordinate.
 */
class AlnSink {

	typedef EList<std::string> StrList;

public:

	explicit AlnSink(
		OutputQueue& oq,
		const StrList& refnames,
		bool quiet) :
		oq_(oq),
		refnames_(refnames),
		quiet_(quiet)
	{ }

	/**
	 * Destroy HitSinkobject;
	 */
	virtual ~AlnSink() { }

	/**
	 * Called when the AlnSink is wrapped by a new AlnSinkWrap.  This helps us
	 * keep track of whether the main lock or any of the per-stream locks will
	 * be contended by multiple threads.
	 */
	void addWrapper() { numWrappers_++; }

	/**
	 * Append a single hit to the given output stream.  If
	 * synchronization is required, append() assumes the caller has
	 * already grabbed the appropriate lock.
	 */
	virtual void append(
		BTString&             o,
		StackedAln&           staln,
		size_t                threadId,
		const Read           *rd1,
		const Read           *rd2,
		const TReadId         rdid,
		AlnRes               *rs1,
		AlnRes               *rs2,
		const AlnSetSumm&     summ,
		const SeedAlSumm&     ssm1,
		const SeedAlSumm&     ssm2,
		const AlnFlags*       flags1,
		const AlnFlags*       flags2,
		const PerReadMetrics& prm,
		const Mapq&           mapq,
		const Scoring&        sc,
		bool                  report2) = 0;

	/**
	 * Report a given batch of hits for the given read or read pair.
	 * Should be called just once per read pair.  Assumes all the
	 * alignments are paired, split between rs1 and rs2.
	 *
	 * The caller hasn't decided which alignments get reported as primary
	 * or secondary; that's up to the routine.  Because the caller might
	 * want to know this, we use the pri1 and pri2 out arguments to
	 * convey this.
	 */
	virtual void reportHits(
		BTString&             o,              // write to this buffer
		StackedAln&           staln,       // StackedAln to write stacked alignment
		size_t                threadId,       // which thread am I?
		const Read           *rd1,            // mate #1
		const Read           *rd2,            // mate #2
		const TReadId         rdid,           // read ID
		const EList<size_t>&  select1,        // random subset of rd1s
		const EList<size_t>*  select2,        // random subset of rd2s
		EList<AlnRes>        *rs1,            // alignments for mate #1
		EList<AlnRes>        *rs2,            // alignments for mate #2
		bool                  maxed,          // true iff -m/-M exceeded
		const AlnSetSumm&     summ,           // summary
		const SeedAlSumm&     ssm1,           // seed alignment summ
		const SeedAlSumm&     ssm2,           // seed alignment summ
		const AlnFlags*       flags1,         // flags for mate #1
		const AlnFlags*       flags2,         // flags for mate #2
		const PerReadMetrics& prm,            // per-read metrics
		const Mapq&           mapq,           // MAPQ generator
		const Scoring&        sc,             // scoring scheme
		bool                  getLock = true) // true iff lock held by caller
	{
		// There are a few scenarios:
		// 1. Read is unpaired, in which case rd2 is NULL
		// 2. Read is paired-end and we're reporting concordant alignments
		// 3. Read is paired-end and we're reporting discordant alignments
		// 4. Read is paired-end and we're reporting unpaired alignments for
		//    both mates
		// 5. Read is paired-end and we're reporting an unpaired alignments for
		//    just one mate or the other
		assert(rd1 != NULL || rd2 != NULL);
		assert(rs1 != NULL || rs2 != NULL);
		AlnFlags flagscp1, flagscp2;
		if(flags1 != NULL) {
			flagscp1 = *flags1;
			flags1 = &flagscp1;
			flagscp1.setPrimary(true);
		}
		if(flags2 != NULL) {
			flagscp2 = *flags2;
			flags2 = &flagscp2;
			flagscp2.setPrimary(true);
		}
		if(select2 != NULL) {
			// Handle case 5
			assert(rd1 != NULL); assert(flags1 != NULL);
			assert(rd2 != NULL); assert(flags2 != NULL);
			assert_gt(select1.size(), 0);
			assert_gt(select2->size(), 0);
			AlnRes* r1pri = ((rs1 != NULL) ? &rs1->get(select1[0]) : NULL);
			AlnRes* r2pri = ((rs2 != NULL) ? &rs2->get((*select2)[0]) : NULL);
			append(o, staln, threadId, rd1, rd2, rdid, r1pri, r2pri, summ,
			       ssm1, ssm2, flags1, flags2, prm, mapq, sc, true);
			flagscp1.setPrimary(false);
			flagscp2.setPrimary(false);
			for(size_t i = 1; i < select1.size(); i++) {
				AlnRes* r1 = ((rs1 != NULL) ? &rs1->get(select1[i]) : NULL);
				append(o, staln, threadId, rd1, rd2, rdid, r1, r2pri, summ,
				       ssm1, ssm2, flags1, flags2, prm, mapq, sc, false);
			}
			for(size_t i = 1; i < select2->size(); i++) {
				AlnRes* r2 = ((rs2 != NULL) ? &rs2->get((*select2)[i]) : NULL);
				append(o, staln, threadId, rd2, rd1, rdid, r2, r1pri, summ,
				       ssm2, ssm1, flags2, flags1, prm, mapq, sc, false);
			}
		} else {
			// Handle cases 1-4
			for(size_t i = 0; i < select1.size(); i++) {
				AlnRes* r1 = ((rs1 != NULL) ? &rs1->get(select1[i]) : NULL);
				AlnRes* r2 = ((rs2 != NULL) ? &rs2->get(select1[i]) : NULL);
				append(o, staln, threadId, rd1, rd2, rdid, r1, r2, summ,
				       ssm1, ssm2, flags1, flags2, prm, mapq, sc, true);
				if(flags1 != NULL) {
					flagscp1.setPrimary(false);
				}
				if(flags2 != NULL) {
					flagscp2.setPrimary(false);
				}
			}
		}
	}

	/**
	 * Report an unaligned read.  Typically we do nothing, but we might
	 * want to print a placeholder when output is chained.
	 */
	virtual void reportUnaligned(
		BTString&             o,              // write to this string
		StackedAln&           staln,          // StackedAln to write stacked alignment
		size_t                threadId,       // which thread am I?
		const Read           *rd1,            // mate #1
		const Read           *rd2,            // mate #2
		const TReadId         rdid,           // read ID
		const AlnSetSumm&     summ,           // summary
		const SeedAlSumm&     ssm1,           // seed alignment summary
		const SeedAlSumm&     ssm2,           // seed alignment summary
		const AlnFlags*       flags1,         // flags for mate #1
		const AlnFlags*       flags2,         // flags for mate #2
		const PerReadMetrics& prm,            // per-read metrics
		const Mapq&           mapq,           // MAPQ calculator
		const Scoring&        sc,             // scoring scheme
		bool                  report2,        // report alns for both mates?
		bool                  getLock = true) // true iff lock held by caller
	{
		append(o, staln, threadId, rd1, rd2, rdid, NULL, NULL, summ,
		       ssm1, ssm2, flags1, flags2, prm, mapq, sc, report2);
	}

	/**
	 * Print summary of how many reads aligned, failed to align and aligned
	 * repetitively.  Write it to stderr.  Optionally write Hadoop counter
	 * updates.
	 */
	void printAlSumm(
		const ReportingMetrics& met,
		size_t repThresh, // threshold for uniqueness, or max if no thresh
		bool discord,     // looked for discordant alignments
		bool mixed,       // looked for unpaired alignments where paired failed?
		bool hadoopOut);  // output Hadoop counters?

	/**
	 * Called when all alignments are complete.  It is assumed that no
	 * synchronization is necessary.
	 */
	void finish(
		size_t repThresh,
		bool discord,
		bool mixed,
		bool hadoopOut)
	{
		// Close output streams
		if(!quiet_) {
			printAlSumm(
				met_,
				repThresh,
				discord,
				mixed,
				hadoopOut);
		}
	}

#ifndef NDEBUG
	/**
	 * Check that hit sink is internally consistent.
	 */
	bool repOk() const { return true; }
#endif
	
	//
	// Related to reporting seed hits
	//

	/**
	 * Given a Read and associated, filled-in SeedResults objects,
	 * print a record summarizing the seed hits.
	 */
	void reportSeedSummary(
		BTString&          o,
		const Read&        rd,
		TReadId            rdid,
		size_t             threadId,
		const SeedResults& rs,
		bool               getLock = true);

	/**
	 * Given a Read, print an empty record (all 0s).
	 */
	void reportEmptySeedSummary(
		BTString&          o,
		const Read&        rd,
		TReadId            rdid,
		size_t             threadId,
		bool               getLock = true);

	/**
	 * Append a batch of unresolved seed alignment results (i.e. seed
	 * alignments where all we know is the reference sequence aligned
	 * to and its SA range, not where it falls in the reference
	 * sequence) to the given output stream in Bowtie's seed-alignment
	 * verbose-mode format.
	 */
	virtual void appendSeedSummary(
		BTString&     o,
		const Read&   rd,
		const TReadId rdid,
		size_t        seedsTried,
		size_t        nonzero,
		size_t        ranges,
		size_t        elts,
		size_t        seedsTriedFw,
		size_t        nonzeroFw,
		size_t        rangesFw,
		size_t        eltsFw,
		size_t        seedsTriedRc,
		size_t        nonzeroRc,
		size_t        rangesRc,
		size_t        eltsRc);

	/**
	 * Merge given metrics in with ours by summing all individual metrics.
	 */
	void mergeMetrics(const ReportingMetrics& met, bool getLock = true) {
		met_.merge(met, getLock);
	}

	/**
	 * Return mutable reference to the shared OutputQueue.
	 */
	OutputQueue& outq() {
		return oq_;
	}

protected:

	OutputQueue&       oq_;           // output queue
	int                numWrappers_;  // # threads owning a wrapper for this HitSink
	const StrList&     refnames_;     // reference names
	bool               quiet_;        // true -> don't print alignment stats at the end
	ReportingMetrics   met_;          // global repository of reporting metrics
};

/**
 * Per-thread hit sink "wrapper" for the MultiSeed aligner.  Encapsulates
 * aspects of the MultiSeed aligner hit sink that are per-thread.  This
 * includes aspects relating to:
 *
 * (a) Enforcement of the reporting policy
 * (b) Tallying of results
 * (c) Storing of results for the previous read in case this allows us to
 *     short-circuit some work for the next read (i.e. if it's identical)
 *
 * PHASED ALIGNMENT ASSUMPTION
 *
 * We make some assumptions about how alignment proceeds when we try to
 * short-circuit work for identical reads.  Specifically, we assume that for
 * each read the aligner proceeds in a series of stages (or perhaps just one
 * stage).  In each stage, the aligner either:
 *
 * (a)  Finds no alignments, or
 * (b)  Finds some alignments and short circuits out of the stage with some
 *      random reporting involved (e.g. in -k and/or -M modes), or
 * (c)  Finds all of the alignments in the stage
 *
 * In the event of (a), the aligner proceeds to the next stage and keeps
 * trying; we can skip the stage entirely for the next read if it's identical.
 * In the event of (b), or (c), the aligner stops and does not proceed to
 * further stages.  In the event of (b1), if the next read is identical we
 * would like to tell the aligner to start again at the beginning of the stage
 * that was short-circuited.
 *
 * In any event, the rs1_/rs2_/rs1u_/rs2u_ fields contain the alignments found
 * in the last alignment stage attempted.
 *
 * HANDLING REPORTING LIMITS
 *
 * The user can specify reporting limits, like -k (specifies number of
 * alignments to report out of those found) and -M (specifies a ceiling s.t. if
 * there are more alignments than the ceiling, read is called repetitive and
 * best found is reported).  Enforcing these limits is straightforward for
 * unpaired alignments: if a new alignment causes us to exceed the -M ceiling,
 * we can stop looking.
 *
 * The case where both paired-end and unpaired alignments are possible is
 * trickier.  Once we have a number of unpaired alignments that exceeds the
 * ceiling, we can stop looking *for unpaired alignments* - but we can't
 * necessarily stop looking for paired-end alignments, since there may yet be
 * more to find.  However, if the input read is not a pair, then we can stop at
 * this point.  If the input read is a pair and we have a number of paired
 * aligments that exceeds the -M ceiling, we can stop looking.
 *
 * CONCORDANT & DISCORDANT, PAIRED & UNPAIRED
 *
 * A note on paired-end alignment: Clearly, if an input read is
 * paired-end and we find either concordant or discordant paired-end
 * alignments for the read, then we would like to tally and report
 * those alignments as such (and not as groups of 2 unpaired
 * alignments).  And if we fail to find any paired-end alignments, but
 * we do find some unpaired alignments for one mate or the other, then
 * we should clearly tally and report those alignments as unpaired
 * alignments (if the user so desires).
 *
 * The situation is murkier when there are no paired-end alignments,
 * but there are unpaired alignments for *both* mates.  In this case,
 * we might want to pick out zero or more pairs of mates and classify
 * those pairs as discordant paired-end alignments.  And we might want
 * to classify the remaining alignments as unpaired.  But how do we
 * pick which pairs if any to call discordant?
 *
 * Because the most obvious use for discordant pairs is for identifying
 * large-scale variation, like rearrangements or large indels, we would
 * usually like to be conservative about what we call a discordant
 * alignment.  If there's a good chance that one or the other of the
 * two mates has a good alignment to another place on the genome, this
 * compromises the evidence for the large-scale variant.  For this
 * reason, Bowtie 2's policy is: if there are no paired-end alignments
 * and there is *exactly one alignment each* for both mates, then the
 * two alignments are paired and treated as a discordant paired-end
 * alignment.  Otherwise, all alignments are treated as unpaired
 * alignments.
 *
 * When both paired and unpaired alignments are discovered by the
 * aligner, only the paired alignments are reported by default.  This
 * is sensible considering relative likelihoods: if a good paired-end
 * alignment is found, it is much more likely that the placement of
 * the two mates implied by that paired alignment is correct than any
 * placement implied by an unpaired alignment.
 *
 * 
 */
class AlnSinkWrap {
public:

	AlnSinkWrap(
		AlnSink& g,                // AlnSink being wrapped
		const ReportingParams& rp, // Parameters governing reporting
		Mapq& mapq,                // Mapq calculator
		size_t threadId) :         // Thread ID
		g_(g),
		rp_(rp),
		threadid_(threadId),
		mapq_(mapq),
		init_(false),   
		maxed1_(false),       // read is pair and we maxed out mate 1 unp alns
		maxed2_(false),       // read is pair and we maxed out mate 2 unp alns
		maxedOverall_(false), // alignments found so far exceed -m/-M ceiling
		bestPair_(std::numeric_limits<TAlScore>::min()),
		best2Pair_(std::numeric_limits<TAlScore>::min()),
		bestUnp1_(std::numeric_limits<TAlScore>::min()),
		best2Unp1_(std::numeric_limits<TAlScore>::min()),
		bestUnp2_(std::numeric_limits<TAlScore>::min()),
		best2Unp2_(std::numeric_limits<TAlScore>::min()),
		rd1_(NULL),    // mate 1
		rd2_(NULL),    // mate 2
		rdid_(std::numeric_limits<TReadId>::max()), // read id
		rs1_(),        // mate 1 alignments for paired-end alignments
		rs2_(),        // mate 2 alignments for paired-end alignments
		rs1u_(),       // mate 1 unpaired alignments
		rs2u_(),       // mate 2 unpaired alignments
		select1_(),    // for selecting random subsets for mate 1
		select2_(),    // for selecting random subsets for mate 2
		st_(rp)        // reporting state - what's left to do?
	{
		assert(rp_.repOk());
	}

	/**
	 * Initialize the wrapper with a new read pair and return an
	 * integer >= -1 indicating which stage the aligner should start
	 * at.  If -1 is returned, the aligner can skip the read entirely.
	 * at.  If .  Checks if the new read pair is identical to the
	 * previous pair.  If it is, then we return the id of the first
	 * stage to run.
	 */
	int nextRead(
		// One of the other of rd1, rd2 will = NULL if read is unpaired
		const Read* rd1,      // new mate #1
		const Read* rd2,      // new mate #2
		TReadId rdid,         // read ID for new pair
		bool qualitiesMatter);// aln policy distinguishes b/t quals?

	/**
	 * Inform global, shared AlnSink object that we're finished with
	 * this read.  The global AlnSink is responsible for updating
	 * counters, creating the output record, and delivering the record
	 * to the appropriate output stream.
	 */
	void finishRead(
		const SeedResults *sr1,         // seed alignment results for mate 1
		const SeedResults *sr2,         // seed alignment results for mate 2
		bool               exhaust1,    // mate 1 exhausted?
		bool               exhaust2,    // mate 2 exhausted?
		bool               nfilt1,      // mate 1 N-filtered?
		bool               nfilt2,      // mate 2 N-filtered?
		bool               scfilt1,     // mate 1 score-filtered?
		bool               scfilt2,     // mate 2 score-filtered?
		bool               lenfilt1,    // mate 1 length-filtered?
		bool               lenfilt2,    // mate 2 length-filtered?
		bool               qcfilt1,     // mate 1 qc-filtered?
		bool               qcfilt2,     // mate 2 qc-filtered?
		bool               sortByScore, // prioritize alignments by score
		RandomSource&      rnd,         // pseudo-random generator
		ReportingMetrics&  met,         // reporting metrics
		const PerReadMetrics& prm,      // per-read metrics
		const Scoring& sc,              // scoring scheme
		bool suppressSeedSummary = true,
		bool suppressAlignments = false);
	
	/**
	 * Called by the aligner when a new unpaired or paired alignment is
	 * discovered in the given stage.  This function checks whether the
	 * addition of this alignment causes the reporting policy to be
	 * violated (by meeting or exceeding the limits set by -k, -m, -M),
	 * in which case true is returned immediately and the aligner is
	 * short circuited.  Otherwise, the alignment is tallied and false
	 * is returned.
	 */
	bool report(
		int stage,
		const AlnRes* rs1,
		const AlnRes* rs2);

#ifndef NDEBUG
	/**
	 * Check that hit sink wrapper is internally consistent.
	 */
	bool repOk() const {
		assert_eq(rs2_.size(), rs1_.size());
		if(rp_.mhitsSet()) {
			assert_gt(rp_.mhits, 0);
			assert_leq((int)rs1_.size(), rp_.mhits+1);
			assert_leq((int)rs2_.size(), rp_.mhits+1);
			assert(readIsPair() || (int)rs1u_.size() <= rp_.mhits+1);
			assert(readIsPair() || (int)rs2u_.size() <= rp_.mhits+1);
		}
		if(init_) {
			assert(rd1_ != NULL);
			assert_neq(std::numeric_limits<TReadId>::max(), rdid_);
		}
		assert_eq(st_.numConcordant() + st_.numDiscordant(), rs1_.size());
		//assert_eq(st_.numUnpaired1(), rs1u_.size());
		//assert_eq(st_.numUnpaired2(), rs2u_.size());
		assert(st_.repOk());
		return true;
	}
#endif
	
	/**
	 * Return true iff no alignments have been reported to this wrapper
	 * since the last call to nextRead().
	 */
	bool empty() const {
		return rs1_.empty() && rs1u_.empty() && rs2u_.empty();
	}
	
	/**
	 * Return true iff we have already encountered a number of alignments that
	 * exceeds the -m/-M ceiling.  TODO: how does this distinguish between
	 * pairs and mates?
	 */
	bool maxed() const {
		return maxedOverall_;
	}
	
	/**
	 * Return true if the current read is paired.
	 */
	bool readIsPair() const {
		return rd1_ != NULL && rd2_ != NULL;
	}
	
	/**
	 * Return true iff nextRead() has been called since the last time
	 * finishRead() was called.
	 */
	bool inited() const { return init_; }

	/**
	 * Return a const ref to the ReportingState object associated with the
	 * AlnSinkWrap.
	 */
	const ReportingState& state() const { return st_; }
	
	/**
	 * Return true iff we're in -M mode.
	 */
	bool Mmode() const {
		return rp_.mhitsSet();
	}
	
	/**
	 * Return true iff the policy is to report all hits.
	 */
	bool allHits() const {
		return rp_.allHits();
	}
	
	/**
	 * Return true iff at least two alignments have been reported so far for an
	 * unpaired read or mate 1.
	 */
	bool hasSecondBestUnp1() const {
		return best2Unp1_ != std::numeric_limits<TAlScore>::min();
	}

	/**
	 * Return true iff at least two alignments have been reported so far for
	 * mate 2.
	 */
	bool hasSecondBestUnp2() const {
		return best2Unp2_ != std::numeric_limits<TAlScore>::min();
	}

	/**
	 * Return true iff at least two paired-end alignments have been reported so
	 * far.
	 */
	bool hasSecondBestPair() const {
		return best2Pair_ != std::numeric_limits<TAlScore>::min();
	}
	
	/**
	 * Get best score observed so far for an unpaired read or mate 1.
	 */
	TAlScore bestUnp1() const {
		return bestUnp1_;
	}

	/**
	 * Get second-best score observed so far for an unpaired read or mate 1.
	 */
	TAlScore secondBestUnp1() const {
		return best2Unp1_;
	}

	/**
	 * Get best score observed so far for mate 2.
	 */
	TAlScore bestUnp2() const {
		return bestUnp2_;
	}

	/**
	 * Get second-best score observed so far for mate 2.
	 */
	TAlScore secondBestUnp2() const {
		return best2Unp2_;
	}

	/**
	 * Get best score observed so far for paired-end read.
	 */
	TAlScore bestPair() const {
		return bestPair_;
	}

	/**
	 * Get second-best score observed so far for paired-end read.
	 */
	TAlScore secondBestPair() const {
		return best2Pair_;
	}

protected:

	/**
	 * Return true iff the read in rd1/rd2 matches the last read handled, which
	 * should still be in rd1_/rd2_.
	 */
	bool sameRead(
		const Read* rd1,
		const Read* rd2,
		bool qualitiesMatter);

	/**
	 * If there is a configuration of unpaired alignments that fits our
	 * criteria for there being one or more discordant alignments, then
	 * shift the discordant alignments over to the rs1_/rs2_ lists, clear the
	 * rs1u_/rs2u_ lists and return true.  Otherwise, return false.
	 */
	bool prepareDiscordants();

	/**
	 * Given that rs is already populated with alignments, consider the
	 * alignment policy and make random selections where necessary.  E.g. if we
	 * found 10 alignments and the policy is -k 2 -m 20, select 2 alignments at
	 * random.  We "select" an alignment by setting the parallel entry in the
	 * 'select' list to true.
	 */
	size_t selectAlnsToReport(
		const EList<AlnRes>& rs,     // alignments to select from
		uint64_t             num,    // number of alignments to select
		EList<size_t>&       select, // list to put results in
		RandomSource&        rnd)
		const;

	/**
	 * rs1 (possibly together with rs2 if reads are paired) are populated with
	 * alignments.  Here we prioritize them according to alignment score, and
	 * some randomness to break ties.  Priorities are returned in the 'select'
	 * list.
	 */
	size_t selectByScore(
		const EList<AlnRes>* rs1,    // alignments to select from (mate 1)
		const EList<AlnRes>* rs2,    // alignments to select from (mate 2, or NULL)
		uint64_t             num,    // number of alignments to select
		EList<size_t>&       select, // prioritized list to put results in
		RandomSource&        rnd)
		const;

	AlnSink&        g_;     // global alignment sink
	ReportingParams rp_;    // reporting parameters: khits, mhits etc
	size_t          threadid_; // thread ID
	Mapq&           mapq_;  // mapq calculator
	bool            init_;  // whether we're initialized w/ read pair
	bool            maxed1_; // true iff # unpaired mate-1 alns reported so far exceeded -m/-M
	bool            maxed2_; // true iff # unpaired mate-2 alns reported so far exceeded -m/-M
	bool            maxedOverall_; // true iff # paired-end alns reported so far exceeded -m/-M
	TAlScore        bestPair_;     // greatest score so far for paired-end
	TAlScore        best2Pair_;    // second-greatest score so far for paired-end
	TAlScore        bestUnp1_;     // greatest score so far for unpaired/mate1
	TAlScore        best2Unp1_;    // second-greatest score so far for unpaired/mate1
	TAlScore        bestUnp2_;     // greatest score so far for mate 2
	TAlScore        best2Unp2_;    // second-greatest score so far for mate 2
	const Read*     rd1_;   // mate #1
	const Read*     rd2_;   // mate #2
	TReadId         rdid_;  // read ID (potentially used for ordering)
	EList<AlnRes>   rs1_;   // paired alignments for mate #1
	EList<AlnRes>   rs2_;   // paired alignments for mate #2
	EList<AlnRes>   rs1u_;  // unpaired alignments for mate #1
	EList<AlnRes>   rs2u_;  // unpaired alignments for mate #2
	EList<size_t>   select1_; // parallel to rs1_/rs2_ - which to report
	EList<size_t>   select2_; // parallel to rs1_/rs2_ - which to report
	ReportingState  st_;      // reporting state - what's left to do?
	
	EList<std::pair<TAlScore, size_t> > selectBuf_;
	BTString obuf_;
	StackedAln staln_;
};

/**
 * An AlnSink concrete subclass for printing SAM alignments.  The user might
 * want to customize SAM output in various ways.  We encapsulate all these
 * customizations, and some of the key printing routines, in the SamConfig
 * class in sam.h/sam.cpp.
 */
class AlnSinkSam : public AlnSink {

	typedef EList<std::string> StrList;

public:

	AlnSinkSam(
		OutputQueue&     oq,           // output queue
		const SamConfig& samc,         // settings & routines for SAM output
		const StrList&   refnames,     // reference names
		bool             quiet) :      // don't print alignment summary at end
		AlnSink(
			oq,
			refnames,
			quiet),
		samc_(samc)
	{ }
	
	virtual ~AlnSinkSam() { }

	/**
	 * Append a single alignment result, which might be paired or
	 * unpaired, to the given output stream in Bowtie's verbose-mode
	 * format.  If the alignment is paired-end, print mate1's alignment
	 * then mate2's alignment.
	 */
	virtual void append(
		BTString&     o,           // write output to this string
		StackedAln&   staln,       // StackedAln to write stacked alignment
		size_t        threadId,    // which thread am I?
		const Read*   rd1,         // mate #1
		const Read*   rd2,         // mate #2
		const TReadId rdid,        // read ID
		AlnRes* rs1,               // alignments for mate #1
		AlnRes* rs2,               // alignments for mate #2
		const AlnSetSumm& summ,    // summary
		const SeedAlSumm& ssm1,    // seed alignment summary
		const SeedAlSumm& ssm2,    // seed alignment summary
		const AlnFlags* flags1,    // flags for mate #1
		const AlnFlags* flags2,    // flags for mate #2
		const PerReadMetrics& prm, // per-read metrics
		const Mapq& mapq,          // MAPQ calculator
		const Scoring& sc,         // scoring scheme
		bool report2)              // report alns for both mates
	{
		assert(rd1 != NULL || rd2 != NULL);
		if(rd1 != NULL) {
			assert(flags1 != NULL);
			appendMate(o, staln, *rd1, rd2, rdid, rs1, rs2, summ, ssm1, ssm2,
			           *flags1, prm, mapq, sc);
		}
		if(rd2 != NULL && report2) {
			assert(flags2 != NULL);
			appendMate(o, staln, *rd2, rd1, rdid, rs2, rs1, summ, ssm2, ssm1,
			           *flags2, prm, mapq, sc);
		}
	}

protected:

	/**
	 * Append a single per-mate alignment result to the given output
	 * stream.  If the alignment is part of a pair, information about
	 * the opposite mate and its alignment are given in rdo/rso.
	 */
	void appendMate(
		BTString&     o,
		StackedAln&   staln,
		const Read&   rd,
		const Read*   rdo,
		const TReadId rdid,
		AlnRes* rs,
		AlnRes* rso,
		const AlnSetSumm& summ,
		const SeedAlSumm& ssm,
		const SeedAlSumm& ssmo,
		const AlnFlags& flags,
		const PerReadMetrics& prm, // per-read metrics
		const Mapq& mapq,          // MAPQ calculator
		const Scoring& sc);        // scoring scheme

	const SamConfig& samc_;    // settings & routines for SAM output
	BTDnaString      dseq_;    // buffer for decoded read sequence
	BTString         dqual_;   // buffer for decoded quality sequence
};

#endif /*ndef ALN_SINK_H_*/
