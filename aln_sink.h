//
//  aln_sink.h
//

#ifndef ALN_SINK_H_
#define ALN_SINK_H_

#include "filebuf.h"
#include "read.h"
#include "unique.h"
#include "refmap.h"
#include "sam.h"
#include "ds.h"

// Forward decl	
class SeedResults;

/**
 * Metrics summarizing the work done by the reporter and summarizing
 * the number of reads that align, that fail to align, and that align
 * non-uniquely.
 */
struct ReportingMetrics {

	ReportingMetrics() { reset(); MUTEX_INIT(lock); }
	
	void reset() {
		init(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	}
	
	void init(
		uint64_t nread_,
		uint64_t npaired_,
		uint64_t nunpaired_,
		uint64_t nconcord_uni_,
		uint64_t nconcord_rep_,
		uint64_t nconcord_0_,
		uint64_t ndiscord_,
		uint64_t nunp_0_uni_,
		uint64_t nunp_0_rep_,
		uint64_t nunp_0_0_,
		uint64_t nunp_rep_uni_,
		uint64_t nunp_rep_rep_,
		uint64_t nunp_rep_0_,
		uint64_t nunp_uni_,
		uint64_t nunp_rep_,
		uint64_t nunp_0_)
	{
		nread        = nread_;
		
		npaired      = npaired_;
		nunpaired    = nunpaired_;
		
		nconcord_uni = nconcord_uni_;
		nconcord_rep = nconcord_rep_;
		nconcord_0   = nconcord_0_;
		
		ndiscord     = ndiscord_;
		
		nunp_0_uni   = nunp_0_uni_;
		nunp_0_rep   = nunp_0_rep_;
		nunp_0_0     = nunp_0_0_;

		nunp_rep_uni = nunp_rep_uni_;
		nunp_rep_rep = nunp_rep_rep_;
		nunp_rep_0   = nunp_rep_0_;

		nunp_uni     = nunp_uni_;
		nunp_rep     = nunp_rep_;
		nunp_0       = nunp_0_;
	}
	
	/**
	 * Merge (add) the counters in the given ReportingMetrics object
	 * into this object.  This is the only safe way to update a
	 * ReportingMetrics shared by multiple threads.
	 */
	void merge(const ReportingMetrics& met, bool getLock = false) {
		ThreadSafe ts(&lock, getLock);
		
		nread        += met.nread;
		
		npaired      += met.npaired;
		nunpaired    += met.nunpaired;
		
		nconcord_uni += met.nconcord_uni;
		nconcord_rep += met.nconcord_rep;
		nconcord_0   += met.nconcord_0;
		
		ndiscord     += met.ndiscord;
		
		nunp_0_uni   += met.nunp_0_uni;
		nunp_0_rep   += met.nunp_0_rep;
		nunp_0_0     += met.nunp_0_0;

		nunp_rep_uni += met.nunp_rep_uni;
		nunp_rep_rep += met.nunp_rep_rep;
		nunp_rep_0   += met.nunp_rep_0;

		nunp_uni     += met.nunp_uni;
		nunp_rep     += met.nunp_rep;
		nunp_0       += met.nunp_0;
	}

	uint64_t  nread;        // # reads
	uint64_t  npaired;      // # pairs
	uint64_t  nunpaired;    // # unpaired reads
	
	// Paired
	
	//  Concordant
	uint64_t  nconcord_uni; // # pairs with unique concordant alns
	uint64_t  nconcord_rep; // # pairs with repetitive concordant alns
	uint64_t  nconcord_0;   // # pairs with 0 concordant alns
	//  Discordant
	uint64_t  ndiscord;     // # pairs with 1 discordant aln
	
	//  Unpaired from failed pairs
	uint64_t  nunp_0_uni;   // # unique from nconcord_0_ - ndiscord_
	uint64_t  nunp_0_rep;   // # repetitive from 
	uint64_t  nunp_0_0;     // # with 0 alignments

	//  Unpaired from repetitive pairs
	uint64_t  nunp_rep_uni; // # pairs with unique concordant alns
	uint64_t  nunp_rep_rep; // # pairs with repetitive concordant alns
	uint64_t  nunp_rep_0;   // # pairs with 0 concordant alns
	
	// Unpaired
	
	uint64_t  nunp_uni;   // # unique from nconcord_0_ - ndiscord_
	uint64_t  nunp_rep;   // # repetitive from 
	uint64_t  nunp_0;     // # with 0 alignments

	MUTEX_T lock;
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
	
	/**
	 * Check that reporting parameters are internally consistent.
	 */
	bool repOk() const {
		assert_geq(khits, 1);
		assert_geq(mhits, 1);
		return true;
	}
	
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
	void boostThresholds(
		float& posmin,
		float& posfrac,
		float& rowmult)
	{
		THitInt mul = mult();
		assert_gt(mul, 0);
		if(mul == std::numeric_limits<THitInt>::max()) {
			// If -a was specified, boost ROWM and POSF so that all hits are
			// tried for al; positions
			posmin =
			posfrac = 
			rowmult = std::numeric_limits<float>::max();
		} else if(mul > 1) {
			// If -k or -M were specified, boost ROWM and POSF so that an
			// appropriately larger number of hits are tried for more
			// positions
			posmin  *= mul;
			posfrac *= mul;
			rowmult *= mul;
		}
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
		EXIT_SHORT_CIRCUIT_m,         // -m exceeded
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
	bool foundDiscordant();

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
		OutFileBuf*        out,
		const Mapq&        mapq,       // mapping quality calculator
		bool               deleteOuts,
		const StrList&     refnames,
		bool               quiet) :
		outs_(),
		outNames_(),
		locks_(),
		deleteOuts_(deleteOuts),
		refnames_(refnames),
		quiet_(quiet),
		mapq_(mapq)
	{
		// Add the default output stream
		outs_.push_back(out);
		// Add its lock
		locks_.resize(1);
		// Initialize locks
		MUTEX_INIT(locks_[0]);
		MUTEX_INIT(mainlock_);
	}

	/**
	 * Destroy HitSinkobject;
	 */
	virtual ~AlnSink() { closeOuts(deleteOuts_); }

	/**
	 * Called when the AlnSink is wrapped by a new AlnSinkWrap.  This helps us
	 * keep track of whether the main lock or any of the per-stream locks will
	 * be contended by multiple threads.
	 */
	void addWrapper() { numWrappers_++; }

	/**
	 * Maps a read id and a reference coordinate (usually corresponding
	 * to the leftmost position on the Watson strand involved in the
	 * alignment) to a stream id used to determine which output stream
	 * to write results to.  If the user has requested that the output
	 * alignments appear in an order corresponding to the input order,
	 * it may be useful to partition by read id for a future sort step.
	 * If the user has requested that alignments appear binned and/or
	 * sorted by chromosome, then it may be useful to partition by
	 * reference coordinate.
	 */
	size_t streamId(TReadId rdid, Coord c) { return 0; }

	/**
	 * Append a single hit to the given output stream.  If
	 * synchronization is required, append() assumes the caller has
	 * already grabbed the appropriate lock.
	 */
	virtual void append(
		OutFileBuf&        o,
		const Read        *rd1,
		const Read        *rd2,
		const TReadId      rdid,
		const AlnRes      *rs1,
		const AlnRes      *rs2,
		const AlnSetSumm&  summ,
		const AlnFlags*    flags1,
		const AlnFlags*    flags2) = 0;

	/**
	 * Report a given batch of hits for the given read pair.  Should be
	 * called just once per read pair.
	 */
	virtual void reportHits(
		const Read          *rd1,            // mate #1
		const Read          *rd2,            // mate #2
		const TReadId        rdid,           // read ID
		const EList<size_t>& select,         // random subset
		const EList<AlnRes> *rs1,            // alignments for mate #1
		const EList<AlnRes> *rs2,            // alignments for mate #2
		bool                 maxed,          // true iff -m/-M exceeded
		const AlnSetSumm&    summ,           // summary
		const AlnFlags*      flags1,         // flags for mate #1
		const AlnFlags*      flags2,         // flags for mate #2
		bool                 getLock = true) // true iff lock held by caller
	{
		assert(rd1 != NULL || rd2 != NULL);
		assert(rs1 != NULL || rs2 != NULL);
		size_t sz = ((rs1 != NULL) ? rs1->size() : rs2->size());
		// Report all hits in the rs1/rs2 lists
		reportHits(
			rd1,
			rd2,
			rdid,
			select,
			rs1,
			rs2,
			maxed,
			0,
			sz,
			summ,
			flags1,
			flags2,
			getLock);
	}

	/**
	 * Report a given batch of hits for the given read pair.  Should be
	 * called just once per read pair.
	 */
	virtual void reportHits(
		const Read          *rd1,            // mate #1
		const Read          *rd2,            // mate #2
		const TReadId        rdid,           // read ID
		const EList<size_t>& select,         // random subset
		const EList<AlnRes> *rs1,            // alignments for mate #1
		const EList<AlnRes> *rs2,            // alignments for mate #2
		bool                 maxed,          // true iff -m/-M exceeded
		size_t               start,          // alignments to report: start
		size_t               end,            // alignments to report: end 
		const AlnSetSumm&    summ,           // summary
		const AlnFlags*      flags1,         // flags for mate #1
		const AlnFlags*      flags2,         // flags for mate #2
		bool                 getLock = true) // true iff lock held by caller
	{
		assert(rd1 != NULL || rd2 != NULL);
		assert(rs1 != NULL || rs2 != NULL);
		assert_geq(end, start);
		ASSERT_ONLY(size_t sz = ((rs1 != NULL) ? rs1->size() : rs2->size()));
		assert_leq(end, sz);
		size_t num = end - start;
		if(num == 0) {
			// Nothing to report
			return;
		}
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
		size_t sel_start = start;
		bool found = false;
		for(size_t i = start; i < end; i++) {
			if(select[i] == 1) {
				sel_start = i;
				found = true;
				break;
			}
		}
		assert(found);
		size_t i = sel_start;
		do {
			// Determine the stream id using the coordinate of the
			// upstream mate
			Coord c = ((rs1 != NULL) ?
				rs1->get(i).refcoord() :
				rs2->get(i).refcoord());
			const AlnRes* r1 = ((rs1 != NULL) ? &rs1->get(i) : NULL);
			const AlnRes* r2 = ((rs2 != NULL) ? &rs2->get(i) : NULL);
			size_t sid = streamId(rdid, c);
			assert_lt(sid, locks_.size());
			{
				ThreadSafe ts(&locks_[sid], getLock);
				append(out(sid), rd1, rd2, rdid, r1, r2, summ, flags1, flags2);
			}
			if(flags1 != NULL) {
				flagscp1.setPrimary(false);
			}
			if(flags2 != NULL) {
				flagscp2.setPrimary(false);
			}
			if(++i == end) {
				i = 0;
			}
		} while(select[i] > 0 && i != sel_start);
	}

	/**
	 * Report a read that aligned more times than allowed by the -m or
	 * -M ceiling.
	 */
	virtual void reportMaxed(
		const Read          *rd1,            // mate #1
		const Read          *rd2,            // mate #2
		const TReadId        rdid,           // read ID
		const EList<AlnRes> *rs1,            // alignments for mate #1
		const EList<AlnRes> *rs2,            // alignments for mate #2
		const AlnSetSumm&    summ,           // summary
		const AlnFlags*      flags1,         // flags for mate #1
		const AlnFlags*      flags2,         // flags for mate #2
		bool                 getLock = true) // true iff lock held by caller
	{
		assert(rd1 != NULL || rd2 != NULL);
		assert(rs1 != NULL || rs2 != NULL);
		size_t sz = ((rs1 != NULL) ? rs1->size() : rs2->size());
		reportMaxed(
			rd1,
			rd2,
			rdid,
			rs1,
			rs2,
			0,
			sz,
			summ,
			flags1,
			flags2,
			getLock);
	}

	/**
	 * Report a read that aligned more times than allowed by the -m or
	 * -M ceiling.
	 */
	virtual void reportMaxed(
		const Read          *rd1,            // mate #1
		const Read          *rd2,            // mate #2
		const TReadId        rdid,           // read ID
		const EList<AlnRes> *rs1,            // alignments for mate #1
		const EList<AlnRes> *rs2,            // alignments for mate #2
		size_t               start,          // alignments to report: start
		size_t               end,            // alignments to report: end 
		const AlnSetSumm&    summ,           // summary
		const AlnFlags*      flags1,         // flags for mate #1
		const AlnFlags*      flags2,         // flags for mate #2
		bool                 getLock = true) // true iff lock held by caller
	{
		assert(rd1 != NULL || rd2 != NULL);
		assert(rs1 != NULL || rs2 != NULL);
		{
			// Determine the stream id using the coordinate of the
			// upstream mate
			Coord c(0, 0, true);
			size_t sid = streamId(rdid, c);
			assert_lt(sid, locks_.size());
			ThreadSafe ts(&locks_[sid], getLock);
			append(out(sid), rd1, rd2, rdid, NULL, NULL, summ, flags1, flags2);
		}
	}

	/**
	 * Report an unaligned read.  Typically we do nothing, but we might
	 * want to print a placeholder when output is chained.
	 */
	virtual void reportUnaligned(
		const Read          *rd1,            // mate #1
		const Read          *rd2,            // mate #2
		const TReadId        rdid,           // read ID
		const AlnSetSumm&    summ,           // summary
		const AlnFlags*      flags1,         // flags for mate #1
		const AlnFlags*      flags2,         // flags for mate #2
		bool                 getLock = true) // true iff lock held by caller
	{
		{
			// Determine the stream id using the coordinate of the
			// upstream mate
			Coord c(0, 0, true);
			size_t sid = streamId(rdid, c);
			assert_lt(sid, locks_.size());
			ThreadSafe ts(&locks_[sid], getLock);
			append(out(sid), rd1, rd2, rdid, NULL, NULL, summ, flags1, flags2);
		}
	}

	/**
	 * Commit a reported hit.
	 */
	virtual void commitHits(
		const Read*          rd1,
		const Read*          rd2,
		const TReadId        rdid,
		const EList<AlnRes>* rs1,
		const EList<AlnRes>* rs2,
		size_t               start,
		size_t               end,
		bool                 getLock = true)
	{
		
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
		closeOuts(false);
		if(!quiet_) {
			printAlSumm(
				met_,
				repThresh,
				discord,
				mixed,
				hadoopOut);
		}
	}

	/**
	 * Returns the output stream associated with the given stream id.
	 * It lazily initializes the output stream first if necessary.
	 */
	OutFileBuf& out(size_t sid) {
		assert_lt(sid, outs_.size());
		if(outs_[sid] == NULL) {
			outs_[sid] = new OutFileBuf(outNames_[sid]);
		}
		assert(outs_[sid] != NULL);
		return *(outs_[sid]);
	}
	
	/**
	 * Check that hit sink is internally consistent.
	 */
	bool repOk() const {
		return true;
	}

	void dumpMaxed(const Read* m1, const Read* m2, TReadId rdid) {
	}
	
	void dumpUnal(const Read* m1, const Read* m2, TReadId rdid) {
	}
	
	void dumpAlign(const Read* m1, const Read* m2, TReadId rdid) {
	}
	
	//
	// Related to reporting seed hits
	//

	/**
	 * Given a Read and associated, filled-in SeedResults objects,
	 * print a record summarizing the seed hits.
	 */
	void reportSeedSummary(
		const Read&        rd,
		TReadId            rdid,
		const SeedResults& rs,
		bool               getLock = true);

	/**
	 * Given a Read, print an empty record (all 0s).
	 */
	void reportEmptySeedSummary(
		const Read&        rd,
		TReadId            rdid,
		bool               getLock = true);

	/**
	 * Append a batch of unresolved seed alignment results (i.e. seed
	 * alignments where all we know is the reference sequence aligned
	 * to and its SA range, not where it falls in the reference
	 * sequence) to the given output stream in Bowtie's seed-alignment
	 * verbose-mode format.
	 */
	virtual void appendSeedSummary(
		OutFileBuf&   o,
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

protected:

	/**
	 * Close (and flush) all OutFileBufs.
	 */
	void closeOuts(bool del) {
		// Flush and close all non-NULL output streams
		for(size_t i = 0; i < outs_.size(); i++) {
			if(outs_[i] != NULL && !outs_[i]->closed()) {
				outs_[i]->close();
			}
			if(del && outs_[i] != NULL) {
				delete outs_[i];
				outs_[i] = NULL;
			}
		}
	}

	// TODO: allow multiple output streams.  Right now, all output goes
	// to outs_[0].
	
	EList<OutFileBuf*> outs_;         // the alignment output stream(s)
	EList<std::string> outNames_;     // the filenames for the alignment output stream(s)
	EList<MUTEX_T>     locks_;        // pthreads mutexes for per-file critical sections
	bool               deleteOuts_;   // Whether to delete elements of outs_ upon exit
	int                numWrappers_;  // # threads owning a wrapper for this HitSink
	MUTEX_T            mainlock_;     // pthreads mutexes for fields of this object
	const StrList&     refnames_; // reference names
	bool               quiet_;        // true -> don't print alignment stats at the end
	const Mapq&        mapq_;         // mapping quality calculator
	
	ReportingMetrics   met_;          // global repository of reporting metrics
};

/**
 * Per-thread hit sink "wrapper" for the MultiSeed aligner.
 * Encapsulates aspects of the MultiSeed aligner hit sink that are
 * particular to a single thread.  This includes aspects relating to:
 *
 * (a) Enforcement of the global policy
 * (b) Tallying of results
 * (c) Storing of results for the previous read in case this allows us
 *     to short-circuit some work for the next read (i.e. if it's
 *     identical)
 *
 * RANDOM ORDER ASSUMPTION
 *
 * In order to short-circuit the alignment process when the -l limit is
 * reached (and there is no -m or -M limit), we must know that the
 * order in which the alignments are being found and reported is
 * reasonably "random".  If the order is not sufficiently random, then
 * we really ought to expend the additional effort needed to capture a
 * large (maybe comprehensive) sample of alignments and then pick
 * randomly from among those.
 *
 * Note that a reasonable definition of random *might* include a
 * preference for better-scoring alignments.  Thus, if the aligner
 * chooses seeds to extend in some order s.t. better-aligning seeds are
 * chosen before worse-aligning seeds, that may be OK.
 *
 * PHASED ALIGNMENT ASSUMPTION
 *
 * We make some assumptions about how alignment proceeds when we try to
 * short-circuit work for identical reads.  Specifically, we assume
 * that for each read the aligner proceeds in a series of stages (or
 * perhaps just one stage).  In each stage, the aligner either:
 *
 * (a)  Finds no alignments, or
 * (b1) Finds some alignments and short circuits out of the stage with
 *      some random reporting involved (e.g. in -k and/or -M modes), or
 * (b2) Finds some alignments and short circuits out of the stage
 *      without any random reporting involved (e.g. in -m mode), or
 * (c)  Finds all of the alignments in the stage
 *
 * In the event of (a), the aligner proceeds to the next stage and
 * keeps trying; we can skip the stage entirely for the next read if
 * it's identical.  In the event of (b1), (b2), or (c), the aligner
 * stops and does not proceed to further stages.  In the event of (b1),
 * if the next read is identical we would like to tell the aligner to
 * start again at the beginning of the stage that was short-circuited.
 * In the event of (b2), if the next read is identical we can skip the
 * read entirely.
 *
 * In any event, the rs1_ and rs2_ fields contain the alignments found
 * in the last alignment stage attempted.
 *
 * HANDLING REPORTING LIMITS
 *
 * The user can specify reporting limits, like -k (specifies number of
 * alignments to report out of pool of those found) and -m (specifies a
 * ceiling s.t. if there are more alignments than the ceiling, read is
 * called repetitive).  Enforcing these limits is straightforward for
 * unpaired alignments: if a new alignment causes us to exceed the -m
 * ceiling, we can stop looking.
 *
 * The case where both paired-end and unpaired alignments are possible
 * is trickier.  Once we have a number of unpaired alignments that
 * exceeds the ceiling, we can stop looking *for unpaired alignments* -
 * but we can't necessarily stop looking for paired-end alignments,
 * since there may yet be more to find.  However, if the input read is
 * not a pair, then we can stop at this point.  If the input read is a
 * pair and we have a number of paired aligments that exceeds the -m
 * ceiling, we can stop looking.
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
		AlnSink& g,                  // AlnSink being wrapped
		const ReportingParams& rp) : // Parameters governing reporting
		g_(g),
		rp_(rp),
		init_(false),   
		maxed1_(false),       // read is pair and we maxed out mate 1 unp alns
		maxed2_(false),       // read is pair and we maxed out mate 2 unp alns
		maxedOverall_(false), // alignments found so far exceed -m/-M ceiling
		best_(std::numeric_limits<THitInt>::max()),
		rd1_(NULL),    // mate 1
		rd2_(NULL),    // mate 2
		rd1buf_(),     // copy of mate 1 Read object
		rd2buf_(),     // copy of mate 2 Read object
		rdid_(std::numeric_limits<TReadId>::max()), // read id
		rs1_(),        // mate 1 alignments for paired-end alignments
		rs2_(),        // mate 2 alignments for paired-end alignments
		rs1u_(),       // mate 1 unpaired alignments
		rs2u_(),       // mate 2 unpaired alignments
		select_(),     // for selecting random subsets
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
		const SeedResults *sr1,
		const SeedResults *sr2,
		bool               exhaust1,
		bool               exhaust2,
		bool               nfilt1,
		bool               nfilt2,
		bool               scfilt1,
		bool               scfilt2,
		bool               lenfilt1,
		bool               lenfilt2,
		bool               qcfilt1,
		bool               qcfilt2,
		RandomSource&      rnd,
		ReportingMetrics&  met,
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

	AlnSink&        g_;     // global alignment sink
	ReportingParams rp_;    // reporting parameters: khits, mhits etc
	bool            init_;  // whether we're initialized w/ read pair
	bool            maxed1_; // true iff # unpaired mate-1 alns reported so far exceeded -m/-M
	bool            maxed2_; // true iff # unpaired mate-2 alns reported so far exceeded -m/-M
	bool            maxedOverall_; // true iff # paired-end alns reported so far exceeded -m/-M
	THitInt         best_;  // greatest score so far
	const Read*     rd1_;   // mate #1
	const Read*     rd2_;   // mate #2
	Read            rd1buf_;// buffer for mate #1
	Read            rd2buf_;// buffer for mate #2
	TReadId         rdid_;  // read ID (potentially used for ordering)
	EList<AlnRes>   rs1_;   // paired alignments for mate #1
	EList<AlnRes>   rs2_;   // paired alignments for mate #2
	EList<AlnRes>   rs1u_;  // unpaired alignments for mate #1
	EList<AlnRes>   rs2u_;  // unpaired alignments for mate #2
	EList<size_t>   select_;    // parallel to rs1_/rs2_ - which to report
	ReportingState  st_;    // reporting state - what's left to do?
};

/**
 * An AlnSink concrete subclass for printing Bowtie verbose-mode alignments.
 */
class AlnSinkVerbose : public AlnSink {

	typedef EList<std::string> StrList;

public:

	AlnSinkVerbose(
		OutFileBuf*        out,        // initial output stream
		const EList<bool>& suppress,   // suppress columns
		const Mapq&        mapq,       // mapping quality calculator
		bool               deleteOuts, // delete output objects upon destruction
		const StrList&     refnames,   // reference names
		bool               quiet,      // don't print alignment summary at end
		int                offBase,    // add to 0-based offsets before printing
		bool               colorSeq,   // color: print color seq, not decoded nucs
		bool               colorQual,  // color: print color quals, not decoded quals
		bool               exEnds,     // exclude ends for decoded colors alns
		bool               printPlaceholders, // print maxs and unals
		bool               printFlags, // print alignment flags a la SAM
		bool               printCost,  // print penalty in extra column
		bool               printParams,// print alignment parameters
		bool               printSeed,  // print pseudo-random seed
		ReferenceMap*      rmap,       // reference coordinate transformation
		bool               fullRef,    // print entire ref name incl whitespace
		int                partition = 0) : // partition size
		AlnSink(
			out,
			mapq,
			deleteOuts,
			refnames,
			quiet),
		suppress_(suppress),
		offBase_(offBase),
		colorSeq_(colorSeq),
		colorQual_(colorQual),
		exEnds_(exEnds),
		printPlaceholders_(printPlaceholders),
		printFlags_(printFlags),
		printCost_(printCost),
		printParams_(printParams),
		printSeed_(printSeed),
		rmap_(rmap),
		fullRef_(fullRef),
		partition_(partition)
	{ }

	/**
	 * Append a single alignment result, which might be paired or
	 * unpaired, to the given output stream in Bowtie's verbose-mode
	 * format.  If the alignment is paired-end, print mate1's alignment
	 * then mate2's alignment.
	 */
	virtual void append(
		OutFileBuf&   o,        // file buffer to write to
		const Read*   rd1,      // mate #1
		const Read*   rd2,      // mate #2
		const TReadId rdid,     // read ID
		const AlnRes* rs1,      // alignments for mate #1
		const AlnRes* rs2,      // alignments for mate #2
		const AlnSetSumm& summ, // summary
		const AlnFlags* flags1, // flags for mate #1
		const AlnFlags* flags2) // flags for mate #2
	{
		assert(rd1 != NULL || rd2 != NULL);
		if(rd1 != NULL) {
			appendMate(o, *rd1, rd2, rdid, rs1, rs2, summ, *flags1);
		}
		if(rd2 != NULL) {
			appendMate(o, *rd2, rd1, rdid, rs2, rs1, summ, *flags2);
		}
	}

protected:

	/**
	 * Append a single per-mate alignment result to the given output
	 * stream.  If the alignment is part of a pair, information about
	 * the opposite mate and its alignment are given in rdo/rso.
	 */
	void appendMate(
		OutFileBuf&   o,
		const Read&   rd,
		const Read*   rdo,
		const TReadId rdid,
		const AlnRes* rs,
		const AlnRes* rso,
		const AlnSetSumm& summ,
		const AlnFlags& flags);

	const EList<bool>& suppress_; // bit mask of columns to suppress
	int           offBase_;    // add this to 0-based reference offsets before printing
	bool          colorSeq_;   // colorspace: print color seq instead of decoded nucs
	bool          colorQual_;  // colorspace: print color quals instead of decoded quals
	bool          exEnds_;     // exclude ends for decoded colorspace alignments
	bool    printPlaceholders_;// print maxs and unals
	bool          printFlags_; // print alignment flags
	bool          printCost_;  // print penalty in extra column
	bool          printParams_;// print alignment parameters
	bool          printSeed_;  // print pseudo-random seed
	ReferenceMap* rmap_;       // reference coordinate transformation
	bool          fullRef_;    // print entire reference name including whitespace
	int           partition_;  // partition size
	
	BTDnaString   dseq_;       // buffer for decoded read sequence
	BTString      dqual_;      // buffer for decoded quality sequence
	
	EList<char>   tmpop_;      // temporary holder for CIGAR ops
	EList<size_t> tmprun_;     // temporary holder for CIGAR runs
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
		OutFileBuf*      out,        // initial output stream
		const Mapq&      mapq,       // mapping quality calculator
		const SamConfig& samc,       // settings & routines for SAM output
		bool             deleteOuts, // delete output objects upon destruction
		const StrList&   refnames,   // reference names
		bool             quiet,      // don't print alignment summary at end
		bool             exEnds,     // exclude ends for decoded colors alns
		ReferenceMap*    rmap) :     // reference coordinate transformation
		AlnSink(
			out,
			mapq,
			deleteOuts,
			refnames,
			quiet),
		samc_(samc),
		exEnds_(exEnds),
		rmap_(rmap)
	{ }

	/**
	 * Append a single alignment result, which might be paired or
	 * unpaired, to the given output stream in Bowtie's verbose-mode
	 * format.  If the alignment is paired-end, print mate1's alignment
	 * then mate2's alignment.
	 */
	virtual void append(
		OutFileBuf&   o,        // file buffer to write to
		const Read*   rd1,      // mate #1
		const Read*   rd2,      // mate #2
		const TReadId rdid,     // read ID
		const AlnRes* rs1,      // alignments for mate #1
		const AlnRes* rs2,      // alignments for mate #2
		const AlnSetSumm& summ, // summary
		const AlnFlags* flags1, // flags for mate #1
		const AlnFlags* flags2) // flags for mate #2
	{
		assert(rd1 != NULL || rd2 != NULL);
		if(rd1 != NULL) {
			assert(flags1 != NULL);
			appendMate(o, *rd1, rd2, rdid, rs1, rs2, summ, *flags1);
		}
		if(rd2 != NULL) {
			assert(flags2 != NULL);
			appendMate(o, *rd2, rd1, rdid, rs2, rs1, summ, *flags2);
		}
	}

protected:

	/**
	 * Append a single per-mate alignment result to the given output
	 * stream.  If the alignment is part of a pair, information about
	 * the opposite mate and its alignment are given in rdo/rso.
	 */
	void appendMate(
		OutFileBuf&   o,
		const Read&   rd,
		const Read*   rdo,
		const TReadId rdid,
		const AlnRes* rs,
		const AlnRes* rso,
		const AlnSetSumm& summ,
		const AlnFlags& flags);

	const SamConfig& samc_;    // settings & routines for SAM output
	bool             exEnds_;  // exclude ends for decoded colorspace alignments
	ReferenceMap*    rmap_;    // reference coordinate transformation
	
	BTDnaString      dseq_;    // buffer for decoded read sequence
	BTString         dqual_;   // buffer for decoded quality sequence
	
	EList<char>      tmpop_;   // temporary holder for CIGAR ops
	EList<size_t>    tmprun_;  // temporary holder for CIGAR runs
};

#endif /*ndef ALN_SINK_H_*/
