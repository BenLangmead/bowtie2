//
//  aln_sink.h
//

#ifndef ALN_SINK_H_
#define ALN_SINK_H_

#include "filebuf.h"
#include "read.h"
#include "read_sink.h"
#include "unique.h"
#include "refmap.h"
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
	
	void reset() { init(0, 0, 0); }
	
	void init(
		uint64_t al_,
		uint64_t unal_,
		uint64_t max_)
	{
		al = al_;
		unal = unal_;
		max = max_;
	}
	
	/**
	 * Merge (add) the counters in the given ReportingMetrics object
	 * into this object.  This is the only safe way to update a
	 * ReportingMetrics shared by multiple threads.
	 */
	void merge(const ReportingMetrics& met, bool getLock = false) {
		ThreadSafe ts(&lock, getLock);
		al   += met.al;
		unal += met.unal;
		max  += met.max;
	}

	uint64_t al;   // # reads w/ >= 1 alignment
	uint64_t unal; // # reads w/ 0 alignments
	uint64_t max;  // # reads w/ more alignments than the -M/-m ceiling
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
		bool msample_)
	{
		init(khits_, mhits_, pengap_, msample_);
	}

	void init(
		THitInt khits_,
		THitInt mhits_,
		THitInt pengap_,
		bool msample_)
	{
		khits = khits_;     // -k (or high if -a)
		mhits = mhits_;     // -m or -M
		pengap = pengap_;
		msample = msample_;
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
	bool mhitsSet() const {
		return mhits < std::numeric_limits<THitInt>::max();
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

public:

	explicit AlnSink(
		OutFileBuf*               out,
		const EList<bool>&        suppress,
		ReadSink*                 readSink,
		const Mapq&               mapq,       // mapping quality calculator
		bool                      deleteOuts,
		const EList<std::string>* refnames,
		bool                      quiet) :
		outs_(),
		outNames_(),
		suppress_(suppress),
		locks_(),
		deleteOuts_(deleteOuts),
		numWrappers_(0),
		numAligned_(0llu),
		numUnaligned_(0llu),
		numMaxed_(0llu),
		numReported_(0llu),
		numReportedPaired_(0llu),
		refnames_(refnames),
		quiet_(quiet),
		readSink_(readSink),
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
		const AlnSetSumm&  summ) = 0;

	/**
	 * Report a given batch of hits for the given read pair.  Should be
	 * called just once per read pair.
	 */
	virtual void reportHits(
		const Read           *rd1,
		const Read           *rd2,
		const TReadId         rdid,
		const EList<bool>&    select,
		const EList<AlnRes>  *rs1,
		const EList<AlnRes>  *rs2,
		bool                  maxed,
		const AlnSetSumm&     summ,
		bool                  getLock = true)
	{
		assert(rd1 != NULL);
		assert(rs1 != NULL);
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
			rs1->size(),
			summ,
			getLock);
	}

	/**
	 * Report a given batch of hits for the given read pair.  Should be
	 * called just once per read pair.
	 */
	virtual void reportHits(
		const Read          *rd1,
		const Read          *rd2,
		const TReadId        rdid,
		const EList<bool>&   select,
		const EList<AlnRes> *rs1,
		const EList<AlnRes> *rs2,
		bool                 maxed,
		size_t               start,
		size_t               end,
		const AlnSetSumm&    summ,
		bool                 getLock = true)
	{
		assert(rd1 != NULL);
		assert(rs1 != NULL);
		assert_geq(end, start);
		assert_leq(end, rs1->size());
		size_t num = end - start;
		if(num == 0) {
			// Nothing to report
			return;
		}
		bool paired = (rs2 != NULL && !rs2->empty() && !rs2->get(0).empty());
		size_t reported = 0;
		for(size_t i = start; i < end; i++) {
			// Skip if it hasn't been selected
			if(!select[i]) continue;
			// Determine the stream id using the coordinate of the
			// upstream mate
			size_t sid = streamId(rdid, rs1->get(i).refcoord());
			assert_lt(sid, locks_.size());
			ThreadSafe ts(&locks_[sid], getLock);
			const AlnRes* r1 = &rs1->get(i);
			const AlnRes* r2 = NULL;
			if(paired) {
				assert_eq(rs1->size(), rs2->size());
				r2 = &rs2->get(i);
			}
			append(out(sid), rd1, rd2, rdid, r1, r2, summ);
			reported++;
		}
		readSink_->dumpAlign(rd1, rd2, rdid);
		{
			ThreadSafe ts(&mainlock_);
			commitHits(rd1, rd2, rdid, rs1, rs2, start, end, getLock);
			numAligned_++;
			if(paired) numReportedPaired_ += reported;
			else       numReported_       += reported;
		}
	}

	/**
	 * Report a read that aligned more times than allowed by the -m or
	 * -M ceiling.
	 */
	virtual void reportMaxed(
		const Read*          rd1,
		const Read*          rd2,
		const TReadId        rdid,
		const EList<AlnRes>* rs1,
		const EList<AlnRes>* rs2,
		const AlnSetSumm&    summ,
		bool                 getLock = true)
	{
		assert(rd1 != NULL);
		assert(rs1 != NULL);
		reportMaxed(rd1, rd2, rdid, rs1, rs2, 0, rs1->size(), summ, getLock);
	}

	/**
	 * Report a read that aligned more times than allowed by the -m or
	 * -M ceiling.
	 */
	virtual void reportMaxed(
		const Read*          rd1,
		const Read*          rd2,
		const TReadId        rdid,
		const EList<AlnRes>* rs1,
		const EList<AlnRes>* rs2,
		size_t               start,
		size_t               end,
		const AlnSetSumm&    summ,
		bool                 getLock = true)
	{
		assert(rd1 != NULL);
		assert(rs1 != NULL);
		readSink_->dumpMaxed(rd1, rd2, rdid);
		ThreadSafe ts(&mainlock_, getLock);
		numMaxed_++;
	}

	/**
	 * Report an unaligned read.  Typically we do nothing, but we might
	 * want to print a placeholder when output is chained.
	 */
	virtual void reportUnaligned(
		const Read*   rd1,
		const Read*   rd2,
		const TReadId rdid,
		bool          getLock = true)
	{
		readSink_->dumpUnal(rd1, rd2, rdid);
		ThreadSafe ts(&mainlock_, getLock);
		numUnaligned_++;
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
		uint64_t al,
		uint64_t un,
		uint64_t mx,
		uint64_t rep,
		uint64_t repp,
		bool sample,
		bool hadoopOut);

	/**
	 * Called when all alignments are complete.  It is assumed that no
	 * synchronization is necessary.
	 */
	void finish(bool hadoopOut, bool sampleMax) {
		// Close output streams
		closeOuts(false);
		if(!quiet_) {
			printAlSumm(
				numAligned_,
				numUnaligned_,
				numMaxed_,
				numReported_,
				numReportedPaired_,
				sampleMax,
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
		if(readSink_ != NULL) readSink_->dumpMaxed(m1, m2, rdid);
	}
	
	void dumpUnal(const Read* m1, const Read* m2, TReadId rdid) {
		if(readSink_ != NULL) readSink_->dumpUnal(m1, m2, rdid);
	}
	
	void dumpAlign(const Read* m1, const Read* m2, TReadId rdid) {
		if(readSink_ != NULL) readSink_->dumpAlign(m1, m2, rdid);
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
	EList<bool>        suppress_;     // suppress columns
	EList<MUTEX_T>     locks_;        // pthreads mutexes for per-file critical sections
	bool               deleteOuts_;   // Whether to delete elements of outs_ upon exit
	int                numWrappers_;  // # threads owning a wrapper for this HitSink
	MUTEX_T            mainlock_;     // pthreads mutexes for fields of this object
	volatile uint64_t  numAligned_;   // # reads with >= 1 alignment
	volatile uint64_t  numUnaligned_; // # reads with no alignments
	volatile uint64_t  numMaxed_;     // # reads with # alignments exceeding -m ceiling
	volatile uint64_t  numReported_;  // # single-ended alignments reported
	volatile uint64_t  numReportedPaired_; // # paired-end alignments reported
	const EList<std::string>* refnames_; // reference names
	bool               quiet_;        // true -> don't print alignment stats at the end
	ReadSink*          readSink_;     // 
	const Mapq&        mapq_;         // mapping quality calculator
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
		short_(false), // memory of whether lastStage_ was short-circuited
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
		lastStage_(-1) // memory of what alignment stages we've made it thru
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
	 * Finish reporting for the read in rd1_ or in rd2_, depending on
	 * how 'one' is set.
	 *
	 * Called by finishRead.
	 */
	void finishGroup(
		bool paired,        // true iff alns being reported are paired
		bool condord,       // true iff paired-end alns are concordant
		bool one,           // true iff unpaired alns are from mate 1
		RandomSource& rnd,  // pseudo-random generator
		uint64_t& al,       // counter to inc for reads that align
		uint64_t& mx,       // counter to inc for reads that align repetitively
		uint64_t& un) const;// counter to inc for reads that don't align
	
	/**
	 * Inform global, shared AlnSink object that we're finished with
	 * this read.  The global AlnSink is responsible for updating
	 * counters, creating the output record, and delivering the record
	 * to the appropriate output stream.
	 */
	void finishRead(
		const SeedResults *sr1,
		const SeedResults *sr2,
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
	 * Caller uses this function to indicate that all stages up to and
	 * including the given stage have been attempted and no alignments
	 * were found.  Note that we assume that the aligner proceeds
	 * through each stage in order.
	 */
	void finishStage(int stage) {
		assert_gt(stage, lastStage_);
		lastStage_ = stage;
	}

	/**
	 * Check that hit sink wrapper is internally consistent.
	 */
	bool repOk() const {
		assert_geq(lastStage_, -1);
		assert_eq(rs2_.size(), rs1_.size());
		if(rp_.mhitsSet()) {
			assert_leq((int)rs1_.size()-1, rp_.mhits);
			assert_leq((int)rs2_.size()-1, rp_.mhits);
			assert(readIsPair() || (int)rs1u_.size()-1 <= rp_.mhits);
			assert(readIsPair() || (int)rs2u_.size()-1 <= rp_.mhits);
		}
		if(init_) {
			assert(rd1_ != NULL);
			assert_neq(std::numeric_limits<TReadId>::max(), rdid_);
		} else {
		
		}
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
	 * Return true iff we have already encountered a number of
	 * alignments that exceeds the -m/-M ceiling.
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

protected:

	/**
	 * If there is a configuration of unpaired alignments that fits our
	 * criteria for there being one or more discordant alignments, then
	 * shift the discordant alignments over to the rs1_/rs2_ lists and
	 * return true.  Otherwise, return false.
	 */
	bool prepareDiscordants() {
		if(rs1u_.size() == 1 && rs2u_.size() == 1) {
			assert(rs1_.empty());
			assert(rs2_.empty());
			rs1_.push_back(rs1u_[0]);
			rs2_.push_back(rs2u_[0]);
			return true;
		}
		return false;
	}

	/**
	 * Given that rd1_/rd2_/rs1_/rs2_ are already populated with
	 * information about the input reads and their respective
	 * alignments, consider the alignment policy and make random
	 * selections where necessary.  E.g. if we found 10 alignments and
	 * the policy is -k 2 -m 20, select 2 alignments at random.
	 */
	void selectAlnsToReport(
		const EList<AlnRes>& rs,     // alignments to select from
		EList<bool>&         select, // list to put results in
		RandomSource&        rnd)
		const
	{
		assert(init_);
		assert(repOk());
		assert(!maxed() || rp_.msample);
		size_t sz = rs.size();
		if(sz < 1) return;
		select.resize(sz);
		if((size_t)rp_.khits < sz) {
			// Select a random offset into the list of alignments
			uint32_t off = rnd.nextU32() % (uint32_t)sz;
			size_t take = rp_.khits;
			// Now take rp_.khits elements starting at that offset,
			// wrapping back to 0 if necessary, and leave the rest.
			for(size_t i = 0; i < sz; i++) {
				off++;
				if(off == sz) off = 0;
				select[off] = i < take;
			}
		} else {
			// Select them all!  No randomness needed.
			for(size_t i = 0; i < sz; i++) {
				select[i] = true;
			}
		}
	}

	AlnSink&        g_;     // global alignment sink
	ReportingParams rp_;    // reporting parameters: khits, mhits etc
	bool            init_;  // whether we're initialized w/ read pair
	bool            maxed1_; // true iff # unpaired mate-1 alns reported so far exceeded -m/-M
	bool            maxed2_; // true iff # unpaired mate-2 alns reported so far exceeded -m/-M
	bool            maxedOverall_; // true iff # paired-end alns reported so far exceeded -m/-M
	bool            short_; // true iff report() returned true already
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
	EList<bool>     select_;    // parallel to rs1_/rs2_ - which to report
	int             lastStage_; // set to the last stage finished
};

/**
 * An AlnSink concrete subclass for printing Bowtie verbose-mode alignments.
 */
class AlnSinkVerbose : public AlnSink {

public:

	AlnSinkVerbose(
		OutFileBuf*               out,        // initial output stream
		const EList<bool>&        suppress,   // suppress columns
		ReadSink*                 readSink,   // read sink
		const Mapq&               mapq,       // mapping quality calculator
		bool                      deleteOuts, // whether to delete output objects upon destruction
		const EList<std::string>* refnames,   // reference names
		bool                      quiet,      // don't print alignment summary at end
		int                       offBase,    // add this to 0-based offsets before printing
		bool                      colorSeq,   // colorspace: print color seq instead of decoded nucs
		bool                      colorQual,  // colorspace: print color quals instead of decoded quals
		bool                      exEnds,     // exclude ends for decoded colorspace alignments
		bool                      printCost,  // print penalty in extra column
		bool                      printParams,// print alignment parameters
		ReferenceMap*             rmap,       // reference coordinate transformation
		bool                      fullRef,    // print entire reference name including whitespace
		int                       partition = 0) : // partition size
		AlnSink(
			out,
			suppress,
			readSink,
			mapq,
			deleteOuts,
			refnames,
			quiet),
		offBase_(offBase),
		colorSeq_(colorSeq),
		colorQual_(colorQual),
		exEnds_(exEnds),
		printCost_(printCost),
		printParams_(printParams),
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
		OutFileBuf&   o,
		const Read*   rd1,
		const Read*   rd2,
		const TReadId rdid,
		const AlnRes* rs1,
		const AlnRes* rs2,
		const AlnSetSumm& summ)
	{
		assert(rd1 != NULL);
		assert(rs1 != NULL);
		appendMate(o, *rd1, rd2, rdid, *rs1, rs2, summ);
		if(rd2 != NULL) {
			assert(rs2 != NULL);
			appendMate(o, *rd2, rd1, rdid, *rs2, rs1, summ);
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
		const AlnRes& rs,
		const AlnRes* rso,
		const AlnSetSumm& summ);

	int           offBase_;    // add this to 0-based reference offsets before printing
	bool          colorSeq_;   // colorspace: print color seq instead of decoded nucs
	bool          colorQual_;  // colorspace: print color quals instead of decoded quals
	bool          exEnds_;     // exclude ends for decoded colorspace alignments
	bool          printCost_;  // print penalty in extra column
	bool          printParams_;// print alignment parameters
	ReferenceMap* rmap_;       // reference coordinate transformation
	bool          fullRef_;    // print entire reference name including whitespace
	int           partition_;  // partition size
	
	BTDnaString   dseq_;       // buffer for decoded read sequence
	BTString      dqual_;      // buffer for decoded quality sequence
};

#endif /*ndef ALN_SINK_H_*/
