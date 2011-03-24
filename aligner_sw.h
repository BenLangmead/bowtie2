/*
 * aligner_sw.h
 *
 * Classes and routines for solving dynamic programming problems in aid of read
 * alignment.  Goals include the ability to handle:
 *
 * - Both nucleotide queries and colorspace queries
 * - Both read alignment, where the query must align end-to-end, and local
 *   alignment, where we seek a high-scoring alignment that need not involve
 *   the entire query.
 * - Situations where: (a) we've found a seed hit and are trying to extend it
 *   into a larger hit, (b) we've found an alignment for one mate of a pair and
 *   are trying to find a nearby alignment for the other mate, (c) we're
 *   aligning against an entire reference sequence.
 * - Caller-specified indicators for what columns of the dynamic programming
 *   matrix we are allowed to start in or end in.
 *
 * TODO:
 *
 * - A slicker way to filter out alignments that violate a ceiling placed on
 *   the number of Ns permitted in the reference portion of the alignment.
 *   Right now we accomplish this by masking out ending columns that correspond
 *   to *ungapped* alignments with too many Ns.  This results in false
 *   positives and false negatives for gapped alignments.  The margin of error
 *   (# of Ns by which we might miscount) is bounded by the number of gaps.
 */

#ifndef ALIGNER_SW_H_
#define ALIGNER_SW_H_

#include <stdint.h>
#include "ds.h"
#include "threading.h"
#include "aligner_seed.h"
#include "reference.h"
#include "random_source.h"
#include "mem_ids.h"
#include "aligner_result.h"
#include "aligner_sw_col.h"
#include "mask.h"

#define QUAL2(d, f) pen_->mm((int)(*rd_)[rdi_ + d], \
                             (int)  rf_ [rfi_ + f], \
						     (int)(*qu_)[rdi_ + d] - 33)
#define QUAL(d)     pen_->mm((int)(*rd_)[rdi_ + d], \
                             (int)(*qu_)[rdi_ + d] - 33)
#define N_SNP_PEN(c) (((int)rf_[rfi_ + c] > 15) ? pen_->n(30) : pen_->penSnp)

/**
 * Key parameters to dynamic programming alignment.
 */
struct SwParams {
	SwParams() { initFromGlobals(); }
	
	/**
	 * Initialize fields of SwParams by copying values of globals (see
	 * search_globals.h).
	 */
	void initFromGlobals();

	int   gapBar;  // no gaps allowed within 'gapBar' positions of either end
	bool  exEnds;  // true -> nucleotide alignment is 1 char shorter than color read
};

/**
 * Encapsulates the result of a dynamic programming alignment, including
 * colorspace alignments.  In our case, the result is a combination of:
 *
 * 1. All the nucleotide edits
 * 2. All the "edits" where an ambiguous reference char is resolved to
 *    an unambiguous char.
 * 3. All the color edits (if applicable)
 * 4. All the color miscalls (if applicable).  This is a subset of 3.
 * 5. The score of the best alginment
 * 6. The score of the second-best alignment
 *
 * Having scores for the best and second-best alignments gives us an
 * idea of where gaps may make reassembly beneficial.
 */
struct SwResult {

	SwResult() :
		alres(),
		sws(0),
		swcups(0),
		swrows(0),
		swskiprows(0),
		swsucc(0),
		swfail(0),
		swbts(0)
	{ }

	/**
	 * Clear all contents.
	 */
	void reset() {
		sws = swcups = swrows = swskiprows = swsucc =
		swfail = swbts = 0;
		alres.reset();
	}
	
	/**
	 * Reverse all edit lists.
	 */
	void reverse() {
		alres.reverseEdits();
	}
	
	/**
	 * Return true iff no result has been installed.
	 */
	bool empty() const {
		return alres.empty();
	}
	
	/**
	 * Check that result is internally consistent.
	 */
	bool repOk() const {
		assert(alres.repOk());
		return true;
	}

	AlnRes alres;
	uint64_t sws;    // # dynamic programming problems solved
	uint64_t swcups; // # dynamic programming cell updates
	uint64_t swrows; // # dynamic programming row updates
	uint64_t swskiprows; // # skipped dynamic programming row updates (b/c no valid alignments can go thru row)
	uint64_t swsucc; // # dynamic programming problems resulting in alignment
	uint64_t swfail; // # dynamic programming problems not resulting in alignment
	uint64_t swbts;  // # dynamic programming backtrace steps
	
	int nup;         // upstream decoded nucleotide; for colorspace reads
	int ndn;         // downstream decoded nucleotide; for colorspace reads
};

/**
 * Encapsulates counters that measure how much work has been done by
 * the dynamic programming driver and aligner.
 */
struct SwMetrics {

	SwMetrics() { reset(); MUTEX_INIT(lock); }
	
	void reset() {
		sws = swcups = swrows = swskiprows = swsucc = swfail = swbts =
		rshit = 0;
	}
	
	void init(
		uint64_t s1,
		uint64_t s2,
		uint64_t s3,
		uint64_t s4,
		uint64_t s5,
		uint64_t s6,
		uint64_t s7,
		uint64_t s8)
	{
		sws        = s1;
		swcups     = s2;
		swrows     = s3;
		swskiprows = s4;
		swsucc     = s5;
		swfail     = s6;
		swbts      = s7;
		rshit      = s8;
	}
	
	/**
	 * Merge (add) the counters in the given SwResult object into this
	 * SwMetrics object.
	 */
	void update(const SwResult& r) {
		sws        += r.sws;
		swcups     += r.swcups;
		swrows     += r.swrows;
		swskiprows += r.swskiprows;
		swsucc     += r.swsucc;
		swfail     += r.swfail;
		swbts      += r.swbts;
	}
	
	/**
	 * Merge (add) the counters in the given SwMetrics object into this
	 * object.  This is the only safe way to update a SwMetrics shared
	 * by multiple threads.
	 */
	void merge(const SwMetrics& r, bool getLock = false) {
		ThreadSafe ts(&lock, getLock);
		sws        += r.sws;
		swcups     += r.swcups;
		swrows     += r.swrows;
		swskiprows += r.swskiprows;
		swsucc     += r.swsucc;
		swfail     += r.swfail;
		swbts      += r.swbts;
	}

	uint64_t sws;    // # dynamic programming problems solved
	uint64_t swcups; // # dynamic programming cell updates
	uint64_t swrows; // # dynamic programming row updates
	uint64_t swskiprows; // # skipped dynamic programming row updates (b/c no valid alignments can go thru row)
	uint64_t swsucc; // # dynamic programming problems resulting in alignment
	uint64_t swfail; // # dynamic programming problems not resulting in alignment
	uint64_t swbts;  // # dynamic programming backtrace steps
	uint64_t rshit;  // # dynamic programming problems avoided b/c seed hit was redundant
	MUTEX_T lock;
};

/**
 * Counters characterizing work done by 
 */
struct SwCounters {
	uint64_t cups;    // cell updates
	uint64_t bts;     // backtracks
	
	/**
	 * Set all counters to 0.
	 */
	void reset() {
		cups = bts = 0;
	}
};

/**
 * Abstract parent class for encapsulating SeedAligner actions.
 */
struct SwAction {
};

/**
 * Abstract parent for a class with a method that gets passed every
 * set of counters for every join attempt.
 */
class SwCounterSink {
public:
	SwCounterSink() { MUTEX_INIT(lock_); }
	virtual ~SwCounterSink() { }
	/**
	 * Grab the lock and call abstract member reportCountersImpl()
	 */
	virtual void reportCounters(const SwCounters& c) {
		ThreadSafe(&this->lock_);
		reportCountersImpl(c);
	}
protected:
	virtual void reportCountersImpl(const SwCounters& c) = 0;
	MUTEX_T lock_;
};

/**
 * Write each per-SW set of counters to an output stream using a
 * simple record-per-line tab-delimited format.
 */
class StreamTabSwCounterSink : public SwCounterSink {
public:
	StreamTabSwCounterSink(std::ostream& os) : SwCounterSink(), os_(os) { }
protected:
	virtual void reportCountersImpl(const SwCounters& c)
	{
		os_ << c.cups << "\t"
			<< c.bts  << "\n"; // avoid 'endl' b/c flush is unnecessary
	}
	std::ostream& os_;
};

/**
 * Abstract parent for a class with a method that gets passed every
 * set of counters for every join attempt.
 */
class SwActionSink {
public:
	SwActionSink() { MUTEX_INIT(lock_); }
	virtual ~SwActionSink() { }
	/**
	 * Grab the lock and call abstract member reportActionsImpl()
	 */
	virtual void reportActions(const EList<SwAction>& as) {
		ThreadSafe(&this->lock_);
		reportActionsImpl(as);
	}
protected:
	virtual void reportActionsImpl(const EList<SwAction>& as) = 0;
	MUTEX_T lock_;
};

/**
 * Write each per-SW set of Actions to an output stream using a
 * simple record-per-line tab-delimited format.
 */
class StreamTabSwActionSink : public SwActionSink {
public:
	StreamTabSwActionSink(std::ostream& os) : SwActionSink(), os_(os) { }
	virtual ~StreamTabSwActionSink() { }
protected:
	virtual void reportActionsImpl(const EList<SwAction>& as)
	{
		for(size_t i = 0; i < as.size(); i++) {
			os_ << "\n"; // avoid 'endl' b/c flush is unnecessary
		}
	}
	std::ostream& os_;
};

enum {
	SW_BT_DIAG,
	SW_BT_REF_OPEN,
	SW_BT_REF_EXTEND,
	SW_BT_READ_OPEN,
	SW_BT_READ_EXTEND
};

/**
 * Encapsulates a bitmask.  The bitmask encodes which backtracking paths out of
 * a cell lie on optimal subpaths.
 */
struct SwNucCellMask {

	/**
	 * Set all flags to 0 (meaning either: there's no way to backtrack from
	 * this cell to an optimal answer, or we haven't set the mask yet)
	 */
	void clear() {
		*((uint8_t*)this) = 0;
	}

	/**
	 * Return true iff the mask is empty.
	 */
	inline bool empty() const {
		return *((uint8_t*)this) == 0;
	}

	/**
	 * Return true iff it's possible to extend a gap in the reference in the
	 * cell below this one.
	 */
	inline bool refExtendPossible() const {
		return rfop || rfex;
	}

	/**
	 * Return true iff it's possible to open a gap in the reference
	 * in the cell below this one (false implies that only extension
	 * is possible).
	 */
	inline bool refOpenPossible() const {
		return diag || rfop || rfex;
	}

	/**
	 * Return true iff it's possible to extend a gap in the read
	 * in the cell to the right of this one.
	 */
	inline bool readExtendPossible() const {
		return rdop || rdex;
	}

	/**
	 * Return true iff it's possible to open a gap in the read in the
	 * cell to the right of this one (false implies that only extension
	 * is possible).
	 */
	inline bool readOpenPossible() const {
		return diag || rdop || rdex;
	}
	
	/**
	 * Return true iff there is >0 possible way to backtrack from this
	 * cell.
	 */
	inline int numPossible() const {
		return diag + rfop + rfex + rdop + rdex;
	}

	/**
	 * Select a path for backtracking from this cell.  If there is a tie among
	 * eligible paths, break it randomly.  Return value is a flag indicating
	 * the backtrack type (see enum defining SW_BT_* above).
	 */
	int randBacktrack(RandomSource& rand);

	uint8_t diag     : 1;
	uint8_t rfop     : 1;
	uint8_t rfex     : 1;
	uint8_t rdop     : 1;
	uint8_t rdex     : 1;
	uint8_t reserved : 3;
};

/**
 * Encapsulates all information needed to encode the optimal subproblem
 * at a cell in a colorspace SW matrix.
 */
struct SwNucCell {

	/**
	 * Clear this cell so that it's ready for updates.
	 */
	void clear() {
		// Initially, best score is invalid
		best = AlignmentScore::INVALID();
		// Initially, there's no way to backtrack from this cell
		mask.clear();
		empty = true;
		ASSERT_ONLY(finalized = false);
	}

	/**
	 * Return true if best is valid.
	 */
	bool valid() const {
		return VALID_AL_SCORE(best);
	}
	
	/**
	 * Caller supplies a current-best and current-second best score and
	 * we update them according to the incoming scores for this cell.
	 */
	bool updateBest(AlignmentScore& bestSc, int penceil) const {
		if(best > bestSc) {
			assert_leq(abs(best.score()), penceil);
			bestSc = best;
			return true;
		}
		return false;
	}

	/**
	 * We finished updating the cell; set empty and finalized
	 * appropriately.
	 */
	bool finalize(int penceil) {
		ASSERT_ONLY(finalized = true);
		assert(empty);
		// Profiling shows cache misses on following line
		if(!mask.empty()) {
			assert(VALID_AL_SCORE(best));
			assert_leq(abs(best.score()), penceil);
			empty = false;
		}
		return !empty;
	}

	// Best incoming score for each 'to' character
	AlignmentScore best;
	// Mask for tied-for-best incoming paths for each 'to' character
	SwNucCellMask mask;
	
	bool empty;
	ASSERT_ONLY(bool finalized);
};

/**
 * SwAligner
 * =========
 *
 * Ensapsulates routines for performing dynamic programming alignments of
 * nucleotide or colorspace reads against reference nucleotides.  In the
 * colorspace case, decoding takes place simultaneously with alignment.
 *
 * The class is stateful.  First the user must call init() to initialize the
 * object with details regarding the dynamic programming problem to be solved.
 * Next, the user calls align() to fill the dynamic programming matrix and
 * calculate summaries describing the solutions.  Finally the user calls 
 * nextAlignment(...), perhaps repeatedly, to populate the SwResult object with
 * the next result.  Results are dispensend in best-to-worst, left-to-right
 * order.
 *
 * There is a design tradeoff between hiding/exposing details of the genome and
 * its strands to the SwAligner.  In a sense, a better design is to hide
 * details such as the id of the reference sequence aligned to, or whether
 * we're aligning the read in its original forward orientation or its reverse
 * complement.  But this means that any alignment results returned by SwAligner
 * have to be extended to include those details before they're useful to the
 * caller.  We opt for messy but expedient - the reference id and orientation
 * of the read are given to SwAligner, remembered, and used to populate
 * SwResults.
 */
class SwAligner {
public:

	SwAligner() :
		inited_(false),
		rfwbuf_(SW_CAT),
		maskst_(SW_CAT),
		masken_(SW_CAT),
		ntab_(SW_CAT),
		ctab_(SW_CAT),
		cursol_(0),
		sols_(SW_CAT)
	{ }

	/**
	 * Initialize with a new alignment problem.
	 */
	void init(
		const BTDnaString& rd, // read sequence
		const BTString& qu,    // read qualities
		size_t rdi,            // offset of first read char to align
		size_t rdf,            // offset of last read char to align
		bool fw,               // true iff read sequence is original fw read
		bool color,            // true iff read is colorspace
		uint32_t refidx,       // id of reference aligned against
		int64_t refoff,        // offset of upstream ref char aligned against
		char *rf,              // reference sequence
		size_t rfi,            // offset of first reference char to align to
		size_t rff,            // offset of last reference char to align to
		size_t width,          // # bands to do (width of parallelogram)
		const EList<bool>* st, // mask indicating which columns we can start in
		const EList<bool>* en, // mask indicating which columns we can end in
		const SwParams& pa,    // params for SW alignment
		const Penalties& pen,  // penalties for edit types
		int penceil)           // penalty ceiling for valid alignments
	{
		int readGaps = pen.maxReadGaps(penceil);
		int refGaps  = pen.maxRefGaps(penceil);
		int maxGaps  = max(readGaps, refGaps);
		int nceil    = (int)pen.nCeil(rd.length());
		assert_geq(readGaps, 0);
		assert_geq(refGaps, 0);
		rd_      = &rd;        // read sequence
		qu_      = &qu;        // read qualities
		rdi_     = rdi;        // offset of first read char to align
		rdf_     = rdf;        // offset of last read char to align
		fw_      = fw;         // true iff read sequence is original fw read
		color_   = color;      // true iff read is colorspace
		refidx_  = refidx;     // id of reference aligned against
		refoff_  = refoff;     // offset of upstream ref char aligned against
		rf_      = rf;         // reference sequence
		rfi_     = rfi;        // offset of first reference char to align to
		rff_     = rff;        // offset of last reference char to align to
		rdgap_   = readGaps;   // max # gaps in read
		rfgap_   = refGaps;    // max # gaps in reference
		maxgap_  = maxGaps;    // max(readGaps, refGaps)
		width_   = width;      // # bands to do (width of parallelogram)
		st_      = st;         // mask indicating which columns we can start in
		en_      = en;         // mask indicating which columns we can end in
		pa_      = &pa;        // params for SW alignment
		pen_     = &pen;       // penalties for edit types
		penceil_ = penceil;    // penalty ceiling for valid alignments
		nceil_   = nceil;      // max # Ns allowed in ref portion of aln
		inited_  = true;       // indicate we're initialized now
	}
	
	/**
	 * Given a read, an alignment orientation, a range of characters in a
	 * referece sequence, and a bit-encoded version of the reference,
	 * execute the corresponding dynamic programming problem.
	 *
	 * Here we expect that the caller has already narrowed down the relevant
	 * portion of the reference (e.g. using a seed hit) and all we do is
	 * banded dynamic programming in the vicinity of that portion.  This is not
	 * the function to call if we are trying to solve the whole alignment
	 * problem with dynamic programming (that is TODO).
	 *
	 * Returns true if an alignment was found, false otherwise.
	 */
	void init(
		const Read& rd,        // read to align
		size_t rdi,            // off of first char in 'rd' to consider
		size_t rdf,            // off of last char (excl) in 'rd' to consider
		bool fw,               // whether to align forward or revcomp read
		bool color,            // colorspace?
		uint32_t refidx,       // reference aligned against
		int64_t rfi,           // off of first character in ref to consider
		int64_t rff,           // off of last char (excl) in ref to consider
		const BitPairReference& refs, // Reference strings
		size_t reflen,         // length of reference sequence
		size_t width,          // # bands to do (width of parallelogram)
		const EList<bool>* st, // mask indicating which columns we can start in
		const EList<bool>* en, // mask indicating which columns we can end in
		const SwParams& pa,    // dynamic programming parameters
		const Penalties& pen,  // penalty scheme
		int penceil);          // penalty ceiling for valid alignments
	
	/**
	 * Align read 'rd' to reference using read & reference information given
	 * last time init() was called.  If the read is colorspace, the decoding is
	 * determined simultaneously with alignment.  Uses dynamic programming.
	 */
	bool align(
		SwResult& res,
		RandomSource& rnd);
	
	/**
	 * Return the next alignment, if there is one.
	 */
	void nextAlignment(SwResult& res);
	
	/**
	 * Return the estimated number of distinct alignments uncovered by the
	 * dynamic programming matrix.
	 */
	size_t numAlignments() const {
		return sols_.size();
	}
	
	/**
	 * Return true iff this SwAligner has been initialized with a dynamic
	 * programming problem.
	 */
	bool inited() const { return inited_; }
	
	/**
	 * Reset, signaling that we're done with this dynamic programming problem
	 * and won't be asking for any more alignments.
	 */
	void reset() { inited_ = false; }

protected:

	/**
	 * Align and simultaneously decode the colorspace read 'rd' to the
	 * reference string 'rf' using dynamic programming.
	 */
	int alignNucleotides(
		SwResult& res,         // store results (edits and scores) here
		RandomSource& rnd);    // pseudo-random generator

	/**
	 * Align the nucleotide read 'rd' to the reference string 'rf' using
	 * dynamic programming.
	 */
	int alignColors(
		SwResult& res,         // store results (edits and scores) here
		RandomSource& rnd);    // pseudo-random generator

	/**
	 * Given the dynamic programming table, trace backwards from the lower
	 * right-hand corner and populate 'decoded' 'nedits', 'cedits'
	 * accordingly.
	 */
	int backtrackColors(
		AlignmentScore dscore, // score we expect to get over backtrack
		SwResult& res,         // store results (edits and scores) here
		int col,               // start in this column
		int lastC,             // character to backtrace from in lower-right corner
		RandomSource& rand);   // pseudo-random generator

	/**
	 * Given the dynamic programming table, trace backwards from the lower
	 * right-hand corner and populate 'decoded' 'nedits', 'cedits'
	 * accordingly.
	 */
	int backtrackNucleotides(
		AlignmentScore dscore, // score we expect to get over backtrack
		SwResult& res,         // store results (edits and scores) here
		int col,               // start in this column
		RandomSource& rand);   // pseudo-random generator

	// Members for updating cells in nucleotide dynamic programming tables

	inline void updateNucHoriz(
		const SwNucCell& lc,
		SwNucCell& dstc,
		int rfm);

	inline void updateNucDiag(
		const SwNucCell& dc,
		SwNucCell& dstc,
		int rdc,
		int rfm,
		int pen);

	inline void updateNucVert(
		const SwNucCell& uc,
		SwNucCell& dstc,
		int rdc);

	// Members for updating cells in colorspace dynamic programming tables

	inline void updateColorHoriz(
		const SwColorCell& lc,
		SwColorCell& dstc,
		int refMask);

	inline void updateColorDiag(
		const SwColorCell& uc,
		SwColorCell& dstc,
		int refMask,
		int prevColor,
		int prevQual);

	inline void updateColorVert(
		const SwColorCell& uc,
		SwColorCell& dstc,
		int prevColor,
		int prevQual);

	const BTDnaString*  rd_;     // read sequence
	const BTString*     qu_;     // read qualities
	size_t              rdi_;    // offset of first read char to align
	size_t              rdf_;    // offset of last read char to align
	bool                fw_;     // true iff read sequence is original fw read
	bool                color_;  // true iff read is colorspace
	uint32_t            refidx_; // id of reference aligned against
	int64_t             refoff_; // offset of upstream ref char aligned against
	char               *rf_;     // reference sequence
	size_t              rfi_;    // offset of first ref char to align to
	size_t              rff_;    // offset of last ref char to align to
	size_t              width_;  // # bands to do (width of parallelogram)
	const EList<bool>*  st_;     // mask indicating which cols we can start in
	const EList<bool>*  en_;     // mask indicating which cols we can end in
	int                 rdgap_;  // max # gaps in read
	int                 rfgap_;  // max # gaps in reference
	int                 maxgap_; // max(rdgap_, rfgap_)
	const SwParams*     pa_;     // params for SW alignment
	const Penalties*    pen_;    // penalties for edit types
	int                 penceil_;// penalty ceiling for valid alignments
	int                 nceil_;  // max # Ns allowed in ref portion of aln

	bool                inited_; // true iff initialized with DP problem
	EList<uint32_t>     rfwbuf_; // buffer for wordized refernece stretches
	EList<bool>         maskst_; // list of bools: which cols can we start in?
	EList<bool>         masken_; // list of bools: which cols can we end in?
	ELList<SwNucCell>   ntab_;   // dynamic programming table for nucleotide SW
	ELList<SwColorCell> ctab_;   // dynamic programming table for colorspace SW
	size_t              cursol_; // idx of next solution to dish out
	EList<size_t>       sols_;   // list of elts in last row ending solutions
};

#endif /*ALIGNER_SW_H_*/
