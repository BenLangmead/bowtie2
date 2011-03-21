/*
 *  aligner_sw.h
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

/**
 * Return 1 if a 2-bit-encoded base ('i') match any bit in the mask
 * ('j') and the mask < 16.  Returns -1 if either the reference or the
 * read character was ambiguous.  Returns 0 if the characters
 * unambiguously mismatch.
 */
static inline int matches(int i, int j) {
	if(j >= 16 || i > 3) return -1;
	return (((1 << i) & j) != 0) ? 1 : 0;
}

static inline ostream& operator<<(ostream& os, const AlignmentScore& o) {
	os << o.score();
	return os;
}

#define QUAL2(d, f) pen.mm((int)rd[rdi + d], (int)rf[rfi + f], (int)qu[rdi + d] - 33)
#define QUAL(d)     pen.mm((int)rd[rdi + d], (int)qu[rdi + d] - 33)
#define N_SNP_PEN(c) (((int)rf[rfi + c] > 15) ? pen.n(30) : pen.penSnp)

/**
 * Key parameters to Smith-Waterman alignment.
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
 * Encapsulates the result of a Smith-Waterman alignment, including
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
};

/**
 * Encapsulates counters that measure how much work has been done by
 * the Smith-Waterman driver and aligner.
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
 * A bitmask encoding which backtracking paths out of a particular cell
 * correspond to optimal subpaths.
 */
struct SwNucCellMask {

	/**
	 * Set all flags to 0, indicating there is no way to backtrack from
	 * this cell to an optimal answer.
	 */
	void clear() {
		*((uint8_t*)this) = 0;
	}

	/**
	 * Return true iff there are no backward paths recorded in this
	 * mask.
	 */
	inline bool empty() const {
		return *((uint8_t*)this) == 0;
	}

	/**
	 * Return true iff it's possible to extend a gap in the reference
	 * in the cell below this one.
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
	 * Select a path for backtracking from this cell.  If there is a
	 * tie among eligible paths, break it randomly.  Return value is
	 * a pair where first = a flag indicating the backtrack type (see
	 * enum defining SW_BT_* above), and second = a selection for
	 * the read character for the next row up.  second should be
	 * ignored if the backtrack type is a gap in the read.
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

	inline void updateHoriz(
		const SwNucCell& lc,
		int rfm,
		const Penalties& pen,
		int nceil,
		int penceil);

	inline void updateDiag(
		const SwNucCell& dc,
		int rdc,
		int rfm,
		int pen,
		const Penalties& pens,
		int nceil,
		int penceil);

	inline void updateVert(
		const SwNucCell& uc,
		int rdc,
		const Penalties& pen,
		int nceil,
		int penceil);

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
 * Ensapsulates routines for performing Smith-Waterman alignments of
 * nucleotide or colorspace reads against reference nucleotides.  In
 * the colorspace case, decoding takes place simultaneously with
 * alignment.
 */
class SwAligner {
public:

	SwAligner() :
		rfbuf_(SW_CAT),
		mayend_(SW_CAT),
		ntab_(SW_CAT),
		ctab_(SW_CAT)
	{ }

	/**
	 * Given a list of reads, execute the corresponding Smith-Waterman
	 * problem across the entire extent of the reference.
	 */
	void alignBatchToFasta(
		const EList<Read>& rds,
		const EList<std::string>& fas,
		ELList<SwResult>& res);
	
	/**
	 * Given a read, an alignment orientation, a range of characters in a
	 * referece sequence, and a bit-encoded version of the reference,
	 * execute the corresponding Smith-Waterman problem.
	 *
	 * Here we expect that the caller has already narrowed down the
	 * relevant portion of the reference (e.g. using a seed hit) and all we
	 * have to do is a banded Smith-Waterman in the vicinity of that
	 * portion.  This is not the function to call if we are trying to solve
	 * the whole alignment problem with SW (see alignToFasta).
	 *
	 * Returns true if an alignment was found, false otherwise.
	 */
	bool alignToBitPairReference(
		const Read& rd,       // read to align
		bool color,           // colorspace?
		size_t rdi,           // offset of first character within 'read' to consider
		size_t rdf,           // offset of last char (exclusive) in 'read' to consider
		bool fw,              // whether to align forward or revcomp read
		uint32_t refidx,      // reference aligned against
		int64_t rfi,          // first reference base to SW align against
		int64_t rff,          // last reference base (exclusive) to SW align against
		const BitPairReference& refs, // Reference strings
		size_t reflen,        // length of reference sequence
		const SwParams& pa,   // Smith-Waterman parameters
		const Penalties& pen, // penalty scheme
		int penceil,          // maximum penalty we can incur for a valid alignment
		SwResult& res,        // results of Smith-Waterman alignment (score, edits, etc)
		EList<SwCounterSink*>* swCounterSinks, // send counter updates to these
		EList<SwActionSink*>* swActionSinks);  // send action-list updates to these

	/**
	 * Align and simultaneously decode the colorspace read 'read' as
	 * aligned against the reference string 'rf' using dynamic
	 * programming.
	 */
	int alignNucleotides(
		const BTDnaString& rd, // read sequence
		const BTString& qu,    // read qualities
		size_t rdi,            // offset of first character within 'read' to consider
		size_t rdf,            // offset of last char (exclusive) in 'read' to consider
		BTString& rf,          // reference sequence, as masks
		size_t rfi,            // offset of first character within 'rf' to consider
		size_t rff,            // offset of last char (exclusive) in 'rf' to consider
		const SwParams& pa,    // parameters governing Smith-Waterman problem
		const Penalties& pen,  // penalties for various edits
		int penceil,           // penalty ceiling for valid alignments
		SwResult& res,         // edits and scores
		RandomSource& rnd);    // pseudo-random generator

	/**
	 * Align the nucleotide read 'read' as aligned against the reference
	 * string 'rf' using dynamic programming.
	 */
	int alignColors(
		const BTDnaString& rd, // read sequence
		const BTString& qu,    // read qualities
		size_t rdi,            // offset of first character within 'read' to consider
		size_t rdf,            // offset of last char (exclusive) in 'read' to consider
		BTString& rf,          // reference sequence, as masks
		size_t rfi,            // offset of first character within 'rf' to consider
		size_t rff,            // offset of last char (exclusive) in 'rf' to consider
		const SwParams& pa,    // parameters governing Smith-Waterman problem
		const Penalties& pen,  // penalties for various edits
		int penceil,           // penalty ceiling for valid alignments
		int& nup,              // upstream decoded nucleotide
		int& ndn,              // downstream decoded nucleotide
		SwResult& res,         // edits and scores
		RandomSource& rnd);    // pseudo-random generator

protected:

	/**
	 * Given the dynamic programming table, trace backwards from the lower
	 * right-hand corner and populate 'decoded' 'nedits', 'cedits'
	 * accordingly.
	 */
	int backtrackColors(
		const BTDnaString& rd, // read sequence
		const BTString& qu,    // read qualities
		size_t rdi,            // offset of first read char to align
		size_t rdf,            // offset of last read char to align
		BTString& rf,          // reference sequence
		size_t rfi,            // offset of first reference char to align to
		size_t rff,            // offset of last reference char to align to
		AlignmentScore decodeScore,   // score we expect to get over backtrack
		int readGaps,          // max # gaps in read
		int refGaps,           // max # gaps in ref
		const SwParams& pa,    // params for SW alignment
		const Penalties& pen,  // penalties for edit types
		int& nup,              // upstream decoded nucleotide
		int& ndn,              // downstream decoded nucleotide
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
		const BTDnaString& rd, // read sequence
		const BTString& qu,    // read qualities
		size_t rdi,            // offset of first read char to align
		size_t rdf,            // offset of last read char to align
		BTString& rf,          // reference sequence
		size_t rfi,            // offset of first reference char to align to
		size_t rff,            // offset of last reference char to align to
		AlignmentScore decodeScore,   // score we expect to get over backtrack
		int readGaps,          // max # gaps in read
		int refGaps,           // max # gaps in ref
		const SwParams& pa,    // params for SW alignment
		const Penalties& pen,  // penalties for edit types
		SwResult& res,         // store results (edits and scores) here
		int col,               // start in this column
		RandomSource& rand);   // pseudo-random generator

	EList<uint32_t>     rfbuf_; // buffer for refernece stretches
	EList<bool>         mayend_;// list of bools: which columns can we end in?
	ELList<SwNucCell>   ntab_;  // dynamic programming table for nucleotide SW
	ELList<SwColorCell> ctab_;  // dynamic programming table for colorspace SW
};

#endif /*ALIGNER_SW_H_*/
