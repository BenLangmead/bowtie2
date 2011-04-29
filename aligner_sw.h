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
#include <iostream>
#include "aligner_sw_common.h"
#include "aligner_sw_nuc.h"
#include "aligner_sw_col.h"
#include "ds.h"
#include "threading.h"
#include "aligner_seed.h"
#include "reference.h"
#include "random_source.h"
#include "mem_ids.h"
#include "aligner_result.h"
#include "mask.h"

#define QUAL2(d, f) sc_->mm((int)(*rd_)[rdi_ + d], \
							(int)  rf_ [rfi_ + f], \
							(int)(*qu_)[rdi_ + d] - 33)
#define QUAL(d)     sc_->mm((int)(*rd_)[rdi_ + d], \
							(int)(*qu_)[rdi_ + d] - 33)
#define N_SNP_PEN(c) (((int)rf_[rfi_ + c] > 15) ? sc_->n(30) : sc_->penSnp)

/**
 * SwAligner
 * =========
 *
 * Ensapsulates facilities for alignment using dynamic programming.  Handles
 * alignment of nucleotide or colorspace reads against known reference
 * nucleotides.  In the colorspace case, decoding takes place simultaneously
 * with alignment.
 *
 * The class is stateful.  First the user must call init() to initialize the
 * object with details regarding the dynamic programming problem to be solved.
 * Next, the user calls align() to fill the dynamic programming matrix and
 * calculate summaries describing the solutions.  Finally the user calls 
 * nextAlignment(...), perhaps repeatedly, to populate the SwResult object with
 * the next result.  Results are dispensend in best-to-worst, left-to-right
 * order.
 *
 * The class expects the read string, quality string, and reference string
 * provided by the caller live at least until the user is finished aligning and
 * obtaining alignments from this object.
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
 *
 * LOCAL VS GLOBAL
 *
 * The dynamic programming aligner supports both local and global alignment,
 * and one option in between.  To implement global alignment, the aligner (a)
 * allows negative scores (i.e. doesn't necessarily clamp them up to 0), (b)
 * checks in rows other than the last row for acceptable solutions, and (c)
 * optionally adds a bonus to the score for matches.
 * 
 * For global alignment, we:
 *
 * (a) Allow negative scores
 * (b) Check only in the last row
 * (c) Either add a bonus for matches or not (doesn't matter)
 *
 * For local alignment, we:
 *
 * (a) Clamp scores to 0
 * (b) Check in any row for a sufficiently high score
 * (c) Add a bonus for matches
 *
 * An in-between solution is to allow alignments to be curtailed on the
 * right-hand side if a better score can be achieved thereby, but not on the
 * left.  For this, we:
 *
 * (a) Allow negative scores
 * (b) Check in any row for a sufficiently high score
 * (c) Either add a bonus for matches or not (doesn't matter)
 *
 * REDUNDANT ALIGNMENTS
 *
 * When are two alignments distinct and when are they redundant (not distinct)?
 * At one extreme, we might say the best alignment from any given dynamic
 * programming problem is redundant with all other alignments from that
 # problem.  At the other extreme, we might say that any two alignments with
 * distinct starting points and edits are distinct.  The former is probably too
 * conservative for mate-finding DP problems.  The latter is certainly too
 * permissive, since two alignments that differ only in how gaps are arranged
 * should not be considered distinct.
 *
 * Some in-between solutions are:
 *
 * (a) If two alignments share an end point on either end, they are redundant.
 *     Otherwise, they are distinct.
 * (b) If two alignments share *both* end points, they are redundant.
 * (c) If two alignments share any cells in the DP table, they are redundant.
 * (d) 2 alignments are redundant if either end within N poss of each other
 * (e) Like (d) but both instead of either
 * (f, g) Like d, e, but where N is tied to maxgaps somehow
 *
 * Why not (a)?  One reason is that it's possible for two alignments to have
 * different start & end positions but share many cells.  Consider alignments 1
 * and 2 below; their end-points are labeled.
 *
 *  1 2
 *  \ \
 *    -\
 *      \
 *       \
 *        \
 *        -\
 *        \ \
 *        1 2
 *
 * 1 and 2 are distinct according to (a) but they share many cells in common.
 *
 * Why not (f, g)?  It fixes the problem with (a) above by forcing the
 * alignments to be spread so far that they can't possibly share diagonal cells
 * in common
 */
class SwAligner {

	typedef std::pair<size_t, size_t> SizeTPair;

	// States that the aligner can be in
	enum {
		STATE_UNINIT,  // init() hasn't been called yet
		STATE_INITED,  // init() has been called, but not align()
		STATE_ALIGNED, // align() has been called
	};

public:

	SwAligner() :
		state_(STATE_UNINIT),
		inited_(false),
		rfwbuf_(DP_CAT),
		ntab_(DP_CAT),
		ctab_(DP_CAT),
		solcols_(DP_CAT),
		nfills_(0),
		ncups_(0),
		nrowups_(0),
		nrowskips_(0)
	{
		SwAligner::EXTREMES.first = std::numeric_limits<size_t>::max();
		SwAligner::EXTREMES.second = std::numeric_limits<size_t>::min();
	}

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
		TRefOff refoff,        // offset of upstream ref char aligned against
		char *rf,              // reference sequence
		size_t rfi,            // offset of first reference char to align to
		size_t rff,            // offset of last reference char to align to
		size_t width,          // # bands to do (width of parallelogram)
		EList<bool>* st,       // mask indicating which columns we can start in
		EList<bool>* en,       // mask indicating which columns we can end in
		const Scoring& sc,     // scoring scheme
		TAlScore minsc,        // minimum score a cell must achieve to have sol
		TAlScore floorsc)      // local-alignment score floor
	{
		size_t readGaps = sc.maxReadGaps(minsc, rd.length());
		size_t refGaps  = sc.maxRefGaps(minsc, rd.length());
		size_t maxGaps  = max(readGaps, refGaps);
		int nceil    = (int)sc.nCeil(rd.length());
		assert_geq(readGaps, 0);
		assert_geq(refGaps, 0);
		state_   = STATE_INITED;
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
		sc_      = &sc;        // scoring scheme
		minsc_   = minsc;      // minimum score a cell must achieve to have sol
		floorsc_ = floorsc;    // local-alignment score floor
		nceil_   = nceil;      // max # Ns allowed in ref portion of aln
		inited_  = true;       // indicate we're initialized now
		solrows_ = SwAligner::EXTREMES; // range of rows with at least 1 sol
		solcols_.clear();      // for each row, cells with >0 sols remaining
		INVALIDATE_SCORE(solbest_); // best score
		solrowbest_.clear();   // best score in each row
		solrowlo_ = sc.rowlo;  // if row >= this, solutions are possible
		soldone_  = true;      // true iff there are no more cells with sols
		nsols_    = 0;         // # cells with acceptable sols so far
		cural_    = 0;         // idx of next alignment to give out
		if(solrowlo_ < 0) {
			solrowlo_ = (int64_t)(dpRows()-1);
		}
		assert_geq(solrowlo_, 0);
		assert(en_ == NULL || en_->size() == width_);
		assert(st_ == NULL || st_->size() == width_);
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
		TRefOff rfi,           // off of first character in ref to consider
		TRefOff rff,           // off of last char (excl) in ref to consider
		const BitPairReference& refs, // Reference strings
		size_t reflen,         // length of reference sequence
		size_t width,          // # bands to do (width of parallelogram)
		EList<bool>* st,       // mask indicating which columns we can start in
		EList<bool>* en,       // mask indicating which columns we can end in
		const Scoring& sc,     // scoring scheme
		TAlScore minsc,        // minimum score a cell must achieve to have sol
		TAlScore floorsc,      // local-alignment score floor
		int nceil);            // max # Ns allowed in reference part of aln
	
	/**
	 * Align read 'rd' to reference using read & reference information given
	 * last time init() was called.  If the read is colorspace, the decoding is
	 * determined simultaneously with alignment.  Uses dynamic programming.
	 */
	bool align(RandomSource& rnd);
	
	/**
	 * Populate the given SwResult with information about the "next best"
	 * alignment if there is one.  If there isn't one, false is returned.  Note
	 * that false might be returned even though a call to done() would have
	 * returned false.
	 *
	 * Which alignment is "next best" depends on the 'rowfirst_' setting.  If
	 * it's true, then alignments in rows further down in the DP matrix get
	 * higher priority than rows higher up.  Within a row, alignments are
	 * ordered according to their score.  If 'rowfirst_' is false, then
	 * priority is given first according to score, then according to row.
	 */
	bool nextAlignment(
		SwResult& res,
		RandomSource& rnd);
	
	/**
	 * Print out an alignment result as an ASCII DP table.
	 */
	void printResultStacked(
		const SwResult& res,
		std::ostream& os)
	{
		res.alres.printStacked(*rd_, os);
	}
	
	/**
	 * Return true iff there are no more solution cells to backtace from.
	 * Note that this may return false in situations where there are actually
	 * no more solutions, but that hasn't been discovered yet.
	 */
	bool done() const {
		return soldone_;
	}

	/**
	 * Return true iff this SwAligner has been initialized with a dynamic
	 * programming problem.
	 */
	inline bool inited() const { return inited_; }
	
	/**
	 * Reset, signaling that we're done with this dynamic programming problem
	 * and won't be asking for any more alignments.
	 */
	inline void reset() { inited_ = false; }

	/**
	 * Do some quick filtering, setting elements of st_ and en_ to false
	 * according to the alignment policy and properties of the read/ref
	 * sequence.
	 */
	inline void filter(size_t nlim) { nfilter(nlim); }
	
	/**
	 * Check that aligner is internally consistent.
	 */
	bool repOk() const {
		assert_gt(dpRows(), 0);
		if(state_ == STATE_ALIGNED) {
			assert(soldone_ || solrows_.second >= solrows_.first);
			if(solrows_.second >= solrows_.first) {
				for(size_t i = solrows_.first; i <= solrows_.second; i++) {
					assert_leq(solrowbest_[i], solbest_);
					assert(!VALID_SCORE(solrowbest_[i]) || 
					       solcols_[i].second >= solcols_[i].first);
				}
			}
		}
		assert(st_ == NULL || st_->size() == width_);
		assert(en_ == NULL || en_->size() == width_);
		return true;
	}
	
	/**
	 * Return the number of alignments given out so far by nextAlignment().
	 */
	size_t numAlignmentsReported() const { return cural_; }

protected:
	
	/**
	 * Set elements of en_ to false if an ungapped alignment extending
	 * diagonally back from the corresponding cell in the last row would
	 * overlap too many Ns (more than nlim).
	 */
	size_t nfilter(size_t nlim);
	
	/**
	 * Return the number of rows that will be in the dynamic programming table.
	 */
	inline size_t dpRows() const {
		assert(inited_);
		return rdf_ - rdi_ + (color_ ? 1 : 0);
	}

	/**
	 * Align and simultaneously decode the colorspace read 'rd' to the
	 * reference string 'rf' using dynamic programming.  Return true iff
	 * zero or more alignments are possible.
	 */
	bool alignNucleotides(RandomSource& rnd);

	/**
	 * Align the nucleotide read 'rd' to the reference string 'rf' using
	 * dynamic programming.  Return true iff zero or more alignments are
	 * possible.
	 */
	bool alignColors(RandomSource& rnd);

	/**
	 * Given the dynamic programming table and a cell (both the table offset
	 * and the reference character), trace backwards from the cell and install
	 * the edits and score/penalty in the appropriate fields of res.  The
	 * RandomSource is used to break ties among equally good ways of tracing
	 * back.
	 *
	 * Note that the subject nucleotide sequence is decoded at the same time as
	 * the alignment is constructed.  So the traceback reveals both the
	 * nucleotide decoding (in ned) and the colorspace error pattern (in ced).
	 * The approach is very similar to the one described in the SHRiMP
	 * paper:
	 *
	 * Rumble SM, Lacroute P, Dalca AV, Fiume M, Sidow A, Brudno M. SHRiMP:
	 * accurate mapping of short color-space reads. PLoS Comput Biol. 2009
	 * May;5(5)
	 *
	 * Whenever we enter a cell, we check whether the read/ref coordinates of
	 * that cell correspond to a cell we traversed constructing a previous
	 * alignment.  If so, we backtrack to the last decision point, mask out the
	 * path that led to the previously observed cell, and continue along a
	 * different path; or, if there are no more paths to try, we give up.
	 *
	 * If an alignment is found, 'off' is set to the alignment's upstream-most
	 * reference character's offset into the chromosome and true is returned.
	 * Otherwise, false is returned.
	 */
	bool backtrackColors(
		TAlScore      dscore,  // score we expect to get over backtrack
		SwResult&     res,     // out: store results (edits and scores) here
		size_t&       off,     // out: leftmost ref char involved in aln
		size_t        row,     // start in this row (w/r/t full matrix)
		size_t        col,     // start in this column (w/r/t full matrix)
		int           lastC,   // cell to backtrace from in lower-right corner
		RandomSource& rand);   // pseudo-random generator

	/**
	 * Given the dynamic programming table and a cell, trace backwards from the
	 * cell and install the edits and score/penalty in the appropriate fields
	 * of res.  The RandomSource is used to break ties among equally good ways
	 * of tracing back.
	 *
	 * Whenever we enter a cell, we check whether the read/ref coordinates of
	 * that cell correspond to a cell we traversed constructing a previous
	 * alignment.  If so, we backtrack to the last decision point, mask out the
	 * path that led to the previously observed cell, and continue along a
	 * different path; or, if there are no more paths to try, we give up.
	 *
	 * If an alignment is found, 'off' is set to the alignment's upstream-most
	 * reference character's offset into the chromosome and true is returned.
	 * Otherwise, false is returned.
	 */
	bool backtrackNucleotides(
		TAlScore      escore,  // score we expect to get over backtrack
		SwResult&     res,     // out: store results (edits and scores) here
		size_t&       off,     // out: leftmost ref char involved in aln
		size_t        row,     // start in this row (w/r/t full matrix)
		size_t        col,     // start in this column (w/r/t full matrix)
		RandomSource& rand);   // pseudo-random generator

	/**
	 * Update the overall best alignment.
	 */
	void nextAlignmentUpdateBest();

	/**
	 * Given a range of columns in a row, update solrowbest_, solcols_, solbest_.
	 * Return true iff there areone or more solution cells in the row.
	 */
	bool nextAlignmentUpdateColorRow(size_t row);

	/**
	 * Given a range of columns in a row, update solrowbest_, solcols_, solbest_.
	 * Return true iff there areone or more solution cells in the row.
	 */
	bool nextAlignmentUpdateNucRow(size_t row);

	/**
	 * Try to report an alignment with the given score from the given row.
	 * After each attempt, update information about the row.
	 */
	bool nextAlignmentTryColorRow(
		SwResult& res,       // install result here
		size_t& off,         // ref offset of alignment
		size_t row,          // row to try a solution from
		TAlScore sc,         // potential solutions must have this score
		RandomSource& rnd);  // pseudo-random generator for backtrack

	/**
	 * Try to report an alignment with the given score from the given row.
	 * After each attempt, update information about the row.
	 */
	bool nextAlignmentTryNucRow(
		SwResult& res,       // install result here
		size_t& off,         // ref offset of alignment
		size_t row,          // row to try a solution from
		TAlScore sc,         // potential solutions must have this score
		RandomSource& rnd);  // pseudo-random generator for backtrack

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
		int pen,
		bool& improved);

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

	const BTDnaString  *rd_;     // read sequence
	const BTString     *qu_;     // read qualities
	size_t              rdi_;    // offset of first read char to align
	size_t              rdf_;    // offset of last read char to align
	bool                fw_;     // true iff read sequence is original fw read
	bool                color_;  // true iff read is colorspace
	uint32_t            refidx_; // id of reference aligned against
	TRefOff             refoff_; // offset of upstream ref char aligned against
	char               *rf_;     // reference sequence
	size_t              rfi_;    // offset of first ref char to align to
	size_t              rff_;    // offset of last ref char to align to (excl)
	size_t              width_;  // # bands to do (width of parallelogram)
	EList<bool>*        st_;     // mask indicating which cols we can start in
	EList<bool>*        en_;     // mask indicating which cols we can end in
	size_t              rdgap_;  // max # gaps in read
	size_t              rfgap_;  // max # gaps in reference
	int                 maxgap_; // max(rdgap_, rfgap_)
	const Scoring      *sc_;     // penalties for edit types
	TAlScore            minsc_;  // penalty ceiling for valid alignments
	TAlScore            floorsc_;// local-alignment score floor
	int                 nceil_;  // max # Ns allowed in ref portion of aln
	bool                monotone_; // true iff scores only go down

	int                 state_;  // state
	bool                inited_; // true iff initialized with DP problem
	EList<uint32_t>     rfwbuf_; // buffer for wordized refernece stretches
	ELList<SwNucCell>   ntab_;   // DP table for nucleotide read
	ELList<SwColorCell> ctab_;   // DP table for colorspace read
	
	EList<DpNucFrame>   btnstack_;// backtrace stack for nucleotides
	EList<DpColFrame>   btcstack_;// backtrace stack for colors
	EList<SizeTPair>    btcells_; // cells involved in current backtrace
	
	SizeTPair           solrows_; // range of rows with at least 1 sol
	EList<SizeTPair>    solcols_; // per row, cell range with >0 sols remaining
	TAlScore            solbest_; // best score
	EList<TAlScore>   solrowbest_;// best score in each row
	int64_t             solrowlo_;// if row >= this, solutions are possible
	bool                soldone_; // true iff there are no more cells with sols
	size_t              nsols_;   // # cells with acceptable sols so far
	
	size_t              cural_;   // index of next alignment to be given
	
	SizeTPair           EXTREMES; // invalid, uninitialized range
	
	// Counters
	uint64_t nfills_;    // table fills
	uint64_t ncups_;     // cell updates
	uint64_t nrowups_;   // row updates
	uint64_t nrowskips_; // row skips
};

#endif /*ALIGNER_SW_H_*/
