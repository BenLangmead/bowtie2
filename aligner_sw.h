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

/**
 *  |-maxgaps-|
 *  ***********oooooooooooooooooooooo    -
 *   ***********ooooooooooooooooooooo    |
 *    ***********oooooooooooooooooooo    |
 *     ***********ooooooooooooooooooo    |
 *      ***********oooooooooooooooooo    |
 *       ***********ooooooooooooooooo read len
 *        ***********oooooooooooooooo    |
 *         ***********ooooooooooooooo    |
 *          ***********oooooooooooooo    |
 *           ***********ooooooooooooo    |
 *            ***********oooooooooooo    -
 *            |-maxgaps-|
 *  |-readlen-|
 *  |-------skip--------|
 */

#ifndef ALIGNER_SW_H_
#define ALIGNER_SW_H_

#define INLINE_CUPS

#include <stdint.h>
#include <iostream>

#ifndef NO_SSE
#include <emmintrin.h>
#endif

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
#include "seed_scan.h"

#ifndef NO_SSE
#include "aligner_swsse.h"
#endif

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
	
	const static size_t ALPHA_SIZE = 5;

public:

	explicit SwAligner(bool sse = true) :
#ifndef NO_SSE
		sse_(sse),
		sseU8fw_(DP_CAT),
		sseU8rc_(DP_CAT),
		sseI16fw_(DP_CAT),
		sseI16rc_(DP_CAT),
		sseI32fw_(DP_CAT),
		sseI32rc_(DP_CAT),
#endif
		state_(STATE_UNINIT),
		initedRead_(false),
		initedRef_(false),
		rfwbuf_(DP_CAT),
		ntab_(DP_CAT),
		ctab_(DP_CAT),
		btnstack_(DP_CAT),
		btcstack_(DP_CAT),
		btcells_(DP_CAT),
		btncand_(DP_CAT),
		btccand_(DP_CAT),
		nfills_(0),
		ncups_(0),
		nrowups_(0),
		nrowskips_(0),
		nskip_(0),
		nsucc_(0),
		nfail_(0),
		nbts_(0)
		ASSERT_ONLY(, cand_tmp_(DP_CAT))
	{
		SwAligner::EXTREMES.first = std::numeric_limits<size_t>::max();
		SwAligner::EXTREMES.second = std::numeric_limits<size_t>::min();
	}
	
	/**
	 * Prepare the dynamic programming driver with a new read and a new scoring
	 * scheme.
	 */
	void initRead(
		const BTDnaString& rdfw, // read sequence for fw read
		const BTDnaString& rdrc, // read sequence for rc read
		const BTString& qufw,    // read qualities for fw read
		const BTString& qurc,    // read qualities for rc read
		size_t rdi,            // offset of first read char to align
		size_t rdf,            // offset of last read char to align
		bool color,            // true iff read is colorspace
		const Scoring& sc,     // scoring scheme
		TAlScore minsc,        // minimum score a cell must achieve to have sol
		TAlScore floorsc);     // local-alignment score floor

	/**
	 * Initialize with a new alignment problem.
	 */
	void initRef(
		bool fw,               // whether to forward or revcomp read aligned
		uint32_t refidx,       // id of reference aligned against
		TRefOff refoff,        // offset of upstream ref char aligned against
		char *rf,              // reference sequence
		size_t rfi,            // offset of first reference char to align to
		size_t rff,            // offset of last reference char to align to
		size_t width,          // # bands to do (width of parallelogram)
		size_t solwidth,       // # rightmost cols where solns can end
		size_t maxgaps,        // max of max # read gaps, max # ref gaps
		size_t truncLeft,      // # cols/diags to truncate from LHS
		EList<bool>* en,       // mask indicating which columns we can end in
		bool extend);          // true iff this is a seed extension

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
	void initRef(
		bool fw,               // whether to forward or revcomp read aligned
		uint32_t refidx,       // reference aligned against
		TRefOff rfi,           // off of first character in ref to consider
		TRefOff rff,           // off of last char (excl) in ref to consider
		const BitPairReference& refs, // Reference strings
		size_t reflen,         // length of reference sequence
		size_t width,          // # bands to do (width of parallelogram)
		size_t solwidth,       // # rightmost cols where solns can end
		size_t maxgaps,        // max of max # read, ref gaps
		size_t truncLeft,      // columns to truncate from left-hand side of rect
		EList<bool>* en,       // mask indicating which columns we can end in
		bool extend,           // true iff this is a seed extension
		SeedScanner *sscan,    // optional seed scanner to feed ref chars to
		size_t  upto,          // count the number of Ns up to this offset
		size_t& nsUpto);       // output: the number of Ns up to 'upto'

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
		assert(initedRead() && initedRef());
		if(color_) {
			return cural_ == btccand_.size();
		} else {
			return cural_ == btncand_.size();
		}
	}

	/**
	 * Return true iff this SwAligner has been initialized with a read to align.
	 */
	inline bool initedRef() const { return initedRef_; }

	/**
	 * Return true iff this SwAligner has been initialized with a reference to
	 * align against.
	 */
	inline bool initedRead() const { return initedRead_; }
	
	/**
	 * Reset, signaling that we're done with this dynamic programming problem
	 * and won't be asking for any more alignments.
	 */
	inline void reset() { initedRef_ = initedRead_ = false; }

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
		//assert(st_ == NULL || st_->size() == solwidth_);
		assert(en_ == NULL || en_->size() == solwidth_);
		if(color_) {
			// Check btccand_
			for(size_t i = 0; i < btccand_.size(); i++) {
				assert_lt(btccand_[i].row, ctab_.size());
				assert(solrowlo_ < 0 || btccand_[i].row >= (size_t)solrowlo_);
				// The SSE aligner, when operating in local (not end-to-end)
				// mode might find solutions ending outside of the
				// parallelogram.
				//assert_lt(btccand_[i].col, ctab_[btccand_[i].row].size());
				assert(btccand_[i].repOk());
				assert_geq(btccand_[i].score, minsc_);
			}
		} else {
			// Check btncand_
			for(size_t i = 0; i < btncand_.size(); i++) {
				assert(sse_ || btncand_[i].row < ntab_.size());
				assert(solrowlo_ < 0 || btncand_[i].row >= (size_t)solrowlo_);
				// The SSE aligner, when operating in local (not end-to-end)
				// mode might find solutions ending outside of the
				// parallelogram.
				//assert_lt(btncand_[i].col, ntab_[btncand_[i].row].size());
				assert(btncand_[i].repOk());
				assert_geq(btncand_[i].score, minsc_);
			}
		}
		return true;
	}
	
	/**
	 * Return the number of alignments given out so far by nextAlignment().
	 */
	size_t numAlignmentsReported() const { return cural_; }

	/**
	 * Merge tallies in the counters related to filling the DP table.
	 */
	void merge(
		SSEMetrics& sseU8ExtendMet,
		SSEMetrics& sseU8MateMet,
		SSEMetrics& sseI16ExtendMet,
		SSEMetrics& sseI16MateMet)
	{
		sseU8ExtendMet.merge(sseU8ExtendMet_);
		sseU8MateMet.merge(sseU8MateMet_);
		sseI16ExtendMet.merge(sseI16ExtendMet_);
		sseI16MateMet.merge(sseI16MateMet_);
	}
	
	/**
	 * Reset all the counters related to filling in the DP table to 0.
	 */
	void resetCounters() {
		sseU8ExtendMet_.reset();
		sseU8MateMet_.reset();
		sseI16ExtendMet_.reset();
		sseI16MateMet_.reset();
	}

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
		assert(initedRead_);
		return rdf_ - rdi_ + (color_ ? 1 : 0);
	}

	/**
	 * Align nucleotides from read 'rd' to the reference string 'rf' using
	 * dynamic programming and normal, serial code.  Return true iff zero or
	 * more alignments are possible.
	 */
	TAlScore alignNucleotides();

#ifndef NO_SSE

	/**
	 * Align nucleotides from read 'rd' to the reference string 'rf' using
	 * vector instructions.  Return the score of the best alignment found, or
	 * the minimum integer if an alignment could not be found.  Flag is set to
	 * 0 if an alignment is found, -1 if no valid alignment is found, or -2 if
	 * the score saturated at any point during alignment.
	 */
	TAlScore alignNucleotidesEnd2EndSseU8(  // unsigned 8-bit elements
		int& flag);
	TAlScore alignNucleotidesLocalSseU8(    // unsigned 8-bit elements
		int& flag);
	TAlScore alignNucleotidesEnd2EndSseI16( // signed 16-bit elements
		int& flag);
	TAlScore alignNucleotidesLocalSseI16(   // signed 16-bit elements
		int& flag);
	
	/**
	 * Build query profile look up tables for the read.  The query profile look
	 * up table is organized as a 1D array indexed by [i][j] where i is the
	 * reference character in the current DP column (0=A, 1=C, etc), and j is
	 * the segment of the query we're currently working on.
	 */
	void buildQueryProfileEnd2EndSseU8(bool fw);
	void buildQueryProfileLocalSseU8(bool fw);

	/**
	 * Build query profile look up tables for the read.  The query profile look
	 * up table is organized as a 1D array indexed by [i][j] where i is the
	 * reference character in the current DP column (0=A, 1=C, etc), and j is
	 * the segment of the query we're currently working on.
	 */
	void buildQueryProfileEnd2EndSseI16(bool fw);
	void buildQueryProfileLocalSseI16(bool fw);
	
	bool gatherCellsNucleotidesLocalSseU8(TAlScore best);
	bool gatherCellsNucleotidesEnd2EndSseU8(TAlScore best);

	bool gatherCellsNucleotidesLocalSseI16(TAlScore best);
	bool gatherCellsNucleotidesEnd2EndSseI16(TAlScore best);

	bool backtraceNucleotidesLocalSseU8(
		TAlScore       escore, // in: expected score
		SwResult&      res,    // out: store results (edits and scores) here
		size_t&        off,    // out: store diagonal projection of origin
		size_t&        nbts,   // out: # backtracks
		size_t         row,    // start in this rectangle row
		size_t         col,    // start in this rectangle column
		RandomSource&  rand);  // random gen, to choose among equal paths

	bool backtraceNucleotidesLocalSseI16(
		TAlScore       escore, // in: expected score
		SwResult&      res,    // out: store results (edits and scores) here
		size_t&        off,    // out: store diagonal projection of origin
		size_t&        nbts,   // out: # backtracks
		size_t         row,    // start in this rectangle row
		size_t         col,    // start in this rectangle column
		RandomSource&  rand);  // random gen, to choose among equal paths

	bool backtraceNucleotidesEnd2EndSseU8(
		TAlScore       escore, // in: expected score
		SwResult&      res,    // out: store results (edits and scores) here
		size_t&        off,    // out: store diagonal projection of origin
		size_t&        nbts,   // out: # backtracks
		size_t         row,    // start in this rectangle row
		size_t         col,    // start in this rectangle column
		RandomSource&  rand);  // random gen, to choose among equal paths

	bool backtraceNucleotidesEnd2EndSseI16(
		TAlScore       escore, // in: expected score
		SwResult&      res,    // out: store results (edits and scores) here
		size_t&        off,    // out: store diagonal projection of origin
		size_t&        nbts,   // out: # backtracks
		size_t         row,    // start in this rectangle row
		size_t         col,    // start in this rectangle column
		RandomSource&  rand);  // random gen, to choose among equal paths

#endif

	/**
	 * Align the nucleotide read 'rd' to the reference string 'rf' using
	 * dynamic programming.  Return true iff zero or more alignments are
	 * possible.
	 */
	TAlScore alignColors();

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
	bool backtraceNucleotides(
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
		int prevQual,
		bool& improved);

	inline void updateColorVert(
		const SwColorCell& uc,
		SwColorCell& dstc,
		int prevColor,
		int prevQual);

	const BTDnaString  *rd_;     // read sequence
	const BTString     *qu_;     // read qualities
	const BTDnaString  *rdfw_;   // read sequence for fw read
	const BTDnaString  *rdrc_;   // read sequence for rc read
	const BTString     *qufw_;   // read qualities for fw read
	const BTString     *qurc_;   // read qualities for rc read
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
	size_t              solwidth_;// # bands ending @ exhaustively scored cells
	size_t              maxgaps_;// max of max # read gaps, max # ref gaps
	size_t              truncLeft_; // # cols/diags to truncate from LHS
	EList<bool>*        en_;     // mask indicating which cols we can end in
	size_t              rdgap_;  // max # gaps in read
	size_t              rfgap_;  // max # gaps in reference
	bool                extend_; // true iff this is a seed-extend problem
	const Scoring      *sc_;     // penalties for edit types
	TAlScore            minsc_;  // penalty ceiling for valid alignments
	TAlScore            floorsc_;// local-alignment score floor
	int                 nceil_;  // max # Ns allowed in ref portion of aln
	bool                monotone_; // true iff scores only go down

#ifndef NO_SSE
	bool                sse_;       // true -> use SSE 128-bit instructs
	bool                sse8succ_;  // whether 8-bit worked
	bool                sse16succ_; // whether 16-bit worked
	SSEData             sseU8fw_;   // buf for fw query, 8-bit score
	SSEData             sseU8rc_;   // buf for rc query, 8-bit score
	SSEData             sseI16fw_;  // buf for fw query, 16-bit score
	SSEData             sseI16rc_;  // buf for rc query, 16-bit score
	bool                sseU8fwBuilt_;   // built fw query profile, 8-bit score
	bool                sseU8rcBuilt_;   // built rc query profile, 8-bit score
	bool                sseI16fwBuilt_;  // built fw query profile, 16-bit score
	bool                sseI16rcBuilt_;  // built rc query profile, 16-bit score
	SSEData             sseI32fw_;  // buf for fw query, 32-bit score
	SSEData             sseI32rc_;  // buf for rc query, 32-bit score
#endif

	SSEMetrics			sseU8ExtendMet_;
	SSEMetrics			sseU8MateMet_;
	SSEMetrics			sseI16ExtendMet_;
	SSEMetrics			sseI16MateMet_;

	int                 state_;  // state
	bool                initedRead_; // true iff initialized with initRead
	bool                initedRef_;  // true iff initialized with initRef
	EList<uint32_t>     rfwbuf_; // buffer for wordized refernece stretches
	ELList<SwNucCell>   ntab_;   // DP table for nucleotide read
	ELList<SwColorCell> ctab_;   // DP table for colorspace read
	
	EList<DpNucFrame>   btnstack_;// backtrace stack for nucleotides
	EList<DpColFrame>   btcstack_;// backtrace stack for colors
	EList<SizeTPair>    btcells_; // cells involved in current backtrace

	int64_t             solrowlo_;// if row >= this, solutions are possible
	EList<DpNucBtCandidate> btncand_; // cells we might backtrace from
	EList<DpColBtCandidate> btccand_; // cells we might backtrace from
	
	size_t              cural_;   // index of next alignment to be given
	
	SizeTPair           EXTREMES; // invalid, uninitialized range
	
	// Holds potential solutions to backtrace from
	uint64_t nfills_;    // table fills
	uint64_t ncups_;     // cell updates
	uint64_t nrowups_;   // row updates
	uint64_t nrowskips_; // row skips
	uint64_t nskip_;     // # fills skipped b/c of SSE
	uint64_t nsucc_;     // # fills with at least 1 solution cell
	uint64_t nfail_;     // # fills with no solution cells
	uint64_t nbts_;      // backtrace steps
	
	ASSERT_ONLY(SStringExpandable<uint32_t> tmp_destU32_);
	ASSERT_ONLY(BTDnaString tmp_editstr_, tmp_refstr_);
	ASSERT_ONLY(EList<DpNucBtCandidate> cand_tmp_);
};

#endif /*ALIGNER_SW_H_*/
