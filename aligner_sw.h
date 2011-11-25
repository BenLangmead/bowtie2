/*
 * Copyright 2011, Ben Langmead <blangmea@jhsph.edu>
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
 * aligner_sw.h
 *
 * Classes and routines for solving dynamic programming problems in aid of read
 * alignment.  Goals include the ability to handle:
 *
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
#include <limits>
#include <emmintrin.h>
#include "aligner_sw_common.h"
#include "aligner_sw_nuc.h"
#include "ds.h"
#include "threading.h"
#include "aligner_seed.h"
#include "reference.h"
#include "random_source.h"
#include "mem_ids.h"
#include "aligner_result.h"
#include "mask.h"
#include "dp_framer.h"
#include "aligner_swsse.h"

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
 * alignment of nucleotide reads against known reference nucleotides.
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

	explicit SwAligner() :
		sseU8fw_(DP_CAT),
		sseU8rc_(DP_CAT),
		sseI16fw_(DP_CAT),
		sseI16rc_(DP_CAT),
		sseI32fw_(DP_CAT),
		sseI32rc_(DP_CAT),
		state_(STATE_UNINIT),
		initedRead_(false),
		readSse16_(false),
		initedRef_(false),
		rfwbuf_(DP_CAT),
		btnstack_(DP_CAT),
		btcells_(DP_CAT),
		btncand_(DP_CAT),
		btncanddone_(DP_CAT),
		btncanddoneSucc_(0),
		btncanddoneFail_(0),
		colstop_(0),
		lastsolcol_(0),
		cural_(0)
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
		size_t rdi,              // offset of first read char to align
		size_t rdf,              // offset of last read char to align
		const Scoring& sc);      // scoring scheme
	
	/**
	 * Initialize with a new alignment problem.
	 */
	void initRef(
		bool fw,               // whether to forward or revcomp read is aligning
		uint32_t refidx,       // id of reference aligned against
		const DPRect& rect,    // DP rectangle
		char *rf,              // reference sequence
		size_t rfi,            // offset of first reference char to align to
		size_t rff,            // offset of last reference char to align to
		const Scoring& sc,     // scoring scheme
		TAlScore minsc,        // minimum score
		bool enable8,          // use 8-bit SSE if possible?
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
		const DPRect& rect,    // DP rectangle
		const BitPairReference& refs, // Reference strings
		size_t reflen,         // length of reference sequence
		const Scoring& sc,     // scoring scheme
		TAlScore minsc,        // minimum alignment score
		bool enable8,          // use 8-bit SSE if possible?
		bool extend,           // true iff this is a seed extension
		size_t  upto,          // count the number of Ns up to this offset
		size_t& nsUpto);       // output: the number of Ns up to 'upto'

	/**
	 * Given a read, an alignment orientation, a range of characters in a
	 * referece sequence, and a bit-encoded version of the reference, set up
	 * and execute the corresponding ungapped alignment problem.  There can
	 * only be one solution.
	 *
	 * The caller has already narrowed down the relevant portion of the
	 * reference using, e.g., the location of a seed hit, or the range of
	 * possible fragment lengths if we're searching for the opposite mate in a
	 * pair.
	 */
	int ungappedAlign(
		const BTDnaString&      rd,     // read sequence (could be RC)
		const BTString&         qu,     // qual sequence (could be rev)
		const Coord&            coord,  // coordinate aligned to
		const BitPairReference& refs,   // Reference strings
		size_t                  reflen, // length of reference sequence
		const Scoring&          sc,     // scoring scheme
		bool                    ohang,  // allow overhang?
		TAlScore                minsc,  // minimum score
		SwResult&               res);   // put alignment result here

	/**
	 * Align read 'rd' to reference using read & reference information given
	 * last time init() was called.  Uses dynamic programming.
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
		TAlScore minsc,
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
		return cural_ == btncand_.size();
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
	 * Check that aligner is internally consistent.
	 */
	bool repOk() const {
		assert_gt(dpRows(), 0);
		// Check btncand_
		for(size_t i = 0; i < btncand_.size(); i++) {
			assert(solrowlo_ < 0 || btncand_[i].row >= (size_t)solrowlo_);
			assert(btncand_[i].repOk());
			assert_geq(btncand_[i].score, minsc_);
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
		SSEMetrics& sseI16MateMet,
		uint64_t&   nbtfiltst,
		uint64_t&   nbtfiltsc,
		uint64_t&   nbtfiltdo)
	{
		sseU8ExtendMet.merge(sseU8ExtendMet_);
		sseU8MateMet.merge(sseU8MateMet_);
		sseI16ExtendMet.merge(sseI16ExtendMet_);
		sseI16MateMet.merge(sseI16MateMet_);
		nbtfiltst += nbtfiltst_;
		nbtfiltsc += nbtfiltsc_;
		nbtfiltdo += nbtfiltdo_;
	}
	
	/**
	 * Reset all the counters related to filling in the DP table to 0.
	 */
	void resetCounters() {
		sseU8ExtendMet_.reset();
		sseU8MateMet_.reset();
		sseI16ExtendMet_.reset();
		sseI16MateMet_.reset();
		nbtfiltst_ = nbtfiltsc_ = nbtfiltdo_ = 0;
	}
	
	/**
	 * Return the size of the DP problem.
	 */
	size_t size() const {
		return dpRows() * (rff_ - rfi_);
	}

protected:
	
	/**
	 * Return the number of rows that will be in the dynamic programming table.
	 */
	inline size_t dpRows() const {
		assert(initedRead_);
		return rdf_ - rdi_;
	}

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

	const BTDnaString  *rd_;     // read sequence
	const BTString     *qu_;     // read qualities
	const BTDnaString  *rdfw_;   // read sequence for fw read
	const BTDnaString  *rdrc_;   // read sequence for rc read
	const BTString     *qufw_;   // read qualities for fw read
	const BTString     *qurc_;   // read qualities for rc read
	size_t              rdi_;    // offset of first read char to align
	size_t              rdf_;    // offset of last read char to align
	bool                fw_;     // true iff read sequence is original fw read
	uint32_t            refidx_; // id of reference aligned against
	const DPRect*       rect_;   // DP rectangle
	char               *rf_;     // reference sequence
	size_t              rfi_;    // offset of first ref char to align to
	size_t              rff_;    // offset of last ref char to align to (excl)
	size_t              rdgap_;  // max # gaps in read
	size_t              rfgap_;  // max # gaps in reference
	bool                enable8_;// enable 8-bit sse
	bool                extend_; // true iff this is a seed-extend problem
	const Scoring      *sc_;     // penalties for edit types
	TAlScore            minsc_;  // penalty ceiling for valid alignments
	int                 nceil_;  // max # Ns allowed in ref portion of aln
	bool                monotone_; // true iff scores only go down

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

	SSEMetrics			sseU8ExtendMet_;
	SSEMetrics			sseU8MateMet_;
	SSEMetrics			sseI16ExtendMet_;
	SSEMetrics			sseI16MateMet_;

	int                 state_;      // state
	bool                initedRead_; // true iff initialized with initRead
	bool                readSse16_;  // true -> sse16 from now on for this read
	bool                initedRef_;  // true iff initialized with initRef
	EList<uint32_t>     rfwbuf_;     // buffer for wordized refernece stretches
	
	EList<DpNucFrame>   btnstack_;// backtrace stack for nucleotides
	EList<SizeTPair>    btcells_; // cells involved in current backtrace

	int64_t             solrowlo_;// if row >= this, solutions are possible
	EList<DpNucBtCandidate> btncand_;     // cells we might backtrace from
	EList<DpNucBtCandidate> btncanddone_; // candidates that we investigated
	size_t              btncanddoneSucc_; // # investigated and succeeded
	size_t              btncanddoneFail_; // # investigated and failed
	
	size_t              colstop_; // bailed on DP loop after this many cols
	size_t              lastsolcol_; // last DP col with valid cell
	size_t              cural_;   // index of next alignment to be given
	
	SizeTPair           EXTREMES; // invalid, uninitialized range
	
	uint64_t nbtfiltst_; // # candidates filtered b/c starting cell was seen
	uint64_t nbtfiltsc_; // # candidates filtered b/c score uninteresting
	uint64_t nbtfiltdo_; // # candidates filtered b/c dominated by other cell
	
	ASSERT_ONLY(SStringExpandable<uint32_t> tmp_destU32_);
	ASSERT_ONLY(BTDnaString tmp_editstr_, tmp_refstr_);
	ASSERT_ONLY(EList<DpNucBtCandidate> cand_tmp_);
};

#endif /*ALIGNER_SW_H_*/
