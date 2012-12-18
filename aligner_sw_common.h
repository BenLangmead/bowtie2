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

#ifndef ALIGNER_SW_COMMON_H_
#define ALIGNER_SW_COMMON_H_

#include "aligner_result.h"

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
		swskip(0),
		swsucc(0),
		swfail(0),
		swbts(0)
	{ }

	/**
	 * Clear all contents.
	 */
	void reset() {
		sws = swcups = swrows = swskiprows = swskip = swsucc =
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
	
#ifndef NDEBUG
	/**
	 * Check that result is internally consistent.
	 */
	bool repOk() const {
		assert(alres.repOk());
		return true;
	}

	/**
	 * Check that result is internally consistent w/r/t read.
	 */
	bool repOk(const Read& rd) const {
		assert(alres.repOk(rd));
		return true;
	}
#endif

	AlnRes alres;
	uint64_t sws;    // # DP problems solved
	uint64_t swcups; // # DP cell updates
	uint64_t swrows; // # DP row updates
	uint64_t swskiprows; // # skipped DP row updates (b/c no valid alignments can go thru row)
	uint64_t swskip; // # DP problems skipped by sse filter
	uint64_t swsucc; // # DP problems resulting in alignment
	uint64_t swfail; // # DP problems not resulting in alignment
	uint64_t swbts;  // # DP backtrace steps
	
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
		sws = swcups = swrows = swskiprows = swskip = swsucc = swfail = swbts =
		sws10 = sws5 = sws3 =
		rshit = ungapsucc = ungapfail = ungapnodec = 0;
		exatts = exranges = exrows = exsucc = exooms = 0;
		mm1atts = mm1ranges = mm1rows = mm1succ = mm1ooms = 0;
		sdatts = sdranges = sdrows = sdsucc = sdooms = 0;
	}
	
	void init(
		uint64_t sws_,
		uint64_t sws10_,
		uint64_t sws5_,
		uint64_t sws3_,
		uint64_t swcups_,
		uint64_t swrows_,
		uint64_t swskiprows_,
		uint64_t swskip_,
		uint64_t swsucc_,
		uint64_t swfail_,
		uint64_t swbts_,
		uint64_t rshit_,
		uint64_t ungapsucc_,
		uint64_t ungapfail_,
		uint64_t ungapnodec_,
		uint64_t exatts_,
		uint64_t exranges_,
		uint64_t exrows_,
		uint64_t exsucc_,
		uint64_t exooms_,
		uint64_t mm1atts_,
		uint64_t mm1ranges_,
		uint64_t mm1rows_,
		uint64_t mm1succ_,
		uint64_t mm1ooms_,
		uint64_t sdatts_,
		uint64_t sdranges_,
		uint64_t sdrows_,
		uint64_t sdsucc_,
		uint64_t sdooms_)
	{
		sws        = sws_;
		sws10      = sws10_;
		sws5       = sws5_;
		sws3       = sws3_;
		swcups     = swcups_;
		swrows     = swrows_;
		swskiprows = swskiprows_;
		swskip     = swskip_;
		swsucc     = swsucc_;
		swfail     = swfail_;
		swbts      = swbts_;
		ungapsucc  = ungapsucc_;
		ungapfail  = ungapfail_;
		ungapnodec = ungapnodec_;
		
		// Exact end-to-end attempts
		exatts     = exatts_;
		exranges   = exranges_;
		exrows     = exrows_;
		exsucc     = exsucc_;
		exooms     = exooms_;

		// 1-mismatch end-to-end attempts
		mm1atts    = mm1atts_;
		mm1ranges  = mm1ranges_;
		mm1rows    = mm1rows_;
		mm1succ    = mm1succ_;
		mm1ooms    = mm1ooms_;
		
		// Seed attempts
		sdatts     = sdatts_;
		sdranges   = sdranges_;
		sdrows     = sdrows_;
		sdsucc     = sdsucc_;
		sdooms     = sdooms_;
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
		swskip     += r.swskip;
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
		sws10      += r.sws10;
		sws5       += r.sws5;
		sws3       += r.sws3;
		swcups     += r.swcups;
		swrows     += r.swrows;
		swskiprows += r.swskiprows;
		swskip     += r.swskip;
		swsucc     += r.swsucc;
		swfail     += r.swfail;
		swbts      += r.swbts;
		rshit      += r.rshit;
		ungapsucc  += r.ungapsucc;
		ungapfail  += r.ungapfail;
		ungapnodec += r.ungapnodec;
		exatts     += r.exatts;
		exranges   += r.exranges;
		exrows     += r.exrows;
		exsucc     += r.exsucc;
		exooms     += r.exooms;
		mm1atts    += r.mm1atts;
		mm1ranges  += r.mm1ranges;
		mm1rows    += r.mm1rows;
		mm1succ    += r.mm1succ;
		mm1ooms    += r.mm1ooms;
		sdatts     += r.sdatts;
		sdranges   += r.sdranges;
		sdrows     += r.sdrows;
		sdsucc     += r.sdsucc;
		sdooms     += r.sdooms;
	}
	
	void tallyGappedDp(size_t readGaps, size_t refGaps) {
		size_t mx = max(readGaps, refGaps);
		if(mx < 10) sws10++;
		if(mx < 5)  sws5++;
		if(mx < 3)  sws3++;
	}

	uint64_t sws;        // # DP problems solved
	uint64_t sws10;      // # DP problems solved where max gaps < 10
	uint64_t sws5;       // # DP problems solved where max gaps < 5
	uint64_t sws3;       // # DP problems solved where max gaps < 3
	uint64_t swcups;     // # DP cell updates
	uint64_t swrows;     // # DP row updates
	uint64_t swskiprows; // # skipped DP rows (b/c no valid alns go thru row)
	uint64_t swskip;     // # DP problems skipped by sse filter
	uint64_t swsucc;     // # DP problems resulting in alignment
	uint64_t swfail;     // # DP problems not resulting in alignment
	uint64_t swbts;      // # DP backtrace steps
	uint64_t rshit;      // # DP problems avoided b/c seed hit was redundant
	uint64_t ungapsucc;  // # DP problems avoided b/c seed hit was redundant
	uint64_t ungapfail;  // # DP problems avoided b/c seed hit was redundant
	uint64_t ungapnodec; // # DP problems avoided b/c seed hit was redundant

	uint64_t exatts;     // total # attempts at exact-hit end-to-end aln
	uint64_t exranges;   // total # ranges returned by exact-hit queries
	uint64_t exrows;     // total # rows returned by exact-hit queries
	uint64_t exsucc;     // exact-hit yielded non-empty result
	uint64_t exooms;     // exact-hit offset memory exhausted
	
	uint64_t mm1atts;    // total # attempts at 1mm end-to-end aln
	uint64_t mm1ranges;  // total # ranges returned by 1mm-hit queries
	uint64_t mm1rows;    // total # rows returned by 1mm-hit queries
	uint64_t mm1succ;    // 1mm-hit yielded non-empty result
	uint64_t mm1ooms;    // 1mm-hit offset memory exhausted

	uint64_t sdatts;     // total # attempts to find seed alignments
	uint64_t sdranges;   // total # seed-alignment ranges found
	uint64_t sdrows;     // total # seed-alignment rows found
	uint64_t sdsucc;     // # times seed alignment yielded >= 1 hit
	uint64_t sdooms;     // # times an OOM occurred during seed alignment
	
	MUTEX_T lock;
};

// The various ways that one might backtrack from a later cell (either oall,
// rdgap or rfgap) to an earlier cell
enum {
	SW_BT_OALL_DIAG,         // from oall cell to oall cell
	SW_BT_OALL_REF_OPEN,     // from oall cell to oall cell
	SW_BT_OALL_READ_OPEN,    // from oall cell to oall cell
	SW_BT_RDGAP_EXTEND,      // from rdgap cell to rdgap cell
	SW_BT_RFGAP_EXTEND       // from rfgap cell to rfgap cell
};

#endif /*def ALIGNER_SW_COMMON_H_*/
