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

#ifndef READ_H_
#define READ_H_

#include <stdint.h>
#include <sys/time.h>
#include "ds.h"
#include "sstring.h"
#include "filebuf.h"
#include "util.h"

typedef uint64_t TReadId;
typedef size_t TReadOff;
typedef int64_t TAlScore;

class HitSet;

/**
 * A buffer for keeping all relevant information about a single read.
 */
struct Read {

	Read() { reset(); }
	
	Read(const char *nm, const char *seq, const char *ql) { init(nm, seq, ql); }

	void reset() {
		rdid = 0;
		endid = 0;
		alts = 0;
		trimmed5 = trimmed3 = 0;
		readOrigBuf.clear();
		patFw.clear();
		patRc.clear();
		qual.clear();
		patFwRev.clear();
		patRcRev.clear();
		qualRev.clear();
		name.clear();
		for(int j = 0; j < 3; j++) {
			altPatFw[j].clear();
			altPatFwRev[j].clear();
			altPatRc[j].clear();
			altPatRcRev[j].clear();
			altQual[j].clear();
			altQualRev[j].clear();
		}
		color = fuzzy = false;
		primer = '?';
		trimc = '?';
		filter = '?';
		seed = 0;
		ns_ = 0;
	}
	
	/**
	 * Finish initializing a new read.
	 */
	void finalize() {
		for(size_t i = 0; i < patFw.length(); i++) {
			if((int)patFw[i] > 3) {
				ns_++;
			}
		}
		constructRevComps();
		constructReverses();
	}

	/**
	 * Simple init function, used for testing.
	 */
	void init(
		const char *nm,
		const char *seq,
		const char *ql)
	{
		reset();
		patFw.installChars(seq);
		qual.install(ql);
		for(size_t i = 0; i < patFw.length(); i++) {
			if((int)patFw[i] > 3) {
				ns_++;
			}
		}
		constructRevComps();
		constructReverses();
		if(nm != NULL) name.install(nm);
	}

	/// Return true iff the read (pair) is empty
	bool empty() const {
		return patFw.empty();
	}

	/// Return length of the read in the buffer
	size_t length() const {
		return patFw.length();
	}
	
	/**
	 * Return the number of Ns in the read.
	 */
	size_t ns() const {
		return ns_;
	}

	/**
	 * Construct reverse complement of the pattern and the fuzzy
	 * alternative patters.  If read is in colorspace, just reverse
	 * them.
	 */
	void constructRevComps() {
		if(color) {
			patRc.installReverse(patFw);
			for(int j = 0; j < alts; j++) {
				altPatRc[j].installReverse(altPatFw[j]);
			}
		} else {
			patRc.installReverseComp(patFw);
			for(int j = 0; j < alts; j++) {
				altPatRc[j].installReverseComp(altPatFw[j]);
			}
		}
	}

	/**
	 * Given patFw, patRc, and qual, construct the *Rev versions in
	 * place.  Assumes constructRevComps() was called previously.
	 */
	void constructReverses() {
		patFwRev.installReverse(patFw);
		patRcRev.installReverse(patRc);
		qualRev.installReverse(qual);
		for(int j = 0; j < alts; j++) {
			altPatFwRev[j].installReverse(altPatFw[j]);
			altPatRcRev[j].installReverse(altPatRc[j]);
			altQualRev[j].installReverse(altQual[j]);
		}
	}

	/**
	 * Append a "/1" or "/2" string onto the end of the name buf if
	 * it's not already there.
	 */
	void fixMateName(int i) {
		assert(i == 1 || i == 2);
		size_t namelen = name.length();
		bool append = false;
		if(namelen < 2) {
			// Name is too short to possibly have /1 or /2 on the end
			append = true;
		} else {
			if(i == 1) {
				// append = true iff mate name does not already end in /1
				append =
					name[namelen-2] != '/' ||
					name[namelen-1] != '1';
			} else {
				// append = true iff mate name does not already end in /2
				append =
					name[namelen-2] != '/' ||
					name[namelen-1] != '2';
			}
		}
		if(append) {
			name.append('/');
			name.append("012"[i]);
		}
	}

	/**
	 * Dump basic information about this read to the given ostream.
	 */
	void dump(std::ostream& os) const {
		using namespace std;
		os << name << ' ';
		if(color) {
			os << patFw.toZBufXForm("0123.");
		} else {
			os << patFw;
		}
		os << ' ';
		// Print out the fuzzy alternative sequences
		for(int j = 0; j < 3; j++) {
			bool started = false;
			if(!altQual[j].empty()) {
				for(size_t i = 0; i < length(); i++) {
					if(altQual[j][i] != '!') {
						started = true;
					}
					if(started) {
						if(altQual[j][i] == '!') {
							os << '-';
						} else {
							if(color) {
								os << "0123."[(int)altPatFw[j][i]];
							} else {
								os << altPatFw[j][i];
							}
						}
					}
				}
			}
			cout << " ";
		}
		os << qual.toZBuf() << " ";
		// Print out the fuzzy alternative quality strings
		for(int j = 0; j < 3; j++) {
			bool started = false;
			if(!altQual[j].empty()) {
				for(size_t i = 0; i < length(); i++) {
					if(altQual[j][i] != '!') {
						started = true;
					}
					if(started) {
						os << altQual[j][i];
					}
				}
			}
			if(j == 2) {
				os << endl;
			} else {
				os << " ";
			}
		}
	}
	
	/**
	 * Check whether two reads are the same in the sense that they will
	 * lead to us finding the same set of alignments.
	 */
	static bool same(
		const BTDnaString& seq1,
		const BTString&    qual1,
		const BTDnaString& seq2,
		const BTString&    qual2,
		bool qualitiesMatter)
	{
		if(seq1.length() != seq2.length()) {
			return false;
		}
		for(size_t i = 0; i < seq1.length(); i++) {
			if(seq1[i] != seq2[i]) return false;
		}
		if(qualitiesMatter) {
			if(qual1.length() != qual2.length()) {
				return false;
			}
			for(size_t i = 0; i < qual1.length(); i++) {
				if(qual1[i] != qual2[i]) return false;
			}
		}
		return true;
	}

	/**
	 * Get the nucleotide and quality value at the given offset from 5' end.
	 * If 'fw' is false, get the reverse complement.
	 */
	std::pair<int, int> get(TReadOff off5p, bool fw) const {
		assert_lt(off5p, length());
		int c = (int)patFw[off5p];
        int q = qual[off5p];
        assert_geq(q, 33);
		return make_pair((!fw && c < 4) ? (c ^ 3) : c, q - 33);
	}
	
	/**
	 * Get the nucleotide at the given offset from 5' end.
	 * If 'fw' is false, get the reverse complement.
	 */
	int getc(TReadOff off5p, bool fw) const {
		assert_lt(off5p, length());
		int c = (int)patFw[off5p];
		return (!fw && c < 4) ? (c ^ 3) : c;
	}
	
	/**
	 * Get the quality value at the given offset from 5' end.
	 */
	int getq(TReadOff off5p) const {
		assert_lt(off5p, length());
        int q = qual[off5p];
        assert_geq(q, 33);
		return q-33;
	}

#ifndef NDEBUG
	/**
	 * Check that read info is internally consistent.
	 */
	bool repOk() const {
		if(patFw.empty()) return true;
		assert_eq(qual.length(), patFw.length());
		return true;
	}
#endif

	BTDnaString patFw;            // forward-strand sequence
	BTDnaString patRc;            // reverse-complement sequence
	BTString    qual;             // quality values

	BTDnaString altPatFw[3];
	BTDnaString altPatRc[3];
	BTString    altQual[3];

	BTDnaString patFwRev;
	BTDnaString patRcRev;
	BTString    qualRev;

	BTDnaString altPatFwRev[3];
	BTDnaString altPatRcRev[3];
	BTString    altQualRev[3];

	// For remembering the exact input text used to define a read
	SStringExpandable<char> readOrigBuf;

	BTString name;      // read name
	TReadId  rdid;      // 0-based id based on pair's offset in read file(s)
	TReadId  endid;     // 0-based id based on pair's offset in read file(s)
	                    // and which mate ("end") this is
	int      mate;      // 0 = single-end, 1 = mate1, 2 = mate2
	uint32_t seed;      // random seed
	size_t   ns_;       // # Ns
	int      alts;      // number of alternatives
	bool     fuzzy;     // whether to employ fuzziness
	bool     color;     // whether read is in color space
	char     primer;    // primer base, for csfasta files
	char     trimc;     // trimmed color, for csfasta files
	char     filter;    // if read format permits filter char, set it here
	int      trimmed5;  // amount actually trimmed off 5' end
	int      trimmed3;  // amount actually trimmed off 3' end
	HitSet  *hitset;    // holds previously-found hits; for chaining
};

/**
 * A string of FmStringOps represent a string of tasks performed by the
 * best-first alignment search.  We model the search as a series of FM ops
 * interspersed with reported alignments.
 */
struct FmStringOp {
	bool alignment;  // true -> found an alignment
	TAlScore pen;    // penalty of the FM op or alignment
	size_t n;        // number of FM ops (only relevant for non-alignment)
};

/**
 * A string that summarizes the progress of an FM-index-assistet best-first
 * search.  Useful for trying to figure out what the aligner is spending its
 * time doing for a given read.
 */
struct FmString {

	/**
	 * Add one or more FM index ops to the op string
	 */
	void add(bool alignment, TAlScore pen, size_t nops) {
		if(ops.empty() || ops.back().pen != pen) {
			ops.expand();
			ops.back().alignment = alignment;
			ops.back().pen = pen;
			ops.back().n = 0;
		}
		ops.back().n++;
	}
	
	/**
	 * Reset FmString to uninitialized state.
	 */
	void reset() {
		pen = std::numeric_limits<TAlScore>::max();
		ops.clear();
	}

	/**
	 * Print a :Z optional field where certain characters (whitespace, colon
	 * and percent) are escaped using % escapes.
	 */
	void print(BTString& o, char *buf) const {
		for(size_t i = 0; i < ops.size(); i++) {
			if(i > 0) {
				o.append(';');
			}
			if(ops[i].alignment) {
				o.append("A,");
				itoa10(ops[i].pen, buf);
				o.append(buf);
			} else {
				o.append("F,");
				itoa10(ops[i].pen, buf); o.append(buf);
				o.append(',');
				itoa10(ops[i].n, buf); o.append(buf);
			}
		}
	}

	TAlScore pen;          // current penalty
	EList<FmStringOp> ops; // op string
};

/**
 * Key per-read metrics.  These are used for thresholds, allowing us to bail
 * for unproductive reads.  They also the basis of what's printed when the user
 * specifies --read-times.
 */
struct PerReadMetrics {

	PerReadMetrics() { reset(); }

	void reset() {
		nExIters =
		nExDps   = nExDpSuccs   = nExDpFails   =
		nMateDps = nMateDpSuccs = nMateDpFails =
		nExUgs   = nExUgSuccs   = nExUgFails   =
		nMateUgs = nMateUgSuccs = nMateUgFails =
		nExEes   = nExEeSuccs   = nExEeFails   =
		nRedundants =
		nEeFmops = nSdFmops = nExFmops =
		nDpFail = nDpFailStreak = nDpLastSucc =
		nUgFail = nUgFailStreak = nUgLastSucc =
		nEeFail = nEeFailStreak = nEeLastSucc =
		nFilt = 0;
		nFtabs = 0;
		nRedSkip = 0;
		nRedFail = 0;
		nRedIns = 0;
		doFmString = false;
		nSeedRanges = nSeedElts = 0;
		nSeedRangesFw = nSeedEltsFw = 0;
		nSeedRangesRc = nSeedEltsRc = 0;
		seedMedian = seedMean = 0;
		bestLtMinscMate1 =
		bestLtMinscMate2 = std::numeric_limits<TAlScore>::min();
		fmString.reset();
	}

	struct timeval  tv_beg; // timer start to measure how long alignment takes
	struct timezone tz_beg; // timer start to measure how long alignment takes

	uint64_t nExIters;      // iterations of seed hit extend loop

	uint64_t nExDps;        // # extend DPs run on this read
	uint64_t nExDpSuccs;    // # extend DPs run on this read
	uint64_t nExDpFails;    // # extend DPs run on this read
	
	uint64_t nExUgs;        // # extend ungapped alignments run on this read
	uint64_t nExUgSuccs;    // # extend ungapped alignments run on this read
	uint64_t nExUgFails;    // # extend ungapped alignments run on this read

	uint64_t nExEes;        // # extend ungapped alignments run on this read
	uint64_t nExEeSuccs;    // # extend ungapped alignments run on this read
	uint64_t nExEeFails;    // # extend ungapped alignments run on this read

	uint64_t nMateDps;      // # mate DPs run on this read
	uint64_t nMateDpSuccs;  // # mate DPs run on this read
	uint64_t nMateDpFails;  // # mate DPs run on this read
	
	uint64_t nMateUgs;      // # mate ungapped alignments run on this read
	uint64_t nMateUgSuccs;  // # mate ungapped alignments run on this read
	uint64_t nMateUgFails;  // # mate ungapped alignments run on this read

	uint64_t nRedundants;   // # redundant seed hits
	
	uint64_t nSeedRanges;   // # BW ranges found for seeds
	uint64_t nSeedElts;     // # BW elements found for seeds

	uint64_t nSeedRangesFw; // # BW ranges found for seeds from fw read
	uint64_t nSeedEltsFw;   // # BW elements found for seeds from fw read

	uint64_t nSeedRangesRc; // # BW ranges found for seeds from fw read
	uint64_t nSeedEltsRc;   // # BW elements found for seeds from fw read
	
	uint64_t seedMedian;    // median seed hit count
	uint64_t seedMean;      // rounded mean seed hit count
	
	uint64_t nEeFmops;      // FM Index ops for end-to-end alignment
	uint64_t nSdFmops;      // FM Index ops used to align seeds
	uint64_t nExFmops;      // FM Index ops used to resolve offsets
	
	uint64_t nFtabs;        // # ftab lookups
	uint64_t nRedSkip;      // # times redundant path was detected and aborted
	uint64_t nRedFail;      // # times a path was deemed non-redundant
	uint64_t nRedIns;       // # times a path was added to redundancy list
	
	uint64_t nDpFail;       // number of dp failures in a row up until now
	uint64_t nDpFailStreak; // longest streak of dp failures
	uint64_t nDpLastSucc;   // index of last dp attempt that succeeded
	
	uint64_t nUgFail;       // number of ungap failures in a row up until now
	uint64_t nUgFailStreak; // longest streak of ungap failures
	uint64_t nUgLastSucc;   // index of last ungap attempt that succeeded

	uint64_t nEeFail;       // number of ungap failures in a row up until now
	uint64_t nEeFailStreak; // longest streak of ungap failures
	uint64_t nEeLastSucc;   // index of last ungap attempt that succeeded
	
	uint64_t nFilt;         // # mates filtered
	
	TAlScore bestLtMinscMate1; // best invalid score observed for mate 1
	TAlScore bestLtMinscMate2; // best invalid score observed for mate 2
	
	// For collecting information to go into an FM string
	bool doFmString;
	FmString fmString;
};

#endif /*READ_H_*/
