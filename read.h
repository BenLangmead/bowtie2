#ifndef READ_H_
#define READ_H_

/*
 *  read.h
 */

#include <stdint.h>
#include "ds.h"
#include "sstring.h"
#include "filebuf.h"

typedef uint64_t TReadId;

class HitSet;

/**
 * A buffer for keeping all relevant information about a single read.
 */
struct Read {
	Read() { reset(); }

	void reset() {
		patid = 0;
		alts = 0;
		trimmed5 = trimmed3 = 0;
		readOrigBuf.clear();
		qualOrigBuf.clear();
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
	}

	/**
	 * Simple init function, used for testing.
	 */
	void init(
		const char *nm,
		const char *seq,
		const char *ql,
		bool col)
	{
		reset();
		color = col;
		patFw.installChars(seq);
		qual.install(ql);
		constructReverses();
		constructRevComps();
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
	 * Check that read info is internally consistent.
	 */
	bool repOk() const {
		if(patFw.empty()) return true;
		assert_eq(qual.length(), patFw.length());
		return true;
	}

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
	// For when qualities are in a separate file
	SStringExpandable<char> qualOrigBuf;

	BTString name;      // read name
	TReadId  patid;     // unique 0-based id based on order in read file(s)
	int      mate;      // 0 = single-end, 1 = mate1, 2 = mate2
	uint32_t seed;      // random seed
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

#endif /*READ_H_*/
