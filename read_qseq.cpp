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

#include "pat.h"

/**
 * Parse a name from fb_ and store in r.  Assume that the next
 * character obtained via fb_.get() is the first character of
 * the sequence and the string stops at the next char upto (could
 * be tab, newline, etc.).
 */
static int parseName(
	Read::TBuf& buf, // buffer w/ raw qseq data
	size_t& cur,     // buffer cursor
	Read& r,         // buffer for mate 1
	int upto)        // stop parsing when we first reach character 'upto'
{
	const size_t buflen = buf.length();
	int c;
	while(cur < buflen) {
		c = buf[cur++];
		assert(c != '\r' && c != '\n');
		if(c == upto) {
			break; // Finished with field
		}
		r.name.append(c);
	}
	if(cur >= buflen) {
		return -1; // Error: buffer ended prematurely
	}
	return (int)r.name.length();
}


/**
 * Read another pattern from a Qseq input file.
 */
pair<bool, int> QseqPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a, unsigned readi)
{
	int c = getc_wrapper();
	while(c >= 0 && (c == '\n' || c == '\r')) {
		c = getc_wrapper();
	}
	EList<Read>& readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	// Read until we run out of input or until we've filled the buffer
	for(; readi < pt.max_buf_ && c >= 0; readi++) {
		readbuf[readi].readOrigBuf.clear();
		while(c >= 0 && c != '\n' && c != '\r') {
			readbuf[readi].readOrigBuf.append(c);
			c = getc_wrapper();
		}
		while(c >= 0 && (c == '\n' || c == '\r')) {
			c = getc_wrapper();
		}
	}
	if (c != EOF) {
		ungetc_wrapper(c);
	}
	return make_pair(c < 0, readi);
}

/**
 *
 */
bool QseqPatternSource::parse(Read& r, Read& rb, TReadId rdid) const {
	// Light parser (nextBatchFromFile) puts unparsed data
	// into Read& r, even when the read is paired.
	assert(r.empty());
	assert(!r.readOrigBuf.empty()); // raw data for read/pair is here
	int c = '\t';
	size_t cur = 0;
	const size_t buflen = r.readOrigBuf.length();
	assert(r.name.empty());
	
	// 1. Machine name
	if(parseName(r.readOrigBuf, cur, r, '\t') == -1) {
		return false;
	}
	r.name.append('_');
	// 2. Run number
	if(parseName(r.readOrigBuf, cur, r, '\t') == -1) {
		return false;
	}
	r.name.append('_');
	// 3. Lane number
	if(parseName(r.readOrigBuf, cur, r, '\t') == -1) {
		return false;
	}
	r.name.append('_');
	// 4. Tile number
	if(parseName(r.readOrigBuf, cur, r, '\t') == -1) {
		return false;
	}
	r.name.append('_');
	// 5. X coordinate of spot
	if(parseName(r.readOrigBuf, cur, r, '\t') == -1) {
		return false;
	}
	r.name.append('_');
	// 6. Y coordinate of spot
	if(parseName(r.readOrigBuf, cur, r, '\t') == -1) {
		return false;
	}
	r.name.append('_');
	// 7. Index
	if(parseName(r.readOrigBuf, cur, r, '\t') == -1) {
		return false;
	}
	r.name.append('/');
	// 8. Mate number
	if(parseName(r.readOrigBuf, cur, r, '\t') == -1) {
		return false;
	}
	if(cur >= buflen) {
		return false; // ended prematurely
	}
	c = r.readOrigBuf[cur++];
	assert(c != '\r' && c != '\n');
	// 9. Sequence & 10. Qualities
	if(c == '\t') {
		// empty sequence & qualities
		c = r.readOrigBuf[cur++];
		assert(c != '\r' && c != '\n');
		assert_eq('\t', c);
		cerr << "Warning: skipping empty QSEQ read with name '" << r.name << "'" << endl;
	} else {
		// 9. Sequence
		int nchar = 0;
		while(c != '\t') {
			if(c == '.') {
				c = 'N';
			}
			if(isalpha(c)) {
				assert_in(toupper(c), "ACGTN");
				if(++nchar > pp_.trim5) {
					assert_neq(0, asc2dnacat[c]);
					r.patFw.append(asc2dna[c]);
				}
			}
			if(cur >= buflen) {
				break;
			}
			c = r.readOrigBuf[cur++];
		}
		if(cur >= buflen) {
			return false; // ended prematurely
		}
		// record amt trimmed from 5' end due to --trim5
		r.trimmed5 = (int)(nchar - r.patFw.length());
		// record amt trimmed from 3' end due to --trim3
		r.trimmed3 = (int)(r.patFw.trimEnd(pp_.trim3));
		
		// 10. Qualities
		assert(r.qual.empty());
		int nqual = 0;
		if (pp_.intQuals) {
			int cur_int = 0;
			while(c != '\t') {
				cur_int *= 10;
				cur_int += (int)(c - '0');
				c = r.readOrigBuf[cur++];
				assert(c != '\r' && c != '\n');
				if(c == ' ' || c == '\t') {
					char cadd = intToPhred33(cur_int, pp_.solexa64);
					cur_int = 0;
					assert_geq(cadd, 33);
					if(++nqual > pp_.trim5) {
						r.qual.append(cadd);
					}
				}
			}
		} else {
			while(cur < buflen) {
				c = r.readOrigBuf[cur++];
				assert(c != '\r' && c != '\n');
				if (c == ' ') {
					wrongQualityFormat(r.name);
					return false;
				} else if(c == '\t') {
					break;
				}
				c = charToPhred33(c, pp_.solexa64, pp_.phred64);
				if(++nqual > r.trimmed5) {
					r.qual.append(c);
				}
			}
			r.qual.trimEnd(r.trimmed3);
			if(r.qual.length() < r.patFw.length()) {
				tooFewQualities(r.name);
				return false;
			} else if(r.qual.length() > r.patFw.length()) {
				tooManyQualities(r.name);
				return false;
			}
		}
	}
	assert_eq('\t', c);

	// 11. Filter flag
	if(cur >= buflen) {
		return false;
	}
	int filt = r.readOrigBuf[cur++];
	r.filter = filt;
	if(filt != '0' && filt != '1') {
		// Bad value for filt
		cerr << "Error: Bad value '" << filt
		     << "' for qseq filter flag" << endl;
		throw 1;
	}
	assert_eq(cur, buflen);
	r.parsed = true;
	if(!rb.parsed && !rb.readOrigBuf.empty()) {
		return parse(rb, r, rdid);
	}
	return true;
}
