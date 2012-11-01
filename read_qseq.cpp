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

#include "pat.h"

/**
 * Parse a name from fb_ and store in r.  Assume that the next
 * character obtained via fb_.get() is the first character of
 * the sequence and the string stops at the next char upto (could
 * be tab, newline, etc.).
 */
int QseqPatternSource::parseName(
	Read& r,      // buffer for mate 1
	Read* r2,     // buffer for mate 2 (NULL if mate2 is read separately)
	bool append,     // true -> append characters, false -> skip them
	bool clearFirst, // clear the name buffer first
	bool warnEmpty,  // emit a warning if nothing was added to the name
	bool useDefault, // if nothing is read, put readCnt_ as a default value
	int upto)        // stop parsing when we first reach character 'upto'
{
	if(clearFirst) {
		if(r2 != NULL) r2->name.clear();
		r.name.clear();
	}
	while(true) {
		int c;
		if((c = fb_.get()) < 0) {
			// EOF reached in the middle of the name
			return -1;
		}
		if(c == '\n' || c == '\r') {
			// EOL reached in the middle of the name
			return -1;
		}
		if(c == upto) {
			// Finished with field
			break;
		}
		if(append) {
			if(r2 != NULL) r2->name.append(c);
			r.name.append(c);
		}
	}
	// Set up a default name if one hasn't been set
	if(r.name.empty() && useDefault && append) {
		char cbuf[20];
		itoa10(readCnt_, cbuf);
		r.name.append(cbuf);
		if(r2 != NULL) r2->name.append(cbuf);
	}
	if(r.name.empty() && warnEmpty) {
		cerr << "Warning: read had an empty name field" << endl;
	}
	return (int)r.name.length();
}

/**
 * Parse a single sequence from fb_ and store in r.  Assume
 * that the next character obtained via fb_.get() is the first
 * character of the sequence and the sequence stops at the next
 * char upto (could be tab, newline, etc.).
 */
int QseqPatternSource::parseSeq(
	Read& r,
	int& charsRead,
	int& trim5,
	char upto)
{
	int begin = 0;
	int c = fb_.get();
	assert(c != upto);
	r.patFw.clear();
	r.color = gColor;
	if(gColor) {
		// NOTE: clearly this is not relevant for Illumina output, but
		// I'm keeping it here in case there's some reason to put SOLiD
		// data in this format in the future.
	
		// This may be a primer character.  If so, keep it in the
		// 'primer' field of the read buf and parse the rest of the
		// read without it.
		c = toupper(c);
		if(asc2dnacat[c] > 0) {
			// First char is a DNA char
			int c2 = toupper(fb_.peek());
			// Second char is a color char
			if(asc2colcat[c2] > 0) {
				r.primer = c;
				r.trimc = c2;
				trim5 += 2; // trim primer and first color
			}
		}
		if(c < 0) { return -1; }
	}
	while(c != upto) {
		if(c == '.') c = 'N';
		if(gColor) {
			if(c >= '0' && c <= '4') c = "ACGTN"[(int)c - '0'];
		}
		if(isalpha(c)) {
			assert_in(toupper(c), "ACGTN");
			if(begin++ >= trim5) {
				assert_neq(0, asc2dnacat[c]);
				r.patFw.append(asc2dna[c]);
			}
			charsRead++;
		}
		if((c = fb_.get()) < 0) {
			return -1;
		}
	}
	r.patFw.trimEnd(gTrim3);
	return (int)r.patFw.length();
}

/**
 * Parse a single quality string from fb_ and store in r.
 * Assume that the next character obtained via fb_.get() is
 * the first character of the quality string and the string stops
 * at the next char upto (could be tab, newline, etc.).
 */
int QseqPatternSource::parseQuals(
	Read& r,
	int charsRead,
	int dstLen,
	int trim5,
	char& c2,
	char upto = '\t',
	char upto2 = -1)
{
	int qualsRead = 0;
	int c = 0;
	if (intQuals_) {
		// Probably not relevant
		char buf[4096];
		while (qualsRead < charsRead) {
			qualToks_.clear();
			if(!tokenizeQualLine(fb_, buf, 4096, qualToks_)) break;
			for (unsigned int j = 0; j < qualToks_.size(); ++j) {
				char c = intToPhred33(atoi(qualToks_[j].c_str()), solQuals_);
				assert_geq(c, 33);
				if (qualsRead >= trim5) {
					r.qual.append(c);
				}
				++qualsRead;
			}
		} // done reading integer quality lines
		if (charsRead > qualsRead) tooFewQualities(r.name);
	} else {
		// Non-integer qualities
		while((qualsRead < dstLen + trim5) && c >= 0) {
			c = fb_.get();
			c2 = c;
			if (c == ' ') wrongQualityFormat(r.name);
			if(c < 0) {
				// EOF occurred in the middle of a read - abort
				return -1;
			}
			if(!isspace(c) && c != upto && (upto2 == -1 || c != upto2)) {
				if (qualsRead >= trim5) {
					c = charToPhred33(c, solQuals_, phred64Quals_);
					assert_geq(c, 33);
					r.qual.append(c);
				}
				qualsRead++;
			} else {
				break;
			}
		}
	}
	if(r.qual.length() < (size_t)dstLen) {
		tooFewQualities(r.name);
	}
	// TODO: How to detect too many qualities??
	r.qual.resize(dstLen);
	while(c != -1 && c != upto && (upto2 == -1 || c != upto2)) {
		c = fb_.get();
		c2 = c;
	}
	return qualsRead;
}

/**
 * Read another pattern from a Qseq input file.
 */
bool QseqPatternSource::read(
	Read& r,
	TReadId& rdid,
	TReadId& endid,
	bool& success,
	bool& done)
{
	r.reset();
	r.color = gColor;
	success = true;
	done = false;
	readCnt_++;
	rdid = endid = readCnt_-1;
	peekOverNewline(fb_);
	fb_.resetLastN();
	// 1. Machine name
	if(parseName(r, NULL, true, true,  true, false, '\t') == -1) BAIL_UNPAIRED();
	assert_neq('\t', fb_.peek());
	r.name.append('_');
	// 2. Run number
	if(parseName(r, NULL, true, false, true, false, '\t') == -1) BAIL_UNPAIRED();
	assert_neq('\t', fb_.peek());
	r.name.append('_');
	// 3. Lane number
	if(parseName(r, NULL, true, false, true, false, '\t') == -1) BAIL_UNPAIRED();
	assert_neq('\t', fb_.peek());
	r.name.append('_');
	// 4. Tile number
	if(parseName(r, NULL, true, false, true, false, '\t') == -1) BAIL_UNPAIRED();
	assert_neq('\t', fb_.peek());
	r.name.append('_');
	// 5. X coordinate of spot
	if(parseName(r, NULL, true, false, true, false, '\t') == -1) BAIL_UNPAIRED();
	assert_neq('\t', fb_.peek());
	r.name.append('_');
	// 6. Y coordinate of spot
	if(parseName(r, NULL, true, false, true, false, '\t') == -1) BAIL_UNPAIRED();
	assert_neq('\t', fb_.peek());
	r.name.append('_');
	// 7. Index
	if(parseName(r, NULL, true, false, true, false, '\t') == -1) BAIL_UNPAIRED();
	assert_neq('\t', fb_.peek());
	r.name.append('/');
	// 8. Mate number
	if(parseName(r, NULL, true, false, true, false, '\t') == -1) BAIL_UNPAIRED();
	// Empty sequence??
	if(fb_.peek() == '\t') {
		// Get tab that separates seq from qual
		ASSERT_ONLY(int c =) fb_.get();
		assert_eq('\t', c);
		assert_eq('\t', fb_.peek());
		// Get tab that separates qual from filter
		ASSERT_ONLY(c =) fb_.get();
		assert_eq('\t', c);
		// Next char is first char of filter flag
		assert_neq('\t', fb_.peek());
		fb_.resetLastN();
		cerr << "Warning: skipping empty QSEQ read with name '" << r.name << "'" << endl;
	} else {
		assert_neq('\t', fb_.peek());
		int charsRead = 0;
		int mytrim5 = gTrim5;
		// 9. Sequence
		int dstLen = parseSeq(r, charsRead, mytrim5, '\t');
		assert_neq('\t', fb_.peek());
		if(dstLen < 0) BAIL_UNPAIRED();
		char ct = 0;
		// 10. Qualities
		if(parseQuals(r, charsRead, dstLen, mytrim5, ct, '\t', -1) < 0) BAIL_UNPAIRED();
		r.trimmed3 = gTrim3;
		r.trimmed5 = mytrim5;
		if(ct != '\t') {
			cerr << "Error: QSEQ with name " << r.name << " did not have tab after qualities" << endl;
			throw 1;
		}
		assert_eq(ct, '\t');
	}
	// 11. Filter flag
	int filt = fb_.get();
	if(filt == -1) BAIL_UNPAIRED();
	r.filter = filt;
	if(filt != '0' && filt != '1') {
		// Bad value for filt
	}
	if(fb_.peek() != -1 && fb_.peek() != '\n') {
		// Bad value right after the filt field
	}
	fb_.get();
	r.readOrigBuf.install(fb_.lastN(), fb_.lastNLen());
	fb_.resetLastN();
	if(r.qual.length() < r.patFw.length()) {
		tooFewQualities(r.name);
	} else if(r.qual.length() > r.patFw.length()) {
		tooManyQualities(r.name);
	}
#ifndef NDEBUG
	assert_eq(r.patFw.length(), r.qual.length());
	for(size_t i = 0; i < r.qual.length(); i++) {
		assert_geq((int)r.qual[i], 33);
	}
#endif
	return true;
}
