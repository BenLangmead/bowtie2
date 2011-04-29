#include <cmath>
#include <iostream>
#include <string>
#include <stdexcept>
#include "sstring.h"

#include "pat.h"
#include "filebuf.h"
#include "formats.h"

using namespace std;

/**
 * Return a new dynamically allocated PatternSource for the given
 * format, using the given list of strings as the filenames to read
 * from or as the sequences themselves (i.e. if -c was used).
 */
PatternSource* PatternSource::patsrcFromStrings(
	const PatternParams& p,
	const EList<string>& qs,
	const EList<string>* qualities)
{
	switch(p.format) {
		case FASTA:       return new FastaPatternSource(qs, qualities, p);
		case FASTA_CONT:  return new FastaContinuousPatternSource(qs, p);
		case RAW:         return new RawPatternSource(qs, p);
		case FASTQ:       return new FastqPatternSource(qs, p);
		case TAB_MATE:    return new TabbedPatternSource(qs, p);
		case CMDLINE:     return new VectorPatternSource(qs, p);
		case INPUT_CHAIN: return new ChainPatternSource(qs, p);
		case RANDOM:      return new RandomPatternSource(p);
		case QSEQ:        return new QseqPatternSource(qs, p);
		default: {
			cerr << "Internal error; bad patsrc format: " << p.format << endl;
			throw 1;
		}
	}
}

/**
 * Given the values for all of the various arguments used to specify
 * the read and quality input, create a list of pattern sources to
 * dispense them.
 */
PairedPatternSource* PairedPatternSource::setupPatternSources(
	const EList<string>& si,   // singles, from argv
	const EList<string>& m1,   // mate1's, from -1 arg
	const EList<string>& m2,   // mate2's, from -2 arg
	const EList<string>& m12,  // both mates on each line, from --12 arg
	const EList<string>& q,    // qualities associated with singles
	const EList<string>& q1,   // qualities associated with m1
	const EList<string>& q2,   // qualities associated with m2
	const PatternParams& p,    // read-in parameters
	bool verbose)              // be talkative?
{
	EList<PatternSource*>* a  = new EList<PatternSource*>();
	EList<PatternSource*>* b  = new EList<PatternSource*>();
	EList<PatternSource*>* ab = new EList<PatternSource*>();
	// Create list of pattern sources for paired reads appearing
	// interleaved in a single file
	for(size_t i = 0; i < m12.size(); i++) {
		const EList<string>* qs = &m12;
		EList<string> tmp;
		if(p.fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmp;
			tmp.push_back(m12[i]);
			assert_eq(1, tmp.size());
		}
		ab->push_back(PatternSource::patsrcFromStrings(p, *qs, NULL));
		if(!p.fileParallel) {
			break;
		}
	}

	// Create list of pattern sources for paired reads
	for(size_t i = 0; i < m1.size(); i++) {
		const EList<string>* qs = &m1;
		const EList<string>* quals = &q1;
		EList<string> tmpSeq;
		EList<string> tmpQual;
		if(p.fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmpSeq;
			tmpSeq.push_back(m1[i]);
			quals = &tmpSeq;
			tmpQual.push_back(q1[i]);
			assert_eq(1, tmpSeq.size());
		}
		if(quals->empty()) quals = NULL;
		a->push_back(PatternSource::patsrcFromStrings(p, *qs, quals));
		if(!p.fileParallel) {
			break;
		}
	}

	// Create list of pattern sources for paired reads
	for(size_t i = 0; i < m2.size(); i++) {
		const EList<string>* qs = &m2;
		const EList<string>* quals = &q2;
		EList<string> tmpSeq;
		EList<string> tmpQual;
		if(p.fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmpSeq;
			tmpSeq.push_back(m2[i]);
			quals = &tmpQual;
			tmpQual.push_back(q2[i]);
			assert_eq(1, tmpSeq.size());
		}
		if(quals->empty()) quals = NULL;
		b->push_back(PatternSource::patsrcFromStrings(p, *qs, quals));
		if(!p.fileParallel) {
			break;
		}
	}
	// All mates/mate files must be paired
	assert_eq(a->size(), b->size());

	// Create list of pattern sources for the unpaired reads
	for(size_t i = 0; i < si.size(); i++) {
		const EList<string>* qs = &si;
		const EList<string>* quals = &q;
		PatternSource* patsrc = NULL;
		EList<string> tmpSeq;
		EList<string> tmpQual;
		if(p.fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmpSeq;
			tmpSeq.push_back(si[i]);
			quals = &tmpQual;
			tmpQual.push_back(q[i]);
			assert_eq(1, tmpSeq.size());
		}
		if(quals->empty()) quals = NULL;
		patsrc = PatternSource::patsrcFromStrings(p, *qs, quals);
		assert(patsrc != NULL);
		a->push_back(patsrc);
		b->push_back(NULL);
		if(!p.fileParallel) {
			break;
		}
	}

	PairedPatternSource *patsrc = NULL;
	if(m12.size() > 0) {
		patsrc = new PairedSoloPatternSource(ab, p);
		for(size_t i = 0; i < a->size(); i++) delete (*a)[i];
		for(size_t i = 0; i < b->size(); i++) delete (*b)[i];
		delete a; delete b;
	} else {
		patsrc = new PairedDualPatternSource(a, b, p);
		for(size_t i = 0; i < ab->size(); i++) delete (*ab)[i];
		delete ab;
	}
	return patsrc;
}

/**
 * Parse a single quality string from fb and store qualities in r.
 * Assume the next character obtained via fb.get() is the first
 * character of the quality string.  When returning, the next
 * character returned by fb.peek() or fb.get() should be the first
 * character of the following line.
 */
int parseQuals(Read& r,
               FileBuf& fb,
               int firstc,
               int readLen,
               int trim3,
               int trim5,
               bool intQuals,
               bool phred64,
               bool solexa64)
{
	int c = firstc;
	assert(c != '\n' && c != '\r');
	r.qual.clear();
	if (intQuals) {
		while (c != '\r' && c != '\n' && c != -1) {
			bool neg = false;
			int num = 0;
			while(!isspace(c) && !fb.eof()) {
				if(c == '-') {
					neg = true;
					assert_eq(num, 0);
				} else {
					if(!isdigit(c)) {
						char buf[2048];
						cerr << "Warning: could not parse quality line:" << endl;
						fb.getPastNewline();
						cerr << fb.copyLastN(buf);
						buf[2047] = '\0';
						cerr << buf;
						throw 1;
					}
					assert(isdigit(c));
					num *= 10;
					num += (c - '0');
				}
				c = fb.get();
			}
			if(neg) num = 0;
			// Phred-33 ASCII encode it and add it to the back of the
			// quality string
			r.qual.append('!' + num);
			// Skip over next stretch of whitespace
			while(c != '\r' && c != '\n' && isspace(c) && !fb.eof()) {
				c = fb.get();
			}
		}
	} else {
		while (c != '\r' && c != '\n' && c != -1) {
			r.qual.append(charToPhred33(c, solexa64, phred64));
			c = fb.get();
			while(c != '\r' && c != '\n' && isspace(c) && !fb.eof()) {
				c = fb.get();
			}
		}
	}
	if ((int)r.qual.length() < readLen-1 ||
	    ((int)r.qual.length() < readLen && !r.color))
	{
		tooFewQualities(r.name);
	}
	r.qual.trimEnd(trim3);
	if(r.qual.length()-trim5 < r.patFw.length()) {
		assert(gColor && r.primer != -1);
		assert_gt(trim5, 0);
		trim5--;
	}
	r.qual.trimBegin(trim5);
	if(r.qual.length() <= 0) return 0;
	assert_eq(r.qual.length(), r.patFw.length());
	while(fb.peek() == '\n' || fb.peek() == '\r') fb.get();
	return (int)r.qual.length();
}

/// Read another pattern from a FASTQ input file
void FastqPatternSource::read(Read& r, TReadId& patid) {
	int c;
	int dstLen = 0;
	r.reset();
	r.color = gColor;
	r.fuzzy = fuzzy_;
	// Pick off the first at
	if(first_) {
		c = fb_.get();
		if(c != '@') {
			c = getOverNewline(fb_);
			if(c < 0) { bail(r); return; }
		}
		if(c != '@') {
			cerr << "Error: reads file does not look like a FASTQ file" << endl;
			throw 1;
		}
		assert_eq('@', c);
		first_ = false;
	}

	// Read to the end of the id line, sticking everything after the '@'
	// into *name
	while(true) {
		c = fb_.get();
		if(c < 0) { bail(r); return; }
		if(c == '\n' || c == '\r') {
			// Break at end of line, after consuming all \r's, \n's
			while(c == '\n' || c == '\r') {
				c = fb_.get();
				if(c < 0) { bail(r); return; }
			}
			break;
		}
		r.name.append(c);
	}
	// fb_ now points just past the first character of a
	// sequence line, and c holds the first character
	int charsRead = 0;
	BTDnaString *sbuf = &r.patFw;
	int dstLens[] = {0, 0, 0, 0};
	int *dstLenCur = &dstLens[0];
	int mytrim5 = gTrim5;
	int altBufIdx = 0;
	if(gColor && c != '+') {
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
				mytrim5 += 2; // trim primer and first color
			}
		}
		if(c < 0) { bail(r); return; }
	}
	int trim5 = 0;
	if(c != '+') {
		trim5 = mytrim5;
		while(c != '+') {
			// Convert color numbers to letters if necessary
			if(gColor) {
				if(c >= '0' && c <= '4') c = "ACGTN"[(int)c - '0'];
				if(c == '.') c = 'N';
			}
			if(fuzzy_ && c == '-') c = 'A';
			if(isalpha(c)) {
				// If it's past the 5'-end trim point
				if(charsRead >= trim5) {
					sbuf->append(asc2dna[c]);
					(*dstLenCur)++;
				}
				charsRead++;
			} else if(fuzzy_ && c == ' ') {
				trim5 = 0; // disable 5' trimming for now
				if(charsRead == 0) {
					c = fb_.get();
					continue;
				}
				charsRead = 0;
				if(altBufIdx >= 3) {
					cerr << "At most 3 alternate sequence strings permitted; offending read: " << r.name << endl;
					throw 1;
				}
				// Move on to the next alternate-sequence buffer
				sbuf = &r.altPatFw[altBufIdx++];
				dstLenCur = &dstLens[altBufIdx];
			}
			c = fb_.get();
			if(c < 0) { bail(r); return; }
		}
		dstLen = dstLens[0];
		charsRead = dstLen + mytrim5;
	}
	// Trim from 3' end
	if(gTrim3 > 0) {
		if((int)r.patFw.length() > gTrim3) {
			r.patFw.resize(r.patFw.length() - gTrim3);
			dstLen -= gTrim3;
			assert_eq((int)r.patFw.length(), dstLen);
		} else {
			// Trimmed the whole read; we won't be using this read,
			// but we proceed anyway so that fb_ is advanced
			// properly
			r.patFw.clear();
			dstLen = 0;
		}
	}
	assert_eq('+', c);

	// Chew up the optional name on the '+' line
	ASSERT_ONLY(int pk =) peekToEndOfLine(fb_);
	if(charsRead == 0) {
		assert_eq('@', pk);
		fb_.get();
		fb_.resetLastN();
		cerr << "Warning: skipping empty FASTQ read with name '" << r.name << "'" << endl;
		return;
	}

	// Now read the qualities
	if (intQuals_) {
		assert(!fuzzy_);
		int qualsRead = 0;
		char buf[4096];
		if(gColor && r.primer != -1) {
			// In case the original quality string is one shorter
			mytrim5--;
		}
		qualToks_.clear();
		tokenizeQualLine(fb_, buf, 4096, qualToks_);
		for(unsigned int j = 0; j < qualToks_.size(); ++j) {
			char c = intToPhred33(atoi(qualToks_[j].c_str()), solQuals_);
			assert_geq(c, 33);
			if (qualsRead >= mytrim5) {
				//if(qualsRead - mytrim5 >= 1024) tooManyQualities(r.name);
				r.qual.append(c);
			}
			++qualsRead;
		} // done reading integer quality lines
		if(gColor && r.primer != -1) mytrim5++;
		r.qual.trimEnd(gTrim3);
		if(r.qual.length() < r.patFw.length()) {
			tooFewQualities(r.name);
		} else if(r.qual.length() > r.patFw.length() + 1) {
			tooManyQualities(r.name);
		}
		if(r.qual.length() == r.patFw.length()+1 && gColor && r.primer != -1) {
			r.qual.remove(0);
		}
		// Trim qualities on 3' end
		if(r.qual.length() > r.patFw.length()) {
			r.qual.resize(r.patFw.length());
			assert_eq((int)r.qual.length(), dstLen);
		}
		peekOverNewline(fb_);
	} else {
		// Non-integer qualities
		altBufIdx = 0;
		trim5 = mytrim5;
		int qualsRead[4] = {0, 0, 0, 0};
		int *qualsReadCur = &qualsRead[0];
		BTString *qbuf = &r.qual;
		if(gColor && r.primer != -1) {
			// In case the original quality string is one shorter
			trim5--;
		}
		while(true) {
			c = fb_.get();
			if (!fuzzy_ && c == ' ') {
				wrongQualityFormat(r.name);
			} else if(c == ' ') {
				trim5 = 0; // disable 5' trimming for now
				if((*qualsReadCur) == 0) continue;
				if(altBufIdx >= 3) {
					cerr << "At most 3 alternate quality strings permitted; offending read: " << r.name << endl;
					throw 1;
				}
				qbuf = &r.altQual[altBufIdx++];
				qualsReadCur = &qualsRead[altBufIdx];
				continue;
			}
			if(c < 0) { bail(r); return; }
			if (c != '\r' && c != '\n') {
				if (*qualsReadCur >= trim5) {
					//size_t off = (*qualsReadCur) - trim5;
					//if(off >= 1024) tooManyQualities(r.name);
					c = charToPhred33(c, solQuals_, phred64Quals_);
					assert_geq(c, 33);
					qbuf->append(c);
				}
				(*qualsReadCur)++;
			} else {
				break;
			}
		}
		qualsRead[0] -= gTrim3;
		r.qual.trimEnd(gTrim3);
		if(r.qual.length() < r.patFw.length()) {
			tooFewQualities(r.name);
		} else if(r.qual.length() > r.patFw.length()+1) {
			tooManyQualities(r.name);
		}
		if(r.qual.length() == r.patFw.length()+1 && gColor && r.primer != -1) {
			r.qual.remove(0);
		}

		if(fuzzy_) {
			// Trim from 3' end of alternate basecall and quality strings
			if(gTrim3 > 0) {
				for(int i = 0; i < 3; i++) {
					assert_eq(r.altQual[i].length(), r.altPatFw[i].length());
					if((int)r.altQual[i].length() > gTrim3) {
						r.altPatFw[i].resize(gTrim3);
						r.altQual[i].resize(gTrim3);
					} else {
						r.altPatFw[i].clear();
						r.altQual[i].clear();
					}
					qualsRead[i+1] = dstLens[i+1] =
						max<int>(0, dstLens[i+1] - gTrim3);
				}
			}
			// Shift to RHS, and install in Strings
			assert_eq(0, r.alts);
			for(int i = 1; i < 4; i++) {
				if(qualsRead[i] == 0) continue;
				if(qualsRead[i] > dstLen) {
					// Shift everybody up
					int shiftAmt = qualsRead[i] - dstLen;
					for(int j = 0; j < dstLen; j++) {
						r.altQual[i-1].set(r.altQual[i-1][j+shiftAmt], j);
						r.altPatFw[i-1].set(r.altPatFw[i-1][j+shiftAmt], j);
					}
					r.altQual[i-1].resize(dstLen);
					r.altPatFw[i-1].resize(dstLen);
				} else if (qualsRead[i] < dstLen) {
					r.altQual[i-1].resize(dstLen);
					r.altPatFw[i-1].resize(dstLen);
					// Shift everybody down
					int shiftAmt = dstLen - qualsRead[i];
					for(int j = dstLen-1; j >= shiftAmt; j--) {
						r.altQual[i-1].set(r.altQual[i-1][j-shiftAmt], j);
						r.altPatFw[i-1].set(r.altPatFw[i-1][j-shiftAmt], j);
					}
					// Fill in unset positions
					for(int j = 0; j < shiftAmt; j++) {
						// '!' - indicates no alternate basecall at
						// this position
						r.altQual[i-1].set(33, j);
					}
				}
				r.alts++;
			}
		}

		if(c == '\r' || c == '\n') {
			c = peekOverNewline(fb_);
		} else {
			c = peekToEndOfLine(fb_);
		}
	}
	r.readOrigBuf.install(fb_.lastN(), fb_.lastNLen());
	fb_.resetLastN();

	c = fb_.get();
	// Should either be at end of file or at beginning of next record
	assert(c == -1 || c == '@');

	// Set up a default name if one hasn't been set
	if(r.name.empty()) {
		char cbuf[20];
		itoa10<TReadId>(readCnt_, cbuf);
		r.name.install(cbuf);
	}
	r.trimmed3 = gTrim3;
	r.trimmed5 = mytrim5;
	readCnt_++;
	patid = readCnt_-1;
}

void wrongQualityFormat(const BTString& read_name) {
	cerr << "Encountered a space parsing the quality string for read " << read_name << endl
	     << "If this is a FASTQ file with integer (non-ASCII-encoded) qualities, please" << endl
	     << "re-run Bowtie with the --integer-quals option.  If this is a FASTQ file with" << endl
	     << "alternate basecall information, please re-run Bowtie with the --fuzzy option." << endl;
	throw 1;
}

void tooFewQualities(const BTString& read_name) {
	cerr << "Too few quality values for read: " << read_name << endl
		 << "\tare you sure this is a FASTQ-int file?" << endl;
	throw 1;
}

void tooManyQualities(const BTString& read_name) {
	cerr << "Reads file contained a pattern with more than 1024 quality values." << endl
		 << "Please truncate reads and quality values and and re-run Bowtie" << endl;
	throw 1;
}

void tooManySeqChars(const BTString& read_name) {
	cerr << "Reads file contained a pattern with more than 1024 sequence characters." << endl
		 << "Please truncate reads and quality values and and re-run Bowtie." << endl
		 << "Offending read: " << read_name << endl;
	throw 1;
}
