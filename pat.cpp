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
	const EList<string>& qs)
{
	switch(p.format) {
		case FASTA:       return new FastaPatternSource(qs, p);
		case FASTA_CONT:  return new FastaContinuousPatternSource(qs, p);
		case RAW:         return new RawPatternSource(qs, p);
		case FASTQ:       return new FastqPatternSource(qs, p);
		case TAB_MATE5:   return new TabbedPatternSource(qs, p, false);
		case TAB_MATE6:   return new TabbedPatternSource(qs, p, true);
		case CMDLINE:     return new VectorPatternSource(qs, p);
		case QSEQ:        return new QseqPatternSource(qs, p);
		default: {
			cerr << "Internal error; bad patsrc format: " << p.format << endl;
			throw 1;
		}
	}
}

/**
 * The main member function for dispensing patterns.
 *
 * Returns true iff a pair was parsed succesfully.
 */
bool PatternSource::nextReadPair(
	Read& ra,
	Read& rb,
	TReadId& rdid,
	TReadId& endid,
	bool& success,
	bool& done,
	bool& paired,
	bool fixName)
{
	// nextPatternImpl does the reading from the ultimate source;
	// it is implemented in concrete subclasses
	success = done = paired = false;
	nextReadPairImpl(ra, rb, rdid, endid, success, done, paired);
	if(success) {
		// Construct reversed versions of fw and rc seqs/quals
		ra.finalize();
		if(!rb.empty()) {
			rb.finalize();
		}
		// Fill in the random-seed field using a combination of
		// information from the user-specified seed and the read
		// sequence, qualities, and name
		ra.seed = genRandSeed(ra.patFw, ra.qual, ra.name, seed_);
		if(!rb.empty()) {
			rb.seed = genRandSeed(rb.patFw, rb.qual, rb.name, seed_);
		}
	}
	return success;
}

/**
 * The main member function for dispensing patterns.
 */
bool PatternSource::nextRead(
	Read& r,
	TReadId& rdid,
	TReadId& endid,
	bool& success,
	bool& done)
{
	// nextPatternImpl does the reading from the ultimate source;
	// it is implemented in concrete subclasses
	nextReadImpl(r, rdid, endid, success, done);
	if(success) {
		// Construct the reversed versions of the fw and rc seqs
		// and quals
		r.finalize();
		// Fill in the random-seed field using a combination of
		// information from the user-specified seed and the read
		// sequence, qualities, and name
		r.seed = genRandSeed(r.patFw, r.qual, r.name, seed_);
	}
	return success;
}

/**
 * Get the next paired or unpaired read from the wrapped
 * PairedPatternSource.
 */
bool WrappedPatternSourcePerThread::nextReadPair(
	bool& success,
	bool& done,
	bool& paired,
	bool fixName)
{
	PatternSourcePerThread::nextReadPair(success, done, paired, fixName);
	ASSERT_ONLY(TReadId lastRdId = rdid_);
	buf1_.reset();
	buf2_.reset();
	patsrc_.nextReadPair(buf1_, buf2_, rdid_, endid_, success, done, paired, fixName);
	assert(!success || rdid_ != lastRdId);
	return success;
}

/**
 * The main member function for dispensing pairs of reads or
 * singleton reads.  Returns true iff ra and rb contain a new
 * pair; returns false if ra contains a new unpaired read.
 */
bool PairedSoloPatternSource::nextReadPair(
	Read& ra,
	Read& rb,
	TReadId& rdid,
	TReadId& endid,
	bool& success,
	bool& done,
	bool& paired,
	bool fixName)
{
	uint32_t cur = cur_;
	success = false;
	while(cur < src_->size()) {
		// Patterns from srca_[cur_] are unpaired
		do {
			(*src_)[cur]->nextReadPair(
				ra, rb, rdid, endid, success, done, paired, fixName);
		} while(!success && !done);
		if(!success) {
			assert(done);
			// If patFw is empty, that's our signal that the
			// input dried up
			lock();
			if(cur + 1 > cur_) cur_++;
			cur = cur_;
			unlock();
			continue; // on to next pair of PatternSources
		}
		assert(success);
		ra.seed = genRandSeed(ra.patFw, ra.qual, ra.name, seed_);
		if(!rb.empty()) {
			rb.seed = genRandSeed(rb.patFw, rb.qual, rb.name, seed_);
			if(fixName) {
				ra.fixMateName(1);
				rb.fixMateName(2);
			}
		}
		ra.rdid = rdid;
		ra.endid = endid;
		if(!rb.empty()) {
			rb.rdid = rdid;
			rb.endid = endid+1;
		}
		ra.mate = 1;
		rb.mate = 2;
		return true; // paired
	}
	assert_leq(cur, src_->size());
	done = (cur == src_->size());
	return false;
}

/**
 * The main member function for dispensing pairs of reads or
 * singleton reads.  Returns true iff ra and rb contain a new
 * pair; returns false if ra contains a new unpaired read.
 */
bool PairedDualPatternSource::nextReadPair(
	Read& ra,
	Read& rb,
	TReadId& rdid,
	TReadId& endid,
	bool& success,
	bool& done,
	bool& paired,
	bool fixName)
{
	// 'cur' indexes the current pair of PatternSources
	uint32_t cur;
	{
		lock();
		cur = cur_;
		unlock();
	}
	success = false;
	done = true;
	while(cur < srca_->size()) {
		if((*srcb_)[cur] == NULL) {
			paired = false;
			// Patterns from srca_ are unpaired
			do {
				(*srca_)[cur]->nextRead(ra, rdid, endid, success, done);
			} while(!success && !done);
			if(!success) {
				assert(done);
				lock();
				if(cur + 1 > cur_) cur_++;
				cur = cur_; // Move on to next PatternSource
				unlock();
				continue; // on to next pair of PatternSources
			}
			ra.rdid = rdid;
			ra.endid = endid;
			ra.mate  = 0;
			return success;
		} else {
			paired = true;
			// Patterns from srca_[cur_] and srcb_[cur_] are paired
			TReadId rdid_a = 0, endid_a = 0;
			TReadId rdid_b = 0, endid_b = 0;
			bool success_a = false, done_a = false;
			bool success_b = false, done_b = false;
			// Lock to ensure that this thread gets parallel reads
			// in the two mate files
			lock();
			do {
				(*srca_)[cur]->nextRead(ra, rdid_a, endid_a, success_a, done_a);
			} while(!success_a && !done_a);
			do {
				(*srcb_)[cur]->nextRead(rb, rdid_b, endid_b, success_b, done_b);
			} while(!success_b && !done_b);
			if(!success_a && success_b) {
				cerr << "Error, fewer reads in file specified with -1 than in file specified with -2" << endl;
				throw 1;
			} else if(!success_a) {
				assert(done_a && done_b);
				if(cur + 1 > cur_) cur_++;
				cur = cur_; // Move on to next PatternSource
				unlock();
				continue; // on to next pair of PatternSources
			} else if(!success_b) {
				cerr << "Error, fewer reads in file specified with -2 than in file specified with -1" << endl;
				throw 1;
			}
			assert_eq(rdid_a, rdid_b);
			//assert_eq(endid_a+1, endid_b);
			assert_eq(success_a, success_b);
			unlock();
			if(fixName) {
				ra.fixMateName(1);
				rb.fixMateName(2);
			}
			rdid = rdid_a;
			endid = endid_a;
			success = success_a;
			done = done_a;
			ra.rdid = rdid;
			ra.endid = endid;
			if(!rb.empty()) {
				rb.rdid = rdid;
				rb.endid = endid+1;
			}
			ra.mate = 1;
			rb.mate = 2;
			return success;
		}
	}
	return success;
}

/**
 * Return the number of reads attempted.
 */
pair<TReadId, TReadId> PairedDualPatternSource::readCnt() const {
	uint64_t rets = 0llu, retp = 0llu;
	for(size_t i = 0; i < srca_->size(); i++) {
		if((*srcb_)[i] == NULL) {
			rets += (*srca_)[i]->readCnt();
		} else {
			assert_eq((*srca_)[i]->readCnt(), (*srcb_)[i]->readCnt());
			retp += (*srca_)[i]->readCnt();
		}
	}
	return make_pair(rets, retp);
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
		ab->push_back(PatternSource::patsrcFromStrings(p, *qs));
		if(!p.fileParallel) {
			break;
		}
	}

	// Create list of pattern sources for paired reads
	for(size_t i = 0; i < m1.size(); i++) {
		const EList<string>* qs = &m1;
		EList<string> tmpSeq;
		EList<string> tmpQual;
		if(p.fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmpSeq;
			tmpSeq.push_back(m1[i]);
			assert_eq(1, tmpSeq.size());
		}
		a->push_back(PatternSource::patsrcFromStrings(p, *qs));
		if(!p.fileParallel) {
			break;
		}
	}

	// Create list of pattern sources for paired reads
	for(size_t i = 0; i < m2.size(); i++) {
		const EList<string>* qs = &m2;
		EList<string> tmpSeq;
		EList<string> tmpQual;
		if(p.fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmpSeq;
			tmpSeq.push_back(m2[i]);
			assert_eq(1, tmpSeq.size());
		}
		b->push_back(PatternSource::patsrcFromStrings(p, *qs));
		if(!p.fileParallel) {
			break;
		}
	}
	// All mates/mate files must be paired
	assert_eq(a->size(), b->size());

	// Create list of pattern sources for the unpaired reads
	for(size_t i = 0; i < si.size(); i++) {
		const EList<string>* qs = &si;
		PatternSource* patsrc = NULL;
		EList<string> tmpSeq;
		EList<string> tmpQual;
		if(p.fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmpSeq;
			tmpSeq.push_back(si[i]);
			assert_eq(1, tmpSeq.size());
		}
		patsrc = PatternSource::patsrcFromStrings(p, *qs);
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

void PairedPatternSource::free_EList_pmembers( const EList<PatternSource*> &elist) {
    for (size_t i = 0; i < elist.size(); i++)
        if (elist[i] != NULL)
            delete elist[i];
}

VectorPatternSource::VectorPatternSource(
	const EList<string>& v,
	const PatternParams& p) :
	PatternSource(p),
	cur_(p.skip),
	skip_(p.skip),
	paired_(false),
	v_(),
	quals_()
{
	for(size_t i = 0; i < v.size(); i++) {
		EList<string> ss;
		tokenize(v[i], ":", ss, 2);
		assert_gt(ss.size(), 0);
		assert_leq(ss.size(), 2);
		// Initialize s
		string s = ss[0];
		int mytrim5 = gTrim5;
		if(gColor && s.length() > 1) {
			// This may be a primer character.  If so, keep it in the
			// 'primer' field of the read buf and parse the rest of the
			// read without it.
			int c = toupper(s[0]);
			if(asc2dnacat[c] > 0) {
				// First char is a DNA char
				int c2 = toupper(s[1]);
				// Second char is a color char
				if(asc2colcat[c2] > 0) {
					mytrim5 += 2; // trim primer and first color
				}
			}
		}
		if(gColor) {
			// Convert '0'-'3' to 'A'-'T'
			for(size_t i = 0; i < s.length(); i++) {
				if(s[i] >= '0' && s[i] <= '4') {
					s[i] = "ACGTN"[(int)s[i] - '0'];
				}
				if(s[i] == '.') s[i] = 'N';
			}
		}
		if(s.length() <= (size_t)(gTrim3 + mytrim5)) {
			// Entire read is trimmed away
			s.clear();
		} else {
			// Trim on 5' (high-quality) end
			if(mytrim5 > 0) {
				s.erase(0, mytrim5);
			}
			// Trim on 3' (low-quality) end
			if(gTrim3 > 0) {
				s.erase(s.length()-gTrim3);
			}
		}
		//  Initialize vq
		string vq;
		if(ss.size() == 2) {
			vq = ss[1];
		}
		// Trim qualities
		if(vq.length() > (size_t)(gTrim3 + mytrim5)) {
			// Trim on 5' (high-quality) end
			if(mytrim5 > 0) {
				vq.erase(0, mytrim5);
			}
			// Trim on 3' (low-quality) end
			if(gTrim3 > 0) {
				vq.erase(vq.length()-gTrim3);
			}
		}
		// Pad quals with Is if necessary; this shouldn't happen
		while(vq.length() < s.length()) {
			vq.push_back('I');
		}
		// Truncate quals to match length of read if necessary;
		// this shouldn't happen
		if(vq.length() > s.length()) {
			vq.erase(s.length());
		}
		assert_eq(vq.length(), s.length());
		v_.expand();
		v_.back().installChars(s);
		quals_.push_back(BTString(vq));
		trimmed3_.push_back(gTrim3);
		trimmed5_.push_back(mytrim5);
		ostringstream os;
		os << (names_.size());
		names_.push_back(BTString(os.str()));
	}
	assert_eq(v_.size(), quals_.size());
}
	
bool VectorPatternSource::nextReadImpl(
	Read& r,
	TReadId& rdid,
	TReadId& endid,
	bool& success,
	bool& done)
{
	// Let Strings begin at the beginning of the respective bufs
	r.reset();
	lock();
	if(cur_ >= v_.size()) {
		unlock();
		// Clear all the Strings, as a signal to the caller that
		// we're out of reads
		r.reset();
		success = false;
		done = true;
		assert(r.empty());
		return false;
	}
	// Copy v_*, quals_* strings into the respective Strings
	r.color = gColor;
	r.patFw  = v_[cur_];
	r.qual = quals_[cur_];
	r.trimmed3 = trimmed3_[cur_];
	r.trimmed5 = trimmed5_[cur_];
	ostringstream os;
	os << cur_;
	r.name = os.str();
	cur_++;
	done = cur_ == v_.size();
	rdid = endid = readCnt_;
	readCnt_++;
	unlock();
	success = true;
	return true;
}
	
/**
 * This is unused, but implementation is given for completeness.
 */
bool VectorPatternSource::nextReadPairImpl(
	Read& ra,
	Read& rb,
	TReadId& rdid,
	TReadId& endid,
	bool& success,
	bool& done,
	bool& paired)
{
	// Let Strings begin at the beginning of the respective bufs
	ra.reset();
	rb.reset();
	paired = true;
	if(!paired_) {
		paired_ = true;
		cur_ <<= 1;
	}
	lock();
	if(cur_ >= v_.size()-1) {
		unlock();
		// Clear all the Strings, as a signal to the caller that
		// we're out of reads
		ra.reset();
		rb.reset();
		assert(ra.empty());
		assert(rb.empty());
		success = false;
		done = true;
		return false;
	}
	// Copy v_*, quals_* strings into the respective Strings
	ra.patFw  = v_[cur_];
	ra.qual = quals_[cur_];
	ra.trimmed3 = trimmed3_[cur_];
	ra.trimmed5 = trimmed5_[cur_];
	cur_++;
	rb.patFw  = v_[cur_];
	rb.qual = quals_[cur_];
	rb.trimmed3 = trimmed3_[cur_];
	rb.trimmed5 = trimmed5_[cur_];
	ostringstream os;
	os << readCnt_;
	ra.name = os.str();
	rb.name = os.str();
	ra.color = rb.color = gColor;
	cur_++;
	done = cur_ >= v_.size()-1;
	rdid = endid = readCnt_;
	readCnt_++;
	unlock();
	success = true;
	return true;
}

/**
 * Parse a single quality string from fb and store qualities in r.
 * Assume the next character obtained via fb.get() is the first
 * character of the quality string.  When returning, the next
 * character returned by fb.peek() or fb.get() should be the first
 * character of the following line.
 */
int parseQuals(
	Read& r,
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

/// Read another pattern from a FASTA input file
bool FastaPatternSource::read(
	Read& r,
	TReadId& rdid,
	TReadId& endid,
	bool& success,
	bool& done)
{
	int c, qc = 0;
	success = true;
	done = false;
	assert(fb_.isOpen());
	r.reset();
	r.color = gColor;
	// Pick off the first carat
	c = fb_.get();
	if(c < 0) {
		bail(r); success = false; done = true; return success;
	}
	while(c == '#' || c == ';' || c == '\r' || c == '\n') {
		c = fb_.peekUptoNewline();
		fb_.resetLastN();
		c = fb_.get();
	}
	assert_eq(1, fb_.lastNLen());

	// Pick off the first carat
	if(first_) {
		if(c != '>') {
			cerr << "Error: reads file does not look like a FASTA file" << endl;
			throw 1;
		}
		first_ = false;
	}
	assert_eq('>', c);
	c = fb_.get(); // get next char after '>'

	// Read to the end of the id line, sticking everything after the '>'
	// into *name
	//bool warning = false;
	while(true) {
		if(c < 0 || qc < 0) {
			bail(r); success = false; done = true; return success;
		}
		if(c == '\n' || c == '\r') {
			// Break at end of line, after consuming all \r's, \n's
			while(c == '\n' || c == '\r') {
				if(fb_.peek() == '>') {
					// Empty sequence
					break;
				}
				c = fb_.get();
				if(c < 0 || qc < 0) {
					bail(r); success = false; done = true; return success;
				}
			}
			break;
		}
		r.name.append(c);
		if(fb_.peek() == '>') {
			// Empty sequence
			break;
		}
		c = fb_.get();
	}
	if(c == '>') {
		// Empty sequences!
		cerr << "Warning: skipping empty FASTA read with name '" << r.name << "'" << endl;
		fb_.resetLastN();
		rdid = endid = readCnt_;
		readCnt_++;
		success = true; done = false; return success;
	}
	assert_neq('>', c);

	// _in now points just past the first character of a sequence
	// line, and c holds the first character
	int begin = 0;
	int mytrim5 = gTrim5;
	if(gColor) {
		// This is the primer character, keep it in the
		// 'primer' field of the read buf and keep parsing
		c = toupper(c);
		if(asc2dnacat[c] > 0) {
			// First char is a DNA char
			int c2 = toupper(fb_.peek());
			if(asc2colcat[c2] > 0) {
				// Second char is a color char
				r.primer = c;
				r.trimc = c2;
				mytrim5 += 2;
			}
		}
		if(c < 0) {
			bail(r); success = false; done = true; return success;
		}
	}
	while(c != '>' && c >= 0) {
		if(gColor) {
			if(c >= '0' && c <= '4') c = "ACGTN"[(int)c - '0'];
			if(c == '.') c = 'N';
		}
		if(asc2dnacat[c] > 0 && begin++ >= mytrim5) {
			r.patFw.append(asc2dna[c]);
			r.qual.append('I');
		}
		if(fb_.peek() == '>') break;
		c = fb_.get();
	}
	r.patFw.trimEnd(gTrim3);
	r.qual.trimEnd(gTrim3);
	r.trimmed3 = gTrim3;
	r.trimmed5 = mytrim5;
	// Set up a default name if one hasn't been set
	if(r.name.empty()) {
		char cbuf[20];
		itoa10<TReadId>(readCnt_, cbuf);
		r.name.install(cbuf);
	}
	assert_gt(r.name.length(), 0);
	r.readOrigBuf.install(fb_.lastN(), fb_.lastNLen());
	fb_.resetLastN();
	rdid = endid = readCnt_;
	readCnt_++;
	return success;
}

/// Read another pattern from a FASTQ input file
bool FastqPatternSource::read(
	Read& r,
	TReadId& rdid,
	TReadId& endid,
	bool& success,
	bool& done)
{
	int c;
	int dstLen = 0;
	success = true;
	done = false;
	r.reset();
	r.color = gColor;
	r.fuzzy = fuzzy_;
	// Pick off the first at
	if(first_) {
		c = fb_.get();
		if(c != '@') {
			c = getOverNewline(fb_);
			if(c < 0) {
				bail(r); success = false; done = true; return success;
			}
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
		if(c < 0) {
			bail(r); success = false; done = true; return success;
		}
		if(c == '\n' || c == '\r') {
			// Break at end of line, after consuming all \r's, \n's
			while(c == '\n' || c == '\r') {
				c = fb_.get();
				if(c < 0) {
					bail(r); success = false; done = true; return success;
				}
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
		if(c < 0) {
			bail(r); success = false; done = true; return success;
		}
	}
	int trim5 = 0;
	if(c != '+') {
		trim5 = mytrim5;
		while(c != '+') {
			// Convert color numbers to letters if necessary
			if(c == '.') c = 'N';
			if(gColor) {
				if(c >= '0' && c <= '4') c = "ACGTN"[(int)c - '0'];
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
			if(c < 0) {
				bail(r); success = false; done = true; return success;
			}
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
		rdid = endid = readCnt_;
		readCnt_++;
		return success;
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
			if(c < 0) {
				break; // let the file end just at the end of a quality line
				//bail(r); success = false; done = true; return success;
			}
			if (c != '\r' && c != '\n') {
				if (*qualsReadCur >= trim5) {
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
	rdid = endid = readCnt_;
	readCnt_++;
	return success;
}

/// Read another pattern from a FASTA input file
bool TabbedPatternSource::read(
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
	// fb_ is about to dish out the first character of the
	// name field
	if(parseName(r, NULL, '\t') == -1) {
		peekOverNewline(fb_); // skip rest of line
		r.reset();
		success = false;
		done = true;
		return false;
	}
	assert_neq('\t', fb_.peek());

	// fb_ is about to dish out the first character of the
	// sequence field
	int charsRead = 0;
	int mytrim5 = gTrim5;
	int dstLen = parseSeq(r, charsRead, mytrim5, '\t');
	assert_neq('\t', fb_.peek());
	if(dstLen < 0) {
		peekOverNewline(fb_); // skip rest of line
		r.reset();
		success = false;
		done = true;
		return false;
	}

	// fb_ is about to dish out the first character of the
	// quality-string field
	char ct = 0;
	if(parseQuals(r, charsRead, dstLen, mytrim5, ct, '\n') < 0) {
		peekOverNewline(fb_); // skip rest of line
		r.reset();
		success = false;
		done = true;
		return false;
	}
	r.trimmed3 = gTrim3;
	r.trimmed5 = mytrim5;
	assert_eq(ct, '\n');
	assert_neq('\n', fb_.peek());
	r.readOrigBuf.install(fb_.lastN(), fb_.lastNLen());
	fb_.resetLastN();
	rdid = endid = readCnt_;
	readCnt_++;
	return true;
}

/// Read another pair of patterns from a FASTA input file
bool TabbedPatternSource::readPair(
	Read& ra,
	Read& rb,
	TReadId& rdid,
	TReadId& endid,
	bool& success,
	bool& done,
	bool& paired)
{
	success = true;
	done = false;
	
	// Skip over initial vertical whitespace
	if(fb_.peek() == '\r' || fb_.peek() == '\n') {
		fb_.peekUptoNewline();
		fb_.resetLastN();
	}
	
	// fb_ is about to dish out the first character of the
	// name field
	int mytrim5_1 = gTrim5;
	if(parseName(ra, &rb, '\t') == -1) {
		peekOverNewline(fb_); // skip rest of line
		ra.reset();
		rb.reset();
		fb_.resetLastN();
		success = false;
		done = true;
		return false;
	}
	assert_neq('\t', fb_.peek());

	// fb_ is about to dish out the first character of the
	// sequence field for the first mate
	int charsRead1 = 0;
	int dstLen1 = parseSeq(ra, charsRead1, mytrim5_1, '\t');
	if(dstLen1 < 0) {
		peekOverNewline(fb_); // skip rest of line
		ra.reset();
		rb.reset();
		fb_.resetLastN();
		success = false;
		done = true;
		return false;
	}
	assert_neq('\t', fb_.peek());

	// fb_ is about to dish out the first character of the
	// quality-string field
	char ct = 0;
	if(parseQuals(ra, charsRead1, dstLen1, mytrim5_1, ct, '\t', '\n') < 0) {
		peekOverNewline(fb_); // skip rest of line
		ra.reset();
		rb.reset();
		fb_.resetLastN();
		success = false;
		done = true;
		return false;
	}
	ra.trimmed3 = gTrim3;
	ra.trimmed5 = mytrim5_1;
	assert(ct == '\t' || ct == '\n' || ct == '\r' || ct == -1);
	if(ct == '\r' || ct == '\n' || ct == -1) {
		// Only had 3 fields prior to newline, so this must be an unpaired read
		rb.reset();
		ra.readOrigBuf.install(fb_.lastN(), fb_.lastNLen());
		fb_.resetLastN();
		success = true;
		done = false;
		paired = false;
		rdid = endid = readCnt_;
		readCnt_++;
		return success;
	}
	paired = true;
	assert_neq('\t', fb_.peek());
	
	// Saw another tab after the third field, so this must be a pair
	if(secondName_) {
		// The second mate has its own name
		if(parseName(rb, NULL, '\t') == -1) {
			peekOverNewline(fb_); // skip rest of line
			ra.reset();
			rb.reset();
			fb_.resetLastN();
			success = false;
			done = true;
			return false;
		}
		assert_neq('\t', fb_.peek());
	}

	// fb_ about to give the first character of the second mate's sequence
	int charsRead2 = 0;
	int mytrim5_2 = gTrim5;
	int dstLen2 = parseSeq(rb, charsRead2, mytrim5_2, '\t');
	if(dstLen2 < 0) {
		peekOverNewline(fb_); // skip rest of line
		ra.reset();
		rb.reset();
		fb_.resetLastN();
		success = false;
		done = true;
		return false;
	}
	assert_neq('\t', fb_.peek());

	// fb_ is about to dish out the first character of the
	// quality-string field
	if(parseQuals(rb, charsRead2, dstLen2, mytrim5_2, ct, '\n') < 0) {
		peekOverNewline(fb_); // skip rest of line
		ra.reset();
		rb.reset();
		fb_.resetLastN();
		success = false;
		done = true;
		return false;
	}
	ra.readOrigBuf.install(fb_.lastN(), fb_.lastNLen());
	fb_.resetLastN();
	rb.trimmed3 = gTrim3;
	rb.trimmed5 = mytrim5_2;
	rdid = endid = readCnt_;
	readCnt_++;
	return true;
}

/**
 * Parse a name from fb_ and store in r.  Assume that the next
 * character obtained via fb_.get() is the first character of
 * the sequence and the string stops at the next char upto (could
 * be tab, newline, etc.).
 */
int TabbedPatternSource::parseName(
	Read& r,
	Read* r2,
	char upto /* = '\t' */)
{
	// Read the name out of the first field
	int c = 0;
	if(r2 != NULL) r2->name.clear();
	r.name.clear();
	while(true) {
		if((c = fb_.get()) < 0) {
			return -1;
		}
		if(c == upto) {
			// Finished with first field
			break;
		}
		if(c == '\n' || c == '\r') {
			return -1;
		}
		if(r2 != NULL) r2->name.append(c);
		r.name.append(c);
	}
	// Set up a default name if one hasn't been set
	if(r.name.empty()) {
		char cbuf[20];
		itoa10<TReadId>(readCnt_, cbuf);
		r.name.install(cbuf);
		if(r2 != NULL) r2->name.install(cbuf);
	}
	return (int)r.name.length();
}

/**
 * Parse a single sequence from fb_ and store in r.  Assume
 * that the next character obtained via fb_.get() is the first
 * character of the sequence and the sequence stops at the next
 * char upto (could be tab, newline, etc.).
 */
int TabbedPatternSource::parseSeq(
	Read& r,
	int& charsRead,
	int& trim5,
	char upto /*= '\t'*/)
{
	int begin = 0;
	int c = fb_.get();
	assert(c != upto);
	r.patFw.clear();
	r.color = gColor;
	if(gColor) {
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
		if(gColor) {
			if(c >= '0' && c <= '4') c = "ACGTN"[(int)c - '0'];
			if(c == '.') c = 'N';
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
int TabbedPatternSource::parseQuals(
	Read& r,
	int charsRead,
	int dstLen,
	int trim5,
	char& c2,
	char upto /*= '\t'*/,
	char upto2 /*= -1*/)
{
	int qualsRead = 0;
	int c = 0;
	if (intQuals_) {
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
		if(qualsRead < dstLen + trim5) {
			tooFewQualities(r.name);
		} else if(qualsRead > dstLen + trim5) {
			tooManyQualities(r.name);
		}
	}
	r.qual.resize(dstLen);
	while(c != upto && (upto2 == -1 || c != upto2) && c != -1) {
		c = fb_.get();
		c2 = c;
	}
	return qualsRead;
}

void wrongQualityFormat(const BTString& read_name) {
	cerr << "Error: Encountered one or more spaces while parsing the quality "
	     << "string for read " << read_name << ".  If this is a FASTQ file "
		 << "with integer (non-ASCII-encoded) qualities, try re-running with "
		 << "the --integer-quals option." << endl;
	throw 1;
}

void tooFewQualities(const BTString& read_name) {
	cerr << "Error: Read " << read_name << " has more read characters than "
		 << "quality values." << endl;
	throw 1;
}

void tooManyQualities(const BTString& read_name) {
	cerr << "Error: Read " << read_name << " has more quality values than read "
		 << "characters." << endl;
	throw 1;
}
