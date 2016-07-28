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
 * Calculate a per-read random seed based on a combination of
 * the read data (incl. sequence, name, quals) and the global
 * seed in '_randSeed'.
 */
static uint32_t genRandSeed(
	const BTDnaString& qry,
	const BTString& qual,
	const BTString& name,
	uint32_t seed)
{
	// Calculate a per-read random seed based on a combination of
	// the read data (incl. sequence, name, quals) and the global
	// seed
	uint32_t rseed = (seed + 101) * 59 * 61 * 67 * 71 * 73 * 79 * 83;
	size_t qlen = qry.length();
	// Throw all the characters of the read into the random seed
	for(size_t i = 0; i < qlen; i++) {
		int p = (int)qry[i];
		assert_leq(p, 4);
		size_t off = ((i & 15) << 1);
		rseed ^= (p << off);
	}
	// Throw all the quality values for the read into the random
	// seed
	for(size_t i = 0; i < qlen; i++) {
		int p = (int)qual[i];
		assert_leq(p, 255);
		size_t off = ((i & 3) << 3);
		rseed ^= (p << off);
	}
	// Throw all the characters in the read name into the random
	// seed
	size_t namelen = name.length();
	for(size_t i = 0; i < namelen; i++) {
		int p = (int)name[i];
		if(p == '/') break;
		assert_leq(p, 255);
		size_t off = ((i & 3) << 3);
		rseed ^= (p << off);
	}
	return rseed;
}

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
 * Once name/sequence/qualities have been parsed for an
 * unpaired read, set all the other key fields of the Read
 * struct.
 */
void PatternSourcePerThread::finalize(Read& ra) {
	ra.mate = 1;
	ra.rdid = buf_.rdid();
	ra.seed = genRandSeed(ra.patFw, ra.qual, ra.name, pp_.seed);
	ra.finalize();
	if(pp_.fixName) {
		ra.fixMateName(1);
	}
}

/**
 * Once name/sequence/qualities have been parsed for a
 * paired-end read, set all the other key fields of the Read
 * structs.
 */
void PatternSourcePerThread::finalizePair(Read& ra, Read& rb) {
	ra.mate = 1;
	rb.mate = 2;
	ra.rdid = rb.rdid = buf_.rdid();
	ra.seed = genRandSeed(ra.patFw, ra.qual, ra.name, pp_.seed);
	rb.seed = genRandSeed(rb.patFw, rb.qual, rb.name, pp_.seed);
	ra.finalize();
	rb.finalize();
	if(pp_.fixName) {
		ra.fixMateName(1);
		rb.fixMateName(2);
	}
}

/**
 * Get the next paired or unpaired read from the wrapped
 * PatternComposer.  Returns a pair of bools; first indicates
 * whether we were successful, second indicates whether we're
 * done.
 */
pair<bool, bool> PatternSourcePerThread::nextReadPair() {
	// Prepare batch
	if(buf_.exhausted()) {
		pair<bool, int> res = nextBatch();
		if(res.first && res.second == 0) {
			return make_pair(false, true);
		}
		last_batch_ = res.first;
		last_batch_size_ = res.second;
		assert_eq(0, buf_.cur_buf_);
	} else {
		buf_.next(); // advance cursor
		assert_gt(buf_.cur_buf_, 0);
	}
	// Parse read/pair
	assert(!buf_.read_a().readOrigBuf.empty());
	assert(buf_.read_a().empty());
	if(!parse(buf_.read_a(), buf_.read_b())) {
		return make_pair(false, false);
	}
	// Finalize read/pair
	if(!buf_.read_b().readOrigBuf.empty()) {
		finalizePair(buf_.read_a(), buf_.read_b());
	} else {
		finalize(buf_.read_a());
	}
	bool this_is_last = buf_.cur_buf_ == last_batch_size_-1;
	return make_pair(true, this_is_last ? last_batch_ : false);
}

/**
 * The main member function for dispensing pairs of reads or
 * singleton reads.  Returns true iff ra and rb contain a new
 * pair; returns false if ra contains a new unpaired read.
 */
pair<bool, int> SoloPatternComposer::nextBatch(PerThreadReadBuf& pt) {
	size_t cur = cur_;
	while(cur < src_->size()) {
		// Patterns from srca_[cur_] are unpaired
		pair<bool, int> res;
		do {
			res = (*src_)[cur]->nextBatch(
				pt,
				true,  // batch A (or pairs)
				true); // grab lock below
		} while(!res.first && res.second == 0);
		if(res.second == 0) {
			ThreadSafe ts(&mutex_m);
			if(cur + 1 > cur_) {
				cur_++;
			}
			cur = cur_;
			continue; // on to next pair of PatternSources
		}
		return res;
	}
	assert_leq(cur, src_->size());
	return make_pair(true, 0);
}

/**
 * The main member function for dispensing pairs of reads or
 * singleton reads.  Returns true iff ra and rb contain a new
 * pair; returns false if ra contains a new unpaired read.
 */
pair<bool, int> DualPatternComposer::nextBatch(PerThreadReadBuf& pt) {
	// 'cur' indexes the current pair of PatternSources
	size_t cur = cur_;
	while(cur < srca_->size()) {
		if((*srcb_)[cur] == NULL) {
			// Patterns from srca_ are unpaired
			pair<bool, int> res = (*srca_)[cur]->nextBatch(
				pt,
				true,  // batch A (or pairs)
				true); // grab lock below
			bool done = res.first;
			if(!done && res.second == 0) {
				ThreadSafe ts(&mutex_m);
				if(cur + 1 > cur_) cur_++;
				cur = cur_; // Move on to next PatternSource
				continue; // on to next pair of PatternSources
			}
			return make_pair(done, res.second);
		} else {
			pair<bool, int> resa, resb;
			// Lock to ensure that this thread gets parallel reads
			// in the two mate files
			{
				ThreadSafe ts(&mutex_m);
				resa = (*srca_)[cur]->nextBatch(
					pt,
					true,   // batch A
					false); // don't grab lock below
				resb = (*srcb_)[cur]->nextBatch(
					pt,
					false,  // batch B
					false); // don't grab lock below
				assert_eq((*srca_)[cur]->readCount(),
				          (*srcb_)[cur]->readCount());
			}
			if(resa.second < resb.second) {
				cerr << "Error, fewer reads in file specified with -1 "
				     << "than in file specified with -2" << endl;
				throw 1;
			} else if(resa.second == 0 && resb.second == 0) {
				ThreadSafe ts(&mutex_m);
				if(cur + 1 > cur_) {
					cur_++;
				}
				cur = cur_; // Move on to next PatternSource
				continue; // on to next pair of PatternSources
			} else if(resb.second < resa.second) {
				cerr << "Error, fewer reads in file specified with -2 "
				     << "than in file specified with -1" << endl;
				throw 1;
			}
			assert_eq(resa.first, resb.first);
			assert_eq(resa.second, resb.second);
			return make_pair(resa.first, resa.second);
		}
	}
	assert_leq(cur, srca_->size());
	return make_pair(true, 0);
}

/**
 * Given the values for all of the various arguments used to specify
 * the read and quality input, create a list of pattern sources to
 * dispense them.
 */
PatternComposer* PatternComposer::setupPatternComposer(
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

	PatternComposer *patsrc = NULL;
	if(m12.size() > 0) {
		patsrc = new SoloPatternComposer(ab, p);
		for(size_t i = 0; i < a->size(); i++) delete (*a)[i];
		for(size_t i = 0; i < b->size(); i++) delete (*b)[i];
		delete a; delete b;
	} else {
		patsrc = new DualPatternComposer(a, b, p);
		for(size_t i = 0; i < ab->size(); i++) delete (*ab)[i];
		delete ab;
	}
	return patsrc;
}

void PatternComposer::free_EList_pmembers( const EList<PatternSource*> &elist) {
    for (size_t i = 0; i < elist.size(); i++)
        if (elist[i] != NULL)
            delete elist[i];
}

/**
 * Fill Read with the sequence, quality and name for the next
 * read in the list of read files.  This function gets called by
 * all the search threads, so we must handle synchronization.
 *
 * What does the return value signal?
 * In the event that we read more data, it should at least signal how many
 * reads were read, and whether we're totally done.  It's debatable whether
 * it should convey anything about the individual read files, like whether
 * we finished one of them.
 */
pair<bool, int> BufferedFilePatternSource::nextBatch(
	PerThreadReadBuf& pt,
	bool batch_a,
	bool lock)
{
	bool done = false;
	int nread = 0;
	
	// synchronization at this level because both reading and manipulation of
	// current file pointer have to be protected
	ThreadSafe ts(&mutex, lock);
	pt.setReadId(readCnt_);
	while(true) { // loop that moves on to next file when needed
		do {
			pair<bool, int> ret = nextBatchFromFile(pt, batch_a);
			done = ret.first;
			nread = ret.second;
		} while(!done && nread == 0); // not sure why this would happen
		if(done && filecur_ < infiles_.size()) { // finished with this file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			filecur_++;
			if(nread == 0) {
				continue;
			}
		}
		break;
	}
	assert_geq(nread, 0);
	readCnt_ += nread;
	return make_pair(done, nread);
}


/**
 * Fill Read with the sequence, quality and name for the next
 * read in the list of read files.  This function gets called by
 * all the search threads, so we must handle synchronization.
 *
 * Returns pair<bool, int> where bool indicates whether we're
 * completely done, and int indicates how many reads were read.
 */
pair<bool, int> CFilePatternSource::nextBatch(
	PerThreadReadBuf& pt,
	bool batch_a,
	bool lock)
{
	bool done = false;
	int nread = 0;
	
	// synchronization at this level because both reading and manipulation of
	// current file pointer have to be protected
	ThreadSafe ts(&mutex, lock);
	pt.setReadId(readCnt_);
	while(true) { // loop that moves on to next file when needed
		do {
			pair<bool, int> ret = nextBatchFromFile(pt, batch_a);
			done = ret.first;
			nread = ret.second;
		} while(!done && nread == 0); // not sure why this would happen
		if(done && filecur_ < infiles_.size()) { // finished with this file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			filecur_++;
			if(nread == 0) {
				continue;
			}
		}
		break;
	}
	assert_geq(nread, 0);
	readCnt_ += nread;
	return make_pair(done, nread);
}

/**
 * Open the next file in the list of input files.
 */
void BufferedFilePatternSource::open() {
	if(fb_.isOpen()) {
		fb_.close();
	}
	while(filecur_ < infiles_.size()) {
		// Open read
		FILE *in;
		if(infiles_[filecur_] == "-") {
			in = stdin;
		} else if((in = fopen(infiles_[filecur_].c_str(), "rb")) == NULL) {
			if(!errs_[filecur_]) {
				cerr << "Warning: Could not open read file \""
					 << infiles_[filecur_].c_str()
					 << "\" for reading; skipping..." << endl;
				errs_[filecur_] = true;
			}
			filecur_++;
			continue;
		}
		fb_.newFile(in);
		return;
	}
	cerr << "Error: No input read files were valid" << endl;
	exit(1);
	return;
}

/**
 * Open the next file in the list of input files.
 */
void CFilePatternSource::open() {
	if(is_open_) {
		is_open_ = false;
		fclose(fp_);
		fp_ = NULL;
	}
	while(filecur_ < infiles_.size()) {
		if(infiles_[filecur_] == "-") {
			fp_ = stdin;
		} else if((fp_ = fopen(infiles_[filecur_].c_str(), "rb")) == NULL) {
			if(!errs_[filecur_]) {
				cerr << "Warning: Could not open read file \""
				<< infiles_[filecur_].c_str()
				<< "\" for reading; skipping..." << endl;
				errs_[filecur_] = true;
			}
			filecur_++;
			continue;
		}
		is_open_ = true;
		setvbuf(fp_, buf_, _IOFBF, 64*1024);
		return;
	}
	cerr << "Error: No input read files were valid" << endl;
	exit(1);
	return;
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
		// TODO: more work-intensive than it should be
		ostringstream os;
		os << (names_.size());
		names_.push_back(BTString(os.str()));
	}
	assert_eq(v_.size(), quals_.size());
}

/**
 * Read next batch.  However, batch concept is not very applicable for this
 * PatternSource where all the info has already been parsed into the fields
 * in the contsructor.  This essentially modifies the pt as though we read
 * in some number of patterns.
 */
pair<bool, int> VectorPatternSource::nextBatch(
	PerThreadReadBuf& pt,
	bool batch_a,
	bool lock)
{
	bool success = true;
	int nread = 0;
	pt.reset();
	ThreadSafe ts(&mutex, lock);
	pt.setReadId(readCnt_);
#if 0
	// TODO: set nread to min of pt.size() and total - cur_
	// TODO: implement something like following function
	pt.install_dummies(nread);
#endif
	readCnt_ += nread;
	return make_pair(success, nread);
}

#if 0
/**
 * This is unused, but implementation is given for completeness.
 */
pair<bool, int> VectorPatternSource::nextBatchPair(PerThreadReadBuf& pt)
{
	bool success = true;
	int nread = 0;
	// Let Strings begin at the beginning of the respective bufs
	ra.reset();
	rb.reset();
	paired = true;
	if(!paired_) {
		paired_ = true;
		cur_ <<= 1;
	}
	ThreadSafe ts(&mutex);
	pt.setReadId(readCnt_);
	if(cur_ >= v_.size()-1) {
		ts.~ThreadSafe();
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
	cur_++;
	done = cur_ >= v_.size()-1;
	rdid = readCnt_;
	readCnt_++;
	success = true;
	return make_pair(success, nread);
}
#endif

/**
 * Parse a single quality string from fb and store qualities in r.
 * Assume the next character obtained via fb.get() is the first
 * character of the quality string.  When returning, the next
 * character returned by fb.peek() or fb.get() should be the first
 * character of the following line.
 */
static int parseQuals(
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
	if ((int)r.qual.length() < readLen) {
		tooFewQualities(r.name);
	}
	r.qual.trimEnd(trim3);
	r.qual.trimBegin(trim5);
	if(r.qual.length() <= 0) return 0;
	assert_eq(r.qual.length(), r.patFw.length());
	while(fb.peek() == '\n' || fb.peek() == '\r') fb.get();
	return (int)r.qual.length();
}

/// Read another pattern from a FASTA input file
pair<bool, int> FastaPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a)
{
	int c;
	EList<Read>& readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	if(first_) {
		c = fb_.get();
		while(c == '\r' || c == '\n') {
			c = fb_.get();
		}
		if(c != '>') {
			cerr << "Error: reads file does not look like a FASTA file" << endl;
			throw 1;
		}
		first_ = false;
	}
	bool done = false;
	size_t readi = 0;
	// Read until we run out of input or until we've filled the buffer
	for(; readi < pt.max_buf_ && !done; readi++) {
		SStringExpandable<char, 1024, 2, 1024>& buf = readbuf[readi].readOrigBuf;
		buf.clear();
		buf.append('>');
		while(true) {
			c = fb_.get();
			done = c < 0;
			if(c < 0 || c == '>') {
				break;
			}
			buf.append(c);
		}
	}
	return make_pair(done, readi);
}

bool FastaPatternSource::parse(Read& r, Read& rb) const {
	// We assume the light parser has put the raw data for the separate ends
	// into separate Read objects.  That doesn't have to be the case, but
	// that's how we've chosen to do it for FastqPatternSource
	assert(!r.readOrigBuf.empty());
	assert(r.empty());
	int c;
	size_t cur = 1;
	
	// Parse read name
	assert(r.name.empty());
	while(true) {
		assert(cur < r.readOrigBuf.length());
		c = r.readOrigBuf[cur++];
		if(c == '\n' || c == '\r') {
			do {
				c = r.readOrigBuf[cur++];
			} while(c == '\n' || c == '\r');
			break;
		}
		r.name.append(c);
	}
	
	// Parse sequence
	int nchar = 0;
	assert(r.patFw.empty());
	while(c != '\n') {
		if(c == '.') {
			c = 'N';
		}
		if(isalpha(c)) {
			// If it's past the 5'-end trim point
			if(nchar++ >= gTrim5) {
				r.patFw.append(asc2dna[c]);
			}
		}
		assert(cur < r.readOrigBuf.length());
		c = r.readOrigBuf[cur++];
	}
	r.trimmed5 = (int)(nchar - r.patFw.length());
	r.trimmed3 = (int)(r.patFw.trimEnd(gTrim3));
	
	for(size_t i = 0; i < r.patFw.length(); i++) {
		r.qual.append('I');
	}

	// Set up a default name if one hasn't been set
	if(r.name.empty()) {
		char cbuf[20];
		itoa10<TReadId>(static_cast<TReadId>(readCnt_), cbuf);
		r.name.install(cbuf);
	}
	if(!rb.readOrigBuf.empty() && rb.patFw.empty()) {
		return parse(rb, r);
	}
	return true;
}

/**
 * "Light" parser.  This is inside the critical section, so the key is to do
 * just enough parsing so that another function downstream (finalize()) can do
 * the rest of the parsing.  Really this function's only job is to stick every
 * for lines worth of the input file into a buffer (r.readOrigBuf).  finalize()
 * then parses the contents of r.readOrigBuf later.
 */
pair<bool, int> FastaContinuousPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a)
{
	throw 1;
	return make_pair(false, 0);
}

/// Read another pattern from a FASTA input file
bool FastaContinuousPatternSource::parse(Read& r, Read& rb) const {
	assert(r.empty());
	assert(rb.empty());
#if 0
	r.reset();
	while(true) {
		int c = fb_.get();
		if(c < 0) { return make_pair(true, 0); }
		if(c == '>') {
			resetForNextFile();
			c = fb_.peek();
			bool sawSpace = false;
			while(c != '\n' && c != '\r') {
				if(!sawSpace) {
					sawSpace = isspace(c);
				}
				if(!sawSpace) {
					nameBuf_.append(c);
				}
				fb_.get();
				c = fb_.peek();
			}
			while(c == '\n' || c == '\r') {
				fb_.get();
				c = fb_.peek();
			}
			nameBuf_.append('_');
		} else {
			int cat = asc2dnacat[c];
			if(cat >= 2) c = 'N';
			if(cat == 0) {
				// Encountered non-DNA, non-IUPAC char; skip it
				continue;
			} else {
				// DNA char
				buf_[bufCur_++] = c;
				if(bufCur_ == 1024) bufCur_ = 0;
				if(eat_ > 0) {
					eat_--;
					// Try to keep readCnt_ aligned with the offset
					// into the reference; that lets us see where
					// the sampling gaps are by looking at the read
					// name
					if(!beginning_) readCnt_++;
					continue;
				}
				for(size_t i = 0; i < length_; i++) {
					if(length_ - i <= bufCur_) {
						c = buf_[bufCur_ - (length_ - i)];
					} else {
						// Rotate
						c = buf_[bufCur_ - (length_ - i) + 1024];
					}
					r.patFw.append(asc2dna[c]);
					r.qual.append('I');
				}
				// Set up a default name if one hasn't been set
				r.name = nameBuf_;
				char cbuf[20];
				itoa10<TReadId>(readCnt_ - subReadCnt_, cbuf);
				r.name.append(cbuf);
				eat_ = freq_-1;
				readCnt_++;
				beginning_ = false;
				rdid = readCnt_-1;
				break;
			}
		}
	}
	return make_pair(false, 1);
#endif
	throw 1;
	return false;
}


/**
 * "Light" parser.  This is inside the critical section, so the key is to do
 * just enough parsing so that another function downstream (finalize()) can do
 * the rest of the parsing.  Really this function's only job is to stick every
 * for lines worth of the input file into a buffer (r.readOrigBuf).  finalize()
 * then parses the contents of r.readOrigBuf later.
 */
pair<bool, int> FastqPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a)
{
	int c;
	EList<Read>& readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	if(first_) {
		c = getc_unlocked(fp_);
		while(c == '\r' || c == '\n') {
			c = getc_unlocked(fp_);
		}
		if(c != '@') {
			cerr << "Error: reads file does not look like a FASTQ file" << endl;
			throw 1;
		}
		first_ = false;
		readbuf[0].readOrigBuf.append('@');
	}
	bool done = false, aborted = false;
	size_t readi = 0;
	// Read until we run out of input or until we've filled the buffer
	for(; readi < pt.max_buf_ && !done; readi++) {
		SStringExpandable<char, 1024, 2, 1024>& buf = readbuf[readi].readOrigBuf;
		assert(readi == 0 || buf.empty());
		int newlines = 4;
		while(newlines) {
			c = getc_unlocked(fp_);
			done = c < 0;
			if(c == '\n' || (done && newlines == 1)) {
				// Saw newline, or EOF that we're
				// interpreting as final newline
				newlines--;
				c = '\n';
			} else if(done) {
				aborted = true; // Unexpected EOF
				break;
			}
			buf.append(c);
		}
	}
	if(aborted) {
		readi--;
	}
	return make_pair(done, readi);
}

/// Read another pattern from a FASTQ input file
bool FastqPatternSource::parse(Read &r, Read& rb) const {
	// We assume the light parser has put the raw data for the separate ends
	// into separate Read objects.  That doesn't have to be the case, but
	// that's how we've chosen to do it for FastqPatternSource
	assert(!r.readOrigBuf.empty());
	assert(r.empty());
	int c;
	size_t cur = 1;

	// Parse read name
	assert(r.name.empty());
	while(true) {
		assert(cur < r.readOrigBuf.length());
		c = r.readOrigBuf[cur++];
		if(c == '\n' || c == '\r') {
			do {
				c = r.readOrigBuf[cur++];
			} while(c == '\n' || c == '\r');
			break;
		}
		r.name.append(c);
	}
	
	// Parse sequence
	int nchar = 0;
	assert(r.patFw.empty());
	while(c != '+') {
		if(c == '.') {
			c = 'N';
		}
		if(isalpha(c)) {
			// If it's past the 5'-end trim point
			if(nchar++ >= gTrim5) {
				r.patFw.append(asc2dna[c]);
			}
		}
		assert(cur < r.readOrigBuf.length());
		c = r.readOrigBuf[cur++];
	}
	r.trimmed5 = (int)(nchar - r.patFw.length());
	r.trimmed3 = (int)(r.patFw.trimEnd(gTrim3));
	
	assert_eq('+', c);
	do {
		assert(cur < r.readOrigBuf.length());
		c = r.readOrigBuf[cur++];
	} while(c != '\n' && c != '\r');
	do {
		assert(cur < r.readOrigBuf.length());
		c = r.readOrigBuf[cur++];
	} while(c == '\n' || c == '\r');
	
	// Now we're on the next non-blank line after the + line
	if(r.patFw.empty()) {
		return true; // done parsing empty read
	}

	assert(r.qual.empty());
	int nqual = 0;
	if (intQuals_) {
		throw 1; // not yet implemented
	} else {
		c = charToPhred33(c, solQuals_, phred64Quals_);
		if(nqual++ >= r.trimmed5) {
			r.qual.append(c);
		}
		while(cur < r.readOrigBuf.length()) {
			c = r.readOrigBuf[cur++];
			if (c == ' ') {
				wrongQualityFormat(r.name);
				return false;
			}
			if(c == '\r' || c == '\n') {
				break;
			}
			c = charToPhred33(c, solQuals_, phred64Quals_);
			if(nqual++ >= r.trimmed5) {
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
	// Set up a default name if one hasn't been set
	if(r.name.empty()) {
		char cbuf[20];
		itoa10<TReadId>(static_cast<TReadId>(readCnt_), cbuf);
		r.name.install(cbuf);
	}
	if(!rb.readOrigBuf.empty() && rb.patFw.empty()) {
		return parse(rb, r);
	}
	return true;
}

/// Read another pattern from a FASTA input file
pair<bool, int> TabbedPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a)
{
	bool success = true;
	int c;
	EList<Read>& readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	if(first_) {
		c = fb_.get();
		while(c == '\r' || c == '\n') {
			c = fb_.get();
		}
		if(c != '>') {
			cerr << "Error: reads file does not look like a FASTQ file" << endl;
			throw 1;
		}
		first_ = false;
	}
	bool done = false;
	size_t readi = 0;
	// Read until we run out of input or until we've filled the buffer
	for(; readi < pt.max_buf_ && !done; readi++) {
		SStringExpandable<char, 1024, 2, 1024>& buf = readbuf[readi].readOrigBuf;
		buf.clear();
		buf.append('>'); // TODO: need to handle first char differently
		while(true) {
			c = fb_.get();
			if(c < 0 || c == '>') {
				done = true;
				break;
			}
			buf.append(c);
		}
	}
	return make_pair(success, readi);
}

bool TabbedPatternSource::parse(Read& r, Read& rb) const {
	r.reset();
#if 0
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
		return make_pair(true, 0);
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
		return make_pair(true, 0);
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
		return make_pair(true, 0);
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
		rdid = readCnt_;
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
			return make_pair(true, 0);;
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
		return return make_pair(true, 0);;
	}
	assert_neq('\t', fb_.peek());

	// fb_ is about to dish out the first character of the
	// quality-string field
	if(parseQuals(rb, charsRead2, dstLen2, mytrim5_2, ct, '\n') < 0) {
		peekOverNewline(fb_); // skip rest of line
		ra.reset();
		rb.reset();
		fb_.resetLastN();
		return make_pair(true, 0);;
	}
	ra.readOrigBuf.install(fb_.lastN(), fb_.lastNLen());
	fb_.resetLastN();
	rb.trimmed3 = gTrim3;
	rb.trimmed5 = mytrim5_2;
	rdid = readCnt_;
	readCnt_++;
	return make_pair(false, 1);
#endif
	cerr << "In TabbedPatternSource.parse()" << endl;
	throw 1;
	return false;
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
		itoa10<TReadId>(static_cast<TReadId>(readCnt_), cbuf);
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
	while(c != upto) {
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
		if (charsRead > qualsRead) {
			tooFewQualities(r.name);
		}
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

/**
 * Light-parse a batch into the given buffer.
 */
pair<bool, int> RawPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a)
{
	bool success = true;
	int c;
	EList<Read>& readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	if(first_) {
		c = fb_.get();
		while(c == '\r' || c == '\n') {
			c = fb_.get();
		}
		if(c != '>') {
			cerr << "Error: reads file does not look like a FASTQ file" << endl;
			throw 1;
		}
		first_ = false;
	}
	bool done = false;
	size_t readi = 0;
	// Read until we run out of input or until we've filled the buffer
	for(; readi < pt.max_buf_ && !done; readi++) {
		SStringExpandable<char, 1024, 2, 1024>& buf = readbuf[readi].readOrigBuf;
		buf.clear();
		buf.append('>'); // TODO: need to handle first char differently
		while(true) {
			c = fb_.get();
			if(c < 0 || c == '>') {
				done = true;
				break;
			}
			buf.append(c);
		}
	}
	return make_pair(success, readi);
}

/// Skip to the end of the current string of newline chars and return
/// the first character after the newline chars, or -1 for EOF
static inline int getOverNewline(FileBuf& in) {
	int c;
	while(isspace(c = in.get()));
	return c;
}

/// Skip to the end of the current line such that the next call to
/// get() returns the first character on the next line
static inline int peekToEndOfLine(FileBuf& in) {
	while(true) {
		int c = in.get(); if(c < 0) return c;
		if(c == '\n' || c == '\r') {
			c = in.peek();
			while(c == '\n' || c == '\r') {
				in.get(); if(c < 0) return c; // consume \r or \n
				c = in.peek();
			}
			// next get() gets first character of next line
			return c;
		}
	}
}

/// Read another pattern from a Raw input file
bool RawPatternSource::parse(Read& r, Read& rb) const {
#if 0
	int c;
	r.reset();
	c = getOverNewline(this->fb_);
	if(c < 0) {
		bail(r);
		return make_pair(true, 0);
	}
	assert(!isspace(c));
	int mytrim5 = gTrim5;
	if(first_) {
		// Check that the first character is sane for a raw file
		int cc = c;
		if(asc2dnacat[cc] == 0) {
			cerr << "Error: reads file does not look like a Raw file" << endl;
			if(c == '>') {
				cerr << "Reads file looks like a FASTA file; please use -f" << endl;
			}
			if(c == '@') {
				cerr << "Reads file looks like a FASTQ file; please use -q" << endl;
			}
			throw 1;
		}
		first_ = false;
	}
	// _in now points just past the first character of a sequence
	// line, and c holds the first character
	int chs = 0;
	while(!isspace(c) && c >= 0) {
		// 5' trimming
		if(isalpha(c) && chs >= mytrim5) {
			//size_t len = chs - mytrim5;
			//if(len >= 1024) tooManyQualities(BTString("(no name)"));
			r.patFw.append(asc2dna[c]);
			r.qual.append('I');
		}
		chs++;
		if(isspace(fb_.peek())) break;
		c = fb_.get();
	}
	// 3' trimming
	r.patFw.trimEnd(gTrim3);
	r.qual.trimEnd(gTrim3);
	c = peekToEndOfLine(fb_);
	r.trimmed3 = gTrim3;
	r.trimmed5 = mytrim5;
	r.readOrigBuf.install(fb_.lastN(), fb_.lastNLen());
	fb_.resetLastN();
	
	// Set up name
	char cbuf[20];
	itoa10<TReadId>(readCnt_, cbuf);
	r.name.install(cbuf);
	readCnt_++;
	
	rdid = readCnt_-1;
	return make_pair(false, 1);
#endif
	cerr << "In RawPatternSource.parse()" << endl;
	throw 1;
	return false;
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
