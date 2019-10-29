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
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <string.h>
#include <fcntl.h>
#include "sstring.h"

#include "pat.h"
#include "filebuf.h"
#include "formats.h"
#include "util.h"
#include "tokenize.h"

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
		rseed ^= ((uint32_t)p << off);
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
		case FASTA:       return new FastaPatternSource(qs, p, p.interleaved);
		case FASTA_CONT:  return new FastaContinuousPatternSource(qs, p);
		case RAW:         return new RawPatternSource(qs, p);
		case FASTQ:       return new FastqPatternSource(qs, p, p.interleaved);
		case BAM:         return new BAMPatternSource(qs, p);
		case TAB_MATE5:   return new TabbedPatternSource(qs, p, false);
		case TAB_MATE6:   return new TabbedPatternSource(qs, p, true);
		case CMDLINE:     return new VectorPatternSource(qs, p);
		case QSEQ:        return new QseqPatternSource(qs, p);
#ifdef USE_SRA
		case SRA_FASTA:
		case SRA_FASTQ:   return new SRAPatternSource(qs, p);
#endif
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
		pair<bool, int> res = nextBatch(); // more parsing needed!
		if(res.first && res.second == 0) {
			return make_pair(false, true);
		}
		last_batch_ = res.first;
		last_batch_size_ = res.second;
		assert_eq(0, buf_.cur_buf_);
	} else {
		buf_.next(); // advance cursor; no parsing or locking needed
		assert_gt(buf_.cur_buf_, 0);
	}
	// Now fully parse read/pair *outside* the critical section
	assert(!buf_.read_a().readOrigBuf.empty());
	assert(buf_.read_a().empty());
	if(!parse(buf_.read_a(), buf_.read_b())) {
		return make_pair(false, false);
	}
	// Finalize read/pair
	if(!buf_.read_b().patFw.empty()) {
		trim(buf_.read_a());
		trim(buf_.read_b());
		finalizePair(buf_.read_a(), buf_.read_b());
	} else {
		trim(buf_.read_a());
		finalize(buf_.read_a());
	}
	bool this_is_last = buf_.cur_buf_ == static_cast<unsigned int>(last_batch_size_-1);
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
			ThreadSafe ts(mutex_m);
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
			if(res.second == 0 && cur < srca_->size() - 1) {
				ThreadSafe ts(mutex_m);
				if(cur + 1 > cur_) cur_++;
				cur = cur_; // Move on to next PatternSource
				continue; // on to next pair of PatternSources
			}
			return make_pair(res.first && cur == srca_->size() - 1, res.second);
		} else {
			pair<bool, int> resa, resb;
			// Lock to ensure that this thread gets parallel reads
			// in the two mate files
			{
				ThreadSafe ts(mutex_m);
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
				ThreadSafe ts(mutex_m);
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
			return make_pair(resa.first && cur == srca_->size() - 1, resa.second);
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
#ifdef USE_SRA
	const EList<string>& sra_accs, // SRA accessions
#endif
	PatternParams& p,    // read-in parameters
	bool verbose)              // be talkative?
{
	EList<PatternSource*>* a  = new EList<PatternSource*>();
	EList<PatternSource*>* b  = new EList<PatternSource*>();
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
		a->push_back(PatternSource::patsrcFromStrings(p, *qs));
		b->push_back(NULL);
		p.interleaved = false;
		if(!p.fileParallel) {
			break;
		}
	}

#ifdef USE_SRA
	for(size_t i = 0; i < sra_accs.size(); i++) {
		const EList<string>* qs = &sra_accs;
		EList<string> tmp;
		if(p.fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmp;
			tmp.push_back(sra_accs[i]);
			assert_eq(1, tmp.size());
		}
		a->push_back(PatternSource::patsrcFromStrings(p, *qs));
		b->push_back(NULL);
		if(!p.fileParallel) {
			break;
		}
	}
#endif

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
	patsrc = new DualPatternComposer(a, b, p);
	return patsrc;
}

void PatternComposer::free_EList_pmembers( const EList<PatternSource*> &elist) {
	for (size_t i = 0; i < elist.size(); i++)
		if (elist[i] != NULL)
			delete elist[i];
}

/**
 * Fill Read with the sequence, quality and name for the next
 * read in the list of read files. This function gets called by
 * all the search threads, so we must handle synchronization.
 *
 * Returns pair<bool, int> where bool indicates whether we're
 * completely done, and int indicates how many reads were read.
 */
pair<bool, int> CFilePatternSource::nextBatchImpl(
	PerThreadReadBuf& pt,
	bool batch_a)
{
	bool done = false;
	unsigned nread = 0;
	pt.setReadId(readCnt_);
	while(true) { // loop that moves on to next file when needed
		do {
			pair<bool, int> ret = nextBatchFromFile(pt, batch_a, nread);
			done = ret.first;
			nread = ret.second;
		} while(!done && nread == 0); // not sure why this would happen
		if(done && filecur_ < infiles_.size()) { // finished with this file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			filecur_++;
			if(nread == 0 || (nread < pt.max_buf_)) {
				continue;
			}
			done = false;
		}
		break;
	}
	assert_geq(nread, 0);
	readCnt_ += nread;
	return make_pair(done, nread);
}

pair<bool, int> CFilePatternSource::nextBatch(
	PerThreadReadBuf& pt,
	bool batch_a,
	bool lock)
{
	if(lock) {
		// synchronization at this level because both reading and manipulation of
		// current file pointer have to be protected
		ThreadSafe ts(mutex);
		return nextBatchImpl(pt, batch_a);
	} else {
		return nextBatchImpl(pt, batch_a);
	}
}

/**
 * Open the next file in the list of input files.
 */
void CFilePatternSource::open() {
	if(is_open_) {
		is_open_ = false;
		if (compressed_) {
			gzclose(zfp_);
			zfp_ = NULL;
		}
		else {
			fclose(fp_);
      			fp_ = NULL;
      		}
	}
	while(filecur_ < infiles_.size()) {
		if(infiles_[filecur_] == "-") {
			// always assume that data from stdin is compressed
			compressed_ = true;
			int fd = dup(fileno(stdin));
			zfp_ = gzdopen(fd, "rb");

			if (zfp_ == NULL) {
				close(fd);
			}
		}
		else {
			const char* filename = infiles_[filecur_].c_str();

			int fd = ::open(filename, O_RDONLY);
			bool is_fifo = false;

#ifndef _WIN32
			struct stat st;
			if (fstat(fd, &st) != 0) {
				perror("stat");
			}

			is_fifo = S_ISFIFO(st.st_mode) != 0;
#endif
			if (pp_.format != BAM && (is_fifo || is_gzipped_file(fd))) {
				zfp_ = gzdopen(fd, "r");
				compressed_ = true;
			} else {
				fp_ = fdopen(fd, "rb");
			}

			if((compressed_ && zfp_ == NULL) || (!compressed_ && fp_ == NULL)) {
				if (fd != -1) {
					close(fd);
				}

				if(!errs_[filecur_]) {
					cerr << "Warning: Could not open read file \""
					     << filename
					     << "\" for reading; skipping..." << endl;
					errs_[filecur_] = true;
				}
				filecur_++;
				continue;
			}
		}
		is_open_ = true;
		if (compressed_) {
#if ZLIB_VERNUM < 0x1235
			cerr << "Warning: gzbuffer added in zlib v1.2.3.5. Unable to change "
				"buffer size from default of 8192." << endl;
#else
			gzbuffer(zfp_, 128*1024);
#endif
		}
		else {
			setvbuf(fp_, buf_, _IOFBF, 64*1024);
		}
		return;
	}
	cerr << "Error: No input read files were valid" << endl;
	exit(1);
	return;
}

/**
 * Constructor for vector pattern source, used when the user has
 * specified the input strings on the command line using the -c
 * option.
 */
VectorPatternSource::VectorPatternSource(
	const EList<string>& seqs,
	const PatternParams& p) :
	PatternSource(p),
	cur_(p.skip),
	skip_(p.skip),
	paired_(false),
	tokbuf_(),
	bufs_()
{
	// Install sequences in buffers, ready for immediate copying in
	// nextBatch().  Formatting of the buffer is just like
	// TabbedPatternSource.
	const size_t seqslen = seqs.size();
	for(size_t i = 0; i < seqslen; i++) {
		tokbuf_.clear();
		tokenize(seqs[i], ":", tokbuf_, 2);
		assert_gt(tokbuf_.size(), 0);
		assert_leq(tokbuf_.size(), 2);
		// Get another buffer ready
		bufs_.expand();
		bufs_.back().clear();
		// Install name
		itoa10<TReadId>(static_cast<TReadId>(i), nametmp_);
		bufs_.back().install(nametmp_);
		bufs_.back().append('\t');
		// Install sequence
		bufs_.back().append(tokbuf_[0].c_str());
		bufs_.back().append('\t');
		// Install qualities
		if(tokbuf_.size() > 1) {
			bufs_.back().append(tokbuf_[1].c_str());
		} else {
			const size_t len = tokbuf_[0].length();
			for(size_t i = 0; i < len; i++) {
				bufs_.back().append('I');
			}
		}
	}
}

/**
 * Read next batch.  However, batch concept is not very applicable for this
 * PatternSource where all the info has already been parsed into the fields
 * in the contsructor.	This essentially modifies the pt as though we read
 * in some number of patterns.
 */
pair<bool, int> VectorPatternSource::nextBatchImpl(
	PerThreadReadBuf& pt,
	bool batch_a)
{
	pt.setReadId(cur_);
	EList<Read>& readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	size_t readi = 0;
	for(; readi < pt.max_buf_ && cur_ < bufs_.size(); readi++, cur_++) {
		readbuf[readi].readOrigBuf = bufs_[cur_];
	}
	readCnt_ += readi;
	return make_pair(cur_ == bufs_.size(), readi);
}

pair<bool, int> VectorPatternSource::nextBatch(
	PerThreadReadBuf& pt,
	bool batch_a,
	bool lock)
{
	if(lock) {
		ThreadSafe ts(mutex);
		return nextBatchImpl(pt, batch_a);
	} else {
		return nextBatchImpl(pt, batch_a);
	}
}

/**
 * Finishes parsing outside the critical section.
 */
bool VectorPatternSource::parse(Read& ra, Read& rb, TReadId rdid) const {
	// Very similar to TabbedPatternSource

	// Light parser (nextBatchFromFile) puts unparsed data
	// into Read& r, even when the read is paired.
	assert(ra.empty());
	assert(!ra.readOrigBuf.empty()); // raw data for read/pair is here
	int c = '\t';
	size_t cur = 0;
	const size_t buflen = ra.readOrigBuf.length();
	
	// Loop over the two ends
	for(int endi = 0; endi < 2 && c == '\t'; endi++) {
		Read& r = ((endi == 0) ? ra : rb);
		assert(r.name.empty());
		// Parse name if (a) this is the first end, or
		// (b) this is tab6
		if(endi < 1 || paired_) {
			// Parse read name
			c = ra.readOrigBuf[cur++];
			while(c != '\t' && cur < buflen) {
				r.name.append(c);
				c = ra.readOrigBuf[cur++];
			}
			assert_eq('\t', c);
			if(cur >= buflen) {
				return false; // record ended prematurely
			}
		} else if(endi > 0) {
			// if this is the second end and we're parsing
			// tab5, copy name from first end
			rb.name = ra.name;
		}

		// Parse sequence
		assert(r.patFw.empty());
		c = ra.readOrigBuf[cur++];
		int nchar = 0;
		while(c != '\t' && cur < buflen) {
			if(isalpha(c)) {
				assert_in(toupper(c), "ACGTN");
				if(nchar++ >= pp_.trim5) {
					assert_neq(0, asc2dnacat[c]);
					r.patFw.append(asc2dna[c]); // ascii to int
				}
			}
			c = ra.readOrigBuf[cur++];
		}
		assert_eq('\t', c);
		if(cur >= buflen) {
			return false; // record ended prematurely
		}
		// record amt trimmed from 5' end due to --trim5
		r.trimmed5 = (int)(nchar - r.patFw.length());
		// record amt trimmed from 3' end due to --trim3
		r.trimmed3 = (int)(r.patFw.trimEnd(pp_.trim3));
		
		// Parse qualities
		assert(r.qual.empty());
		c = ra.readOrigBuf[cur++];
		int nqual = 0;
		while(c != '\t' && c != '\n' && c != '\r') {
			if(c == ' ') {
				wrongQualityFormat(r.name);
				return false;
			}
			char cadd = charToPhred33(c, false, false);
			if(++nqual > pp_.trim5) {
				r.qual.append(cadd);
			}
			if(cur >= buflen) break;
			c = ra.readOrigBuf[cur++];
		}
		if(nchar > nqual) {
			tooFewQualities(r.name);
			return false;
		} else if(nqual > nchar) {
			tooManyQualities(r.name);
			return false;
		}
		r.qual.trimEnd(pp_.trim3);
		assert(c == '\t' || c == '\n' || c == '\r' || cur >= buflen);
		assert_eq(r.patFw.length(), r.qual.length());
	}
	ra.parsed = true;
	if(!rb.parsed && !rb.readOrigBuf.empty()) {
		return parse(rb, ra, rdid);
	}
	return true;
}

/**
 * Light-parse a FASTA batch into the given buffer.
 */
pair<bool, int> FastaPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a, unsigned readi)
{
	int c;
	EList<Read>* readbuf = batch_a ? &pt.bufa_ : &pt.bufb_;
	if(first_) {
		c = getc_wrapper();
		if (c == EOF) {
			return make_pair(true, 0);
		}
		while(c == '\r' || c == '\n') {
			c = getc_wrapper();
		}
		if(c != '>') {
			cerr << "Error: reads file does not look like a FASTA file" << endl;
			throw 1;
		}
		first_ = false;
	}
	bool done = false;
	// Read until we run out of input or until we've filled the buffer
	while (readi < pt.max_buf_ && !done) {
		Read::TBuf& buf = (*readbuf)[readi].readOrigBuf;
		buf.clear();
		buf.append('>');
		while(true) {
			c = getc_wrapper();
			if(c < 0 || c == '>') {
				done = c < 0;
				break;
			}
			buf.append(c);
		}
		if (interleaved_) {
			// alternate between read buffers
			batch_a = !batch_a;
			readbuf = batch_a ? &pt.bufa_ : &pt.bufb_;
			// increment read counter after each pair gets read
			readi = batch_a ? readi+1 : readi;
		} else {
			readi++;
		}

	}
	// Immediate EOF case
	if(done && (*readbuf)[readi-1].readOrigBuf.length() == 1) {
		readi--;
	}
	return make_pair(done, readi);
}

/**
 * Finalize FASTA parsing outside critical section.
 */
bool FastaPatternSource::parse(Read& r, Read& rb, TReadId rdid) const {
	// We assume the light parser has put the raw data for the separate ends
	// into separate Read objects.	That doesn't have to be the case, but
	// that's how we've chosen to do it for FastqPatternSource
	assert(!r.readOrigBuf.empty());
	assert(r.empty());
	int c = -1;
	size_t cur = 1;
	const size_t buflen = r.readOrigBuf.length();

	// Parse read name
	assert(r.name.empty());
	while(cur < buflen) {
		c = r.readOrigBuf[cur++];
		if(c == '\n' || c == '\r') {
			do {
				c = r.readOrigBuf[cur++];
			} while((c == '\n' || c == '\r') && cur < buflen);
			break;
		}
		r.name.append(c);
	}
	if(cur >= buflen) {
		return false; // FASTA ended prematurely
	}

	// Parse sequence
	int nchar = 0;
	assert(r.patFw.empty());
	assert(c != '\n' && c != '\r');
	assert_lt(cur, buflen);
	while(cur < buflen) {
		if(c == '.') {
			c = 'N';
		}
		if(isalpha(c)) {
			// If it's past the 5'-end trim point
			if(nchar++ >= pp_.trim5) {
				r.patFw.append(asc2dna[c]);
			}
		}
		assert_lt(cur, buflen);
		c = r.readOrigBuf[cur++];
		if ((c == '\n' || c == '\r')
				&& cur < buflen
				&& r.readOrigBuf[cur] != '>') {
			c = r.readOrigBuf[cur++];
		}
	}
	r.trimmed5 = (int)(nchar - r.patFw.length());
	r.trimmed3 = (int)(r.patFw.trimEnd(pp_.trim3));

	for(size_t i = 0; i < r.patFw.length(); i++) {
		r.qual.append('I');
	}

	// Set up a default name if one hasn't been set
	if(r.name.empty()) {
		char cbuf[20];
		itoa10<TReadId>(static_cast<TReadId>(rdid), cbuf);
		r.name.install(cbuf);
	}
	r.parsed = true;
	if(!rb.parsed && !rb.readOrigBuf.empty()) {
		return parse(rb, r, rdid);
	}
	return true;
}

/**
 * Light-parse a FASTA-continuous batch into the given buffer.
 * This is trickier for FASTA-continuous than for other formats,
 * for several reasons:
 *
 * 1. Reads are substrings of a longer FASTA input string
 * 2. Reads may overlap w/r/t the longer FASTA string
 * 3. Read names depend on the most recently observed FASTA
 *	  record name
 */
pair<bool, int> FastaContinuousPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a, unsigned readi)
{
	int c = -1;
	EList<Read>& readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	while(readi < pt.max_buf_) {
		c = getc_wrapper();
		if(c < 0) {
			break;
		}
		if(c == '>') {
			resetForNextFile();
			c = getc_wrapper();
			bool sawSpace = false;
			while(c != '\n' && c != '\r') {
				if(!sawSpace) {
					sawSpace = isspace(c);
				}
				if(!sawSpace) {
					name_prefix_buf_.append(c);
				}
				c = getc_wrapper();
			}
			while(c == '\n' || c == '\r') {
				c = getc_wrapper();
			}
			if(c < 0) {
				break;
			}
			name_prefix_buf_.append('_');
		}
		int cat = asc2dnacat[c];
		if(cat >= 2) c = 'N';
		if(cat == 0) {
			// Non-DNA, non-IUPAC char; skip
			continue;
		} else {
			// DNA char
			buf_[bufCur_++] = c;
			if(bufCur_ == 1024) {
				bufCur_ = 0; // wrap around circular buf
			}
			if(eat_ > 0) {
				eat_--;
				// Try to keep cur_ aligned with the offset
				// into the reference; that lets us see where
				// the sampling gaps are by looking at the read
				// name
				if(!beginning_) {
					cur_++;
				}
				continue;
			}
			// install name
			readbuf[readi].readOrigBuf = name_prefix_buf_;
			itoa10<TReadId>(cur_ - last_, name_int_buf_);
			readbuf[readi].readOrigBuf.append(name_int_buf_);
			readbuf[readi].readOrigBuf.append('\t');
			// install sequence
			for(size_t i = 0; i < length_; i++) {
				if(length_ - i <= bufCur_) {
					c = buf_[bufCur_ - (length_ - i)];
				} else {
					// Rotate
					c = buf_[bufCur_ - (length_ - i) + 1024];
				}
				readbuf[readi].readOrigBuf.append(c);
			}
			eat_ = freq_-1;
			cur_++;
			beginning_ = false;
			readi++;
		}
	}
	return make_pair(c < 0, readi);
}

/**
 * Finalize FASTA-continuous parsing outside critical section.
 */
bool FastaContinuousPatternSource::parse(
	Read& ra,
	Read& rb,
	TReadId rdid) const
{
	// Light parser (nextBatchFromFile) puts unparsed data
	// into Read& r, even when the read is paired.
	assert(ra.empty());
	assert(rb.empty());
	assert(!ra.readOrigBuf.empty()); // raw data for read/pair is here
	assert(rb.readOrigBuf.empty());
	int c = '\t';
	size_t cur = 0;
	const size_t buflen = ra.readOrigBuf.length();
	
	// Parse read name
	c = ra.readOrigBuf[cur++];
	while(c != '\t' && cur < buflen) {
		ra.name.append(c);
		c = ra.readOrigBuf[cur++];
	}
	assert_eq('\t', c);
	if(cur >= buflen) {
		return false; // record ended prematurely
	}

	// Parse sequence
	assert(ra.patFw.empty());
	int nchar = 0;
	while(cur < buflen) {
		c = ra.readOrigBuf[cur++];
		if(isalpha(c)) {
			assert_in(toupper(c), "ACGTN");
			if(nchar++ >= pp_.trim5) {
				assert_neq(0, asc2dnacat[c]);
				ra.patFw.append(asc2dna[c]); // ascii to int
			}
		}
	}
	// record amt trimmed from 5' end due to --trim5
	ra.trimmed5 = (int)(nchar - ra.patFw.length());
	// record amt trimmed from 3' end due to --trim3
	ra.trimmed3 = (int)(ra.patFw.trimEnd(pp_.trim3));
	
	// Make fake qualities
	assert(ra.qual.empty());
	const size_t len = ra.patFw.length();
	for(size_t i = 0; i < len; i++) {
		ra.qual.append('I');
	}
	return true;
}


/**
 * "Light" parser. This is inside the critical section, so the key is to do
 * just enough parsing so that another function downstream (finalize()) can do
 * the rest of the parsing.  Really this function's only job is to stick every
 * for lines worth of the input file into a buffer (r.readOrigBuf).  finalize()
 * then parses the contents of r.readOrigBuf later.
 */
pair<bool, int> FastqPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a, unsigned readi)
{
	int c = -1;
	EList<Read>* readbuf = batch_a ? &pt.bufa_ : &pt.bufb_;
	if(first_) {
		c = getc_wrapper();
		if (c == EOF) {
			return make_pair(true, 0);
		}
		while(c == '\r' || c == '\n') {
			c = getc_wrapper();
		}
		if(c != '@') {
			cerr << "Error: reads file does not look like a FASTQ file" << endl;
			throw 1;
		}
		first_ = false;
		(*readbuf)[readi].readOrigBuf.append('@');
	}

	bool done = false, aborted = false;
	// Read until we run out of input or until we've filled the buffer
	while (readi < pt.max_buf_ && !done) {
		Read::TBuf& buf = (*readbuf)[readi].readOrigBuf;
		int newlines = 4;
		while(newlines) {
			c = getc_wrapper();
			done = c < 0;
			if(c == '\n' || (done && newlines == 1)) {
				// Saw newline, or EOF that we're
				// interpreting as final newline
				newlines--;
				c = '\n';
			} else if(done) {
				// account for newline at the end of the file
				if (newlines == 4) {
					newlines = 0;
				}
				else {
					aborted = true; // Unexpected EOF
				}
				break;
			}
			buf.append(c);
		}
		if (c > 0) {
			if (interleaved_) {
				// alternate between read buffers
				batch_a = !batch_a;
				readbuf = batch_a ? &pt.bufa_ : &pt.bufb_;
				// increment read counter after each pair gets read
				readi = batch_a ? readi+1 : readi;
			}
			else {
				readi++;
			}
		}
	}
	if(aborted) {
		readi--;
	}
	return make_pair(done, readi);
}

/**
 * Finalize FASTQ parsing outside critical section.
 */
bool FastqPatternSource::parse(Read &r, Read& rb, TReadId rdid) const {
	// We assume the light parser has put the raw data for the separate ends
	// into separate Read objects. That doesn't have to be the case, but
	// that's how we've chosen to do it for FastqPatternSource
	assert(!r.readOrigBuf.empty());
	assert(r.empty());
	int c;
	size_t cur = 1;
	const size_t buflen = r.readOrigBuf.length();

	// Parse read name
	assert(r.name.empty());
	while(true) {
		assert_lt(cur, buflen);
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
			if(nchar++ >= pp_.trim5) {
				r.patFw.append(asc2dna[c]);
			}
		}
		assert(cur < r.readOrigBuf.length());
		c = r.readOrigBuf[cur++];
	}
	r.trimmed5 = (int)(nchar - r.patFw.length());
	r.trimmed3 = (int)(r.patFw.trimEnd(pp_.trim3));
	
	assert_eq('+', c);
	do {
		assert(cur < r.readOrigBuf.length());
		c = r.readOrigBuf[cur++];
	} while(c != '\n' && c != '\r');
	while(cur < buflen && (c == '\n' || c == '\r')) {
		c = r.readOrigBuf[cur++];
	}
	
	assert(r.qual.empty());
	if(nchar > 0) {
		int nqual = 0;
		if (pp_.intQuals) {
			int cur_int = 0;
			while(c != '\t' && c != '\n' && c != '\r') {
				cur_int *= 10;
				cur_int += (int)(c - '0');
				c = r.readOrigBuf[cur++];
				if(c == ' ' || c == '\t' || c == '\n' || c == '\r') {
					char cadd = intToPhred33(cur_int, pp_.solexa64);
					cur_int = 0;
					assert_geq(cadd, 33);
					if(++nqual > pp_.trim5) {
						r.qual.append(cadd);
					}
				}
			}
		} else {
			c = charToPhred33(c, pp_.solexa64, pp_.phred64);
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
				c = charToPhred33(c, pp_.solexa64, pp_.phred64);
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
	}
	// Set up a default name if one hasn't been set
	if(r.name.empty()) {
		char cbuf[20];
		itoa10<TReadId>(static_cast<TReadId>(rdid), cbuf);
		r.name.install(cbuf);
	}
	r.parsed = true;
	if(!rb.parsed && !rb.readOrigBuf.empty()) {
		return parse(rb, r, rdid);
	}
	return true;
}

const int BAMPatternSource::offset[] = {
	0,   //refID
	4,   //pos
	8,   //l_read_name
	9,   //mapq
	10,  //bin
	12,  //n_cigar_op
	14,  //flag
	16,  //l_seq
	20,  //next_refID
	24,  //next_pos
	28,  //tlen
	32,  //read_name
};

const uint8_t BAMPatternSource::EOF_MARKER[] = {
	0x1f,  0x8b,  0x08,  0x04,  0x00,  0x00,  0x00,  0x00,  0x00,  0xff,
	0x06,  0x00,  0x42,  0x43,  0x02,  0x00,  0x1b,  0x00,  0x03,  0x00,
	0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00
};

bool BAMPatternSource::parse_bam_header() {
	char magic[4];

	if (zread(magic, 4) != 4 || strncmp(magic, "BAM\001", 4) != 0) {
		std::cerr << "This file is not a BAM file" << std::endl;
		return false;
	}

	int32_t l_text = 0;
	int32_t n_ref  = 0;
	char* data = NULL;
	int32_t size = 0;

	if (zread(&l_text, sizeof(l_text)) != sizeof(l_text)) {
		return false;
	}

	size = l_text + 1;
	data = new char[size];
	if (data == NULL) {
		return false;
	}
	if (zread(data, l_text) != l_text) {
		delete[] data;
		return false;
	}
	if (zread(&n_ref, sizeof(n_ref)) != sizeof(n_ref)) {
		delete[] data;
		return false;
	}

	for (int i = 0; i != n_ref; i++) {
		int32_t l_name, l_ref;
		if (zread(&l_name, sizeof(l_name)) != sizeof(l_name)) {
			delete[] data;
			return false;
		}
		if (l_name > size) {
			size = l_name;
			delete[] data;
			data = new char[size];
		}
		if (zread(data, l_name) != l_name) {
			delete[] data;
			return false;
		}
		if (zread(&l_ref, sizeof(l_ref)) != sizeof(l_ref)) {
			delete[] data;
			return false;
		}
	}

	delete[] data;

	return true;
}

uint16_t BAMPatternSource::nextBGZFBlockFromFile(BGZF& b) {
	size_t r = fread(&b.hdr, sizeof(b.hdr), 1, fp_);
	if ((r == 0 && ferror_unlocked(fp_)) || feof_unlocked(fp_)) {
		return 0;
	}

	uint8_t *extra = new uint8_t[b.hdr.xlen];
	r = fread(extra, b.hdr.xlen, 1, fp_);
	if (r == 0 && ferror_unlocked(fp_)) {
		return 0;
	}

	if (memcmp(EOF_MARKER, &b.hdr, sizeof(b.hdr)) == 0
			&& memcmp(EOF_MARKER + sizeof(b.hdr), extra, 6 /* sizeof BAM subfield */) == 0)
	{
		delete[] extra;
		fclose(fp_);
		return 0;
	}

	uint16_t bsize = 0;
	for (uint16_t i = 0; i < b.hdr.xlen;) {
		if (extra[0] == 66 && extra[1] == 67) {
			bsize = *((uint16_t *)(extra + 4));
			bsize -= (b.hdr.xlen + 19);
			break;
		}
		i = i + 2;
		uint16_t sub_field_len = *((uint16_t *)(extra + 2));
		i = i + 2 + sub_field_len;
	}

	delete[] extra;

	if (bsize == 0) {
		return 0;
	}

	r = fread(b.cdata, bsize, 1, fp_);
	if (r == 0 && ferror_unlocked(fp_)) {
		return 0;
	}

	r = fread(&b.ftr, sizeof b.ftr, 1, fp_);
	if (r == 0 && ferror_unlocked(fp_)) {
		return 0;
	}

	return bsize;
}

std::pair<bool, int> BAMPatternSource::nextBatch(PerThreadReadBuf& pt, bool batch_a, bool lock) {
	uint16_t cdata_len;

	unsigned nread = 0;
	bool done = false;

	do {
		if (bam_batch_indexes_[pt.tid_] >= bam_batches_[pt.tid_].size()) {
			BGZF& block = blocks_[pt.tid_];
			std::vector<uint8_t>& batch = bam_batches_[pt.tid_];
			if (lock) {
				ThreadSafe ts(mutex);
				if (first_) {
					nextBGZFBlockFromFile(block);
					first_ = false;
				}
				cdata_len = nextBGZFBlockFromFile(block);
			} else {
				cdata_len = nextBGZFBlockFromFile(block);
			}

			if (cdata_len == 0) {
				if (lock) {
					ThreadSafe ts(orphan_mates_mutex_);
					get_orphaned_pairs(pt.bufa_, pt.bufb_, pt.max_buf_, nread);
				} else {
					get_orphaned_pairs(pt.bufa_, pt.bufb_, pt.max_buf_, nread);
				}

				done = nread == 0;
				break;
			}

			bam_batch_indexes_[pt.tid_] = 0;

			batch.resize(block.ftr.isize);

			int ret_code = decompress_bgzf_block(&batch[0], block.ftr.isize, block.cdata, cdata_len);

			if (ret_code != Z_OK) {
				return make_pair(true, 0);
			}

			uLong crc = crc32(0L, Z_NULL, 0);
			crc = crc32(crc, &batch[0], batch.size());
			assert(crc == block.ftr.crc32);
		}

		std::pair<bool, int> ret = get_alignments(pt, batch_a, nread, lock);

		done = ret.first;
	} while (!done && nread < pt.max_buf_);

	if (lock) {
		ThreadSafe ts(mutex);
		pt.setReadId(readCnt_);
		readCnt_ += nread;
	} else {
		pt.setReadId(readCnt_);
		readCnt_ += nread;
	}

	return make_pair(done, nread);
}

std::pair<bool, int> BAMPatternSource::get_alignments(PerThreadReadBuf& pt, bool batch_a, unsigned& readi, bool lock) {
	size_t& i = bam_batch_indexes_[pt.tid_];
	bool done = false;
	bool read1 = true;

	while (readi < pt.max_buf_) {
		if (i >= bam_batches_[pt.tid_].size()) {
			return make_pair(false, readi);
		}

		uint16_t flag;
		int32_t block_size;
		EList<Read>& readbuf = pp_.align_paired_reads && !read1 ? pt.bufb_ : pt.bufa_;

		memcpy(&block_size, &bam_batches_[pt.tid_][0] + i, sizeof(block_size));
		if (block_size <= 0) {
			return make_pair(done, readi);
		}
		i += sizeof(block_size);

		memcpy(&flag, &bam_batches_[pt.tid_][0] + i + offset[BAMField::flag], sizeof(flag));
		if (!pp_.align_paired_reads && ((flag & 0x40) != 0 || (flag & 0x80) != 0)) {
			readbuf[readi].readOrigBuf.clear();
			i += block_size;
			continue;
		}

		if (pp_.align_paired_reads && ((flag & 0x40) == 0 && (flag & 0x80) == 0)) {
			readbuf[readi].readOrigBuf.clear();
			i += block_size;
			continue;
		}

		if (pp_.align_paired_reads && (((flag & 0x40) != 0
				&& i + block_size == bam_batches_[pt.tid_].size())
				|| ((flag & 0x80) != 0 && i == sizeof(block_size))))
		{
			if (lock) {
				ThreadSafe ts(orphan_mates_mutex_);
				store_orphan_mate(&bam_batches_[pt.tid_][0] + i, block_size);
				i += block_size;
				get_orphaned_pairs(pt.bufa_, pt.bufb_, pt.max_buf_, readi);
			} else {
				store_orphan_mate(&bam_batches_[pt.tid_][0] + i, block_size);
				i += block_size;
				get_orphaned_pairs(pt.bufa_, pt.bufb_, pt.max_buf_, readi);
			}

		} else {
			readbuf[readi].readOrigBuf.resize(block_size);

			memcpy(readbuf[readi].readOrigBuf.wbuf(), &bam_batches_[pt.tid_][0] + i, block_size);
			i += block_size;

			read1 = !read1;
			readi = (pp_.align_paired_reads
					 && pt.bufb_[readi].readOrigBuf.length() == 0) ? readi : readi + 1;
		}
	}

	return make_pair(done, readi);
}

void BAMPatternSource::store_orphan_mate(const uint8_t* r, size_t read_len) {
	uint8_t flag;
	memcpy(&flag, r + offset[BAMField::flag], sizeof(flag));

	std::vector<orphan_mate_t>&
		orphan_mates = (flag & 0x40) != 0 ? orphan_mate1s : orphan_mate2s;

	size_t i;
	for (i = 0; i < orphan_mates.size() && !orphan_mates[i].empty(); ++i) ;

	if (i == orphan_mates.size()) {
		orphan_mates.resize(orphan_mates.size() * 2);
	}

	orphan_mate_t& mate = orphan_mates[i];

	if (mate.data == NULL || mate.cap < read_len) {
		mate.data = new uint8_t[read_len];
		mate.cap = read_len;
	}

	if (mate.cap < read_len) {
		mate.data = new uint8_t[read_len];
		mate.cap = read_len;
	}

	mate.size = read_len;

	memcpy(mate.data, r, read_len);
}

int BAMPatternSource::compare_read_names(const void* m1, const void* m2) {
	const orphan_mate_t* mate1 = (const orphan_mate_t*)m1;
	const orphan_mate_t* mate2 = (const orphan_mate_t*)m2;

	const char* r1 = (const char *)(mate1->data + offset[BAMField::read_name]);
	const char* r2 = (const char *)(mate2->data + offset[BAMField::read_name]);

	return strcmp(r1, r2);
}

bool BAMPatternSource::compare_read_names2(const orphan_mate_t& m1, const orphan_mate_t& m2) {
	if (m1.empty()) {
		return false;
	}

	if (m2.empty()) {
		return true;
	}

	const char* r1 = (const char *)(m1.data + offset[BAMField::read_name]);
	const char* r2 = (const char *)(m2.data + offset[BAMField::read_name]);

	return strcmp(r1, r2);
}

void BAMPatternSource::get_orphaned_pairs(EList<Read>& buf_a, EList<Read>& buf_b, const size_t max_buf, unsigned& readi) {
	std::sort(orphan_mate1s.begin(), orphan_mate1s.end(), &BAMPatternSource::compare_read_names2);
	std::sort(orphan_mate2s.begin(), orphan_mate2s.end(), &BAMPatternSource::compare_read_names2);

	size_t lim1, lim2;
	for (lim1 = 0; lim1 < orphan_mate1s.size() && !orphan_mate1s[lim1].empty(); lim1++) ;
	for (lim2 = 0; lim2 < orphan_mate2s.size() && !orphan_mate2s[lim2].empty(); lim2++) ;

	if (lim2 == 0) {
		return;
	}

	for (size_t i = 0; i < lim1 && readi < max_buf; ++i) {
		orphan_mate_t* mate1 = &orphan_mate1s[i];
		orphan_mate_t* mate2 = (orphan_mate_t*)bsearch(mate1, &orphan_mate2s[0], lim2,
		                            sizeof(orphan_mate2s[0]), compare_read_names);

		if (mate2 != NULL) {
			Read& ra = buf_a[readi];
			Read& rb = buf_b[readi];

			ra.readOrigBuf.resize(mate1->size);
			rb.readOrigBuf.resize(mate2->size);

			memcpy(ra.readOrigBuf.wbuf(), mate1->data, mate1->size);
			memcpy(rb.readOrigBuf.wbuf(), mate2->data, mate2->size);

			mate1->reset();
			mate2->reset();

			readi++;
		}
	}
}

int BAMPatternSource::decompress_bgzf_block(uint8_t *dst, size_t dst_len, uint8_t *src, size_t src_len) {
	z_stream strm;

	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;

	strm.avail_in = src_len;
	strm.next_in = src;
	strm.avail_out = dst_len;
	strm.next_out = dst;

	int ret  = inflateInit2(&strm, -8);
	if (ret != Z_OK) {
		return ret;
	}

	ret = inflate(&strm, Z_FINISH);
	if (ret != Z_STREAM_END) {
		return ret;
	}

	return inflateEnd(&strm);
}

bool BAMPatternSource::parse(Read& ra, Read& rb, TReadId rdid) const {
	uint8_t l_read_name;
	int32_t l_seq;
	uint16_t n_cigar_op;
	const char* buf = ra.readOrigBuf.buf();
	int block_size = ra.readOrigBuf.length();

	memcpy(&l_read_name, buf + offset[BAMField::l_read_name], sizeof(l_read_name));
	memcpy(&n_cigar_op, buf + offset[BAMField::n_cigar_op], sizeof(n_cigar_op));
	memcpy(&l_seq, buf + offset[BAMField::l_seq], sizeof(l_seq));

	int off = offset[BAMField::read_name];
	ra.name.install(buf + off, l_read_name-1);
	off += (l_read_name + sizeof(uint32_t) * n_cigar_op);
	const char* seq = buf + off;
	off += (l_seq+1)/2;
	const char* qual = buf + off;
	for (int i = 0; i < l_seq; i++) {
		if (i < pp_.trim5) {
			ra.trimmed5 += 1;
		} else {
			ra.qual.append(qual[i] + 33);
			int base = "=ACMGRSVTWYHKDBN"[static_cast<uint8_t>(seq[i/2]) >> 4*(1-(i%2)) & 0xf];
			ra.patFw.append(asc2dna[base]);
		}
	}
	ra.trimmed3 = (int)(ra.patFw.trimEnd(pp_.trim3));
	ra.qual.trimEnd(ra.trimmed3);

	if (pp_.preserve_tags) {
		off += l_seq;
		ra.preservedOptFlags.install(buf + off, block_size - off);
	}

	ra.parsed = true;
	if (!rb.parsed && rb.readOrigBuf.length() != 0) {
		return parse(rb, ra, rdid);
	}

	return true;
}

/**
 * Light-parse a batch of tabbed-format reads into given buffer.
 */
pair<bool, int> TabbedPatternSource::nextBatchFromFile(
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
		while(c >= 0 && (c == '\n' || c == '\r') && readi < pt.max_buf_ - 1) {
			c = getc_wrapper();
		}
	}
	return make_pair(c < 0, readi);
}

/**
 * Finalize tabbed parsing outside critical section.
 */
bool TabbedPatternSource::parse(Read& ra, Read& rb, TReadId rdid) const {
	// Light parser (nextBatchFromFile) puts unparsed data
	// into Read& r, even when the read is paired.
	assert(ra.empty());
	assert(rb.empty());
	assert(!ra.readOrigBuf.empty()); // raw data for read/pair is here
	assert(rb.readOrigBuf.empty());
	int c = '\t';
	size_t cur = 0;
	const size_t buflen = ra.readOrigBuf.length();
	
	// Loop over the two ends
	for(int endi = 0; endi < 2 && c == '\t'; endi++) {
		Read& r = ((endi == 0) ? ra : rb);
		assert(r.name.empty());
		// Parse name if (a) this is the first end, or
		// (b) this is tab6
		if(endi < 1 || secondName_) {
			// Parse read name
			c = ra.readOrigBuf[cur++];
			while(c != '\t' && cur < buflen) {
				r.name.append(c);
				c = ra.readOrigBuf[cur++];
			}
			assert_eq('\t', c);
			if(cur >= buflen) {
				return false; // record ended prematurely
			}
		} else if(endi > 0) {
			// if this is the second end and we're parsing
			// tab5, copy name from first end
			rb.name = ra.name;
		}

		// Parse sequence
		assert(r.patFw.empty());
		c = ra.readOrigBuf[cur++];
		int nchar = 0;
		while(c != '\t' && cur < buflen) {
			if(isalpha(c)) {
				assert_in(toupper(c), "ACGTN");
				if(nchar++ >= pp_.trim5) {
					assert_neq(0, asc2dnacat[c]);
					r.patFw.append(asc2dna[c]); // ascii to int
				}
			}
			c = ra.readOrigBuf[cur++];
		}
		assert_eq('\t', c);
		if(cur >= buflen) {
			return false; // record ended prematurely
		}
		// record amt trimmed from 5' end due to --trim5
		r.trimmed5 = (int)(nchar - r.patFw.length());
		// record amt trimmed from 3' end due to --trim3
		r.trimmed3 = (int)(r.patFw.trimEnd(pp_.trim3));
		
		// Parse qualities
		assert(r.qual.empty());
		c = ra.readOrigBuf[cur++];
		int nqual = 0;
		if (pp_.intQuals) {
			int cur_int = 0;
			while(c != '\t' && c != '\n' && c != '\r' && cur < buflen) {
				cur_int *= 10;
				cur_int += (int)(c - '0');
				c = ra.readOrigBuf[cur++];
				if(c == ' ' || c == '\t' || c == '\n' || c == '\r') {
					char cadd = intToPhred33(cur_int, pp_.solexa64);
					cur_int = 0;
					assert_geq(cadd, 33);
					if(++nqual > pp_.trim5) {
						r.qual.append(cadd);
					}
				}
			}
		} else {
			while(c != '\t' && c != '\n' && c != '\r') {
				if(c == ' ') {
					wrongQualityFormat(r.name);
					return false;
				}
				char cadd = charToPhred33(c, pp_.solexa64, pp_.phred64);
				if(++nqual > pp_.trim5) {
					r.qual.append(cadd);
				}
				if(cur >= buflen) break;
				c = ra.readOrigBuf[cur++];
			}
		}
		if(nchar > nqual) {
			tooFewQualities(r.name);
			return false;
		} else if(nqual > nchar) {
			tooManyQualities(r.name);
			return false;
		}
		r.qual.trimEnd(pp_.trim3);
		assert(c == '\t' || c == '\n' || c == '\r' || cur >= buflen);
		assert_eq(r.patFw.length(), r.qual.length());
	}
	return true;
}

/**
 * Light-parse a batch of raw-format reads into given buffer.
 */
pair<bool, int> RawPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a,
    unsigned readi)
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
	// incase a valid character is consumed between batches
	if (c >= 0 && c != '\n' && c != '\r') {
		ungetc_wrapper(c);
	}
	return make_pair(c < 0, readi);
}

/**
 * Finalize raw parsing outside critical section.
 */
bool RawPatternSource::parse(Read& r, Read& rb, TReadId rdid) const {
	assert(r.empty());
	assert(!r.readOrigBuf.empty()); // raw data for read/pair is here
	int c = '\n';
	size_t cur = 0;
	const size_t buflen = r.readOrigBuf.length();

	// Parse sequence
	assert(r.patFw.empty());
	int nchar = 0;
	while(cur < buflen) {
		c = r.readOrigBuf[cur++];
		assert(c != '\r' && c != '\n');
		if(isalpha(c)) {
			assert_in(toupper(c), "ACGTN");
			if(nchar++ >= pp_.trim5) {
				assert_neq(0, asc2dnacat[c]);
				r.patFw.append(asc2dna[c]); // ascii to int
			}
		}
	}
	assert_eq(cur, buflen);
	// record amt trimmed from 5' end due to --trim5
	r.trimmed5 = (int)(nchar - r.patFw.length());
	// record amt trimmed from 3' end due to --trim3
	r.trimmed3 = (int)(r.patFw.trimEnd(pp_.trim3));
	
	// Give the name field a dummy value
	char cbuf[20];
	itoa10<TReadId>(rdid, cbuf);
	r.name.install(cbuf);
	
	// Give the base qualities dummy values
	assert(r.qual.empty());
	const size_t len = r.patFw.length();
	for(size_t i = 0; i < len; i++) {
		r.qual.append('I');
	}
	r.parsed = true;
	if(!rb.parsed && !rb.readOrigBuf.empty()) {
		return parse(rb, r, rdid);
	}
	return true;
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

#ifdef USE_SRA

std::pair<bool, int> SRAPatternSource::nextBatchImpl(
	PerThreadReadBuf& pt,
	bool batch_a)
{
	EList<Read>& readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	size_t readi = 0;
	bool done = false;

	for(; readi < pt.max_buf_; readi++) {
		if(!sra_its_[pt.tid_]->nextRead() || !sra_its_[pt.tid_]->nextFragment()) {
			done = true;
			break;
		}
		const ngs::StringRef rname = sra_its_[pt.tid_]->getReadId();
		const ngs::StringRef ra_seq = sra_its_[pt.tid_]->getFragmentBases();
		const ngs::StringRef ra_qual = sra_its_[pt.tid_]->getFragmentQualities();
		readbuf[readi].readOrigBuf.install(rname.data(), rname.size());
		readbuf[readi].readOrigBuf.append('\t');
		readbuf[readi].readOrigBuf.append(ra_seq.data(), ra_seq.size());
		readbuf[readi].readOrigBuf.append('\t');
		readbuf[readi].readOrigBuf.append(ra_qual.data(), ra_qual.size());
		if(sra_its_[pt.tid_]->nextFragment()) {
			const ngs::StringRef rb_seq = sra_its_[pt.tid_]->getFragmentBases();
			const ngs::StringRef rb_qual = sra_its_[pt.tid_]->getFragmentQualities();
			readbuf[readi].readOrigBuf.append('\t');
			readbuf[readi].readOrigBuf.append(rb_seq.data(), rb_seq.size());
			readbuf[readi].readOrigBuf.append('\t');
			readbuf[readi].readOrigBuf.append(rb_qual.data(), rb_qual.size());
		}
		readbuf[readi].readOrigBuf.append('\n');
	}

	pt.setReadId(readCnt_);

	{
		ThreadSafe ts(mutex);
		readCnt_ += readi;
	}

	return make_pair(done, readi);
}

/**
 * TODO: need to think about whether this can be done in a sensible way
 */
std::pair<bool, int> SRAPatternSource::nextBatch(
	PerThreadReadBuf& pt,
	bool batch_a,
	bool lock)
{
	if(lock) {
		// synchronization at this level because both reading and manipulation of
		// current file pointer have to be protected
		return nextBatchImpl(pt, batch_a);
	} else {
		return nextBatchImpl(pt, batch_a);
	}
}

/**
 * Finalize tabbed parsing outside critical section.
 */
bool SRAPatternSource::parse(Read& ra, Read& rb, TReadId rdid) const {
	// Light parser (nextBatchFromFile) puts unparsed data
	// into Read& r, even when the read is paired.
	assert(ra.empty());
	assert(rb.empty());
	assert(!ra.readOrigBuf.empty()); // raw data for read/pair is here
	assert(rb.readOrigBuf.empty());
	int c = '\t';
	size_t cur = 0;
	const size_t buflen = ra.readOrigBuf.length();

	// Loop over the two ends
	for(int endi = 0; endi < 2 && c == '\t'; endi++) {
		Read& r = ((endi == 0) ? ra : rb);
		assert(r.name.empty());
		// Parse name if (a) this is the first end, or
		// (b) this is tab6
		if(endi < 1) {
			// Parse read name
			c = ra.readOrigBuf[cur++];
			while(c != '\t' && cur < buflen) {
				r.name.append(c);
				c = ra.readOrigBuf[cur++];
			}
			assert_eq('\t', c);
			if(cur >= buflen) {
				return false; // record ended prematurely
			}
		} else if(endi > 0) {
			// if this is the second end and we're parsing
			// tab5, copy name from first end
			rb.name = ra.name;
		}

		// Parse sequence
		assert(r.patFw.empty());
		c = ra.readOrigBuf[cur++];
		int nchar = 0;
		while(c != '\t' && cur < buflen) {
			if(isalpha(c)) {
				assert_in(toupper(c), "ACGTN");
				if(nchar++ >= pp_.trim5) {
					assert_neq(0, asc2dnacat[c]);
					r.patFw.append(asc2dna[c]); // ascii to int
				}
			}
			c = ra.readOrigBuf[cur++];
		}
		assert_eq('\t', c);
		if(cur >= buflen) {
			return false; // record ended prematurely
		}
		// record amt trimmed from 5' end due to --trim5
		r.trimmed5 = (int)(nchar - r.patFw.length());
		// record amt trimmed from 3' end due to --trim3
		r.trimmed3 = (int)(r.patFw.trimEnd(pp_.trim3));

		// Parse qualities
		assert(r.qual.empty());
		c = ra.readOrigBuf[cur++];
		int nqual = 0;
		if (pp_.intQuals) {
			int cur_int = 0;
			while(c != '\t' && c != '\n' && c != '\r' && cur < buflen) {
				cur_int *= 10;
				cur_int += (int)(c - '0');
				c = ra.readOrigBuf[cur++];
				if(c == ' ' || c == '\t' || c == '\n' || c == '\r') {
					char cadd = intToPhred33(cur_int, pp_.solexa64);
					cur_int = 0;
					assert_geq(cadd, 33);
					if(++nqual > pp_.trim5) {
						r.qual.append(cadd);
					}
				}
			}
		} else {
			while(c != '\t' && c != '\n' && c != '\r') {
				if(c == ' ') {
					wrongQualityFormat(r.name);
					return false;
				}
				char cadd = charToPhred33(c, pp_.solexa64, pp_.phred64);
				if(++nqual > pp_.trim5) {
					r.qual.append(cadd);
				}
				if(cur >= buflen) break;
				c = ra.readOrigBuf[cur++];
			}
		}
		if(nchar > nqual) {
			tooFewQualities(r.name);
			return false;
		} else if(nqual > nchar) {
			tooManyQualities(r.name);
			return false;
		}
		r.qual.trimEnd(pp_.trim3);
		assert(c == '\t' || c == '\n' || c == '\r' || cur >= buflen);
		assert_eq(r.patFw.length(), r.qual.length());
	}
	return true;
}

void SRAPatternSource::open() {
	const string& sra_acc = sra_accs_[sra_acc_cur_];
	string version = "bowtie2";
	ncbi::NGS::setAppVersionString(version);
	assert(!sra_acc.empty());
	try {
		// open requested accession using SRA implementation of the API
		ngs::ReadCollection sra_run = ncbi::NGS::openReadCollection(sra_acc);

		// compute window to iterate through
		size_t MAX_ROW = sra_run.getReadCount();
		pp_.upto -= pp_.skip;

		if (pp_.upto <= MAX_ROW) {
			MAX_ROW = pp_.upto;
		}
		if(MAX_ROW < 0) {
			return;
		}

		size_t window_size = MAX_ROW / sra_its_.size();
		size_t remainder = MAX_ROW % sra_its_.size();
		size_t i = 0, start = 1;

		if (pp_.skip > 0) {
			start = pp_.skip + 1;
			readCnt_ = pp_.skip;
		}

		while (i < sra_its_.size()) {
			sra_its_[i] = new ngs::ReadIterator(sra_run.getReadRange(start, window_size, ngs::Read::all));
			assert(sra_its_[i] != NULL);

			i++;
			start += window_size;
			if (i == sra_its_.size() - 1) {
				window_size += remainder;
			}
		}

	} catch(...) {
		cerr << "Warning: Could not access \"" << sra_acc << "\" for reading; skipping..." << endl;
	}
}

#endif
