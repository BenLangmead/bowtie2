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

#ifndef PAT_H_
#define PAT_H_

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>
#include <cstring>
#include <ctype.h>
#include <fstream>
#include "alphabet.h"
#include "assert_helpers.h"
#include "tokenize.h"
#include "random_source.h"
#include "spinlock.h"
#include "threading.h"
#include "filebuf.h"
#include "qual.h"
#include "search_globals.h"
#include "sstring.h"
#include "ds.h"
#include "read.h"
#include "util.h"

/**
 * Classes and routines for reading reads from various input sources.
 */

using namespace std;

/**
 * Calculate a per-read random seed based on a combination of
 * the read data (incl. sequence, name, quals) and the global
 * seed in '_randSeed'.
 */
static inline uint32_t genRandSeed(const BTDnaString& qry,
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
 * Parameters affecting how reads and read in.
 */
struct PatternParams {
	PatternParams(
		int format_,
		bool fileParallel_,
		uint32_t seed_,
		bool useSpinlock_,
		bool solexa64_,
		bool phred64_,
		bool intQuals_,
		bool fuzzy_,
		int sampleLen_,
		int sampleFreq_,
		uint32_t skip_) :
		format(format_),
		fileParallel(fileParallel_),
		seed(seed_),
		useSpinlock(useSpinlock_),
		solexa64(solexa64_),
		phred64(phred64_),
		intQuals(intQuals_),
		fuzzy(fuzzy_),
		sampleLen(sampleLen_),
		sampleFreq(sampleFreq_),
		skip(skip_) { }

	int format;           // file format
	bool fileParallel;    // true -> wrap files with separate PairedPatternSources
	uint32_t seed;        // pseudo-random seed
	bool useSpinlock;     // use spin locks instead of pthreads
	bool solexa64;        // true -> qualities are on solexa64 scale
	bool phred64;         // true -> qualities are on phred64 scale
	bool intQuals;        // true -> qualities are space-separated numbers
	bool fuzzy;           // true -> try to parse fuzzy fastq
	int sampleLen;        // length of sampled reads for FastaContinuous...
	int sampleFreq;       // frequency of sampled reads for FastaContinuous...
	uint32_t skip;        // skip the first 'skip' patterns
};

/**
 * Encapsulates a synchronized source of patterns; usually a file.
 * Optionally reverses reads and quality strings before returning them,
 * though that is usually more efficiently done by the concrete
 * subclass.  Concrete subclasses should delimit critical sections with
 * calls to lock() and unlock().
 */
class PatternSource {

public:

	PatternSource(const PatternParams& p) :
		seed_(p.seed),
		readCnt_(0),
		numWrappers_(0),
		doLocking_(true),
		useSpinlock_(p.useSpinlock),
		lock_()
	{
		MUTEX_INIT(lock_);
	}

	virtual ~PatternSource() { }

	/**
	 * Call this whenever this PatternSource is wrapped by a new
	 * WrappedPatternSourcePerThread.  This helps us keep track of
	 * whether locks will be contended.
	 */
	void addWrapper() {
		lock();
		numWrappers_++;
		unlock();
	}
	
	/**
	 * The main member function for dispensing patterns.
	 *
	 * Returns true iff a pair was parsed succesfully.
	 */
	virtual bool nextReadPair(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired,
		bool fixName);

	/**
	 * The main member function for dispensing patterns.
	 */
	virtual bool nextRead(
		Read& r,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done);

	/**
	 * Implementation to be provided by concrete subclasses.  An
	 * implementation for this member is only relevant for formats that
	 * can read in a pair of reads in a single transaction with a
	 * single input source.  If paired-end input is given as a pair of
	 * parallel files, this member should throw an error and exit.
	 */
	virtual bool nextReadPairImpl(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired) = 0;

	/**
	 * Implementation to be provided by concrete subclasses.  An
	 * implementation for this member is only relevant for formats
	 * where individual input sources look like single-end-read
	 * sources, e.g., formats where paired-end reads are specified in
	 * parallel read files.
	 */
	virtual bool nextReadImpl(
		Read& r,
		TReadId& rdid, 
		TReadId& endid, 
		bool& success,
		bool& done) = 0;

	/// Reset state to start over again with the first read
	virtual void reset() { readCnt_ = 0; }

	/**
	 * Concrete subclasses call lock() to enter a critical region.
	 * What constitutes a critical region depends on the subclass.
	 */
	void lock() {
		if(!doLocking_) return; // no contention
#ifdef USE_SPINLOCK
		if(useSpinlock_) {
			// User can ask to use the normal pthreads lock even if
			// spinlocks are compiled in.
			spinlock_.Enter();
		} else {
#endif
			MUTEX_LOCK(lock_);
#ifdef USE_SPINLOCK
		}
#endif
	}

	/**
	 * Concrete subclasses call unlock() to exit a critical region
	 * What constitutes a critical region depends on the subclass.
	 */
	void unlock() {
		if(!doLocking_) return; // no contention
#ifdef USE_SPINLOCK
		if(useSpinlock_) {
			// User can ask to use the normal pthreads lock even if
			// spinlocks are compiled in.
			spinlock_.Leave();
		} else {
#endif
			MUTEX_UNLOCK(lock_);
#ifdef USE_SPINLOCK
		}
#endif
	}

	/**
	 * Return a new dynamically allocated PatternSource for the given
	 * format, using the given list of strings as the filenames to read
	 * from or as the sequences themselves (i.e. if -c was used).
	 */
	static PatternSource* patsrcFromStrings(
		const PatternParams& p,
		const EList<string>& qs);

	/**
	 * Return the number of reads attempted.
	 */
	TReadId readCnt() const { return readCnt_ - 1; }

protected:

	uint32_t seed_;

	/// The number of reads read by this PatternSource
	TReadId readCnt_;

	int numWrappers_;      /// # threads that own a wrapper for this PatternSource
	bool doLocking_;       /// override whether to lock (true = don't override)
	/// User can ask to use the normal pthreads-style lock even if
	/// spinlocks is enabled and compiled in.  This is sometimes better
	/// if we expect bad I/O latency on some reads.
	bool useSpinlock_;
#ifdef USE_SPINLOCK
	SpinLock spinlock_;
#endif
	MUTEX_T lock_; /// mutex for locking critical regions
};

/**
 * Abstract parent class for synhconized sources of paired-end reads
 * (and possibly also single-end reads).
 */
class PairedPatternSource {
public:
	PairedPatternSource(const PatternParams& p) : seed_(p.seed) {
		MUTEX_INIT(lock_);
	}
	virtual ~PairedPatternSource() { }

	virtual void addWrapper() = 0;
	virtual void reset() = 0;
	
	virtual bool nextReadPair(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired,
		bool fixName) = 0;
	
	virtual pair<TReadId, TReadId> readCnt() const = 0;

	/**
	 * Lock this PairedPatternSource, usually because one of its shared
	 * fields is being updated.
	 */
	void lock() {
#ifdef USE_SPINLOCK
		spinlock_.Enter();
#else
		MUTEX_LOCK(lock_);
#endif
	}

	/**
	 * Unlock this PairedPatternSource.
	 */
	void unlock() {
#ifdef USE_SPINLOCK
		spinlock_.Leave();
#else
		MUTEX_UNLOCK(lock_);
#endif
	}

	/**
	 * Given the values for all of the various arguments used to specify
	 * the read and quality input, create a list of pattern sources to
	 * dispense them.
	 */
	static PairedPatternSource* setupPatternSources(
		const EList<string>& si,    // singles, from argv
		const EList<string>& m1,    // mate1's, from -1 arg
		const EList<string>& m2,    // mate2's, from -2 arg
		const EList<string>& m12,   // both mates on each line, from --12 arg
		const EList<string>& q,     // qualities associated with singles
		const EList<string>& q1,    // qualities associated with m1
		const EList<string>& q2,    // qualities associated with m2
		const PatternParams& p,     // read-in params
		bool verbose);              // be talkative?

protected:

#ifdef USE_SPINLOCK
	SpinLock spinlock_;
#endif
	MUTEX_T lock_; /// mutex for locking critical regions
	uint32_t seed_;
};

/**
 * Encapsulates a synchronized source of both paired-end reads and
 * unpaired reads, where the paired-end must come from parallel files.
 */
class PairedSoloPatternSource : public PairedPatternSource {

public:

	PairedSoloPatternSource(
		const EList<PatternSource*>* src,
		const PatternParams& p) :
		PairedPatternSource(p),
		cur_(0),
		src_(src)
	{
		assert(src_ != NULL);
		for(size_t i = 0; i < src_->size(); i++) {
			assert((*src_)[i] != NULL);
		}
	}

	virtual ~PairedSoloPatternSource() { delete src_; }

	/**
	 * Call this whenever this PairedPatternSource is wrapped by a new
	 * WrappedPatternSourcePerThread.  This helps us keep track of
	 * whether locks within PatternSources will be contended.
	 */
	virtual void addWrapper() {
		for(size_t i = 0; i < src_->size(); i++) {
			(*src_)[i]->addWrapper();
		}
	}

	/**
	 * Reset this object and all the PatternSources under it so that
	 * the next call to nextReadPair gets the very first read pair.
	 */
	virtual void reset() {
		for(size_t i = 0; i < src_->size(); i++) {
			(*src_)[i]->reset();
		}
		cur_ = 0;
	}

	/**
	 * The main member function for dispensing pairs of reads or
	 * singleton reads.  Returns true iff ra and rb contain a new
	 * pair; returns false if ra contains a new unpaired read.
	 */
	virtual bool nextReadPair(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired,
		bool fixName);

	/**
	 * Return the number of reads attempted.
	 */
	virtual pair<TReadId, TReadId> readCnt() const {
		uint64_t ret = 0llu;
		for(size_t i = 0; i < src_->size(); i++) ret += (*src_)[i]->readCnt();
		return make_pair(ret, 0llu);
	}

protected:

	volatile uint32_t cur_; // current element in parallel srca_, srcb_ vectors
	const EList<PatternSource*>* src_; /// PatternSources for paired-end reads
};

/**
 * Encapsulates a synchronized source of both paired-end reads and
 * unpaired reads, where the paired-end must come from parallel files.
 */
class PairedDualPatternSource : public PairedPatternSource {

public:

	PairedDualPatternSource(
		const EList<PatternSource*>* srca,
		const EList<PatternSource*>* srcb,
		const PatternParams& p) :
		PairedPatternSource(p), cur_(0), srca_(srca), srcb_(srcb)
	{
		assert(srca_ != NULL);
		assert(srcb_ != NULL);
		// srca_ and srcb_ must be parallel
		assert_eq(srca_->size(), srcb_->size());
		for(size_t i = 0; i < srca_->size(); i++) {
			// Can't have NULL first-mate sources.  Second-mate sources
			// can be NULL, in the case when the corresponding first-
			// mate source is unpaired.
			assert((*srca_)[i] != NULL);
			for(size_t j = 0; j < srcb_->size(); j++) {
				assert_neq((*srca_)[i], (*srcb_)[j]);
			}
		}
	}

	virtual ~PairedDualPatternSource() {
		delete srca_;
		delete srcb_;
	}

	/**
	 * Call this whenever this PairedPatternSource is wrapped by a new
	 * WrappedPatternSourcePerThread.  This helps us keep track of
	 * whether locks within PatternSources will be contended.
	 */
	virtual void addWrapper() {
		for(size_t i = 0; i < srca_->size(); i++) {
			(*srca_)[i]->addWrapper();
			if((*srcb_)[i] != NULL) {
				(*srcb_)[i]->addWrapper();
			}
		}
	}

	/**
	 * Reset this object and all the PatternSources under it so that
	 * the next call to nextReadPair gets the very first read pair.
	 */
	virtual void reset() {
		for(size_t i = 0; i < srca_->size(); i++) {
			(*srca_)[i]->reset();
			if((*srcb_)[i] != NULL) {
				(*srcb_)[i]->reset();
			}
		}
		cur_ = 0;
	}

	/**
	 * The main member function for dispensing pairs of reads or
	 * singleton reads.  Returns true iff ra and rb contain a new
	 * pair; returns false if ra contains a new unpaired read.
	 */
	virtual bool nextReadPair(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired,
		bool fixName);
	
	/**
	 * Return the number of reads attempted.
	 */
	virtual pair<TReadId, TReadId> readCnt() const;

protected:

	volatile uint32_t cur_; // current element in parallel srca_, srcb_ vectors
	const EList<PatternSource*>* srca_; /// PatternSources for 1st mates and/or unpaired reads
	const EList<PatternSource*>* srcb_; /// PatternSources for 2nd mates
};

/**
 * Encapsulates a single thread's interaction with the PatternSource.
 * Most notably, this class holds the buffers into which the
 * PatterSource will write sequences.  This class is *not* threadsafe
 * - it doesn't need to be since there's one per thread.  PatternSource
 * is thread-safe.
 */
class PatternSourcePerThread {

public:

	PatternSourcePerThread() :
		buf1_(), buf2_(), rdid_(0xffffffff), endid_(0xffffffff) { }

	virtual ~PatternSourcePerThread() { }

	/**
	 * Read the next read pair.
	 */
	virtual bool nextReadPair(
		bool& success,
		bool& done,
		bool& paired,
		bool fixName)
	{
		return success;
	}

	Read& bufa()             { return buf1_;    }	
	Read& bufb()             { return buf2_;    }
	const Read& bufa() const { return buf1_;    }
	const Read& bufb() const { return buf2_;    }

	TReadId       rdid()  const { return rdid_;  }
	TReadId       endid() const { return endid_; }
	virtual void  reset()       { rdid_ = endid_ = 0xffffffff;  }
	
	/**
	 * Return the length of mate 1 or mate 2.
	 */
	size_t length(int mate) const {
		return (mate == 1) ? buf1_.length() : buf2_.length();
	}

protected:

	Read  buf1_;    // read buffer for mate a
	Read  buf2_;    // read buffer for mate b
	TReadId rdid_;  // index of read just read
	TReadId endid_; // index of read just read
};

/**
 * Abstract parent factory for PatternSourcePerThreads.
 */
class PatternSourcePerThreadFactory {
public:
	virtual ~PatternSourcePerThreadFactory() { }
	virtual PatternSourcePerThread* create() const = 0;
	virtual EList<PatternSourcePerThread*>* create(uint32_t n) const = 0;

	/// Free memory associated with a pattern source
	virtual void destroy(PatternSourcePerThread* patsrc) const {
		assert(patsrc != NULL);
		// Free the PatternSourcePerThread
		delete patsrc;
	}

	/// Free memory associated with a pattern source list
	virtual void destroy(EList<PatternSourcePerThread*>* patsrcs) const {
		assert(patsrcs != NULL);
		// Free all of the PatternSourcePerThreads
		for(size_t i = 0; i < patsrcs->size(); i++) {
			if((*patsrcs)[i] != NULL) {
				delete (*patsrcs)[i];
				(*patsrcs)[i] = NULL;
			}
		}
		// Free the vector
		delete patsrcs;
	}
};

/**
 * A per-thread wrapper for a PairedPatternSource.
 */
class WrappedPatternSourcePerThread : public PatternSourcePerThread {
public:
	WrappedPatternSourcePerThread(PairedPatternSource& __patsrc) :
		patsrc_(__patsrc)
	{
		patsrc_.addWrapper();
	}

	/**
	 * Get the next paired or unpaired read from the wrapped
	 * PairedPatternSource.
	 */
	virtual bool nextReadPair(
		bool& success,
		bool& done,
		bool& paired,
		bool fixName);

private:

	/// Container for obtaining paired reads from PatternSources
	PairedPatternSource& patsrc_;
};

/**
 * Abstract parent factory for PatternSourcePerThreads.
 */
class WrappedPatternSourcePerThreadFactory : public PatternSourcePerThreadFactory {
public:
	WrappedPatternSourcePerThreadFactory(PairedPatternSource& patsrc) :
		patsrc_(patsrc) { }

	/**
	 * Create a new heap-allocated WrappedPatternSourcePerThreads.
	 */
	virtual PatternSourcePerThread* create() const {
		return new WrappedPatternSourcePerThread(patsrc_);
	}

	/**
	 * Create a new heap-allocated vector of heap-allocated
	 * WrappedPatternSourcePerThreads.
	 */
	virtual EList<PatternSourcePerThread*>* create(uint32_t n) const {
		EList<PatternSourcePerThread*>* v = new EList<PatternSourcePerThread*>;
		for(size_t i = 0; i < n; i++) {
			v->push_back(new WrappedPatternSourcePerThread(patsrc_));
			assert(v->back() != NULL);
		}
		return v;
	}

private:
	/// Container for obtaining paired reads from PatternSources
	PairedPatternSource& patsrc_;
};

/// Skip to the end of the current string of newline chars and return
/// the first character after the newline chars, or -1 for EOF
static inline int getOverNewline(FileBuf& in) {
	int c;
	while(isspace(c = in.get()));
	return c;
}

/// Skip to the end of the current string of newline chars such that
/// the next call to get() returns the first character after the
/// whitespace
static inline int peekOverNewline(FileBuf& in) {
	while(true) {
		int c = in.peek();
		if(c != '\r' && c != '\n') {
			return c;
		}
		in.get();
	}
}

/// Skip to the end of the current line; return the first character
/// of the next line or -1 for EOF
static inline int getToEndOfLine(FileBuf& in) {
	while(true) {
		int c = in.get(); if(c < 0) return -1;
		if(c == '\n' || c == '\r') {
			while(c == '\n' || c == '\r') {
				c = in.get(); if(c < 0) return -1;
			}
			// c now holds first character of next line
			return c;
		}
	}
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

extern void wrongQualityFormat(const BTString& read_name);
extern void tooFewQualities(const BTString& read_name);
extern void tooManyQualities(const BTString& read_name);

/**
 * Encapsulates a source of patterns which is an in-memory vector.
 */
class VectorPatternSource : public PatternSource {

public:

	VectorPatternSource(
		const EList<string>& v,
		const PatternParams& p);
	
	virtual ~VectorPatternSource() { }
	
	virtual bool nextReadImpl(
		Read& r,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done);
	
	/**
	 * This is unused, but implementation is given for completeness.
	 */
	virtual bool nextReadPairImpl(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired);
	
	virtual void reset() {
		PatternSource::reset();
		cur_ = skip_;
		paired_ = false;
	}
	
private:

	size_t cur_;
	uint32_t skip_;
	bool paired_;
	EList<BTDnaString> v_;  // forward sequences
	EList<BTString> quals_; // forward qualities
	EList<BTString> names_; // names
	EList<int> trimmed3_;   // names
	EList<int> trimmed5_;   // names
};

/**
 *
 */
class BufferedFilePatternSource : public PatternSource {
public:
	BufferedFilePatternSource(
		const EList<string>& infiles,
		const PatternParams& p) :
		PatternSource(p),
		infiles_(infiles),
		filecur_(0),
		fb_(),
		skip_(p.skip),
		first_(true)
	{
		assert_gt(infiles.size(), 0);
		errs_.resize(infiles_.size());
		errs_.fill(0, infiles_.size(), false);
		assert(!fb_.isOpen());
		open(); // open first file in the list
		filecur_++;
	}

	virtual ~BufferedFilePatternSource() {
		if(fb_.isOpen()) fb_.close();
	}

	/**
	 * Fill Read with the sequence, quality and name for the next
	 * read in the list of read files.  This function gets called by
	 * all the search threads, so we must handle synchronization.
	 */
	virtual bool nextReadImpl(
		Read& r,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done)
	{
		// We'll be manipulating our file handle/filecur_ state
		lock();
		while(true) {
			do { read(r, rdid, endid, success, done); }
			while(!success && !done);
			if(!success && filecur_ < infiles_.size()) {
				assert(done);
				open();
				resetForNextFile(); // reset state to handle a fresh file
				filecur_++;
				continue;
			}
			break;
		}
		assert(r.repOk());
		// Leaving critical region
		unlock();
		return success;
	}
	
	/**
	 *
	 */
	virtual bool nextReadPairImpl(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired)
	{
		// We'll be manipulating our file handle/filecur_ state
		lock();
		while(true) {
			do { readPair(ra, rb, rdid, endid, success, done, paired); }
			while(!success && !done);
			if(!success && filecur_ < infiles_.size()) {
				assert(done);
				open();
				resetForNextFile(); // reset state to handle a fresh file
				filecur_++;
				continue;
			}
			break;
		}
		assert(ra.repOk());
		assert(rb.repOk());
		// Leaving critical region
		unlock();
		return success;
	}
	
	/**
	 * Reset state so that we read start reading again from the
	 * beginning of the first file.  Should only be called by the
	 * master thread.
	 */
	virtual void reset() {
		PatternSource::reset();
		filecur_ = 0,
		open();
		filecur_++;
	}

protected:

	/// Read another pattern from the input file; this is overridden
	/// to deal with specific file formats
	virtual bool read(
		Read& r,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done) = 0;
	
	/// Read another pattern pair from the input file; this is
	/// overridden to deal with specific file formats
	virtual bool readPair(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired) = 0;
	
	/// Reset state to handle a fresh file
	virtual void resetForNextFile() { }
	
	void open() {
		if(fb_.isOpen()) fb_.close();
		while(filecur_ < infiles_.size()) {
			// Open read
			FILE *in;
			if(infiles_[filecur_] == "-") {
				in = stdin;
			} else if((in = fopen(infiles_[filecur_].c_str(), "rb")) == NULL) {
				if(!errs_[filecur_]) {
					cerr << "Warning: Could not open read file \"" << infiles_[filecur_].c_str() << "\" for reading; skipping..." << endl;
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
	
	EList<string> infiles_;  // filenames for read files
	EList<bool> errs_;       // whether we've already printed an error for each file
	size_t filecur_;         // index into infiles_ of next file to read
	FileBuf fb_;             // read file currently being read from
	TReadId skip_;           // number of reads to skip
	bool first_;
};

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
	bool solexa64);

/**
 * Synchronized concrete pattern source for a list of FASTA or CSFASTA
 * (if color = true) files.
 */
class FastaPatternSource : public BufferedFilePatternSource {
public:
	FastaPatternSource(const EList<string>& infiles,
	                   const PatternParams& p) :
		BufferedFilePatternSource(infiles, p),
		first_(true), solexa64_(p.solexa64), phred64_(p.phred64), intQuals_(p.intQuals)
	{ }
	virtual void reset() {
		first_ = true;
		BufferedFilePatternSource::reset();
	}
protected:
	/**
	 * Scan to the next FASTA record (starting with >) and return the first
	 * character of the record (which will always be >).
	 */
	static int skipToNextFastaRecord(FileBuf& in) {
		int c;
		while((c = in.get()) != '>') {
			if(in.eof()) return -1;
		}
		return c;
	}

	/// Called when we have to bail without having parsed a read.
	void bail(Read& r) {
		r.reset();
		fb_.resetLastN();
	}

	/// Read another pattern from a FASTA input file
	virtual bool read(
		Read& r,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done);
	
	/// Read another pair of patterns from a FASTA input file
	virtual bool readPair(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired)
	{
		// (For now, we shouldn't ever be here)
		cerr << "In FastaPatternSource.readPair()" << endl;
		throw 1;
		return false;
	}
	
	virtual void resetForNextFile() {
		first_ = true;
	}
	
private:
	bool first_;
	bool solexa64_;
	bool phred64_;
	bool intQuals_;
};


/**
 * Tokenize a line of space-separated integer quality values.
 */
static inline bool tokenizeQualLine(
	FileBuf& filebuf,
	char *buf,
	size_t buflen,
	EList<string>& toks)
{
	size_t rd = filebuf.gets(buf, buflen);
	if(rd == 0) return false;
	assert(NULL == strrchr(buf, '\n'));
	tokenize(string(buf), " ", toks);
	return true;
}

/**
 * Synchronized concrete pattern source for a list of files with tab-
 * delimited name, seq, qual fields (or, for paired-end reads,
 * basename, seq1, qual1, seq2, qual2).
 */
class TabbedPatternSource : public BufferedFilePatternSource {

public:

	TabbedPatternSource(
		const EList<string>& infiles,
		const PatternParams& p,
		bool  secondName) :
		BufferedFilePatternSource(infiles, p),
		solQuals_(p.solexa64),
		phred64Quals_(p.phred64),
		intQuals_(p.intQuals),
		secondName_(secondName) { }

protected:

	/// Read another pattern from a FASTA input file
	virtual bool read(
		Read& r,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done);

	/// Read another pair of patterns from a FASTA input file
	virtual bool readPair(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired);
	
private:

	/**
	 * Parse a name from fb_ and store in r.  Assume that the next
	 * character obtained via fb_.get() is the first character of
	 * the sequence and the string stops at the next char upto (could
	 * be tab, newline, etc.).
	 */
	int parseName(Read& r, Read* r2, char upto = '\t');

	/**
	 * Parse a single sequence from fb_ and store in r.  Assume
	 * that the next character obtained via fb_.get() is the first
	 * character of the sequence and the sequence stops at the next
	 * char upto (could be tab, newline, etc.).
	 */
	int parseSeq(Read& r, int& charsRead, int& trim5, char upto = '\t');

	/**
	 * Parse a single quality string from fb_ and store in r.
	 * Assume that the next character obtained via fb_.get() is
	 * the first character of the quality string and the string stops
	 * at the next char upto (could be tab, newline, etc.).
	 */
	int parseQuals(Read& r, int charsRead, int dstLen, int trim5,
	               char& c2, char upto = '\t', char upto2 = -1);

	bool solQuals_;
	bool phred64Quals_;
	bool intQuals_;
	EList<string> qualToks_;
	bool secondName_;
};

/**
 * Synchronized concrete pattern source for Illumina Qseq files.  In
 * Qseq files, each read appears on a separate line and the tab-
 * delimited fields are:
 *
 * 1. Machine name
 * 2. Run number
 * 3. Lane number
 * 4. Tile number
 * 5. X coordinate of spot
 * 6. Y coordinate of spot
 * 7. Index: "Index sequence or 0. For no indexing, or for a file that
 *    has not been demultiplexed yet, this field should have a value of
 *    0."
 * 8. Read number: 1 for unpaired, 1 or 2 for paired
 * 9. Sequence
 * 10. Quality
 * 11. Filter: 1 = passed, 0 = didn't
 */
class QseqPatternSource : public BufferedFilePatternSource {

public:

	QseqPatternSource(
		const EList<string>& infiles,
	    const PatternParams& p) :
		BufferedFilePatternSource(infiles, p),
		solQuals_(p.solexa64),
		phred64Quals_(p.phred64),
		intQuals_(p.intQuals) { }

protected:

#define BAIL_UNPAIRED() { \
	peekOverNewline(fb_); \
	r.reset(); \
	success = false; \
	done = true; \
	return success; \
}

	/**
	 * Parse a name from fb_ and store in r.  Assume that the next
	 * character obtained via fb_.get() is the first character of
	 * the sequence and the string stops at the next char upto (could
	 * be tab, newline, etc.).
	 */
	int parseName(
		Read& r,      // buffer for mate 1
		Read* r2,     // buffer for mate 2 (NULL if mate2 is read separately)
		bool append,     // true -> append characters, false -> skip them
		bool clearFirst, // clear the name buffer first
		bool warnEmpty,  // emit a warning if nothing was added to the name
		bool useDefault, // if nothing is read, put readCnt_ as a default value
		int upto);       // stop parsing when we first reach character 'upto'

	/**
	 * Parse a single sequence from fb_ and store in r.  Assume
	 * that the next character obtained via fb_.get() is the first
	 * character of the sequence and the sequence stops at the next
	 * char upto (could be tab, newline, etc.).
	 */
	int parseSeq(
		Read& r,      // buffer for read
		int& charsRead,
		int& trim5,
		char upto);

	/**
	 * Parse a single quality string from fb_ and store in r.
	 * Assume that the next character obtained via fb_.get() is
	 * the first character of the quality string and the string stops
	 * at the next char upto (could be tab, newline, etc.).
	 */
	int parseQuals(
		Read& r,      // buffer for read
		int charsRead,
		int dstLen,
		int trim5,
		char& c2,
		char upto,
		char upto2);

	/**
	 * Read another pattern from a Qseq input file.
	 */
	virtual bool read(
		Read& r,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done);

	/**
	 * Read a pair of patterns from 1 Qseq file.  Note: this is never used.
	 */
	virtual bool readPair(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired)
	{
		// (For now, we shouldn't ever be here)
		cerr << "In QseqPatternSource.readPair()" << endl;
		throw 1;
		return false;
	}

	bool solQuals_;
	bool phred64Quals_;
	bool intQuals_;
	EList<string> qualToks_;
};

/**
 * Synchronized concrete pattern source for a list of FASTA files where
 * reads need to be extracted from long continuous sequences.
 */
class FastaContinuousPatternSource : public BufferedFilePatternSource {
public:
	FastaContinuousPatternSource(const EList<string>& infiles, const PatternParams& p) :
		BufferedFilePatternSource(infiles, p),
		length_(p.sampleLen), freq_(p.sampleFreq),
		eat_(length_-1), beginning_(true),
		bufCur_(0), subReadCnt_(0llu)
	{
		resetForNextFile();
	}

	virtual void reset() {
		BufferedFilePatternSource::reset();
		resetForNextFile();
	}

protected:

	/// Read another pattern from a FASTA input file
	virtual bool read(
		Read& r,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done)
	{
		success = true;
		done = false;
		r.reset();
		while(true) {
			r.color = gColor;
			int c = fb_.get();
			if(c < 0) { success = false; done = true; return success; }
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
					rdid = endid = readCnt_-1;
					break;
				}
			}
		}
		return true;
	}
	
	/// Shouldn't ever be here; it's not sensible to obtain read pairs
	// from a continuous input.
	virtual bool readPair(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired)
	{
		cerr << "In FastaContinuousPatternSource.readPair()" << endl;
		throw 1;
		return false;
	}

	/**
	 * Reset state to be read for the next file.
	 */
	virtual void resetForNextFile() {
		eat_ = length_-1;
		beginning_ = true;
		bufCur_ = 0;
		nameBuf_.clear();
		subReadCnt_ = readCnt_;
	}

private:
	size_t length_;     /// length of reads to generate
	size_t freq_;       /// frequency to sample reads
	size_t eat_;        /// number of characters we need to skip before
	                    /// we have flushed all of the ambiguous or
	                    /// non-existent characters out of our read
	                    /// window
	bool beginning_;    /// skipping over the first read length?
	char buf_[1024];    /// read buffer
	BTString nameBuf_;  /// read buffer for name of fasta record being
	                    /// split into mers
	size_t bufCur_;     /// buffer cursor; points to where we should
	                    /// insert the next character
	uint64_t subReadCnt_;/// number to subtract from readCnt_ to get
	                    /// the pat id to output (so it resets to 0 for
	                    /// each new sequence)
};

/**
 * Read a FASTQ-format file.
 * See: http://maq.sourceforge.net/fastq.shtml
 */
class FastqPatternSource : public BufferedFilePatternSource {

public:

	FastqPatternSource(const EList<string>& infiles, const PatternParams& p) :
		BufferedFilePatternSource(infiles, p),
		first_(true),
		solQuals_(p.solexa64),
		phred64Quals_(p.phred64),
		intQuals_(p.intQuals),
		fuzzy_(p.fuzzy)
	{ }
	
	virtual void reset() {
		first_ = true;
		fb_.resetLastN();
		BufferedFilePatternSource::reset();
	}
	
protected:

	/**
	 * Scan to the next FASTQ record (starting with @) and return the first
	 * character of the record (which will always be @).  Since the quality
	 * line may start with @, we keep scanning until we've seen a line
	 * beginning with @ where the line two lines back began with +.
	 */
	static int skipToNextFastqRecord(FileBuf& in, bool sawPlus) {
		int line = 0;
		int plusLine = -1;
		int c = in.get();
		int firstc = c;
		while(true) {
			if(line > 20) {
				// If we couldn't find our desired '@' in the first 20
				// lines, it's time to give up
				if(firstc == '>') {
					// That firstc is '>' may be a hint that this is
					// actually a FASTA file, so return it intact
					return '>';
				}
				// Return an error
				return -1;
			}
			if(c == -1) return -1;
			if(c == '\n') {
				c = in.get();
				if(c == '@' && sawPlus && plusLine == (line-2)) {
					return '@';
				}
				else if(c == '+') {
					// Saw a '+' at the beginning of a line; remember where
					// we saw it
					sawPlus = true;
					plusLine = line;
				}
				else if(c == -1) {
					return -1;
				}
				line++;
			}
			c = in.get();
		}
	}

	/// Read another pattern from a FASTQ input file
	virtual bool read(
		Read& r,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done);
	
	/// Read another read pair from a FASTQ input file
	virtual bool readPair(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired)
	{
		// (For now, we shouldn't ever be here)
		cerr << "In FastqPatternSource.readPair()" << endl;
		throw 1;
		return false;
	}
	
	virtual void resetForNextFile() {
		first_ = true;
	}
	
private:

	/**
	 * Do things we need to do if we have to bail in the middle of a
	 * read, usually because we reached the end of the input without
	 * finishing.
	 */
	void bail(Read& r) {
		r.patFw.clear();
		fb_.resetLastN();
	}

	bool first_;
	bool solQuals_;
	bool phred64Quals_;
	bool intQuals_;
	bool fuzzy_;
	EList<string> qualToks_;
};

/**
 * Read a Raw-format file (one sequence per line).  No quality strings
 * allowed.  All qualities are assumed to be 'I' (40 on the Phred-33
 * scale).
 */
class RawPatternSource : public BufferedFilePatternSource {

public:

	RawPatternSource(const EList<string>& infiles, const PatternParams& p) :
		BufferedFilePatternSource(infiles, p), first_(true) { }

	virtual void reset() {
		first_ = true;
		BufferedFilePatternSource::reset();
	}

protected:

	/// Read another pattern from a Raw input file
	virtual bool read(
		Read& r,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done)
	{
		int c;
		success = true;
		done = false;
		r.reset();
		c = getOverNewline(this->fb_);
		if(c < 0) {
			bail(r); success = false; done = true; return success;
		}
		assert(!isspace(c));
		r.color = gColor;
		int mytrim5 = gTrim5;
		if(first_) {
			// Check that the first character is sane for a raw file
			int cc = c;
			if(gColor) {
				if(cc >= '0' && cc <= '4') cc = "ACGTN"[(int)cc - '0'];
				if(cc == '.') cc = 'N';
			}
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
					mytrim5 += 2; // trim primer and first color
				}
			}
			if(c < 0) {
				bail(r); success = false; done = true; return success;
			}
		}
		// _in now points just past the first character of a sequence
		// line, and c holds the first character
		int chs = 0;
		while(!isspace(c) && c >= 0) {
			if(gColor) {
				if(c >= '0' && c <= '4') c = "ACGTN"[(int)c - '0'];
				if(c == '.') c = 'N';
			}
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

		rdid = endid = readCnt_-1;
		return success;
	}
	
	/// Read another read pair from a FASTQ input file
	virtual bool readPair(
		Read& ra,
		Read& rb,
		TReadId& rdid,
		TReadId& endid,
		bool& success,
		bool& done,
		bool& paired)
	{
		// (For now, we shouldn't ever be here)
		cerr << "In RawPatternSource.readPair()" << endl;
		throw 1;
		return false;
	}
	
	virtual void resetForNextFile() {
		first_ = true;
	}
	
private:

	/**
	 * Do things we need to do if we have to bail in the middle of a
	 * read, usually because we reached the end of the input without
	 * finishing.
	 */
	void bail(Read& r) {
		r.patFw.clear();
		fb_.resetLastN();
	}
	
	bool first_;
};

#endif /*PAT_H_*/
