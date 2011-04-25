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
#include "hit_set.h"
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
		const char *dumpfile_,
		bool solexa64_,
		bool phred64_,
		bool intQuals_,
		bool fuzzy_,
		int numRandom_,
		int lenRandom_,
		int sampleLen_,
		int sampleFreq_,
		uint32_t skip_) :
		format(format_),
		fileParallel(fileParallel_),
		seed(seed_),
		useSpinlock(useSpinlock_),
		dumpfile(dumpfile_),
		solexa64(solexa64_),
		phred64(phred64_),
		intQuals(intQuals_),
		fuzzy(fuzzy_),
		numRandom(numRandom_),
		lenRandom(lenRandom_),
		sampleLen(sampleLen_),
		sampleFreq(sampleFreq_),
		skip(skip_) { }

	int format;           // file format
	bool fileParallel;    // true -> wrap files with separate PairedPatternSources
	uint32_t seed;        // pseudo-random seed
	bool useSpinlock;     // use spin locks instead of pthreads
	const char *dumpfile; // file name of file to dump to
	bool solexa64;        // true -> qualities are on solexa64 scale
	bool phred64;         // true -> qualities are on phred64 scale
	bool intQuals;        // true -> qualities are space-separated numbers
	bool fuzzy;           // true -> try to parse fuzzy fastq
	int numRandom;        // number of randomly-generated reads to make
	int lenRandom;        // for randomly-generated reads, length they should be 
	int sampleLen;        // length of sampled reads for FastaContinuous...
	int sampleFreq;       // frequency of sampled reads for FastaContinuous...
	uint32_t skip;        // skip the first 'skip' patterns
};

/**
 * Encapsulates a synchronized source of patterns; usually a file.
 * Handles dumping patterns to a logfile (useful for debugging).  Also
 * optionally reverses reads and quality strings before returning them,
 * though that is usually more efficiently done by the concrete
 * subclass.  Concrete subclasses should delimit critical sections with
 * calls to lock() and unlock().
 */
class PatternSource {
public:
	PatternSource(const PatternParams& p) :
		seed_(p.seed),
		readCnt_(0),
		dumpfile_(p.dumpfile),
		numWrappers_(0),
		doLocking_(true),
		useSpinlock_(p.useSpinlock),
		lock_()
	{
		// Open dumpfile, if specified
		if(dumpfile_ != NULL) {
			out_.open(dumpfile_, ios_base::out);
			if(!out_.good()) {
				cerr << "Could not open pattern dump file \"" << dumpfile_ << "\" for writing" << endl;
				throw 1;
			}
		}
		MUTEX_INIT(lock_);
	}

	virtual ~PatternSource() { }

	/**
	 * Call this whenever this PatternSource is wrapped by a new
	 * WrappedPatternSourcePerThread.  This helps us keep track of
	 * whether locks will be contended.
	 */
	void addWrapper() {
		numWrappers_++;
	}

	/**
	 * The main member function for dispensing patterns.
	 */
	virtual void nextReadPair(
		Read& ra,
		Read& rb,
		TReadId& patid,
		bool fixName)
	{
		// nextPatternImpl does the reading from the ultimate source;
		// it is implemented in concrete subclasses
		nextReadPairImpl(ra, rb, patid);
		if(!ra.empty()) {
			// Construct the reversed versions of the fw and rc seqs
			// and quals
			ra.constructRevComps();
			ra.constructReverses();
			if(!rb.empty()) {
				rb.constructRevComps();
				rb.constructReverses();
			}
			// Fill in the random-seed field using a combination of
			// information from the user-specified seed and the read
			// sequence, qualities, and name
			ra.seed = genRandSeed(ra.patFw, ra.qual, ra.name, seed_);
			if(!rb.empty()) {
				rb.seed = genRandSeed(rb.patFw, rb.qual, rb.name, seed_);
			}
			// Output it, if desired
			if(dumpfile_ != NULL) {
				dumpBuf(ra);
				if(!rb.empty()) {
					dumpBuf(rb);
				}
			}
			if(gVerbose) {
				cout << "Parsed mate 1: "; ra.dump(cout);
				cout << "Parsed mate 2: "; rb.dump(cout);
			}
		}
	}

	/**
	 * The main member function for dispensing patterns.
	 */
	virtual void nextRead(Read& r, TReadId& patid) {
		// nextPatternImpl does the reading from the ultimate source;
		// it is implemented in concrete subclasses
		nextReadImpl(r, patid);
		if(!r.empty()) {
			// Construct the reversed versions of the fw and rc seqs
			// and quals
			r.constructRevComps();
			r.constructReverses();
			// Fill in the random-seed field using a combination of
			// information from the user-specified seed and the read
			// sequence, qualities, and name
			r.seed = genRandSeed(r.patFw, r.qual, r.name, seed_);
			// Output it, if desired
			if(dumpfile_ != NULL) {
				dumpBuf(r);
			}
			if(gVerbose) {
				cout << "Parsed read: "; r.dump(cout);
			}
		}
	}

	/**
	 * Implementation to be provided by concrete subclasses.  An
	 * implementation for this member is only relevant for formats that
	 * can read in a pair of reads in a single transaction with a
	 * single input source.  If paired-end input is given as a pair of
	 * parallel files, this member should throw an error and exit.
	 */
	virtual void nextReadPairImpl(Read& ra, Read& rb, TReadId& patid) = 0;

	/**
	 * Implementation to be provided by concrete subclasses.  An
	 * implementation for this member is only relevant for formats
	 * where individual input sources look like single-end-read
	 * sources, e.g., formats where paired-end reads are specified in
	 * parallel read files.
	 */
	virtual void nextReadImpl(Read& r, TReadId& patid) = 0;

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
		const EList<string>& qs,
		const EList<string>* qualities);

	/**
	 * Return the number of reads attempted.
	 */
	TReadId readCnt() const { return readCnt_ - 1; }

protected:

	/**
	 * Dump the contents of the Read to the dump file.
	 */
	void dumpBuf(const Read& r) {
		assert(dumpfile_ != NULL);
		BTString empty("(empty)");
		dump(out_, r.patFw,
		     r.qual.empty()    ? empty : r.qual,
		     r.name.empty()    ? empty : r.name);
		dump(out_, r.patRc,
		     r.qualRev.empty() ? empty : r.qualRev,
		     r.name.empty()    ? empty : r.name);
	}

	/**
	 * Default format for dumping a read to an output stream.  Concrete
	 * subclasses might want to do something fancier.
	 */
	virtual void dump(ostream& out,
	                  const BTDnaString& seq,
	                  const BTString& qual,
	                  const BTString& name)
	{
		out << name << ": " << seq << " " << qual << endl;
	}

	uint32_t seed_;

	/// The number of reads read by this PatternSource
	TReadId readCnt_;

	const char *dumpfile_; /// dump patterns to this file before returning them
	ofstream out_;         /// output stream for dumpfile
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
		TReadId& patid,
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

	PairedSoloPatternSource(const EList<PatternSource*>* src,
	                        const PatternParams& p) :
		PairedPatternSource(p), cur_(0), src_(src)
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
		TReadId& patid,
		bool fixName)
	{
		uint32_t cur = cur_;
		while(cur < src_->size()) {
			// Patterns from srca_[cur_] are unpaired
			(*src_)[cur]->nextReadPair(ra, rb, patid, fixName);
			if(ra.patFw.empty()) {
				// If patFw is empty, that's our signal that the
				// input dried up
				lock();
				if(cur + 1 > cur_) cur_++;
				cur = cur_;
				unlock();
				continue; // on to next pair of PatternSources
			}
			ra.seed = genRandSeed(ra.patFw, ra.qual, ra.name, seed_);
			if(!rb.empty()) {
				rb.seed = genRandSeed(rb.patFw, rb.qual, rb.name, seed_);
				if(fixName) {
					ra.fixMateName(1);
					rb.fixMateName(2);
				}
			}
			ra.patid = patid;
			ra.mate  = 1;
			rb.mate  = 2;
			return true; // paired
		}
		return false;
	}

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

	PairedDualPatternSource(const EList<PatternSource*>* srca,
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
		TReadId& patid,
		bool fixName)
	{
		uint32_t cur = cur_;
		while(cur < srca_->size()) {
			if((*srcb_)[cur] == NULL) {
				// Patterns from srca_[cur_] are unpaired
				(*srca_)[cur]->nextRead(ra, patid);
				if(ra.patFw.empty()) {
					// If patFw is empty, that's our signal that the
					// input dried up
					lock();
					if(cur + 1 > cur_) cur_++;
					cur = cur_;
					unlock();
					continue; // on to next pair of PatternSources
				}
				ra.patid = patid;
				ra.mate  = 0;
				return false; // unpaired
			} else {
				// Patterns from srca_[cur_] and srcb_[cur_] are paired
				TReadId patid_a = 0;
				TReadId patid_b = 0;
				// Lock to ensure that this thread gets parallel reads
				// in the two mate files
				lock();
				(*srca_)[cur]->nextRead(ra, patid_a);
				(*srcb_)[cur]->nextRead(rb, patid_b);
				bool cont = false;
				// Did the pair obtained fail to match up?
				while(patid_a != patid_b) {
					// Is either input exhausted?  If so, bail.
					if(ra.patFw.empty() || rb.patFw.empty()) {
						ra.patFw.clear();
						if(cur + 1 > cur_) cur_++;
						cur = cur_;
						cont = true;
						break;
					}
					if(patid_a < patid_b) {
						(*srca_)[cur]->nextRead(ra, patid_a);
						if(fixName) {
							ra.fixMateName(1);
						}
					} else {
						(*srcb_)[cur]->nextRead(rb, patid_b);
						if(fixName) {
							rb.fixMateName(2);
						}
					}
				}
				unlock();
				if(cont) continue; // on to next pair of PatternSources
				if(fixName) {
					ra.fixMateName(1);
					rb.fixMateName(2);
				}
				if(ra.patFw.empty()) {
					// If patFw is empty, that's our signal that the
					// input dried up
					lock();
					if(cur + 1 > cur_) cur_++;
					cur = cur_;
					unlock();
					continue; // on to next pair of PatternSources
				}
				assert_eq(patid_a, patid_b);
				patid = patid_a;
				ra.patid = patid;
				rb.patid = patid;
				ra.mate  = 1;
				rb.mate  = 2;
				return true; // paired
			}
		}
		return false;
	}

	/**
	 * Return the number of reads attempted.
	 */
	virtual pair<TReadId, TReadId> readCnt() const {
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
		buf1_(), buf2_(), patid_(0xffffffff) { }

	virtual ~PatternSourcePerThread() { }

	/**
	 * Read the next read pair.
	 */
	virtual void nextReadPair(bool fixName) { }

	Read& bufa()        { return buf1_;         }
	Read& bufb()        { return buf2_;         }

	TReadId       patid() const { return patid_;        }
	virtual void  reset()       { patid_ = 0xffffffff;  }
	bool          empty() const { return buf1_.empty(); }
	
	size_t length(int mate) const {
		return (mate == 1) ? buf1_.length() : buf2_.length();
	}

	/**
	 * Return true iff the buffers jointly contain a paired-end read.
	 */
	bool paired() {
		bool ret = !buf2_.empty();
		assert(!ret || !empty());
		return ret;
	}

protected:
	Read  buf1_;    // read buffer for mate a
	Read  buf2_;    // read buffer for mate b
	TReadId patid_; // index of read just read
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
	virtual void nextReadPair(bool fixName) {
		PatternSourcePerThread::nextReadPair(fixName);
		ASSERT_ONLY(TReadId lastPatid = patid_);
		buf1_.reset();
		buf2_.reset();
		patsrc_.nextReadPair(buf1_, buf2_, patid_, fixName);
		assert(buf1_.empty() || patid_ != lastPatid);
	}

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

/**
 * A synchronized pattern source that simply returns random reads
 * without reading from the disk or storing lists of reads in memory.
 * Reads are generated with a RandomSource.
 */
class RandomPatternSource : public PatternSource {
public:
	RandomPatternSource(const PatternParams& p) :
		PatternSource(p), numReads_(p.numRandom), length_(p.lenRandom), seed_(p.seed)
	{
		//if(length_ > 1024) {
		//	cerr << "Read length for RandomPatternSource may not exceed 1024; got " << length_ << endl;
		//	throw 1;
		//}
		rand_.init(seed_);
	}

	/** Get the next random read and set patid */
	virtual void nextReadImpl(Read& r, TReadId& patid) {
		// Begin critical section
		lock();
		if(readCnt_ >= numReads_) {
			r.reset();
			unlock();
			return;
		}
		uint32_t ra = rand_.nextU32();
		patid = readCnt_;
		readCnt_++;
		unlock();
		fillRandomRead(r, ra, length_, patid);
	}

	/** Get the next random read and set patid */
	virtual void nextReadPairImpl(Read& ra, Read& rb, TReadId& patid) {
		// Begin critical section
		lock();
		if(readCnt_ >= numReads_) {
			ra.reset();
			rb.reset();
			unlock();
			return;
		}
		uint32_t rna = rand_.nextU32();
		uint32_t rnb = rand_.nextU32();
		patid = readCnt_;
		readCnt_++;
		unlock();
		fillRandomRead(ra, rna, length_, patid);
		fillRandomRead(rb, rnb, length_, patid);
	}

	/** */
	static void fillRandomRead(Read& r,
	                           uint32_t ra,
	                           int length,
	                           TReadId patid)
	{
		// End critical section
		for(int i = 0; i < length; i++) {
			ra = RandomSource::nextU32(ra) >> 8;
			r.patFw.append((char)(ra & 3));
			char c = 'I' - ((ra >> 2) & 31);
			r.qual.append(c);
		}
		char cbuf[20];
		itoa10<TReadId>(patid, cbuf);
		r.name.install(cbuf);
	}

	/** Reset the pattern source to the beginning */
	virtual void reset() {
		PatternSource::reset();
		// reset pseudo-random generator; next string of calls to
		// nextU32() will return same pseudo-randoms as the last
		rand_.init(seed_);
	}
private:
	TReadId      numReads_; /// number of reads to dish out
	int          length_;   /// length of reads
	uint32_t     seed_;     /// seed for pseudo-randoms
	RandomSource rand_;     /// pseudo-random generator
};

/**
 * A version of PatternSourcePerThread that dishes out random patterns
 * without any synchronization.
 */
class RandomPatternSourcePerThread : public PatternSourcePerThread {
public:
	RandomPatternSourcePerThread(TReadId numreads,
	                             int length,
	                             int numthreads,
	                             int thread) :
		PatternSourcePerThread(),
		numreads_(numreads),
		length_(length),
		numthreads_(numthreads),
		thread_(thread)
	{
		patid_ = thread_;
		//if(length_ > 1024) {
		//	cerr << "Read length for RandomPatternSourcePerThread may not exceed 1024; got " << length_ << endl;
		//	throw 1;
		//}
		rand_.init(thread_);
	}

	virtual void nextReadPair(bool fixName) {
		PatternSourcePerThread::nextReadPair(fixName);
		if(patid_ >= numreads_) {
			buf1_.reset();
			buf2_.reset();
			return;
		}
		RandomPatternSource::fillRandomRead(
			buf1_, rand_.nextU32(), length_, patid_);
		RandomPatternSource::fillRandomRead(
			buf2_, rand_.nextU32(), length_, patid_);
		patid_ += numthreads_;
	}

	virtual void reset() {
		PatternSourcePerThread::reset();
		patid_ = thread_;
		rand_.init(thread_);
	}

private:
	TReadId      numreads_;
	int          length_;
	int          numthreads_;
	int          thread_;
	RandomSource rand_;
};

/**
 * Abstract parent factory for PatternSourcePerThreads.
 */
class RandomPatternSourcePerThreadFactory : public PatternSourcePerThreadFactory {
public:
	RandomPatternSourcePerThreadFactory(
	        uint32_t numreads,
	        int length,
	        int numthreads,
	        int thread) :
	        numreads_(numreads),
	        length_(length),
	        numthreads_(numthreads),
	        thread_(thread) { }

	/**
	 * Create a new heap-allocated WrappedPatternSourcePerThreads.
	 */
	virtual PatternSourcePerThread* create() const {
		return new RandomPatternSourcePerThread(
			numreads_, length_, numthreads_, thread_);
	}

	/**
	 * Create a new heap-allocated vector of heap-allocated
	 * WrappedPatternSourcePerThreads.
	 */
	virtual EList<PatternSourcePerThread*>* create(uint32_t n) const {
		EList<PatternSourcePerThread*>* v = new EList<PatternSourcePerThread*>;
		for(size_t i = 0; i < n; i++) {
			v->push_back(new RandomPatternSourcePerThread(
				numreads_, length_, numthreads_, thread_));
			assert(v->back() != NULL);
		}
		return v;
	}

private:
	TReadId numreads_;
	int length_;
	int numthreads_;
	int thread_;
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
extern void tooManySeqChars(const BTString& read_name);

/**
 * Encapsulates a source of patterns which is an in-memory vector.
 */
class VectorPatternSource : public PatternSource {
public:
	VectorPatternSource(const EList<string>& v, const PatternParams& p) :
		PatternSource(p),
		cur_(p.skip), skip_(p.skip), paired_(false), v_(), quals_()
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
				continue;
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
	virtual ~VectorPatternSource() { }
	virtual void nextReadImpl(Read& r, TReadId& patid) {
		// Let Strings begin at the beginning of the respective bufs
		r.reset();
		lock();
		if(cur_ >= v_.size()) {
			unlock();
			// Clear all the Strings, as a signal to the caller that
			// we're out of reads
			r.reset();
			assert(r.empty());
			return;
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
		readCnt_++;
		patid = readCnt_;
		unlock();
	}
	/**
	 * This is unused, but implementation is given for completeness.
	 */
	virtual void nextReadPairImpl(Read& ra, Read& rb, TReadId& patid) {
		// Let Strings begin at the beginning of the respective bufs
		ra.reset();
		rb.reset();
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
			return;
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
		readCnt_++;
		patid = readCnt_;
		unlock();
	}
	virtual void reset() {
		PatternSource::reset();
		cur_ = skip_;
		paired_ = false;
	}
private:
	size_t cur_;
	uint32_t skip_;
	bool paired_;
	EList<BTDnaString> v_;     // forward sequences
	EList<BTString> quals_; // forward qualities
	EList<BTString> names_; // names
	EList<int> trimmed3_; // names
	EList<int> trimmed5_; // names
};

/**
 *
 */
class BufferedFilePatternSource : public PatternSource {
public:
	BufferedFilePatternSource(const EList<string>& infiles,
	                          const EList<string>* qinfiles,
	                          const PatternParams& p) :
		PatternSource(p),
		infiles_(infiles),
		filecur_(0),
		fb_(),
		qfb_(),
		skip_(p.skip),
		first_(true)
	{
		qinfiles_.clear();
		if(qinfiles != NULL) qinfiles_ = *qinfiles;
		assert_gt(infiles.size(), 0);
		errs_.resize(infiles_.size());
		errs_.fill(0, infiles_.size(), false);
		if(qinfiles_.size() > 0 &&
		   qinfiles_.size() != infiles_.size())
		{
			cerr << "Error: Different numbers of input FASTA/quality files ("
			     << infiles_.size() << "/" << qinfiles_.size() << ")" << endl;
			throw 1;
		}
		assert(!fb_.isOpen());
		assert(!qfb_.isOpen());
		open(); // open first file in the list
		filecur_++;
	}

	virtual ~BufferedFilePatternSource() {
		if(fb_.isOpen()) fb_.close();
		if(qfb_.isOpen()) {
			assert_gt(qinfiles_.size(), 0);
			qfb_.close();
		}
	}

	/**
	 * Fill Read with the sequence, quality and name for the next
	 * read in the list of read files.  This function gets called by
	 * all the search threads, so we must handle synchronization.
	 */
	virtual void nextReadImpl(Read& r, TReadId& patid) {
		// We are entering a critical region, because we're
		// manipulating our file handle and filecur_ state
		lock();
		bool notDone = true;
		do {
			read(r, patid);
			// Try again if r is empty (indicating an error) and input
			// is not yet exhausted, OR if we have more reads to skip
			// over
			notDone = r.patFw.empty() && !fb_.eof();
		} while(notDone || (!fb_.eof() && patid < skip_));
		if(patid < skip_) {
			unlock();
			r.reset();
			assert(r.patFw.empty());
			return;
		}
		if(first_ && r.patFw.empty() /* && !quiet_ */) {
			// No reads could be extracted from the first _infile
			cerr << "Warning: Could not find any reads in \"" << infiles_[0] << "\"" << endl;
		}
		first_ = false;
		while(r.patFw.empty() && filecur_ < infiles_.size()) {
			// Open next file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			do {
				read(r, patid);
			} while((r.patFw.empty() && !fb_.eof()));
			assert_geq(patid, skip_);
			if(r.patFw.empty() /*&& !quiet_ */) {
				// No reads could be extracted from this _infile
				cerr << "Warning: Could not find any reads in \"" << infiles_[filecur_] << "\"" << endl;
			}
			filecur_++;
		}
		assert(r.repOk());
		// Leaving critical region
		unlock();
		// If r.patFw is empty, then the caller knows that we are
		// finished with the reads
	}
	
	/**
	 *
	 */
	virtual void nextReadPairImpl(Read& ra, Read& rb, TReadId& patid) {
		// We are entering a critical region, because we're
		// manipulating our file handle and filecur_ state
		lock();
		bool notDone = true;
		do {
			readPair(ra, rb, patid);
			// Try again if ra is empty (indicating an error) and input
			// is not yet exhausted, OR if we have more reads to skip
			// over
			notDone = ra.patFw.empty() && !fb_.eof();
		} while(notDone || (!fb_.eof() && patid < skip_));
		if(patid < skip_) {
			unlock();
			ra.reset();
			rb.reset();
			assert(ra.patFw.empty());
			return;
		}
		if(first_ && ra.patFw.empty() /*&& !quiet_*/) {
			// No reads could be extracted from the first _infile
			cerr << "Warning: Could not find any read pairs in \"" << infiles_[0] << "\"" << endl;
		}
		first_ = false;
		while(ra.patFw.empty() && filecur_ < infiles_.size()) {
			// Open next file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			do {
				readPair(ra, rb, patid);
			} while((ra.patFw.empty() && !fb_.eof()));
			assert_geq(patid, skip_);
			if(ra.patFw.empty() /*&& !quiet_*/) {
				// No reads could be extracted from this _infile
				cerr << "Warning: Could not find any reads in \"" << infiles_[filecur_] << "\"" << endl;
			}
			filecur_++;
		}
		// Leaving critical region
		unlock();
		// If ra.patFw is empty, then the caller knows that we are
		// finished with the reads
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
	virtual void read(Read& r, TReadId& patid) = 0;
	/// Read another pattern pair from the input file; this is
	/// overridden to deal with specific file formats
	virtual void readPair(Read& ra, Read& rb, TReadId& patid) = 0;
	/// Reset state to handle a fresh file
	virtual void resetForNextFile() { }
	void open() {
		if(fb_.isOpen()) fb_.close();
		if(qfb_.isOpen()) qfb_.close();
		while(filecur_ < infiles_.size()) {
			// Open read
			FILE *in;
			if(infiles_[filecur_] == "-") {
				in = stdin;
			} else if((in = fopen(infiles_[filecur_].c_str(), "rb")) == NULL) {
				if(!errs_[filecur_]) {
					cerr << "Warning: Could not open read file \"" << infiles_[filecur_] << "\" for reading; skipping..." << endl;
					errs_[filecur_] = true;
				}
				filecur_++;
				continue;
			}
			fb_.newFile(in);
			// Open quality
			if(!qinfiles_.empty()) {
				FILE *in;
				if(qinfiles_[filecur_] == "-") {
					in = stdin;
				} else if((in = fopen(qinfiles_[filecur_].c_str(), "rb")) == NULL) {
					if(!errs_[filecur_]) {
						cerr << "Warning: Could not open quality file \"" << qinfiles_[filecur_] << "\" for reading; skipping..." << endl;
						errs_[filecur_] = true;
					}
					filecur_++;
					continue;
				}
				qfb_.newFile(in);
			}
			return;
		}
		throw 1;
	}
	EList<string> infiles_; /// filenames for read files
	EList<string> qinfiles_; /// filenames for quality files
	EList<bool> errs_; /// whether we've already printed an error for each file
	size_t filecur_;   /// index into infiles_ of next file to read
	FileBuf fb_;  /// read file currently being read from
	FileBuf qfb_; /// quality file currently being read from
	TReadId skip_;     /// number of reads to skip
	bool first_;
};

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
               bool solexa64);

/**
 * Synchronized concrete pattern source for a list of FASTA or CSFASTA
 * (if color = true) files.
 */
class FastaPatternSource : public BufferedFilePatternSource {
public:
	FastaPatternSource(const EList<string>& infiles,
	                   const EList<string>* qinfiles,
	                   const PatternParams& p) :
		BufferedFilePatternSource(infiles, qinfiles, p),
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
		qfb_.resetLastN();
	}

	/// Read another pattern from a FASTA input file
	virtual void read(Read& r, TReadId& patid) {
		int c, qc = 0;
		bool doquals = qinfiles_.size() > 0;
		assert(!doquals || qfb_.isOpen());
		assert(fb_.isOpen());
		r.reset();
		r.color = gColor;
		// Pick off the first carat
		c = fb_.get();
		if(c < 0) { bail(r); return; }
		while(c == '#' || c == ';') {
			c = fb_.peekUptoNewline();
			fb_.resetLastN();
			c = fb_.get();
		}
		assert_eq(1, fb_.lastNLen());
		if(doquals) {
			qc = qfb_.get();
			if(qc < 0) { bail(r); return; }
			while(qc == '#' || qc == ';') {
				qc = qfb_.peekUptoNewline();
				qfb_.resetLastN();
				qc = qfb_.get();
			}
			assert_eq(1, qfb_.lastNLen());
		}

		// Pick off the first carat
		if(first_) {
			if(c != '>') {
				cerr << "Error: reads file does not look like a FASTA file" << endl;
				throw 1;
			}
			if(doquals && qc != '>') {
				cerr << "Error: quality file does not look like a FASTA quality file" << endl;
				throw 1;
			}
			first_ = false;
		}
		assert_eq('>', c);
		if(doquals) assert_eq('>', qc);
		c = fb_.get(); // get next char after '>'
		if(doquals) qc = qfb_.get();

		// Read to the end of the id line, sticking everything after the '>'
		// into *name
		bool warning = false;
		while(true) {
			if(c < 0 || qc < 0) { bail(r); return; }
			if(c == '\n' || c == '\r') {
				// Break at end of line, after consuming all \r's, \n's
				while(c == '\n' || c == '\r') {
					if(doquals && c != qc) {
						cerr << "Warning: one or more mismatched read names between FASTA and quality files" << endl;
						warning = true;
					}
					if(fb_.peek() == '>') {
						// Empty sequence
						break;
					}
					c = fb_.get();
					if(doquals) qc = qfb_.get();
					if(c < 0 || qc < 0) { bail(r); return; }
				}
				break;
			}
			if(doquals && c != qc) {
				cerr << "Warning: one or more mismatched read names between FASTA and quality files" << endl;
				warning = true;
			}
			r.name.append(c);
			if(fb_.peek() == '>') {
				// Empty sequence
				break;
			}
			c = fb_.get();
			if(doquals) qc = qfb_.get();
		}
		if(c == '>') {
			// Empty sequences!
			cerr << "Warning: skipping empty FASTA read with name '" << r.name << "'" << endl;
			fb_.resetLastN();
			return;
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
			if(c < 0) { bail(r); return; }
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
		if(doquals) {
			parseQuals(r, qfb_, qc, (int)(r.patFw.length() + r.trimmed3 + r.trimmed5),
			           r.trimmed3, r.trimmed5, intQuals_, phred64_,
			           solexa64_);
		}
		// Set up a default name if one hasn't been set
		if(r.name.empty()) {
			char cbuf[20];
			itoa10<TReadId>(readCnt_, cbuf);
			r.name.install(cbuf);
		}
		assert_gt(r.name.length(), 0);
		readCnt_++;
		patid = readCnt_-1;
		r.readOrigBuf.install(fb_.lastN(), fb_.lastNLen());
		fb_.resetLastN();
		if(doquals) {
			r.qualOrigBuf.install(qfb_.lastN(), qfb_.lastNLen());
			qfb_.resetLastN();
			if(false) {
				cout << "Name: " << r.name << endl
					 << " Seq: " << r.patFw << " (" << r.patFw.length() << ")" << endl
					 << "Qual: " << r.qual  << " (" << r.qual.length() << ")" << endl
					 << "Orig seq:" << endl;
				cout << r.readOrigBuf.toZBuf();
				cout << "Orig qual:" << endl;
				cout << r.qualOrigBuf.toZBuf();
				cout << endl;
			}
		}
	}
	/// Read another pair of patterns from a FASTA input file
	virtual void readPair(Read& ra, Read& rb, TReadId& patid) {
		// (For now, we shouldn't ever be here)
		cerr << "In FastaPatternSource.readPair()" << endl;
		throw 1;
		read(ra, patid);
		readCnt_--;
		read(rb, patid);
	}
	virtual void resetForNextFile() {
		first_ = true;
	}
	virtual void dump(ostream& out,
	                  const BTDnaString& seq,
	                  const BTString& qual,
	                  const BTString& name)
	{
		out << ">" << name << endl << seq << endl;
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
static inline bool tokenizeQualLine(FileBuf& filebuf, char *buf, size_t buflen, EList<string>& toks) {
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
	TabbedPatternSource(const EList<string>& infiles,
	                    const PatternParams& p) :
		BufferedFilePatternSource(infiles, NULL, p),
		solQuals_(p.solexa64),
		phred64Quals_(p.phred64),
		intQuals_(p.intQuals) { }

protected:

	/// Read another pattern from a FASTA input file
	virtual void read(Read& r, TReadId& patid) {
		r.reset();
		r.color = gColor;
		// fb_ is about to dish out the first character of the
		// name field
		if(parseName(r, NULL, '\t') == -1) {
			peekOverNewline(fb_); // skip rest of line
			r.reset();
			return;
		}
		assert_neq('\t', fb_.peek());

		// fb_ is about to dish out the first character of the
		// sequence field
		int charsRead = 0;
		int mytrim5 = gTrim5;
		int dstLen = parseSeq(r, charsRead, mytrim5, '\t');
		assert_neq('\t', fb_.peek());
		if(dstLen <= 0) {
			peekOverNewline(fb_); // skip rest of line
			r.reset();
			return;
		}

		// fb_ is about to dish out the first character of the
		// quality-string field
		char ct = 0;
		if(parseQuals(r, charsRead, dstLen, mytrim5, ct, '\n') <= 0) {
			peekOverNewline(fb_); // skip rest of line
			r.reset();
			return;
		}
		r.trimmed3 = gTrim3;
		r.trimmed5 = mytrim5;
		assert_eq(ct, '\n');
		assert_neq('\n', fb_.peek());
		r.readOrigBuf.install(fb_.lastN(), fb_.lastNLen());
		fb_.resetLastN();
		// The last character read in parseQuals should have been a
		// '\n'

		readCnt_++;
		patid = readCnt_-1;
	}

	/// Read another pair of patterns from a FASTA input file
	virtual void readPair(Read& ra, Read& rb, TReadId& patid) {
		// fb_ is about to dish out the first character of the
		// name field
		int mytrim5_1 = gTrim5;
		if(parseName(ra, &rb, '\t') == -1) {
			peekOverNewline(fb_); // skip rest of line
			ra.reset();
			rb.reset();
			fb_.resetLastN();
			return;
		}
		assert_neq('\t', fb_.peek());

		// fb_ is about to dish out the first character of the
		// sequence field for the first mate
		int charsRead1 = 0;
		int dstLen1 = parseSeq(ra, charsRead1, mytrim5_1, '\t');
		if(dstLen1 <= -1) {
			peekOverNewline(fb_); // skip rest of line
			ra.reset();
			rb.reset();
			fb_.resetLastN();
			return;
		}
		assert_neq('\t', fb_.peek());

		// fb_ is about to dish out the first character of the
		// quality-string field
		char ct = 0;
		if(parseQuals(ra, charsRead1, dstLen1, mytrim5_1, ct, '\t', '\n') <= 0) {
			peekOverNewline(fb_); // skip rest of line
			ra.reset();
			rb.reset();
			fb_.resetLastN();
			return;
		}
		ra.trimmed3 = gTrim3;
		ra.trimmed5 = mytrim5_1;
		assert(ct == '\t' || ct == '\n');
		if(ct == '\n') {
			rb.reset();
			peekOverNewline(fb_);
			ra.readOrigBuf.install(fb_.lastN(), fb_.lastNLen());
			fb_.resetLastN();
			readCnt_++;
			patid = readCnt_-1;
			return;
		}
		assert_neq('\t', fb_.peek());

		// fb_ is about to dish out the first character of the
		// sequence field for the second mate
		int charsRead2 = 0;
		int mytrim5_2 = gTrim5;
		int dstLen2 = parseSeq(rb, charsRead2, mytrim5_2, '\t');
		if(dstLen2 <= 0) {
			peekOverNewline(fb_); // skip rest of line
			ra.reset();
			rb.reset();
			fb_.resetLastN();
			return;
		}
		assert_neq('\t', fb_.peek());

		// fb_ is about to dish out the first character of the
		// quality-string field
		if(parseQuals(rb, charsRead2, dstLen2, mytrim5_2, ct, '\n') <= 0) {
			peekOverNewline(fb_); // skip rest of line
			ra.reset();
			rb.reset();
			fb_.resetLastN();
			return;
		}
		assert_eq('\n', ct);
		if(fb_.peek() == '\n') {
			assert(false);
		}
		peekOverNewline(fb_);
		ra.readOrigBuf.install(fb_.lastN(), fb_.lastNLen());
		fb_.resetLastN();

		rb.trimmed3 = gTrim3;
		rb.trimmed5 = mytrim5_2;

		// The last character read in parseQuals should have been a
		// '\n'

		readCnt_++;
		patid = readCnt_-1;
	}

	/**
	 * Dump a FASTQ-style record for the read.
	 */
	virtual void dump(ostream& out,
	                  const BTDnaString& seq,
	                  const BTString& qual,
	                  const BTString& name)
	{
		out << "@" << name << endl << seq << endl
		    << "+" << endl << qual << endl;
	}
private:

	/**
	 * Parse a name from fb_ and store in r.  Assume that the next
	 * character obtained via fb_.get() is the first character of
	 * the sequence and the string stops at the next char upto (could
	 * be tab, newline, etc.).
	 */
	int parseName(Read& r, Read* r2, char upto = '\t') {
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
	int parseSeq(Read& r, int& charsRead, int& trim5, char upto = '\t') {
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
	int parseQuals(Read& r, int charsRead, int dstLen, int trim5,
	               char& c2, char upto = '\t', char upto2 = -1)
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
						//size_t off = qualsRead - trim5;
						//if(off >= 1024) tooManyQualities(r.name);
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
						//size_t off = qualsRead - trim5;
						//if(off >= 1024) tooManyQualities(r.name);
						c = charToPhred33(c, solQuals_, phred64Quals_);
						assert_geq(c, 33);
						r.qual.append(c);
					}
					qualsRead++;
				} else {
					break;
				}
			}
			if(qualsRead != dstLen + trim5) {
				assert(false);
			}
		}
		r.qual.resize(dstLen);
		while(c != upto && (upto2 == -1 || c != upto2)) {
			c = fb_.get();
			c2 = c;
		}
		return qualsRead;
	}

	bool solQuals_;
	bool phred64Quals_;
	bool intQuals_;
	EList<string> qualToks_;
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
		BufferedFilePatternSource(infiles, NULL, p),
		solQuals_(p.solexa64),
		phred64Quals_(p.phred64),
		intQuals_(p.intQuals) { }

protected:

#define BAIL_UNPAIRED() { \
	peekOverNewline(fb_); \
	r.reset(); \
	return; \
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
	virtual void read(Read& r, TReadId& patid);

	/**
	 * Read a pair of patterns from 1 Qseq file.  Note: this is never used.
	 */
	virtual void readPair(Read& ra, Read& rb, TReadId& patid) {
		// (For now, we shouldn't ever be here)
		cerr << "In QseqPatternSource.readPair()" << endl;
		throw 1;
		read(ra, patid);
		readCnt_--;
		read(rb, patid);
	}

	/**
	 * Dump a FASTQ-style record for the read.
	 */
	virtual void dump(ostream& out,
	                  const BTDnaString& seq,
	                  const BTString& qual,
	                  const BTString& name)
	{
		out << "@" << name << endl << seq << endl
		    << "+" << endl << qual << endl;
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
		BufferedFilePatternSource(infiles, NULL, p),
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
	virtual void read(Read& r, TReadId& patid) {
		r.reset();
		while(true) {
			r.color = gColor;
			int c = fb_.get();
			if(c < 0) {
				r.patFw.clear();
				return;
			}
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
						// into the reference; that let's us see where
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
					patid = readCnt_-1;
					break;
				}
			}
		}
	}
	/// Shouldn't ever be here; it's not sensible to obtain read pairs
	// from a continuous input.
	virtual void readPair(Read& ra, Read& rb, TReadId& patid) {
		cerr << "In FastaContinuousPatternSource.readPair()" << endl;
		throw 1;
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
		BufferedFilePatternSource(infiles, NULL, p),
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
	virtual void read(Read& r, TReadId& patid) {
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
			if(c < 0) { bail(r); return; }
		}
		int trim5 = mytrim5;
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
	/// Read another read pair from a FASTQ input file
	virtual void readPair(Read& ra, Read& rb, TReadId& patid) {
		// (For now, we shouldn't ever be here)
		cerr << "In FastqPatternSource.readPair()" << endl;
		throw 1;
		read(ra, patid);
		readCnt_--;
		read(rb, patid);
	}
	virtual void resetForNextFile() {
		first_ = true;
	}
	virtual void dump(ostream& out,
	                  const BTDnaString& seq,
	                  const BTString& qual,
	                  const BTString& name)
	{
		out << "@" << name << endl << seq << endl << "+" << endl << qual << endl;
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
		BufferedFilePatternSource(infiles, NULL, p), first_(true) { }
	virtual void reset() {
		first_ = true;
		BufferedFilePatternSource::reset();
	}
protected:
	/// Read another pattern from a Raw input file
	virtual void read(Read& r, TReadId& patid) {
		int c;
		r.reset();
		c = getOverNewline(this->fb_);
		if(c < 0) { bail(r); return; }
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
			if(c < 0) { bail(r); return; }
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

		patid = readCnt_-1;
	}
	/// Read another read pair from a FASTQ input file
	virtual void readPair(Read& ra, Read& rb, TReadId& patid) {
		// (For now, we shouldn't ever be here)
		cerr << "In RawPatternSource.readPair()" << endl;
		throw 1;
		read(ra, patid);
		readCnt_--;
		read(rb, patid);
	}
	virtual void resetForNextFile() {
		first_ = true;
	}
	virtual void dump(ostream& out,
	                  const BTDnaString& seq,
	                  const BTString& qual,
	                  const BTString& name)
	{
		out << seq << endl;
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

/**
 * Read a Raw-format file (one sequence per line).  No quality strings
 * allowed.  All qualities are assumed to be 'I' (40 on the Phred-33
 * scale).
 */
class ChainPatternSource : public BufferedFilePatternSource {
public:
	ChainPatternSource(const EList<string>& infiles, const PatternParams& p) :
	BufferedFilePatternSource(infiles, NULL, p) { }

protected:

	/// Read another pattern from a Raw input file
	virtual void read(Read& r, TReadId& patid) {
		r.reset();
		fb_.peek();
		if(fb_.eof()) {
			fb_.resetLastN();
			return;
		}
		do {
			r.hitset->deserialize(fb_);
		} while(!r.hitset->initialized() && !fb_.eof());
		if(!r.hitset->initialized()) {
			fb_.resetLastN();
			r.patFw.clear();
			return;
		}
		// Now copy the name/sequence/quals into r.name/r.patFw/r.qualFw
		r.hitset->name = r.name;
		r.hitset->seq  = r.patFw;
		r.hitset->qual = r.qual;

		r.readOrigBuf.install(fb_.lastN(), fb_.lastNLen());
		fb_.resetLastN();

		readCnt_++;
		patid = readCnt_-1;
	}

	/// Read another read pair
	virtual void readPair(Read& ra, Read& rb, TReadId& patid) {
		// (For now, we shouldn't ever be here)
		cerr << "In ChainPatternSource.readPair()" << endl;
		throw 1;
	}
};

#endif /*PAT_H_*/
