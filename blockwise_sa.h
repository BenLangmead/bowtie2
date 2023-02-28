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

#ifndef BLOCKWISE_SA_H_
#define BLOCKWISE_SA_H_

#include <stdint.h>
#include <stdlib.h>
#include <future>
#include <iostream>
#include <sstream>
#include <thread>
#include <memory>
#include <stdexcept>

#ifdef USE_SAIS
#ifdef BOWTIE_64BIT_INDEX
#include <libsais64.h>
#else
#include <libsais.h>
#endif
#endif

#include "assert_helpers.h"
#include "diff_sample.h"
#include "multikey_qsort.h"
#include "random_source.h"
#include "binary_sa_search.h"
#include "zbox.h"
#include "alphabet.h"
#include "timer.h"
#include "ds.h"
#include "mem_ids.h"
#include "word_io.h"
#include "threadpool.h"

using namespace std;

// Helpers for printing verbose messages

#ifndef VMSG_NL
#define VMSG_NL(...)				\
	if(this->verbose()) {			\
		stringstream tmp;		\
		tmp << __VA_ARGS__ << endl;	\
		this->verbose(tmp.str());	\
	}
#endif

#ifndef VMSG
#define VMSG(...)				\
	if(this->verbose()) {			\
		stringstream tmp;		\
		tmp << __VA_ARGS__;		\
		this->verbose(tmp.str());	\
	}
#endif

/**
 * Abstract parent class for blockwise suffix-array building schemes.
 */
template<typename TStr>
class BlockwiseSA {
public:
	BlockwiseSA(const TStr& __text,
	            TIndexOffU __bucketSz,
	            bool __sanityCheck = false,
	            bool __passMemExc = false,
	            bool __verbose = false,
	            ostream& __logger = cout) :
		_text(__text),
		_bucketSz(max<TIndexOffU>(__bucketSz, 2u)),
		_sanityCheck(__sanityCheck),
		_passMemExc(__passMemExc),
		_verbose(__verbose),
		_itrBucket(EBWTB_CAT),
		_itrBucketPos(OFF_MASK),
		_itrPushedBackSuffix(OFF_MASK),
		_logger(__logger)
		{ }

	virtual ~BlockwiseSA() { }

	/**
	 * Get the next suffix; compute the next bucket if necessary.
	 */
	virtual TIndexOffU nextSuffix() = 0;

	/**
	 * Return true iff the next call to nextSuffix will succeed.
	 */
	bool hasMoreSuffixes() {
		if(_itrPushedBackSuffix != OFF_MASK) return true;
		try {
			_itrPushedBackSuffix = nextSuffix();
		} catch(out_of_range& e) {
			assert_eq(OFF_MASK, _itrPushedBackSuffix);
			return false;
		}
		return true;
	}

	/**
	 * Reset the suffix iterator so that the next call to nextSuffix()
	 * returns the lexicographically-first suffix.
	 */
	void resetSuffixItr() {
		_itrBucket.clear();
		_itrBucketPos = OFF_MASK;
		_itrPushedBackSuffix = OFF_MASK;
		reset();
		assert(suffixItrIsReset());
	}

	/**
	 * Returns true iff the next call to nextSuffix() returns the
	 * lexicographically-first suffix.
	 */
	bool suffixItrIsReset() {
		return _itrBucket.size()    == 0 &&
			_itrBucketPos        == OFF_MASK &&
			_itrPushedBackSuffix == OFF_MASK &&
			isReset();
	}

	const TStr& text()  const { return _text; }
	TIndexOffU bucketSz() const { return _bucketSz; }
	bool sanityCheck()  const { return _sanityCheck; }
	bool verbose()      const { return _verbose; }
	ostream& log()      const { return _logger; }
	size_t size()       const { return _text.length()+1; }

protected:
	/// Reset back to the first block
	virtual void reset() = 0;
	/// Return true iff reset to the first block
	virtual bool isReset() = 0;

	/**
	 * Grab the next block of sorted suffixes.  The block is guaranteed
	 * to have at most _bucketSz elements.
	 */
	virtual void nextBlock(int cur_block, int tid = 0) = 0;
	/// Return true iff more blocks are available
	virtual bool hasMoreBlocks() const = 0;
	/// Optionally output a verbose message
	void verbose(const string& s) const {
		if(this->verbose()) {
			this->log() << s.c_str();
			this->log().flush();
		}
	}

	const TStr&      _text;        /// original string
	const TIndexOffU   _bucketSz;    /// target maximum bucket size
	const bool       _sanityCheck; /// whether to perform sanity checks
	const bool       _passMemExc;  /// true -> pass on memory exceptions
	const bool       _verbose;     /// be talkative
	EList<TIndexOffU>  _itrBucket;   /// current bucket
	TIndexOffU         _itrBucketPos;/// offset into current bucket
	TIndexOffU         _itrPushedBackSuffix; /// temporary slot for lookahead
	ostream&         _logger;      /// write log messages here
};

/**
 * Abstract parent class for a blockwise suffix array builder that
 * always doles out blocks in lexicographical order.
 */
template<typename TStr>
class InorderBlockwiseSA : public BlockwiseSA<TStr> {
public:
	InorderBlockwiseSA(const TStr& __text,
	                   TIndexOffU __bucketSz,
	                   bool __sanityCheck = false,
			   bool __passMemExc = false,
	                   bool __verbose = false,
	                   ostream& __logger = cout) :
		BlockwiseSA<TStr>(__text, __bucketSz, __sanityCheck, __passMemExc, __verbose, __logger)
		{ }
};
#ifdef USE_SAIS
template<typename TStr>
class SAISBlockwiseSA : public InorderBlockwiseSA<TStr> {
public:
	SAISBlockwiseSA(TStr& __text, TIndexOffU __bucketSz, int __nthreads) :
		InorderBlockwiseSA<TStr>(__text, __bucketSz, false, false, false, cout), _i(0), _nthreads(__nthreads)
		{
			__text.wbuf()[__text.length()] = (char)127; // $ is larger than any character in the suffix array
			_suffixes.resize(__text.length() + 1);
#ifdef BOWTIE_64BIT_INDEX
			libsais64_omp((const uint8_t *)__text.buf(), (int64_t *)(_suffixes.ptr()), (int64_t)_suffixes.size(), 0, NULL, _nthreads);

#else
			libsais_omp((const uint8_t *)__text.buf(), (int32_t *)(_suffixes.ptr()), (int32_t)_suffixes.size(), 0, NULL, _nthreads);

#endif
			_suffixes[__text.length()] = __text.length();
		}


	inline TIndexOffU nextSuffix() {
		return _suffixes[_i++];
	}

	void reset() {
		_i = 0;
	}

	bool isReset() {
		return _i == 0;
	}

	void nextBlock(int a, int b) {
		return;
	}

	bool hasMoreBlocks() const {
		return _i == _suffixes.size();
	}

private:
	TIndexOffU _i;
	EList<TIndexOffU> _suffixes;
	int _nthreads;
};
#endif
/**
 * Build the SA a block at a time according to the scheme outlined in
 * Karkkainen's "Fast BWT" paper.
 */
template<typename TStr>
class KarkkainenBlockwiseSA : public InorderBlockwiseSA<TStr> {
public:
	typedef DifferenceCoverSample<TStr> TDC;

	KarkkainenBlockwiseSA(const TStr& __text,
			      TIndexOffU __bucketSz,
			      int __nthreads,
			      thread_pool& pool,
			      uint32_t __dcV,
			      uint32_t __seed = 0,
			      bool __sanityCheck = false,
			      bool __passMemExc = false,
			      bool __verbose = false,
			      string base_fname = "",
			      ostream& __logger = cout) :
		InorderBlockwiseSA<TStr>(__text, __bucketSz, __sanityCheck, __passMemExc, __verbose, __logger),
		_sampleSuffs(EBWTB_CAT),
		_nthreads(__nthreads),
		_pool(pool),
		_itrBucketIdx(0),
		_cur(0),
		_dcV(__dcV),
		_dc(EBWTB_CAT),
		_built(false),
		_base_fname(base_fname),
		_bigEndian(currentlyBigEndian()),
		_done(NULL)
		{ _randomSrc.init(__seed); reset(); }

	~KarkkainenBlockwiseSA() throw() {}

	/**
	 * Allocate an amount of memory that simulates the peak memory
	 * usage of the DifferenceCoverSample with the given text and v.
	 * Throws bad_alloc if it's not going to fit in memory.  Returns
	 * the approximate number of bytes the Cover takes at all times.
	 */
	static size_t simulateAllocs(const TStr& text, TIndexOffU bucketSz) {
		size_t len = text.length();
		// _sampleSuffs and _itrBucket are in memory at the peak
		size_t bsz = bucketSz;
		size_t sssz = len / max<TIndexOffU>(bucketSz-1, 1);
		AutoArray<TIndexOffU> tmp(bsz + sssz + (1024 * 1024 /*out of caution*/), EBWT_CAT);
		return bsz;
	}


	/**
	 * Get the next suffix; compute the next bucket if necessary.
	 */
	virtual TIndexOffU nextSuffix()
		{
			// Launch threads if not
			if(this->_nthreads > 1) {
				if(_tparams.size() == 0) {
					_done = new volatile bool[_sampleSuffs.size() + 1];
					for (size_t i = 0; i < _sampleSuffs.size() + 1; i++) {
						_done[i] = false;
					}
					_itrBuckets.resize(this->_nthreads);
					_tparams.resize(this->_nthreads);
					for(int tid = 0; tid < this->_nthreads; tid++) {
						_tparams[tid].first = this;
						_tparams[tid].second = tid;
						if (tid == _nthreads - 1)
							nextBlock_Worker((void *)&_tparams[tid]);
						else
							_pool.submit(nextBlock_Worker((void *)&_tparams[tid]));
					}
				}
			}
			if(this->_itrPushedBackSuffix != OFF_MASK) {
				TIndexOffU tmp = this->_itrPushedBackSuffix;
				this->_itrPushedBackSuffix = OFF_MASK;
				return tmp;
			}
			while(this->_itrBucketPos >= this->_itrBucket.size() ||
			      this->_itrBucket.size() == 0)
			{
				if(!hasMoreBlocks()) {
					throw out_of_range("No more suffixes");
				}
				if(this->_nthreads == 1) {
					nextBlock((int)_cur);
					_cur++;
				} else {
					while(!_done[this->_itrBucketIdx]) {
						SLEEP(1);
					}
					// Read suffixes from a file
					std::ostringstream number; number << this->_itrBucketIdx;
					const string fname = _base_fname + "." + number.str() + ".sa";
					ifstream sa_file(fname.c_str(), ios::binary);
					if(!sa_file.good()) {
						cerr << "Could not open file for reading a suffix array: \"" << fname << "\"" << endl;
						throw 1;
					}
					size_t numSAs = readU<TIndexOffU>(sa_file, false /* do not endian swap */);
					this->_itrBucket.resizeExact(numSAs);
					for(size_t i = 0; i < numSAs; i++) {
						this->_itrBucket[i] = readU<TIndexOffU>(sa_file, false);
					}
					sa_file.close();
					std::remove(fname.c_str());
				}
				this->_itrBucketIdx++;
				this->_itrBucketPos = 0;
			}
			return this->_itrBucket[this->_itrBucketPos++];
		}

	/// Defined in blockwise_sa.cpp
	virtual void nextBlock(int cur_block, int tid = 0);

	/// Defined in blockwise_sa.cpp
	virtual void qsort(EList<TIndexOffU>& bucket);

	/// Return true iff more blocks are available
	virtual bool hasMoreBlocks() const {
		return this->_itrBucketIdx <= _sampleSuffs.size();
	}

	/// Return the difference-cover period
	uint32_t dcV() const { return _dcV; }

//TBB requires a Functor to be passed to the thread group
//hence the nested class

	class nextBlock_Worker {
		void *vp;

	public:
		nextBlock_Worker(const nextBlock_Worker& W): vp(W.vp) {};
		nextBlock_Worker(void *vp_):vp(vp_) {};
		void operator()() const {
			pair<KarkkainenBlockwiseSA*, int> param = *(pair<KarkkainenBlockwiseSA*, int>*)vp;
			KarkkainenBlockwiseSA* sa = param.first;
			int tid = param.second;
			while(true) {
				size_t cur = 0;
				{
					ThreadSafe ts(sa->_mutex);

					cur = sa->_cur;
					if(cur > sa->_sampleSuffs.size()) break;
					sa->_cur++;
				}
				sa->nextBlock((int)cur, tid);
				// Write suffixes into a file
				std::ostringstream number; number << cur;
				const string fname = sa->_base_fname + "." + number.str() + ".sa";
				ofstream sa_file(fname.c_str(), ios::binary);
				if(!sa_file.good()) {
					cerr << "Could not open file for writing a reference graph: \"" << fname << "\"" << endl;
					throw 1;
				}
				const EList<TIndexOffU>& bucket = sa->_itrBuckets[tid];
				writeU<TIndexOffU>(sa_file, (TIndexOffU)bucket.size(), sa->_bigEndian);
				for(size_t i = 0; i < bucket.size(); i++) {
					writeU<TIndexOffU>(sa_file, bucket[i], sa->_bigEndian);
				}
				sa_file.close();
				sa->_itrBuckets[tid].clear();
				sa->_done[cur] = true;
			}
		}
	};

protected:

	/**
	 * Initialize the state of the blockwise suffix sort.  If the
	 * difference cover sample and the sample set have not yet been
	 * built, build them.  Then reset the block cursor to point to
	 * the first block.
	 */
	virtual void reset() {
		if(!_built) {
			build();
		}
		assert(_built);
		_cur = 0;
	}

	/// Return true iff we're about to dole out the first bucket
	virtual bool isReset() {
		return _cur == 0;
	}

private:

	/**
	 * Calculate the difference-cover sample and sample suffixes.
	 */
	void build() {
		// Calculate difference-cover sample
		assert(_dc.get() == NULL);
		if(_dcV != 0) {
			_dc.init(new TDC(this->text(), _dcV, this->verbose(), this->sanityCheck()));
			_dc.get()->build(_pool, this->_nthreads);
		}
		// Calculate sample suffixes
		if(this->bucketSz() <= this->text().length()) {
			VMSG_NL("Building samples");
			buildSamples();
		} else {
			VMSG_NL("Skipping building samples since text length " <<
				this->text().length() << " is less than bucket size: " <<
				this->bucketSz());
		}
		_built = true;
	}

	/**
	 * Calculate the lcp between two suffixes using the difference
	 * cover as a tie-breaker.  If the tie-breaker is employed, then
	 * the calculated lcp may be an underestimate.
	 *
	 * Defined in blockwise_sa.cpp
	 */
	inline bool tieBreakingLcp(TIndexOffU aOff,
	                           TIndexOffU bOff,
	                           TIndexOffU& lcp,
	                           bool& lcpIsSoft);

	/**
	 * Compare two suffixes using the difference-cover sample.
	 */
	inline bool suffixCmp(TIndexOffU cmp,
	                      TIndexOffU i,
	                      int64_t& j,
	                      int64_t& k,
	                      bool& kSoft,
	                      const EList<TIndexOffU>& z);

	void buildSamples();

	EList<TIndexOffU>  _sampleSuffs; /// sample suffixes
	int                _nthreads;    /// # of threads
	thread_pool& _pool;
	TIndexOffU         _itrBucketIdx;
	TIndexOffU         _cur;         /// offset to 1st elt of next block
	const uint32_t   _dcV;         /// difference-cover periodicity
	PtrWrap<TDC>     _dc;          /// queryable difference-cover data
	bool             _built;       /// whether samples/DC have been built
	RandomSource     _randomSrc;   /// source of pseudo-randoms

	MUTEX_T                 _mutex;       /// synchronization of output message
	string                  _base_fname;  /// base file name for storing SA blocks
	bool                    _bigEndian;   /// bigEndian?
	EList<std::thread*> _threads;     /// thread list
	EList<pair<KarkkainenBlockwiseSA*, int> > _tparams;
	ELList<TIndexOffU>      _itrBuckets;  /// buckets
	volatile bool *_done;        /// is a block processed?
};



/**
 * Qsort the set of suffixes whose offsets are in 'bucket'.
 */
template<typename TStr>
inline void KarkkainenBlockwiseSA<TStr>::qsort(EList<TIndexOffU>& bucket) {
	const TStr& t = this->text();
	TIndexOffU *s = bucket.ptr();
	size_t slen = bucket.size();
	TIndexOffU len = (TIndexOffU)t.length();
	if(_dc.get() != NULL) {
		// Use the difference cover as a tie-breaker if we have it
		VMSG_NL("  (Using difference cover)");
		// Extract the 'host' array because it's faster to work
		// with than the EList<> container
		const uint8_t *host = (const uint8_t *)t.buf();
		assert(_dc.get() != NULL);
		mkeyQSortSufDcU8(t, host, len, s, slen, *_dc.get(), 4,
		                 this->verbose(), this->sanityCheck());
	} else {
		VMSG_NL("  (Not using difference cover)");
		// We don't have a difference cover - just do a normal
		// suffix sort
		mkeyQSortSuf(t, s, slen, 4,
		             this->verbose(), this->sanityCheck());
	}
}

/**
 * Qsort the set of suffixes whose offsets are in 'bucket'.  This
 * specialization for packed strings does not attempt to extract and
 * operate directly on the host string; the fact that the string is
 * packed means that the array cannot be sorted directly.
 */
template<>
inline void KarkkainenBlockwiseSA<S2bDnaString>::qsort(
	EList<TIndexOffU>& bucket)
{
	const S2bDnaString& t = this->text();
	TIndexOffU *s = bucket.ptr();
	size_t slen = bucket.size();
	size_t len = t.length();
	if(_dc.get() != NULL) {
		// Use the difference cover as a tie-breaker if we have it
		VMSG_NL("  (Using difference cover)");
		// Can't use the text's 'host' array because the backing
		// store for the packed string is not one-char-per-elt.
		mkeyQSortSufDcU8(t, t, len, s, slen, *_dc.get(), 4,
		                 this->verbose(), this->sanityCheck());
	} else {
		VMSG_NL("  (Not using difference cover)");
		// We don't have a difference cover - just do a normal
		// suffix sort
		mkeyQSortSuf(t, s, slen, 4,
		             this->verbose(), this->sanityCheck());
	}
}

template<typename TStr>
struct BinarySortingParam {
	const TStr*              t;
	const EList<TIndexOffU>* sampleSuffs;
	EList<TIndexOffU>        bucketSzs;
	EList<TIndexOffU>        bucketReps;
	size_t                   begin;
	size_t                   end;
};

template<typename TStr>
class  BinarySorting_worker {
        void *vp;

public:
	BinarySorting_worker(const BinarySorting_worker& W): vp(W.vp) {};
	BinarySorting_worker(void *vp_):vp(vp_) {};
	void operator()() const {
		BinarySortingParam<TStr>* param = (BinarySortingParam<TStr>*)vp;
		const TStr& t = *(param->t);
		size_t len = t.length();
		const EList<TIndexOffU>& sampleSuffs = *(param->sampleSuffs);
		EList<TIndexOffU>& bucketSzs = param->bucketSzs;
		EList<TIndexOffU>& bucketReps = param->bucketReps;
		ASSERT_ONLY(size_t numBuckets = bucketSzs.size());
		size_t begin = param->begin;
		size_t end = param->end;
		// Iterate through every suffix in the text, determine which
		// bucket it falls into by doing a binary search across the
		// sorted list of samples, and increment a counter associated
		// with that bucket.  Also, keep one representative for each
		// bucket so that we can split it later.  We loop in ten
		// stretches so that we can print out a helpful progress
		// message.  (This step can take a long time.)
		for(TIndexOffU i = (TIndexOffU)begin; i < end && i < len; i++) {
			TIndexOffU r = binarySASearch(t, i, sampleSuffs);
			if(r == std::numeric_limits<TIndexOffU>::max()) continue; // r was one of the samples
			assert_lt(r, numBuckets);
			bucketSzs[r]++;
			assert_lt(bucketSzs[r], len);
			if(bucketReps[r] == OFF_MASK || (i & 100) == 0) {
				bucketReps[r] = i; // clobbers previous one, but that's OK
			}
		}
	}

};

/**
 * Select a set of bucket-delineating sample suffixes such that no
 * bucket is greater than the requested upper limit.  Some care is
 * taken to make each bucket's size close to the limit without
 * going over.
 */
template<typename TStr>
void KarkkainenBlockwiseSA<TStr>::buildSamples() {
	const TStr& t = this->text();
	TIndexOffU bsz = this->bucketSz()-1; // subtract 1 to leave room for sample
	size_t len = this->text().length();
	// Prepare _sampleSuffs array
	_sampleSuffs.clear();
	TIndexOffU numSamples = (TIndexOffU)((len/bsz)+1)<<1; // ~len/bsz x 2
	assert_gt(numSamples, 0);
	VMSG_NL("Reserving space for " << numSamples << " sample suffixes");
	if(this->_passMemExc) {
		_sampleSuffs.resizeExact(numSamples);
		// Randomly generate samples.  Allow duplicates for now.
		VMSG_NL("Generating random suffixes");
		for(size_t i = 0; i < numSamples; i++) {
#ifdef BOWTIE_64BIT_INDEX
			_sampleSuffs[i] = (TIndexOffU)(_randomSrc.nextU64() % len);
#else
			_sampleSuffs[i] = (TIndexOffU)(_randomSrc.nextU32() % len);
#endif
		}
	} else {
		try {
			_sampleSuffs.resizeExact(numSamples);
			// Randomly generate samples.  Allow duplicates for now.
			VMSG_NL("Generating random suffixes");
			for(size_t i = 0; i < numSamples; i++) {
#ifdef BOWTIE_64BIT_INDEX
				_sampleSuffs[i] = (TIndexOffU)(_randomSrc.nextU64() % len);
#else
				_sampleSuffs[i] = (TIndexOffU)(_randomSrc.nextU32() % len);
#endif
			}
		} catch(bad_alloc &e) {
			if(this->_passMemExc) {
				throw e; // rethrow immediately
			} else {
				cerr << "Could not allocate sample suffix container of " << (numSamples * OFF_SIZE) << " bytes." << endl
				     << "Please try using a smaller number of blocks by specifying a larger --bmax or" << endl
				     << "a smaller --bmaxdivn" << endl;
				throw 1;
			}
		}
	}
	// Remove duplicates; very important to do this before the call to
	// mkeyQSortSuf so that it doesn't try to calculate lexicographical
	// relationships between very long, identical strings, which takes
	// an extremely long time in general, and causes the stack to grow
	// linearly with the size of the input
	{
		Timer timer(cout, "QSorting sample offsets, eliminating duplicates time: ", this->verbose());
		VMSG_NL("QSorting " << _sampleSuffs.size() << " sample offsets, eliminating duplicates");
		_sampleSuffs.sort();
		size_t sslen = _sampleSuffs.size();
		for(size_t i = 0; i < sslen-1; i++) {
			if(_sampleSuffs[i] == _sampleSuffs[i+1]) {
				_sampleSuffs.erase(i--);
				sslen--;
			}
		}
	}
	// Multikey quicksort the samples
	{
		Timer timer(cout, "  Multikey QSorting samples time: ", this->verbose());
		VMSG_NL("Multikey QSorting " << _sampleSuffs.size() << " samples");
		this->qsort(_sampleSuffs);
	}
	// Calculate bucket sizes
	VMSG_NL("Calculating bucket sizes");
	int limit = 5;
	// Iterate until all buckets are less than
	while(--limit >= 0) {
		TIndexOffU numBuckets = (TIndexOffU)_sampleSuffs.size()+1;
		std::vector<std::future<void> > threads(_pool.size());
		EList<BinarySortingParam<TStr> > tparams;
		tparams.resize(this->_nthreads);
		for(int tid = 0; tid < this->_nthreads; tid++) {
			// Calculate bucket sizes by doing a binary search for each
			// suffix and noting where it lands
			try {
				// Allocate and initialize containers for holding bucket
				// sizes and representatives.
				tparams[tid].bucketSzs.resizeExact(numBuckets);
				tparams[tid].bucketReps.resizeExact(numBuckets);
				tparams[tid].bucketSzs.fillZero();
				tparams[tid].bucketReps.fill(OFF_MASK);
			} catch(bad_alloc &e) {
				if(this->_passMemExc) {
					throw e; // rethrow immediately
				} else {
					cerr << "Could not allocate sizes, representatives (" << ((numBuckets*8)>>10) << " KB) for blocks." << endl
					     << "Please try using a smaller number of blocks by specifying a larger --bmax or a" << endl
					     << "smaller --bmaxdivn." << endl;
					throw 1;
				}
			}
			tparams[tid].t = &t;
			tparams[tid].sampleSuffs = &_sampleSuffs;
			tparams[tid].begin = (tid == 0 ? 0 : len / this->_nthreads * tid);
			tparams[tid].end = (tid + 1 == this->_nthreads ? len : len / this->_nthreads * (tid + 1));
			if(this->_nthreads == 1 || tid == _nthreads - 1) {
				BinarySorting_worker<TStr>((void*)&tparams[tid])();
			} else {
				threads[tid] = _pool.submit(BinarySorting_worker<TStr>(((void*)&tparams[tid])));
			}
		}

		if(this->_nthreads > 1) {
			for (int tid = 0; tid < _pool.size(); tid++) {
				threads[tid].get();
			}
		}
		EList<TIndexOffU>& bucketSzs = tparams[0].bucketSzs;
		EList<TIndexOffU>& bucketReps = tparams[0].bucketReps;
		for(int tid = 1; tid < this->_nthreads; tid++) {
			for(size_t j = 0; j < numBuckets; j++) {
				bucketSzs[j] += tparams[tid].bucketSzs[j];
				if(bucketReps[j] == OFF_MASK) {
					bucketReps[j] = tparams[tid].bucketReps[j];
				}
			}
		}
		// Check for large buckets and mergeable pairs of small buckets
		// and split/merge as necessary
		TIndexOff added = 0;
		TIndexOff merged = 0;
		assert_eq(bucketSzs.size(), numBuckets);
		assert_eq(bucketReps.size(), numBuckets);
		{
			Timer timer(cout, "  Splitting and merging time: ", this->verbose());
			VMSG_NL("Splitting and merging");
			for(TIndexOffU i = 0; i < numBuckets; i++) {
				TIndexOffU mergedSz = bsz + 1;
				assert(bucketSzs[(size_t)i] == 0 || bucketReps[(size_t)i] != OFF_MASK);
				if(i < numBuckets-1) {
					mergedSz = bucketSzs[(size_t)i] + bucketSzs[(size_t)i+1] + 1;
				}
				// Merge?
				if(mergedSz <= bsz) {
					bucketSzs[(size_t)i+1] += (bucketSzs[(size_t)i]+1);
					// The following may look strange, but it's necessary
					// to ensure that the merged bucket has a representative
					bucketReps[(size_t)i+1] = _sampleSuffs[(size_t)i+added];
					_sampleSuffs.erase((size_t)i+added);
					bucketSzs.erase((size_t)i);
					bucketReps.erase((size_t)i);
					i--; // might go to -1 but ++ will overflow back to 0
					numBuckets--;
					merged++;
					assert_eq(numBuckets, _sampleSuffs.size()+1-added);
					assert_eq(numBuckets, bucketSzs.size());
				}
				// Split?
				else if(bucketSzs[(size_t)i] > bsz) {
					// Add an additional sample from the bucketReps[]
					// set accumulated in the binarySASearch loop; this
					// effectively splits the bucket
					_sampleSuffs.insert(bucketReps[(size_t)i], (TIndexOffU)(i + (added++)));
				}
			}
		}
		if(added == 0) {
			//if(this->verbose()) {
			//	cout << "Final bucket sizes:" << endl;
			//	cout << "  (begin): " << bucketSzs[0] << " (" << (int)(bsz - bucketSzs[0]) << ")" << endl;
			//	for(uint32_t i = 1; i < numBuckets; i++) {
			//		cout << "  " << bucketSzs[i] << " (" << (int)(bsz - bucketSzs[i]) << ")" << endl;
			//	}
			//}
			break;
		}
		// Otherwise, continue until no more buckets need to be
		// split
		VMSG_NL("Split " << added << ", merged " << merged << "; iterating...");
	}
	// Do *not* force a do-over
	//	if(limit == 0) {
	//		VMSG_NL("Iterated too many times; trying again...");
	//		buildSamples();
	//	}
	VMSG_NL("Avg bucket size: " << ((double)(len-_sampleSuffs.size()) / (_sampleSuffs.size()+1)) << " (target: " << bsz << ")");
}

/**
 * Do a simple LCP calculation on two strings.
 */
template<typename T> inline
static TIndexOffU suffixLcp(const T& t, TIndexOffU aOff, TIndexOffU bOff) {
	TIndexOffU c = 0;
	size_t len = t.length();
	assert_leq(aOff, len);
	assert_leq(bOff, len);
	while(aOff + c < len && bOff + c < len && t[aOff + c] == t[bOff + c]) c++;
	return c;
}

/**
 * Calculate the lcp between two suffixes using the difference
 * cover as a tie-breaker.  If the tie-breaker is employed, then
 * the calculated lcp may be an underestimate.  If the tie-breaker is
 * employed, lcpIsSoft will be set to true (otherwise, false).
 */
template<typename TStr> inline
bool KarkkainenBlockwiseSA<TStr>::tieBreakingLcp(TIndexOffU aOff,
                                                 TIndexOffU bOff,
                                                 TIndexOffU& lcp,
                                                 bool& lcpIsSoft)
{
	const TStr& t = this->text();
	TIndexOffU c = 0;
	TIndexOffU tlen = (TIndexOffU)t.length();
	assert_leq(aOff, tlen);
	assert_leq(bOff, tlen);
	assert(_dc.get() != NULL);
	uint32_t dcDist = _dc.get()->tieBreakOff(aOff, bOff);
	lcpIsSoft = false; // hard until proven soft
	while(c < dcDist &&    // we haven't hit the tie breaker
	      c < tlen-aOff && // we haven't fallen off of LHS suffix
	      c < tlen-bOff && // we haven't fallen off of RHS suffix
	      t[aOff+c] == t[bOff+c]) // we haven't hit a mismatch
		c++;
	lcp = c;
	if(c == tlen-aOff) {
		// Fell off LHS (a), a is greater
		return false;
	} else if(c == tlen-bOff) {
		// Fell off RHS (b), b is greater
		return true;
	} else if(c == dcDist) {
		// Hit a tie-breaker element
		lcpIsSoft = true;
		assert_neq(dcDist, 0xffffffff);
		return _dc.get()->breakTie(aOff+c, bOff+c) < 0;
	} else {
		assert_neq(t[aOff+c], t[bOff+c]);
		return t[aOff+c] < t[bOff+c];
	}
}

/**
 * Lookup a suffix LCP in the given z array; if the element is not
 * filled in then calculate it from scratch.
 */
template<typename T>
static TIndexOffU lookupSuffixZ(
	const T& t,
	TIndexOffU zOff,
	TIndexOffU off,
	const EList<TIndexOffU>& z)
{
	if(zOff < z.size()) {
		TIndexOffU ret = z[zOff];
		assert_eq(ret, suffixLcp(t, off + zOff, off));
		return ret;
	}
	assert_leq(off + zOff, t.length());
	return suffixLcp(t, off + zOff, off);
}

/**
 * true -> i < cmp
 * false -> i > cmp
 */
template<typename TStr> inline
bool KarkkainenBlockwiseSA<TStr>::suffixCmp(
	TIndexOffU cmp,
	TIndexOffU i,
	int64_t& j,
	int64_t& k,
	bool& kSoft,
	const EList<TIndexOffU>& z)
{
	const TStr& t = this->text();
	TIndexOffU len = (TIndexOffU)t.length();
	// i is not covered by any previous match
	TIndexOffU l;
	if((int64_t)i > k) {
		k = i; // so that i + lHi == kHi
		l = 0; // erase any previous l
		kSoft = false;
		// To be extended
	}
	// i is covered by a previous match
	else /* i <= k */ {
		assert_gt((int64_t)i, j);
		TIndexOffU zIdx = (TIndexOffU)(i-j);
		assert_leq(zIdx, len-cmp);
		if(zIdx < _dcV || _dc.get() == NULL) {
			// Go as far as the Z-box says
			l = lookupSuffixZ(t, zIdx, cmp, z);
			if(i + l > len) {
				l = len-i;
			}
			assert_leq(i + l, len);
			// Possibly to be extended
		} else {
			// But we're past the point of no-more-Z-boxes
			bool ret = tieBreakingLcp(i, cmp, l, kSoft);
			// Sanity-check tie-breaker
			if(this->sanityCheck()) {
				if(ret) assert(sstr_suf_lt(t, i, t, cmp, false));
				else    assert(sstr_suf_gt(t, i, t, cmp, false));
			}
			j = i;
			k = i + l;
			if(this->sanityCheck()) {
				if(kSoft) { assert_leq(l, suffixLcp(t, i, cmp)); }
				else      { assert_eq (l, suffixLcp(t, i, cmp)); }
			}
			return ret;
		}
	}

	// Z box extends exactly as far as previous match (or there
	// is neither a Z box nor a previous match)
	if((int64_t)(i + l) == k) {
		// Extend
		while(l < len-cmp && k < (int64_t)len && t[(size_t)(cmp+l)] == t[(size_t)k]) {
			k++; l++;
		}
		j = i; // update furthest-extending LHS
		kSoft = false;
		assert_eq(l, suffixLcp(t, i, cmp));
	}
	// Z box extends further than previous match
	else if((int64_t)(i + l) > k) {
		l = (TIndexOffU)(k - i); // point to just after previous match
		j = i; // update furthest-extending LHS
		if(kSoft) {
			while(l < len-cmp && k < (int64_t)len && t[(size_t)(cmp+l)] == t[(size_t)k]) {
				k++; l++;
			}
			kSoft = false;
			assert_eq(l, suffixLcp(t, i, cmp));
		} else assert_eq(l, suffixLcp(t, i, cmp));
	}

	// Check that calculated lcp matches actual lcp
	if(this->sanityCheck()) {
		if(!kSoft) {
			// l should exactly match lcp
			assert_eq(l, suffixLcp(t, i, cmp));
		} else {
			// l is an underestimate of LCP
			assert_leq(l, suffixLcp(t, i, cmp));
		}
	}
	assert_leq(l+i, len);
	assert_leq(l, len-cmp);

	// i and cmp should not be the same suffix
	assert(l != len-cmp || i+l != len);

	// Now we're ready to do a comparison on the next char
	if(l+i != len && (
		   l == len-cmp || // departure from paper algorithm:
		   // falling off pattern implies
		   // pattern is *greater* in our case
		   t[i + l] < t[cmp + l]))
	{
		// Case 2: Text suffix is less than upper sample suffix
#ifndef NDEBUG
		if(this->sanityCheck()) {
			assert(sstr_suf_lt(t, i, t, cmp, false));
		}
#endif
		return true; // suffix at i is less than suffix at cmp
	}
	else {
		// Case 3: Text suffix is greater than upper sample suffix
#ifndef NDEBUG
		if(this->sanityCheck()) {
			assert(sstr_suf_gt(t, i, t, cmp, false));
		}
#endif
		return false; // suffix at i is less than suffix at cmp
	}
}

/**
 * Retrieve the next block.  This is the most performance-critical part
 * of the blockwise suffix sorting process.
 */
template<typename TStr>
void KarkkainenBlockwiseSA<TStr>::nextBlock(int cur_block, int tid) {
#ifndef NDEBUG
	if(this->_nthreads > 1) {
		assert_lt(tid, (int)this->_itrBuckets.size());
	}
#endif
	EList<TIndexOffU>& bucket = (this->_nthreads > 1 ? this->_itrBuckets[tid] : this->_itrBucket);
	{
		ThreadSafe ts(_mutex);
		VMSG_NL("Getting block " << (cur_block+1) << " of " << _sampleSuffs.size()+1);
	}
	assert(_built);
	assert_gt(_dcV, 3);
	assert_leq(cur_block, (int)_sampleSuffs.size());
	const TStr& t = this->text();
	TIndexOffU len = (TIndexOffU)t.length();
	// Set up the bucket
	bucket.clear();
	TIndexOffU lo = OFF_MASK, hi = OFF_MASK;
	if(_sampleSuffs.size() == 0) {
		// Special case: if _sampleSuffs is 0, then multikey-quicksort
		// everything
		{
			ThreadSafe ts(_mutex);
			VMSG_NL("  No samples; assembling all-inclusive block");
		}
		assert_eq(0, cur_block);
		try {
			if(bucket.capacity() < this->bucketSz()) {
				bucket.reserveExact(len+1);
			}
			bucket.resize(len);
			for(TIndexOffU i = 0; i < len; i++) {
				bucket[i] = i;
			}
		} catch(bad_alloc &e) {
			if(this->_passMemExc) {
				throw e; // rethrow immediately
			} else {
				cerr << "Could not allocate a master suffix-array block of " << ((len+1) * 4) << " bytes" << endl
				     << "Please try using a larger number of blocks by specifying a smaller --bmax or" << endl
				     << "a larger --bmaxdivn" << endl;
				throw 1;
			}
		}
	} else {
		try {
			{
				ThreadSafe ts(_mutex);
				VMSG_NL("  Reserving size (" << this->bucketSz() << ") for bucket " << (cur_block+1));
			}
			// BTL: Add a +100 fudge factor; there seem to be instances
			// where a bucket ends up having one more elt than bucketSz()
			if(bucket.size() < this->bucketSz()+100) {
				bucket.reserveExact(this->bucketSz()+100);
			}
		} catch(bad_alloc &e) {
			if(this->_passMemExc) {
				throw e; // rethrow immediately
			} else {
				cerr << "Could not allocate a suffix-array block of " << ((this->bucketSz()+1) * 4) << " bytes" << endl;
				cerr << "Please try using a larger number of blocks by specifying a smaller --bmax or" << endl
				     << "a larger --bmaxdivn" << endl;
				throw 1;
			}
		}
		// Select upper and lower bounds from _sampleSuffs[] and
		// calculate the Z array up to the difference-cover periodicity
		// for both.  Be careful about first/last buckets.
		EList<TIndexOffU> zLo(EBWTB_CAT), zHi(EBWTB_CAT);
		assert_geq(cur_block, 0);
		assert_leq((size_t)cur_block, _sampleSuffs.size());
		bool first = (cur_block == 0);
		bool last  = ((size_t)cur_block == _sampleSuffs.size());
		try {
			// Timer timer(cout, "  Calculating Z arrays time: ", this->verbose());
			{
				ThreadSafe ts(_mutex);
				VMSG_NL("  Calculating Z arrays for bucket " << (cur_block+1));
			}
			if(!last) {
				// Not the last bucket
				assert_lt(cur_block, (int)_sampleSuffs.size());
				hi = _sampleSuffs[cur_block];
				zHi.resizeExact(_dcV);
				zHi.fillZero();
				assert_eq(zHi[0], 0);
				calcZ(t, hi, zHi, this->verbose(), this->sanityCheck());
			}
			if(!first) {
				// Not the first bucket
				assert_gt(cur_block, 0);
				assert_leq(cur_block, (int)_sampleSuffs.size());
				lo = _sampleSuffs[cur_block-1];
				zLo.resizeExact(_dcV);
				zLo.fillZero();
				assert_gt(_dcV, 3);
				assert_eq(zLo[0], 0);
				calcZ(t, lo, zLo, this->verbose(), this->sanityCheck());
			}
		} catch(bad_alloc &e) {
			if(this->_passMemExc) {
				throw e; // rethrow immediately
			} else {
				cerr << "Could not allocate a z-array of " << (_dcV * 4) << " bytes" << endl;
				cerr << "Please try using a larger number of blocks by specifying a smaller --bmax or" << endl
				     << "a larger --bmaxdivn" << endl;
				throw 1;
			}
		}

		// This is the most critical loop in the algorithm; this is where
		// we iterate over all suffixes in the text and pick out those that
		// fall into the current bucket.
		//
		// This loop is based on the SMALLERSUFFIXES function outlined on
		// p7 of the "Fast BWT" paper
		//
		int64_t kHi = -1, kLo = -1;
		int64_t jHi = -1, jLo = -1;
		bool kHiSoft = false, kLoSoft = false;
		assert_eq(0, bucket.size());
		{
			// Timer timer(cout, "  Block accumulator loop time: ", this->verbose());
			{
				ThreadSafe ts(_mutex);
				VMSG_NL("  Entering block accumulator loop for bucket " << (cur_block+1) << ":");
			}
			TIndexOffU lenDiv10 = (len + 9) / 10;
			for(TIndexOffU iten = 0, ten = 0; iten < len; iten += lenDiv10, ten++) {
				TIndexOffU itenNext = iten + lenDiv10;
				{
					ThreadSafe ts(_mutex);
					if(ten > 0) VMSG_NL("  bucket " << (cur_block+1) << ": " << (ten * 10) << "%");
				}
				for(TIndexOffU i = iten; i < itenNext && i < len; i++) {
					assert_lt(jLo, (int64_t)i); assert_lt(jHi, (int64_t)i);
					// Advance the upper-bound comparison by one character
					if(i == hi || i == lo) continue; // equal to one of the bookends
					if(hi != OFF_MASK && !suffixCmp(hi, i, jHi, kHi, kHiSoft, zHi)) {
						continue; // not in the bucket
					}
					if(lo != OFF_MASK && suffixCmp(lo, i, jLo, kLo, kLoSoft, zLo)) {
						continue; // not in the bucket
					}
					// In the bucket! - add it
					assert_lt(i, len);
					try {
						bucket.push_back(i);
					} catch(bad_alloc &e) {
						cerr << "Could not append element to block of " << ((bucket.size()) * OFF_SIZE) << " bytes" << endl;
						if(this->_passMemExc) {
							throw e; // rethrow immediately
						} else {
							cerr << "Please try using a larger number of blocks by specifying a smaller --bmax or" << endl
							     << "a larger --bmaxdivn" << endl;
							throw 1;
						}
					}
					// Not necessarily true; we allow overflowing buckets
					// since we can't guarantee that a good set of sample
					// suffixes can be found in a reasonable amount of time
					//assert_lt(bucket.size(), this->bucketSz());
				}
			} // end loop over all suffixes of t
			{
				ThreadSafe ts(_mutex);
				VMSG_NL("  bucket " << (cur_block+1) << ": 100%");
			}
		}
	} // end else clause of if(_sampleSuffs.size() == 0)
	// Sort the bucket
	if(bucket.size() > 0) {
		Timer timer(cout, "  Sorting block time: ", this->verbose());
		{
			ThreadSafe ts(_mutex);
			VMSG_NL("  Sorting block of length " << bucket.size() << " for bucket " << (cur_block+1));
		}
		this->qsort(bucket);
	}
	if(hi != OFF_MASK) {
		// Not the final bucket; throw in the sample on the RHS
		bucket.push_back(hi);
	} else {
		// Final bucket; throw in $ suffix
		bucket.push_back(len);
	}
	{
		ThreadSafe ts(_mutex);
		VMSG_NL("Returning block of " << bucket.size() << " for bucket " << (cur_block+1));
	}
}

#endif /*BLOCKWISE_SA_H_*/
