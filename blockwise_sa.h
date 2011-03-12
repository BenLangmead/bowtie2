#ifndef BLOCKWISE_SA_H_
#define BLOCKWISE_SA_H_

#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <seqan/sequence.h>
#include <seqan/index.h>
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

using namespace std;
using namespace seqan;

// Helpers for printing verbose messages

#ifndef VMSG_NL
#define VMSG_NL(args...) \
if(this->verbose()) { \
	stringstream tmp; \
	tmp << args << endl; \
	this->verbose(tmp.str()); \
}
#endif

#ifndef VMSG
#define VMSG(args...) \
if(this->verbose()) { \
	stringstream tmp; \
	tmp << args; \
	this->verbose(tmp.str()); \
}
#endif

/**
 * Abstract parent class for blockwise suffix-array building schemes.
 */
template<typename TStr>
class BlockwiseSA {
public:
	BlockwiseSA(const TStr& __text,
	            uint32_t __bucketSz,
	            bool __sanityCheck = false,
	            bool __passMemExc = false,
	            bool __verbose = false,
	            ostream& __logger = cout) :
	_text(__text),
	_bucketSz(max<uint32_t>(__bucketSz, 2u)),
	_sanityCheck(__sanityCheck),
	_passMemExc(__passMemExc),
	_verbose(__verbose),
	_itrBucket(),
	_itrBucketPos(0xffffffff),
	_itrPushedBackSuffix(0xffffffff),
	_logger(__logger)
	{ }

	virtual ~BlockwiseSA() { }

	/**
	 * Get the next suffix; compute the next bucket if necessary.
	 */
	uint32_t nextSuffix() {
		if(_itrPushedBackSuffix != 0xffffffff) {
			uint32_t tmp = _itrPushedBackSuffix;
			_itrPushedBackSuffix = 0xffffffff;
			return tmp;
		}
		while(_itrBucketPos >= length(_itrBucket) ||
		      length(_itrBucket) == 0)
		{
			if(!hasMoreBlocks()) {
				throw out_of_range("No more suffixes");
			}
			nextBlock();
			_itrBucketPos = 0;
		}
		return _itrBucket[_itrBucketPos++];
	}

	/**
	 * Return true iff the next call to nextSuffix will succeed.
	 */
	bool hasMoreSuffixes() {
		if(_itrPushedBackSuffix != 0xffffffff) return true;
		try {
			_itrPushedBackSuffix = nextSuffix();
		} catch(out_of_range& e) {
			assert_eq(0xffffffff, _itrPushedBackSuffix);
			return false;
		}
		return true;
	}

	/**
	 * Reset the suffix iterator so that the next call to nextSuffix()
	 * returns the lexicographically-first suffix.
	 */
	void resetSuffixItr() {
		clear(_itrBucket);
		_itrBucketPos = 0xffffffff;
		_itrPushedBackSuffix = 0xffffffff;
		reset();
		assert(suffixItrIsReset());
	}

	/**
	 * Returns true iff the next call to nextSuffix() returns the
	 * lexicographically-first suffix.
	 */
	bool suffixItrIsReset() {
		return length(_itrBucket)   == 0 &&
		       _itrBucketPos        == 0xffffffff &&
		       _itrPushedBackSuffix == 0xffffffff &&
		       isReset();
	}

	const TStr& text()  const { return _text; }
	uint32_t bucketSz() const { return _bucketSz; }
	bool sanityCheck()  const { return _sanityCheck; }
	bool verbose()      const { return _verbose; }
	ostream& log()      const { return _logger; }
	uint32_t size()     const { return length(_text)+1; }

protected:
	/// Reset back to the first block
	virtual void reset() = 0;
	/// Return true iff reset to the first block
	virtual bool isReset() = 0;

	/**
	 * Grab the next block of sorted suffixes.  The block is guaranteed
	 * to have at most _bucketSz elements.
	 */
	virtual void nextBlock() = 0;
	/// Return true iff more blocks are available
	virtual bool hasMoreBlocks() const = 0;
	/// Optionally output a verbose message
	void verbose(const string& s) const {
		if(this->verbose()) {
			this->log() << s;
			this->log().flush();
		}
	}

	const TStr&      _text;        /// original string
	const uint32_t   _bucketSz;    /// target maximum bucket size
	const bool       _sanityCheck; /// whether to perform sanity checks
	const bool       _passMemExc;  /// true -> pass on memory exceptions
	const bool       _verbose;     /// be talkative
	String<uint32_t> _itrBucket;   /// current bucket
	uint32_t         _itrBucketPos;/// offset into current bucket
	uint32_t         _itrPushedBackSuffix; /// temporary slot for lookahead
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
	                   uint32_t __bucketSz,
	                   bool __sanityCheck = false,
	   	               bool __passMemExc = false,
	                   bool __verbose = false,
	                   ostream& __logger = cout) :
	BlockwiseSA<TStr>(__text, __bucketSz, __sanityCheck, __passMemExc, __verbose, __logger)
	{ }
};

/**
 * Build the SA a block at a time according to the scheme outlined in
 * Karkkainen's "Fast BWT" paper.
 */
template<typename TStr>
class KarkkainenBlockwiseSA : public InorderBlockwiseSA<TStr> {
public:
	typedef DifferenceCoverSample<TStr> TDC;

	KarkkainenBlockwiseSA(const TStr& __text,
	                      uint32_t __bucketSz,
	                      uint32_t __dcV,
	                      uint32_t __seed = 0,
	      	              bool __sanityCheck = false,
	   	                  bool __passMemExc = false,
	      	              bool __verbose = false,
	      	              ostream& __logger = cout) :
	InorderBlockwiseSA<TStr>(__text, __bucketSz, __sanityCheck, __passMemExc, __verbose, __logger),
	_sampleSuffs(), _cur(0), _dcV(__dcV), _dc(), _built(false)
	{ _randomSrc.init(__seed); reset(); }

	~KarkkainenBlockwiseSA() { }

	/**
	 * Allocate an amount of memory that simulates the peak memory
	 * usage of the DifferenceCoverSample with the given text and v.
	 * Throws bad_alloc if it's not going to fit in memory.  Returns
	 * the approximate number of bytes the Cover takes at all times.
	 */
	static size_t simulateAllocs(const TStr& text, uint32_t bucketSz) {
		size_t len = length(text);
		// _sampleSuffs and _itrBucket are in memory at the peak
		size_t bsz = bucketSz;
		size_t sssz = len / max<uint32_t>(bucketSz-1, 1);
		AutoArray<uint32_t> tmp(bsz + sssz + (1024 * 1024 /*out of caution*/), EBWT_CAT);
		return bsz;
	}

	/// Defined in blockwise_sa.cpp
	virtual void nextBlock();

	/// Defined in blockwise_sa.cpp
	virtual void qsort(String<uint32_t>& bucket);

	/// Return true iff more blocks are available
	virtual bool hasMoreBlocks() const {
		return _cur <= length(_sampleSuffs);
	}

	/// Return the difference-cover period
	uint32_t dcV() const { return _dcV; }

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
			_dc.get()->build();
		}
		// Calculate sample suffixes
		if(this->bucketSz() <= length(this->text())) {
			VMSG_NL("Building samples");
			buildSamples();
		} else {
			VMSG_NL("Skipping building samples since text length " <<
			        length(this->text()) << " is less than bucket size: " <<
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
	inline bool tieBreakingLcp(uint32_t aOff,
	                           uint32_t bOff,
	                           uint32_t& lcp,
	                           bool& lcpIsSoft);

	/**
	 * Compare two suffixes using the difference-cover sample.
	 */
	inline bool suffixCmp(uint32_t cmp,
	                      uint32_t i,
	                      int64_t& j,
	                      int64_t& k,
	                      bool& kSoft,
	                      const String<uint32_t>& z);

	void buildSamples();

	String<uint32_t> _sampleSuffs; /// sample suffixes
	uint32_t         _cur;         /// offset to 1st elt of next block
	const uint32_t   _dcV;         /// difference-cover periodicity
	PtrWrap<TDC>     _dc;          /// queryable difference-cover data
	bool             _built;       /// whether samples/DC have been built
	RandomSource     _randomSrc;   /// source of pseudo-randoms
};

/**
 * Qsort the set of suffixes whose offsets are in 'bucket'.
 */
template<typename TStr>
inline void KarkkainenBlockwiseSA<TStr>::qsort(String<uint32_t>& bucket) {
	typedef typename Value<TStr>::Type TAlphabet;
	const TStr& t = this->text();
	uint32_t *s = begin(bucket);
	uint32_t slen = seqan::length(bucket);
	uint32_t len = seqan::length(t);
	if(_dc.get() != NULL) {
		// Use the difference cover as a tie-breaker if we have it
		VMSG_NL("  (Using difference cover)");
		// Extract the 'host' array because it's faster to work
		// with than the String<> container
		uint8_t *host = (uint8_t*)t.data_begin;
		assert(_dc.get() != NULL);
		mkeyQSortSufDcU8(t, host, len, s, slen, *_dc.get(),
		                 ValueSize<TAlphabet>::VALUE,
		                 this->verbose(), this->sanityCheck());
	} else {
		VMSG_NL("  (Not using difference cover)");
		// We don't have a difference cover - just do a normal
		// suffix sort
		mkeyQSortSuf(t, s, slen, ValueSize<TAlphabet>::VALUE,
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
inline void KarkkainenBlockwiseSA<String<Dna, Packed<> > >::qsort(String<uint32_t>& bucket) {
	const String<Dna, Packed<> >& t = this->text();
	uint32_t *s = begin(bucket);
	uint32_t slen = seqan::length(bucket);
	uint32_t len = seqan::length(t);
	if(_dc.get() != NULL) {
		// Use the difference cover as a tie-breaker if we have it
		VMSG_NL("  (Using difference cover)");
		// Can't use the text's 'host' array because the backing
		// store for the packed string is not one-char-per-elt.
		mkeyQSortSufDcU8(t, t, len, s, slen, *_dc.get(),
		                 ValueSize<Dna>::VALUE,
		                 this->verbose(), this->sanityCheck());
	} else {
		VMSG_NL("  (Not using difference cover)");
		// We don't have a difference cover - just do a normal
		// suffix sort
		mkeyQSortSuf(t, s, slen, ValueSize<Dna>::VALUE,
		             this->verbose(), this->sanityCheck());
	}
}

/**
 * Select a set of bucket-delineating sample suffixes such that no
 * bucket is greater than the requested upper limit.  Some care is
 * taken to make each bucket's size close to the limit without
 * going over.
 */
template<typename TStr>
void KarkkainenBlockwiseSA<TStr>::buildSamples() {
	typedef typename Value<TStr>::Type TAlphabet;
	const TStr& t = this->text();
	uint32_t bsz = this->bucketSz()-1; // subtract 1 to leave room for sample
	uint32_t len = length(this->text());
	// Prepare _sampleSuffs array
	clear(_sampleSuffs);
	uint32_t numSamples = ((len/bsz)+1)<<1; // ~len/bsz x 2
	assert_gt(numSamples, 0);
	VMSG_NL("Reserving space for " << numSamples << " sample suffixes");
	if(this->_passMemExc) {
		reserve(_sampleSuffs, numSamples, Exact());
		// Randomly generate samples.  Allow duplicates for now.
		VMSG_NL("Generating random suffixes");
		for(size_t i = 0; i < numSamples; i++) {
			appendValue(_sampleSuffs, _randomSrc.nextU32() % len);
		}
	} else {
		try {
			reserve(_sampleSuffs, numSamples, Exact());
			// Randomly generate samples.  Allow duplicates for now.
			VMSG_NL("Generating random suffixes");
			for(size_t i = 0; i < numSamples; i++) {
				appendValue(_sampleSuffs, _randomSrc.nextU32() % len);
			}
		} catch(bad_alloc &e) {
			if(this->_passMemExc) {
				throw e; // rethrow immediately
			} else {
				cerr << "Could not allocate sample suffix container of " << (numSamples * 4) << " bytes." << endl
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
		VMSG_NL("QSorting " << length(_sampleSuffs) << " sample offsets, eliminating duplicates");
		sort(begin(_sampleSuffs), end(_sampleSuffs));
		size_t sslen = length(_sampleSuffs);
		for(size_t i = 0; i < sslen-1; i++) {
			if(_sampleSuffs[i] == _sampleSuffs[i+1]) {
				erase(_sampleSuffs, i--);
				sslen--;
			}
		}
	}
	// Multikey quicksort the samples
	{
		Timer timer(cout, "  Multikey QSorting samples time: ", this->verbose());
		VMSG_NL("Multikey QSorting " << length(_sampleSuffs) << " samples");
		this->qsort(_sampleSuffs);
	}
	// Calculate bucket sizes
	VMSG_NL("Calculating bucket sizes");
	int limit = 5;
	// Iterate until all buckets are less than
	while(--limit >= 0) {
		// Calculate bucket sizes by doing a binary search for each
		// suffix and noting where it lands
		uint32_t numBuckets = length(_sampleSuffs)+1;
		String<uint32_t> bucketSzs; // holds computed bucket sizes
		String<uint32_t> bucketReps; // holds 1 member of each bucket (for splitting)
		try {
			// Allocate and initialize containers for holding bucket
			// sizes and representatives.
			fill(bucketSzs, numBuckets, 0, Exact());
			fill(bucketReps, numBuckets, 0xffffffff, Exact());
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
		// Iterate through every suffix in the text, determine which
		// bucket it falls into by doing a binary search across the
		// sorted list of samples, and increment a counter associated
		// with that bucket.  Also, keep one representative for each
		// bucket so that we can split it later.  We loop in ten
		// stretches so that we can print out a helpful progress
		// message.  (This step can take a long time.)
		{
			VMSG_NL("  Binary sorting into buckets");
			Timer timer(cout, "  Binary sorting into buckets time: ", this->verbose());
			uint32_t lenDiv10 = (len + 9) / 10;
			for(uint32_t iten = 0, ten = 0; iten < len; iten += lenDiv10, ten++) {
				uint32_t itenNext = iten + lenDiv10;
				if(ten > 0) VMSG_NL("  " << (ten * 10) << "%");
				for(uint32_t i = iten; i < itenNext && i < len; i++) {
					uint32_t r = binarySASearch(t, i, _sampleSuffs);
					if(r == 0xffffffff) continue; // r was one of the samples
					assert_lt(r, numBuckets);
					bucketSzs[r]++;
					assert_lt(bucketSzs[r], len);
					if(bucketReps[r] == 0xffffffff ||
					   (_randomSrc.nextU32() & 100) == 0)
					{
						bucketReps[r] = i; // clobbers previous one, but that's OK
					}
				}
			}
			VMSG_NL("  100%");
		}
		// Check for large buckets and mergeable pairs of small buckets
		// and split/merge as necessary
		int added = 0;
		int merged = 0;
		assert_eq(length(bucketSzs), numBuckets);
		assert_eq(length(bucketReps), numBuckets);
		{
			Timer timer(cout, "  Splitting and merging time: ", this->verbose());
			VMSG_NL("Splitting and merging");
			for(int64_t i = 0; i < numBuckets; i++) {
				uint32_t mergedSz = bsz + 1;
				assert(bucketSzs[i] == 0 || bucketReps[i] != 0xffffffff);
				if(i < (int64_t)numBuckets-1) {
					mergedSz = bucketSzs[i] + bucketSzs[i+1] + 1;
				}
				// Merge?
				if(mergedSz <= bsz) {
					bucketSzs[i+1] += (bucketSzs[i]+1);
					// The following may look strange, but it's necessary
					// to ensure that the merged bucket has a representative
					bucketReps[i+1] = _sampleSuffs[i+added];
					erase(_sampleSuffs, i+added);
					erase(bucketSzs, i);
					erase(bucketReps, i);
					i--; // might go to -1 but ++ will overflow back to 0
					numBuckets--;
					merged++;
					assert_eq(numBuckets, length(_sampleSuffs)+1-added);
					assert_eq(numBuckets, length(bucketSzs));
				}
				// Split?
				else if(bucketSzs[i] > bsz) {
					// Add an additional sample from the bucketReps[]
					// set accumulated in the binarySASearch loop; this
					// effectively splits the bucket
					insertValue(_sampleSuffs, i + (added++), bucketReps[i]);
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
	VMSG_NL("Avg bucket size: " << ((float)(len-length(_sampleSuffs)) / (length(_sampleSuffs)+1)) << " (target: " << bsz << ")");
}

/**
 * Do a simple LCP calculation on two strings.
 */
template<typename T> inline
static uint32_t suffixLcp(const T& t, uint32_t aOff, uint32_t bOff) {
	uint32_t c = 0;
	size_t len = length(t);
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
bool KarkkainenBlockwiseSA<TStr>::tieBreakingLcp(uint32_t aOff,
                                                 uint32_t bOff,
                                                 uint32_t& lcp,
                                                 bool& lcpIsSoft)
{
	const TStr& t = this->text();
	uint32_t c = 0;
	uint32_t tlen = length(t);
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
static uint32_t lookupSuffixZ(const T& t,
                              uint32_t zOff,
                              uint32_t off,
                              const String<uint32_t>& z)
{
	if(zOff < length(z)) {
		uint32_t ret = z[zOff];
		assert_eq(ret, suffixLcp(t, off + zOff, off));
		return ret;
	}
	assert_leq(off + zOff, length(t));
	return suffixLcp(t, off + zOff, off);
}

/**
 * true -> i < cmp
 * false -> i > cmp
 */
template<typename TStr> inline
bool KarkkainenBlockwiseSA<TStr>::suffixCmp(uint32_t cmp,
                                            uint32_t i,
                                            int64_t& j,
                                            int64_t& k,
                                            bool& kSoft,
                                            const String<uint32_t>& z)
{
	const TStr& t = this->text();
	uint32_t len = length(t);
	// i is not covered by any previous match
	uint32_t l;
	if(i > k) {
		k = i; // so that i + lHi == kHi
		l = 0; // erase any previous l
		kSoft = false;
		// To be extended
	}
	// i is covered by a previous match
	else /* i <= k */ {
		assert_gt((int64_t)i, j);
		uint32_t zIdx = i-j;
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
				if(ret) assert(dollarLt(suffix(t, i), suffix(t, cmp)));
				else    assert(dollarGt(suffix(t, i), suffix(t, cmp)));
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
	if(i + l == k) {
		// Extend
		while(l < len-cmp && k < len && t[cmp+l] == t[k]) {
			k++; l++;
		}
		j = i; // update furthest-extending LHS
		kSoft = false;
		assert_eq(l, suffixLcp(t, i, cmp));
	}
	// Z box extends further than previous match
	else if(i + l > k) {
		l = k - i; // point to just after previous match
		j = i; // update furthest-extending LHS
		if(kSoft) {
			while(l < len-cmp && k < len && t[cmp+l] == t[k]) {
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
		if(this->sanityCheck()) assert(dollarLt(suffix(t, i), suffix(t, cmp)));
		return true; // suffix at i is less than suffix at cmp
	}
	else {
		// Case 3: Text suffix is greater than upper sample suffix
		if(this->sanityCheck()) assert(dollarGt(suffix(t, i), suffix(t, cmp)));
		return false; // suffix at i is less than suffix at cmp
	}
}

/**
 * Retrieve the next block.  This is the most performance-critical part
 * of the blockwise suffix sorting process.
 */
template<typename TStr>
void KarkkainenBlockwiseSA<TStr>::nextBlock() {
	typedef typename Value<TStr>::Type TAlphabet;
	String<uint32_t>& bucket = this->_itrBucket;
	VMSG_NL("Getting block " << (_cur+1) << " of " << length(_sampleSuffs)+1);
	assert(_built);
	assert_gt(_dcV, 3);
	assert_leq(_cur, length(_sampleSuffs));
	const TStr& t = this->text();
	uint32_t len = length(t);
	// Set up the bucket
	clear(bucket);
	uint32_t lo = 0xffffffff, hi = 0xffffffff;
	if(length(_sampleSuffs) == 0) {
		// Special case: if _sampleSuffs is 0, then multikey-quicksort
		// everything
		VMSG_NL("  No samples; assembling all-inclusive block");
		assert_eq(0, _cur);
		try {
			if(capacity(bucket) < this->bucketSz()) {
				reserve(bucket, len+1, Exact());
			}
			for(uint32_t i = 0; i < len; i++) append(bucket, i);
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
			VMSG_NL("  Reserving size (" << this->bucketSz() << ") for bucket");
			// BTL: Add a +100 fudge factor; there seem to be instances
			// where a bucket ends up having one more elt than bucketSz()
			if(capacity(bucket) < this->bucketSz()+100) {
				reserve(bucket, this->bucketSz()+100, Exact());
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
		String<uint32_t> zLo, zHi;
		assert_geq(_cur, 0);
		assert_leq(_cur, length(_sampleSuffs));
		bool first = (_cur == 0);
		bool last  = (_cur == length(_sampleSuffs));
		try {
			Timer timer(cout, "  Calculating Z arrays time: ", this->verbose());
			VMSG_NL("  Calculating Z arrays");
			if(!last) {
				// Not the last bucket
				assert_lt(_cur, length(_sampleSuffs));
				hi = _sampleSuffs[_cur];
				fill(zHi, _dcV, 0, Exact());
				assert_eq(zHi[0], 0);
				calcZ(t, hi, zHi, this->verbose(), this->sanityCheck());
			}
			if(!first) {
				// Not the first bucket
				assert_gt(_cur, 0);
				assert_leq(_cur, length(_sampleSuffs));
				lo = _sampleSuffs[_cur-1];
				fill(zLo, _dcV, 0, Exact());
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
		assert_eq(0, length(bucket));
		{
			Timer timer(cout, "  Block accumulator loop time: ", this->verbose());
			VMSG_NL("  Entering block accumulator loop:");
			uint32_t lenDiv10 = (len + 9) / 10;
			for(uint32_t iten = 0, ten = 0; iten < len; iten += lenDiv10, ten++) {
			uint32_t itenNext = iten + lenDiv10;
			if(ten > 0) VMSG_NL("  " << (ten * 10) << "%");
			for(uint32_t i = iten; i < itenNext && i < len; i++) {
				assert_lt(jLo, i); assert_lt(jHi, i);
				// Advance the upper-bound comparison by one character
				if(i == hi || i == lo) continue; // equal to one of the bookends
				if(hi != 0xffffffff && !suffixCmp(hi, i, jHi, kHi, kHiSoft, zHi)) {
					continue; // not in the bucket
				}
				if(lo != 0xffffffff && suffixCmp(lo, i, jLo, kLo, kLoSoft, zLo)) {
					continue; // not in the bucket
				}
				// In the bucket! - add it
				assert_lt(i, len);
				try {
					appendValue(bucket, i);
				} catch(bad_alloc &e) {
					if(this->_passMemExc) {
						throw e; // rethrow immediately
					} else {
						cerr << "Could not append element to block of " << ((length(bucket)) * 4) << " bytes" << endl;
						cerr << "Please try using a larger number of blocks by specifying a smaller --bmax or" << endl
						     << "a larger --bmaxdivn" << endl;
						throw 1;
					}
				}
				// Not necessarily true; we allow overflowing buckets
				// since we can't guarantee that a good set of sample
				// suffixes can be found in a reasonable amount of time
				//assert_lt(length(bucket), this->bucketSz());
			}
			} // end loop over all suffixes of t
			VMSG_NL("  100%");
		}
	} // end else clause of if(length(_sampleSuffs) == 0)
	// Sort the bucket
	if(length(bucket) > 0) {
		Timer timer(cout, "  Sorting block time: ", this->verbose());
		VMSG_NL("  Sorting block of length " << length(bucket));
		this->qsort(bucket);
	}
	if(hi != 0xffffffff) {
		// Not the final bucket; throw in the sample on the RHS
		appendValue(bucket, hi);
	} else {
		// Final bucket; throw in $ suffix
		appendValue(bucket, len);
	}
	VMSG_NL("Returning block of " << length(bucket));
	_cur++; // advance to next bucket
}

#endif /*BLOCKWISE_SA_H_*/
