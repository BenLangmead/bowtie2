/**
 * seed_scan.h
 *
 * Somewhat generic routines for searching for any of a set of up to 64-bit
 * words in a longer string of 2-bit characters.  Uses a simple N-bitpair
 * filter to find candidates for extension to the left.  N might be 4 or 5; it
 * can't be so high that a lookup table with 2^N entries is impractical.
 * Scanner then checks the last up-to-32 characters seen to see if they match
 * any of the sequences in a list (list = binfull of queries with the seed
 * sequence).  There is one such list for every possible filter byte (for a
 * total of 256).
 *
 * 
 *
 *              -> Scan direction ->
 *
 *  CATTACGGCTAGCGTAGTGCTAGCTAGCTGCAGCATTATAAGCGGCTA
 *                                    |--|
 *                                    ATTA (rotating 1-byte buffer)
 *        |------------------------------|
 *        GGCTAGCGTAGTGCTAGCTAGCTGCAGCATTA
 *
 * When a hit is found 
 */

#ifndef SEED_SCAN_H_
#define SEED_SCAN_H_

#include <stdint.h>
#include <utility>
#include "ds.h"

#define BP_SUFFIX(seq, bps) \
	(seq & ~(0xfffffffffffffffflu << (bps << 1)))

#define BP_SUFFIX_EQ(seq1, seq2, bps) \
	(BP_SUFFIX(seq1, bps) == BP_SUFFIX(seq2, bps))

/**
 * Encapsulates something we're looking for (i.e. a uint64_t) along with an
 * identifier used to report at hit to the user.
 */
class SeedScanQuery {
public:

	typedef std::pair<uint32_t,uint32_t> U32Pair;

	SeedScanQuery(
		U32Pair id,
		uint64_t seq,
		size_t len)
	{
		init(id, seq, len);
	}

	SeedScanQuery() { reset(); }

	void reset() {
		id_ = std::make_pair(0, 0);
		seq_ = 0;
		len_ = 0;
		assert(!inited());
	}

	void init(U32Pair id, uint64_t seq, size_t len) {
		assert_gt(len, 0);
		id_ = id;
		seq_ = seq;
		len_ = len;
		assert(inited());
	}
	
	/**
	 * Return true iff this object has been initialized with a sequence and
	 * valid length.
	 */
	bool inited() const {
		return len_ != 0;
	}
	
	inline U32Pair  id()  const { return id_; }
	inline uint64_t seq() const { return seq_; }
	inline size_t   len() const { return len_; }

protected:
	U32Pair  id_;
	uint64_t seq_;
	size_t   len_;
};

/**
 * Encapsulates an occurrence of a query sequence in the reference stream.
 */
class SeedScanHit {
public:

	typedef std::pair<uint32_t,uint32_t> U32Pair;

	inline SeedScanHit() { reset(); }

	inline SeedScanHit(U32Pair id, uint32_t off) { init(id, off); }
	
	/**
	 * Set to uninitialized state.
	 */
	void reset() {
		id_ = std::make_pair(0xffffffff, 0);
		off_ = 0;
	}
	
	/**
	 * Return true iff hit is initialized.
	 */
	bool inited() const {
		return id_.first != 0xffffffff;
	}
	
	/**
	 * Initialize to given parameters.
	 */
	inline void init(U32Pair id, uint32_t off) {
		id_ = id;
		off_ = off;
	}

	inline U32Pair  id()  const { return id_; }
	inline uint32_t off() const { return off_; }

protected:
	U32Pair  id_;  // id of query sequence that hit
	uint32_t off_; // offset into the stream where it hit
};

/**
 * The table that (a) holds all the seeds, and (b), 
 */
class SeedScanTable {

public:

	typedef std::pair<uint32_t,uint32_t> U32Pair;

	SeedScanTable() : qrys_() {
		reset();
	}

	SeedScanTable(size_t keybps) : qrys_() {
		init(keybps);
	}
	
	/**
	 * Reset the table to be empty and uninitialized.
	 */
	void reset() {
		qrys_.clear();
		first_ = true;
		assert(empty());
		assert(!inited());
	}
	
	/**
	 * Initialize with a potentially new key length.  This determines the
	 * number of bins.  Resize qrys_ to the appropriate number of bins.
	 */
	void init(size_t keybps) {
		keybps_ = keybps;
		qrys_.resize(1 << (keybps_ << 1));
		for(size_t i = 0; i < qrys_.size(); i++) {
			qrys_[i].clear();
		}
		first_ = true;
		assert(inited());
	}
	
	/**
	 * Return true iff this table has been initialized and bins are ready to
	 * accept queries.
	 */
	bool inited() const {
		return !qrys_.empty();
	}
	
	/**
	 * Check that table is internally consistent.
	 */
	bool repOk() const {
		assert(!inited() || qrys_.size() == (1ul << (keybps_ << 1ul)));
		return true;
	}
	
	/**
	 * Add a new query and place it in the appropriate bin.
	 */
	void add(U32Pair id, uint64_t seq, size_t len) {
		assert(empty() || inited());
		add(SeedScanQuery(id, seq, len));
	}

	/**
	 * Add a new query and place it in the appropriate bin.
	 */
	void add(const SeedScanQuery& qry) {
		assert(empty() || inited());
		if(first_) {
			len_ = qry.len();
			first_ = false;
		} else {
			assert_eq(len_, qry.len());
		}
		size_t key = (size_t)BP_SUFFIX(qry.seq(), keybps_);
		// TODO: check if a query has been added twice?
		qrys_[key].push_back(qry);
	}
	
	/**
	 * Check if the current buffer contains one or more query strings.  If so,
	 * add each as a hit to the 'hits' list.
	 */
	inline void query(
		uint64_t buf,
		uint32_t off,
		EList<SeedScanHit>& hits) const
	{
		size_t key = (size_t)BP_SUFFIX(buf, keybps_);
		const size_t binsz = qrys_[key].size();
		if(binsz == 0) {
			return;
		}
		buf = BP_SUFFIX(buf, len());
		for(size_t i = 0; i < binsz; i++) {
			if(buf == qrys_[key][i].seq()) {
				// Hit!
				hits.push_back(SeedScanHit(qrys_[key][i].id(), off));
			}
		}
	}
	
	inline size_t len()   const { return len_; }
	inline bool   empty() const { return first_; }
	
protected:

	bool first_;                 // true -> at least query already added
	size_t len_;                 // length of all queries
	size_t keybps_;              // number of bitpairs to use as bin id
	ELList<SeedScanQuery> qrys_; // queries divided into bin by key
	
};

/**
 * User provides scanner with successive characters from a text; scanner
 * accumualtes a list of hits for its query sequences, already installed in the
 * SeedScanTable.
 */
class SeedScanner {
public:

	typedef std::pair<uint32_t,uint32_t> U32Pair;

	SeedScanner() { }

	/**
	 * Give the next consecutive character to the scanner.  If the character
	 * triggers any hits, they are stored in the hits_ list.
	 */
	inline void nextChar(int c) {
		assert(tab_ != NULL);
		if(tab_->empty()) {
			return;
		}
		// Add it to the buffer
		assert_range(0, 4, c);
		if(c == 4) {
			lastUnmatchable_ = 0;
			c = 0;
		}
		assert_range(0, 3, c);
		buf_ <<= 2;
		buf_ |= c;
		if(lastUnmatchable_ >= tab_->len()) {
			tab_->query(buf_, off_, hits_);
		}
		lastUnmatchable_++;
		off_++;
	}
	
	/**
	 * Initialize with new table in preparation for new character stream.
	 */
	void init(const SeedScanTable& tab) {
		tab_ = &tab;
		buf_ = off_ = 0;
		hits_.clear();
		lastUnmatchable_ = 1;
	}
	
	/**
	 * Return true iff this scanner is currently initialized with a problem.
	 */
	bool inited() const {
		return tab_ != NULL;
	}
	
	/**
	 * Return true iff there are no queries to scan for.
	 */
	bool empty() const {
		return tab_ == NULL || tab_->empty();
	}
	
	/**
	 * Return the list of hits so far.
	 */
	const EList<SeedScanHit>& hits() const {
		return hits_;
	}

protected:
	uint64_t buf_;             // rotating buffer of recent characters
	uint32_t off_;             // offset w/r/t last time init() was called
	EList<SeedScanHit> hits_;  // hits so far
	const SeedScanTable *tab_; // table with queries
	size_t lastUnmatchable_;   // how many chars ago was last N?
};

#endif /*ndef SEED_SCAN_H_*/
