/*
 * range_cache.h
 *
 * Classes that encapsulate the caching of
 */

#ifndef RANGE_CACHE_H_
#define RANGE_CACHE_H_

#include <stdint.h>
#include <utility>
#include <iostream>
#include <stdexcept>
#include <map>
#include "ebwt.h"
#include "row_chaser.h"

#define RANGE_NOT_SET 0xffffffff
#define RANGE_CACHE_BAD_ALLOC 0xffffffff

/**
 * Manages a pool of memory used exclusively for range cache entries.
 * This manager is allocate-only; it exists mainly so that we can avoid
 * lots of new[]s and delete[]s.
 *
 * A given stretch of words may be one of two types: a cache entry, or
 * a cache entry wrapper.  A cache entry has a length and a list of
 * already-resolved reference positions.  A cache entry wrapper has a
 * pointer to a cache entry for a different range, along with an
 * integer indicating how many "jumps" to the left that range is from
 * the one that owns the wrapper.
 */
class RangeCacheMemPool {
public:
	RangeCacheMemPool(uint32_t lim /* max cache size in bytes */) :
		lim_(lim >> 2 /* convert to words */), occ_(0), buf_(NULL),
		closed_(false)
	{
		if(lim_ > 0) {
			try {
				buf_ = new uint32_t[lim_];
				if(buf_ == NULL) throw std::bad_alloc();
			} catch(std::bad_alloc& e) {
				cerr << "Allocation error allocating " << lim
					 << " words of range-cache memory" << endl;
				throw 1;
			}
			assert(buf_ != NULL);
			// Fill with 1s to signal that these elements are
			// uninitialized
			memset(buf_, 0xff, lim_ << 2 /* convert back to bytes */);
		}
	}

	~RangeCacheMemPool() {
		// Release all word memory!
		if(lim_ > 0) delete[] buf_;
	}

	/**
	 * Allocate numElts elements from the word pool.
	 */
	uint32_t alloc(uint32_t numElts) {
		assert_gt(numElts, 0);
		assert_leq(occ_, lim_);
		if(occ_ + numElts > lim_ || numElts >= 0x80000000) {
			return RANGE_CACHE_BAD_ALLOC;
		}
		assert_gt(lim_, 0);
		uint32_t ret = occ_;
		assert(allocs_.find(ret) == allocs_.end());
		ASSERT_ONLY(allocs_.insert(ret));
		// Clear the first elt so that we don't think there's already
		// something there
#ifndef NDEBUG
		for(size_t i = 0; i < numElts; i++) {
			assert_eq(0xffffffff, buf_[occ_ + i]);
		}
#endif
		buf_[occ_] = 0;
		occ_ += numElts;
		assert_leq(occ_, lim_);
		if(lim_ - occ_ < 10) {
			// No more room - don't try anymore
			closed_ = true;
		}
		return ret;
	}

	/**
	 * Turn a pool-array index into a pointer; check that it doesn't
	 * fall outside the pool first.
	 */
	inline uint32_t *get(uint32_t off) {
		assert_gt(lim_, 0);
		assert_lt(off, lim_);
		assert(allocs_.find(off) != allocs_.end());
		uint32_t *ret = buf_ + off;
		assert_neq(0x80000000, ret[0]);
		assert_neq(0xffffffff, ret[0]);
		return ret;
	}

	/**
	 * Return true iff there's no more room in the cache.
	 */
	inline bool closed() {
		return closed_;
	}

private:
	uint32_t lim_;  /// limit on number of 32-bit words to dish out in total
	uint32_t occ_;  /// number of occupied words
	uint32_t *buf_; /// buffer of 32-bit words
	bool closed_;   ///
#ifndef NDEBUG
	std::set<uint32_t> allocs_; // elements allocated
#endif
};

/**
 * A view to a range of cached reference positions.
 */
class RangeCacheEntry {

	typedef std::pair<uint32_t,uint32_t> U32Pair;

public:
	/**
	 *
	 */
	RangeCacheEntry(bool sanity = false) :
		top_(0xffffffff), jumps_(0), len_(0), ents_(NULL), ebwt_(NULL),
		sanity_(sanity)
	{ }

	/**
	 * Create a new RangeCacheEntry from the data in the pool at 'ents'.
	 */
	RangeCacheEntry(RangeCacheMemPool& pool, uint32_t top,
	                uint32_t ent, Ebwt* ebwt, bool sanity = false) :
	    sanity_(sanity)
	{
		init(pool, top, ent, ebwt);
	}

	/**
	 * Initialize a RangeCacheEntry from the data in the pool at 'ents'.
	 */
	void init(RangeCacheMemPool& pool, uint32_t top, uint32_t ent, Ebwt* ebwt) {
		assert(ebwt != NULL);
		top_ = top;
		ebwt_ = ebwt;
		uint32_t *ents = pool.get(ent);
		assert_neq(0x80000000, ents[0]);
		// Is hi bit set?
		if((ents[0] & 0x80000000) != 0) {
			// If so, the target is a wrapper and the non-hi bits
			// contain the # jumps
			jumps_ = (ents[0] & ~0x80000000);
			assert_gt(jumps_, 0);
			assert_leq(jumps_, ebwt_->_eh._len);
			// Get the target entry
			uint32_t *dest = pool.get(ents[1]);
			// Get the length from the target entry
			len_ = dest[0];
			assert_leq(top_ + len_, ebwt_->_eh._len);
			assert_gt(len_, 0);
			assert_leq(len_, ebwt_->_eh._len);
			// Get the pointer to the entries themselves
			ents_ = dest + 1;
		} else {
			// Not a wrapper, so there are no jumps
			jumps_ = 0;
			// Get the length from the target entry
			len_  = ents[0];
			assert_leq(top_ + len_, ebwt_->_eh._len);
			assert_gt(len_, 0);
			assert_leq(len_, ebwt_->_eh._len);
			// Get the pointer to the entries themselves
			ents_ = ents + 1;
		}
		assert(sanityCheckEnts());
	}

	/**
	 * Initialize a wrapper with given number of jumps and given target
	 * entry index.
	 */
	void init(RangeCacheMemPool& pool, uint32_t top, uint32_t jumps,
	          uint32_t ent, Ebwt* ebwt)
	{
		assert(ebwt != NULL);
		ebwt_ = ebwt;
		top_ = top;
		jumps_ = jumps;
		uint32_t *ents = pool.get(ent);
		// Must not be a wrapper
		assert_eq(0, ents[0] & 0x80000000);
		// Get the length from the target entry
		len_ = ents[0];
		assert_gt(len_, 0);
		assert_leq(len_, ebwt_->_eh._len);
		// Get the pointer to the entries themselves
		ents_ = ents + 1;
		assert_leq(top_ + len_, ebwt_->_eh._len);
		assert(sanityCheckEnts());
	}

	uint32_t len() const   {
		assert(ents_ != NULL);
		assert(ebwt_ != NULL);
		return len_;
	}

	uint32_t jumps() const {
		assert(ents_ != NULL);
		assert(ebwt_ != NULL);
		return jumps_;
	}

	/**
	 *
	 */
	void reset() {
		ents_ = NULL;
	}

	/**
	 * Return true iff this object represents a valid cache entry.
	 */
	bool valid() const {
		return ents_ != NULL;
	}

	Ebwt *ebwt() {
		return ebwt_;
	}

	/**
	 * Install a result obtained by a client of this cache; be sure to
	 * adjust for how many jumps down the tunnel the cache entry is
	 * situated.
	 */
	void install(uint32_t elt, uint32_t val) {
		if(ents_ == NULL) {
			// This is not a valid cache entry; do nothing
			return;
		}
		assert(ents_ != NULL);
		assert(ebwt_ != NULL);
		assert_leq(jumps_, val);
		assert_neq(0xffffffff, val);
		assert_leq(top_ + len_, ebwt_->_eh._len);
		if(elt < len_) {
			val -= jumps_;
			if(verbose_) cout << "Installed reference offset: " << (top_ + elt) << endl;
			ASSERT_ONLY(uint32_t sanity = RowChaser::toFlatRefOff(ebwt_, 1, top_ + elt));
			assert_eq(sanity, val);
#ifndef NDEBUG
			for(size_t i = 0; i < len_; i++) {
				if(i == elt) continue;
				assert_neq(val, ents_[i]);
			}
#endif
			ents_[elt] = val;
		} else {
			// ignore install request
			if(verbose_) cout << "Fell off end of cache entry for install: " << (top_ + elt) << endl;
		}
	}

	/**
	 * Get an element from the cache, adjusted for tunnel jumps.
	 */
	inline uint32_t get(uint32_t elt) const {
		if(ents_ == NULL) {
			// This is not a valid cache entry; do nothing
			return RANGE_NOT_SET;
		}
		assert(ents_ != NULL);
		assert(ebwt_ != NULL);
		assert_leq(top_ + len_, ebwt_->_eh._len);
		if(elt < len_ && ents_[elt] != RANGE_NOT_SET) {
			if(verbose_) cout << "Retrieved result from cache: " << (top_ + elt) << endl;
			uint32_t ret = ents_[elt] + jumps_;
			ASSERT_ONLY(uint32_t sanity = RowChaser::toFlatRefOff(ebwt_, 1, top_ + elt));
			assert_eq(sanity, ret);
			return ret;
		} else {
			if(verbose_) cout << "Cache entry not set: " << (top_ + elt) << endl;
			return RANGE_NOT_SET;
		}
	}

	/**
	 * Check that len_ and the ents_ array both make sense.
	 */
	static bool sanityCheckEnts(uint32_t len, uint32_t *ents, Ebwt* ebwt) {
		assert_gt(len, 0);
		assert_leq(len, ebwt->_eh._len);
		if(len < 10) {
			for(size_t i = 0; i < len; i++) {
				if(ents[i] == 0xffffffff) continue;
				assert_leq(ents[i], ebwt->_eh._len);
				for(size_t j = i+1; j < len; j++) {
					if(ents[j] == 0xffffffff) continue;
					assert_neq(ents[i], ents[j]);
				}
			}
		} else {
			std::set<uint32_t> seen;
			for(size_t i = 0; i < len; i++) {
				if(ents[i] == 0xffffffff) continue;
				assert(seen.find(ents[i]) == seen.end());
				seen.insert(ents[i]);
			}
		}
		return true;
	}

	/**
	 * Check that len_ and the ents_ array both make sense.
	 */
	bool sanityCheckEnts() {
		return RangeCacheEntry::sanityCheckEnts(len_, ents_, ebwt_);
	}

private:

	uint32_t top_;   /// top pointer for this range
	uint32_t jumps_; /// how many tunnel-jumps it is away from the requester
	uint32_t len_;   /// # of entries in cache entry
	uint32_t *ents_; /// ptr to entries, which are flat offs within joined ref
	//U32Pair *ents_;  /// pointer to entries, which are tidx,toff pairs
	Ebwt    *ebwt_; /// index that alignments are in
	bool     verbose_; /// be talkative?
	bool     sanity_;  /// do consistency checks?
};

/**
 *
 */
class RangeCache {

	typedef EList<uint32_t> TU32Vec;
	typedef std::map<uint32_t, uint32_t> TMap;
	typedef std::map<uint32_t, uint32_t>::iterator TMapItr;

public:
	RangeCache(uint32_t lim, Ebwt* ebwt) :
		lim_(lim), map_(), pool_(lim), closed_(false), ebwt_(ebwt), sanity_(true) { }

	/**
	 * Given top and bot offsets, retrieve the canonical cache entry
	 * that best covers that range.  The cache entry may not directly
	 * correspond to the top offset provided, rather, it might be an
	 * entry that lies "at the end of the tunnel" when top and bot are
	 * walked backward.
	 */
	bool lookup(uint32_t top, uint32_t bot, RangeCacheEntry& ent) {
		if(ebwt_ == NULL || lim_ == 0) return false;
		assert_gt(bot, top);
		ent.reset();
		TMapItr itr = map_.find(top);
		if(itr == map_.end()) {
			// No cache entry for the given 'top' offset
			if(closed_) {
				return false; // failed to get cache entry
			} else {
				if(pool_.closed()) {
					closed_ = true;
					return false; // failed to get cache entry
				}
			}
			// Use the tunnel
			bool ret = tunnel(top, bot, ent);
			return ret;
		} else {
			// There is a cache entry for the given 'top' offset
			uint32_t ret = itr->second;
			ent.init(pool_, top, ret, ebwt_);
			return true; // success
		}
	}

	/**
	 * Exhaustively check all entries linked to from map_ to ensure
	 * they're well-formed.
	 */
	bool repOk() {
#ifndef NDEBUG
		for(TMapItr itr = map_.begin(); itr != map_.end(); itr++) {
			uint32_t top = itr->first;
			uint32_t idx = itr->second;
			uint32_t jumps = 0;
			assert_leq(top, ebwt_->_eh._len);
			uint32_t *ents = pool_.get(idx);
			if((ents[0] & 0x80000000) != 0) {
				jumps = ents[0] & ~0x80000000;
				assert_leq(jumps, ebwt_->_eh._len);
				idx = ents[1];
				ents = pool_.get(idx);
			}
			uint32_t len = ents[0];
			assert_leq(top + len, ebwt_->_eh._len);
			RangeCacheEntry::sanityCheckEnts(len, ents + 1, ebwt_);
		}
#endif
		return true;
	}

protected:

	/**
	 * Tunnel through to the first range that 1) includes all the same
	 * suffixes (though longer) as the given range, and 2) has a cache
	 * entry for it.
	 */
	bool tunnel(uint32_t top, uint32_t bot, RangeCacheEntry& ent) {
		assert_gt(bot, top);
		tops_.clear();
		const uint32_t spread = bot - top;
		SideLocus tloc, bloc;
		SideLocus::initFromTopBot(top, bot, ebwt_->_eh, ebwt_->ebwt(), tloc, bloc);
		uint32_t newtop = top, newbot = bot;
		uint32_t jumps = 0;
		// Walk left through the tunnel
		while(true) {
			if(ebwt_->rowL(tloc) != ebwt_->rowL(bloc)) {
				// Different characters at top and bot positions of
				// BWT; this means that the calls to mapLF below are
				// guaranteed to yield rows in two different character-
				// sections of the BWT.
				break;
			}
			// Advance top and bot
			newtop = ebwt_->mapLF(tloc);
			newbot = ebwt_->mapLF(bloc);
			assert_geq(newbot, newtop);
			assert_leq(newbot - newtop, spread);
			// If the new spread is the same as the old spread, we can
			// be confident that the new range includes all of the same
			// suffixes as the last range (though longer by 1 char)
			if((newbot - newtop) == spread) {
				// Check if newtop is already cached
				TMapItr itr = map_.find(newtop);
				jumps++;
				if(itr != map_.end()) {
					// This range, which is further to the left in the
					// same tunnel as the query range, has a cache
					// entry already, so use that
					uint32_t idx = itr->second;
					uint32_t *ents = pool_.get(idx);
					if((ents[0] & 0x80000000) != 0) {
						// The cache entry we found was a wrapper; make
						// a new wrapper that points to that wrapper's
						// target, with the appropriate number of jumps
						jumps += (ents[0] & ~0x80000000);
						idx = ents[1];
					}
					// Allocate a new wrapper
					uint32_t newentIdx = pool_.alloc(2);
					if(newentIdx != RANGE_CACHE_BAD_ALLOC) {
						// We successfully made a new wrapper entry;
						// now populate it and install it in map_
						uint32_t *newent = pool_.get(newentIdx); // get ptr to it
						assert_eq(0, newent[0]);
						newent[0] = 0x80000000 | jumps; // set jumps
						newent[1] = idx;                // set target
						assert(map_.find(top) == map_.end());
						map_[top] = newentIdx;
						if(sanity_) assert(repOk());
					}
					// Initialize the entry
					ent.init(pool_, top, jumps, idx, ebwt_);
					return true;
				}
				// Save this range
				tops_.push_back(newtop);
				SideLocus::initFromTopBot(newtop, newbot, ebwt_->_eh, ebwt_->ebwt(), tloc, bloc);
			} else {
				// Not all the suffixes were preserved, so we can't
				// link the source range's cached result to this
				// range's cached results
				break;
			}
			assert_eq(jumps, tops_.size());
		}
		assert_eq(jumps, tops_.size());
		// Try to create a new cache entry for the leftmost range in
		// the tunnel (which might be the query range)
		uint32_t newentIdx = pool_.alloc(spread + 1);
		if(newentIdx != RANGE_CACHE_BAD_ALLOC) {
			// Successfully allocated new range cache entry; install it
			uint32_t *newent = pool_.get(newentIdx);
			assert_eq(0, newent[0]);
			// Store cache-range length in first word
			newent[0] = spread;
			assert_lt(newent[0], 0x80000000);
			assert_eq(spread, newent[0]);
			uint32_t entTop = top;
			uint32_t jumps = 0;
			if(tops_.size() > 0) {
				entTop = tops_.back();
				jumps = tops_.size();
			}
			// Cache the entry for the end of the tunnel
			assert(map_.find(entTop) == map_.end());
			map_[entTop] = newentIdx;
			if(sanity_) assert(repOk());
			ent.init(pool_, entTop, jumps, newentIdx, ebwt_);
			assert_eq(spread, newent[0]);
			if(jumps > 0) {
				assert_neq(entTop, top);
				// Cache a wrapper entry for the query range (if possible)
				uint32_t wrapentIdx = pool_.alloc(2);
				if(wrapentIdx != RANGE_CACHE_BAD_ALLOC) {
					uint32_t *wrapent = pool_.get(wrapentIdx);
					assert_eq(0, wrapent[0]);
					wrapent[0] = 0x80000000 | jumps;
					wrapent[1] = newentIdx;
					assert(map_.find(top) == map_.end());
					map_[top] = wrapentIdx;
					if(sanity_) assert(repOk());
				}
			}
			return true;
		} else {
			// Could not allocate new range cache entry
			return false;
		}
	}

	uint32_t lim_;           /// Total number of key/val bytes to keep in cache
	TMap map_;               ///
	RangeCacheMemPool pool_; /// Memory pool
	bool closed_;            /// Out of space; no new entries
	Ebwt* ebwt_;             /// Index that alignments are in
	bool sanity_;
	TU32Vec tops_;
};

#endif /* RANGE_CACHE_H_ */
