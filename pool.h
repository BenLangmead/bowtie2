/*
 * Copyright 2011, Ben Langmead <blangmea@jhsph.edu>
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

#ifndef POOL_H_
#define POOL_H_

#include <iostream>
#include <stdexcept>
#include <string.h>
#include <stdlib.h>
#include "bitset.h"
#include "log.h"
#include "ds.h"

/**
 * Very simple allocator for fixed-size chunks of memory.  Chunk size
 * is set at construction time.  Heap memory is only allocated at
 * construction and deallocated at destruction.
 */
class ChunkPool {
public:
	/**
	 * Initialize a new pool with an initial size of about 'bytes'
	 * bytes.  Exit with an error message if we can't allocate it.
	 */
	ChunkPool(uint32_t chunkSz, uint32_t totSz, bool verbose_) :
		verbose(verbose_), patid(0), pool_(NULL), cur_(0),
		chunkSz_(chunkSz), totSz_(totSz), lim_(totSz/chunkSz),
		bits_(lim_), exhaustCrash_(false),
		lastSkippedRead_(0xffffffff), readName_(NULL)
	{
		assert_gt(lim_, 0);
		try {
			if((pool_ = new int8_t[totSz_]) == NULL) {
				throw std::bad_alloc();
			}
		} catch(std::bad_alloc& e) {
			std::cerr << "Error: Could not allocate ChunkPool of "
			          << totSz << " bytes" << std::endl;
			exhausted();
			throw 1; // Exit if we haven't already
		}
	}

	/**
	 * Delete all the pools.
	 */
	~ChunkPool() {
		if(pool_ != NULL) delete[] pool_;
	}

	/**
	 * Reset the pool, freeing all arrays that had been given out.
	 */
	void reset(const BTString* name, TReadId patid_) {
		patid = patid_;
		readName_ = name;
		cur_ = 0;
		bits_.clear();
		assert_eq(0, bits_.test(0));
	}

	/**
	 * Return our current position.
	 */
	uint32_t pos() {
		return cur_;
	}

	/**
	 * Return our current position.
	 */
	uint32_t remaining() {
		assert_geq(lim_, cur_);
		return lim_ - cur_;
	}

	/**
	 * Allocate a single T from the pool.
	 */
	void* alloc() {
		assert_lt(cur_, lim_);
		uint32_t cur = cur_;
		while(bits_.test(cur)) {
			cur++;
			if(cur >= lim_) {
				cur = 0;
			}
			if(cur == cur_) {
				// Wrapped all the way around without finding a free
				// chunk
				return NULL;
			}
		}
		void * ptr = (void *)(&pool_[cur * chunkSz_]);
		assert(!bits_.test(cur));
		bits_.set(cur);
		assert(bits_.test(cur));
		if(verbose) {
			stringstream ss;
			ss << patid << ": Allocating chunk with offset: " << cur;
			glog.msg(ss.str());
		}
		cur_ = cur;
		return ptr;
	}

	/**
	 *
	 */
	void free(void *ptr) {
		uint32_t off = (uint32_t)((int8_t*)ptr - pool_);
		assert_eq(0, off % chunkSz_);
		off /= chunkSz_;
		if(verbose) {
			stringstream ss;
			ss << patid << ": Freeing chunk with offset: " << cur_;
			glog.msg(ss.str());
		}
		bits_.clear(off);
	}

	/**
	 *
	 */
	uint32_t chunkSize() const {
		return chunkSz_;
	}

	/**
	 *
	 */
	uint32_t totalSize() const {
		return totSz_;
	}

	/**
	 * Utility function to call when memory has been exhausted.
	 * Currently just prints a friendly message and quits.
	 */
	void exhausted() {
		if(patid != lastSkippedRead_) {
			if(!exhaustCrash_) {
				std::cerr << "Warning: ";
			}
			std::cerr << "Exhausted best-first chunk memory for read " << (*readName_) << " (patid " << patid << "); skipping read" << std::endl;
		}
		if(exhaustCrash_) {
			std::cerr << "Please try specifying a larger --chunkmbs <int> (default is 32)" << std::endl;
			throw 1;
		}
		lastSkippedRead_ = patid;
	}

	bool verbose;
	TReadId patid;

protected:

	int8_t*  pool_; /// the memory pools
	uint32_t cur_;  /// index of next free element of pool_
	const uint32_t chunkSz_;
	const uint32_t totSz_;
	uint32_t lim_;  /// # elements held in pool_
	FixedBitset2 bits_;
	bool exhaustCrash_; /// abort hard when memory's exhausted?
	TReadId lastSkippedRead_;
	const BTString* readName_;
};

/**
 * Class for managing a pool of memory from which items of type T
 * (which must have a default constructor) are allocated.  Does not
 * support freeing or resizing individual items - just allocation and
 * then freeing all items at once.
 */
template<typename T>
class AllocOnlyPool {
	typedef std::pair<uint32_t, uint32_t> U32Pair;
public:
	/**
	 * Initialize a new pool with an initial size of about 'bytes'
	 * bytes.  Exit with an error message if we can't allocate it.
	 */
	AllocOnlyPool(ChunkPool* pool, const char *name) :
		pool_(pool), name_(name), curPool_(0), cur_(0)
	{
		assert(pool != NULL);
		lim_ = pool->chunkSize() / sizeof(T);
		assert_gt(lim_, 0);
		assert_gt(lim_, 1024);
	}

	/**
	 * Reset the pool, freeing all arrays that had been given out.
	 */
	void reset() {
		pools_.clear();
		lastCurInPool_.clear();
		cur_ = 0;
		curPool_ = 0;
	}

	/**
	 * Allocate a single T from the pool.
	 */
	T* alloc() {
		if(!lazyInit()) return NULL;
		if(cur_ + 1 >= lim_) {
			if(!allocNextPool()) return NULL;
		}
		cur_ ++;
		return &pools_[curPool_][cur_ - 1];
	}

	/**
	 * Allocate a single T from the pool and clear it.
	 */
	T* allocC() {
		T* t = alloc();
		if(t != NULL) {
			memset(t, 0, sizeof(T));
		}
		return t;
	}

	/**
	 * Allocate an array of Ts from the pool.
	 */
	T* alloc(uint32_t num) {
		if(!lazyInit()) return NULL;
		if(cur_ + num >= lim_) {
			if(!allocNextPool()) return NULL;
		}
		cur_ += num;
		return &pools_[curPool_][cur_ - num];
	}

	/**
	 * Allocate an array of Ts and clear them.
	 */
	T* allocC(uint32_t num) {
		T* t = alloc(num);
		if(t != NULL) {
			memset(t, 0, sizeof(T) * num);
		}
		return t;
	}

	/**
	 * Return the current pool.
	 */
	uint32_t curPool() const {
		return curPool_;
	}

	/**
	 * Return the current position within the current pool.
	 */
	uint32_t cur() const {
		return cur_;
	}

	/**
	 * Free a pointer allocated from this pool.  Fow, we only know how
	 * to free the topmost element.
	 */
	void free(T* t) {
		assert(t != NULL);
		if(pool_->verbose) {
			stringstream ss;
			ss << pool_->patid << ": Freeing a " << name_;
			glog.msg(ss.str());
		}
		if(cur_ > 0 && t == &pools_[curPool_][cur_-1]) {
			cur_--;
			ASSERT_ONLY(memset(&pools_[curPool_][cur_], 0, sizeof(T)));
			if(cur_ == 0 && curPool_ > 0) {
				rewindPool();
			}
		}
	}

	/**
	 * Free an array of pointers allocated from this pool.  For now, we
	 * only know how to free the topmost array.
	 */
	bool free(T* t, uint32_t num) {
		assert(t != NULL);
		if(pool_->verbose) {
			stringstream ss;
			ss << pool_->patid << ": Freeing " << num << " " << name_ << "s";
			glog.msg(ss.str());
		}
		if(num <= cur_ && t == &pools_[curPool_][cur_ - num]) {
			cur_ -= num;
			ASSERT_ONLY(memset(&pools_[curPool_][cur_], 0, num * sizeof(T)));
			if(cur_ == 0 && curPool_ > 0) {
				rewindPool();
			}
			return true; // deallocated
		}
		return false; // didn't deallocate
	}

	/**
	 * Return a unique (with respect to every other object allocated
	 * from this pool) identifier for the last object that was just
	 * allocated.
	 */
	uint32_t lastId() const {
		return (curPool_ << 16) | cur_;
	}

#ifndef NDEBUG
	bool empty() const {
		assert(pools_.empty());
		assert_eq(0, cur_);
		assert_eq(0, curPool_);
		return true;
	}
#endif

protected:

	bool allocNextPool() {
		assert_eq(curPool_+1, pools_.size());
		T *pool;
		try {
			if((pool = (T*)pool_->alloc()) == NULL) {
				throw std::bad_alloc();
			}
		} catch(std::bad_alloc& e) {
			//std::cerr << "Error: Could not allocate " << name_ << " pool #" << (curPool_+1) << " of " << (lim_ * sizeof(T)) << " bytes" << std::endl;
			pool_->exhausted();
			return false;
		}
		ASSERT_ONLY(memset(pool, 0, lim_ * sizeof(T)));
		pools_.push_back(pool);
		lastCurInPool_.push_back(cur_);
		curPool_++;
		cur_ = 0;
		return true;
	}

	bool lazyInit() {
		if(cur_ == 0 && pools_.empty()) {
			T *pool;
			try {
				if((pool = (T*)pool_->alloc()) == NULL) {
					throw std::bad_alloc();
				}
			} catch(std::bad_alloc& e) {
				//std::cerr << "Error: Could not allocate " << name_ << " pool #1" << std::endl;
				pool_->exhausted();
				return false;
			}
			ASSERT_ONLY(memset(pool, 0, lim_ * sizeof(T)));
			pools_.push_back(pool);
			assert_eq(1, pools_.size());
		}
		assert(!pools_.empty());
		return true;
	}

	void rewindPool() {
		assert_eq(curPool_+1, pools_.size());
		assert_eq(curPool_, lastCurInPool_.size());
		if(pool_->verbose) {
			stringstream ss;
			ss << pool_->patid << ": Freeing a " << name_ << " pool";
			glog.msg(ss.str());
		}
		pool_->free(pools_.back());
		pools_.pop_back();
		curPool_--;
		assert_gt(lastCurInPool_.size(), 0);
		cur_ = lastCurInPool_.back();
		lastCurInPool_.pop_back();
	}

	ChunkPool*      pool_;
	const char     *name_;
	EList<T*>      pools_; /// the memory pools
	uint32_t        curPool_; /// pool we're current allocating from
	EList<uint32_t> lastCurInPool_;
	uint32_t        lim_;  /// # elements held in pool_
	uint32_t        cur_;  /// index of next free element of pool_
};

#endif /* POOL_H_ */

