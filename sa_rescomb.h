/**
 * sa_rescomb.h
 *
 * Routines for combining results from a Burrows-Wheeler-based offset resolver,
 * where each resolution is associated with a row in the Burrows-Wheeler matrix,
 * and a reference-scanning offset resolver, where resolutions are not
 * associated with any particular row in the Burrows-Wheeler matrix.  The idea
 * is that even though the latter result cannot be associated with a BW row per
 * se, we might immediately or eventually be able to associated it with a row
 * by combining results from both resolves and applying process of elimination.
 *
 * A trivial but common example is if the reference-scanning resolver resolves
 * an offset for a range that has just one element.  In which case, we know
 * that the offset must be associated with that (only) element.
 *
 * When a new offset can be resolved by combining information from the two
 * resolvers, it is set directly in the cache's salist.
 */

#ifndef SA_RESCOMB_H_
#define SA_RESCOMB_H_

#include "ds.h"
#include "aligner_cache.h"

/**
 * For a given seed-hit range, keep track of resolutions in the salist as well
 * as resolutions ascertained by reference scanning.  When it's possible to
 * resolve a new offset in the salist, do so.
 */
class SAResolveCombiner {

public:

	SAResolveCombiner() { reset(); }

	SAResolveCombiner(SATuple& tup) { init(tup); }
	
	/**
	 * Reset, leaving combiner uninitialized.
	 */
	void reset() {
		inited_ = false;
	}
	
	/**
	 * Initialize the combiner with a new range in the salist.
	 */
	void init(SATuple& tup) {
		satup_ = &tup;
		refscan_.clear();
		inited_ = true;
	}
	
	/**
	 * Return true iff combiner has been initialized with a range in the salist.
	 */
	bool inited() const { return inited_; }
	
	/**
	 * Given the contents of satup_ and refscan_, see if any additional
	 * elements of satup_.offs() can be resolved.  Returns true iff all
	 * remaining elements of satup_.offs() are now filled (or were filled to
	 * begin with).
	 */
	bool tryResolving(uint64_t& nresolved);
	
	/**
	 * Try adding a new resolved offset for the sa range resolved using
	 * reference scanning.  If it's the same as a previously resolved offset,
	 * ignore it and return false.  Otherwise add it and return true.
	 */
	bool addRefscan(uint32_t off);
	
	/**
	 * The caller just resolved a new offset in the salist; update the salist
	 * and then call and return the result of tryResolving().
	 */
	bool setSalist(size_t saoff, uint32_t off, uint64_t& nresolved) {
		assert(inited_);
		assert(satup_ != NULL);
		satup_->offs[saoff] = off;
		return tryResolving(nresolved);
	}
	
	/**
	 * Return the SATuple.
	 */
	const SATuple& satup() const {
		return *satup_;
	}

protected:

	bool            inited_;  // true iff init() was called since last reset()
	EList<uint32_t> refscan_; // resolved offsets deduced from ref scanning
	SATuple*        satup_;   // tuple describing BW range
	EList<bool>     found_;   // indicators for which
};

#endif /*SA_RESCOMB_H_*/
