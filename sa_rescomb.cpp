/**
 * sa_rescomb.cpp
 */

#include "sa_rescomb.h"

using namespace std;

/**
 * Given the contents of satup_ and refscan_, see if any additional
 * elements of satup_->offs can be resolved.  Returns true iff all
 * remaining elements of satup_->offs are now filled (or were filled to
 * begin with).
 */
bool SAResolveCombiner::tryResolving() {
	if(satup_->offs.size() == 1) {
		// Just one element in range; easy to check if it's been resolved
		if(satup_->offs[0] != 0xffffffff) {
			return true; // already filled in
		}
		if(refscan_.empty()) {
			return false; // not filled in, nothing to add
		}
		assert_eq(1, refscan_.size());
		assert_neq(0xffffffff, refscan_[0]);
		// Resolved!
		satup_->offs[0] = refscan_[0];
		return true; // filled everything (1 thing) in
	} else {
		// More than one element
		// Count how many elements are yet to be resolved in satup_
		size_t needResolving = 0;
		for(size_t i = 0; i < satup_->offs.size(); i++) {
			if(satup_->offs[i] == 0xffffffff) {
				needResolving++;
			}
		}
		if(needResolving == 0) {
			return true; // already filled in
		}
		// Count how many elements are in refscan_ but not in satup_
		found_.resize(refscan_.size());
		size_t nfound = 0;
		for(size_t i = 0; i < refscan_.size(); i++) {
			found_[i] = true;
			nfound++;
			for(size_t j = 0; j < satup_->offs.size(); j++) {
				if(satup_->offs[j] == refscan_[i]) {
					assert(found_[i]);
					found_[i] = false;
					nfound--;
					break;
				}
			}
		}
		assert_leq(nfound, needResolving);
		if(nfound == needResolving) {
			// We can fill in the as-yet-unresolved elements of the salist
			size_t nexti = 0;
			ASSERT_ONLY(size_t added = 0);
			for(size_t i = 0; i < refscan_.size(); i++) {
				if(found_[i]) {
					// Put it in the next open salist slot
					assert_neq(0xffffffff, refscan_[i]);
					while(satup_->offs[nexti] != 0xffffffff) {
						nexti++;
						assert_lt(nexti, satup_->offs.size());
					}
					satup_->offs[nexti] = refscan_[i];
					ASSERT_ONLY(added++);
					nexti++;
					assert_leq(nexti, satup_->offs.size());
				}
			}
			assert_eq(nfound, added);
			return true; // filled everything in
		}
		return false;
	}
}

/**
 * Try adding a new resolved offset for the sa range resolved using
 * reference scanning.  If it's the same as a previously resolved offset,
 * ignore it and return false.  Otherwise add it and return true.
 */
bool SAResolveCombiner::addRefscan(uint32_t off) {
	for(size_t i = 0; i < refscan_.size(); i++) {
		if(refscan_[i] == off) {
			return false;
		}
	}
	refscan_.push_back(off);
	return true;
}

#ifdef SA_RESCOMB_MAIN

#include <iostream>

int main(void) {
	cerr << "Case 1 (1 element already resolved) ... ";
	{
		SAResolveCombiner comb;
		SAKey k;
		TSlice s;
		Pool pool(10 * 1024, 1024);
		PList<uint32_t, CACHE_PAGE_SZ> list;
		for(size_t i = 0; i < 1; i++) {
			list.add(pool, 0); // already resolved
		}
		assert_eq(1, list.size());
		PListSlice<uint32_t, CACHE_PAGE_SZ> slice(list, 0, 1);
		SATuple sat(k, 0, slice);
		comb.init(sat);
		assert(comb.tryResolving());
		for(size_t i = 0; i < 1; i++) {
			assert_eq(i, list.get(i));
		}
	}
	cerr << "PASSED" << endl;

	cerr << "Case 2 (3 elements already resolved) ... ";
	{
		SAResolveCombiner comb;
		SAKey k;
		TSlice s;
		Pool pool(10 * 1024, 1024);
		PList<uint32_t, CACHE_PAGE_SZ> list;
		for(size_t i = 0; i < 3; i++) {
			list.add(pool, (uint32_t)i); // already resolved
		}
		assert_eq(3, list.size());
		PListSlice<uint32_t, CACHE_PAGE_SZ> slice(list, 0, 3);
		SATuple sat(k, 0, slice);
		comb.init(sat);
		assert(comb.tryResolving());
		for(size_t i = 0; i < 3; i++) {
			assert_eq(i, list.get(i));
		}
	}
	cerr << "PASSED" << endl;

	cerr << "Case 3 (2 elements already resolved, 1 then resolved) ... ";
	{
		SAResolveCombiner comb;
		SAKey k;
		TSlice s;
		Pool pool(10 * 1024, 1024);
		PList<uint32_t, CACHE_PAGE_SZ> list;
		for(size_t i = 0; i < 3; i++) {
			if(i < 2) {
				list.add(pool, (uint32_t)i); // already resolved
			} else {
				list.add(pool, 0xffffffff); // already resolved
			}
		}
		assert_eq(3, list.size());
		PListSlice<uint32_t, CACHE_PAGE_SZ> slice(list, 0, 3);
		SATuple sat(k, 0, slice);
		comb.init(sat);
		assert(!comb.tryResolving());
		comb.addRefscan(2);
		assert(comb.tryResolving());
		for(size_t i = 0; i < 3; i++) {
			assert_eq(i, list.get(i));
		}
	}
	cerr << "PASSED" << endl;

	cerr << "Case 4 (2 elts pre-resolved, 1 resolved after redundants) ... ";
	{
		SAResolveCombiner comb;
		SAKey k;
		TSlice s;
		Pool pool(10 * 1024, 1024);
		PList<uint32_t, CACHE_PAGE_SZ> list;
		for(size_t i = 0; i < 3; i++) {
			if(i < 2) {
				list.add(pool, (uint32_t)i); // already resolved
			} else {
				list.add(pool, 0xffffffff); // already resolved
			}
		}
		assert_eq(3, list.size());
		PListSlice<uint32_t, CACHE_PAGE_SZ> slice(list, 0, 3);
		SATuple sat(k, 0, slice);
		comb.init(sat);
		assert(!comb.tryResolving());
		comb.addRefscan(1); // already there
		assert(!comb.tryResolving());
		comb.addRefscan(0); // already there
		assert(!comb.tryResolving());
		comb.addRefscan(2);
		assert(comb.tryResolving());
		for(size_t i = 0; i < 3; i++) {
			assert_eq(i, list.get(i));
		}
	}
	cerr << "PASSED" << endl;

	cerr << "Case 5 (2 elts pre-resolved, 1 via salist after redundants) ... ";
	{
		SAResolveCombiner comb;
		SAKey k;
		TSlice s;
		Pool pool(10 * 1024, 1024);
		PList<uint32_t, CACHE_PAGE_SZ> list;
		for(size_t i = 0; i < 3; i++) {
			if(i < 2) {
				list.add(pool, (uint32_t)i); // already resolved
			} else {
				list.add(pool, 0xffffffff); // already resolved
			}
		}
		assert_eq(3, list.size());
		PListSlice<uint32_t, CACHE_PAGE_SZ> slice(list, 0, 3);
		SATuple sat(k, 0, slice);
		comb.init(sat);
		assert(!comb.tryResolving());
		comb.addRefscan(1); // already there
		assert(!comb.tryResolving());
		comb.addRefscan(0); // already there
		assert(!comb.tryResolving());
		assert(comb.setSalist(2, 2));
		assert(comb.tryResolving());
		for(size_t i = 0; i < 3; i++) {
			assert_eq(i, list.get(i));
		}
	}
	cerr << "PASSED" << endl;

	cerr << "Case 6 (10 elements, resolved via refscan) ... ";
	{
		SAResolveCombiner comb;
		SAKey k;
		TSlice s;
		Pool pool(10 * 1024, 1024);
		PList<uint32_t, CACHE_PAGE_SZ> list;
		for(size_t i = 0; i < 20; i++) {
			list.add(pool, 0xffffffff);
		}
		assert_eq(20, list.size());
		PListSlice<uint32_t, CACHE_PAGE_SZ> slice(list, 1, 10);
		SATuple sat(k, 0, slice);
		comb.init(sat);
		comb.addRefscan(0);
		assert(!comb.tryResolving());
		comb.addRefscan(1);
		assert(!comb.tryResolving());
		comb.addRefscan(2);
		assert(!comb.tryResolving());
		comb.addRefscan(3);
		assert(!comb.tryResolving());
		comb.addRefscan(4);
		assert(!comb.tryResolving());
		comb.addRefscan(5);
		assert(!comb.tryResolving());
		comb.addRefscan(6);
		assert(!comb.tryResolving());
		comb.addRefscan(7);
		assert(!comb.tryResolving());
		comb.addRefscan(8);
		assert(!comb.tryResolving());
		comb.addRefscan(9);
		assert(comb.tryResolving());
		for(size_t i = 0; i < 20; i++) {
			if(i >= 1 && i < 11) {
				assert_eq(i-1, list.get(i));
			} else {
				assert_eq(0xffffffff, list.get(i));
			}
		}
	}
	cerr << "PASSED" << endl;

	cerr << "Case 7 (10 elements, resolved via salist) ... ";
	{
		SAResolveCombiner comb;
		SAKey k;
		TSlice s;
		Pool pool(10 * 1024, 1024);
		PList<uint32_t, CACHE_PAGE_SZ> list;
		for(size_t i = 0; i < 20; i++) {
			list.add(pool, 0xffffffff);
		}
		assert_eq(20, list.size());
		PListSlice<uint32_t, CACHE_PAGE_SZ> slice(list, 1, 10);
		SATuple sat(k, 0, slice);
		comb.init(sat);
		assert(!comb.setSalist(0, 0));
		assert(!comb.setSalist(1, 1));
		assert(!comb.setSalist(2, 2));
		assert(!comb.setSalist(3, 3));
		assert(!comb.setSalist(4, 4));
		assert(!comb.setSalist(5, 5));
		assert(!comb.setSalist(6, 6));
		assert(!comb.setSalist(7, 7));
		assert(!comb.setSalist(8, 8));
		assert(comb.setSalist(9, 9));
		for(size_t i = 0; i < 20; i++) {
			if(i >= 1 && i < 11) {
				assert_eq(i-1, list.get(i));
			} else {
				assert_eq(0xffffffff, list.get(i));
			}
		}
	}
	cerr << "PASSED" << endl;
}

#endif /*def SA_RESCOMB_MAIN*/

