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

#ifndef BINARY_SA_SEARCH_H_
#define BINARY_SA_SEARCH_H_

#include <stdint.h>
#include <iostream>
#include <limits>
#include "alphabet.h"
#include "assert_helpers.h"
#include "ds.h"
#include "btypes.h"

/**
 * Do a binary search using the suffix of 'host' beginning at offset
 * 'qry' as the query and 'sa' as an already-lexicographically-sorted
 * list of suffixes of host.  'sa' may be all suffixes of host or just
 * a subset.  Returns the index in sa of the smallest suffix of host
 * that is larger than qry, or length(sa) if all suffixes of host are
 * less than qry.
 *
 * We use the Manber and Myers optimization of maintaining a pair of
 * counters for the longest lcp observed so far on the left- and right-
 * hand sides and using the min of the two as a way of skipping over
 * characters at the beginning of a new round.
 *
 * Returns maximum value if the query suffix matches an element of sa.
 */
template<typename TStr, typename TSufElt> inline
TIndexOffU binarySASearch(
	const TStr& host,
	TIndexOffU qry,
	const EList<TSufElt>& sa)
{
	TIndexOffU lLcp = 0, rLcp = 0; // greatest observed LCPs on left and right
	TIndexOffU l = 0, r = (TIndexOffU)sa.size()+1; // binary-search window
	TIndexOffU hostLen = (TIndexOffU)host.length();
	while(true) {
		assert_gt(r, l);
		TIndexOffU m = (l+r) >> 1;
		if(m == l) {
			// Binary-search window has closed: we have an answer
			if(m > 0 && sa[m-1] == qry) {
				return std::numeric_limits<TIndexOffU>::max(); // qry matches
			}
			assert_leq(m, sa.size());
			return m; // Return index of right-hand suffix
		}
		assert_gt(m, 0);
		TIndexOffU suf = sa[m-1];
		if(suf == qry) {
			return std::numeric_limits<TIndexOffU>::max(); // query matches an elt of sa
		}
		TIndexOffU lcp = min(lLcp, rLcp);
#ifndef NDEBUG
		if(sstr_suf_upto_neq(host, qry, host, suf, lcp)) {
			assert(0);
		}
#endif
		// Keep advancing lcp, but stop when query mismatches host or
		// when the counter falls off either the query or the suffix
		while(suf+lcp < hostLen && qry+lcp < hostLen && host[suf+lcp] == host[qry+lcp]) {
			lcp++;
		}
		// Fell off the end of either the query or the sa elt?
		bool fell = (suf+lcp == hostLen || qry+lcp == hostLen);
		if((fell && qry+lcp == hostLen) || (!fell && host[suf+lcp] < host[qry+lcp])) {
			// Query is greater than sa elt
			l = m;                 // update left bound
			lLcp = max(lLcp, lcp); // update left lcp
		}
		else if((fell && suf+lcp == hostLen) || (!fell && host[suf+lcp] > host[qry+lcp])) {
			// Query is less than sa elt
			r = m;                 // update right bound
			rLcp = max(rLcp, lcp); // update right lcp
		} else {
			assert(false); // Must be one or the other!
		}
	}
	// Shouldn't get here
	assert(false);
	return std::numeric_limits<TIndexOffU>::max();
}

#endif /*BINARY_SA_SEARCH_H_*/
