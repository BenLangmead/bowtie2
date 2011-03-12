#ifndef BINARY_SA_SEARCH_H_
#define BINARY_SA_SEARCH_H_

#include <stdint.h>
#include <iostream>
#include <seqan/sequence.h>
#include "alphabet.h"
#include "assert_helpers.h"

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
 * Returns 0xffffffff if the query suffix matches an element of sa.
 */
template<typename TStr, typename TSufElt> inline
uint32_t binarySASearch(const TStr& host,
                        uint32_t qry,
                        const String<TSufElt>& sa)
{
	uint32_t lLcp = 0, rLcp = 0; // greatest observed LCPs on left and right
	uint32_t l = 0, r = length(sa)+1; // binary-search window
	uint32_t hostLen = length(host);
	while(true) {
		assert_gt(r, l);
		uint32_t m = (l+r) >> 1;
		if(m == l) {
			// Binary-search window has closed: we have an answer
			if(m > 0 && sa[m-1] == qry) return 0xffffffff; // qry matches
			assert_leq(m, length(sa));
			return m; // Return index of right-hand suffix
		}
		assert_gt(m, 0);
		uint32_t suf = sa[m-1];
		if(suf == qry) return 0xffffffff; // query matches an elt of sa
		uint32_t lcp = min(lLcp, rLcp);
#ifndef NDEBUG
		if(prefix(suffix(host, qry), lcp) != prefix(suffix(host, suf), lcp)) {
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
	return 0xffffffff;
}

#endif /*BINARY_SA_SEARCH_H_*/
