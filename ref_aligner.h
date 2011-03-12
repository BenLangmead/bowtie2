/*
 * ref_aligner.h
 */

#ifndef REF_ALIGNER_H_
#define REF_ALIGNER_H_

#include <stdint.h>
#include <iostream>
#include <stdexcept>
#include "alphabet.h"
#include "range.h"
#include "reference.h"
#include "search_globals.h"
#include "sstring.h"

// Let the reference-aligner buffer size be 16K by default.  If more
// room is required, a new buffer must be allocated from the heap.
const static int REF_ALIGNER_BUFSZ = 16 * 1024;

struct RefAlignerHit {

	RefAlignerHit() : off(), stratum(), cost(), edits(16) { }

	void clear() {
		off = stratum = cost = 0;
		edits.clear();
	}

	uint32_t off;     // 0-based offset into reference
	int      stratum; // stratum
	uint16_t cost;    // cost (incl. stratum)
	EList<Edit> edits;
};

/**
 * Abstract parent class for classes that look for alignments by
 * matching against the reference sequence directly.  This is useful
 * both for sanity-checking results from the Bowtie index and for
 * finding mates when the reference location of the opposite mate is
 * known.
 */
template<typename TStr>
class RefAligner {

	typedef BTDnaString TDna5Str;
	typedef BTString TCharStr;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:
	RefAligner(uint32_t seedLen = 0,
	           uint32_t qualMax = 0xffffffff) :
		seedLen_(seedLen),
		qualMax_(qualMax), refbuf_(buf_),
		refbufSz_(REF_ALIGNER_BUFSZ), freeRefbuf_(false)
		{ }

	/**
	 * Free the reference-space alignment buffer if this object
	 * allocated it.
	 */
	virtual ~RefAligner() {
		if(freeRefbuf_) {
			delete[] refbuf_;
		}
	}

	/**
	 * Find one alignment of qry:quals in the range begin-end in
	 * reference string ref.  Store the alignment details in range.
	 */
	virtual void find(uint32_t numToFind,
	                  const uint32_t tidx,
	                  const BitPairReference *refs,
	                  const TDna5Str& qry,
	                  const TCharStr& quals,
	                  uint32_t begin,
	                  uint32_t end,
	                  EList<RefAlignerHit>& results,
	                  TSetPairs* pairs = NULL,
	                  uint32_t aoff = 0xffffffff,
	                  bool seedOnLeft = false)
	{
		assert_gt(numToFind, 0);
		assert_gt(end, begin);
		uint32_t spread = end - begin + (gColor ? 1 : 0);
		uint32_t spreadPlus = spread + 12;
		// Make sure the buffer is large enough to accommodate the spread
		if(spreadPlus > this->refbufSz_) {
			this->newBuf(spreadPlus);
		}
		// Read in the relevant stretch of the reference string
		int offset = refs->getStretch(this->refbuf_, tidx, begin, spread);
		uint8_t *buf = ((uint8_t*)this->refbuf_) + offset;
		if(gColor) {
			// Colorize buffer
			for(size_t i = 0; i < (end-begin); i++) {
				assert_leq((int)buf[i], 4);
				buf[i] = dinuc2color[(int)buf[i]][(int)buf[i+1]];
			}
		}
		// Look for alignments
		anchor64Find(numToFind, tidx, buf, qry, quals, begin,
		             end, results, pairs, aoff, seedOnLeft);
	}

	/**
	 * Find one alignment of qry:quals in the range begin-end in
	 * reference string ref.  Store the alignment details in range.
	 * Uses a combination of the anchor bases and 64-bit arithmetic to
	 * find anchors quickly.
	 */
	virtual void anchor64Find(uint32_t numToFind,
	                uint32_t tidx,
	                uint8_t* ref,
	                const TDna5Str& qry,
	                const TCharStr& quals,
	                uint32_t begin,
	                uint32_t end,
	                EList<RefAlignerHit>& results,
	                TSetPairs* pairs = NULL,
	                uint32_t aoff = 0xffffffff,
	                bool seedOnLeft = false) = 0;

	/**
	 * Set a new reference-sequence buffer.
	 */
	void setBuf(uint32_t *newbuf, uint32_t newsz) {
		if(freeRefbuf_) {
			delete[] refbuf_;
			freeRefbuf_ = false;
		}
		refbuf_ = newbuf;
		refbufSz_ = newsz;
	}

	/**
	 * Set a new reference-sequence buffer.
	 */
	void newBuf(uint32_t newsz) {
		if(freeRefbuf_) {
			delete[] refbuf_;
		}
		try {
			refbuf_ = new uint32_t[(newsz + 3) / 4];
			if(refbuf_ == NULL) throw std::bad_alloc();
		} catch(std::bad_alloc& e) {
			cerr << "Error: Could not allocate reference-space alignment buffer of " << newsz << "B" << endl;
			throw 1;
		}
		refbufSz_ = newsz;
		freeRefbuf_ = true;
	}

protected:
	uint32_t  seedLen_;   /// length of seed region for read
	uint32_t  qualMax_;   /// maximum sum of quality penalties
	uint32_t *refbuf_;    /// pointer to current reference buffer
	uint32_t  refbufSz_;  /// size of current reference buffer
	uint32_t  buf_[REF_ALIGNER_BUFSZ / 4]; /// built-in reference buffer (may be superseded)
	bool      freeRefbuf_; /// whether refbuf_ points to something we should delete
};

/**
 * Concrete RefAligner for finding nearby exact hits given an anchor
 * hit.
 */
template<typename TStr>
class ExactRefAligner : public RefAligner<TStr> {

	typedef BTDnaString TDna5Str;
	typedef BTString TCharStr;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:

	virtual ~ExactRefAligner() { }

protected:
	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	void naiveFind(uint32_t numToFind,
	               uint32_t tidx,
	               uint8_t* ref,
	               const TDna5Str& qry,
	               const TCharStr& quals,
	               uint32_t begin,
	               uint32_t end,
	               EList<RefAlignerHit>& results,
	               TSetPairs* pairs,
	               uint32_t aoff,
	               bool seedOnLeft)
	{
		assert_gt(numToFind, 0);
		const uint32_t qlen = qry.length();
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
		uint32_t qend = end - qlen;
		uint32_t lim = qend - begin;
		uint32_t halfway = begin + (lim >> 1);
		bool hi = false;
		for(uint32_t i = 1; i <= lim+1; i++) {
			uint32_t ri;  // leftmost position in candidate alignment
			uint32_t rir; // same, minus begin; for indexing into ref[]
			if(hi) {
				ri = halfway + (i >> 1); rir = ri - begin;
				assert_leq(ri, qend);
			} else {
				ri = halfway - (i >> 1); rir = ri - begin;
				assert_geq(ri, begin);
			}
			hi = !hi;
			// Do the naive comparison
			bool match = true;
			for(uint32_t j = 0; j < qlen; j++) {
#if 0
				// Count Ns in the reference as mismatches
				const int q = (int)qry[j];
				const int r = (int)ref[rir + j];
				assert_leq(q, 4);
				assert_leq(r, 4);
				if(q == 4 || r == 4 || q != r) {
					// Mismatch
					match = false;
					break;
				}
				// Match; continue
#else
				// Disallow alignments that involve an N in the
				// reference
				const int r = (int)ref[rir + j];
				if(r & 4) {
					// N in reference; bail
					match = false;
					break;
				}
				const int q = (int)qry[j];
				assert_leq(q, 4);
				assert_lt(r, 4);
				if(q != r) {
					// Mismatch
					match = false;
					break;
				}
				// Match; continue
#endif
			}
			if(match) {
				results.expand();
				RefAlignerHit& result = results.back();
				result.clear();
				result.off = ri;
			}
		}
		return; // no match
	}

	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	virtual void anchor64Find(uint32_t numToFind,
	                uint32_t tidx,
	                uint8_t *ref,
	                const TDna5Str& qry,
	                const TCharStr& quals,
	                uint32_t begin,
	                uint32_t end,
	                EList<RefAlignerHit>& results,
	                TSetPairs* pairs,
	                uint32_t aoff, // offset of anchor mate
	                bool seedOnLeft)
	{
		assert_gt(numToFind, 0);
		ASSERT_ONLY(const uint32_t resultsISz = results.size());
		ASSERT_ONLY(uint32_t duplicates = 0);
		ASSERT_ONLY(uint32_t r2i = 0);
		const uint32_t qlen = qry.length();
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
#ifndef NDEBUG
		// Get all naive hits
		EList<RefAlignerHit> r2;
		naiveFind(numToFind, tidx, ref, qry, quals, begin, end,
		          r2, pairs, aoff, seedOnLeft);
#endif
		const uint32_t anchorBitPairs = min<int>(qlen, 32);
		// anchorOverhang = # read bases not included in the anchor
		const uint32_t anchorOverhang = qlen <= 32 ? 0 : qlen - 32;
		const uint32_t lim = end - qlen - begin;
		const uint32_t halfway = begin + (lim >> 1);
		uint64_t anchor = 0llu;
		uint64_t buffw = 0llu;
		// Set up a mask that we'll apply to the two bufs every round
		// to discard bits that were rotated out of the anchor area
		uint64_t clearMask = 0xffffffffffffffffllu;
		bool useMask = false;
		if(anchorBitPairs < 32) {
			clearMask >>= ((32-anchorBitPairs) << 1);
			useMask = true;
		}
		const int lhsShift = ((anchorBitPairs - 1) << 1);
		// Build the contents of the 'anchor' dword and the initial
		// contents of the 'buffw' dword.  If there are fewer than 32
		// anchorBitPairs, the content will be packed into the least
		// significant bits of the word.
		uint32_t skipLeftToRights = 0;
		uint32_t skipRightToLefts = 0;
		for(uint32_t i = 0; i < anchorBitPairs; i++) {
			int c = (int)qry[i]; // next query character
			assert_leq(c, 4);
			if(c & 4) {
				assert_eq(r2.size(), results.size() - resultsISz);
				return; // can't match if query has Ns
			}
			int r = (int)ref[halfway - begin + i]; // next reference character
			if(r & 4) {
				r = 0;
				// The right-to-left direction absorbs the candidate
				// alignment based at 'halfway'; so skipLeftToRights is
				// i, not i+1
				skipLeftToRights = max(skipLeftToRights, i);
				skipRightToLefts = max(skipRightToLefts, anchorBitPairs - i);
			}
			assert_lt(r, 4);
			assert_lt(c, 4);
			anchor  = ((anchor  << 2llu) | c);
			buffw = ((buffw << 2llu) | r);
		}
		// Check whether read is disqualified by Ns outside of the anchor
		// region
		for(uint32_t i = anchorBitPairs; i < qlen; i++) {
			if((int)qry[i] == 4) {
				assert_eq(r2.size(), results.size() - resultsISz);
				return; // can't match if query has Ns
			}
		}
		uint64_t bufbw = buffw;
		// We're moving the right-hand edge of the anchor along until
		// it's 'anchorOverhang' chars from the end of the target region.
		// Note that we're not making a 3'/5' distinction here; if we
		// were, we might need to make the 'anchorOverhang' adjustment on
		// the left end of the range rather than the right.
		bool hi = false;
		uint32_t riHi  = halfway;
		uint32_t rirHi = halfway - begin;
		uint32_t rirHiAnchor = rirHi + anchorBitPairs - 1;
		uint32_t riLo  = halfway + 1;
		uint32_t rirLo = halfway - begin + 1;
		for(uint32_t i = 1; i <= lim + 1; i++) {
			int r;       // new reference char
			assert_lt(skipLeftToRights, qlen);
			assert_leq(skipRightToLefts, qlen);
			if(hi) {
				hi = false;
				// Moving left-to-right
				riHi++; rirHi++; rirHiAnchor++;
				r = (int)ref[rirHiAnchor];
				if(r & 4) {
					r = 0;
					skipLeftToRights = anchorBitPairs;
				}
				assert_lt(r, 4);
				// Bring in new base pair at the least significant
				// position
				buffw = ((buffw << 2llu) | r);
				if(useMask) buffw &= clearMask;
				if(skipLeftToRights > 0) {
					skipLeftToRights--;
					continue;
				}
				if(buffw != anchor) {
					continue;
				}
			} else {
				hi = true;
				// Moving right-to-left
				riLo--; rirLo--;
				r = (int)ref[rirLo];
				if(r & 4) {
					r = 0;
					skipRightToLefts = qlen;
				}
				assert_lt(r, 4);
				if(i >= 2) {
					bufbw >>= 2llu;
					// Bring in new base pair at the most significant
					// position
					bufbw |= ((uint64_t)r << lhsShift);
				}
				if(skipRightToLefts > 0) {
					skipRightToLefts--;
					continue;
				}
				if(bufbw != anchor) {
					continue;
				}
			}
			// Seed hit!
			bool foundHit = true;
			uint32_t ri = hi ? riLo : riHi;
			uint32_t rir = hi ? rirLo : rirHi;
			if(anchorOverhang > 0) {
				// Does the non-anchor part of the alignment (the
				// "overhang") ruin it?
				bool skipCandidate = false;
				for(uint32_t j = 0; j < anchorOverhang; j++) {
					assert_lt(ri + anchorBitPairs + j, end);
					int rc = (int)ref[rir + anchorBitPairs + j];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							skipRightToLefts = anchorOverhang - j - 1;
						} else {
							// Left-to-right
							skipLeftToRights = anchorBitPairs + j;
						}
						skipCandidate = true;
						break; // Skip this candidate
					}
					if((int)qry[32 + j] != rc) {
						// Yes, overhang ruins it
						foundHit = false;
						break;
					}
				}
				if(skipCandidate) continue;
			}
			if(foundHit) {
				if(pairs != NULL) {
					TU64Pair p;
					if(ri < aoff) {
						// By convention, the upstream mate's
						// coordinates go in the 'first' field
						p.first  = ((uint64_t)tidx << 32) | (uint64_t)ri;
						p.second = ((uint64_t)tidx << 32) | (uint64_t)aoff;
					} else {
						p.first  = ((uint64_t)tidx << 32) | (uint64_t)aoff;
						p.second = ((uint64_t)tidx << 32) | (uint64_t)ri;
					}
					if(pairs->find(p) != pairs->end()) {
						// We already found this hit!  Continue.
						ASSERT_ONLY(duplicates++);
						ASSERT_ONLY(r2i++);
						continue;
					} else {
						// Record this hit
						pairs->insert(p);
					}
				}
				assert_lt(r2i, r2.size());
				assert_eq(r2[r2i].off, ri);
				results.expand();
				RefAlignerHit& result = results.back();
				result.clear();
				result.stratum = result.cost = 0;
				result.off = ri;
				ASSERT_ONLY(r2i++);
				if(--numToFind == 0) return;
			} else {
				// Keep scanning
			}
		}
		assert_leq(duplicates, r2.size());
		assert_geq(r2.size() - duplicates, results.size() - resultsISz);
		return; // no hit
	}
};

/**
 * Defined in ref_aligner.cpp.  Maps an octet representing the XOR of
 * two two-bit-per-base-encoded DNA sequences to the number of bases
 * that mismatch between the two.
 */
extern unsigned char u8toMms[];

/**
 * Concrete RefAligner for finding nearby 1-mismatch hits given an
 * anchor hit.
 */
template<typename TStr>
class OneMMRefAligner : public RefAligner<TStr> {

	typedef BTDnaString TDna5Str;
	typedef BTString TCharStr;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:

	virtual ~OneMMRefAligner() { }

protected:
	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	void naiveFind(uint32_t numToFind,
	               uint32_t tidx,
	               uint8_t* ref,
	               const TDna5Str& qry,
	               const TCharStr& quals,
	               uint32_t begin,
	               uint32_t end,
	               EList<RefAlignerHit>& results,
	               TSetPairs* pairs,
	               uint32_t aoff,
	               bool seedOnLeft)
	{
		assert_gt(numToFind, 0);
		const uint32_t qlen = qry.length();
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
		uint32_t qend = end - qlen;
		uint32_t lim = qend - begin;
		uint32_t halfway = begin + (lim >> 1);
		bool hi = false;
		for(uint32_t i = 1; i <= lim+1; i++) {
			uint32_t ri;  // leftmost position in candidate alignment
			uint32_t rir; // same, minus begin; for indexing into ref[]
			if(hi) {
				ri = halfway + (i >> 1); rir = ri - begin;
				assert_leq(ri, qend);
			} else {
				ri = halfway - (i >> 1); rir = ri - begin;
				assert_geq(ri, begin);
			}
			hi = !hi;
			// Do the naive comparison
			bool match = true;
			int refc = -1;
			int readc = -1;
			uint32_t mmOff = 0xffffffff;
			int mms = 0;
			for(uint32_t j = 0; j < qlen; j++) {
#if 0
				// Count Ns in the reference as mismatches
				const int q = (int)qry[j];
				const int r = (int)ref[rir + j];
				assert_leq(q, 4);
				assert_leq(r, 4);
				if(q == 4 || r == 4 || q != r) {
#else
				// Disallow alignments that involve an N in the
				// reference
				const int r = (int)ref[rir + j];
				if(r & 4) {
					// N in reference; bail
					match = false;
					break;
				}
				const int q = (int)qry[j];
				assert_leq(q, 4);
				assert_lt(r, 4);
				if(q != r) {
#endif
					// Mismatch!
					if(++mms > 1) {
						// Too many; reject this alignment
						match = false;
						break;
					} else {
						// First one; remember offset and ref char
						refc = "ACGT"[r];
						readc = "ACGTN"[q];
						mmOff = j;
					}
				}
	 		}
			if(match) {
				assert_leq(mms, 1);
				results.expand();
				RefAlignerHit& result = results.back();
				result.clear();
				result.stratum = mms;
				result.cost = 0;
				if(mms == 1) {
					assert_lt(mmOff, qlen);
					assert_in(refc, "ACGT");
					assert_in(readc, "ACGTN");
					result.edits.push_back(Edit(mmOff, refc, readc, EDIT_TYPE_MM));
				}
				result.off = ri;
			}
		}
		return;
	}

	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	virtual void anchor64Find(uint32_t numToFind,
	                uint32_t tidx,
	                uint8_t* ref,
	                const TDna5Str& qry,
	                const TCharStr& quals,
	                uint32_t begin,
	                uint32_t end,
	                EList<RefAlignerHit>& results,
	                TSetPairs* pairs = NULL,
	                uint32_t aoff = 0xffffffff,
	                bool seedOnLeft = false)
	{
		assert_gt(numToFind, 0);
		ASSERT_ONLY(const uint32_t resultsISz = results.size());
		ASSERT_ONLY(uint32_t duplicates = 0);
		ASSERT_ONLY(uint32_t r2i = 0);
		const uint32_t qlen = qry.length();
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
#ifndef NDEBUG
		// Get results from the naive matcher for sanity-checking
		EList<RefAlignerHit> r2;
		naiveFind(numToFind, tidx, ref, qry, quals, begin, end, r2,
		          pairs, aoff, seedOnLeft);
#endif
		const uint32_t anchorBitPairs = min<int>(qlen, 32);
		const int lhsShift = ((anchorBitPairs - 1) << 1);
		const uint32_t anchorCushion  = 32 - anchorBitPairs;
		// anchorOverhang = # read bases not included in the anchor
		const uint32_t anchorOverhang = (qlen <= 32 ? 0 : (qlen - 32));
		const uint32_t lim = end - qlen - begin;
		const uint32_t halfway = begin + (lim >> 1);
		uint64_t anchor = 0llu;
		uint64_t buffw = 0llu;  // rotating ref sequence buffer
		// OR the 'diff' buffer with this so that we can always count
		// 'N's as mismatches
		uint64_t diffMask = 0llu;
		// Set up a mask that we'll apply to the two bufs every round
		// to discard bits that were rotated out of the anchor area
		uint64_t clearMask = 0xffffffffffffffffllu;
		bool useMask = false;
		if(anchorBitPairs < 32) {
			clearMask >>= ((32-anchorBitPairs) << 1);
			useMask = true;
		}
		int nsInAnchor = 0;
		int nPos = -1;
		uint32_t skipLeftToRights = 0;
		uint32_t skipRightToLefts = 0;
		// Construct the 'anchor' 64-bit buffer so that it holds all of
		// the first 'anchorBitPairs' bit pairs of the query.
		for(uint32_t i = 0; i < anchorBitPairs; i++) {
			int c = (int)qry[i]; // next query character
			int r = (int)ref[halfway - begin + i]; // next reference character
			if(r & 4) {
				r = 0;
				// The right-to-left direction absorbs the candidate
				// alignment based at 'halfway'; so skipLeftToRights is
				// i, not i+1
				skipLeftToRights = max(skipLeftToRights, i);
				skipRightToLefts = max(skipRightToLefts, anchorBitPairs - i);
			}
			assert_leq(c, 4);
			assert_lt(r, 4);
			// Special case: query has an 'N'
			if(c == 4) {
				if(++nsInAnchor > 1) {
					// More than one 'N' in the anchor region; can't
					// possibly have a 1-mismatch hit anywhere
					assert_eq(r2.size(), results.size() - resultsISz);
					return;   // can't match if query has Ns
				}
				nPos = (int)i;
				// Make it look like an 'A' in the anchor
				c = 0;
				diffMask = (diffMask << 2llu) | 1llu;
			} else {
				diffMask <<= 2llu;
			}
			anchor  = ((anchor  << 2llu) | c);
			buffw = ((buffw << 2llu) | r);
		}
		// Check whether read is disqualified by Ns outside of the anchor
		// region
		for(uint32_t i = anchorBitPairs; i < qlen; i++) {
			if((int)qry[i] == 4) {
				if(++nsInAnchor > 1) {
					assert_eq(r2.size(), results.size() - resultsISz);
					return; // can't match if query has Ns
				}
			}
		}
		uint64_t bufbw = buffw;
		// We're moving the right-hand edge of the anchor along until
		// it's 'anchorOverhang' chars from the end of the target region.
		// Note that we're not making a 3'/5' distinction here; if we
		// were, we might need to make the 'anchorOverhang' adjustment on
		// the left end of the range rather than the right.
		bool hi = false;
		uint32_t riHi  = halfway;
		uint32_t rirHi = halfway - begin;
		uint32_t rirHiAnchor = rirHi + anchorBitPairs - 1;
		uint32_t riLo  = halfway + 1;
		uint32_t rirLo = halfway - begin + 1;
		for(uint32_t i = 1; i <= lim + 1; i++) {
			int r;       // new reference char
			uint64_t diff;
			assert_lt(skipLeftToRights, qlen);
			assert_leq(skipRightToLefts, qlen);
			if(hi) {
				hi = false;
				// Moving left-to-right
				riHi++;
				rirHi++;
				rirHiAnchor++;
				r = (int)ref[rirHiAnchor];
				if(r & 4) {
					r = 0;
					skipLeftToRights = anchorBitPairs;
				}
				assert_lt(r, 4);
				// Bring in new base pair at the least significant
				// position
				buffw = ((buffw << 2llu) | r);
				if(useMask) buffw &= clearMask;
				if(skipLeftToRights > 0) {
					skipLeftToRights--;
					continue;
				}
				diff = (buffw ^ anchor) | diffMask;
			} else {
				hi = true;
				// Moving right-to-left
				riLo--;
				rirLo--;
				r = (int)ref[rirLo];
				if(r & 4) {
					r = 0;
					skipRightToLefts = qlen;
				}
				assert_lt(r, 4);
				if(i >= 2) {
					bufbw >>= 2llu;
					// Bring in new base pair at the most significant
					// position
					bufbw |= ((uint64_t)r << lhsShift);
				}
				if(skipRightToLefts > 0) {
					skipRightToLefts--;
					continue;
				}
				diff = (bufbw ^ anchor) | diffMask;
			}
			if((diff & 0xffffffff00000000llu) &&
			   (diff & 0x00000000ffffffffllu)) continue;
			uint32_t ri  = hi ? riLo  : riHi;
			uint32_t rir = hi ? rirLo : rirHi;
			// Could use pop count
			uint8_t *diff8 = reinterpret_cast<uint8_t*>(&diff);
			// As a first cut, see if there are too many mismatches in
			// the first and last parts of the anchor
			uint32_t diffs = u8toMms[(int)diff8[0]] + u8toMms[(int)diff8[7]];
			if(diffs > 1) continue;
			diffs += u8toMms[(int)diff8[1]] +
			         u8toMms[(int)diff8[2]] +
			         u8toMms[(int)diff8[3]] +
			         u8toMms[(int)diff8[4]] +
			         u8toMms[(int)diff8[5]] +
			         u8toMms[(int)diff8[6]];
			uint32_t mmpos = 0xffffffff;
			int refc = -1;
			int readc = -1;
			if(diffs > 1) {
				// Too many differences
				continue;
			} else if(diffs == 1 && nPos != -1) {
				// There was one difference, but there was also one N,
				// so we already know where the difference is
				mmpos = nPos;
				refc = "ACGT"[(int)ref[rir + nPos]];
				readc = 'N';
			} else if(diffs == 1) {
				// Figure out which position mismatched
				mmpos = 31;
				if((diff & 0xffffffffllu) == 0) { diff >>= 32llu; mmpos -= 16; }
				assert_neq(0, diff);
				if((diff & 0xffffllu) == 0)     { diff >>= 16llu; mmpos -=  8; }
				assert_neq(0, diff);
				if((diff & 0xffllu) == 0)       { diff >>= 8llu;  mmpos -=  4; }
				assert_neq(0, diff);
				if((diff & 0xfllu) == 0)        { diff >>= 4llu;  mmpos -=  2; }
				assert_neq(0, diff);
				if((diff & 0x3llu) == 0)        { mmpos--; }
				assert_neq(0, diff);
				assert_geq(mmpos, 0);
				assert_lt(mmpos, 32);
				readc = "ACGT"[(anchor >> (62-mmpos*2)) & 3];
				if((diffMask >> (62-mmpos*2)) & 3) {
					readc = 'N';
				}
				mmpos -= anchorCushion;
				refc = "ACGT"[(int)ref[rir + mmpos]];
				assert_neq(refc, readc);
			}
			// Now extend the anchor into a longer alignment
			bool foundHit = true;
			if(anchorOverhang > 0) {
				assert_leq(ri + anchorBitPairs + anchorOverhang, end);
				bool skipCandidate = false;
				for(uint32_t j = 0; j < anchorOverhang; j++) {
					int rc = (int)ref[rir + 32 + j];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							skipRightToLefts = anchorOverhang - j - 1;
						} else {
							// Left-to-right
							skipLeftToRights = anchorBitPairs + j;
						}
						skipCandidate = true; // Skip this candidate
						break;
					}
					int qc = (int)qry[32 + j];
					if(qc != rc) {
						if(++diffs > 1) {
							foundHit = false;
							break;
						} else {
							mmpos = 32 + j;
							refc = "ACGT"[rc];
							readc = "ACGTN"[qc];
						}
					}
				}
				if(skipCandidate) continue;
			}
			if(!foundHit) continue;
			if(pairs != NULL) {
				TU64Pair p;
				if(ri < aoff) {
					// By convention, the upstream mate's
					// coordinates go in the 'first' field
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.second = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				} else {
					p.second = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				}
				if(pairs->find(p) != pairs->end()) {
					// We already found this hit!  Continue.
					ASSERT_ONLY(duplicates++);
					ASSERT_ONLY(r2i++);
					continue;
				} else {
					// Record this hit
					pairs->insert(p);
				}
			}
			assert_leq(diffs, 1);
			assert_lt(r2i, r2.size());
			assert_eq(r2[r2i].off, ri);
			assert_eq(diffs, r2[r2i].edits.size());
			results.expand();
			RefAlignerHit& result = results.back();
			result.clear();
			result.stratum = diffs;
			result.cost = 0;
			assert_eq(0, result.edits.size());
			if(diffs == 1) {
				assert_neq(mmpos, 0xffffffff);
				assert_eq(mmpos, r2[r2i].edits[0].pos);
				assert_in(refc, "ACGT");
				assert_in(readc, "ACGTN");
				assert_eq((unsigned)refc, r2[r2i].edits[0].chr);
				result.edits.push_back(Edit(mmpos, refc, readc, EDIT_TYPE_MM));
			}
			ASSERT_ONLY(r2i++);
			result.off = ri;
			if(--numToFind == 0) return;
		}
		assert_leq(duplicates, r2.size());
		assert_geq(r2.size() - duplicates, results.size() - resultsISz);
		return; // no hit
	}
};

/**
 * Concrete RefAligner for finding nearby 2-mismatch hits given an
 * anchor hit.
 */
template<typename TStr>
class TwoMMRefAligner : public RefAligner<TStr> {

	typedef BTDnaString TDna5Str;
	typedef BTString TCharStr;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:

	virtual ~TwoMMRefAligner() { }

protected:
	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	void naiveFind(uint32_t numToFind,
	               uint32_t tidx,
	               uint8_t* ref,
	               const TDna5Str& qry,
	               const TCharStr& quals,
	               uint32_t begin,
	               uint32_t end,
	               EList<RefAlignerHit>& results,
	               TSetPairs* pairs,
	               uint32_t aoff,
	               bool seedOnLeft)
	{
		assert_gt(numToFind, 0);
		const uint32_t qlen = qry.length();
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
		uint32_t qend = end - qlen;
		uint32_t lim = qend - begin;
		uint32_t halfway = begin + (lim >> 1);
		bool hi = false;
		for(uint32_t i = 1; i <= lim+1; i++) {
			uint32_t ri;  // leftmost position in candidate alignment
			uint32_t rir; // same, minus begin; for indexing into ref[]
			if(hi) {
				ri = halfway + (i >> 1); rir = ri - begin;
				assert_leq(ri, qend);
			} else {
				ri = halfway - (i >> 1); rir = ri - begin;
				assert_geq(ri, begin);
			}
			hi = !hi;
			// Do the naive comparison
			bool match = true;
			int refc1 = -1;
			int readc1 = -1;
			uint32_t mmOff1 = 0xffffffff;
			int refc2 = -1;
			int readc2 = -1;
			uint32_t mmOff2 = 0xffffffff;
			int mms = 0;
			for(uint32_t j = 0; j < qlen; j++) {
#if 0
				// Count Ns in the reference as mismatches
				const int q = (int)qry[j];
				const int r = (int)ref[rir + j];
				assert_leq(q, 4);
				assert_leq(r, 4);
				if(q == 4 || r == 4 || q != r) {
#else
				// Disallow alignments that involve an N in the
				// reference
				const int r = (int)ref[rir + j];
				if(r & 4) {
					// N in reference; bail
					match = false;
					break;
				}
				const int q = (int)qry[j];
				assert_leq(q, 4);
				assert_lt(r, 4);
				if(q != r) {
#endif
					// Mismatch!
					if(++mms > 2) {
						// Too many; reject this alignment
						match = false;
						break;
					} else if(mms == 2) {
						// Second one; remember offset and ref char
						refc2 = "ACGTN"[r];
						readc2 = "ACGTN"[q];
						mmOff2 = j;
					} else {
						assert_eq(1, mms);
						// First one; remember offset and ref char
						refc1 = "ACGT"[r];
						readc1 = "ACGTN"[q];
						mmOff1 = j;
					}
				}
			}
			if(match) {
				assert_leq(mms, 2);
				results.expand();
				RefAlignerHit& result = results.back();
				result.clear();
				result.stratum = mms;
				result.cost = 0;
				assert_eq(0, result.edits.size());
				if(mms > 0) {
					assert_lt(mmOff1, qlen);
					assert_in(refc1, "ACGT");
					assert_in(readc1, "ACGTN");
					result.edits.push_back(Edit(mmOff1, refc1, readc1, EDIT_TYPE_MM));
					if(mms > 1) {
						assert_eq(2, mms);
						assert_lt(mmOff2, qlen);
						assert_in(refc2, "ACGT");
						assert_in(readc2, "ACGTN");
						result.edits.push_back(Edit(mmOff2, refc2, readc2, EDIT_TYPE_MM));
					}
				}
				result.off = ri;
			}
		}
		return;
	}

	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	virtual void anchor64Find(uint32_t numToFind,
					uint32_t tidx,
					uint8_t* ref,
					const TDna5Str& qry,
					const TCharStr& quals,
					uint32_t begin,
					uint32_t end,
					EList<RefAlignerHit>& results,
					TSetPairs* pairs = NULL,
					uint32_t aoff = 0xffffffff,
					bool seedOnLeft = false)
	{
		assert_gt(numToFind, 0);
		ASSERT_ONLY(const uint32_t resultsISz = results.size());
		ASSERT_ONLY(uint32_t duplicates = 0);
		ASSERT_ONLY(uint32_t r2i = 0);
		const uint32_t qlen = qry.length();
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
#ifndef NDEBUG
		// Get results from the naive matcher for sanity-checking
		EList<RefAlignerHit> r2;
		naiveFind(numToFind, tidx, ref, qry, quals, begin, end, r2,
		          pairs, aoff, seedOnLeft);
#endif
		const uint32_t anchorBitPairs = min<int>(qlen, 32);
		const int lhsShift = ((anchorBitPairs - 1) << 1);
		const uint32_t anchorCushion  = 32 - anchorBitPairs;
		// anchorOverhang = # read bases not included in the anchor
		const uint32_t anchorOverhang = (qlen <= 32 ? 0 : (qlen - 32));
		const uint32_t lim = end - qlen - begin;
		const uint32_t halfway = begin + (lim >> 1);
		uint64_t anchor = 0llu;
		uint64_t buffw = 0llu;  // rotating ref sequence buffer
		// OR the 'diff' buffer with this so that we can always count
		// 'N's as mismatches
		uint64_t diffMask = 0llu;
		// Set up a mask that we'll apply to the two bufs every round
		// to discard bits that were rotated out of the anchor area
		uint64_t clearMask = 0xffffffffffffffffllu;
		bool useMask = false;
		if(anchorBitPairs < 32) {
			clearMask >>= ((32-anchorBitPairs) << 1);
			useMask = true;
		}
		int nsInAnchor = 0;
		uint32_t nPoss = 0;
		int nPos1 = -1;
		int nPos2 = -1;
		uint32_t skipLeftToRights = 0;
		uint32_t skipRightToLefts = 0;
		// Construct the 'anchor' 64-bit buffer so that it holds all of
		// the first 'anchorBitPairs' bit pairs of the query.
		for(uint32_t i = 0; i < anchorBitPairs; i++) {
			int c = (int)qry[i]; // next query character
			int r = (int)ref[halfway - begin + i]; // next reference character
			if(r & 4) {
				r = 0;
				// The right-to-left direction absorbs the candidate
				// alignment based at 'halfway'; so skipLeftToRights is
				// i, not i+1
				skipLeftToRights = max(skipLeftToRights, i);
				skipRightToLefts = max(skipRightToLefts, anchorBitPairs - i);
			}
			assert_leq(c, 4);
			assert_lt(r, 4);
			// Special case: query has an 'N'
			if(c == 4) {
				if(++nsInAnchor > 2) {
					// More than two 'N's in the anchor region; can't
					// possibly have a 2-mismatch hit anywhere
					return;   // can't match if query has Ns
				} else if(nsInAnchor == 2) {
					nPos2 = (int)i;
					nPoss++;
				} else {
					assert_eq(1, nsInAnchor);
					nPos1 = (int)i;
					nPoss++;
				}
				// Make it look like an 'A' in the anchor
				c = 0;
				diffMask = (diffMask << 2llu) | 1llu;
			} else {
				diffMask <<= 2llu;
			}
			anchor  = ((anchor  << 2llu) | c);
			buffw = ((buffw << 2llu) | r);
		}
		assert_leq(nPoss, 2);
		// Check whether read is disqualified by Ns outside of the anchor
		// region
		for(uint32_t i = anchorBitPairs; i < qlen; i++) {
			if((int)qry[i] == 4) {
				if(++nsInAnchor > 2) {
					return; // can't match if query has Ns
				}
			}
		}
		uint64_t bufbw = buffw;
		// We're moving the right-hand edge of the anchor along until
		// it's 'anchorOverhang' chars from the end of the target region.
		// Note that we're not making a 3'/5' distinction here; if we
		// were, we might need to make the 'anchorOverhang' adjustment on
		// the left end of the range rather than the right.
		bool hi = false;
		uint32_t riHi  = halfway;
		uint32_t rirHi = halfway - begin;
		uint32_t rirHiAnchor = rirHi + anchorBitPairs - 1;
		uint32_t riLo  = halfway + 1;
		uint32_t rirLo = halfway - begin + 1;
		uint32_t i;
		for(i = 1; i <= lim + 1; i++) {
			int r;       // new reference char
			uint64_t diff;
			assert_lt(skipLeftToRights, qlen);
			assert_leq(skipRightToLefts, qlen);
			if(hi) {
				hi = false;
				// Moving left-to-right
				riHi++;
				rirHi++;
				rirHiAnchor++;
				r = (int)ref[rirHiAnchor];
				if(r & 4) {
					r = 0;
					skipLeftToRights = anchorBitPairs;
				}
				assert_lt(r, 4);
				// Bring in new base pair at the least significant
				// position
				buffw = ((buffw << 2llu) | r);
				if(useMask) buffw &= clearMask;
				if(skipLeftToRights > 0) {
					skipLeftToRights--;
					continue;
				}
				diff = (buffw ^ anchor) | diffMask;
			} else {
				hi = true;
				// Moving right-to-left
				riLo--;
				rirLo--;
				r = (int)ref[rirLo];
				if(r & 4) {
					r = 0;
					skipRightToLefts = qlen;
				}
				assert_lt(r, 4);
				if(i >= 2) {
					bufbw >>= 2llu;
					// Bring in new base pair at the most significant
					// position
					bufbw |= ((uint64_t)r << lhsShift);
				}
				if(skipRightToLefts > 0) {
					skipRightToLefts--;
					continue;
				}
				diff = (bufbw ^ anchor) | diffMask;
			}
			if((diff & 0xfffff00000000000llu) &&
			   (diff & 0x00000ffffff00000llu) &&
			   (diff & 0x00000000000fffffllu)) continue;
			uint32_t ri  = hi ? riLo  : riHi;
			uint32_t rir = hi ? rirLo : rirHi;
			// Could use pop count
			uint8_t *diff8 = reinterpret_cast<uint8_t*>(&diff);
			// As a first cut, see if there are too many mismatches in
			// the first and last parts of the anchor
			uint32_t diffs = u8toMms[(int)diff8[0]] + u8toMms[(int)diff8[7]];
			if(diffs > 2) continue;
			diffs += u8toMms[(int)diff8[1]] +
					 u8toMms[(int)diff8[2]] +
					 u8toMms[(int)diff8[3]] +
					 u8toMms[(int)diff8[4]] +
					 u8toMms[(int)diff8[5]] +
					 u8toMms[(int)diff8[6]];
			uint32_t mmpos1 = 0xffffffff;
			int refc1 = -1;
			int readc1 = -1;
			uint32_t mmpos2 = 0xffffffff;
			int refc2 = -1;
			int readc2 = -1;
			if(diffs > 2) {
				// Too many differences
				continue;
			} else if(nPoss > 1 && diffs == nPoss) {
				// There was one difference, but there was also one N,
				// so we already know where the difference is
				mmpos1 = nPos1;
				refc1 = "ACGT"[(int)ref[rir + nPos1]];
				readc1 = 'N';
				if(nPoss == 2) {
					mmpos2 = nPos2;
					refc2 = "ACGT"[(int)ref[rir + nPos2]];
					readc2 = 'N';
				}
			} else if(diffs > 0) {
				// Figure out which position mismatched
				uint64_t diff2 = diff;
				mmpos1 = 31;
				if((diff & 0xffffffffllu) == 0) { diff >>= 32llu; mmpos1 -= 16; }
				assert_neq(0, diff);
				if((diff & 0xffffllu) == 0)     { diff >>= 16llu; mmpos1 -=  8; }
				assert_neq(0, diff);
				if((diff & 0xffllu) == 0)       { diff >>= 8llu;  mmpos1 -=  4; }
				assert_neq(0, diff);
				if((diff & 0xfllu) == 0)        { diff >>= 4llu;  mmpos1 -=  2; }
				assert_neq(0, diff);
				if((diff & 0x3llu) == 0)        { mmpos1--; }
				assert_neq(0, diff);
				assert_geq(mmpos1, 0);
				assert_lt(mmpos1, 32);
				readc1 = "ACGT"[(anchor >> (62-mmpos1*2)) & 3];
				if((diffMask >> (62-mmpos1*2)) & 3) {
					readc1 = 'N';
				}
				mmpos1 -= anchorCushion;
				refc1 = "ACGT"[(int)ref[rir + mmpos1]];
				assert_neq(readc1, refc1);
				if(diffs > 1) {
					// Figure out the second mismatched position
					ASSERT_ONLY(uint64_t origDiff2 = diff2);
					diff2 &= ~(0xc000000000000000llu >> (uint64_t)((mmpos1+anchorCushion) << 1));
					assert_neq(diff2, origDiff2);
					mmpos2 = 31;
					if((diff2 & 0xffffffffllu) == 0) { diff2 >>= 32llu; mmpos2 -= 16; }
					assert_neq(0, diff2);
					if((diff2 & 0xffffllu) == 0)     { diff2 >>= 16llu; mmpos2 -=  8; }
					assert_neq(0, diff2);
					if((diff2 & 0xffllu) == 0)       { diff2 >>= 8llu;  mmpos2 -=  4; }
					assert_neq(0, diff2);
					if((diff2 & 0xfllu) == 0)        { diff2 >>= 4llu;  mmpos2 -=  2; }
					assert_neq(0, diff2);
					if((diff2 & 0x3llu) == 0)        { mmpos2--; }
					assert_neq(0, diff2);
					assert_geq(mmpos2, 0);
					assert_lt(mmpos2, 32);
					readc2 = "ACGT"[(anchor >> (62-mmpos2*2)) & 3];
					if((diffMask >> (62-mmpos2*2)) & 3) {
						readc2 = 'N';
					}
					mmpos2 -= anchorCushion;
					assert_neq(mmpos1, mmpos2);
					refc2 = "ACGT"[(int)ref[rir + mmpos2]];
					assert_neq(readc2, refc2);
					if(mmpos2 < mmpos1) {
						uint32_t mmtmp = mmpos1;
						mmpos1 = mmpos2;
						mmpos2 = mmtmp;
						int refctmp = refc1;
						refc1 = refc2;
						refc2 = refctmp;
						int readctmp = readc1;
						readc1 = readc2;
						readc2 = readctmp;
					}
					assert_lt(mmpos1, mmpos2);
				}
			}
			// Now extend the anchor into a longer alignment
			bool foundHit = true;
			if(anchorOverhang > 0) {
				assert_leq(ri + anchorBitPairs + anchorOverhang, end);
				bool skipCandidate = false;
				for(uint32_t j = 0; j < anchorOverhang; j++) {
					int rc = (int)ref[rir + 32 + j];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							skipRightToLefts = anchorOverhang - j - 1;
						} else {
							// Left-to-right
							skipLeftToRights = anchorBitPairs + j;
						}
						skipCandidate = true; // Skip this candidate
						break;
					}
					int qc = (int)qry[32 + j];
					if(qc != rc) {
						if(++diffs > 2) {
							foundHit = false;
							break;
						} else if(diffs == 2) {
							mmpos2 = 32 + j;
							refc2 = "ACGT"[rc];
							readc2 = "ACGTN"[qc];
						} else {
							assert_eq(1, diffs);
							mmpos1 = 32 + j;
							refc1 = "ACGT"[rc];
							readc1 = "ACGTN"[qc];
						}
					}
				}
				if(skipCandidate) continue;
			}
			if(!foundHit) continue;
			if(pairs != NULL) {
				TU64Pair p;
				if(ri < aoff) {
					// By convention, the upstream mate's
					// coordinates go in the 'first' field
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.second = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				} else {
					p.second = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				}
				if(pairs->find(p) != pairs->end()) {
					// We already found this hit!  Continue.
					ASSERT_ONLY(duplicates++);
					ASSERT_ONLY(r2i++);
					continue;
				} else {
					// Record this hit
					pairs->insert(p);
				}
			}
			assert_leq(diffs, 2);
			assert_lt(r2i, r2.size());
			assert_eq(r2[r2i].off, ri);
			assert_eq(diffs, r2[r2i].edits.size());
			results.expand();
			RefAlignerHit& result = results.back();
			result.clear();
			result.stratum = diffs;
			result.cost = 0;
			assert_eq(0, result.edits.size());
			if(diffs > 0) {
				assert_neq(mmpos1, 0xffffffff);
				assert_eq(mmpos1, r2[r2i].edits[0].pos);
				assert_in(refc1, "ACGT");
				assert_in(readc1, "ACGTN");
				assert_eq((unsigned)refc1, r2[r2i].edits[0].chr);
				result.edits.push_back(Edit(mmpos1, refc1, readc1, EDIT_TYPE_MM));
				if(diffs > 1) {
					assert_neq(mmpos2, 0xffffffff);
					assert_eq(mmpos2, r2[r2i].edits[1].pos);
					assert_in(refc2, "ACGT");
					assert_in(readc2, "ACGTN");
					assert_eq((unsigned)refc2, r2[r2i].edits[1].chr);
					result.edits.push_back(Edit(mmpos2, refc2, readc2, EDIT_TYPE_MM));
				}
			}
			ASSERT_ONLY(r2i++);
			result.off = ri;
			if(--numToFind == 0) return;
		}
		assert_leq(duplicates, r2.size());
		assert_geq(r2.size() - duplicates, results.size() - resultsISz);
		return; // no hit
	}
};

/**
 * Concrete RefAligner for finding nearby 2-mismatch hits given an
 * anchor hit.
 */
template<typename TStr>
class ThreeMMRefAligner : public RefAligner<TStr> {

	typedef BTDnaString TDna5Str;
	typedef BTString TCharStr;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:

	virtual ~ThreeMMRefAligner() { }

protected:
	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	void naiveFind(uint32_t numToFind,
	               uint32_t tidx,
	               uint8_t* ref,
	               const TDna5Str& qry,
	               const TCharStr& quals,
	               uint32_t begin,
	               uint32_t end,
	               EList<RefAlignerHit>& results,
	               TSetPairs* pairs,
	               uint32_t aoff,
	               bool seedOnLeft)
	{
		assert_gt(numToFind, 0);
		const uint32_t qlen = qry.length();
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
		uint32_t qend = end - qlen;
		uint32_t lim = qend - begin;
		uint32_t halfway = begin + (lim >> 1);
		bool hi = false;
		for(uint32_t i = 1; i <= lim+1; i++) {
			uint32_t ri;  // leftmost position in candidate alignment
			uint32_t rir; // same, minus begin; for indexing into ref[]
			if(hi) {
				ri = halfway + (i >> 1); rir = ri - begin;
				assert_leq(ri, qend);
			} else {
				ri = halfway - (i >> 1); rir = ri - begin;
				assert_geq(ri, begin);
			}
			hi = !hi;
			// Do the naive comparison
			bool match = true;
			int refc1 = -1;
			int readc1 = -1;
			uint32_t mmOff1 = 0xffffffff;
			int refc2 = -1;
			int readc2 = -1;
			uint32_t mmOff2 = 0xffffffff;
			int refc3 = -1;
			int readc3 = -1;
			uint32_t mmOff3 = 0xffffffff;
			int mms = 0;
			for(uint32_t j = 0; j < qlen; j++) {
#if 0
				// Count Ns in the reference as mismatches
				const int q = (int)qry[j];
				const int r = (int)ref[rir + j];
				assert_leq(q, 4);
				assert_leq(r, 4);
				if(q == 4 || r == 4 || q != r) {
#else
				// Disallow alignments that involve an N in the
				// reference
				const int r = (int)ref[rir + j];
				if(r & 4) {
					// N in reference; bail
					match = false;
					break;
				}
				const int q = (int)qry[j];
				assert_leq(q, 4);
				assert_lt(r, 4);
				if(q != r) {
#endif
					// Mismatch!
					if(++mms > 3) {
						// Too many; reject this alignment
						match = false;
						break;
					} else if(mms == 3) {
						// Second one; remember offset and ref char
						refc3 = "ACGTN"[r];
						readc3 = "ACGTN"[q];
						mmOff3 = j;
					} else if(mms == 2) {
						// Second one; remember offset and ref char
						refc2 = "ACGTN"[r];
						readc2 = "ACGTN"[q];
						mmOff2 = j;
					} else {
						assert_eq(1, mms);
						// First one; remember offset and ref char
						refc1 = "ACGT"[r];
						readc1 = "ACGTN"[q];
						mmOff1 = j;
					}
				}
			}
			if(match) {
				assert_leq(mms, 3);
				results.expand();
				RefAlignerHit& result = results.back();
				result.clear();
				result.stratum = mms;
				result.cost = 0;
				assert_eq(0, result.edits.size());
				if(mms > 0) {
					assert_in(refc1, "ACGT");
					assert_in(readc1, "ACGTN");
					assert_lt(mmOff1, qlen);
					result.edits.push_back(Edit(mmOff1, refc1, readc1, EDIT_TYPE_MM));
					if(mms > 1) {
						assert_in(refc2, "ACGT");
						assert_in(readc2, "ACGTN");
						assert_lt(mmOff2, qlen);
						result.edits.push_back(Edit(mmOff2, refc2, readc2, EDIT_TYPE_MM));
						if(mms > 2) {
							assert_in(refc3, "ACGT");
							assert_in(readc3, "ACGTN");
							assert_eq(3, mms);
							assert_lt(mmOff3, qlen);
							result.edits.push_back(Edit(mmOff3, refc3, readc3, EDIT_TYPE_MM));
						}
					}
				}
				result.off = ri;
			}
		}
		return;
	}

	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	virtual void anchor64Find(uint32_t numToFind,
					uint32_t tidx,
					uint8_t* ref,
					const TDna5Str& qry,
					const TCharStr& quals,
					uint32_t begin,
					uint32_t end,
					EList<RefAlignerHit>& results,
					TSetPairs* pairs = NULL,
					uint32_t aoff = 0xffffffff,
					bool seedOnLeft = false)
	{
		assert_gt(numToFind, 0);
		ASSERT_ONLY(const uint32_t resultsISz = results.size());
		ASSERT_ONLY(uint32_t duplicates = 0);
		ASSERT_ONLY(uint32_t r2i = 0);
		const uint32_t qlen = qry.length();
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
#ifndef NDEBUG
		// Get results from the naive matcher for sanity-checking
		EList<RefAlignerHit> r2;
		naiveFind(numToFind, tidx, ref, qry, quals, begin, end, r2,
		          pairs, aoff, seedOnLeft);
#endif
		const uint32_t anchorBitPairs = min<int>(qlen, 32);
		const int lhsShift = ((anchorBitPairs - 1) << 1);
		const uint32_t anchorCushion  = 32 - anchorBitPairs;
		// anchorOverhang = # read bases not included in the anchor
		const uint32_t anchorOverhang = (qlen <= 32 ? 0 : (qlen - 32));
		const uint32_t lim = end - qlen - begin;
		const uint32_t halfway = begin + (lim >> 1);
		uint64_t anchor = 0llu;
		uint64_t buffw = 0llu;  // rotating ref sequence buffer
		// OR the 'diff' buffer with this so that we can always count
		// 'N's as mismatches
		uint64_t diffMask = 0llu;
		// Set up a mask that we'll apply to the two bufs every round
		// to discard bits that were rotated out of the anchor area
		uint64_t clearMask = 0xffffffffffffffffllu;
		bool useMask = false;
		if(anchorBitPairs < 32) {
			clearMask >>= ((32-anchorBitPairs) << 1);
			useMask = true;
		}
		int nsInAnchor = 0;
		uint32_t nPoss = 0;
		int nPos1 = -1;
		int nPos2 = -1;
		int nPos3 = -1;
		uint32_t skipLeftToRights = 0;
		uint32_t skipRightToLefts = 0;
		// Construct the 'anchor' 64-bit buffer so that it holds all of
		// the first 'anchorBitPairs' bit pairs of the query.
		for(uint32_t i = 0; i < anchorBitPairs; i++) {
			int c = (int)qry[i]; // next query character
			int r = (int)ref[halfway - begin + i]; // next reference character
			if(r & 4) {
				r = 0;
				// The right-to-left direction absorbs the candidate
				// alignment based at 'halfway'; so skipLeftToRights is
				// i, not i+1
				skipLeftToRights = max(skipLeftToRights, i);
				skipRightToLefts = max(skipRightToLefts, anchorBitPairs - i);
			}
			assert_leq(c, 4);
			assert_lt(r, 4);
			// Special case: query has an 'N'
			if(c == 4) {
				if(++nsInAnchor > 3) {
					// More than two 'N's in the anchor region; can't
					// possibly have a 2-mismatch hit anywhere
					assert_eq(r2.size(), results.size() - resultsISz);
					return;   // can't match if query has Ns
				} else if(nsInAnchor == 3) {
					nPos3 = (int)i;
					nPoss++;
				} else if(nsInAnchor == 2) {
					nPos2 = (int)i;
					nPoss++;
				} else {
					assert_eq(1, nsInAnchor);
					nPos1 = (int)i;
					nPoss++;
				}
				// Make it look like an 'A' in the anchor
				c = 0;
				diffMask = (diffMask << 2llu) | 1llu;
			} else {
				diffMask <<= 2llu;
			}
			anchor  = ((anchor  << 2llu) | c);
			buffw = ((buffw << 2llu) | r);
		}
		assert_leq(nPoss, 3);
		// Check whether read is disqualified by Ns outside of the anchor
		// region
		for(uint32_t i = anchorBitPairs; i < qlen; i++) {
			if((int)qry[i] == 4) {
				if(++nsInAnchor > 3) {
					assert_eq(r2.size(), results.size() - resultsISz);
					return; // can't match if query has Ns
				}
			}
		}
		uint64_t bufbw = buffw;
		// We're moving the right-hand edge of the anchor along until
		// it's 'anchorOverhang' chars from the end of the target region.
		// Note that we're not making a 3'/5' distinction here; if we
		// were, we might need to make the 'anchorOverhang' adjustment on
		// the left end of the range rather than the right.
		bool hi = false;
		uint32_t riHi  = halfway;
		uint32_t rirHi = halfway - begin;
		uint32_t rirHiAnchor = rirHi + anchorBitPairs - 1;
		uint32_t riLo  = halfway + 1;
		uint32_t rirLo = halfway - begin + 1;
		for(uint32_t i = 1; i <= lim + 1; i++) {
			int r;       // new reference char
			uint64_t diff;
			assert_lt(skipLeftToRights, qlen);
			assert_leq(skipRightToLefts, qlen);
			if(hi) {
				hi = false;
				// Moving left-to-right
				riHi++;
				rirHi++;
				rirHiAnchor++;
				r = (int)ref[rirHiAnchor];
				if(r & 4) {
					r = 0;
					skipLeftToRights = anchorBitPairs;
				}
				assert_lt(r, 4);
				// Bring in new base pair at the least significant
				// position
				buffw = ((buffw << 2llu) | r);
				if(useMask) buffw &= clearMask;
				if(skipLeftToRights > 0) {
					skipLeftToRights--;
					continue;
				}
				diff = (buffw ^ anchor) | diffMask;
			} else {
				hi = true;
				// Moving right-to-left
				riLo--;
				rirLo--;
				r = (int)ref[rirLo];
				if(r & 4) {
					r = 0;
					skipRightToLefts = qlen;
				}
				assert_lt(r, 4);
				if(i >= 2) {
					bufbw >>= 2llu;
					// Bring in new base pair at the most significant
					// position
					bufbw |= ((uint64_t)r << lhsShift);
				}
				if(skipRightToLefts > 0) {
					skipRightToLefts--;
					continue;
				}
				diff = (bufbw ^ anchor) | diffMask;
			}
			if((diff & 0xffff000000000000llu) &&
			   (diff & 0x0000ffff00000000llu) &&
			   (diff & 0x00000000ffff0000llu) &&
			   (diff & 0x000000000000ffffllu)) continue;
			uint32_t ri  = hi ? riLo  : riHi;
			uint32_t rir = hi ? rirLo : rirHi;
			// Could use pop count
			uint8_t *diff8 = reinterpret_cast<uint8_t*>(&diff);
			// As a first cut, see if there are too many mismatches in
			// the first and last parts of the anchor
			uint32_t diffs = u8toMms[(int)diff8[0]] + u8toMms[(int)diff8[7]];
			if(diffs > 3) continue;
			diffs += u8toMms[(int)diff8[1]] +
					 u8toMms[(int)diff8[2]] +
					 u8toMms[(int)diff8[3]] +
					 u8toMms[(int)diff8[4]] +
					 u8toMms[(int)diff8[5]] +
					 u8toMms[(int)diff8[6]];
			uint32_t mmpos1 = 0xffffffff;
			int refc1 = -1;
			int readc1 = -1;
			uint32_t mmpos2 = 0xffffffff;
			int refc2 = -1;
			int readc2 = -1;
			uint32_t mmpos3 = 0xffffffff;
			int refc3 = -1;
			int readc3 = -1;
			if(diffs > 3) {
				// Too many differences
				continue;
			} else if(nPoss > 1 && diffs == nPoss) {
				// There was one difference, but there was also one N,
				// so we already know where the difference is
				mmpos1 = nPos1;
				refc1 = "ACGT"[(int)ref[rir + nPos1]];
				readc1 = 'N';
				if(nPoss > 1) {
					mmpos2 = nPos2;
					refc2 = "ACGT"[(int)ref[rir + nPos2]];
					readc2 = 'N';
					if(nPoss > 2) {
						mmpos3 = nPos3;
						refc3 = "ACGT"[(int)ref[rir + nPos3]];
						readc3 = 'N';
					}
				}
			} else if(diffs > 0) {
				// Figure out which position mismatched
				uint64_t diff2 = diff;
				mmpos1 = 31;
				if((diff & 0xffffffffllu) == 0) { diff >>= 32llu; mmpos1 -= 16; }
				assert_neq(0, diff);
				if((diff & 0xffffllu) == 0)     { diff >>= 16llu; mmpos1 -=  8; }
				assert_neq(0, diff);
				if((diff & 0xffllu) == 0)       { diff >>= 8llu;  mmpos1 -=  4; }
				assert_neq(0, diff);
				if((diff & 0xfllu) == 0)        { diff >>= 4llu;  mmpos1 -=  2; }
				assert_neq(0, diff);
				if((diff & 0x3llu) == 0)        { mmpos1--; }
				assert_neq(0, diff);
				assert_geq(mmpos1, 0);
				assert_lt(mmpos1, 32);
				readc1 = "ACGT"[(anchor >> (62-mmpos1*2)) & 3];
				if((diffMask >> (62-mmpos1*2)) & 3) {
					readc1 = 'N';
				}
				mmpos1 -= anchorCushion;
				refc1 = "ACGT"[(int)ref[rir + mmpos1]];
				assert_neq(readc1, refc1);
				if(diffs > 1) {
					// Figure out the second mismatched position
					diff2 &= ~(0xc000000000000000llu >> (uint64_t)((mmpos1 + anchorCushion) << 1));
					uint64_t diff3 = diff2;
					mmpos2 = 31;
					if((diff2 & 0xffffffffllu) == 0) { diff2 >>= 32llu; mmpos2 -= 16; }
					assert_neq(0, diff2);
					if((diff2 & 0xffffllu) == 0)     { diff2 >>= 16llu; mmpos2 -=  8; }
					assert_neq(0, diff2);
					if((diff2 & 0xffllu) == 0)       { diff2 >>= 8llu;  mmpos2 -=  4; }
					assert_neq(0, diff2);
					if((diff2 & 0xfllu) == 0)        { diff2 >>= 4llu;  mmpos2 -=  2; }
					assert_neq(0, diff2);
					if((diff2 & 0x3llu) == 0)        { mmpos2--; }
					assert_neq(0, diff2);
					readc2 = "ACGT"[(anchor >> (62-mmpos2*2)) & 3];
					if((diffMask >> (62-mmpos2*2)) & 3) {
						readc2 = 'N';
					}
					mmpos2 -= anchorCushion;
					assert_geq(mmpos2, 0);
					assert_lt(mmpos2, 32);
					assert_neq(mmpos1, mmpos2);
					refc2 = "ACGT"[(int)ref[rir + mmpos2]];
					assert_neq(readc2, refc2);
					uint32_t mmpos2orig = mmpos2;
					if(mmpos2 < mmpos1) {
						uint32_t mmtmp = mmpos1;
						mmpos1 = mmpos2;
						mmpos2 = mmtmp;
						int refctmp = refc1;
						refc1 = refc2;
						refc2 = refctmp;
						int readctmp = readc1;
						readc1 = readc2;
						readc2 = readctmp;
					}
					assert_lt(mmpos1, mmpos2);
					if(diffs > 2) {
						// Figure out the second mismatched position
						diff3 &= ~(0xc000000000000000llu >> (uint64_t)((mmpos2orig + anchorCushion) << 1));
						mmpos3 = 31;
						if((diff3 & 0xffffffffllu) == 0) { diff3 >>= 32llu; mmpos3 -= 16; }
						assert_neq(0, diff3);
						if((diff3 & 0xffffllu) == 0)     { diff3 >>= 16llu; mmpos3 -=  8; }
						assert_neq(0, diff3);
						if((diff3 & 0xffllu) == 0)       { diff3 >>= 8llu;  mmpos3 -=  4; }
						assert_neq(0, diff3);
						if((diff3 & 0xfllu) == 0)        { diff3 >>= 4llu;  mmpos3 -=  2; }
						assert_neq(0, diff3);
						if((diff3 & 0x3llu) == 0)        { mmpos3--; }
						assert_neq(0, diff3);
						readc3 = "ACGT"[(anchor >> (62-mmpos3*2)) & 3];
						if((diffMask >> (62-mmpos3*2)) & 3) {
							readc3 = 'N';
						}
						mmpos3 -= anchorCushion;
						assert_geq(mmpos3, 0);
						assert_lt(mmpos3, 32);
						assert_neq(mmpos1, mmpos3);
						assert_neq(mmpos2, mmpos3);
						refc3 = "ACGT"[(int)ref[rir + mmpos3]];
						assert_neq(readc3, refc3);
						if(mmpos3 < mmpos1) {
							uint32_t mmtmp = mmpos1;
							mmpos1 = mmpos3;
							mmpos3 = mmpos2;
							mmpos2 = mmtmp;
							int refctmp = refc1;
							refc1 = refc3;
							refc3 = refc2;
							refc2 = refctmp;
							int readctmp = readc1;
							readc1 = readc3;
							readc3 = readc2;
							readc2 = readctmp;
						} else if(mmpos3 < mmpos2) {
							uint32_t mmtmp = mmpos2;
							mmpos2 = mmpos3;
							mmpos3 = mmtmp;
							int refctmp = refc2;
							refc2 = refc3;
							refc3 = refctmp;
							int readctmp = readc2;
							readc2 = readc3;
							readc3 = readctmp;
						}
						assert_lt(mmpos1, mmpos2);
						assert_lt(mmpos2, mmpos3);
					}
				}
			}
			// Now extend the anchor into a longer alignment
			bool foundHit = true;
			if(anchorOverhang > 0) {
				assert_leq(ri + anchorBitPairs + anchorOverhang, end);
				bool skipCandidate = false;
				for(uint32_t j = 0; j < anchorOverhang; j++) {
					int rc = (int)ref[rir + 32 + j];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							skipRightToLefts = anchorOverhang - j - 1;
						} else {
							// Left-to-right
							skipLeftToRights = anchorBitPairs + j;
						}
						skipCandidate = true; // Skip this candidate
						break;
					}
					int qc = (int)qry[32 + j] ;
					if(qc != rc) {
						if(++diffs > 3) {
							foundHit = false;
							break;
						} else if(diffs == 3) {
							mmpos3 = 32 + j;
							refc3 = "ACGT"[rc];
							readc3 = "ACGTN"[qc];
						} else if(diffs == 2) {
							mmpos2 = 32 + j;
							refc2 = "ACGT"[rc];
							readc2 = "ACGTN"[qc];
						} else {
							assert_eq(1, diffs);
							mmpos1 = 32 + j;
							refc1 = "ACGT"[rc];
							readc1 = "ACGTN"[qc];
						}
					}
				}
				if(skipCandidate) continue;
			}
			if(!foundHit) continue;
			if(pairs != NULL) {
				TU64Pair p;
				if(ri < aoff) {
					// By convention, the upstream mate's
					// coordinates go in the 'first' field
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.second = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				} else {
					p.second = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				}
				if(pairs->find(p) != pairs->end()) {
					// We already found this hit!  Continue.
					ASSERT_ONLY(duplicates++);
					ASSERT_ONLY(r2i++);
					continue;
				} else {
					// Record this hit
					pairs->insert(p);
				}
			}
			assert_leq(diffs, 3);
			assert_lt(r2i, r2.size());
			assert_eq(r2[r2i].off, ri);
			assert_eq(diffs, r2[r2i].edits.size());
			results.expand();
			RefAlignerHit& result = results.back();
			result.clear();
			result.stratum = diffs;
			result.cost = 0;
			assert_eq(0, result.edits.size());
			if(diffs > 0) {
				assert_neq(mmpos1, 0xffffffff);
				assert_eq(mmpos1, r2[r2i].edits[0].pos);
				assert_in(refc1, "ACGT");
				assert_in(readc1, "ACGTN");
				assert_eq((unsigned)refc1, r2[r2i].edits[0].chr);
				result.edits.push_back(Edit(mmpos1, refc1, readc1, EDIT_TYPE_MM));
				if(diffs > 1) {
					assert_neq(mmpos2, 0xffffffff);
					assert_eq(mmpos2, r2[r2i].edits[1].pos);
					assert_in(refc2, "ACGT");
					assert_in(readc2, "ACGTN");
					assert_eq((unsigned)refc2, r2[r2i].edits[1].chr);
					result.edits.push_back(Edit(mmpos2, refc2, readc2, EDIT_TYPE_MM));
					if(diffs > 2) {
						assert_neq(mmpos3, 0xffffffff);
						assert_eq(mmpos3, r2[r2i].edits[2].pos);
						assert_in(refc3, "ACGT");
						assert_in(readc3, "ACGTN");
						assert_eq((unsigned)refc3, r2[r2i].edits[2].chr);
						result.edits.push_back(Edit(mmpos3, refc3, readc3, EDIT_TYPE_MM));
					}
				}
			}
			ASSERT_ONLY(r2i++);
			result.off = ri;
			if(--numToFind == 0) return;
		}
		assert_leq(duplicates, r2.size());
		assert_geq(r2.size() - duplicates, results.size() - resultsISz);
		return; // no hit
	}
};

/**
 * Concrete RefAligner for finding nearby 1-mismatch hits given an
 * anchor hit.
 */
template<typename TStr>
class Seed0RefAligner : public RefAligner<TStr> {

	typedef BTDnaString TDna5Str;
	typedef BTString TCharStr;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:

	Seed0RefAligner(uint32_t seedLen, uint32_t qualMax) :
		RefAligner<TStr>(seedLen, qualMax) { }

	virtual ~Seed0RefAligner() { }

protected:
	/**
	 * This schematic shows the roles played by the begin, qbegin, end,
	 * qend, halfway, slen, qlen, and lim variables:
	 *
	 * seedOnLeft == true:
	 *
	 * |<                   lim                   >|<     qlen       >|
	 *  --------------------------------------------------------------
	 * |                     | slen | qlen-slen |  | slen | qlen-slen |
	 *  --------------------------------------------------------------
	 * ^                     ^                     ^                  ^
	 * begin & qbegin     halfway                qend               end
	 *
	 * seedOnLeft == false:
	 *
	 * |<     qlen       >|<                   lim                   >|
	 *  --------------------------------------------------------------
	 * | qlen-slen | slen |  | qlen-slen | slen |                     |
	 *  --------------------------------------------------------------
	 * ^                  ^                     ^                     ^
	 * begin            qbegin             halfway           qend & end
	 *
	 * Note that, for seeds longer than 32 base-pairs, the seed is
	 * further subdivided into a 32-bit anchor and a seed overhang of
	 * length > 0.
	 */
	void naiveFind(uint32_t numToFind,
				   uint32_t tidx,
				   uint8_t* ref,
				   const TDna5Str& qry,
				   const TCharStr& quals,
				   uint32_t begin,
				   uint32_t end,
				   EList<RefAlignerHit>& results,
				   TSetPairs* pairs,
				   uint32_t aoff,
				   bool seedOnLeft)
	{
		assert_gt(numToFind, 0);
		assert_gt(end, begin);
		const uint32_t qlen = qry.length();
		assert_gt(qlen, 0);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(this->seedLen_, 0);
		const uint32_t slen = min(qlen, this->seedLen_);
		uint32_t qend = end;
		uint32_t qbegin = begin;
		// If the seed is on the left-hand side of the alignment, then
		// leave a gap at the right-hand side of the interval;
		// otherwise, do the opposite
		if(seedOnLeft) {
			// Leave gap on right-hand side of the interval
			qend -= qlen;
		} else {
			// Leave gap on left-hand side of the interval
			qbegin += qlen;
		}
		// lim = number of alignments to try
		const uint32_t lim = qend - qbegin;
		// halfway = position in the reference to start at (and then
		// we work our way out to the right and to the left).
		const uint32_t halfway = qbegin + (lim >> 1);
		nonSeedEdits.clear();
		bool hi = false;
		for(uint32_t i = 1; i <= lim+1; i++) {
			uint32_t ri;  // leftmost position in candidate alignment
			uint32_t rir; // same, minus begin; for indexing into ref[]
			if(hi) {
				ri = halfway + (i >> 1); rir = ri - begin;
				assert_leq(ri, qend);
			} else {
				ri = halfway - (i >> 1); rir = ri - begin;
				assert_geq(ri, begin);
			}
			hi = !hi;
			// Do the naive comparison
			bool match = true;
			int mms = 0;
			unsigned int ham = 0;
			nonSeedEdits.clear();
			// Walk through each position of the alignment
			for(uint32_t jj = 0; jj < qlen; jj++) {
				uint32_t j = jj;
				if(!seedOnLeft) {
					// If seed is on the right, scan right-to-left
					j = qlen - jj - 1;
				} else {
					// Go left-to-right
				}
				uint32_t rirj = rir + j;
				if(!seedOnLeft) {
					assert_geq(rir, jj);
					rirj = rir - jj - 1;
				}
#if 0
				// Count Ns in the reference as mismatches
				const int q = (int)qry[j];
				const int r = (int)ref[rirj];
				assert_leq(q, 4);
				assert_leq(r, 4);
				if(q == 4 || r == 4 || q != r) {
#else
				// Disallow alignments that involve an N in the
				// reference
				const int r = (int)ref[rirj];
				if(r & 4) {
					// N in reference; bail on this alignment
					match = false;
					break;
				}
				const int q = (int)qry[j];
				assert_leq(q, 4);
				assert_lt(r, 4);
				if(q != r) {
#endif
					// Mismatch!
					if(jj < slen) {
						// More than one mismatch in the anchor; reject
						match = false;
						break;
					}
					uint8_t qual = phredcToPhredq(quals[j]);
					ham += mmPenalty(!gNoMaqRound, qual);
					if(ham > this->qualMax_) {
						// Exceeded quality ceiling; reject
						match = false;
						break;
					} else {
						// Legal mismatch outside of the anchor; record it
						mms++;
						nonSeedEdits.push_back(Edit(j, "ACGT"[r], "ACGTN"[q], EDIT_TYPE_MM));
						assert_leq(nonSeedEdits.size(), (size_t)mms);
					}
				}
			}
			if(match) {
				results.expand();
				RefAlignerHit& result = results.back();
				result.clear();
				assert_eq(0, result.edits.size());
				if(mms >= 1) {
					// Be careful to add edits in left-to-right order
					// with respect to the read/alignment
					const size_t nonSeedEditsSz = nonSeedEdits.size();
					if(nonSeedEditsSz > 0) {
						if(seedOnLeft) {
							for(size_t k = 0; k < nonSeedEditsSz; k++) {
								result.edits.push_back(nonSeedEdits[k]);
							}
						} else {
							for(size_t k = nonSeedEditsSz; k > 0; k--) {
								result.edits.push_back(nonSeedEdits[k-1]);
							}
						}
					}
					result.cost = ham;
				}
				assert_eq((size_t)mms, result.edits.size());
				if(seedOnLeft) {
					result.off = ri;
				} else {
					result.off = ri - qlen;
				}
			}
		}
		return;
	}

	/**
	 * This schematic shows the roles played by the begin, qbegin, end,
	 * qend, halfway, slen, qlen, and lim variables:
	 *
	 * seedOnLeft == true:
	 *
	 * |<                   lim                   >|<     qlen       >|
	 *  --------------------------------------------------------------
	 * |                     | slen | qlen-slen |  | slen | qlen-slen |
	 *  --------------------------------------------------------------
	 * ^                     ^                     ^                  ^
	 * begin & qbegin     halfway                qend               end
	 *
	 * seedOnLeft == false:
	 *
	 *             |<                   lim                   >|
	 *  --------------------------------------------------------------
	 * | qlen-slen |         | qlen-slen | slen |              | slen |
	 *  --------------------------------------------------------------
	 * ^           ^                     ^                     ^      ^
	 * begin       qbegin             halfway                qend   end
	 *
	 * Note that, for seeds longer than 32 base-pairs, the seed is
	 * further subdivided into a 32-bit anchor and a seed overhang of
	 * length > 0.
	 */
	virtual void anchor64Find(uint32_t numToFind,
					uint32_t tidx,
					uint8_t* ref,
					const TDna5Str& qry,
					const TCharStr& quals,
					uint32_t begin,
					uint32_t end,
					EList<RefAlignerHit>& results,
					TSetPairs* pairs = NULL,
					uint32_t aoff = 0xffffffff,
					bool seedOnLeft = false)
	{
		assert_gt(numToFind, 0);
		ASSERT_ONLY(const uint32_t resultsISz = results.size());
		ASSERT_ONLY(uint32_t duplicates = 0);
		ASSERT_ONLY(uint32_t r2i = 0);
		const uint32_t qlen = qry.length();
		assert_gt(qlen, 0);
		assert_gt(end, begin);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(this->seedLen_, 0);
		uint32_t slen = min(qlen, this->seedLen_);
#ifndef NDEBUG
		// Get results from the naive matcher for sanity-checking
		EList<RefAlignerHit> r2;
		naiveFind(numToFind, tidx, ref, qry, quals, begin, end, r2,
		          pairs, aoff, seedOnLeft);
#endif
		const uint32_t anchorBitPairs = min<int>(slen, 32);
		const int lhsShift = ((anchorBitPairs - 1) << 1);
		ASSERT_ONLY(const uint32_t anchorCushion  = 32 - anchorBitPairs);
		// seedAnchorOverhang = # seed bases not included in the anchor
		const uint32_t seedAnchorOverhang = (slen <= anchorBitPairs ? 0 : (slen - anchorBitPairs));
		// seedAnchorOverhang = # seed bases not included in the anchor
		const uint32_t readSeedOverhang = (slen == qlen ? 0 : (qlen - slen));
		assert(anchorCushion == 0 || seedAnchorOverhang == 0);
		assert_eq(qlen, readSeedOverhang + slen);
		uint32_t qend = end;
		uint32_t qbegin = begin;
		if(seedOnLeft) {
			// Leave read-sized gap on right-hand side of the interval
			qend -= qlen;
		} else {
			// Leave seed-sized gap on right-hand side and
			// non-seed-sized gap on the left-hand side
			qbegin += readSeedOverhang;
			qend -= slen;
		}
		// lim = # possible alignments in the range
		const uint32_t lim = qend - qbegin;
		// halfway = point on the genome to radiate out from
		const uint32_t halfway = qbegin + (lim >> 1);
		uint64_t anchor = 0llu;
		uint64_t buffw = 0llu;  // rotating ref sequence buffer
		// Set up a mask that we'll apply to the two bufs every round
		// to discard bits that were rotated out of the anchor area
		uint64_t clearMask = 0xffffffffffffffffllu;
		bool useMask = false;
		if(anchorBitPairs < 32) {
			clearMask >>= ((32-anchorBitPairs) << 1);
			useMask = true;
		}
		uint32_t skipLeftToRights = 0;
		uint32_t skipRightToLefts = 0;
		const uint32_t halfwayRi = halfway - begin;
		// Construct the 'anchor' 64-bit buffer so that it holds all of
		// the first 'anchorBitPairs' bit pairs of the query.
		for(uint32_t ii = 0; ii < anchorBitPairs; ii++) {
			uint32_t i = ii;
			if(!seedOnLeft) {
				// Fill in the anchor using characters from the right-
				// hand side of the query (but take the characters in
				// left-to-right order)
				i = qlen - slen + ii;
			}
			int c = (int)qry[i]; // next query character
			int r = (int)ref[halfwayRi + ii]; // next reference character
			if(r & 4) {
				// The reference character is an N; to mimic the
				// behavior of BW alignment, we have to skip all
				// alignments that involve an N in the reference.  Set
				// the skip* variables accordingly.
				r = 0;
				uint32_t lrSkips = ii;
				uint32_t rlSkips = qlen - ii;
				if(!seedOnLeft && readSeedOverhang) {
					lrSkips += readSeedOverhang;
					assert_geq(rlSkips, readSeedOverhang);
					rlSkips -= readSeedOverhang;
				}
				// The right-to-left direction absorbs the candidate
				// alignment based at 'halfway'
				skipLeftToRights = max(skipLeftToRights, lrSkips);
				skipRightToLefts = max(skipRightToLefts, rlSkips);
				assert_leq(skipLeftToRights, qlen);
				assert_leq(skipRightToLefts, qlen);
			}
			assert_leq(c, 4);
			assert_lt(r, 4);
			// Special case: query has an 'N'
			if(c == 4) {
				// One or more 'N's in the anchor region; can't
				// possibly have a 0-mismatch hit anywhere
				assert_eq(r2.size(), results.size() - resultsISz);
				return;   // can't match if query has Ns
			}
			anchor = ((anchor << 2llu) | c);
			buffw = ((buffw << 2llu) | r);
		}
		// Check whether read is disqualified by Ns inside the seed
		// region but outside the anchor region
		if(seedAnchorOverhang) {
			assert_lt(anchorBitPairs, slen);
			for(uint32_t ii = anchorBitPairs; ii < slen; ii++) {
				uint32_t i = ii;
				if(!seedOnLeft) {
					i = qlen - slen + ii;
				}
				if((int)qry[i] == 4) {
					assert_eq(r2.size(), results.size() - resultsISz);
					return; // can't match if query has Ns
				}
			}
		} else {
			assert_eq(anchorBitPairs, slen);
		}
		uint64_t bufbw = buffw;
		// Slide the anchor out in either direction, alternating
		// between right-to-left and left-to-right shifts, until all of
		// the positions from qbegin to qend have been covered.
		bool hi = false;
		uint32_t riHi  = halfway;
		uint32_t rirHi = halfway - begin;
		uint32_t rirHiAnchor = rirHi + anchorBitPairs - 1;
		uint32_t riLo  = halfway + 1;
		uint32_t rirLo = halfway - begin + 1;
		uint32_t lrSkips = anchorBitPairs;
		uint32_t rlSkips = qlen;
		if(!seedOnLeft && readSeedOverhang) {
			lrSkips += readSeedOverhang;
			assert_geq(rlSkips, readSeedOverhang);
			rlSkips -= readSeedOverhang;
		}
		for(uint32_t i = 1; i <= lim + 1; i++) {
			int r;       // new reference char
			uint64_t diff;
			assert_leq(skipLeftToRights, qlen);
			assert_leq(skipRightToLefts, qlen);
			if(hi) {
				hi = false;
				// Moving left-to-right
				riHi++;
				rirHi++;
				rirHiAnchor++;
				r = (int)ref[rirHiAnchor];
				if(r & 4) {
					r = 0;
					skipLeftToRights = lrSkips;
				}
				assert_lt(r, 4);
				// Bring in new base pair at the least significant
				// position
				buffw = ((buffw << 2llu) | r);
				if(useMask) buffw &= clearMask;
				if(skipLeftToRights > 0) {
					skipLeftToRights--;
					continue;
				}
				diff = buffw ^ anchor;
			} else {
				hi = true;
				// Moving right-to-left
				riLo--;
				rirLo--;
				r = (int)ref[rirLo];
				if(r & 4) {
					r = 0;
					skipRightToLefts = rlSkips;
				}
				assert_lt(r, 4);
				if(i >= 2) {
					bufbw >>= 2llu;
					// Bring in new base pair at the most significant
					// position
					bufbw |= ((uint64_t)r << lhsShift);
				}
				if(skipRightToLefts > 0) {
					skipRightToLefts--;
					continue;
				}
				diff = bufbw ^ anchor;
			}
			if(diff) continue;
			uint32_t ri  = hi ? riLo  : riHi;
			uint32_t rir = hi ? rirLo : rirHi;
			unsigned int ham = 0;
			// If the seed is longer than the anchor, then scan the
			// rest of the seed characters
			bool foundHit = true;
			if(seedAnchorOverhang) {
				for(uint32_t j = 0; j < seedAnchorOverhang; j++) {
					int rc = (int)ref[rir + anchorBitPairs + j];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							// Skip out of the seedAnchorOverhang
							assert_eq(0, skipRightToLefts);
							skipRightToLefts = seedAnchorOverhang - j - 1;
							if(seedOnLeft) {
								// ...and skip out of the rest of the read
								skipRightToLefts += readSeedOverhang;
							}
						} else {
							// Left-to-right
							// Skip out of the seedAnchorOverhang
							assert_eq(0, skipLeftToRights);
							skipLeftToRights = anchorBitPairs + j;
							if(!seedOnLeft) {
								// ...and skip out of the rest of the read
								skipLeftToRights += readSeedOverhang;
							}
						}
						foundHit = false; // Skip this candidate
						break;
					}
					uint32_t qoff = anchorBitPairs + j;
					if(!seedOnLeft) {
						qoff += readSeedOverhang;
					}
					assert_lt(qoff, qlen);
					if((int)qry[qoff] != rc) {
						foundHit = false;
						break;
					}
				}
				if(!foundHit) continue;
			}
			// If the read is longer than the seed, then scan the rest
			// of the read characters; mismatches no longer count
			// toward the stratum or the 1-mm limit.
			// Vectors for holding edit information
			nonSeedEdits.clear();
			int mms = 0; // start counting total mismatches
			if((qlen - slen) > 0) {
				// Going left-to-right
				for(uint32_t j = 0; j < readSeedOverhang; j++) {
					uint32_t roff = rir + slen + j;
					uint32_t qoff = slen + j;
					if(!seedOnLeft) {
						assert_geq(roff, qlen);
						roff -= qlen;
						qoff = j;
					}
					int rc = (int)ref[roff];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							// Skip what's left of the readSeedOverhang
							skipRightToLefts = readSeedOverhang - j - 1;
							if(!seedOnLeft) {
								// ...and skip the seed if it's on the right
								skipRightToLefts += slen;
							}
						} else {
							// Left-to-right
							// Skip what we've matched of the overhang
							skipLeftToRights = j;
							if(seedOnLeft) {
								// ...and skip the seed if it's on the left
								skipLeftToRights += slen;
							}
						}
						foundHit = false; // Skip this candidate
						break;
					}
					int qc = (int)qry[qoff];
					if(qc != rc) {
						// Calculate quality of mismatched base
						char q = phredcToPhredq(quals[qoff]);
						ham += mmPenalty(!gNoMaqRound, q);
						if(ham > this->qualMax_) {
							// Exceeded quality limit
							foundHit = false;
							break;
						}
						// Legal mismatch outside of the anchor; record it
						mms++;
						nonSeedEdits.push_back(Edit(qoff, "ACGT"[rc], "ACGTN"[qc], EDIT_TYPE_MM));
						assert_leq(nonSeedEdits.size(), (size_t)mms);
					}
				}
				if(!foundHit) continue;
			}
			assert(foundHit);
			// Adjust ri if seed is on the right-hand side
			if(!seedOnLeft) {
				ri -= readSeedOverhang;
				rir -= readSeedOverhang;
			}
			if(pairs != NULL) {
				TU64Pair p;
				if(ri < aoff) {
					// By convention, the upstream mate's
					// coordinates go in the 'first' field
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.second = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				} else {
					p.second = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				}
				if(pairs->find(p) != pairs->end()) {
					// We already found this hit!  Continue.
					ASSERT_ONLY(duplicates++);
					ASSERT_ONLY(r2i++);
					continue;
				} else {
					// Record this hit
					pairs->insert(p);
				}
			}
			if(gVerbose) {
				cout << "About to report seed0:" << endl;
				cout << "  ";
				for(size_t i = 0; i < qlen; i++) {
					cout << (char)qry[i];
				}
				cout << endl;
				cout << "  ";
				for(size_t i = 0; i < qlen; i++) {
					cout << "ACGT"[ref[rir+i]];
				}
				cout << endl;
			}
			assert_lt(r2i, r2.size());
			assert_eq(r2[r2i].off, ri);
			assert_eq((size_t)mms, r2[r2i].edits.size());
			results.expand();
			RefAlignerHit& result = results.back();
			result.clear();
			result.stratum = 0;
			result.cost = ham;
			assert_eq(0, result.edits.size());
			if(mms > 0) {
				ASSERT_ONLY(size_t mmcur = 0);
				const size_t nonSeedEditsSz = nonSeedEdits.size();
				for(size_t i = 0; i < nonSeedEditsSz; i++) {
					assert(nonSeedEdits[i].initialized());
					assert_lt(mmcur, (size_t)mms);
					assert_eq(nonSeedEdits[i].pos, r2[r2i].edits[mmcur].pos);
					result.edits.push_back(nonSeedEdits[i]);
					assert_eq(nonSeedEdits[i].chr, r2[r2i].edits[mmcur].chr);
					ASSERT_ONLY(mmcur++);
				}
				assert_eq(mmcur, r2[r2i].edits.size());
			}
			assert_eq((size_t)mms, result.edits.size());
			ASSERT_ONLY(r2i++);
			result.off = ri;
			if(--numToFind == 0) return;
		}
		assert_leq(duplicates, r2.size());
		assert_geq(r2.size() - duplicates, results.size() - resultsISz);
		return; // no hit
	}

  private:

	// Vectors for holding edit information
	EList<Edit> nonSeedEdits;

};

/**
 * Concrete RefAligner for finding nearby 1-mismatch hits given an
 * anchor hit.
 */
template<typename TStr>
class Seed1RefAligner : public RefAligner<TStr> {

	typedef BTDnaString TDna5Str;
	typedef BTString TCharStr;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:

	Seed1RefAligner(uint32_t seedLen, uint32_t qualMax) :
		RefAligner<TStr>(seedLen, qualMax) { }

	virtual ~Seed1RefAligner() { }

protected:
	/**
	 * This schematic shows the roles played by the begin, qbegin, end,
	 * qend, halfway, slen, qlen, and lim variables:
	 *
	 * seedOnLeft == true:
	 *
	 * |<                   lim                   >|<     qlen       >|
	 *  --------------------------------------------------------------
	 * |                     | slen | qlen-slen |  | slen | qlen-slen |
	 *  --------------------------------------------------------------
	 * ^                     ^                     ^                  ^
	 * begin & qbegin     halfway                qend               end
	 *
	 * seedOnLeft == false:
	 *
	 * |<     qlen       >|<                   lim                   >|
	 *  --------------------------------------------------------------
	 * | qlen-slen | slen |  | qlen-slen | slen |                     |
	 *  --------------------------------------------------------------
	 * ^                  ^                     ^                     ^
	 * begin            qbegin             halfway           qend & end
	 *
	 * Note that, for seeds longer than 32 base-pairs, the seed is
	 * further subdivided into a 32-bit anchor and a seed overhang of
	 * length > 0.
	 */
	void naiveFind(uint32_t numToFind,
				   uint32_t tidx,
				   uint8_t* ref,
				   const TDna5Str& qry,
				   const TCharStr& quals,
				   uint32_t begin,
				   uint32_t end,
				   EList<RefAlignerHit>& results,
				   TSetPairs* pairs,
				   uint32_t aoff,
				   bool seedOnLeft)
	{
		assert_gt(numToFind, 0);
		assert_gt(end, begin);
		const uint32_t qlen = qry.length();
		assert_gt(qlen, 0);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(this->seedLen_, 0);
		const uint32_t slen = min(qlen, this->seedLen_);
		uint32_t qend = end;
		uint32_t qbegin = begin;
		// If the seed is on the left-hand side of the alignment, then
		// leave a gap at the right-hand side of the interval;
		// otherwise, do the opposite
		if(seedOnLeft) {
			// Leave gap on right-hand side of the interval
			qend -= qlen;
		} else {
			// Leave gap on left-hand side of the interval
			qbegin += qlen;
		}
		// lim = number of alignments to try
		const uint32_t lim = qend - qbegin;
		// halfway = position in the reference to start at (and then
		// we work our way out to the right and to the left).
		const uint32_t halfway = qbegin + (lim >> 1);
		// Vectors for holding edit information
		nonSeedEdits.clear();
		bool hi = false;
		for(uint32_t i = 1; i <= lim+1; i++) {
			uint32_t ri;  // leftmost position in candidate alignment
			uint32_t rir; // same, minus begin; for indexing into ref[]
			if(hi) {
				ri = halfway + (i >> 1); rir = ri - begin;
				assert_leq(ri, qend);
			} else {
				ri = halfway - (i >> 1); rir = ri - begin;
				assert_geq(ri, begin);
			}
			hi = !hi;
			// Do the naive comparison
			bool match = true;
			int refc = -1;
			int readc = -1;
			uint32_t mmOff = 0xffffffff;
			int mms = 0;
			int seedMms = 0;
			unsigned int ham = 0;
			nonSeedEdits.clear();
			// Walk through each position of the alignment
			for(uint32_t jj = 0; jj < qlen; jj++) {
				uint32_t j = jj;
				if(!seedOnLeft) {
					// If seed is on the right, scan right-to-left
					j = qlen - jj - 1;
				} else {
					// Go left-to-right
				}
				uint32_t rirj = rir + j;
				if(!seedOnLeft) {
					assert_geq(rir, jj);
					rirj = rir - jj - 1;
				}
#if 0
				// Count Ns in the reference as mismatches
				const int q = (int)qry[j];
				const int r = (int)ref[rirj];
				assert_leq(q, 4);
				assert_leq(r, 4);
				if(q == 4 || r == 4 || q != r) {
#else
				// Disallow alignments that involve an N in the
				// reference
				const int r = (int)ref[rirj];
				if(r & 4) {
					// N in reference; bail on this alignment
					match = false;
					break;
				}
				const int q = (int)qry[j];
				assert_leq(q, 4);
				assert_lt(r, 4);
				if(q != r) {
#endif
					// Mismatch!
					mms++;
					if(mms > 1 && jj < slen) {
						// More than one mismatch in the anchor; reject
						match = false;
						break;
					}
					uint8_t qual = phredcToPhredq(quals[j]);
					ham += mmPenalty(!gNoMaqRound, qual);
					if(ham > this->qualMax_) {
						// Exceeded quality ceiling; reject
						match = false;
						break;
					} else if(jj < slen) {
						// First mismatch in the anchor; remember offset
						// and ref char
						refc = "ACGT"[r];
						readc = "ACGTN"[q];
						mmOff = j;
						seedMms = 1;
					} else {
						// Legal mismatch outside of the anchor; record it
						nonSeedEdits.push_back(Edit(j, "ACGT"[r], "ACGTN"[q], EDIT_TYPE_MM));
						assert_leq(nonSeedEdits.size(), (size_t)mms);
					}
				}
			}
			if(match) {
				assert_leq(seedMms, mms);
				results.expand();
				RefAlignerHit& result = results.back();
				result.clear();
				result.stratum = seedMms;
				result.cost = ham;
				assert_eq(0, result.edits.size());
				if(mms >= 1) {
					// Be careful to add edits in left-to-right order
					// with respect to the read/alignment
					if(seedOnLeft && seedMms) {
						assert_lt(mmOff, qlen);
						assert_in(refc, "ACGT");
						assert_in(readc, "ACGTN");
						result.edits.push_back(Edit(mmOff, refc, readc, EDIT_TYPE_MM));
					}
					const size_t nonSeedEditsSz = nonSeedEdits.size();
					if(nonSeedEditsSz > 0) {
						if(seedOnLeft) {
							for(size_t k = 0; k < nonSeedEditsSz; k++) {
								result.edits.push_back(nonSeedEdits[k]);
							}
						} else {
							for(size_t k = nonSeedEditsSz; k > 0; k--) {
								result.edits.push_back(nonSeedEdits[k-1]);
							}
						}
					}
					if(!seedOnLeft && seedMms) {
						assert_lt(mmOff, qlen);
						assert_in(refc, "ACGT");
						assert_in(readc, "ACGTN");
						result.edits.push_back(Edit(mmOff, refc, readc, EDIT_TYPE_MM));
					}
				}
				assert_eq((size_t)mms, result.edits.size());
				if(seedOnLeft) {
					result.off = ri;
				} else {
					result.off = ri - qlen;
				}
			}
		}
		return;
	}

	/**
	 * This schematic shows the roles played by the begin, qbegin, end,
	 * qend, halfway, slen, qlen, and lim variables:
	 *
	 * seedOnLeft == true:
	 *
	 * |<                   lim                   >|<     qlen       >|
	 *  --------------------------------------------------------------
	 * |                     | slen | qlen-slen |  | slen | qlen-slen |
	 *  --------------------------------------------------------------
	 * ^                     ^                     ^                  ^
	 * begin & qbegin     halfway                qend               end
	 *
	 * seedOnLeft == false:
	 *
	 *             |<                   lim                   >|
	 *  --------------------------------------------------------------
	 * | qlen-slen |         | qlen-slen | slen |              | slen |
	 *  --------------------------------------------------------------
	 * ^           ^                     ^                     ^      ^
	 * begin       qbegin             halfway                qend   end
	 *
	 * Note that, for seeds longer than 32 base-pairs, the seed is
	 * further subdivided into a 32-bit anchor and a seed overhang of
	 * length > 0.
	 */
	virtual void anchor64Find(uint32_t numToFind,
					uint32_t tidx,
					uint8_t* ref,
					const TDna5Str& qry,
					const TCharStr& quals,
					uint32_t begin,
					uint32_t end,
					EList<RefAlignerHit>& results,
					TSetPairs* pairs = NULL,
					uint32_t aoff = 0xffffffff,
					bool seedOnLeft = false)
	{
		assert_gt(numToFind, 0);
		ASSERT_ONLY(const uint32_t resultsISz = results.size());
		ASSERT_ONLY(uint32_t duplicates = 0);
		ASSERT_ONLY(uint32_t r2i = 0);
		const uint32_t qlen = qry.length();
		assert_gt(qlen, 0);
		assert_gt(end, begin);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(this->seedLen_, 0);
		uint32_t slen = min(qlen, this->seedLen_);
#ifndef NDEBUG
		// Get results from the naive matcher for sanity-checking
		EList<RefAlignerHit> r2;
		naiveFind(numToFind, tidx, ref, qry, quals, begin, end, r2,
		          pairs, aoff, seedOnLeft);
#endif
		const uint32_t anchorBitPairs = min<int>(slen, 32);
		const int lhsShift = ((anchorBitPairs - 1) << 1);
		const uint32_t anchorCushion  = 32 - anchorBitPairs;
		// seedAnchorOverhang = # seed bases not included in the anchor
		const uint32_t seedAnchorOverhang = (slen <= anchorBitPairs ? 0 : (slen - anchorBitPairs));
		// seedAnchorOverhang = # seed bases not included in the anchor
		const uint32_t readSeedOverhang = (slen == qlen ? 0 : (qlen - slen));
		assert(anchorCushion == 0 || seedAnchorOverhang == 0);
		assert_eq(qlen, readSeedOverhang + slen);
		uint32_t qend = end;
		uint32_t qbegin = begin;
		if(seedOnLeft) {
			// Leave read-sized gap on right-hand side of the interval
			qend -= qlen;
		} else {
			// Leave seed-sized gap on right-hand side and
			// non-seed-sized gap on the left-hand side
			qbegin += readSeedOverhang;
			qend -= slen;
		}
		// lim = # possible alignments in the range
		const uint32_t lim = qend - qbegin;
		// halfway = point on the genome to radiate out from
		const uint32_t halfway = qbegin + (lim >> 1);
		uint64_t anchor = 0llu;
		uint64_t buffw = 0llu;  // rotating ref sequence buffer
		// OR the 'diff' buffer with this so that we can always count
		// 'N's as mismatches
		uint64_t diffMask = 0llu;
		// Set up a mask that we'll apply to the two bufs every round
		// to discard bits that were rotated out of the anchor area
		uint64_t clearMask = 0xffffffffffffffffllu;
		bool useMask = false;
		if(anchorBitPairs < 32) {
			clearMask >>= ((32-anchorBitPairs) << 1);
			useMask = true;
		}
		int nsInSeed = 0;
		int nPos = -1;
		uint32_t skipLeftToRights = 0;
		uint32_t skipRightToLefts = 0;
		const uint32_t halfwayRi = halfway - begin;
		// Construct the 'anchor' 64-bit buffer so that it holds all of
		// the first 'anchorBitPairs' bit pairs of the query.
		for(uint32_t ii = 0; ii < anchorBitPairs; ii++) {
			uint32_t i = ii;
			if(!seedOnLeft) {
				// Fill in the anchor using characters from the right-
				// hand side of the query (but take the characters in
				// left-to-right order)
				i = qlen - slen + ii;
			}
			int c = (int)qry[i]; // next query character
			int r = (int)ref[halfwayRi + ii]; // next reference character
			if(r & 4) {
				// The reference character is an N; to mimic the
				// behavior of BW alignment, we have to skip all
				// alignments that involve an N in the reference.  Set
				// the skip* variables accordingly.
				r = 0;
				uint32_t lrSkips = ii;
				uint32_t rlSkips = qlen - ii;
				if(!seedOnLeft && readSeedOverhang) {
					lrSkips += readSeedOverhang;
					assert_geq(rlSkips, readSeedOverhang);
					rlSkips -= readSeedOverhang;
				}
				// The right-to-left direction absorbs the candidate
				// alignment based at 'halfway'
				skipLeftToRights = max(skipLeftToRights, lrSkips);
				skipRightToLefts = max(skipRightToLefts, rlSkips);
				assert_leq(skipLeftToRights, qlen);
				assert_leq(skipRightToLefts, qlen);
			}
			assert_leq(c, 4);
			assert_lt(r, 4);
			// Special case: query has an 'N'
			if(c == 4) {
				if(++nsInSeed > 1) {
					// More than one 'N' in the anchor region; can't
					// possibly have a 1-mismatch hit anywhere
					assert_eq(r2.size(), results.size() - resultsISz);
					return;   // can't match if query has Ns
				}
				nPos = (int)ii; // w/r/t LHS of anchor
				// Make it look like an 'A' in the anchor
				c = 0;
				diffMask = (diffMask << 2llu) | 1llu;
			} else {
				diffMask <<= 2llu;
			}
			anchor = ((anchor << 2llu) | c);
			buffw = ((buffw << 2llu) | r);
		}
		// Check whether read is disqualified by Ns inside the seed
		// region but outside the anchor region
		if(seedAnchorOverhang) {
			assert_lt(anchorBitPairs, slen);
			for(uint32_t ii = anchorBitPairs; ii < slen; ii++) {
				uint32_t i = ii;
				if(!seedOnLeft) {
					i = qlen - slen + ii;
				}
				if((int)qry[i] == 4) {
					if(++nsInSeed > 1) {
						assert_eq(r2.size(), results.size() - resultsISz);
						return; // can't match if query has Ns
					}
				}
			}
		} else {
			assert_eq(anchorBitPairs, slen);
		}
		uint64_t bufbw = buffw;
		// Slide the anchor out in either direction, alternating
		// between right-to-left and left-to-right shifts, until all of
		// the positions from qbegin to qend have been covered.
		bool hi = false;
		uint32_t riHi  = halfway;
		uint32_t rirHi = halfway - begin;
		uint32_t rirHiAnchor = rirHi + anchorBitPairs - 1;
		uint32_t riLo  = halfway + 1;
		uint32_t rirLo = halfway - begin + 1;
		uint32_t lrSkips = anchorBitPairs;
		uint32_t rlSkips = qlen;
		if(!seedOnLeft && readSeedOverhang) {
			lrSkips += readSeedOverhang;
			assert_geq(rlSkips, readSeedOverhang);
			rlSkips -= readSeedOverhang;
		}
		for(uint32_t i = 1; i <= lim + 1; i++) {
			int r;       // new reference char
			uint64_t diff;
			assert_leq(skipLeftToRights, qlen);
			assert_leq(skipRightToLefts, qlen);
			if(hi) {
				hi = false;
				// Moving left-to-right
				riHi++;
				rirHi++;
				rirHiAnchor++;
				r = (int)ref[rirHiAnchor];
				if(r & 4) {
					r = 0;
					skipLeftToRights = lrSkips;
				}
				assert_lt(r, 4);
				// Bring in new base pair at the least significant
				// position
				buffw = ((buffw << 2llu) | r);
				if(useMask) buffw &= clearMask;
				if(skipLeftToRights > 0) {
					skipLeftToRights--;
					continue;
				}
				diff = (buffw ^ anchor) | diffMask;
			} else {
				hi = true;
				// Moving right-to-left
				riLo--;
				rirLo--;
				r = (int)ref[rirLo];
				if(r & 4) {
					r = 0;
					skipRightToLefts = rlSkips;
				}
				assert_lt(r, 4);
				if(i >= 2) {
					bufbw >>= 2llu;
					// Bring in new base pair at the most significant
					// position
					bufbw |= ((uint64_t)r << lhsShift);
				}
				if(skipRightToLefts > 0) {
					skipRightToLefts--;
					continue;
				}
				diff = (bufbw ^ anchor) | diffMask;
			}
			if((diff & 0xffffffff00000000llu) &&
			   (diff & 0x00000000ffffffffllu)) continue;
			uint32_t ri  = hi ? riLo  : riHi;
			uint32_t rir = hi ? rirLo : rirHi;
			// Could use pop count
			uint8_t *diff8 = reinterpret_cast<uint8_t*>(&diff);
			// As a first cut, see if there are too many mismatches in
			// the first and last parts of the anchor
			uint32_t diffs = u8toMms[(int)diff8[0]] + u8toMms[(int)diff8[7]];
			if(diffs > 1) continue;
			diffs += u8toMms[(int)diff8[1]] +
					 u8toMms[(int)diff8[2]] +
					 u8toMms[(int)diff8[3]] +
					 u8toMms[(int)diff8[4]] +
					 u8toMms[(int)diff8[5]] +
					 u8toMms[(int)diff8[6]];
			uint32_t mmpos = 0xffffffff;
			int refc = -1;
			int readc = -1;
			unsigned int ham = 0;
			if(diffs > 1) {
				// Too many differences in the seed; stop
				continue;
			} else if(diffs == 1 && nPos != -1) {
				// There was one difference, but there was also one N,
				// so we already know where the difference is
				mmpos = nPos;
				assert_lt(mmpos, anchorBitPairs);
				refc = "ACGT"[(int)ref[rir + nPos]];
				readc = 'N';
				if(!seedOnLeft) {
					mmpos += readSeedOverhang;
				}
				char q = quals[mmpos];
				ham += mmPenalty(!gNoMaqRound, phredcToPhredq(q));
				if(ham > this->qualMax_) {
					// Exceeded quality limit
					continue;
				}
			} else if(diffs == 1) {
				// Figure out which position mismatched
				mmpos = 31;
				if((diff & 0xffffffffllu) == 0) { diff >>= 32llu; mmpos -= 16; }
				assert_neq(0, diff);
				if((diff & 0xffffllu) == 0)     { diff >>= 16llu; mmpos -=  8; }
				assert_neq(0, diff);
				if((diff & 0xffllu) == 0)       { diff >>= 8llu;  mmpos -=  4; }
				assert_neq(0, diff);
				if((diff & 0xfllu) == 0)        { diff >>= 4llu;  mmpos -=  2; }
				assert_neq(0, diff);
				if((diff & 0x3llu) == 0)        { mmpos--; }
				assert_neq(0, diff);
				assert_geq(mmpos, 0);
				assert_lt(mmpos, 32);
				readc = "ACGT"[(anchor >> (62-mmpos*2)) & 3];
				if((diffMask >> (62-mmpos*2)) & 3) {
					readc = 'N';
				}
				mmpos -= anchorCushion;
				assert_lt(mmpos, anchorBitPairs);
				refc = "ACGT"[(int)ref[rir + mmpos]];
				assert_neq(refc, readc);
				if(!seedOnLeft) {
					mmpos += readSeedOverhang;
				}
				char q = quals[mmpos];
				ham += mmPenalty(!gNoMaqRound, phredcToPhredq(q));
				if(ham > this->qualMax_) {
					// Exceeded quality limit
					continue;
				}
			}
			// If the seed is longer than the anchor, then scan the
			// rest of the seed characters
			bool foundHit = true;
			if(seedAnchorOverhang) {
				for(uint32_t j = 0; j < seedAnchorOverhang; j++) {
					int rc = (int)ref[rir + anchorBitPairs + j];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							// Skip out of the seedAnchorOverhang
							assert_eq(0, skipRightToLefts);
							skipRightToLefts = seedAnchorOverhang - j - 1;
							if(seedOnLeft) {
								// ...and skip out of the rest of the read
								skipRightToLefts += readSeedOverhang;
							}
						} else {
							// Left-to-right
							// Skip out of the seedAnchorOverhang
							assert_eq(0, skipLeftToRights);
							skipLeftToRights = anchorBitPairs + j;
							if(!seedOnLeft) {
								// ...and skip out of the rest of the read
								skipLeftToRights += readSeedOverhang;
							}
						}
						foundHit = false; // Skip this candidate
						break;
					}
					uint32_t qoff = anchorBitPairs + j;
					if(!seedOnLeft) {
						qoff += readSeedOverhang;
					}
					assert_lt(qoff, qlen);
					int qc = (int)qry[qoff];
					if(qc != rc) {
						if(++diffs > 1) {
							foundHit = false;
							break;
						} else {
							assert_eq(0xffffffff, mmpos);
							mmpos = qoff;
							assert_eq(-1, refc);
							refc = "ACGT"[rc];
							readc = "ACGTN"[qc];
							char q = phredcToPhredq(quals[qoff]);
							ham += mmPenalty(!gNoMaqRound, q);
							if(ham > this->qualMax_) {
								// Exceeded quality limit
								foundHit = false;
								break;
							}
						}
					}
				}
				if(!foundHit) continue;
			}
			// If the read is longer than the seed, then scan the rest
			// of the read characters; mismatches no longer count
			// toward the stratum or the 1-mm limit.
			// Vectors for holding edit information
			nonSeedEdits.clear();
			int mms = diffs; // start counting total mismatches
			if((qlen - slen) > 0) {
				// Going left-to-right
				for(uint32_t j = 0; j < readSeedOverhang; j++) {
					uint32_t roff = rir + slen + j;
					uint32_t qoff = slen + j;
					if(!seedOnLeft) {
						assert_geq(roff, qlen);
						roff -= qlen;
						qoff = j;
					}
					int rc = (int)ref[roff];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							// Skip what's left of the readSeedOverhang
							skipRightToLefts = readSeedOverhang - j - 1;
							if(!seedOnLeft) {
								// ...and skip the seed if it's on the right
								skipRightToLefts += slen;
							}
						} else {
							// Left-to-right
							// Skip what we've matched of the overhang
							skipLeftToRights = j;
							if(seedOnLeft) {
								// ...and skip the seed if it's on the left
								skipLeftToRights += slen;
							}
						}
						foundHit = false; // Skip this candidate
						break;
					}
					int qc = (int)qry[qoff];
					if(qc != rc) {
						// Calculate quality of mismatched base
						char q = phredcToPhredq(quals[qoff]);
						ham += mmPenalty(!gNoMaqRound, q);
						if(ham > this->qualMax_) {
							// Exceeded quality limit
							foundHit = false;
							break;
						}
						// Legal mismatch outside of the anchor; record it
						mms++;
						nonSeedEdits.push_back(Edit(qoff, "ACGT"[rc], "ACGTN"[qc], EDIT_TYPE_MM));
						assert_leq(nonSeedEdits.size(), (size_t)mms);
					}
				}
				if(!foundHit) continue;
			}
			assert(foundHit);
			// Adjust ri if seed is on the right-hand side
			if(!seedOnLeft) {
				ri -= readSeedOverhang;
				rir -= readSeedOverhang;
			}
			if(pairs != NULL) {
				TU64Pair p;
				if(ri < aoff) {
					// By convention, the upstream mate's
					// coordinates go in the 'first' field
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.second = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				} else {
					p.second = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				}
				if(pairs->find(p) != pairs->end()) {
					// We already found this hit!  Continue.
					ASSERT_ONLY(duplicates++);
					ASSERT_ONLY(r2i++);
					continue;
				} else {
					// Record this hit
					pairs->insert(p);
				}
			}
			if(gVerbose) {
				cout << "About to report:" << endl;
				cout << "  ";
				for(size_t i = 0; i < qlen; i++) {
					cout << (char)qry[i];
				}
				cout << endl;
				cout << "  ";
				for(size_t i = 0; i < qlen; i++) {
					cout << "ACGT"[ref[rir+i]];
				}
				cout << endl;
			}
			assert_leq(diffs, 1);
			assert_geq((size_t)mms, diffs);
			assert_lt(r2i, r2.size());
			assert_eq(r2[r2i].off, ri);
			assert_eq((size_t)mms, r2[r2i].edits.size());
			results.expand();
			RefAlignerHit& result = results.back();
			result.clear();
			result.stratum = diffs;
			result.cost = ham;
			assert_eq(0, result.edits.size());
			if(mms > 0) {
				ASSERT_ONLY(size_t mmcur = 0);
				if(seedOnLeft && diffs > 0) {
					assert_neq(mmpos, 0xffffffff);
					assert_lt(mmpos, qlen);
					assert_lt(mmcur, (size_t)mms);
					assert_eq(mmpos, r2[r2i].edits[mmcur].pos);
					assert_in(refc, "ACGT");
					assert_in(readc, "ACGTN");
					assert_eq((unsigned)refc, r2[r2i].edits[mmcur].chr);
					ASSERT_ONLY(mmcur++);
					result.edits.push_back(Edit(mmpos, refc, readc, EDIT_TYPE_MM));
				}
				const size_t nonSeedEditsSz = nonSeedEdits.size();
				for(size_t i = 0; i < nonSeedEditsSz; i++) {
					assert(nonSeedEdits[i].initialized());
					assert_lt(mmcur, (size_t)mms);
					assert_eq(nonSeedEdits[i].pos, r2[r2i].edits[mmcur].pos);
					result.edits.push_back(nonSeedEdits[i]);
					assert_eq(nonSeedEdits[i].chr, r2[r2i].edits[mmcur].chr);
					ASSERT_ONLY(mmcur++);
				}
				if(!seedOnLeft && diffs > 0) {
					assert_neq(mmpos, 0xffffffff);
					assert_lt(mmpos, qlen);
					assert_lt(mmcur, (size_t)mms);
					assert_eq(mmpos, r2[r2i].edits[mmcur].pos);
					assert_in(refc, "ACGT");
					assert_in(readc, "ACGTN");
					assert_eq((unsigned)refc, r2[r2i].edits[mmcur].chr);
					ASSERT_ONLY(mmcur++);
					result.edits.push_back(Edit(mmpos, refc, readc, EDIT_TYPE_MM));
				}
				assert_eq(mmcur, r2[r2i].edits.size());
			}
			assert_eq((size_t)mms, result.edits.size());
			ASSERT_ONLY(r2i++);
			result.off = ri;
			if(--numToFind == 0) return;
		}
		assert_leq(duplicates, r2.size());
		assert_geq(r2.size() - duplicates, results.size() - resultsISz);
		return; // no hit
	}

  private:

	EList<Edit> nonSeedEdits;
};

/**
 * Concrete RefAligner for finding nearby 2-mismatch hits given an
 * anchor hit.
 */
template<typename TStr>
class Seed2RefAligner : public RefAligner<TStr> {

	typedef BTDnaString TDna5Str;
	typedef BTString TCharStr;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:

	Seed2RefAligner(uint32_t seedLen, uint32_t qualMax) :
		RefAligner<TStr>(seedLen, qualMax) { }

	virtual ~Seed2RefAligner() { }

protected:
	/**
	 * This schematic shows the roles played by the begin, qbegin, end,
	 * qend, halfway, slen, qlen, and lim variables:
	 *
	 * seedOnLeft == true:
	 *
	 * |<                   lim                   >|<     qlen       >|
	 *  --------------------------------------------------------------
	 * |                     | slen | qlen-slen |  | slen | qlen-slen |
	 *  --------------------------------------------------------------
	 * ^                     ^                     ^                  ^
	 * begin & qbegin     halfway                qend               end
	 *
	 * seedOnLeft == false:
	 *
	 * |<     qlen       >|<                   lim                   >|
	 *  --------------------------------------------------------------
	 * | qlen-slen | slen |  | qlen-slen | slen |                     |
	 *  --------------------------------------------------------------
	 * ^                  ^                     ^                     ^
	 * begin            qbegin             halfway           qend & end
	 *
	 * Note that, for seeds longer than 32 base-pairs, the seed is
	 * further subdivided into a 32-bit anchor and a seed overhang of
	 * length > 0.
	 */
	void naiveFind(uint32_t numToFind,
				   uint32_t tidx,
				   uint8_t* ref,
				   const TDna5Str& qry,
				   const TCharStr& quals,
				   uint32_t begin,
				   uint32_t end,
				   EList<RefAlignerHit>& results,
				   TSetPairs* pairs,
				   uint32_t aoff,
				   bool seedOnLeft)
	{
		assert_gt(numToFind, 0);
		assert_gt(end, begin);
		const uint32_t qlen = qry.length();
		assert_gt(qlen, 0);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(this->seedLen_, 0);
		const uint32_t slen = min(qlen, this->seedLen_);
		uint32_t qend = end;
		uint32_t qbegin = begin;
		// If the seed is on the left-hand side of the alignment, then
		// leave a gap at the right-hand side of the interval;
		// otherwise, do the opposite
		if(seedOnLeft) {
			// Leave gap on right-hand side of the interval
			qend -= qlen;
		} else {
			// Leave gap on left-hand side of the interval
			qbegin += qlen;
		}
		// lim = number of alignments to try
		const uint32_t lim = qend - qbegin;
		// halfway = position in the reference to start at (and then
		// we work our way out to the right and to the left).
		const uint32_t halfway = qbegin + (lim >> 1);
		// Vectors for holding edit information
		nonSeedEdits.clear();
		bool hi = false;
		for(uint32_t i = 1; i <= lim+1; i++) {
			uint32_t ri;  // leftmost position in candidate alignment
			uint32_t rir; // same, minus begin; for indexing into ref[]
			if(hi) {
				ri = halfway + (i >> 1); rir = ri - begin;
				assert_leq(ri, qend);
			} else {
				ri = halfway - (i >> 1); rir = ri - begin;
				assert_geq(ri, begin);
			}
			hi = !hi;
			// Do the naive comparison
			bool match = true;
			int refc1 = -1;
			int readc1 = -1;
			uint32_t mmOff1 = 0xffffffff;
			int refc2 = -1;
			int readc2 = -1;
			uint32_t mmOff2 = 0xffffffff;
			int mms = 0;
			int seedMms = 0;
			unsigned int ham = 0;
			nonSeedEdits.clear();
			// Walk through each position of the alignment
			for(uint32_t jj = 0; jj < qlen; jj++) {
				uint32_t j = jj;
				if(!seedOnLeft) {
					// If seed is on the right, scan right-to-left
					j = qlen - jj - 1;
				} else {
					// Go left-to-right
				}
				uint32_t rirj = rir + j;
				if(!seedOnLeft) {
					assert_geq(rir, jj);
					rirj = rir - jj - 1;
				}
#if 0
				// Count Ns in the reference as mismatches
				const int q = (int)qry[j];
				const int r = (int)ref[rirj];
				assert_leq(q, 4);
				assert_leq(r, 4);
				if(q == 4 || r == 4 || q != r) {
#else
				// Disallow alignments that involve an N in the
				// reference
				const int r = (int)ref[rirj];
				if(r & 4) {
					// N in reference; bail on this alignment
					match = false;
					break;
				}
				const int q = (int)qry[j];
				assert_leq(q, 4);
				assert_lt(r, 4);
				if(q != r) {
#endif
					// Mismatch!
					mms++;
					if(mms > 2 && jj < slen) {
						// More than one mismatch in the anchor; reject
						match = false;
						break;
					}
					uint8_t qual = phredcToPhredq(quals[j]);
					ham += mmPenalty(!gNoMaqRound, qual);
					if(ham > this->qualMax_) {
						// Exceeded quality ceiling; reject
						match = false;
						break;
					} else if(mms == 1 && jj < slen) {
						// First mismatch in the anchor; remember offset
						// and ref char
						refc1 = "ACGT"[r];
						readc1 = "ACGTN"[q];
						mmOff1 = j;
						seedMms = 1;
					} else if(mms == 2 && jj < slen) {
						// Second mismatch in the anchor; remember offset
						// and ref char
						refc2 = "ACGT"[r];
						readc2 = "ACGTN"[q];
						mmOff2 = j;
						seedMms = 2;
					} else {
						// Legal mismatch outside of the anchor; record it
						nonSeedEdits.push_back(Edit(j, "ACGT"[r], "ACGTN"[q], EDIT_TYPE_MM));
						assert_leq(nonSeedEdits.size(), (size_t)mms);
					}
				}
			}
			if(match) {
				assert_leq(seedMms, mms);
				results.expand();
				RefAlignerHit& result = results.back();
				result.clear();
				result.stratum = seedMms;
				result.cost = ham;
				assert_eq(0, result.edits.size());
				if(mms >= 1) {
					// Be careful to add edits in left-to-right order
					// with respect to the read/alignment
					if(seedOnLeft && seedMms) {
						assert_lt(mmOff1, qlen);
						assert_in(refc1, "ACGT");
						assert_in(readc1, "ACGTN");
						result.edits.push_back(Edit(mmOff1, refc1, readc1, EDIT_TYPE_MM));
						if(seedMms > 1) {
							assert_lt(mmOff1, mmOff2);
							assert_lt(mmOff2, qlen);
							assert_in(refc2, "ACGT");
							assert_in(readc2, "ACGTN");
							result.edits.push_back(Edit(mmOff2, refc2, readc2, EDIT_TYPE_MM));
						}
					}
					const size_t nonSeedEditsSz = nonSeedEdits.size();
					if(nonSeedEditsSz > 0) {
						if(seedOnLeft) {
							for(size_t k = 0; k < nonSeedEditsSz; k++) {
								result.edits.push_back(nonSeedEdits[k]);
							}
						} else {
							for(size_t k = nonSeedEditsSz; k > 0; k--) {
								result.edits.push_back(nonSeedEdits[k-1]);
							}
						}
					}
					if(!seedOnLeft && seedMms) {
						if(seedMms > 1) {
							assert_lt(mmOff2, mmOff1);
							assert_lt(mmOff2, qlen);
							assert_in(refc2, "ACGT");
							assert_in(readc2, "ACGTN");
							result.edits.push_back(Edit(mmOff2, refc2, readc2, EDIT_TYPE_MM));
						}
						assert_lt(mmOff1, qlen);
						assert_in(refc1, "ACGT");
						assert_in(readc1, "ACGTN");
						result.edits.push_back(Edit(mmOff1, refc1, readc1, EDIT_TYPE_MM));
					}
				}
				assert_eq((size_t)mms, result.edits.size());
				if(seedOnLeft) {
					result.off = ri;
				} else {
					result.off = ri - qlen;
				}
			}
		}
		return;
	}

	/**
	 * This schematic shows the roles played by the begin, qbegin, end,
	 * qend, halfway, slen, qlen, and lim variables:
	 *
	 * seedOnLeft == true:
	 *
	 * |<                   lim                   >|<     qlen       >|
	 *  --------------------------------------------------------------
	 * |                     | slen | qlen-slen |  | slen | qlen-slen |
	 *  --------------------------------------------------------------
	 * ^                     ^                     ^                  ^
	 * begin & qbegin     halfway                qend               end
	 *
	 * seedOnLeft == false:
	 *
	 *             |<                   lim                   >|
	 *  --------------------------------------------------------------
	 * | qlen-slen |         | qlen-slen | slen |              | slen |
	 *  --------------------------------------------------------------
	 * ^           ^                     ^                     ^      ^
	 * begin       qbegin             halfway                qend   end
	 *
	 * Note that, for seeds longer than 32 base-pairs, the seed is
	 * further subdivided into a 32-bit anchor and a seed overhang of
	 * length > 0.
	 */
	virtual void anchor64Find(uint32_t numToFind,
					uint32_t tidx,
					uint8_t* ref,
					const TDna5Str& qry,
					const TCharStr& quals,
					uint32_t begin,
					uint32_t end,
					EList<RefAlignerHit>& results,
					TSetPairs* pairs = NULL,
					uint32_t aoff = 0xffffffff,
					bool seedOnLeft = false)
	{
		assert_gt(numToFind, 0);
		ASSERT_ONLY(const uint32_t resultsISz = results.size());
		ASSERT_ONLY(uint32_t duplicates = 0);
		ASSERT_ONLY(uint32_t r2i = 0);
		const uint32_t qlen = qry.length();
		assert_gt(qlen, 0);
		assert_gt(end, begin);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(this->seedLen_, 0);
		uint32_t slen = min(qlen, this->seedLen_);
#ifndef NDEBUG
		// Get results from the naive matcher for sanity-checking
		EList<RefAlignerHit> r2;
		naiveFind(numToFind, tidx, ref, qry, quals, begin, end, r2,
		          pairs, aoff, seedOnLeft);
#endif
		const uint32_t anchorBitPairs = min<int>(slen, 32);
		const int lhsShift = ((anchorBitPairs - 1) << 1);
		const uint32_t anchorCushion  = 32 - anchorBitPairs;
		// seedAnchorOverhang = # seed bases not included in the anchor
		const uint32_t seedAnchorOverhang = (slen <= anchorBitPairs ? 0 : (slen - anchorBitPairs));
		// seedAnchorOverhang = # seed bases not included in the anchor
		const uint32_t readSeedOverhang = (slen == qlen ? 0 : (qlen - slen));
		assert(anchorCushion == 0 || seedAnchorOverhang == 0);
		assert_eq(qlen, readSeedOverhang + slen);
		uint32_t qend = end;
		uint32_t qbegin = begin;
		if(seedOnLeft) {
			// Leave read-sized gap on right-hand side of the interval
			qend -= qlen;
		} else {
			// Leave seed-sized gap on right-hand side and
			// non-seed-sized gap on the left-hand side
			qbegin += readSeedOverhang;
			qend -= slen;
		}
		// lim = # possible alignments in the range
		const uint32_t lim = qend - qbegin;
		// halfway = point on the genome to radiate out from
		const uint32_t halfway = qbegin + (lim >> 1);
		uint64_t anchor = 0llu;
		uint64_t buffw = 0llu;  // rotating ref sequence buffer
		// OR the 'diff' buffer with this so that we can always count
		// 'N's as mismatches
		uint64_t diffMask = 0llu;
		// Set up a mask that we'll apply to the two bufs every round
		// to discard bits that were rotated out of the anchor area
		uint64_t clearMask = 0xffffffffffffffffllu;
		bool useMask = false;
		if(anchorBitPairs < 32) {
			clearMask >>= ((32-anchorBitPairs) << 1);
			useMask = true;
		}
		int nsInSeed = 0;
		uint32_t nPoss = 0;
		int nPos1 = -1;
		int nPos2 = -1;
		uint32_t skipLeftToRights = 0;
		uint32_t skipRightToLefts = 0;
		const uint32_t halfwayRi = halfway - begin;
		assert_leq(anchorBitPairs, slen);
		// Construct the 'anchor' 64-bit buffer so that it holds all of
		// the first 'anchorBitPairs' bit pairs of the query.
		for(uint32_t ii = 0; ii < anchorBitPairs; ii++) {
			uint32_t i = ii;
			if(!seedOnLeft) {
				// Fill in the anchor using characters from the seed
				// portion of the read, starting at the left.  Note
				// that we're subtracting by slen rather than
				// anchorBitPairs because we want the seed anchor
				// overhang to be on the right-hand side
				i = qlen - slen + ii;
			}
			int c = (int)qry[i]; // next query character
			int r = (int)ref[halfwayRi + ii]; // next reference character
			if(r & 4) {
				// The reference character is an N; to mimic the
				// behavior of BW alignment, we have to skip all
				// alignments that involve an N in the reference.  Set
				// the skip* variables accordingly.
				r = 0;
				uint32_t lrSkips = ii;
				uint32_t rlSkips = qlen - ii;
				if(!seedOnLeft && readSeedOverhang) {
					lrSkips += readSeedOverhang;
					assert_geq(rlSkips, readSeedOverhang);
					rlSkips -= readSeedOverhang;
				}
				// The right-to-left direction absorbs the candidate
				// alignment based at 'halfway'
				skipLeftToRights = max(skipLeftToRights, lrSkips);
				skipRightToLefts = max(skipRightToLefts, rlSkips);
				assert_leq(skipLeftToRights, qlen);
				assert_leq(skipRightToLefts, qlen);
			}
			assert_leq(c, 4);
			assert_lt(r, 4);
			// Special case: query has an 'N'
			if(c == 4) {
				if(++nsInSeed > 2) {
					// More than one 'N' in the anchor region; can't
					// possibly have a 1-mismatch hit anywhere
					assert_eq(r2.size(), results.size() - resultsISz);
					return;   // can't match if query has Ns
				}
				if(nsInSeed == 1) {
					nPos1 = (int)ii; // w/r/t LHS of anchor
				} else {
					assert_eq(2, nsInSeed);
					nPos2 = (int)ii; // w/r/t LHS of anchor
					assert_gt(nPos2, nPos1);
				}
				// Make it look like an 'A' in the anchor
				c = 0;
				diffMask = (diffMask << 2llu) | 1llu;
			} else {
				diffMask <<= 2llu;
			}
			anchor = ((anchor << 2llu) | c);
			buffw = ((buffw << 2llu) | r);
		}
		// Check whether read is disqualified by Ns inside the seed
		// region but outside the anchor region
		if(seedAnchorOverhang) {
			assert_lt(anchorBitPairs, slen);
			for(uint32_t ii = anchorBitPairs; ii < slen; ii++) {
				uint32_t i = ii;
				if(!seedOnLeft) {
					i = qlen - slen + ii;
				}
				if((int)qry[i] == 4) {
					if(++nsInSeed > 2) {
						assert_eq(r2.size(), results.size() - resultsISz);
						return; // can't match if query has Ns
					}
				}
			}
		} else {
			assert_eq(anchorBitPairs, slen);
		}
		uint64_t bufbw = buffw;
		// Slide the anchor out in either direction, alternating
		// between right-to-left and left-to-right shifts, until all of
		// the positions from qbegin to qend have been covered.
		bool hi = false;
		uint32_t riHi  = halfway;
		uint32_t rirHi = halfway - begin;
		uint32_t rirHiAnchor = rirHi + anchorBitPairs - 1;
		uint32_t riLo  = halfway + 1;
		uint32_t rirLo = halfway - begin + 1;
		uint32_t lrSkips = anchorBitPairs;
		uint32_t rlSkips = qlen;
		if(!seedOnLeft && readSeedOverhang) {
			lrSkips += readSeedOverhang;
			assert_geq(rlSkips, readSeedOverhang);
			rlSkips -= readSeedOverhang;
		}
		for(uint32_t i = 1; i <= lim + 1; i++) {
			int r;       // new reference char
			uint64_t diff;
			assert_leq(skipLeftToRights, qlen);
			assert_leq(skipRightToLefts, qlen);
			if(hi) {
				hi = false;
				// Moving left-to-right
				riHi++;
				rirHi++;
				rirHiAnchor++;
				r = (int)ref[rirHiAnchor];
				if(r & 4) {
					r = 0;
					skipLeftToRights = lrSkips;
				}
				assert_lt(r, 4);
				// Bring in new base pair at the least significant
				// position
				buffw = ((buffw << 2llu) | r);
				if(useMask) buffw &= clearMask;
				if(skipLeftToRights > 0) {
					skipLeftToRights--;
					continue;
				}
				diff = (buffw ^ anchor) | diffMask;
			} else {
				hi = true;
				// Moving right-to-left
				riLo--;
				rirLo--;
				r = (int)ref[rirLo];
				if(r & 4) {
					r = 0;
					skipRightToLefts = rlSkips;
				}
				assert_lt(r, 4);
				if(i >= 2) {
					bufbw >>= 2llu;
					// Bring in new base pair at the most significant
					// position
					bufbw |= ((uint64_t)r << lhsShift);
				}
				if(skipRightToLefts > 0) {
					skipRightToLefts--;
					continue;
				}
				diff = (bufbw ^ anchor) | diffMask;
			}
			if((diff & 0xf00f00f00f00f00fllu) &&
			   (diff & 0x0f00f00f00f00f00llu) &&
			   (diff & 0x00f00f00f00f00f0llu)) continue;
			if((diff & 0xc30c30c30c30c30cllu) &&
			   (diff & 0x30c30c30c30c30c3llu) &&
			   (diff & 0x0c30c30c30c30c30llu)) continue;
			uint32_t ri  = hi ? riLo  : riHi;
			uint32_t rir = hi ? rirLo : rirHi;
			// Could use pop count
			uint8_t *diff8 = reinterpret_cast<uint8_t*>(&diff);
			// As a first cut, see if there are too many mismatches in
			// the first and last parts of the anchor
			uint32_t diffs = u8toMms[(int)diff8[0]] + u8toMms[(int)diff8[7]];
			if(diffs > 2) continue;
			diffs += u8toMms[(int)diff8[1]] +
					 u8toMms[(int)diff8[2]] +
					 u8toMms[(int)diff8[3]] +
					 u8toMms[(int)diff8[4]] +
					 u8toMms[(int)diff8[5]] +
					 u8toMms[(int)diff8[6]];
			uint32_t mmpos1 = 0xffffffff;
			int refc1 = -1;
			int readc1 = -1;
			uint32_t mmpos2 = 0xffffffff;
			int refc2 = -1;
			int readc2 = -1;
			unsigned int ham = 0;
			if(diffs > 2) {
				// Too many differences in the seed; stop
				continue;
			} else if(nPoss > 1 && diffs == nPoss) {
				// There was one difference, but there was also one N,
				// so we already know where the difference is
				mmpos1 = nPos1;
				refc1 = "ACGT"[(int)ref[rir + nPos1]];
				readc1 = 'N';
				if(!seedOnLeft) {
					mmpos1 += readSeedOverhang;
				}
				char q = quals[mmpos1];
				ham += mmPenalty(!gNoMaqRound, phredcToPhredq(q));
				if(ham > this->qualMax_) {
					// Exceeded quality limit
					continue;
				}
				if(nPoss == 2) {
					mmpos2 = nPos2;
					refc2 = "ACGT"[(int)ref[rir + nPos2]];
					readc2 = 'N';
					if(!seedOnLeft) {
						mmpos2 += readSeedOverhang;
					}
					q = quals[mmpos2];
					ham += mmPenalty(!gNoMaqRound, phredcToPhredq(q));
					if(ham > this->qualMax_) {
						// Exceeded quality limit
						continue;
					}
				}

			} else if(diffs > 0) {
				// Figure out which position mismatched
				uint64_t diff2 = diff;
				mmpos1 = 31;
				if((diff & 0xffffffffllu) == 0) { diff >>= 32llu; mmpos1 -= 16; }
				assert_neq(0, diff);
				if((diff & 0xffffllu) == 0)     { diff >>= 16llu; mmpos1 -=  8; }
				assert_neq(0, diff);
				if((diff & 0xffllu) == 0)       { diff >>= 8llu;  mmpos1 -=  4; }
				assert_neq(0, diff);
				if((diff & 0xfllu) == 0)        { diff >>= 4llu;  mmpos1 -=  2; }
				assert_neq(0, diff);
				if((diff & 0x3llu) == 0)        { mmpos1--; }
				assert_neq(0, diff);
				assert_geq(mmpos1, 0);
				assert_lt(mmpos1, 32);
				readc1 = "ACGT"[(anchor >> (62-mmpos1*2)) & 3];
				if((diffMask >> (62-mmpos1*2)) & 3) {
					readc1 = 'N';
				}
				uint32_t savedMmpos1 = mmpos1;
				mmpos1 -= anchorCushion;
				assert_lt(mmpos1, anchorBitPairs);
				refc1 = "ACGT"[(int)ref[rir + mmpos1]];
				assert_neq(readc1, refc1);
				if(!seedOnLeft) {
					mmpos1 += readSeedOverhang;
				}
				char q = quals[mmpos1];
				ham += mmPenalty(!gNoMaqRound, phredcToPhredq(q));
				if(ham > this->qualMax_) {
					// Exceeded quality limit
					continue;
				}
				if(diffs > 1) {
					// Figure out the second mismatched position
					ASSERT_ONLY(uint64_t origDiff2 = diff2);
					diff2 &= ~(0xc000000000000000llu >> (uint64_t)((savedMmpos1) << 1));
					assert_neq(diff2, origDiff2);
					mmpos2 = 31;
					if((diff2 & 0xffffffffllu) == 0) { diff2 >>= 32llu; mmpos2 -= 16; }
					assert_neq(0, diff2);
					if((diff2 & 0xffffllu) == 0)     { diff2 >>= 16llu; mmpos2 -=  8; }
					assert_neq(0, diff2);
					if((diff2 & 0xffllu) == 0)       { diff2 >>= 8llu;  mmpos2 -=  4; }
					assert_neq(0, diff2);
					if((diff2 & 0xfllu) == 0)        { diff2 >>= 4llu;  mmpos2 -=  2; }
					assert_neq(0, diff2);
					if((diff2 & 0x3llu) == 0)        { mmpos2--; }
					assert_neq(0, diff2);
					assert_geq(mmpos2, 0);
					assert_lt(mmpos2, 32);
					readc2 = "ACGT"[(anchor >> (62-mmpos2*2)) & 3];
					if((diffMask >> (62-mmpos2*2)) & 3) {
						readc2 = 'N';
					}
					mmpos2 -= anchorCushion;
					assert_neq(mmpos1, mmpos2);
					refc2 = "ACGT"[(int)ref[rir + mmpos2]];
					assert_neq(readc2, refc2);
					if(!seedOnLeft) {
						mmpos2 += readSeedOverhang;
					}
					q = quals[mmpos2];
					ham += mmPenalty(!gNoMaqRound, phredcToPhredq(q));
					if(ham > this->qualMax_) {
						// Exceeded quality limit
						continue;
					}
					if(mmpos2 < mmpos1) {
						uint32_t mmtmp = mmpos1;
						mmpos1 = mmpos2;
						mmpos2 = mmtmp;
						int refctmp = refc1;
						refc1 = refc2;
						refc2 = refctmp;
						int readctmp = readc1;
						readc1 = readc2;
						readc2 = readctmp;
					}
					assert_lt(mmpos1, mmpos2);
				}
			}
			// If the seed is longer than the anchor, then scan the
			// rest of the seed characters
			bool foundHit = true;
			if(seedAnchorOverhang) {
				for(uint32_t j = 0; j < seedAnchorOverhang; j++) {
					int rc = (int)ref[rir + anchorBitPairs + j];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							// Skip out of the seedAnchorOverhang
							assert_eq(0, skipRightToLefts);
							skipRightToLefts = seedAnchorOverhang - j - 1;
							if(seedOnLeft) {
								// ...and skip out of the rest of the read
								skipRightToLefts += readSeedOverhang;
							}
						} else {
							// Left-to-right
							// Skip out of the seedAnchorOverhang
							assert_eq(0, skipLeftToRights);
							skipLeftToRights = anchorBitPairs + j;
							if(!seedOnLeft) {
								// ...and skip out of the rest of the read
								skipLeftToRights += readSeedOverhang;
							}
						}
						foundHit = false; // Skip this candidate
						break;
					}
					uint32_t qoff = anchorBitPairs + j;
					if(!seedOnLeft) {
						qoff += readSeedOverhang;
					}
					assert_lt(qoff, qlen);
					int qc = (int)qry[qoff];
					if(qc != rc) {
						diffs++;
						if(diffs > 2) {
							foundHit = false;
							break;
						} else if(diffs == 2) {
							assert_eq(0xffffffff, mmpos2);
							mmpos2 = qoff;
							assert_eq(-1, refc2);
							refc2 = "ACGT"[(int)ref[rir + anchorBitPairs + j]];
							readc2 = "ACGTN"[qc];
						} else {
							assert_eq(1, diffs);
							assert_eq(0xffffffff, mmpos1);
							mmpos1 = qoff;
							assert_eq(-1, refc1);
							refc1 = "ACGT"[(int)ref[rir + anchorBitPairs + j]];
							readc1 = "ACGTN"[qc];
						}
						char q = phredcToPhredq(quals[qoff]);
						ham += mmPenalty(!gNoMaqRound, q);
						if(ham > this->qualMax_) {
							// Exceeded quality limit
							foundHit = false;
							break;
						}
					}
				}
				if(!foundHit) continue;
			}
			// If the read is longer than the seed, then scan the rest
			// of the read characters; mismatches no longer count
			// toward the stratum or the 1-mm limit.
			// Vectors for holding edit information
			nonSeedEdits.clear();
			int mms = diffs; // start counting total mismatches
			if((qlen - slen) > 0) {
				// Going left-to-right
				for(uint32_t j = 0; j < readSeedOverhang; j++) {
					uint32_t roff = rir + slen + j;
					uint32_t qoff = slen + j;
					if(!seedOnLeft) {
						assert_geq(roff, qlen);
						roff -= qlen;
						qoff = j;
					}
					int rc = (int)ref[roff];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							// Skip what's left of the readSeedOverhang
							skipRightToLefts = readSeedOverhang - j - 1;
							if(!seedOnLeft) {
								// ...and skip the seed if it's on the right
								skipRightToLefts += slen;
							}
						} else {
							// Left-to-right
							// Skip what we've matched of the overhang
							skipLeftToRights = j;
							if(seedOnLeft) {
								// ...and skip the seed if it's on the left
								skipLeftToRights += slen;
							}
						}
						foundHit = false; // Skip this candidate
						break;
					}
					int qc = (int)qry[qoff];
					if(qc != rc) {
						// Calculate quality of mismatched base
						char q = phredcToPhredq(quals[qoff]);
						ham += mmPenalty(!gNoMaqRound, q);
						if(ham > this->qualMax_) {
							// Exceeded quality limit
							foundHit = false;
							break;
						}
						// Legal mismatch outside of the anchor; record it
						mms++;
						nonSeedEdits.push_back(Edit(qoff, "ACGT"[rc], "ACGTN"[qc], EDIT_TYPE_MM));
						assert_leq(nonSeedEdits.size(), (size_t)mms);
					}
				}
				if(!foundHit) continue;
			}
			assert(foundHit);
			// Adjust ri if seed is on the right-hand side
			if(!seedOnLeft) {
				ri -= readSeedOverhang;
				rir -= readSeedOverhang;
			}
			if(pairs != NULL) {
				TU64Pair p;
				if(ri < aoff) {
					// By convention, the upstream mate's
					// coordinates go in the 'first' field
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.second = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				} else {
					p.second = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				}
				if(pairs->find(p) != pairs->end()) {
					// We already found this hit!  Continue.
					ASSERT_ONLY(duplicates++);
					ASSERT_ONLY(r2i++);
					continue;
				} else {
					// Record this hit
					pairs->insert(p);
				}
			}
			if(gVerbose) {
				cout << "About to report:" << endl;
				cout << "  ";
				for(size_t i = 0; i < qlen; i++) {
					cout << (char)qry[i];
				}
				cout << endl;
				cout << "  ";
				for(size_t i = 0; i < qlen; i++) {
					cout << "ACGT"[ref[rir+i]];
				}
				cout << endl;
			}
			assert_leq(diffs, 2);
			assert_geq((size_t)mms, diffs);
			assert_lt(r2i, r2.size());
			assert_eq(r2[r2i].off, ri);
			results.expand();
			RefAlignerHit& result = results.back();
			result.clear();
			assert_eq((size_t)mms, r2[r2i].edits.size());
			result.stratum = diffs;
			result.cost = ham;
			if(mms > 0) {
				ASSERT_ONLY(size_t mmcur = 0);
				if(seedOnLeft && diffs > 0) {
					assert_neq(mmpos1, 0xffffffff);
					assert_lt(mmpos1, qlen);
					assert_lt(mmcur, (size_t)mms);
					assert_eq(mmpos1, r2[r2i].edits[mmcur].pos);
					assert_in(refc1, "ACGT");
					assert_in(readc1, "ACGTN");
					assert_eq((unsigned)refc1, r2[r2i].edits[mmcur].chr);
					ASSERT_ONLY(mmcur++);
					result.edits.push_back(Edit(mmpos1, refc1, readc1, EDIT_TYPE_MM));
					if(diffs > 1) {
						assert_eq(2, diffs);
						assert_neq(mmpos2, 0xffffffff);
						assert_lt(mmpos2, qlen);
						assert_lt(mmcur, (size_t)mms);
						assert_eq(mmpos2, r2[r2i].edits[mmcur].pos);
						assert_in(refc2, "ACGT");
						assert_in(readc2, "ACGTN");
						assert_eq((unsigned)refc2, r2[r2i].edits[mmcur].chr);
						ASSERT_ONLY(mmcur++);
						result.edits.push_back(Edit(mmpos2, refc2, readc2, EDIT_TYPE_MM));
					}
				}
				const size_t nonSeedEditsSz = nonSeedEdits.size();
				for(size_t i = 0; i < nonSeedEditsSz; i++) {
					assert(nonSeedEdits[i].initialized());
					assert_lt(mmcur, (size_t)mms);
					assert_eq(nonSeedEdits[i].pos, r2[r2i].edits[mmcur].pos);
					result.edits.push_back(nonSeedEdits[i]);
					assert_eq(nonSeedEdits[i].chr, r2[r2i].edits[mmcur].chr);
					ASSERT_ONLY(mmcur++);
				}
				if(!seedOnLeft && diffs > 0) {
					assert_neq(mmpos1, 0xffffffff);
					assert_lt(mmpos1, qlen);
					assert_lt(mmcur, (size_t)mms);
					assert_eq(mmpos1, r2[r2i].edits[mmcur].pos);
					assert_in(refc1, "ACGT");
					assert_in(readc1, "ACGTN");
					assert_eq((unsigned)refc1, r2[r2i].edits[mmcur].chr);
					ASSERT_ONLY(mmcur++);
					result.edits.push_back(Edit(mmpos1, refc1, readc1, EDIT_TYPE_MM));
					if(diffs > 1) {
						assert_eq(2, diffs);
						assert_neq(mmpos2, 0xffffffff);
						assert_lt(mmpos2, qlen);
						assert_lt(mmcur, (size_t)mms);
						assert_eq(mmpos2, r2[r2i].edits[mmcur].pos);
						assert_in(refc2, "ACGT");
						assert_in(readc2, "ACGTN");
						assert_eq((unsigned)refc2, r2[r2i].edits[mmcur].chr);
						ASSERT_ONLY(mmcur++);
						result.edits.push_back(Edit(mmpos2, refc2, readc2, EDIT_TYPE_MM));
					}
				}
				assert_eq(mmcur, r2[r2i].edits.size());
			}
			assert_eq((size_t)mms, result.edits.size());
			ASSERT_ONLY(r2i++);
			result.off = ri;
			if(--numToFind == 0) return;
		}
		assert_leq(duplicates, r2.size());
		assert_geq(r2.size() - duplicates, results.size() - resultsISz);
		return; // no hit
	}

  private:

	EList<Edit> nonSeedEdits;
};

/**
 * Concrete RefAligner for finding nearby 3-mismatch hits given an
 * anchor hit.
 */
template<typename TStr>
class Seed3RefAligner : public RefAligner<TStr> {

	typedef BTDnaString TDna5Str;
	typedef BTString TCharStr;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:

	Seed3RefAligner(uint32_t seedLen, uint32_t qualMax) :
		RefAligner<TStr>(seedLen, qualMax) { }

	virtual ~Seed3RefAligner() { }

protected:
	/**
	 * This schematic shows the roles played by the begin, qbegin, end,
	 * qend, halfway, slen, qlen, and lim variables:
	 *
	 * seedOnLeft == true:
	 *
	 * |<                   lim                   >|<     qlen       >|
	 *  --------------------------------------------------------------
	 * |                     | slen | qlen-slen |  | slen | qlen-slen |
	 *  --------------------------------------------------------------
	 * ^                     ^                     ^                  ^
	 * begin & qbegin     halfway                qend               end
	 *
	 * seedOnLeft == false:
	 *
	 * |<     qlen       >|<                   lim                   >|
	 *  --------------------------------------------------------------
	 * | qlen-slen | slen |  | qlen-slen | slen |                     |
	 *  --------------------------------------------------------------
	 * ^                  ^                     ^                     ^
	 * begin            qbegin             halfway           qend & end
	 *
	 * Note that, for seeds longer than 32 base-pairs, the seed is
	 * further subdivided into a 32-bit anchor and a seed overhang of
	 * length > 0.
	 */
	void naiveFind(uint32_t numToFind,
				   uint32_t tidx,
				   uint8_t* ref,
				   const TDna5Str& qry,
				   const TCharStr& quals,
				   uint32_t begin,
				   uint32_t end,
				   EList<RefAlignerHit>& results,
				   TSetPairs* pairs,
				   uint32_t aoff,
				   bool seedOnLeft)
	{
		assert_gt(numToFind, 0);
		assert_gt(end, begin);
		const uint32_t qlen = qry.length();
		assert_gt(qlen, 0);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(this->seedLen_, 0);
		const uint32_t slen = min(qlen, this->seedLen_);
		uint32_t qend = end;
		uint32_t qbegin = begin;
		// If the seed is on the left-hand side of the alignment, then
		// leave a gap at the right-hand side of the interval;
		// otherwise, do the opposite
		if(seedOnLeft) {
			// Leave gap on right-hand side of the interval
			qend -= qlen;
		} else {
			// Leave gap on left-hand side of the interval
			qbegin += qlen;
		}
		// lim = number of alignments to try
		const uint32_t lim = qend - qbegin;
		// halfway = position in the reference to start at (and then
		// we work our way out to the right and to the left).
		const uint32_t halfway = qbegin + (lim >> 1);
		// Vectors for holding edit information
		nonSeedEdits.clear();
		bool hi = false;
		for(uint32_t i = 1; i <= lim+1; i++) {
			uint32_t ri;  // leftmost position in candidate alignment
			uint32_t rir; // same, minus begin; for indexing into ref[]
			if(hi) {
				ri = halfway + (i >> 1); rir = ri - begin;
				assert_leq(ri, qend);
			} else {
				ri = halfway - (i >> 1); rir = ri - begin;
				assert_geq(ri, begin);
			}
			hi = !hi;
			// Do the naive comparison
			bool match = true;
			int refc1 = -1;
			int readc1 = -1;
			uint32_t mmOff1 = 0xffffffff;
			int refc2 = -1;
			int readc2 = -1;
			uint32_t mmOff2 = 0xffffffff;
			int refc3 = -1;
			int readc3 = -1;
			uint32_t mmOff3 = 0xffffffff;
			int mms = 0;
			int seedMms = 0;
			unsigned int ham = 0;
			nonSeedEdits.clear();
			// Walk through each position of the alignment
			for(uint32_t jj = 0; jj < qlen; jj++) {
				uint32_t j = jj;
				if(!seedOnLeft) {
					// If seed is on the right, scan right-to-left
					j = qlen - jj - 1;
				} else {
					// Go left-to-right
				}
				uint32_t rirj = rir + j;
				if(!seedOnLeft) {
					assert_geq(rir, jj);
					rirj = rir - jj - 1;
				}
#if 0
				// Count Ns in the reference as mismatches
				const int q = (int)qry[j];
				const int r = (int)ref[rirj];
				assert_leq(q, 4);
				assert_leq(r, 4);
				if(q == 4 || r == 4 || q != r) {
#else
				// Disallow alignments that involve an N in the
				// reference
				const int r = (int)ref[rirj];
				if(r & 4) {
					// N in reference; bail on this alignment
					match = false;
					break;
				}
				const int q = (int)qry[j];
				assert_leq(q, 4);
				assert_lt(r, 4);
				if(q != r) {
#endif
					// Mismatch!
					mms++;
					if(mms > 3 && jj < slen) {
						// More than one mismatch in the anchor; reject
						match = false;
						break;
					}
					uint8_t qual = phredcToPhredq(quals[j]);
					ham += mmPenalty(!gNoMaqRound, qual);
					if(ham > this->qualMax_) {
						// Exceeded quality ceiling; reject
						match = false;
						break;
					} else if(mms == 1 && jj < slen) {
						// First mismatch in the anchor; remember offset
						// and ref char
						refc1 = "ACGT"[r];
						readc1 = "ACGTN"[q];
						mmOff1 = j;
						seedMms = 1;
					} else if(mms == 2 && jj < slen) {
						// Second mismatch in the anchor; remember offset
						// and ref char
						refc2 = "ACGT"[r];
						readc2 = "ACGTN"[q];
						mmOff2 = j;
						seedMms = 2;
					} else if(mms == 3 && jj < slen) {
						// Third mismatch in the anchor; remember offset
						// and ref char
						refc3 = "ACGT"[r];
						readc3 = "ACGTN"[q];
						mmOff3 = j;
						seedMms = 3;
					} else {
						// Legal mismatch outside of the anchor; record it
						nonSeedEdits.push_back(Edit(j, "ACGT"[r], "ACGTN"[q], EDIT_TYPE_MM));
						assert_leq(nonSeedEdits.size(), (size_t)mms);
					}
				}
			}
			if(match) {
				results.expand();
				RefAlignerHit& result = results.back();
				result.clear();
				assert_leq(seedMms, mms);
				result.stratum = seedMms;
				result.cost = ham;
				if(mms >= 1) {
					// Be careful to add edits in left-to-right order
					// with respect to the read/alignment
					if(seedOnLeft && seedMms) {
						assert_lt(mmOff1, qlen);
						assert_in(refc1, "ACGT");
						assert_in(readc1, "ACGTN");
						result.edits.push_back(Edit(mmOff1, refc1, readc1, EDIT_TYPE_MM));
						if(seedMms > 1) {
							assert_lt(mmOff1, mmOff2);
							assert_lt(mmOff2, qlen);
							assert_in(refc2, "ACGT");
							assert_in(readc2, "ACGTN");
							result.edits.push_back(Edit(mmOff2, refc2, readc2, EDIT_TYPE_MM));
							if(seedMms > 2) {
								assert_lt(mmOff2, mmOff3);
								assert_lt(mmOff3, qlen);
								assert_in(refc3, "ACGT");
								assert_in(readc3, "ACGTN");
								result.edits.push_back(Edit(mmOff3, refc3, readc3, EDIT_TYPE_MM));
							}
						}
					}
					const size_t nonSeedEditsSz = nonSeedEdits.size();
					if(nonSeedEditsSz > 0) {
						if(seedOnLeft) {
							for(size_t k = 0; k < nonSeedEditsSz; k++) {
								result.edits.push_back(nonSeedEdits[k]);
							}
						} else {
							for(size_t k = nonSeedEditsSz; k > 0; k--) {
								result.edits.push_back(nonSeedEdits[k-1]);
							}
						}
					}
					if(!seedOnLeft && seedMms) {
						if(seedMms > 1) {
							if(seedMms > 2) {
								assert_lt(mmOff3, mmOff2);
								assert_lt(mmOff3, qlen);
								assert_in(refc3, "ACGT");
								assert_in(readc3, "ACGTN");
								result.edits.push_back(Edit(mmOff3, refc3, readc3, EDIT_TYPE_MM));
							}
							assert_lt(mmOff2, mmOff1);
							assert_lt(mmOff2, qlen);
							assert_in(refc2, "ACGT");
							assert_in(readc2, "ACGTN");
							result.edits.push_back(Edit(mmOff2, refc2, readc2, EDIT_TYPE_MM));
						}
						assert_lt(mmOff1, qlen);
						assert_in(refc1, "ACGT");
						assert_in(readc1, "ACGTN");
						result.edits.push_back(Edit(mmOff1, refc1, readc1, EDIT_TYPE_MM));
					}
				}
				assert_eq((size_t)mms, result.edits.size());
				if(seedOnLeft) {
					result.off = ri;
				} else {
					result.off = ri - qlen;
				}
			}
		}
		return;
	}

	/**
	 * This schematic shows the roles played by the begin, qbegin, end,
	 * qend, halfway, slen, qlen, and lim variables:
	 *
	 * seedOnLeft == true:
	 *
	 * |<                   lim                   >|<     qlen       >|
	 *  --------------------------------------------------------------
	 * |                     | slen | qlen-slen |  | slen | qlen-slen |
	 *  --------------------------------------------------------------
	 * ^                     ^                     ^                  ^
	 * begin & qbegin     halfway                qend               end
	 *
	 * seedOnLeft == false:
	 *
	 *             |<                   lim                   >|
	 *  --------------------------------------------------------------
	 * | qlen-slen |         | qlen-slen | slen |              | slen |
	 *  --------------------------------------------------------------
	 * ^           ^                     ^                     ^      ^
	 * begin       qbegin             halfway                qend   end
	 *
	 * Note that, for seeds longer than 32 base-pairs, the seed is
	 * further subdivided into a 32-bit anchor and a seed overhang of
	 * length > 0.
	 */
	virtual void anchor64Find(uint32_t numToFind,
					uint32_t tidx,
					uint8_t* ref,
					const TDna5Str& qry,
					const TCharStr& quals,
					uint32_t begin,
					uint32_t end,
					EList<RefAlignerHit>& results,
					TSetPairs* pairs = NULL,
					uint32_t aoff = 0xffffffff,
					bool seedOnLeft = false)
	{
		assert_gt(numToFind, 0);
		ASSERT_ONLY(const uint32_t resultsISz = results.size());
		ASSERT_ONLY(uint32_t duplicates = 0);
		ASSERT_ONLY(uint32_t r2i = 0);
		const uint32_t qlen = qry.length();
		assert_gt(qlen, 0);
		assert_gt(end, begin);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(this->seedLen_, 0);
		uint32_t slen = min(qlen, this->seedLen_);
#ifndef NDEBUG
		// Get results from the naive matcher for sanity-checking
		EList<RefAlignerHit> r2;
		naiveFind(numToFind, tidx, ref, qry, quals, begin, end, r2,
		          pairs, aoff, seedOnLeft);
#endif
		const uint32_t anchorBitPairs = min<int>(slen, 32);
		const int lhsShift = ((anchorBitPairs - 1) << 1);
		const uint32_t anchorCushion  = 32 - anchorBitPairs;
		// seedAnchorOverhang = # seed bases not included in the anchor
		const uint32_t seedAnchorOverhang = (slen <= anchorBitPairs ? 0 : (slen - anchorBitPairs));
		// seedAnchorOverhang = # seed bases not included in the anchor
		const uint32_t readSeedOverhang = (slen == qlen ? 0 : (qlen - slen));
		assert(anchorCushion == 0 || seedAnchorOverhang == 0);
		assert_eq(qlen, readSeedOverhang + slen);
		uint32_t qend = end;
		uint32_t qbegin = begin;
		if(seedOnLeft) {
			// Leave read-sized gap on right-hand side of the interval
			qend -= qlen;
		} else {
			// Leave seed-sized gap on right-hand side and
			// non-seed-sized gap on the left-hand side
			qbegin += readSeedOverhang;
			qend -= slen;
		}
		// lim = # possible alignments in the range
		const uint32_t lim = qend - qbegin;
		// halfway = point on the genome to radiate out from
		const uint32_t halfway = qbegin + (lim >> 1);
		uint64_t anchor = 0llu;
		uint64_t buffw = 0llu;  // rotating ref sequence buffer
		// OR the 'diff' buffer with this so that we can always count
		// 'N's as mismatches
		uint64_t diffMask = 0llu;
		// Set up a mask that we'll apply to the two bufs every round
		// to discard bits that were rotated out of the anchor area
		uint64_t clearMask = 0xffffffffffffffffllu;
		bool useMask = false;
		if(anchorBitPairs < 32) {
			clearMask >>= ((32-anchorBitPairs) << 1);
			useMask = true;
		}
		int nsInSeed = 0;
		uint32_t nPoss = 0;
		int nPos1 = -1;
		int nPos2 = -1;
		int nPos3 = -1;
		uint32_t skipLeftToRights = 0;
		uint32_t skipRightToLefts = 0;
		const uint32_t halfwayRi = halfway - begin;
		// Construct the 'anchor' 64-bit buffer so that it holds all of
		// the first 'anchorBitPairs' bit pairs of the query.
		for(uint32_t ii = 0; ii < anchorBitPairs; ii++) {
			uint32_t i = ii;
			if(!seedOnLeft) {
				// Fill in the anchor using characters from the right-
				// hand side of the query (but take the characters in
				// left-to-right order).  Be sure to subtract slen from
				// qlen; not anchorBitPairs from qlen.  We want the
				// characters in the seedAnchorOverhang region to be to
				// the right of the characters in the anchor.
				i = qlen - slen + ii;
			}
			int c = (int)qry[i]; // next query character
			int r = (int)ref[halfwayRi + ii]; // next reference character
			if(r & 4) {
				// The reference character is an N; to mimic the
				// behavior of BW alignment, we have to skip all
				// alignments that involve an N in the reference.  Set
				// the skip* variables accordingly.
				r = 0;
				uint32_t lrSkips = ii;
				uint32_t rlSkips = qlen - ii;
				if(!seedOnLeft && readSeedOverhang) {
					lrSkips += readSeedOverhang;
					assert_geq(rlSkips, readSeedOverhang);
					rlSkips -= readSeedOverhang;
				}
				// The right-to-left direction absorbs the candidate
				// alignment based at 'halfway'
				skipLeftToRights = max(skipLeftToRights, lrSkips);
				skipRightToLefts = max(skipRightToLefts, rlSkips);
				assert_leq(skipLeftToRights, qlen);
				assert_leq(skipRightToLefts, qlen);
			}
			assert_leq(c, 4);
			assert_lt(r, 4);
			// Special case: query has an 'N'
			if(c == 4) {
				if(++nsInSeed > 3) {
					// More than one 'N' in the anchor region; can't
					// possibly have a 1-mismatch hit anywhere
					assert_eq(r2.size(), results.size() - resultsISz);
					return;   // can't match if query has Ns
				}
				if(nsInSeed == 1) {
					nPos1 = (int)ii; // w/r/t LHS of anchor
				} else if(nsInSeed == 2) {
					nPos2 = (int)ii; // w/r/t LHS of anchor
					assert_gt(nPos2, nPos1);
				} else {
					assert_eq(3, nsInSeed);
					nPos3 = (int)ii; // w/r/t LHS of anchor
					assert_gt(nPos3, nPos2);
				}
				// Make it look like an 'A' in the anchor
				c = 0;
				diffMask = (diffMask << 2llu) | 1llu;
			} else {
				diffMask <<= 2llu;
			}
			anchor = ((anchor << 2llu) | c);
			buffw = ((buffw << 2llu) | r);
		}
		// Check whether read is disqualified by Ns inside the seed
		// region but outside the anchor region
		if(seedAnchorOverhang) {
			assert_lt(anchorBitPairs, slen);
			for(uint32_t ii = anchorBitPairs; ii < slen; ii++) {
				uint32_t i = ii;
				if(!seedOnLeft) {
					i = qlen - slen + ii;
				}
				if((int)qry[i] == 4) {
					if(++nsInSeed > 3) {
						assert_eq(r2.size(), results.size() - resultsISz);
						return; // can't match if query has Ns
					}
				}
			}
		} else {
			assert_eq(anchorBitPairs, slen);
		}
		uint64_t bufbw = buffw;
		// Slide the anchor out in either direction, alternating
		// between right-to-left and left-to-right shifts, until all of
		// the positions from qbegin to qend have been covered.
		bool hi = false;
		uint32_t riHi  = halfway;
		uint32_t rirHi = halfway - begin;
		uint32_t rirHiAnchor = rirHi + anchorBitPairs - 1;
		uint32_t riLo  = halfway + 1;
		uint32_t rirLo = halfway - begin + 1;
		uint32_t lrSkips = anchorBitPairs;
		uint32_t rlSkips = qlen;
		if(!seedOnLeft && readSeedOverhang) {
			lrSkips += readSeedOverhang;
			assert_geq(rlSkips, readSeedOverhang);
			rlSkips -= readSeedOverhang;
		}
		for(uint32_t i = 1; i <= lim + 1; i++) {
			int r;       // new reference char
			uint64_t diff;
			assert_leq(skipLeftToRights, qlen);
			assert_leq(skipRightToLefts, qlen);
			if(hi) {
				hi = false;
				// Moving left-to-right
				riHi++;
				rirHi++;
				rirHiAnchor++;
				r = (int)ref[rirHiAnchor];
				if(r & 4) {
					r = 0;
					skipLeftToRights = lrSkips;
				}
				assert_lt(r, 4);
				// Bring in new base pair at the least significant
				// position
				buffw = ((buffw << 2llu) | r);
				if(useMask) buffw &= clearMask;
				if(skipLeftToRights > 0) {
					skipLeftToRights--;
					continue;
				}
				diff = (buffw ^ anchor) | diffMask;
			} else {
				hi = true;
				// Moving right-to-left
				riLo--;
				rirLo--;
				r = (int)ref[rirLo];
				if(r & 4) {
					r = 0;
					skipRightToLefts = rlSkips;
				}
				assert_lt(r, 4);
				if(i >= 2) {
					bufbw >>= 2llu;
					// Bring in new base pair at the most significant
					// position
					bufbw |= ((uint64_t)r << lhsShift);
				}
				if(skipRightToLefts > 0) {
					skipRightToLefts--;
					continue;
				}
				diff = (bufbw ^ anchor) | diffMask;
			}
			if((diff & 0xf000f000f000f000llu) &&
			   (diff & 0x0f000f000f000f00llu) &&
			   (diff & 0x00f000f000f000f0llu) &&
			   (diff & 0x000f000f000f000fllu)) continue;
			if((diff & 0xc003c003c003c003llu) &&
			   (diff & 0x3c003c003c003c00llu) &&
			   (diff & 0x03c003c003c003c0llu) &&
			   (diff & 0x003c003c003c003cllu)) continue;
			uint32_t ri  = hi ? riLo  : riHi;
			uint32_t rir = hi ? rirLo : rirHi;
			// Could use pop count
			uint8_t *diff8 = reinterpret_cast<uint8_t*>(&diff);
			// As a first cut, see if there are too many mismatches in
			// the first and last parts of the anchor
			uint32_t diffs = u8toMms[(int)diff8[0]] + u8toMms[(int)diff8[7]];
			if(diffs > 3) continue;
			diffs += u8toMms[(int)diff8[1]] +
					 u8toMms[(int)diff8[2]] +
					 u8toMms[(int)diff8[3]] +
					 u8toMms[(int)diff8[4]] +
					 u8toMms[(int)diff8[5]] +
					 u8toMms[(int)diff8[6]];
			uint32_t mmpos1 = 0xffffffff;
			int refc1 = -1;
			int readc1 = -1;
			uint32_t mmpos2 = 0xffffffff;
			int refc2 = -1;
			int readc2 = -1;
			uint32_t mmpos3 = 0xffffffff;
			int refc3 = -1;
			int readc3 = -1;
			unsigned int ham = 0;
			if(diffs > 3) {
				// Too many differences in the seed; stop
				continue;
			} else if(nPoss > 1 && diffs == nPoss) {
				// There was one difference, but there was also one N,
				// so we already know where the difference is
				mmpos1 = nPos1;
				refc1 = "ACGT"[(int)ref[rir + nPos1]];
				readc1 = 'N'; // read char was an N
				if(!seedOnLeft) {
					mmpos1 += readSeedOverhang;
				}
				char q = quals[mmpos1];
				ham += mmPenalty(!gNoMaqRound, phredcToPhredq(q));
				if(ham > this->qualMax_) {
					// Exceeded quality limit
					continue;
				}
				if(nPoss > 1) {
					mmpos2 = nPos2;
					readc2 = 'N'; // read char was an N
					refc2 = "ACGT"[(int)ref[rir + nPos2]];
					if(!seedOnLeft) {
						mmpos2 += readSeedOverhang;
					}
					q = quals[mmpos2];
					ham += mmPenalty(!gNoMaqRound, phredcToPhredq(q));
					if(ham > this->qualMax_) {
						// Exceeded quality limit
						continue;
					}
					if(nPoss > 3) {
						mmpos3 = nPos3;
						readc3 = 'N'; // read char was an N
						refc3 = "ACGT"[(int)ref[rir + nPos3]];
						if(!seedOnLeft) {
							mmpos3 += readSeedOverhang;
						}
						q = quals[mmpos3];
						ham += mmPenalty(!gNoMaqRound, phredcToPhredq(q));
						if(ham > this->qualMax_) {
							// Exceeded quality limit
							continue;
						}
					}
				}
			} else if(diffs > 0) {
				// Figure out which position mismatched
				uint64_t diff2 = diff;
				mmpos1 = 31;
				if((diff & 0xffffffffllu) == 0) { diff >>= 32llu; mmpos1 -= 16; }
				assert_neq(0, diff);
				if((diff & 0xffffllu) == 0)     { diff >>= 16llu; mmpos1 -=  8; }
				assert_neq(0, diff);
				if((diff & 0xffllu) == 0)       { diff >>= 8llu;  mmpos1 -=  4; }
				assert_neq(0, diff);
				if((diff & 0xfllu) == 0)        { diff >>= 4llu;  mmpos1 -=  2; }
				assert_neq(0, diff);
				if((diff & 0x3llu) == 0)        { mmpos1--; }
				assert_neq(0, diff);
				assert_geq(mmpos1, 0);
				assert_lt(mmpos1, 32);
				readc1 = "ACGT"[(anchor >> (62-mmpos1*2)) & 3];
				if((diffMask >> (62-mmpos1*2)) & 3) {
					readc1 = 'N';
				}
				uint32_t savedMmpos1 = mmpos1;
				mmpos1 -= anchorCushion;
				assert_lt(mmpos1, anchorBitPairs);
				refc1 = "ACGT"[(int)ref[rir + mmpos1]];
				assert_neq(readc1, refc1);
				if(!seedOnLeft) {
					mmpos1 += readSeedOverhang;
				}
				char q = quals[mmpos1];
				ham += mmPenalty(!gNoMaqRound, phredcToPhredq(q));
				if(ham > this->qualMax_) {
					// Exceeded quality limit
					continue;
				}
				if(diffs > 1) {
					// Figure out the second mismatched position
					ASSERT_ONLY(uint64_t origDiff2 = diff2);
					diff2 &= ~(0xc000000000000000llu >> (uint64_t)((savedMmpos1) << 1));
					uint64_t diff3 = diff2;
					assert_neq(diff2, origDiff2);
					mmpos2 = 31;
					if((diff2 & 0xffffffffllu) == 0) { diff2 >>= 32llu; mmpos2 -= 16; }
					assert_neq(0, diff2);
					if((diff2 & 0xffffllu) == 0)     { diff2 >>= 16llu; mmpos2 -=  8; }
					assert_neq(0, diff2);
					if((diff2 & 0xffllu) == 0)       { diff2 >>= 8llu;  mmpos2 -=  4; }
					assert_neq(0, diff2);
					if((diff2 & 0xfllu) == 0)        { diff2 >>= 4llu;  mmpos2 -=  2; }
					assert_neq(0, diff2);
					if((diff2 & 0x3llu) == 0)        { mmpos2--; }
					assert_neq(0, diff2);
					assert_geq(mmpos2, 0);
					assert_lt(mmpos2, 32);
					readc2 = "ACGT"[(anchor >> (62-mmpos2*2)) & 3];
					if((diffMask >> (62-mmpos2*2)) & 3) {
						readc2 = 'N';
					}
					uint32_t savedMmpos2 = mmpos2;
					mmpos2 -= anchorCushion;
					assert_neq(mmpos1, mmpos2);
					refc2 = "ACGT"[(int)ref[rir + mmpos2]];
					assert_neq(readc2, refc2);
					if(!seedOnLeft) {
						mmpos2 += readSeedOverhang;
					}
					q = quals[mmpos2];
					ham += mmPenalty(!gNoMaqRound, phredcToPhredq(q));
					if(ham > this->qualMax_) {
						// Exceeded quality limit
						continue;
					}
					if(mmpos2 < mmpos1) {
						uint32_t mmtmp = mmpos1;
						mmpos1 = mmpos2;
						mmpos2 = mmtmp;
						int refctmp = refc1;
						refc1 = refc2;
						refc2 = refctmp;
						int readctmp = readc1;
						readc1 = readc2;
						readc2 = readctmp;
					}
					assert_lt(mmpos1, mmpos2);
					if(diffs > 2) {
						// Figure out the second mismatched position
						ASSERT_ONLY(uint32_t origDiff3 = diff3);
						diff3 &= ~(0xc000000000000000llu >> (uint64_t)((savedMmpos2) << 1));
						assert_neq(diff3, origDiff3);
						mmpos3 = 31;
						if((diff3 & 0xffffffffllu) == 0) { diff3 >>= 32llu; mmpos3 -= 16; }
						assert_neq(0, diff3);
						if((diff3 & 0xffffllu) == 0)     { diff3 >>= 16llu; mmpos3 -=  8; }
						assert_neq(0, diff3);
						if((diff3 & 0xffllu) == 0)       { diff3 >>= 8llu;  mmpos3 -=  4; }
						assert_neq(0, diff3);
						if((diff3 & 0xfllu) == 0)        { diff3 >>= 4llu;  mmpos3 -=  2; }
						assert_neq(0, diff3);
						if((diff3 & 0x3llu) == 0)        { mmpos3--; }
						assert_neq(0, diff3);
						assert_geq(mmpos3, 0);
						assert_lt(mmpos3, 32);
						readc3 = "ACGT"[(anchor >> (62-mmpos3*2)) & 3];
						if((diffMask >> (62-mmpos3*2)) & 3) {
							readc3 = 'N';
						}
						mmpos3 -= anchorCushion;
						assert_neq(mmpos2, mmpos3);
						assert_neq(mmpos1, mmpos3);
						refc3 = "ACGT"[(int)ref[rir + mmpos3]];
						assert_neq(readc3, refc3);
						if(!seedOnLeft) {
							mmpos3 += readSeedOverhang;
						}
						q = quals[mmpos3];
						ham += mmPenalty(!gNoMaqRound, phredcToPhredq(q));
						if(ham > this->qualMax_) {
							// Exceeded quality limit
							continue;
						}
						if(mmpos3 < mmpos1) {
							uint32_t mmtmp = mmpos1;
							mmpos1 = mmpos3;
							mmpos3 = mmpos2;
							mmpos2 = mmtmp;
							int refctmp = refc1;
							refc1 = refc3;
							refc3 = refc2;
							refc2 = refctmp;
							int readctmp = readc1;
							readc1 = readc3;
							readc3 = readc2;
							readc2 = readctmp;
						} else if(mmpos3 < mmpos2) {
							uint32_t mmtmp = mmpos2;
							mmpos2 = mmpos3;
							mmpos3 = mmtmp;
							int refctmp = refc2;
							refc2 = refc3;
							refc3 = refctmp;
							int readctmp = readc2;
							readc2 = readc3;
							readc3 = readctmp;
						}
						assert_lt(mmpos1, mmpos2);
						assert_lt(mmpos2, mmpos3);
					}
				}
			}
			// If the seed is longer than the anchor, then scan the
			// rest of the seed characters
			bool foundHit = true;
			if(seedAnchorOverhang) {
				for(uint32_t j = 0; j < seedAnchorOverhang; j++) {
					int rc = (int)ref[rir + anchorBitPairs + j];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							// Skip out of the seedAnchorOverhang
							assert_eq(0, skipRightToLefts);
							skipRightToLefts = seedAnchorOverhang - j - 1;
							if(seedOnLeft) {
								// ...and skip out of the rest of the read
								skipRightToLefts += readSeedOverhang;
							}
						} else {
							// Left-to-right
							// Skip out of the seedAnchorOverhang
							assert_eq(0, skipLeftToRights);
							skipLeftToRights = anchorBitPairs + j;
							if(!seedOnLeft) {
								// ...and skip out of the rest of the read
								skipLeftToRights += readSeedOverhang;
							}
						}
						foundHit = false; // Skip this candidate
						break;
					}
					uint32_t qoff = anchorBitPairs + j;
					if(!seedOnLeft) {
						qoff += readSeedOverhang;
					}
					assert_lt(qoff, qlen);
					int qc = (int)qry[qoff];
					if(qc != rc) {
						diffs++;
						if(diffs > 3) {
							foundHit = false;
							break;
						} else if(diffs == 3) {
							assert_eq(0xffffffff, mmpos3);
							mmpos3 = qoff;
							assert_eq(-1, refc3);
							refc3 = "ACGT"[(int)ref[rir + anchorBitPairs + j]];
							readc3 = "ACGTN"[qc];
						} else if(diffs == 2) {
							assert_eq(0xffffffff, mmpos2);
							mmpos2 = qoff;
							assert_eq(-1, refc2);
							refc2 = "ACGT"[(int)ref[rir + anchorBitPairs + j]];
							readc2 = "ACGTN"[qc];
						} else {
							assert_eq(1, diffs);
							assert_eq(0xffffffff, mmpos1);
							mmpos1 = qoff;
							assert_eq(-1, refc1);
							refc1 = "ACGT"[(int)ref[rir + anchorBitPairs + j]];
							readc1 = "ACGTN"[qc];
						}
						char q = phredcToPhredq(quals[qoff]);
						ham += mmPenalty(!gNoMaqRound, q);
						if(ham > this->qualMax_) {
							// Exceeded quality limit
							foundHit = false;
							break;
						}
					}
				}
				if(!foundHit) continue;
			}
			// If the read is longer than the seed, then scan the rest
			// of the read characters; mismatches no longer count
			// toward the stratum or the 1-mm limit.
			// Vectors for holding edit information
			nonSeedEdits.clear();
			int mms = diffs; // start counting total mismatches
			if((qlen - slen) > 0) {
				// Going left-to-right
				for(uint32_t j = 0; j < readSeedOverhang; j++) {
					uint32_t roff = rir + slen + j;
					uint32_t qoff = slen + j;
					if(!seedOnLeft) {
						assert_geq(roff, qlen);
						roff -= qlen;
						qoff = j;
					}
					int rc = (int)ref[roff];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							// Skip what's left of the readSeedOverhang
							skipRightToLefts = readSeedOverhang - j - 1;
							if(!seedOnLeft) {
								// ...and skip the seed if it's on the right
								skipRightToLefts += slen;
							}
						} else {
							// Left-to-right
							// Skip what we've matched of the overhang
							skipLeftToRights = j;
							if(seedOnLeft) {
								// ...and skip the seed if it's on the left
								skipLeftToRights += slen;
							}
						}
						foundHit = false; // Skip this candidate
						break;
					}
					int qc = (int)qry[qoff];
					if(qc != rc) {
						// Calculate quality of mismatched base
						char q = phredcToPhredq(quals[qoff]);
						ham += mmPenalty(!gNoMaqRound, q);
						if(ham > this->qualMax_) {
							// Exceeded quality limit
							foundHit = false;
							break;
						}
						// Legal mismatch outside of the anchor; record it
						mms++;
						nonSeedEdits.push_back(Edit(qoff, "ACGT"[rc], "ACGTN"[qc], EDIT_TYPE_MM));
						assert_leq(nonSeedEdits.size(), (size_t)mms);
					}
				}
				if(!foundHit) continue;
			}
			assert(foundHit);
			// Adjust ri if seed is on the right-hand side
			if(!seedOnLeft) {
				ri -= readSeedOverhang;
				rir -= readSeedOverhang;
			}
			if(pairs != NULL) {
				TU64Pair p;
				if(ri < aoff) {
					// By convention, the upstream mate's
					// coordinates go in the 'first' field
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.second = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				} else {
					p.second = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				}
				if(pairs->find(p) != pairs->end()) {
					// We already found this hit!  Continue.
					ASSERT_ONLY(duplicates++);
					ASSERT_ONLY(r2i++);
					continue;
				} else {
					// Record this hit
					pairs->insert(p);
				}
			}
			if(gVerbose) {
				cout << "About to report:" << endl;
				cout << "  ";
				for(size_t i = 0; i < qlen; i++) {
					cout << (char)qry[i];
				}
				cout << endl;
				cout << "  ";
				for(size_t i = 0; i < qlen; i++) {
					cout << "ACGT"[ref[rir+i]];
				}
				cout << endl;
			}
			assert_leq(diffs, 3);
			assert_geq((size_t)mms, diffs);
			assert_lt(r2i, r2.size());
			assert_eq(r2[r2i].off, ri);
			results.expand();
			RefAlignerHit& result = results.back();
			result.clear();
			assert_eq((size_t)mms, r2[r2i].edits.size());
			result.stratum = diffs;
			result.cost = ham;
			if(mms > 0) {
				ASSERT_ONLY(size_t mmcur = 0);
				if(seedOnLeft && diffs > 0) {
					assert_neq(mmpos1, 0xffffffff);
					assert_lt(mmpos1, qlen);
					assert_lt(mmcur, (size_t)mms);
					assert_eq(mmpos1, r2[r2i].edits[mmcur].pos);
					assert_in(refc1, "ACGT");
					assert_in(readc1, "ACGTN");
					assert_eq((unsigned)refc1, r2[r2i].edits[mmcur].chr);
					ASSERT_ONLY(mmcur++);
					result.edits.push_back(Edit(mmpos1, refc1, readc1, EDIT_TYPE_MM));
					if(diffs > 1) {
						assert_neq(mmpos2, 0xffffffff);
						assert_lt(mmpos2, qlen);
						assert_lt(mmcur, (size_t)mms);
						assert_eq(mmpos2, r2[r2i].edits[mmcur].pos);
						assert_in(refc2, "ACGT");
						assert_in(readc2, "ACGTN");
						assert_eq((unsigned)refc2, r2[r2i].edits[mmcur].chr);
						ASSERT_ONLY(mmcur++);
						result.edits.push_back(Edit(mmpos2, refc2, readc2, EDIT_TYPE_MM));
						if(diffs > 2) {
							assert_eq(3, diffs);
							assert_neq(mmpos3, 0xffffffff);
							assert_lt(mmpos3, qlen);
							assert_lt(mmcur, (size_t)mms);
							assert_eq(mmpos3, r2[r2i].edits[mmcur].pos);
							assert_in(refc3, "ACGT");
							assert_in(readc3, "ACGTN");
							assert_eq((unsigned)refc3, r2[r2i].edits[mmcur].chr);
							ASSERT_ONLY(mmcur++);
							result.edits.push_back(Edit(mmpos3, refc3, readc3, EDIT_TYPE_MM));
						}
					}
				}
				const size_t nonSeedEditsSz = nonSeedEdits.size();
				for(size_t i = 0; i < nonSeedEditsSz; i++) {
					assert(nonSeedEdits[i].initialized());
					assert_lt(mmcur, (size_t)mms);
					assert_eq(nonSeedEdits[i].pos, r2[r2i].edits[mmcur].pos);
					result.edits.push_back(nonSeedEdits[i]);
					assert_eq(nonSeedEdits[i].chr, r2[r2i].edits[mmcur].chr);
					ASSERT_ONLY(mmcur++);
				}
				if(!seedOnLeft && diffs > 0) {
					assert_neq(mmpos1, 0xffffffff);
					assert_lt(mmpos1, qlen);
					assert_lt(mmcur, (size_t)mms);
					assert_eq(mmpos1, r2[r2i].edits[mmcur].pos);
					assert_in(refc1, "ACGT");
					assert_in(readc1, "ACGTN");
					assert_eq((unsigned)refc1, r2[r2i].edits[mmcur].chr);
					ASSERT_ONLY(mmcur++);
					result.edits.push_back(Edit(mmpos1, refc1, readc1, EDIT_TYPE_MM));
					if(diffs > 1) {
						assert_neq(mmpos2, 0xffffffff);
						assert_lt(mmpos2, qlen);
						assert_lt(mmcur, (size_t)mms);
						assert_eq(mmpos2, r2[r2i].edits[mmcur].pos);
						assert_in(refc2, "ACGT");
						assert_in(readc2, "ACGTN");
						assert_eq((unsigned)refc2, r2[r2i].edits[mmcur].chr);
						ASSERT_ONLY(mmcur++);
						result.edits.push_back(Edit(mmpos2, refc2, readc2, EDIT_TYPE_MM));
						if(diffs > 2) {
							assert_eq(3, diffs);
							assert_neq(mmpos3, 0xffffffff);
							assert_lt(mmpos3, qlen);
							assert_lt(mmcur, (size_t)mms);
							assert_eq(mmpos3, r2[r2i].edits[mmcur].pos);
							assert_in(refc3, "ACGT");
							assert_in(readc3, "ACGTN");
							assert_eq((unsigned)refc3, r2[r2i].edits[mmcur].chr);
							ASSERT_ONLY(mmcur++);
							result.edits.push_back(Edit(mmpos3, refc3, readc3, EDIT_TYPE_MM));
						}
					}
				}
				assert_eq(mmcur, r2[r2i].edits.size());
			}
			assert_eq((size_t)mms, result.edits.size());
			ASSERT_ONLY(r2i++);
			result.off = ri;
			if(--numToFind == 0) return;
		}
		assert_leq(duplicates, r2.size());
		assert_geq(r2.size() - duplicates, results.size() - resultsISz);
		return; // no hit
	}

  private:

	EList<Edit> nonSeedEdits;
};

#endif /* REF_ALIGNER_H_ */
