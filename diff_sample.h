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

#ifndef DIFF_SAMPLE_H_
#define DIFF_SAMPLE_H_

#include <stdint.h>
#include <string.h>
#include "assert_helpers.h"
#include "multikey_qsort.h"
#include "timer.h"
#include "ds.h"
#include "mem_ids.h"
#include "ls.h"

using namespace std;

#ifndef VMSG_NL
#define VMSG_NL(...) \
if(this->verbose()) { \
	stringstream tmp; \
	tmp << __VA_ARGS__ << endl; \
	this->verbose(tmp.str()); \
}
#endif

#ifndef VMSG
#define VMSG(...) \
if(this->verbose()) { \
	stringstream tmp; \
	tmp << __VA_ARGS__; \
	this->verbose(tmp.str()); \
}
#endif

/**
 * Routines for calculating, sanity-checking, and dispensing difference
 * cover samples to clients.
 */

/**
 *
 */
struct sampleEntry {
	uint32_t maxV;
	uint32_t numSamples;
	uint32_t samples[128];
};

/// Array of Colbourn and Ling calculated difference covers up to
/// r = 16 (maxV = 5953)
extern struct sampleEntry clDCs[16];
extern bool clDCs_calced; /// have clDCs been calculated?

/**
 * Check that the given difference cover 'ds' actually covers all
 * differences for a periodicity of v.
 */
template<typename T>
static bool dcRepOk(T v, EList<T>& ds) {
	// diffs[] records all the differences observed
	AutoArray<bool> covered(v, EBWT_CAT);
	for(T i = 1; i < v; i++) {
		covered[i] = false;
	}
	for(T di = T(); di < ds.size(); di++) {
		for(T dj = di+1; dj < ds.size(); dj++) {
			assert_lt(ds[di], ds[dj]);
			T d1 = (ds[dj] - ds[di]);
			T d2 = (ds[di] + v - ds[dj]);
			assert_lt(d1, v);
			assert_lt(d2, v);
			covered[d1] = true;
			covered[d2] = true;
		}
	}
	bool ok = true;
	for(T i = 1; i < v; i++) {
		if(covered[i] == false) {
			ok = false;
			break;
		}
	}
	return ok;
}

/**
 * Return true iff each element of ts (with length 'limit') is greater
 * than the last.
 */
template<typename T>
static bool increasing(T* ts, size_t limit) {
	for(size_t i = 0; i < limit-1; i++) {
		if(ts[i+1] <= ts[i]) return false;
	}
	return true;
}

/**
 * Return true iff the given difference cover covers difference 'diff'
 * mod 'v'.
 */
template<typename T>
static inline bool hasDifference(T *ds, T d, T v, T diff) {
	// diffs[] records all the differences observed
	for(T di = T(); di < d; di++) {
		for(T dj = di+1; dj < d; dj++) {
			assert_lt(ds[di], ds[dj]);
			T d1 = (ds[dj] - ds[di]);
			T d2 = (ds[di] + v - ds[dj]);
			assert_lt(d1, v);
			assert_lt(d2, v);
			if(d1 == diff || d2 == diff) return true;
		}
	}
	return false;
}

/**
 * Exhaustively calculate optimal difference cover samples for v = 4,
 * 8, 16, 32, 64, 128, 256 and store results in p2DCs[]
 */
template<typename T>
void calcExhaustiveDC(T i, bool verbose = false, bool sanityCheck = false) {
	T v = i;
	AutoArray<bool> diffs(v, EBWT_CAT);
	// v is the target period
	T ld = (T)ceil(sqrt(v));
	// ud is the upper bound on |D|
	T ud = v / 2;
	// for all possible |D|s
	bool ok = true;
	T *ds = NULL;
	T d;
	for(d = ld; d <= ud+1; d++) {
		// for all possible |D| samples
		AutoArray<T> ds(d, EBWT_CAT);
		for(T j = 0; j < d; j++) {
			ds[j] = j;
		}
		assert(increasing(ds, d));
		while(true) {
			// reset diffs[]
			for(T t = 1; t < v; t++) {
				diffs[t] = false;
			}
			T diffCnt = 0;
			// diffs[] records all the differences observed
			for(T di = 0; di < d; di++) {
				for(T dj = di+1; dj < d; dj++) {
					assert_lt(ds[di], ds[dj]);
					T d1 = (ds[dj] - ds[di]);
					T d2 = (ds[di] + v - ds[dj]);
					assert_lt(d1, v);
					assert_lt(d2, v);
					assert_gt(d1, 0);
					assert_gt(d2, 0);
					if(!diffs[d1]) diffCnt++; diffs[d1] = true;
					if(!diffs[d2]) diffCnt++; diffs[d2] = true;
				}
			}
			// Do we observe all possible differences (except 0)
			ok = diffCnt == v-1;
			if(ok) {
				// Yes, all differences are covered
				break;
			} else {
				// Advance ds
				// (Following is commented out because it turns out
				// it's slow)
				// Find a missing difference
				//uint32_t missing = 0xffffffff;
				//for(uint32_t t = 1; t < v; t++) {
				//	if(diffs[t] == false) {
				//		missing = diffs[t];
				//		break;
				//	}
				//}
				//assert_neq(missing, 0xffffffff);
				assert(increasing(ds, d));
				bool advanced = false;
				bool keepGoing = false;
				do {
					keepGoing = false;
					for(T bd = d-1; bd > 1; bd--) {
						T dif = (d-1)-bd;
						if(ds[bd] < v-1-dif) {
							ds[bd]++;
							assert_neq(0, ds[bd]);
							// Reset subsequent ones
							for(T bdi = bd+1; bdi < d; bdi++) {
								assert_eq(0, ds[bdi]);
								ds[bdi] = ds[bdi-1]+1;
								assert_gt(ds[bdi], ds[bdi-1]);
							}
							assert(increasing(ds, d));
							// (Following is commented out because
							// it turns out it's slow)
							// See if the new DC has the missing value
							//if(!hasDifference(ds, d, v, missing)) {
							//	keepGoing = true;
							//	break;
							//}
							advanced = true;
							break;
						} else {
							ds[bd] = 0;
							// keep going
						}
					}
				} while(keepGoing);
				// No solution for this |D|
				if(!advanced) break;
				assert(increasing(ds, d));
			}
		} // next sample assignment
		if(ok) {
			break;
		}
	} // next |D|
	assert(ok);
	cout << "Did exhaustive v=" << v << " |D|=" << d << endl;
	cout << "  ";
	for(T i = 0; i < d; i++) {
		cout << ds[i];
		if(i < d-1) cout << ",";
	}
	cout << endl;
}

/**
 * Routune for calculating the elements of clDCs up to r = 16 using the
 * technique of Colbourn and Ling.
 *
 * See http://citeseer.ist.psu.edu/211575.html
 */
template <typename T>
void calcColbournAndLingDCs(bool verbose = false, bool sanityCheck = false) {
	for(T r = 0; r < 16; r++) {
		T maxv = 24*r*r + 36*r + 13; // Corollary 2.3
		T numsamp = 6*r + 4;
		clDCs[r].maxV = maxv;
		clDCs[r].numSamples = numsamp;
		memset(clDCs[r].samples, 0, 4 * 128);
		T i;
		// clDCs[r].samples[0] = 0;
		// Fill in the 1^r part of the B series
		for(i = 1; i < r+1; i++) {
			clDCs[r].samples[i] = clDCs[r].samples[i-1] + 1;
		}
		// Fill in the (r + 1)^1 part
		clDCs[r].samples[r+1] = clDCs[r].samples[r] + r + 1;
		// Fill in the (2r + 1)^r part
		for(i = r+2; i < r+2+r; i++) {
			clDCs[r].samples[i] = clDCs[r].samples[i-1] + 2*r + 1;
		}
		// Fill in the (4r + 3)^(2r + 1) part
		for(i = r+2+r; i < r+2+r+2*r+1; i++) {
			clDCs[r].samples[i] = clDCs[r].samples[i-1] + 4*r + 3;
		}
		// Fill in the (2r + 2)^(r + 1) part
		for(i = r+2+r+2*r+1; i < r+2+r+2*r+1+r+1; i++) {
			clDCs[r].samples[i] = clDCs[r].samples[i-1] + 2*r + 2;
		}
		// Fill in the last 1^r part
		for(i = r+2+r+2*r+1+r+1; i < r+2+r+2*r+1+r+1+r; i++) {
			clDCs[r].samples[i] = clDCs[r].samples[i-1] + 1;
		}
		assert_eq(i, numsamp);
		assert_lt(i, 128);
		if(sanityCheck) {
			// diffs[] records all the differences observed
			AutoArray<bool> diffs(maxv, EBWT_CAT);
			for(T i = 0; i < numsamp; i++) {
				for(T j = i+1; j < numsamp; j++) {
					T d1 = (clDCs[r].samples[j] - clDCs[r].samples[i]);
					T d2 = (clDCs[r].samples[i] + maxv - clDCs[r].samples[j]);
					assert_lt(d1, maxv);
					assert_lt(d2, maxv);
					diffs[d1] = true;
					diffs[d2] = true;
				}
			}
			// Should have observed all possible differences (except 0)
			for(T i = 1; i < maxv; i++) {
				if(diffs[i] == false) cout << r << ", " << i << endl;
				assert(diffs[i] == true);
			}
		}
	}
	clDCs_calced = true;
}

/**
 * A precalculated list of difference covers.
 */
extern uint32_t dc0to64[65][10];

/**
 * Get a difference cover for the requested periodicity v.
 */
template <typename T>
static EList<T> getDiffCover(
	T v,
	bool verbose = false,
	bool sanityCheck = false)
{
	assert_gt(v, 2);
	EList<T> ret;
	ret.clear();
	// Can we look it up in our hardcoded array?
	if(v <= 64 && dc0to64[v][0] == 0xffffffff) {
		if(verbose) cout << "v in hardcoded area, but hardcoded entry was all-fs" << endl;
		return ret;
	} else if(v <= 64) {
		ret.push_back(0);
		for(size_t i = 0; i < 10; i++) {
			if(dc0to64[v][i] == 0) break;
			ret.push_back(dc0to64[v][i]);
		}
		if(sanityCheck) assert(dcRepOk(v, ret));
		return ret;
	}

	// Can we look it up in our calcColbournAndLingDCs array?
	if(!clDCs_calced) {
		calcColbournAndLingDCs<uint32_t>(verbose, sanityCheck);
		assert(clDCs_calced);
	}
	for(size_t i = 0; i < 16; i++) {
		if(v <= clDCs[i].maxV) {
			for(size_t j = 0; j < clDCs[i].numSamples; j++) {
				T s = clDCs[i].samples[j];
				if(s >= v) {
					s %= v;
					for(size_t k = 0; k < ret.size(); k++) {
						if(s == ret[k]) break;
						if(s < ret[k]) {
							ret.insert(s, k);
							break;
						}
					}
				} else {
					ret.push_back(s % v);
				}
			}
			if(sanityCheck) assert(dcRepOk(v, ret));
			return ret;
		}
	}
	cerr << "Error: Could not find a difference cover sample for v=" << v << endl;
	throw 1;
}

/**
 * Calculate and return a delta map based on the given difference cover
 * and periodicity v.
 */
template <typename T>
static EList<T> getDeltaMap(T v, const EList<T>& dc) {
	// Declare anchor-map-related items
	EList<T> amap;
	size_t amapEnts = 1;
	amap.resizeExact((size_t)v);
	amap.fill(0xffffffff);
	amap[0] = 0;
	// Print out difference cover (and optionally calculate
	// anchor map)
	for(size_t i = 0; i < dc.size(); i++) {
		for(size_t j = i+1; j < dc.size(); j++) {
			assert_gt(dc[j], dc[i]);
			T diffLeft  = dc[j] - dc[i];
			T diffRight = dc[i] + v - dc[j];
			assert_lt(diffLeft, v);
			assert_lt(diffRight, v);
			if(amap[diffLeft] == 0xffffffff) {
				amap[diffLeft] = dc[i];
				amapEnts++;
			}
			if(amap[diffRight] == 0xffffffff) {
				amap[diffRight] = dc[j];
				amapEnts++;
			}
		}
	}
	return amap;
}

/**
 * Return population count (count of all bits set to 1) of i.
 */
template<typename T>
static unsigned int popCount(T i) {
	unsigned int cnt = 0;
	for(size_t j = 0; j < sizeof(T)*8; j++) {
		if(i & 1) cnt++;
		i >>= 1;
	}
	return cnt;
}

/**
 * Calculate log-base-2 of i
 */
template<typename T>
static unsigned int myLog2(T i) {
	assert_eq(1, popCount(i)); // must be power of 2
	for(size_t j = 0; j < sizeof(T)*8; j++) {
		if(i & 1) return (int)j;
		i >>= 1;
	}
	assert(false);
	return 0xffffffff;
}

/**
 *
 */
template<typename TStr>
class DifferenceCoverSample {
public:

	DifferenceCoverSample(const TStr& __text,
	                      uint32_t __v,
	                      bool __verbose = false,
	                      bool __sanity = false,
	                      ostream& __logger = cout) :
		_text(__text),
		_v(__v),
		_verbose(__verbose),
		_sanity(__sanity),
		_ds(getDiffCover(_v, _verbose, _sanity)),
		_dmap(getDeltaMap(_v, _ds)),
		_d((uint32_t)_ds.size()),
		_doffs(),
		_isaPrime(),
		_dInv(),
		_log2v(myLog2(_v)),
		_vmask(0xffffffff << _log2v),
		_logger(__logger)
	{
		assert_gt(_d, 0);
		assert_eq(1, popCount(_v)); // must be power of 2
		// Build map from d's to idx's
		_dInv.resizeExact((size_t)v());
		_dInv.fill(0xffffffff);
		uint32_t lim = (uint32_t)_ds.size();
		for(uint32_t i = 0; i < lim; i++) {
			_dInv[_ds[i]] = i;
		}
	}
	
	/**
	 * Allocate an amount of memory that simulates the peak memory
	 * usage of the DifferenceCoverSample with the given text and v.
	 * Throws bad_alloc if it's not going to fit in memory.  Returns
	 * the approximate number of bytes the Cover takes at all times.
	 */
	static size_t simulateAllocs(const TStr& text, uint32_t v) {
		EList<uint32_t> ds(getDiffCover(v, false /*verbose*/, false /*sanity*/));
		size_t len = text.length();
		size_t sPrimeSz = (len / v) * ds.size();
		// sPrime, sPrimeOrder, _isaPrime all exist in memory at
		// once and that's the peak
		AutoArray<uint32_t> aa(sPrimeSz * 3 + (1024 * 1024 /*out of caution*/), EBWT_CAT);
		return sPrimeSz * 4; // sPrime array
	}

	uint32_t v() const                   { return _v; }
	uint32_t log2v() const               { return _log2v; }
	uint32_t vmask() const               { return _vmask; }
	uint32_t modv(uint32_t i) const      { return i & ~_vmask; }
	uint32_t divv(uint32_t i) const      { return i >> _log2v; }
	uint32_t d() const                   { return _d; }
	bool verbose() const                 { return _verbose; }
	bool sanityCheck() const             { return _sanity; }
	const TStr& text() const             { return _text; }
	const EList<uint32_t>& ds() const    { return _ds; }
	const EList<uint32_t>& dmap() const  { return _dmap; }
	ostream& log() const                 { return _logger; }

	void     build();
	uint32_t tieBreakOff(uint32_t i, uint32_t j) const;
	int64_t  breakTie(uint32_t i, uint32_t j) const;
	bool     isCovered(uint32_t i) const;
	uint32_t rank(uint32_t i) const;

	/**
	 * Print out the suffix array such that every sample offset has its
	 * rank filled in and every non-sample offset is shown as '-'.
	 */
	void print(ostream& out) {
		for(size_t i = 0; i < _text.length(); i++) {
			if(isCovered(i)) {
				out << rank(i);
			} else {
				out << "-";
			}
			if(i < _text.length()-1) {
				out << ",";
			}
		}
		out << endl;
	}

private:

	void doBuiltSanityCheck() const;
	void buildSPrime(EList<uint32_t>& sPrime, size_t padding);

	bool built() const {
		return _isaPrime.size() > 0;
	}

	void verbose(const string& s) const {
		if(this->verbose()) {
			this->log() << s;
			this->log().flush();
		}
	}

	const TStr&      _text;     // text to sample
	uint32_t         _v;        // periodicity of sample
	bool             _verbose;  //
	bool             _sanity;   //
	EList<uint32_t>  _ds;       // samples: idx -> d
	EList<uint32_t>  _dmap;     // delta map
	uint32_t         _d;        // |D| - size of sample
	EList<uint32_t>  _doffs;    // offsets into sPrime/isaPrime for each d idx
	EList<uint32_t>  _isaPrime; // ISA' array
	EList<uint32_t>  _dInv;     // Map from d -> idx
	uint32_t         _log2v;
	uint32_t         _vmask;
	ostream&         _logger;
};

/**
 * Sanity-check the difference cover by first inverting _isaPrime then
 * checking that each successive suffix really is less than the next.
 */
template <typename TStr>
void DifferenceCoverSample<TStr>::doBuiltSanityCheck() const {
	uint32_t v = this->v();
	assert(built());
	VMSG_NL("  Doing sanity check");
	uint32_t added = 0;
	EList<uint32_t> sorted;
	sorted.resizeExact(_isaPrime.size());
	sorted.fill(0xffffffff);
	for(size_t di = 0; di < this->d(); di++) {
		uint32_t d = _ds[di];
		size_t i = 0;
		for(size_t doi = _doffs[di]; doi < _doffs[di+1]; doi++, i++) {
			assert_eq(0xffffffff, sorted[_isaPrime[doi]]);
			// Maps the offset of the suffix to its rank
			sorted[_isaPrime[doi]] = (uint32_t)(v*i + d);
			added++;
		}
	}
	assert_eq(added, _isaPrime.size());
#ifndef NDEBUG
	for(size_t i = 0; i < sorted.size()-1; i++) {
		assert(sstr_suf_lt(this->text(), sorted[i], this->text(), sorted[i+1], false));
	}
#endif
}

/**
 * Build the s' array by sampling suffixes (suffix offsets, actually)
 * from t according to the difference-cover sample and pack them into
 * an array of machine words in the order dictated by the "mu" mapping
 * described in Burkhardt.
 *
 * Also builds _doffs map.
 */
template <typename TStr>
void DifferenceCoverSample<TStr>::buildSPrime(
	EList<uint32_t>& sPrime,
	size_t padding)
{
	const TStr& t = this->text();
	const EList<uint32_t>& ds = this->ds();
	uint32_t tlen = (uint32_t)t.length();
	uint32_t v = this->v();
	uint32_t d = this->d();
	assert_gt(v, 2);
	assert_lt(d, v);
	// Record where each d section should begin in sPrime
	uint32_t tlenDivV = this->divv(tlen);
	uint32_t tlenModV = this->modv(tlen);
	uint32_t sPrimeSz = 0;
	assert(_doffs.empty());
	_doffs.resizeExact((size_t)d+1);
	for(uint32_t di = 0; di < d; di++) {
		// mu mapping
		uint32_t sz = tlenDivV + ((ds[di] <= tlenModV) ? 1 : 0);
		assert_geq(sz, 0);
		_doffs[di] = sPrimeSz;
		sPrimeSz += sz;
	}
	_doffs[d] = sPrimeSz;
#ifndef NDEBUG
	if(tlenDivV > 0) {
		for(size_t i = 0; i < d; i++) {
			assert_gt(_doffs[i+1], _doffs[i]);
			uint32_t diff = _doffs[i+1] - _doffs[i];
			assert(diff == tlenDivV || diff == tlenDivV+1);
		}
	}
#endif
	assert_eq(_doffs.size(), d+1);
	// Size sPrime appropriately
	sPrime.resizeExact((size_t)sPrimeSz + padding);
	sPrime.fill(0xffffffff);
	// Slot suffixes from text into sPrime according to the mu
	// mapping; where the mapping would leave a blank, insert a 0
	uint32_t added = 0;
	uint32_t i = 0;
	for(uint32_t ti = 0; ti <= tlen; ti += v) {
		for(uint32_t di = 0; di < d; di++) {
			uint32_t tti = ti + ds[di];
			if(tti > tlen) break;
			uint32_t spi = _doffs[di] + i;
			assert_lt(spi, _doffs[di+1]);
			assert_leq(tti, tlen);
			assert_lt(spi, sPrimeSz);
			assert_eq(0xffffffff, sPrime[spi]);
			sPrime[spi] = tti; added++;
		}
		i++;
	}
	assert_eq(added, sPrimeSz);
}

/**
 * Return true iff suffixes with offsets suf1 and suf2 out of host
 * string 'host' are identical up to depth 'v'.
 */
template <typename TStr>
static inline bool suffixSameUpTo(
	const TStr& host,
	uint32_t suf1,
	uint32_t suf2,
	uint32_t v)
{
	for(uint32_t i = 0; i < v; i++) {
		bool endSuf1 = suf1+i >= host.length();
		bool endSuf2 = suf2+i >= host.length();
		if((endSuf1 && !endSuf2) || (!endSuf1 && endSuf2)) return false;
		if(endSuf1 && endSuf2) return true;
		if(host[suf1+i] != host[suf2+i]) return false;
	}
	return true;
}

/**
 * Calculates a ranking of all suffixes in the sample and stores them,
 * packed according to the mu mapping, in _isaPrime.
 */
template <typename TStr>
void DifferenceCoverSample<TStr>::build() {
	// Local names for relevant types
	VMSG_NL("Building DifferenceCoverSample");
	// Local names for relevant data
	const TStr& t = this->text();
	uint32_t v = this->v();
	assert_gt(v, 2);
	// Build s'
	EList<uint32_t> sPrime;
	// Need to allocate 2 extra elements at the end of the sPrime and _isaPrime
	// arrays.  One element that's less than all others, and another that acts
	// as needed padding for the Larsson-Sadakane sorting code.
	size_t padding = 1;
	VMSG_NL("  Building sPrime");
	buildSPrime(sPrime, padding);
	size_t sPrimeSz = sPrime.size() - padding;
	assert_gt(sPrime.size(), padding);
	assert_leq(sPrime.size(), t.length() + padding + 1);
	uint32_t nextRank = 0;
	{
		VMSG_NL("  Building sPrimeOrder");
		EList<uint32_t> sPrimeOrder;
		sPrimeOrder.resizeExact(sPrimeSz);
		for(uint32_t i = 0; i < sPrimeSz; i++) {
			sPrimeOrder[i] = i;
		}
		// sPrime now holds suffix-offsets for DC samples.
		{
			Timer timer(cout, "  V-Sorting samples time: ", this->verbose());
			VMSG_NL("  V-Sorting samples");
			// Extract backing-store array from sPrime and sPrimeOrder;
			// the mkeyQSortSuf2 routine works on the array for maximum
			// efficiency
			uint32_t *sPrimeArr = (uint32_t*)sPrime.ptr();
			assert_eq(sPrimeArr[0], sPrime[0]);
			assert_eq(sPrimeArr[sPrimeSz-1], sPrime[sPrimeSz-1]);
			uint32_t *sPrimeOrderArr = (uint32_t*)sPrimeOrder.ptr();
			assert_eq(sPrimeOrderArr[0], sPrimeOrder[0]);
			assert_eq(sPrimeOrderArr[sPrimeSz-1], sPrimeOrder[sPrimeSz-1]);
			// Sort sample suffixes up to the vth character using a
			// multikey quicksort.  Sort time is proportional to the
			// number of samples times v.  It isn't quadratic.
			// sPrimeOrder is passed in as a swapping partner for
			// sPrimeArr, i.e., every time the multikey qsort swaps
			// elements in sPrime, it swaps the same elements in
			// sPrimeOrder too.  This allows us to easily reconstruct
			// what the sort did.
			mkeyQSortSuf2(t, sPrimeArr, sPrimeSz, sPrimeOrderArr, 4,
			              this->verbose(), this->sanityCheck(), v);
			// Make sure sPrime and sPrimeOrder are consistent with
			// their respective backing-store arrays
			assert_eq(sPrimeArr[0], sPrime[0]);
			assert_eq(sPrimeArr[sPrimeSz-1], sPrime[sPrimeSz-1]);
			assert_eq(sPrimeOrderArr[0], sPrimeOrder[0]);
			assert_eq(sPrimeOrderArr[sPrimeSz-1], sPrimeOrder[sPrimeSz-1]);
		}
		// Now assign the ranking implied by the sorted sPrime/sPrimeOrder
		// arrays back into sPrime.
		VMSG_NL("  Allocating rank array");
		_isaPrime.resizeExact(sPrime.size());
		ASSERT_ONLY(_isaPrime.fill(0xffffffff));
		assert_gt(_isaPrime.size(), 0);
		{
			Timer timer(cout, "  Ranking v-sort output time: ", this->verbose());
			VMSG_NL("  Ranking v-sort output");
			for(size_t i = 0; i < sPrimeSz-1; i++) {
				// Place the appropriate ranking
				_isaPrime[sPrimeOrder[i]] = nextRank;
				// If sPrime[i] and sPrime[i+1] are identical up to v, then we
				// should give the next suffix the same rank
				if(!suffixSameUpTo(t, sPrime[i], sPrime[i+1], v)) nextRank++;
			}
			_isaPrime[sPrimeOrder[sPrimeSz-1]] = nextRank; // finish off
#ifndef NDEBUG
			for(size_t i = 0; i < sPrimeSz; i++) {
				assert_neq(0xffffffff, _isaPrime[i]);
				assert_lt(_isaPrime[i], sPrimeSz);
			}
#endif
		}
		// sPrimeOrder is destroyed
		// All the information we need is now in _isaPrime
	}
	_isaPrime[_isaPrime.size()-1] = (uint32_t)sPrimeSz;
	sPrime[sPrime.size()-1] = (uint32_t)sPrimeSz;
	// _isaPrime[_isaPrime.size()-1] and sPrime[sPrime.size()-1] are just
	// spacer for the Larsson-Sadakane routine to use
	{
		Timer timer(cout, "  Invoking Larsson-Sadakane on ranks time: ", this->verbose());
		VMSG_NL("  Invoking Larsson-Sadakane on ranks");
		if(sPrime.size() >= 0x10000000) {
			cerr << "Error; sPrime array has so many elements that it can't be converted to a signed array without overflow." << endl;
			throw 1;
		}
		LarssonSadakane<int> ls;
		ls.suffixsort(
			(int*)_isaPrime.ptr(),
			(int*)sPrime.ptr(),
			(int)sPrimeSz,
			(int)sPrime.size(),
			0);
	}
	// chop off final character of _isaPrime
	_isaPrime.resizeExact(sPrimeSz);
	for(size_t i = 0; i < _isaPrime.size(); i++) {
		_isaPrime[i]--;
	}
#ifndef NDEBUG
	for(size_t i = 0; i < sPrimeSz-1; i++) {
		assert_lt(_isaPrime[i], sPrimeSz);
		assert(i == 0 || _isaPrime[i] != _isaPrime[i-1]);
	}
#endif
	VMSG_NL("  Sanity-checking and returning");
	if(this->sanityCheck()) doBuiltSanityCheck();
}

/**
 * Return true iff index i within the text is covered by the difference
 * cover sample.  Allow i to be off the end of the text; simplifies
 * logic elsewhere.
 */
template <typename TStr>
bool DifferenceCoverSample<TStr>::isCovered(uint32_t i) const {
	assert(built());
	uint32_t modi = this->modv(i);
	assert_lt(modi, _dInv.size());
	return _dInv[modi] != 0xffffffff;
}

/**
 * Given a text offset that's covered, return its lexicographical rank
 * among the sample suffixes.
 */
template <typename TStr>
uint32_t DifferenceCoverSample<TStr>::rank(uint32_t i) const {
	assert(built());
	assert_lt(i, this->text().length());
	uint32_t imodv = this->modv(i);
	assert_neq(0xffffffff, _dInv[imodv]); // must be in the sample
	uint32_t ioff = this->divv(i);
	assert_lt(ioff, _doffs[_dInv[imodv]+1] - _doffs[_dInv[imodv]]);
	uint32_t isaIIdx = _doffs[_dInv[imodv]] + ioff;
	assert_lt(isaIIdx, _isaPrime.size());
	uint32_t isaPrimeI = _isaPrime[isaIIdx];
	assert_leq(isaPrimeI, _isaPrime.size());
	return isaPrimeI;
}

/**
 * Return: < 0 if suffix i is lexicographically less than suffix j; > 0
 * if suffix j is lexicographically greater.
 */
template <typename TStr>
int64_t DifferenceCoverSample<TStr>::breakTie(uint32_t i, uint32_t j) const {
	assert(built());
	assert_neq(i, j);
	assert_lt(i, this->text().length());
	assert_lt(j, this->text().length());
	uint32_t imodv = this->modv(i);
	uint32_t jmodv = this->modv(j);
	assert_neq(0xffffffff, _dInv[imodv]); // must be in the sample
	assert_neq(0xffffffff, _dInv[jmodv]); // must be in the sample
	uint32_t dimodv = _dInv[imodv];
	uint32_t djmodv = _dInv[jmodv];
	uint32_t ioff = this->divv(i);
	uint32_t joff = this->divv(j);
	assert_lt(dimodv+1, _doffs.size());
	assert_lt(djmodv+1, _doffs.size());
	// assert_lt: expected (32024) < (0)
	assert_lt(ioff, _doffs[dimodv+1] - _doffs[dimodv]);
	assert_lt(joff, _doffs[djmodv+1] - _doffs[djmodv]);
	uint32_t isaIIdx = _doffs[dimodv] + ioff;
	uint32_t isaJIdx = _doffs[djmodv] + joff;
	assert_lt(isaIIdx, _isaPrime.size());
	assert_lt(isaJIdx, _isaPrime.size());
	assert_neq(isaIIdx, isaJIdx); // ranks must be unique
	uint32_t isaPrimeI = _isaPrime[isaIIdx];
	uint32_t isaPrimeJ = _isaPrime[isaJIdx];
	assert_neq(isaPrimeI, isaPrimeJ); // ranks must be unique
	assert_leq(isaPrimeI, _isaPrime.size());
	assert_leq(isaPrimeJ, _isaPrime.size());
	return (int64_t)isaPrimeI - (int64_t)isaPrimeJ;
}

/**
 * Given i, j, return the number of additional characters that need to
 * be compared before the difference cover can break the tie.
 */
template <typename TStr>
uint32_t DifferenceCoverSample<TStr>::tieBreakOff(uint32_t i, uint32_t j) const {
	const TStr& t = this->text();
	const EList<uint32_t>& dmap = this->dmap();
	assert(built());
	// It's actually convenient to allow this, but we're permitted to
	// return nonsense in that case
	if(t[i] != t[j]) return 0xffffffff;
	//assert_eq(t[i], t[j]); // if they're unequal, there's no tie to break
	uint32_t v = this->v();
	assert_neq(i, j);
	assert_lt(i, t.length());
	assert_lt(j, t.length());
	uint32_t imod = this->modv(i);
	uint32_t jmod = this->modv(j);
	uint32_t diffLeft = (jmod >= imod)? (jmod - imod) : (jmod + v - imod);
	uint32_t diffRight = (imod >= jmod)? (imod - jmod) : (imod + v - jmod);
	assert_lt(diffLeft, dmap.size());
	assert_lt(diffRight, dmap.size());
	uint32_t destLeft = dmap[diffLeft];   // offset where i needs to be
	uint32_t destRight = dmap[diffRight]; // offset where i needs to be
	assert(isCovered(destLeft));
	assert(isCovered(destLeft+diffLeft));
	assert(isCovered(destRight));
	assert(isCovered(destRight+diffRight));
	assert_lt(destLeft, v);
	assert_lt(destRight, v);
	uint32_t deltaLeft = (destLeft >= imod)? (destLeft - imod) : (destLeft + v - imod);
	if(deltaLeft == v) deltaLeft = 0;
	uint32_t deltaRight = (destRight >= jmod)? (destRight - jmod) : (destRight + v - jmod);
	if(deltaRight == v) deltaRight = 0;
	assert_lt(deltaLeft, v);
	assert_lt(deltaRight, v);
	assert(isCovered(i+deltaLeft));
	assert(isCovered(j+deltaLeft));
	assert(isCovered(i+deltaRight));
	assert(isCovered(j+deltaRight));
	return min(deltaLeft, deltaRight);
}

#endif /*DIFF_SAMPLE_H_*/
