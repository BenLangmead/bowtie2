#ifndef ZBOX_H_
#define ZBOX_H_

/**
 * Fill z with Z-box information for s.  String z will not be resized
 * and will only be filled up to its size cap.  This is the linear-time
 * algorithm from Gusfield.  An optional sanity-check uses a naive
 * algorithm to double-check results.
 */
template<typename T>
void calcZ(const T& s,
           uint32_t off,
           String<uint32_t>& z,
           bool verbose = false,
           bool sanityCheck = false)
{
	size_t lCur = 0, rCur = 0;
	size_t zlen = length(z);
	size_t slen = length(s);
	assert_gt(zlen, 0);
	assert_eq(z[0], 0);
	//assert_leq(zlen, slen);
	for (size_t k = 1; k < zlen && k+off < slen; k++) {
		assert_lt(lCur, k);
		assert(z[lCur] == 0 || z[lCur] == rCur - lCur + 1);
		if(k > rCur) {
			// compare starting at k with prefix starting at 0
			size_t ki = k;
			while(off+ki < length(s) && s[off+ki] == s[off+ki-k]) ki++;
			z[k] = ki - k;
			assert_lt(off+z[k], slen);
			if(z[k] > 0) {
				lCur = k;
				rCur = k + z[k] - 1;
			}
		} else {
			// position k is contained in a Z-box
			size_t betaLen = rCur - k + 1;
			size_t kPrime = k - lCur;
			assert_eq(s[off+k], s[off+kPrime]);
			if(z[kPrime] < betaLen) {
				z[k] = z[kPrime];
				assert_lt(off+z[k], slen);
				// lCur, rCur unchanged
			} else if (z[kPrime] > 0) {
				int q = 0;
				while (off+q+rCur+1 < length(s) && s[off+q+rCur+1] == s[off+betaLen+q]) q++;
				z[k] = betaLen + q;
				assert_lt(off+z[k], slen);
				rCur = rCur + q;
				assert_geq(k, lCur);
				lCur = k;
			} else {
				z[k] = 0;
				assert_lt(off+z[k], slen);
				// lCur, rCur unchanged
			}
		}
	}
#ifndef NDEBUG
	if(sanityCheck) {
		// Recalculate Z-boxes using naive quadratic-time algorithm and
		// compare to linear-time result
		assert_eq(0, z[0]);
		for(size_t i = 1; i < length(z); i++) {
			size_t j;
			for(j = i; off+j < length(s); j++) {
				if(s[off+j] != s[off+j-i]) break;
			}
			assert_eq(j-i, z[i]);
		}
	}
#endif
}

#endif /*ZBOX_H_*/
