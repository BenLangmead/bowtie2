#ifndef STR_UTIL_H_
#define STR_UTIL_H_

#include <string>

/**
 * Given a string, return an int hash for it.
 */
static inline int
hash_string(const std::string& s) {
	int ret = 0;
	int a = 63689;
	int b = 378551;
	for(size_t i = 0; i < s.length(); i++) {
		ret = (ret * a) + (int)s[i];
		if(a == 0) {
			a += b;
		} else {
			a *= b;
		}
		if(a == 0) {
			a += b;
		}
	}
	return ret;
}

#endif /* STR_UTIL_H_ */
