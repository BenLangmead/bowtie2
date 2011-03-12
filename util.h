/*
 * util.h
 */

#ifndef UTIL_H_
#define UTIL_H_

/**
 * C++ version char* style "itoa":
 */
template<typename T>
char* itoa10(const T& value, char* result) {
	// Check that base is valid
	char* out = result;
	T quotient = value;
	if(quotient < 0) quotient = -quotient;
	// Now write each digit from most to least significant
	do {
		*out = "0123456789"[quotient % 10];
		++out;
		quotient /= 10;
	} while (quotient > 0);
	// Only apply negative sign for base 10
	if (value < 0) *out++ = '-';
	std::reverse( result, out );
	*out = 0; // terminator
	return out;
}

#endif /*ndef UTIL_H_*/
