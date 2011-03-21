//
//  ls.cpp
//

#ifdef MAIN_LS

#include <string.h>
#include <iostream>
#include "sstring.h"
#include "ls.h"
#include "ds.h"

using namespace std;

int main(void) {
	cerr << "Test LarssonSadakana for int...";
	{
		typedef int T;
		const char *t = "banana";
		EList<T> sa;
		EList<T> isa;
		for(size_t i = 0; i < strlen(t); i++) {
			isa.push_back(t[i]);
		}
		isa.push_back(0); // disregarded
		sa.resize(isa.size());
		LarssonSadakane<T> ls;
		ls.suffixsort(isa.ptr(), sa.ptr(), (T)sa.size()-1, 'z', 0);
		assert_eq((T)'a', t[sa[1]]); assert_eq(5, sa[1]);
		assert_eq((T)'a', t[sa[2]]); assert_eq(3, sa[2]);
		assert_eq((T)'a', t[sa[3]]); assert_eq(1, sa[3]);
		assert_eq((T)'b', t[sa[4]]); assert_eq(0, sa[4]);
		assert_eq((T)'n', t[sa[5]]); assert_eq(4, sa[5]);
		assert_eq((T)'n', t[sa[6]]); assert_eq(2, sa[6]);
	}
	cerr << "PASSED" << endl;

	cerr << "Test LarssonSadakana for uint32_t...";
	{
		typedef uint32_t T;
		const char *t = "banana";
		EList<T> sa;
		EList<T> isa;
		for(size_t i = 0; i < strlen(t); i++) {
			isa.push_back(t[i]);
		}
		isa.push_back(0); // disregarded
		sa.resize(isa.size());
		LarssonSadakane<int> ls;
		ls.suffixsort(
			(int*)isa.ptr(),
			(int*)sa.ptr(),
			(int)sa.size()-1,
			'z',
			0);
		assert_eq((T)'a', t[sa[1]]); assert_eq(5, sa[1]);
		assert_eq((T)'a', t[sa[2]]); assert_eq(3, sa[2]);
		assert_eq((T)'a', t[sa[3]]); assert_eq(1, sa[3]);
		assert_eq((T)'b', t[sa[4]]); assert_eq(0, sa[4]);
		assert_eq((T)'n', t[sa[5]]); assert_eq(4, sa[5]);
		assert_eq((T)'n', t[sa[6]]); assert_eq(2, sa[6]);
	}
	cerr << "PASSED" << endl;

	cerr << "Last elt is < or > others ...";
	{
		{
		typedef int T;
		const char *t = "aaa";
		EList<T> sa;
		EList<T> isa;
		for(size_t i = 0; i < strlen(t); i++) {
			isa.push_back(t[i]);
		}
		isa.push_back(0); // disregarded
		sa.resize(isa.size());
		LarssonSadakane<T> ls;
		ls.suffixsort(isa.ptr(), sa.ptr(), (T)sa.size()-1, 'z', 0);
		assert_eq(3, sa[0]);
		assert_eq(2, sa[1]);
		assert_eq(1, sa[2]);
		assert_eq(0, sa[3]);
		}

		{
		typedef int T;
		const char *t = "aaa";
		EList<T> sa;
		EList<T> isa;
		for(size_t i = 0; i < strlen(t); i++) {
			isa.push_back(t[i]);
		}
		isa.push_back('y'); // doesn't matter if this is > others
		sa.resize(isa.size());
		LarssonSadakane<T> ls;
		ls.suffixsort(isa.ptr(), sa.ptr(), (T)sa.size()-1, 'z', 0);
		assert_eq(3, sa[0]);
		assert_eq(2, sa[1]);
		assert_eq(1, sa[2]);
		assert_eq(0, sa[3]);
		}
		
		{
		typedef int T;
		const char *t = "aaa";
		EList<T> sa;
		EList<T> isa;
		for(size_t i = 0; i < strlen(t); i++) {
			isa.push_back(t[i]);
		}
		isa.push_back('y'); // breaks ties
		isa.push_back(0);   // disregarded
		sa.resize(isa.size());
		LarssonSadakane<T> ls;
		ls.suffixsort(isa.ptr(), sa.ptr(), (T)sa.size()-1, 'z', 0);
		assert_eq(4, sa[0]);
		assert_eq(0, sa[1]);
		assert_eq(1, sa[2]);
		assert_eq(2, sa[3]);
		assert_eq(3, sa[4]);
		}
		
	}
	cerr << "PASSED" << endl;
}

#endif
