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

#ifdef MAIN_SSTRING

#include <string.h>
#include <iostream>
#include "ds.h"
#include "sstring.h"

using namespace std;

int main(void) {
	cerr << "Test inter-class comparison operators...";
	{
		SString<int> s(2);
		s.set('a', 0);
		s.set('b', 1);
		assert(sstr_eq(s, (const char *)"ab"));
		assert(!sstr_neq(s, (const char *)"ab"));
		assert(!sstr_lt(s, (const char *)"ab"));
		assert(!sstr_gt(s, (const char *)"ab"));
		assert(sstr_leq(s, (const char *)"ab"));
		assert(sstr_geq(s, (const char *)"ab"));
		
		SStringExpandable<int> s2;
		s2.append('a');
		s2.append('b');
		assert(sstr_eq(s, s2));
		assert(sstr_eq(s2, (const char *)"ab"));
		assert(!sstr_neq(s, s2));
		assert(!sstr_neq(s2, (const char *)"ab"));
		assert(!sstr_lt(s, s2));
		assert(!sstr_lt(s2, (const char *)"ab"));
		assert(!sstr_gt(s, s2));
		assert(!sstr_gt(s2, (const char *)"ab"));
		assert(sstr_leq(s, s2));
		assert(sstr_leq(s2, (const char *)"ab"));
		assert(sstr_geq(s, s2));
		assert(sstr_geq(s2, (const char *)"ab"));

		SStringFixed<int, 12> s3;
		s3.append('a');
		s3.append('b');
		assert(sstr_eq(s, s3));
		assert(sstr_eq(s2, s3));
		assert(sstr_eq(s3, (const char *)"ab"));
		assert(!sstr_neq(s, s3));
		assert(!sstr_neq(s2, s3));
		assert(!sstr_neq(s3, (const char *)"ab"));
		assert(!sstr_lt(s, s3));
		assert(!sstr_lt(s2, s3));
		assert(!sstr_lt(s3, (const char *)"ab"));
		assert(!sstr_gt(s, s3));
		assert(!sstr_gt(s2, s3));
		assert(!sstr_gt(s3, (const char *)"ab"));
		assert(sstr_geq(s, s3));
		assert(sstr_geq(s2, s3));
		assert(sstr_geq(s3, (const char *)"ab"));
		assert(sstr_leq(s, s3));
		assert(sstr_leq(s2, s3));
		assert(sstr_leq(s3, (const char *)"ab"));
	}
	cerr << "PASSED" << endl;
	
	cerr << "Test flag for whether to consider end-of-word < other chars ...";
	{
		SString<char> ss("String");
		SString<char> sl("String1");
		assert(sstr_lt(ss, sl));
		assert(sstr_gt(ss, sl, false));
		assert(sstr_leq(ss, sl));
		assert(sstr_geq(ss, sl, false));
	}
	cerr << "PASSED" << endl;
	
	cerr << "Test toZBuf and toZBufXForm ...";
	{
		SString<uint32_t> s(10);
		for(int i = 0; i < 10; i++) {
			s[i] = (uint32_t)i;
		}
		assert(strcmp(s.toZBufXForm("0123456789"), "0123456789") == 0);
	}
	cerr << "PASSED" << endl;

	cerr << "Test S2bDnaString ...";
	{
		const char *str =
			"ACGTACGTAC" "ACGTACGTAC" "ACGTACGTAC"
			"ACGTACGTAC" "ACGTACGTAC" "ACGTACGTAC";
		const char *gs =
			"GGGGGGGGGG" "GGGGGGGGGG" "GGGGGGGGGG"
			"GGGGGGGGGG" "GGGGGGGGGG" "GGGGGGGGGG";
		for(size_t i = 0; i < 60; i++) {
			S2bDnaString s(str, i, true);
			S2bDnaString sr;
			BTDnaString s2(str, i, true);
			assert(sstr_eq(s, s2));
			if(i >= 10) {
				BTDnaString s3;
				s.windowGetDna(s3, true, false, 3, 4);
				assert(sstr_eq(s3.toZBuf(), (const char*)"TACG"));
				s.windowGetDna(s3, false, false, 3, 4);
				assert(sstr_eq(s3.toZBuf(), (const char*)"CGTA"));
				assert_eq('A', s.toChar(0));
				assert_eq('G', s.toChar(2));
				assert_eq('A', s.toChar(4));
				assert_eq('G', s.toChar(6));
				assert_eq('A', s.toChar(8));
				
				s.reverseWindow(1, 8);
				s2.reverseWindow(1, 8);
				
				assert_eq('A', s.toChar(1));
				assert_eq('T', s.toChar(2));
				assert_eq('G', s.toChar(3));
				assert_eq('C', s.toChar(4));
				assert_eq('A', s.toChar(5));
				assert_eq('T', s.toChar(6));
				assert_eq('G', s.toChar(7));
				assert_eq('C', s.toChar(8));
				assert(sstr_eq(s, s2));

				s.reverseWindow(1, 8);
				s2.reverseWindow(1, 8);
				assert(sstr_eq(s, s2));
			}
			if(i > 1) {
				s.reverse();
				sr.installReverseChars(str, i);
				s2.reverse();
				assert(sstr_eq(s, s2));
				assert(sstr_eq(sr, s2));
				s.reverse();
				sr.reverse();
				assert(sstr_neq(s, s2));
				assert(sstr_neq(sr, s2));
				s.fill(2);
				s2.reverse();
				assert(sstr_leq(s, gs));
				assert(sstr_gt(s, s2));
				assert(sstr_gt(s, sr));
				s2.fill(2);
				sr.fill(2);
				assert(sstr_eq(s, s2));
				assert(sstr_eq(s, sr));
			}
		}
		S2bDnaString s(str, true);
		S2bDnaString sr;
		BTDnaString s2(str, true);
		assert(sstr_eq(s2.toZBuf(), str));
		assert(sstr_eq(s, s2));
		s.reverse();
		sr.installReverseChars(str);
		s2.reverse();
		assert(sstr_eq(s, s2));
		assert(sstr_eq(sr, s2));
		s.reverse();
		sr.reverse();
		assert(sstr_neq(s, s2));
		assert(sstr_neq(sr, s2));
	}
	cerr << "PASSED" << endl;

	cerr << "Test operator=() ...";
	{
		S2bDnaString s;
		s.installChars(string("gtcagtca"));
		assert(sstr_eq(s.toZBuf(), (const char *)"GTCAGTCA"));
	}
	cerr << "PASSED" << endl;
	
	cerr << "Conversions from string ...";
	{
		SStringExpandable<char> se(string("hello"));
		EList<SStringExpandable<char> > sel;
		sel.push_back(SStringExpandable<char>(string("hello")));
	}
	cerr << "PASSED" << endl;
	
	cerr << "PASSED" << endl;
}

#endif /*def MAIN_SSTRING*/
