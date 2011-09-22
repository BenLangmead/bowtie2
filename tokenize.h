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

#ifndef TOKENIZE_H_
#define TOKENIZE_H_

#include <string>
#include <sstream>

using namespace std;

/**
 * Split string s according to given delimiters.  Mostly borrowed
 * from C++ Programming HOWTO 7.3.
 */
template<typename T>
static inline void tokenize(const string& s, const string& delims,
                            T& ss, size_t max = 9999)
{
	//string::size_type lastPos = s.find_first_not_of(delims, 0);
	string::size_type lastPos = 0;
	string::size_type pos = s.find_first_of(delims, lastPos);
	while (string::npos != pos || string::npos != lastPos) {
		ss.push_back(s.substr(lastPos, pos - lastPos));
		lastPos = s.find_first_not_of(delims, pos);
		pos = s.find_first_of(delims, lastPos);
		if(ss.size() == (max - 1)) {
			pos = string::npos;
		}
	}
}

template<typename T>
static inline void tokenize(const std::string& s, char delim, T& ss) {
	std::string token;
	std::istringstream iss(s);
	while(getline(iss, token, delim)) {
		ss.push_back(token);
	}
}

#endif /*TOKENIZE_H_*/
