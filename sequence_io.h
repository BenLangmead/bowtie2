#ifndef SEQUENCE_IO_H_
#define SEQUENCE_IO_H_

#include <string>
#include <stdexcept>
#include <fstream>
#include <stdio.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include "assert_helpers.h"
#include "ds.h"

using namespace std;
using namespace seqan;

/**
 * Read a sequence file of the given format and alphabet type.  Store
 * all of the extracted sequences in vector ss.  Note that SeqAn's
 * policy for when it encounters characters not from the specified
 * alphabet is to convert them to the lexicographically smallest
 * character in the alphabet.
 */
template <typename TStr, typename TFile>
static void readSequenceFile(const std::string& infile,
                             EList<TStr>& ss,
                             int64_t& baseCutoff, // limit for total bases
                             int seqCutoff = -1,  // limit for sequences
                             bool reverse = false)
{
	typedef typename Value<TStr>::Type TVal;
	static char buf[256 * 1024]; // fairly large input buffer
	if(baseCutoff <= 0) return;
	FILE *in = fopen(infile.c_str(), "r");
	if(in == NULL) {
		cerr << "Could not open sequence file" << endl;
		throw 1;
	}
	// Associate large input buffer with FILE *in
	if(setvbuf(in, buf, _IOFBF, 256 * 1024) != 0) {
		cerr << "Could not create input buffer for sequence file" << endl;
		throw 1;
	}
	// Read entries using SeqAn
	int cnt = 0;
	while(!feof(in)) {
		while(true) {
			ss.push_back(TStr()); // add a new empty string to the end
			// Fill the new empty string with the next sequence from
			// the file.  SeqAn allocates just enough mem for it (at
			// the expense of lots of file seeks, which can add up)
			seqan::read(in, ss.back(), TFile()); 
			if(seqan::empty(ss.back())) {
				ss.pop_back();
				break;
			}
			// Enforce the base cutoff
			if((int64_t)length(ss.back()) > baseCutoff) {
				resize(ss.back(), baseCutoff);
				baseCutoff = 0;
			} else {
				baseCutoff -= length(ss.back());
			}
			// Reverse the newly-read sequence in-place if desired
			if(reverse) {
				size_t len = length(ss.back());
				for(size_t i = 0; i < len/2; i++) {
					TVal t = ss.back()[i];
					ss.back()[i] = ss.back()[len-i-1];
					ss.back()[len-i-1] = t;
				}
			}
			#ifndef NDEBUG
			// Sanity check that all (int) values are in range
			for(size_t i = 0; i < length(ss.back()); i++) {
				assert_lt(ss.back()[i], (int)(ValueSize<TVal>::VALUE));
				assert_geq(ss.back()[i], 0);
			}
			#endif
			cnt++;
			// Enforce the sequence cutoff
			if(seqCutoff != -1 && cnt >= seqCutoff) {
				fclose(in);
				return;
			}
		}
	}
	fclose(in);
}

/**
 * Read a set of sequence files of the given format and alphabet type.
 * Store all of the extracted sequences in vector ss.
 */
template <typename TStr, typename TFile>
static void readSequenceFiles(const EList<std::string>& infiles,
                              EList<TStr>& ss,
                              int64_t& baseCutoff,
                              int seqCutoff = -1,
                              bool reverse = false)
{
	for(size_t i = 0; i < infiles.size() && baseCutoff > 0; i++) {
		readSequenceFile<TStr,TFile>(infiles[i], ss, baseCutoff, seqCutoff, reverse);
		if(baseCutoff <= 0) break;
	}
}

/**
 * Read a set of sequence files of the given format and alphabet type.
 * Store all of the extracted sequences in vector ss.
 */
template <typename TStr, typename TFile>
static void readSequenceFiles(const EList<std::string>& infiles,
                              EList<TStr>& ss,
                              int seqCutoff = -1,
                              bool reverse = false)
{
	int64_t i = 0xffffffffll;
	readSequenceFiles<TStr,TFile>(infiles, ss, i, seqCutoff, reverse);
}

/**
 * Parse a comma-delimited list of strings of type T into a vector. 
 */
template <typename T>
void readSequenceString(const std::string& s,
                        EList<T>& ss,
                        int64_t& baseCutoff,
                        int seqCutoff = -1,
                        bool reverse = false)
{
	// Split string s using comma as a delimiter.  Borrowed from C++
	// Programming HOWTO 7.3
	std::string::size_type lastPos = s.find_first_not_of(",", 0);
	std::string::size_type pos = s.find_first_of(",", lastPos);
	while (baseCutoff > 0 && (std::string::npos != pos || std::string::npos != lastPos)) {
		string stmp = s.substr(lastPos, pos - lastPos);
		if((int64_t)stmp.length() < baseCutoff) {
			baseCutoff -= stmp.length();
		} else {
			stmp = stmp.substr(0, baseCutoff);
			baseCutoff = 0;
		}
		if(reverse) {
			size_t len = stmp.length();
			for(size_t i = 0; i < len/2; i++) {
				char tmp = stmp[i];
				stmp[i] = stmp[len-i-1];
				stmp[len-i-1] = tmp;
			}
			ss.push_back(T(stmp.c_str()));
		} else {
			ss.push_back(T(stmp.c_str()));
		}
		if(seqCutoff != -1 && ss.size() >= (size_t)seqCutoff) {
			return;
		}
	    lastPos = s.find_first_not_of(",", pos);
	    pos = s.find_first_of(",", lastPos);
	}
}

/**
 * Parse a comma-delimited list of strings of type T into a vector. 
 * Doesn't require callee to supply a baseCutoff.
 */
template <typename T>
void readSequenceString(const std::string& s,
                        EList<T>& ss,
                        int seqCutoff = -1,
                        bool reverse = false)
{
	int64_t i = 0xffffffffll;
	readSequenceString(s, ss, i, seqCutoff, reverse);
}

#endif /*SEQUENCE_IO_H_*/
