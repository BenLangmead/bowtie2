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

#ifndef SEQUENCE_IO_H_
#define SEQUENCE_IO_H_

#include <string>
#include <stdexcept>
#include <fstream>
#include <stdio.h>
#include "assert_helpers.h"
#include "ds.h"
#include "filebuf.h"
#include "sstring.h"

using namespace std;

/**
 * Parse the fasta file 'infile'.  Store 
 */
template<typename TFnStr>
static void parseFastaLens(
	const TFnStr&  infile,   // filename
	EList<size_t>& namelens, // destination for fasta name lengths
	EList<size_t>& seqlens)  // destination for fasta sequence lengths
{
	FILE *in = fopen(sstr_to_cstr(infile), "r");
	if(in == NULL) {
		cerr << "Could not open sequence file" << endl;
		throw 1;
	}
	FileBuf fb(in);
	while(!fb.eof()) {
		namelens.expand(); namelens.back() = 0;
		seqlens.expand();  seqlens.back() = 0;
		fb.parseFastaRecordLength(namelens.back(), seqlens.back());
		if(seqlens.back() == 0) {
			// Couldn't read a record.  We're probably done with this file.
			namelens.pop_back();
			seqlens.pop_back();
			continue;
		}
	}
	fb.close();
}

/**
 * Parse the fasta file 'infile'.  Store each name record in 'names', each
 * sequence record  in 'seqs', and the lengths of each 
 */
template<typename TFnStr, typename TNameStr, typename TSeqStr>
static void parseFasta(
	const TFnStr&    infile,   // filename
	EList<TNameStr>& names,    // destination for fasta names
	EList<size_t>&   namelens, // destination for fasta name lengths
	EList<TSeqStr>&  seqs,     // destination for fasta sequences
	EList<size_t>&   seqlens)  // destination for fasta sequence lengths
{
	assert_eq(namelens.size(), seqlens.size());
	assert_eq(names.size(),    namelens.size());
	assert_eq(seqs.size(),     seqlens.size());
	size_t cur = namelens.size();
	parseFastaLens(infile, namelens, seqlens);
	FILE *in = fopen(sstr_to_cstr(infile), "r");
	if(in == NULL) {
		cerr << "Could not open sequence file" << endl;
		throw 1;
	}
	FileBuf fb(in);
	while(!fb.eof()) {
		// Add a new empty record to the end
		names.expand();
		seqs.expand();
		names.back() = new char[namelens[cur]+1];
		seqs.back() = new char[seqlens[cur]+1];
		fb.parseFastaRecord(names.back(), seqs.back());
		if(seqs.back().empty()) {
			// Couldn't read a record.  We're probably done with this file.
			names.pop_back();
			seqs.pop_back();
			continue;
		}
	}
	fb.close();
}

/**
 * Read a set of FASTA sequence files of the given format and alphabet type.
 * Store all of the extracted sequences in vector ss.
 */
template <typename TFnStr, typename TNameStr, typename TSeqStr>
static void parseFastas(
	const EList<TFnStr>& infiles, // filenames
	EList<TNameStr>& names,    // destination for fasta names
	EList<size_t>&   namelens, // destination for fasta name lengths
	EList<TSeqStr>&  seqs,     // destination for fasta sequences
	EList<size_t>&   seqlens)  // destination for fasta sequence lengths
{
	for(size_t i = 0; i < infiles.size(); i++) {
		parseFasta<TFnStr, TNameStr, TSeqStr>(
			infiles[i],
			names,
			namelens,
			seqs,
			seqlens);
	}
}

#endif /*SEQUENCE_IO_H_*/
