/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
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

#ifndef ALPHABETS_H_
#define ALPHABETS_H_

#include <stdexcept>
#include <string>
#include <sstream>
#include <stdint.h>
#include "assert_helpers.h"

using namespace std;

/// Convert an ascii char to a DNA category.  Categories are:
/// 0 -> invalid
/// 1 -> unambiguous a, c, g or t
/// 2 -> ambiguous
/// 3 -> unmatchable
extern uint8_t asc2dnacat[];
/// Convert masks to ambiguous nucleotides
extern char mask2dna[];
/// Convert ambiguous ASCII nuceleotide to mask
extern uint8_t asc2dnamask[];
/// Convert mask to # of alternative in the mask
extern int mask2popcnt[];
/// Convert an ascii char to a 2-bit base: 0=A, 1=C, 2=G, 3=T, 4=N
extern uint8_t asc2dna[];
/// Convert an ascii char representing a base or a color to a 2-bit
/// code: 0=A,0; 1=C,1; 2=G,2; 3=T,3; 4=N,.
extern uint8_t asc2dnaOrCol[];
/// Convert a pair of DNA masks to a color mask
extern uint8_t dnamasks2colormask[16][16];

/// Convert an ascii char to a color category.  Categories are:
/// 0 -> invalid
/// 1 -> unambiguous 0, 1, 2 or 3
/// 2 -> ambiguous (not applicable for colors)
/// 3 -> unmatchable
extern uint8_t asc2colcat[];
/// Convert an ascii char to a 2-bit base: 0=A, 1=C, 2=G, 3=T, 4=N
extern uint8_t asc2col[];
/// Convert an ascii char to its DNA complement, including IUPACs
extern char asc2dnacomp[];

/// Convert a pair of 2-bit (and 4=N) encoded DNA bases to a color
extern uint8_t dinuc2color[5][5];
/// Convert a 2-bit nucleotide (and 4=N) and a color to the
/// corresponding 2-bit nucleotide
extern uint8_t nuccol2nuc[5][5];
/// Convert a 4-bit mask into an IUPAC code
extern char mask2iupac[16];

/// Convert an ascii color to an ascii dna char
extern char col2dna[];
/// Convert an ascii dna to a color char
extern char dna2col[];
/// Convert an ascii dna to a color char
extern const char* dna2colstr[];

/// Convert bit encoded DNA char to its complement
extern int dnacomp[5];

/// String of all DNA and IUPAC characters
extern const char *iupacs;

/// Map from masks to their reverse-complement masks
extern int maskcomp[16];

/**
 * Return true iff c is a Dna character.
 */
static inline bool isDna(char c) {
	return asc2dnacat[(int)c] > 0;
}

/**
 * Return true iff c is a color character.
 */
static inline bool isColor(char c) {
	return asc2colcat[(int)c] > 0;
}

/**
 * Return true iff c is an ambiguous Dna character.
 */
static inline bool isAmbigNuc(char c) {
	return asc2dnacat[(int)c] == 2;
}

/**
 * Return true iff c is an ambiguous color character.
 */
static inline bool isAmbigColor(char c) {
	return asc2colcat[(int)c] == 2;
}

/**
 * Return true iff c is an ambiguous character.
 */
static inline bool isAmbig(char c, bool color) {
	return (color ? asc2colcat[(int)c] : asc2dnacat[(int)c]) == 2;
}

/**
 * Return true iff c is an unambiguous DNA character.
 */
static inline bool isUnambigNuc(char c) {
	return asc2dnacat[(int)c] == 1;
}

/**
 * Return the DNA complement of the given ASCII char.
 */
static inline char comp(char c) {
	switch(c) {
	case 'a': return 't';
	case 'A': return 'T';
	case 'c': return 'g';
	case 'C': return 'G';
	case 'g': return 'c';
	case 'G': return 'C';
	case 't': return 'a';
	case 'T': return 'A';
	default: return c;
	}
}

/**
 * Return the reverse complement of a bit-encoded nucleotide.
 */
static inline int compDna(int c) {
	assert_leq(c, 4);
	return dnacomp[c];
}

/**
 * Return true iff c is an unambiguous Dna character.
 */
static inline bool isUnambigDna(char c) {
	return asc2dnacat[(int)c] == 1;
}

/**
 * Return true iff c is an unambiguous color character (0,1,2,3).
 */
static inline bool isUnambigColor(char c) {
	return asc2colcat[(int)c] == 1;
}

/// Convert a pair of 2-bit (and 4=N) encoded DNA bases to a color
extern uint8_t dinuc2color[5][5];

/**
 * Decode a not-necessarily-ambiguous nucleotide.
 */
static inline void decodeNuc(char c , int& num, int *alts) {
	switch(c) {
	case 'A': alts[0] = 0; num = 1; break;
	case 'C': alts[0] = 1; num = 1; break;
	case 'G': alts[0] = 2; num = 1; break;
	case 'T': alts[0] = 3; num = 1; break;
	case 'M': alts[0] = 0; alts[1] = 1; num = 2; break;
	case 'R': alts[0] = 0; alts[1] = 2; num = 2; break;
	case 'W': alts[0] = 0; alts[1] = 3; num = 2; break;
	case 'S': alts[0] = 1; alts[1] = 2; num = 2; break;
	case 'Y': alts[0] = 1; alts[1] = 3; num = 2; break;
	case 'K': alts[0] = 2; alts[1] = 3; num = 2; break;
	case 'V': alts[0] = 0; alts[1] = 1; alts[2] = 2; num = 3; break;
	case 'H': alts[0] = 0; alts[1] = 1; alts[2] = 3; num = 3; break;
	case 'D': alts[0] = 0; alts[1] = 2; alts[2] = 3; num = 3; break;
	case 'B': alts[0] = 1; alts[1] = 2; alts[2] = 3; num = 3; break;
	case 'N': alts[0] = 0; alts[1] = 1; alts[2] = 2; alts[3] = 3; num = 4; break;
	default: {
		std::cerr << "Bad IUPAC code: " << c << ", (int: " << (int)c << ")" << std::endl;
		throw std::runtime_error("");
	}
	}
}

extern void setIupacsCat(uint8_t cat);

#endif /*ALPHABETS_H_*/
