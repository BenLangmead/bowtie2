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

#include <stdint.h>
#include <cassert>
#include <string>
#include "alphabet.h"

using namespace std;

/**
 * Mapping from ASCII characters to DNA categories:
 *
 * 0 = invalid - error
 * 1 = DNA
 * 2 = IUPAC (ambiguous DNA)
 * 3 = not an error, but unmatchable; alignments containing this
 *     character are invalid
 */
uint8_t asc2dnacat[] = {
	/*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0,
	       /*                                        - */
	/*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  64 */ 0, 1, 2, 1, 2, 0, 0, 1, 2, 0, 0, 2, 0, 2, 2, 0,
	       /*    A  B  C  D        G  H        K     M  N */
	/*  80 */ 0, 0, 2, 2, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
	       /*       R  S  T     V  W  X  Y */
	/*  96 */ 0, 1, 2, 1, 2, 0, 0, 1, 2, 0, 0, 2, 0, 2, 2, 0,
	       /*    a  b  c  d        g  h        k     m  n */
	/* 112 */ 0, 0, 2, 2, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
	       /*       r  s  t     v  w  x  y */
	/* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

// 5-bit pop count
int mask2popcnt[] = {
	0, 1, 1, 2, 1, 2, 2, 3,
	1, 2, 2, 3, 2, 3, 3, 4,
	1, 2, 2, 3, 2, 3, 3, 4,
	2, 3, 3, 4, 3, 4, 4, 5
};

/**
 * Mapping from masks to ASCII characters for ambiguous nucleotides.
 */
char mask2dna[] = {
	'?', // 0
	'A', // 1
	'C', // 2
	'M', // 3
	'G', // 4
	'R', // 5
	'S', // 6
	'V', // 7
	'T', // 8
	'W', // 9
	'Y', // 10
	'H', // 11
	'K', // 12
	'D', // 13
	'B', // 14
	'N', // 15 (inclusive N)
	'N'  // 16 (exclusive N)
};

/**
 * Mapping from ASCII characters for ambiguous nucleotides into masks:
 */
uint8_t asc2dnamask[] = {
	/*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  64 */ 0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0,
	       /*    A  B  C  D        G  H        K     M  N */
	/*  80 */ 0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0,
	       /*       R  S  T     V  W     Y */
	/*  96 */ 0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0,
	       /*    a  b  c  d        g  h        k     m  n */
	/* 112 */ 0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0,
	       /*       r  s  t     v  w     y */
	/* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

/**
 * Set the category for all IUPAC codes.  By default they're in
 * category 2 (IUPAC), but sometimes we'd like to put them in category
 * 3 (unmatchable), for example.
 */
void setIupacsCat(uint8_t cat) {
	assert(cat < 4);
	asc2dnacat[(int)'B'] = asc2dnacat[(int)'b'] =
	asc2dnacat[(int)'D'] = asc2dnacat[(int)'d'] =
	asc2dnacat[(int)'H'] = asc2dnacat[(int)'h'] =
	asc2dnacat[(int)'K'] = asc2dnacat[(int)'k'] =
	asc2dnacat[(int)'M'] = asc2dnacat[(int)'m'] =
	asc2dnacat[(int)'N'] = asc2dnacat[(int)'n'] =
	asc2dnacat[(int)'R'] = asc2dnacat[(int)'r'] =
	asc2dnacat[(int)'S'] = asc2dnacat[(int)'s'] =
	asc2dnacat[(int)'V'] = asc2dnacat[(int)'v'] =
	asc2dnacat[(int)'W'] = asc2dnacat[(int)'w'] =
	asc2dnacat[(int)'X'] = asc2dnacat[(int)'x'] =
	asc2dnacat[(int)'Y'] = asc2dnacat[(int)'y'] = cat;
}

/// For converting from ASCII to the Dna5 code where A=0, C=1, G=2,
/// T=3, N=4
/// According to the manual all the other characters, including
/// IUPAC codes are being converted to N
uint8_t asc2dna[] = {
	/*   0 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/*  16 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/*  32 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/*                                               - */
	/*  48 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/*  64 */ 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	       /*    A  B  C  D        G  H        K     M  N */
	/*  80 */ 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	       /*       R  S  T  U  V  W     Y */
	/*  96 */ 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	       /*    a  b  c  d        g  h        k     m  n */
	/* 112 */ 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	       /*       r  s  t  u  v  w     y */
	/* 128 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 144 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 160 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 176 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 192 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 208 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 224 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 240 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};

/// For converting from ASCII to the Dna5 code where A=0, C=1, G=2,
/// T=3, N=4
uint8_t asc2col[] = {
	/*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 0,
	       /*                                        -  . */
	/*  48 */ 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	       /* 0  1  2  3 */
	/*  64 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  80 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  96 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 112 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

/**
 * Convert a pair of nucleotides to a color.
 */
uint8_t dinuc2color[5][5] = {
	/* A */ {0, 1, 2, 3, 4},
	/* C */ {1, 0, 3, 2, 4},
	/* G */ {2, 3, 0, 1, 4},
	/* T */ {3, 2, 1, 0, 4},
	/* N */ {4, 4, 4, 4, 4}
};

/// Convert bit encoded DNA char to its complement
int dnacomp[5] = {
	3, 2, 1, 0, 4
};

const char *iupacs = "!ACMGRSVTWYHKDBN!acmgrsvtwyhkdbn";

char mask2iupac[16] = {
	-1,
	'A', // 0001
	'C', // 0010
	'M', // 0011
	'G', // 0100
	'R', // 0101
	'S', // 0110
	'V', // 0111
	'T', // 1000
	'W', // 1001
	'Y', // 1010
	'H', // 1011
	'K', // 1100
	'D', // 1101
	'B', // 1110
	'N', // 1111
};

int maskcomp[16] = {
	0,  // 0000 (!) -> 0000 (!)
	8,  // 0001 (A) -> 1000 (T)
	4,  // 0010 (C) -> 0100 (G)
	12, // 0011 (M) -> 1100 (K)
	2,  // 0100 (G) -> 0010 (C)
	10, // 0101 (R) -> 1010 (Y)
	6,  // 0110 (S) -> 0110 (S)
	14, // 0111 (V) -> 1110 (B)
	1,  // 1000 (T) -> 0001 (A)
	9,  // 1001 (W) -> 1001 (W)
	5,  // 1010 (Y) -> 0101 (R)
	13, // 1011 (H) -> 1101 (D)
	3,  // 1100 (K) -> 0011 (M)
	11, // 1101 (D) -> 1011 (H)
	7,  // 1110 (B) -> 0111 (V)
	15, // 1111 (N) -> 1111 (N)
};

