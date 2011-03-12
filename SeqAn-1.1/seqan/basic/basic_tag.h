 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: basic_tag.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_TAG_H
#define SEQAN_HEADER_BASIC_TAG_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

/**
.Tag.DotDrawing
..summary:Switch to trigger drawing in dot format.
..value.DotDrawing:Graphs in dot format.
*/

struct DotDrawing_;
typedef Tag<DotDrawing_> const DotDrawing;


/**
.Tag.HammingDistance
..summary:Switch to trigger Hamming distance, which is a measure of character substitutions.
*/

/**
.Tag.LevenshteinDistance
..summary:Switch to trigger Levenshtein distance, which is a measure of edit operations (character substitutions, deletions or insertions).
*/

struct _HammingDistance;
struct _LevenshteinDistance;

typedef Tag<_HammingDistance>		HammingDistance;
typedef Tag<_LevenshteinDistance>	LevenshteinDistance;
typedef Tag<_LevenshteinDistance>	EditDistance;


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Alignment: Tags
//////////////////////////////////////////////////////////////////////////////
//Sollte eigentlich nach align/, aber da jetzt ja so viele
//alignment algorithmen in graph/ gelandet sind...

/**
.Tag.Global Alignment Algorithms:
..summary:Global alignment algorithm used by globalAlignment.
..see:Function.globalAlignment
..see:Tag.Local Alignment Algorithms
*/

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.NeedlemanWunsch:
	Dynamic programming algorithm for alignments by Needleman and Wunsch.
*/

struct NeedlemanWunsch_;
typedef Tag<NeedlemanWunsch_> const NeedlemanWunsch;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.Gotoh:
	Gotoh's affine gap cost alignment algorithm.
*/
struct Gotoh_;
typedef Tag<Gotoh_> const Gotoh;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.MyersBitVector:
	Myers' bit vector alignment algorithm for edit distance.
	Note that this algorithm does not returns the alignment itself, but only computes the score.
*/
struct MyersBitVector_;
typedef Tag<MyersBitVector_> const MyersBitVector;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.MyersHirschberg:
	Myers' bit vector algorithm for edit distance combined with Hirschberg's linear space alignment algorithm.
*/
struct MyersHirschberg_;
typedef Tag<MyersHirschberg_> const MyersHirschberg;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.Hirschberg:
	Hirschberg's linear space global alignment algorithm.
*/
struct Hirschberg_;
typedef Tag<Hirschberg_> const Hirschberg;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms:
..summary:Local alignment algorithm used by localAlignment.
..see:Function.localAlignment
*/

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms.value.SmithWaterman:
	Triggers a Smith Waterman local alignment algorithm.
*/
struct SmithWaterman_;
typedef Tag<SmithWaterman_> const SmithWaterman;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms.value.SmithWatermanClump:
	Local alignment algorithm with "declumping" by Waterman and Eggert.
*/
struct SmithWatermanClump_;
typedef Tag<SmithWatermanClump_> const SmithWatermanClump;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms.value.SmithWatermanIsland:
	Local alignment algorithm finding all "islands".
*/
struct SmithWatermanIsland_;
typedef Tag<SmithWatermanIsland_> const SmithWatermanIsland;

//////////////////////////////////////////////////////////////////////////////

/*DISABLED
.Tag.RNA Folding Algorithms.value.Nussinov:
	Nussinov style RNA folding algorithm
*/
struct Nussinov_;
typedef Tag<Nussinov_> const Nussinov;

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
