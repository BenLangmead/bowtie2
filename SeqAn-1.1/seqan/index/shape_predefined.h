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
  $Id: shape_predefined.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SHAPE_PREDEFINED_H
#define SEQAN_HEADER_SHAPE_PREDEFINED_H

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// some predefined gapped shapes


	//////////////////////////////////////////////////////////////////////////////
	// Single seed of
	// B.Ma and J.Tromp and M.Li, 
	// "PatternHunter: faster and more sensitive homology search"
	// Bioinformatics 18, 2002
	//
	// weight:11 
	// length:18
	// 
	// shape:
	// 111010010100110111

	typedef FixedGappedShape< 
		HardwiredShape< 1, 1, 2, 3, 2, 3, 1, 2, 1, 1 > 
	> ShapePatternHunter;



	//////////////////////////////////////////////////////////////////////////////
	// Multiple seeds of
	// L.Ilie and S.Ilie, "Fast Computation of Good Multiple Spaced Seeds"
	// WABI, 2007
	//
	// weight:9 
	// length:15 
	//
	// shapes:
	// 111010100100111
	// 110100110011101
	// 111010001011011
	//
	// sensitivity:
	// 65% 0.747975		70% 0.897741
	// 75% 0.973134		80% 0.996226

	typedef FixedGappedShape< 
		HardwiredShape< 1, 1, 2, 2, 3, 3, 1, 1 > 
	> ShapeIlie_9_15_1;

	typedef FixedGappedShape< 
		HardwiredShape< 1, 2, 3, 1, 3, 1, 1, 2 > 
	> ShapeIlie_9_15_2;

	typedef FixedGappedShape< 
		HardwiredShape< 1, 1, 2, 4, 2, 1, 2, 1 > 
	> ShapeIlie_9_15_3;



	//////////////////////////////////////////////////////////////////////////////
	// Multiple seeds of
	// L.Ilie and S.Ilie, "Fast Computation of Good Multiple Spaced Seeds"
	// WABI 2007
	//
	// weight:9 
	// length:13..23 
	//
	// shapes:
	// 1110110100111
	// 11010000110010111
	// 11100010010000101011
	//
	// sensitivity:
	// 65% 0.767413		70% 0.910949
	// 75% 0.978558		80% 0.997357

	typedef FixedGappedShape< 
		HardwiredShape< 1, 1, 2, 1, 2, 3, 1, 1 > 
	> ShapeIlie_9_1323_1;

	typedef FixedGappedShape< 
		HardwiredShape< 1, 2, 5, 1, 3, 2, 1, 1 > 
	> ShapeIlie_9_1323_2;

	typedef FixedGappedShape< 
		HardwiredShape< 1, 1, 4, 3, 5, 2, 2, 1 > 
	> ShapeIlie_9_1323_3;


}	// namespace seqan

#endif
