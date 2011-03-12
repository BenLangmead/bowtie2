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
  $Id: index_manual_forwards.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_MANUAL_FORWARDS_H 
#define SEQAN_HEADER_INDEX_MANUAL_FORWARDS_H 

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

//////////////////////////////////////////////////////////////////////////////
// CLASSES
//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN {

	struct _Fibre_Text;		// Original text. Can be a String or a StringSet
	struct _Fibre_RawText;	// Concatenation of the strings above
	struct _Fibre_SA;		// suffix array (of raw text with virtual $-delimiters) with Pair entries
	struct _Fibre_RawSA;	// suffix array with integer entries
	struct _Fibre_SAE;		// suffix array reordered in a b-tree
	struct _Fibre_LCP;		// lcp table of raw text
	struct _Fibre_LCPE;		// lcp interval tree
	struct _Fibre_ChildTab;	// childtab (Kurtz et al.) of raw text
	struct _Fibre_BWT;		// burrows wheeler table of raw text

	typedef Tag<_Fibre_Text> const		Fibre_Text;
	typedef Tag<_Fibre_RawText> const	Fibre_RawText;
	typedef Tag<_Fibre_SA> const		Fibre_SA;
	typedef Tag<_Fibre_RawSA> const		Fibre_RawSA;
	typedef Tag<_Fibre_SAE> const		Fibre_SAE;
	typedef Tag<_Fibre_LCP> const		Fibre_LCP;
	typedef Tag<_Fibre_LCPE> const		Fibre_LCPE;
	typedef Tag<_Fibre_ChildTab> const	Fibre_ChildTab;
	typedef Tag<_Fibre_BWT> const		Fibre_BWT;

}

#endif

