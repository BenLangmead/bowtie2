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
  $Id: sequence_forwards.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEQUENCE_FORWARDS_H 
#define SEQAN_HEADER_SEQUENCE_FORWARDS_H 

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN 
{
    
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData> 
void read(TFile & file, TData & data);       	// "projects/library/seqan/file/file_format_raw.h"(307)

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData> 
void write(TFile & file, TData & data);       	// "projects/library/seqan/file/file_format_raw.h"(327)

template <typename TFile, typename TData> 
void write(TFile & file, TData const & data);   // "projects/library/seqan/file/file_format_raw.h"(335)


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif

