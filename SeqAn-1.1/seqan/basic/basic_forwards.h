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
  $Id: basic_forwards.h,v 1.1 2008/08/25 16:20:02 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_FORWARD2_H
#define SEQAN_HEADER_BASIC_FORWARD2_H

//forward declarations (make GCC 4.x happy)

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// basic_transport.h::assign

template <typename TTarget, typename TSource>
inline void
assign(TTarget & target,
	   TSource & source);

template <typename TTarget, typename TSource>
inline void
assign(TTarget & target,
	   TSource const & source);

//////////////////////////////////////////////////////////////////////////////
// string_pointer.h::assignValue

template <typename TValue, typename TPos>
inline void
assignValue(TValue * me,
			TPos pos, 
			TValue const & _value);

//////////////////////////////////////////////////////////////////////////////
// string_pointer.h::moveValue

template <typename TValue, typename TPos>
inline void
moveValue(TValue * me,
			TPos pos, 
			TValue const & _value);

//////////////////////////////////////////////////////////////////////////////
// string_pointer.h::value

template <typename TValue, typename TPos>
inline TValue &
value(TValue * me,
	  TPos pos);

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
