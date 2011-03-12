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
  $Id: basic_compare.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_COMPARE_H
#define SEQAN_HEADER_BASIC_COMPARE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename TLeft, typename TRight>
struct CompareType;

template <typename TLeft, typename TRight>
struct CompareType<TLeft const, TRight>
{
	typedef typename CompareType<TLeft, TRight>::Type const Type;
};
template <typename TLeft, typename TRight>
struct CompareType<TLeft, TRight const>
{
	typedef typename CompareType<TLeft, TRight>::Type const Type;
};
template <typename TLeft, typename TRight>
struct CompareType<TLeft const, TRight const>
{
	typedef typename CompareType<TLeft, TRight>::Type const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename _T> inline
bool lexLess(const _T& _Left, const _T& _Right)
{	// return lexicographical _Left < _Right
	typedef typename _MakeUnsigned<_T>::Type TUnsigned;
    return (TUnsigned)_Left < (TUnsigned)_Right;
}

//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
