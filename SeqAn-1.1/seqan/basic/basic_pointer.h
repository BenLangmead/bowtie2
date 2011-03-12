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
  $Id: basic_pointer.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_POINTER_H
#define SEQAN_HEADER_BASIC_POINTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Adaption.char array

template <typename TValue>
struct Value< TValue * >
{
	typedef TValue Type;
};
template <typename TValue>
struct Value< TValue * const>
{
	typedef TValue Type;
};

//The next two metafunctions dont work in VC++ due to a compiler bug.
//(the default implementation in common_type.h is called instead)
//work-around: convert arrays to pointers.
template <typename TValue, size_t SIZE>
struct Value< TValue [SIZE] >
{
	typedef TValue Type;
};
template <typename TValue, size_t SIZE>
struct Value< TValue const [SIZE] >
{
	typedef TValue const Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Adaption.char array

template <typename TValue>
struct Iterator< TValue *, Standard>
{
	typedef TValue * Type;
};
template <typename TValue>
struct Iterator< TValue * const, Standard>
{
	typedef TValue * Type;
};

//____________________________________________________________________________

template <typename TValue, size_t SIZE>
struct Iterator< TValue [SIZE], Standard>:
	Iterator<TValue *, Standard>
{
};
template <typename TValue, size_t SIZE>
struct Iterator< TValue const [SIZE], Standard>:
	Iterator<TValue const *, Standard>
{
};

template <typename TValue, size_t SIZE>
struct Iterator< TValue [SIZE], Rooted>:
	Iterator<TValue *, Rooted>
{
};
template <typename TValue, size_t SIZE>
struct Iterator< TValue const [SIZE], Rooted>:
	Iterator<TValue const *, Rooted>
{
};




//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
