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
  $Id: basic_iterator_simple.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_ITERATOR_SIMPLE_H
#define SEQAN_HEADER_BASIC_ITERATOR_SIMPLE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Iter
//////////////////////////////////////////////////////////////////////////////

struct SimpleIterator;

/**
.Spec.SimpleIterator:
..cat:Iterators
..summary:A simple iterator.
..signature:Iter<TContainer, SimpleIterator>
..param.TContainer:Type of the container that can be iterated.
...metafunction:Metafunction.Container
..general:Class.Iter
*/
template <typename TContainer>
class Iter<TContainer, SimpleIterator>
{
public:
	typedef typename Value<TContainer>::Type TValue;
	TValue * data_ptr;

	Iter()
	{
	}
	Iter(Iter const & other_):
		data_ptr(other_.data_ptr)
	{
	}
	Iter(TValue * other_data_ptr):
		data_ptr(other_data_ptr)
	{
	}
	template <typename TContainer2>
	Iter(Iter<TContainer2, SimpleIterator> const & other_):
		data_ptr(other_.data_ptr)
	{
	}
	~Iter()
	{
	}
	Iter const &
	operator = (Iter const & other_)
	{
		this->data_ptr = other_.data_ptr;
		return *this;
	}
	Iter const &
	operator = (TValue * other_data_ptr)
	{
		data_ptr = other_data_ptr;
		return *this;
	}
	template <typename TContainer2>
	Iter const &
	operator = (Iter<TContainer2, SimpleIterator> const & other_)
	{
		this->data_ptr = other_.data_ptr;
		return *this;
	}

	operator TValue * ()
	{
		return data_ptr;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Iterator_Default_Imp<T, Standard>
{
	typedef typename Value<T>::Type * Type;
//	typedef Iter<T, SimpleIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TContainer2>
inline typename Position<Iter<TContainer, SimpleIterator> const>::Type 
position(Iter<TContainer, SimpleIterator> const & me,
		 TContainer2 const & cont)
{
SEQAN_CHECKPOINT
	return me.data_ptr - begin(cont);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
