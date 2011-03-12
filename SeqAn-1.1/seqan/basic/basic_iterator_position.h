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
  $Id: basic_iterator_position.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_ITERATOR_POSITION_H
#define SEQAN_HEADER_BASIC_ITERATOR_POSITION_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Tag

struct PositionIterator;

//////////////////////////////////////////////////////////////////////////////
// Position Iterator
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Position Iterator:
..cat:Iterators
..general:Class.Iter
..summary:Adapts @Metafunction.Position.position@ to @Concept.Rooted Iterator.iterator@.
..signature:Iter<TContainer, PositionIterator>
..param.TContainer:Type of the container.
...metafunction:Metafunction.Container
..remarks
...text:Position Iterators provide the concept @Concept.Rooted Iterator@.
..see:Metafunction.Position
*/

template <typename TContainer>
class Iter<TContainer, PositionIterator>
{
public:
	typedef typename Position<TContainer>::Type TPosition;

	typename _Pointer<TContainer>::Type data_container;
	TPosition data_position;
//____________________________________________________________________________

public:
/**
.Memfunc.PositionIterator#Iter:
..class:Spec.Position Iterator
..summary:Constructor
..signature:Iter()
..signature:Iter(iter)
..signature:Iter(container [, position])
..param.iter:Another position iterator object.
..param.container:The corresponding container object.
...metafunction:Metafunction.Container
..param.position:A position in $container$. (optional)
...metafunction:Metafunction.Position
...remarks.text:If this argument is omitted, the adaptor iterator is initialized to the @Function.beginPosition.begin position@ of $container$.
*/
	Iter()
	{
SEQAN_CHECKPOINT
	}
	Iter(typename _Parameter<TContainer>::Type container_, TPosition position_ = 0):
		data_container(_toPointer(container_)),
		data_position(position_)
	{
SEQAN_CHECKPOINT
	}
	Iter(Iter const & other_):
		data_container(other_.data_container),
		data_position(other_.data_position)
	{
SEQAN_CHECKPOINT
	}
	template <typename TContainer2, typename TSpec2>
	Iter(Iter<TContainer2, TSpec2> const & other_)
	{
SEQAN_CHECKPOINT
		assign(*this, other_);
	}
	~Iter()
	{
SEQAN_CHECKPOINT
	}
	Iter const & 
	operator = (Iter const & other_)
	{
SEQAN_CHECKPOINT
		data_container = other_.data_container;
		data_position = other_.data_position;
		return *this;
	}
//____________________________________________________________________________

	friend inline typename _Parameter<TContainer>::Type 
	container(Iter & me)
	{
SEQAN_CHECKPOINT
		return _toParameter<TContainer>(me.data_container);
	}
	friend inline typename _Parameter<TContainer>::Type 
	container(Iter const & me)
	{
SEQAN_CHECKPOINT
		return _toParameter<TContainer>(me.data_container);
	}
//____________________________________________________________________________

	friend inline void
	setContainer(Iter & me,	typename _Parameter<TContainer>::Type container_)
	{
SEQAN_CHECKPOINT
		typename Position<Iter>::Type pos = position(me);
		me.data_container = _toPointer(container_);
		setPosition(me, pos);
	}

//____________________________________________________________________________

	friend inline TPosition &
	position(Iter & me)
	{
SEQAN_CHECKPOINT
		return me.data_position;
	}
	friend inline TPosition const &
	position(Iter const & me)
	{
SEQAN_CHECKPOINT
		return me.data_position;
	}
//____________________________________________________________________________

	friend inline void
	setPosition(Iter & me, TPosition position_)
	{
SEQAN_CHECKPOINT
		me.data_position = position_;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// value
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline typename Reference<Iter<TContainer, PositionIterator> >::Type 
value(Iter<TContainer, PositionIterator> & me)
{
SEQAN_CHECKPOINT
	return value(container(me), position(me));
}
template <typename TContainer>
inline typename Reference<Iter<TContainer, PositionIterator> >::Type 
value(Iter<TContainer, PositionIterator> const & me)
{
SEQAN_CHECKPOINT
	return value(container(me), position(me));
}

/////////////////////////////////////////////////////////////////////////////
// assignValue
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TValue>
inline void
assignValue(Iter<TContainer, PositionIterator> & me,
			TValue _value)
{
SEQAN_CHECKPOINT
	assignValue(container(me), position(me), _value);
}
template <typename TContainer, typename TValue>
inline void
assignValue(Iter<TContainer, PositionIterator> const & me,
			TValue _value)
{
SEQAN_CHECKPOINT
	assignValue(container(me), position(me), _value);
}

/////////////////////////////////////////////////////////////////////////////
// moveValue
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TValue>
inline void
moveValue(Iter<TContainer, PositionIterator> & me,
		  TValue _value)
{
SEQAN_CHECKPOINT
	moveValue(container(me), position(me), _value);
}
template <typename TContainer, typename TValue>
inline void
moveValue(Iter<TContainer, PositionIterator> const & me,
		  TValue _value)
{
SEQAN_CHECKPOINT
	moveValue(container(me), position(me), _value);
}

//////////////////////////////////////////////////////////////////////////////
// operator ==
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator == (Iter<TContainer, PositionIterator> const & left,
			 Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) == position(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator !=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator != (Iter<TContainer, PositionIterator> const & left,
			 Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) != position(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator < / >
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator < (Iter<TContainer, PositionIterator> const & left,
			Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) < position(right);
}

template <typename TContainer>
inline bool 
operator > (Iter<TContainer, PositionIterator> const & left,
			Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) > position(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator <= / >=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator <= (Iter<TContainer, PositionIterator> const & left,
			 Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) <= position(right);
}

template <typename TContainer>
inline bool 
operator >= (Iter<TContainer, PositionIterator> const & left,
			 Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) >= position(right);
}

//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline void
goNext(Iter<TContainer, PositionIterator> & me)
{
SEQAN_CHECKPOINT
	setPosition(me, position(me) + 1);
}

//////////////////////////////////////////////////////////////////////////////
// goPrevious
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline void
goPrevious(Iter<TContainer, PositionIterator> & me)
{
SEQAN_CHECKPOINT
	setPosition(me, position(me) - 1);
}

//////////////////////////////////////////////////////////////////////////////
// operator +
//////////////////////////////////////////////////////////////////////////////
template <typename TContainer, typename TIntegral>
inline Iter<TContainer, PositionIterator>  
operator + (Iter<TContainer, PositionIterator> const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, PositionIterator>(container(left), position(left) + right);
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, PositionIterator>  
operator + (Iter<TContainer, PositionIterator> const & left,
			int right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, PositionIterator>(container(left), position(left) + right);
}

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, PositionIterator>  
operator + (TIntegral left,
			Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, PositionIterator>(container(right), position(right) + left);
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, PositionIterator>  
operator + (int left,
			Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, PositionIterator>(container(right), position(right) + left);
}

//////////////////////////////////////////////////////////////////////////////
// operator +=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, PositionIterator> &
operator += (Iter<TContainer, PositionIterator> & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	setPosition(left, position(left) + right);
	return left;
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, PositionIterator> &
operator += (Iter<TContainer, PositionIterator> & left,
			 int right)
{
SEQAN_CHECKPOINT
	setPosition(left, position(left) + right);
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// operator -
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, PositionIterator>  
operator - (Iter<TContainer, PositionIterator> const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, PositionIterator>(container(left), position(left) - right);
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, PositionIterator>  
operator - (Iter<TContainer, PositionIterator> const & left,
			int right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, PositionIterator>(container(left), position(left) - right);
}

//____________________________________________________________________________

template <typename TContainer>
inline typename Difference<TContainer>::Type  
operator - (Iter<TContainer, PositionIterator> const & left,
			Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) - position(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator -=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, PositionIterator> &
operator -= (Iter<TContainer, PositionIterator> & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	setPosition(left, position(left) - right);
	return left;
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, PositionIterator> &
operator -= (Iter<TContainer, PositionIterator> & left,
			 int right)
{
SEQAN_CHECKPOINT
	setPosition(left, position(left) - right);
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// assign (Conversion)
//////////////////////////////////////////////////////////////////////////////

template <typename TTargetContainer, typename TSource>
inline void
assign(Iter<TTargetContainer, PositionIterator> & target,
	   TSource const & source)
{
SEQAN_CHECKPOINT
	target.data_container = container(source);
	target.data_position = position(source);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
