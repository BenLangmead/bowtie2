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
  $Id: basic_iterator_adapt_std.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_ITERATOR_ADAPT_STD_H
#define SEQAN_HEADER_BASIC_ITERATOR_ADAPT_STD_H

//////////////////////////////////////////////////////////////////////////////

//adapt SeqAn iterator to std
namespace std
{
	template<typename TContainer, typename TSpec>
	struct iterator_traits<seqan::Iter<TContainer, TSpec> >
	{
		typedef ::seqan::Iter<TContainer, TSpec> TIter;

		typedef random_access_iterator_tag iterator_category;
		typedef typename ::seqan::Value<TIter>::Type value_type;
		typedef typename ::seqan::Difference<TIter>::Type difference_type;
		typedef typename ::seqan::Value<TIter>::Type * pointer;
		typedef typename ::seqan::Reference<TIter>::Type reference;
	};
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

//helper Metafunction

/* This simple, general implementation cannot be used due to strange VC++ 2003 behavior

template <typename TStdContainer>
struct StdContainerIterator
{
	typedef typename TStdContainer::iterator Type;
};

template <typename TStdContainer>
struct StdContainerIterator<TStdContainer const>
{
	typedef typename TStdContainer::const_iterator Type;
};
*/

//we use this instead: specialize StdContainerIterator for each std-container
template <typename TStdContainer>
struct StdContainerIterator
{
	typedef void * Type; //dummy, just to make VC++ 2003 happy
};

template <typename TChar, typename TCharTraits, typename TAlloc>
struct StdContainerIterator< ::std::basic_string<TChar, TCharTraits, TAlloc> >
{
	typedef ::std::basic_string<TChar, TCharTraits, TAlloc> TContainer;
	typedef typename TContainer::iterator Type;
};
template <typename TChar, typename TCharTraits, typename TAlloc>
struct StdContainerIterator< ::std::basic_string<TChar, TCharTraits, TAlloc> const>
{
	typedef ::std::basic_string<TChar, TCharTraits, TAlloc> TContainer;
	typedef typename TContainer::const_iterator Type;
};

//////////////////////////////////////////////////////////////////////////////
//adapt std iterator to SeqAn


struct StdIteratorAdaptor;

template <typename TContainer>
class Iter<TContainer, StdIteratorAdaptor>
{
public:
	typedef typename StdContainerIterator<TContainer>::Type TIterator;
	TIterator data_iterator;

	Iter() {}
	Iter(Iter const & other_): data_iterator(other_.data_iterator) {}
	Iter(TIterator const & iter_): data_iterator(iter_) {}
	Iter(TContainer const & cont_): data_iterator(begin(cont_)) {}

	Iter const & operator = (Iter const & other_)
	{
		data_iterator = other_.data_iterator;
		return *this;
	}
	Iter const & operator = (TIterator const & iter_)
	{
		data_iterator = iter_;
		return *this;
	}

	operator TIterator &()
	{
		return data_iterator;
	}

	~Iter() {}
};

//////////////////////////////////////////////////////////////////////////////
// hostIterator
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline typename StdContainerIterator<TContainer>::Type &
hostIterator(Iter<TContainer, StdIteratorAdaptor> & me)
{
	return me.data_iterator;
}
template <typename TContainer>
inline typename StdContainerIterator<TContainer>::Type const &
hostIterator(Iter<TContainer, StdIteratorAdaptor> const & me)
{
	return me.data_iterator;
}

//////////////////////////////////////////////////////////////////////////////
// value
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline typename Reference<Iter<TContainer, StdIteratorAdaptor> >::Type 
value(Iter<TContainer, StdIteratorAdaptor> & me)
{
	return *(me.data_iterator);
}
template <typename TContainer>
inline typename Reference<Iter<TContainer, StdIteratorAdaptor> const>::Type 
value(Iter<TContainer, StdIteratorAdaptor> const & me)
{
	return *(me.data_iterator);
}


/////////////////////////////////////////////////////////////////////////////
// assignValue
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TValue>
inline void 
assignValue(Iter<TContainer, StdIteratorAdaptor> & me,
			TValue & val)
{
	*(me.data_iterator) = val;
}
template <typename TContainer, typename TValue>
inline void 
assignValue(Iter<TContainer, StdIteratorAdaptor> & me,
			TValue const & val)
{
	*(me.data_iterator) = val;
}

/////////////////////////////////////////////////////////////////////////////
// moveValue
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TValue>
inline void 
moveValue(Iter<TContainer, StdIteratorAdaptor> & me,
		  TValue & val)
{
	move(*(me.data_iterator), val);
}
template <typename TContainer, typename TValue>
inline void 
moveValue(Iter<TContainer, StdIteratorAdaptor> & me,
		  TValue const & val)
{
	move(*(me.data_iterator), val);
}

//////////////////////////////////////////////////////////////////////////////
// operator ==
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator == (Iter<TContainer, StdIteratorAdaptor> const & left,
			 Iter<TContainer, StdIteratorAdaptor> const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) == hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator !=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator != (Iter<TContainer, StdIteratorAdaptor> const & left,
			 Iter<TContainer, StdIteratorAdaptor> const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) != hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator <
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator < (Iter<TContainer, StdIteratorAdaptor> const & left,
			Iter<TContainer, StdIteratorAdaptor> const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) < hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator >
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator > (Iter<TContainer, StdIteratorAdaptor> const & left,
			Iter<TContainer, StdIteratorAdaptor> const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) > hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator <=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator <= (Iter<TContainer, StdIteratorAdaptor> const & left,
			Iter<TContainer, StdIteratorAdaptor> const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) <= hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator >=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator >= (Iter<TContainer, StdIteratorAdaptor> const & left,
			Iter<TContainer, StdIteratorAdaptor> const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) >= hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline void
goNext(Iter<TContainer, StdIteratorAdaptor> & me)
{
SEQAN_CHECKPOINT
	goNext(hostIterator(me));
}

//////////////////////////////////////////////////////////////////////////////
// goPrevious
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline void
goPrevious(Iter<TContainer, StdIteratorAdaptor> & me)
{
SEQAN_CHECKPOINT
	goPrevious(hostIterator(me));
}

//////////////////////////////////////////////////////////////////////////////
// operator +
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, StdIteratorAdaptor>  
operator + (Iter<TContainer, StdIteratorAdaptor> const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, StdIteratorAdaptor>(hostIterator(left) + right);
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, StdIteratorAdaptor>  
operator + (Iter<TContainer, StdIteratorAdaptor> const & left,
			int right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, StdIteratorAdaptor>(hostIterator(left) + right);
}

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, StdIteratorAdaptor>  
operator + (TIntegral left,
			Iter<TContainer, StdIteratorAdaptor> const & right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, StdIteratorAdaptor>(hostIterator(right) + left);
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, StdIteratorAdaptor>  
operator + (int left,
			Iter<TContainer, StdIteratorAdaptor> const & right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, StdIteratorAdaptor>(hostIterator(right) + left);
}

//////////////////////////////////////////////////////////////////////////////
// operator +=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, StdIteratorAdaptor> &
operator += (Iter<TContainer, StdIteratorAdaptor> & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	hostIterator(left) += right;
	return left;
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, StdIteratorAdaptor> &
operator += (Iter<TContainer, StdIteratorAdaptor> & left,
			 int right)
{
SEQAN_CHECKPOINT
	hostIterator(left) += right;
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// operator -
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, StdIteratorAdaptor>  
operator - (Iter<TContainer, StdIteratorAdaptor> const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, StdIteratorAdaptor>(hostIterator(left) - right);
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, StdIteratorAdaptor>  
operator - (Iter<TContainer, StdIteratorAdaptor> const & left,
			int right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, StdIteratorAdaptor>(hostIterator(left) - right);
}

//____________________________________________________________________________

template <typename TContainer>
inline typename Difference<Iter<TContainer, StdIteratorAdaptor> >::Type  
operator - (Iter<TContainer, StdIteratorAdaptor> const & left,
			Iter<TContainer, StdIteratorAdaptor> const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) - hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator -=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, StdIteratorAdaptor> &
operator -= (Iter<TContainer, StdIteratorAdaptor> & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	hostIterator(left) -= right;
	return left;
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, StdIteratorAdaptor> &
operator -= (Iter<TContainer, StdIteratorAdaptor> & left,
			 int right)
{
SEQAN_CHECKPOINT
	hostIterator(left) -= right;
	return left;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// assign (Conversion)
//////////////////////////////////////////////////////////////////////////////

template <typename TTargetContainer, typename TSource>
inline void
assign(Iter<TTargetContainer, StdIteratorAdaptor> & target,
	   TSource const & source)
{
SEQAN_CHECKPOINT
	target.data_iterator = begin(container(source)) + position(source);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
