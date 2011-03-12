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
  $Id: std_string.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_STD_STRING_H
#define SEQAN_HEADER_STD_STRING_H


//Adaption for ::std::basic_string

#include <string>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
/**
.Adaption."std::basic_string":
..summary:Standard library string class.
*/



//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IsContiguous.param.T.type:Adaption.std::basic_string

template <typename  TChar, typename TCharTraits, typename TAlloc>
struct IsContiguous< ::std::basic_string<TChar, TCharTraits, TAlloc> >
{
    enum { VALUE = true };
};

template <typename  TChar, typename TCharTraits, typename TAlloc>
struct IsContiguous< ::std::basic_string<TChar, TCharTraits, TAlloc> const>
{
    enum { VALUE = true };
};

///.Metafunction.Value.param.T.type:Adaption.std::basic_string
template <typename TChar, typename TCharTraits, typename TAlloc>
struct Value< ::std::basic_string<TChar, TCharTraits, TAlloc> >
{
	typedef typename ::std::basic_string<TChar, TCharTraits, TAlloc>::value_type Type;
};
template <typename TChar, typename TCharTraits, typename TAlloc>
struct Value< ::std::basic_string<TChar, TCharTraits, TAlloc> const>
{
	typedef typename ::std::basic_string<TChar, TCharTraits, TAlloc>::value_type Type;
};

///.Metafunction.GetValue.param.T.type:Adaption.std::basic_string
template <typename TChar, typename TCharTraits, typename TAlloc>
struct GetValue< ::std::basic_string<TChar, TCharTraits, TAlloc> >
{
	typedef typename ::std::basic_string<TChar, TCharTraits, TAlloc>::reference Type;
};
template <typename TChar, typename TCharTraits, typename TAlloc>
struct GetValue< ::std::basic_string<TChar, TCharTraits, TAlloc> const>
{
	typedef typename ::std::basic_string<TChar, TCharTraits, TAlloc>::const_reference Type;
};

//???GetValue<vector<bool> > ist bool

//____________________________________________________________________________

///.Metafunction.Iterator.param.T.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc>, Rooted>
{
	typedef ::std::basic_string<TChar, TCharTraits, TAlloc> TString;
	typedef Iter<TString, StdIteratorAdaptor> TIterator;
	typedef Iter<TString, AdaptorIterator<TIterator> > Type;
};
template <typename TChar, typename TCharTraits, typename TAlloc>
struct Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc> const, Rooted>
{
	typedef ::std::basic_string<TChar, TCharTraits, TAlloc> const TString;
	typedef Iter<TString, StdIteratorAdaptor> TIterator;
	typedef Iter<TString, AdaptorIterator<TIterator> > Type;
};


template <typename TChar, typename TCharTraits, typename TAlloc>
struct Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc>, Standard >
{
	typedef Iter< ::std::basic_string<TChar, TCharTraits, TAlloc>, StdIteratorAdaptor > Type;
};
template <typename TChar, typename TCharTraits, typename TAlloc>
struct Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>
{
	typedef Iter< ::std::basic_string<TChar, TCharTraits, TAlloc> const, StdIteratorAdaptor > Type;
};

//____________________________________________________________________________

///.Metafunction.Position.param.T.type:Adaption.std::basic_string
template <typename TChar, typename TCharTraits, typename TAlloc>
struct Position< ::std::basic_string<TChar, TCharTraits, TAlloc> >
{
	typedef typename ::std::basic_string<TChar, TCharTraits, TAlloc>::size_type Type;
};

//____________________________________________________________________________

///.Metafunction.Size.param.T.type:Adaption.std::basic_string
template <typename TChar, typename TCharTraits, typename TAlloc>
struct Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >
{
	typedef typename ::std::basic_string<TChar, TCharTraits, TAlloc>::size_type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Size.param.T.type:Adaption.std::basic_string
template <typename TChar, typename TCharTraits, typename TAlloc>
struct DefaultOverflowImplicit< ::std::basic_string<TChar, TCharTraits, TAlloc> >
{
	typedef Generous Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Function.id.param.object.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
inline void const * 
id(::std::basic_string<TChar, TCharTraits, TAlloc> const & me)
{
SEQAN_CHECKPOINT
	return & *end(me, Standard());
}

//////////////////////////////////////////////////////////////////////////////

///.Function.begin.param.object.type:Adaption.std::basic_string
template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc>, Standard>::Type 
begin(::std::basic_string<TChar, TCharTraits, TAlloc> & me,
	  Standard)
{
SEQAN_CHECKPOINT
	return typename Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc>, Standard>::Type(me.begin());
}
template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>::Type 
begin(::std::basic_string<TChar, TCharTraits, TAlloc> const & me,
	  Standard)
{
SEQAN_CHECKPOINT
	return typename Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>::Type(me.begin());
}

//////////////////////////////////////////////////////////////////////////////

///.Function.end.param.object.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc>, Standard>::Type 
end(::std::basic_string<TChar, TCharTraits, TAlloc> & me,
	Standard)
{
SEQAN_CHECKPOINT
	return typename Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc>, Standard>::Type(me.end());
}
template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>::Type 
end(::std::basic_string<TChar, TCharTraits, TAlloc> const & me,
	Standard)
{
SEQAN_CHECKPOINT
	return typename Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>::Type(me.end());
}

//////////////////////////////////////////////////////////////////////////////

///.Function.value.param.container.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TPos>
inline typename GetValue< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type
value(::std::basic_string<TChar, TCharTraits, TAlloc> & me, 
	  TPos pos)
{
SEQAN_CHECKPOINT
	return me[pos];
} 
template <typename TChar, typename TCharTraits, typename TAlloc, typename TPos>
inline typename GetValue< ::std::basic_string<TChar, TCharTraits, TAlloc> const>::Type
value(::std::basic_string<TChar, TCharTraits, TAlloc> const & me, 
	  TPos pos)
{
SEQAN_CHECKPOINT
	return me[pos];
}

//////////////////////////////////////////////////////////////////////////////

///.Function.length.param.object.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type
length(::std::basic_string<TChar, TCharTraits, TAlloc> const & me)
{
SEQAN_CHECKPOINT
	return me.length();
}

//////////////////////////////////////////////////////////////////////////////

///.Function.capacity.param.object.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type
capacity(::std::basic_string<TChar, TCharTraits, TAlloc> const & me)
{
SEQAN_CHECKPOINT
	return me.capacity();
}

//////////////////////////////////////////////////////////////////////////////

///.Function.empty.param.object.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
inline bool
empty(::std::basic_string<TChar, TCharTraits, TAlloc> const & me)
{
SEQAN_CHECKPOINT
	return me.empty();
}

//////////////////////////////////////////////////////////////////////////////

///.Function.clear.param.object.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
inline void
clear(::std::basic_string<TChar, TCharTraits, TAlloc> & me)
{
SEQAN_CHECKPOINT
	me.clear();
}

//////////////////////////////////////////////////////////////////////////////
//assign to ::std::basic_string

///.Function.assign.param.target.type:Adaption.std::basic_string
///.Function.assign.param.source.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
assign(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource & source)
{
SEQAN_CHECKPOINT
	assign(target, source, Generous());
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
assign(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource const & source)
{
SEQAN_CHECKPOINT
	assign(target, source, Generous());
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TSize>
inline void 
assign(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource & source,
	   TSize limit)
{
SEQAN_CHECKPOINT
	assign(target, source, limit, Generous());
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TSize>
inline void 
assign(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource const & source,
	   TSize limit)
{
SEQAN_CHECKPOINT
	assign(target, source, limit, Generous());
}

//____________________________________________________________________________

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
assign(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource & source,
	   Generous)
{
SEQAN_CHECKPOINT
	target.assign(begin(source, Standard()), end(source, Standard()));
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
assign(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource const & source,
	   Generous)
{
SEQAN_CHECKPOINT
	target.assign(begin(source, Standard()), end(source, Standard()));
}


template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
assign_std_string_Generous_impl(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
								TSource & source,
								typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit)
{
SEQAN_CHECKPOINT
	typename Iterator<TSource const, Standard>::Type source_begin = begin(source, Standard());
	typename Size<TSource const>::Type source_length = length(source);
	if (source_length > limit)
	{
		source_length = limit;
	}
	target.assign(source_begin, source_begin + source_length);
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
assign(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource & source,
	   typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
	   Generous)
{
SEQAN_CHECKPOINT
	assign_std_string_Generous_impl(target, source, limit);
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
assign(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource const & source,
	   typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
	   Generous)
{
SEQAN_CHECKPOINT
	assign_std_string_Generous_impl(target, source, limit);
}

//____________________________________________________________________________

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
assign(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource & source,
	   Limit)
{
SEQAN_CHECKPOINT
	assign(target, source, target.capacity(), Generous());
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
assign(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource const & source,
	   Limit)
{
SEQAN_CHECKPOINT
	assign(target, source, target.capacity(), Generous());
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
assign(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource & source,
	   typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
	   Limit)
{
SEQAN_CHECKPOINT
	if (limit > target.capacity()) 
	{
		limit = target.capacity();
	}

	assign(target, source, limit, Generous());
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
assign(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource const & source,
	   typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
	   Limit)
{
SEQAN_CHECKPOINT
	if (limit > target.capacity()) 
	{
		limit = target.capacity();
	}

	assign(target, source, limit, Generous());
}

//////////////////////////////////////////////////////////////////////////////
//append to ::std::basic_string

///.Function.append.param.target.type:Adaption.std::basic_string
///.Function.append.param.source.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
append(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource const & source,
	   Generous)
{
SEQAN_CHECKPOINT
	target.append(begin(source, Standard()), end(source, Standard()));
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
append(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource const & source,
	   typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
	   Generous)
{
SEQAN_CHECKPOINT
	typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type target_length = target.length();
	if (target_length > limit)
	{
		target.resize(limit);
	}
	else
	{
		limit -= target_length;
		typename Iterator<TSource const, Standard>::Type source_begin = begin(source, Standard());
		typename Size<TSource const>::Type source_length = length(source);
		if (source_length > limit)
		{
			source_length = limit;
		}
		target.append(source_begin, source_begin + source_length);
	}
}

//____________________________________________________________________________

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
append(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource const & source,
	   Limit)
{
SEQAN_CHECKPOINT
	append(target, source, target.capacity(), Generous());
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
append(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
	   TSource const & source,
	   typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
	   Limit)
{
SEQAN_CHECKPOINT
	if (limit > target.capacity()) 
	{
		limit = target.capacity();
	}

	append(target, source, limit, Generous());
}

//////////////////////////////////////////////////////////////////////////////
///.Function.appendValue.param.target.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TValue, typename TTag>
inline void
appendValue(::std::basic_string<TChar, TCharTraits, TAlloc> & me, 
			TValue const & _value,
			TTag)
{
SEQAN_CHECKPOINT
	me.push_back(_value);
} 

template <typename TChar, typename TCharTraits, typename TAlloc, typename TValue>
inline void
appendValue(::std::basic_string<TChar, TCharTraits, TAlloc> & me, 
			TValue const & _value,
			Limit)
{
SEQAN_CHECKPOINT
	if (capacity(me) > length(me)) me.push_back(_value);
} 

//////////////////////////////////////////////////////////////////////////////
//replace to ::std::basic_string

///.Function.replace.param.target.type:Adaption.std::basic_string
///.Function.replace.param.source.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
replace(::std::basic_string<TChar, TCharTraits, TAlloc> & target,
		typename Position< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_begin,
		typename Position< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_end,
		TSource const & source,
		Generous)
{
SEQAN_CHECKPOINT
	target.replace(target.begin() + pos_begin, target.begin() + pos_end, begin(source, Standard()), end(source, Standard()));
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
replace(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
		typename Position< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_begin,
		typename Position< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_end,
		TSource const & source,
		typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
		Generous)
{
SEQAN_CHECKPOINT
	if (pos_begin >= limit)
	{
		target.resize(limit);
	}
	else
	{
		typename Iterator<TSource const, Standard>::Type source_begin = begin(source, Standard());
		typename Size<TSource const>::Type source_length = length(source);
		typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_mid = pos_begin + source_length;
		if (pos_mid > limit)
		{
			target.replace(target.begin() + pos_begin, target.begin() + limit, source_begin, source_begin + limit - pos_begin);
			target.resize(limit);
		}
		else
		{
			target.replace(target.begin() + pos_begin, target.begin() + pos_end, source_begin, end(source, Standard()));
			if (target.length() > limit)
			{
				target.resize(limit);
			}
		}
	}
}

//____________________________________________________________________________

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
replace(::std::basic_string<TChar, TCharTraits, TAlloc> & target,
		typename Position< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_begin,
		typename Position< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_end,
		TSource const & source,
		Limit)
{
SEQAN_CHECKPOINT
	replace(target, pos_begin, pos_end, source, target.capacity(), Generous());
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void 
replace(::std::basic_string<TChar, TCharTraits, TAlloc> & target, 
		typename Position< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_begin,
		typename Position< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_end,
		TSource const & source,
		typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
		Limit)
{
SEQAN_CHECKPOINT
	if (limit > target.capacity()) 
	{
		limit = target.capacity();
	}

	replace(target, pos_begin, pos_end, source, limit, Generous());
}

//////////////////////////////////////////////////////////////////////////////
// handling of iterators as begin and end

/*
template<typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TExpand>
inline void 
replace(::std::basic_string<TChar, TCharTraits, TAlloc> & target,
		typename Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc>, Rooted>::Type pos_begin,
		typename Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc>, Rooted>::Type pos_end,
		TSource & source,
		Tag<TExpand> const tag)
{
	replace(target, position(pos_begin), position(pos_end), source, tag);
}

template<typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TExpand>
inline void 
replace(::std::basic_string<TChar, TCharTraits, TAlloc> & target,
		typename Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc>, Rooted>::Type pos_begin,
		typename Iterator< ::std::basic_string<TChar, TCharTraits, TAlloc>, Rooted>::Type pos_end,
		TSource & source,
		typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
		Tag<TExpand> const tag)
{
	replace(target,  position(pos_begin),  position(pos_end), source, tag);
}
*/

//////////////////////////////////////////////////////////////////////////////

///.Function.reserve.param.object.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TExpand>
inline typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type 
reserve(
	::std::basic_string<TChar, TCharTraits, TAlloc> & seq, 
	typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type new_capacity,
	Tag<TExpand> const &)
{
SEQAN_CHECKPOINT
    seq.reserve(new_capacity);
    if (new_capacity < seq.capacity())
    {
        return seq.capacity();
    }
	return new_capacity;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.resize.param.object.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSize, typename TExpand>
inline typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type 
resize(
	::std::basic_string<TChar, TCharTraits, TAlloc> & me,
	TSize new_length,
	Tag<TExpand> const &)
{
SEQAN_CHECKPOINT
    me.resize(new_length);
	return me.length();
}

//////////////////////////////////////////////////////////////////////////////

///.Function.fill.param.object.type:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSize, typename TExpand>
inline typename Size< ::std::basic_string<TChar, TCharTraits, TAlloc> >::Type 
fill(
	::std::basic_string<TChar, TCharTraits, TAlloc> & me,
	TSize new_length,
	TChar const & val,
	Tag<TExpand> const &)
{
SEQAN_CHECKPOINT
    me.resize(new_length, val);
	return me.length();
}


/* (veraltet)
//////////////////////////////////////////////////////////////////////////////
// Iterator Handling
//////////////////////////////////////////////////////////////////////////////

//??????
template <typename TChar, typename TContainer>
struct GetValue< ::__gnu_cxx::__normal_iterator<TChar, TContainer> >
{
	typedef typename ::__gnu_cxx::__normal_iterator<TChar, TContainer>::reference Type;
};
*/

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
