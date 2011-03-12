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
  $Id: segment_suffix.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEGMENT_SUFFIX_H
#define SEQAN_HEADER_SEGMENT_SUFFIX_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// SuffixSegment
//////////////////////////////////////////////////////////////////////////////


/**
.Spec.SuffixSegment:
..cat:Segments
..summary:End part segment of a sequence.
..general:Class.Segment
..signature:Segment<THost, SuffixSegment>
..param.THost:Type of the whole sequence.
...text:Instances of $Segment<THost, SuffixSegment>$ are suffixes of $THost$ objects.
...remarks:Use @Metafunction.Host@ to get the host type for a given class.
..remarks.note:Since the appropriate segment type depends on the host sequence type, 
	it is recommended to use the metafunction @Metafunction.Suffix@ instead of explicitely 
	choose a specialization of @Class.Segment@.
..see:Spec.InfixSegment
..see:Metafunction.Suffix
*/

struct SuffixSegment;

template <typename THost_>
class Segment<THost_, SuffixSegment>
{
protected:
	typedef typename Host<Segment>::Type THost;

	typename _Pointer<THost>::Type data_host;
	typename Position<THost>::Type data_begin_position;

//____________________________________________________________________________

public:

/**
.Memfunc.SuffixSegment#Segment:
..class:Spec.SuffixSegment
..summary:Constructor
..signature:Segment<THost, SuffixSegment> ()
..signature:Segment<THost, SuffixSegment> (suffix)
..signature:Segment<THost, SuffixSegment> (host [, begin])
..param.suffix:Other suffix object. (copy constructor)
..param.host:The whole sequence.
..param.begin:Position in $host$ of the first item in segment. (optional)
...default:$0$
...type:Metafunction.Position.$Position<THost>::Type$
...type:Metafunction.Iterator.$Iterator<THost>::Type$
..remarks:
...text:A Segment object cannot work without a host. If the object is default constructed,
the host must be set by @Function.setHost@ before the segment can be used.
...text:If a segment object is constructed by the copy constructor, the
members of the new constructed object are set to the same values as the members in the
source object; the host object is not modified.
Note that this is a special case, since all other copy operations result in changes 
of the host object.
...text:$begin$ must be a valid position/iterator in $host$.
If $begin$ is omitted, the suffix segment corresponding to
the whole sequence $host$ is constructed.
This is the same segment that is returned by @Function.goBegin@.
*/
	Segment():
		data_begin_position(0)
	{
SEQAN_CHECKPOINT
	}

	Segment(THost & _host):
		data_host(& _host),
		data_begin_position(0)
	{
SEQAN_CHECKPOINT
	}

	Segment(typename _Parameter<THost>::Type _host, typename Position<THost>::Type _begin_index):
		data_host(_toPointer(_host)),
		data_begin_position(_begin_index)
	{
SEQAN_CHECKPOINT
	}
/*
	Segment(typename _Parameter<THost>::Type _host, typename Iterator<THost, Rooted>::Type _begin):
		data_host(_toPointer(_host)),
		data_begin_position(position(_begin))
	{
SEQAN_CHECKPOINT
	}
*/
	Segment(typename _Parameter<THost>::Type _host, typename Iterator<THost, Standard>::Type _begin):
		data_host(_toPointer(_host)),
		data_begin_position(position(_begin, _host))
	{
SEQAN_CHECKPOINT
	}
/*
	Segment(Segment const & _other):
		data_host(_other.data_host),
		data_begin_position(_other.data_begin_position)
	{
SEQAN_CHECKPOINT
	}
*/
	template <typename THost2, typename TSpec2>
	Segment(Segment<THost2, TSpec2> const & _other):
		data_host(_toPointer(host(_other))),
		data_begin_position(beginPosition(_other))
	{
SEQAN_CHECKPOINT
	}

	~ Segment() 
	{
SEQAN_CHECKPOINT
	}

	template <typename TSource>
	inline Segment & 
	operator = (TSource const & source)
	{
		assign(*this, source);
		return *this;
	}
	inline Segment & 
	operator = (Segment const & source)
	{
		assign(*this, source);
		return *this;
	}
//____________________________________________________________________________

public:

	friend inline typename _Parameter<THost>::Type 
	host(Segment & me)
	{
SEQAN_CHECKPOINT
		return _toParameter<THost>(me.data_host);
	}

	friend inline typename _Parameter<THost>::Type 
	host(Segment const & me)
	{
SEQAN_CHECKPOINT
		return _toParameter<THost>(me.data_host);
	}

//____________________________________________________________________________

	friend inline void 
	setHost(Segment & me, typename _Parameter<THost>::Type _host)
	{
SEQAN_CHECKPOINT
		me.data_host = _toPointer(_host);
	}

//____________________________________________________________________________

	template <typename TPos>
	inline typename Reference<Segment>::Type
	operator [] (TPos pos)
	{
SEQAN_CHECKPOINT
		return value(*this, pos);
	}

	template <typename TPos>
	inline typename Reference<Segment const>::Type 
	operator [] (TPos pos) const
	{
SEQAN_CHECKPOINT
		return value(*this, pos);
	}

//____________________________________________________________________________

	friend inline typename Iterator<Segment, Standard>::Type 
	begin(Segment & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return begin(host(me), Standard()) + me.data_begin_position;
	}
	friend inline typename Iterator<Segment const, Standard>::Type 
	begin(Segment const & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return begin(host(me), Standard()) + me.data_begin_position;
	}

//____________________________________________________________________________

	friend inline typename Position<Segment const>::Type 
	beginPosition(Segment const & me)
	{
SEQAN_CHECKPOINT
		return me.data_begin_position;
	}
	friend inline typename Position<Segment>::Type 
	beginPosition(Segment & me)
	{
SEQAN_CHECKPOINT
		return me.data_begin_position;
	}
//____________________________________________________________________________

	template <typename TIterator>
	friend inline void 
	setBegin(Segment & me, TIterator new_begin)
	{
SEQAN_CHECKPOINT
		me.data_begin_position = new_begin - begin(host(me));//, Standard());
	}

	friend inline void 
	setBegin(typename Iterator<Segment, Rooted>::Type new_begin)
	{
SEQAN_CHECKPOINT
		container(new_begin).data_begin_position = hostIterator(new_begin) - begin(host(container(new_begin)));//, Standard());
	}

//____________________________________________________________________________

	template <typename TPosition>
	friend inline void 
	setBeginPosition(Segment & me, TPosition new_begin)
	{
SEQAN_CHECKPOINT
		me.data_begin_position = new_begin;
	}

//____________________________________________________________________________

	friend inline typename Iterator<Segment, Standard>::Type 
	end(Segment & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return end(host(me), Standard());
	}
	friend inline typename Iterator<Segment const, Standard>::Type 
	end(Segment const & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return end(host(me), Standard());
	}

//____________________________________________________________________________


	friend inline typename Position<Segment>::Type 
	endPosition(Segment & me)
	{
SEQAN_CHECKPOINT
		return length(host(me));
	}

	friend inline typename Position<Segment const>::Type 
	endPosition(Segment const & me)
	{
SEQAN_CHECKPOINT
		return length(host(me));
	}

//____________________________________________________________________________

	template <typename TIterator>
	friend inline void
	setEnd(Segment &, TIterator)
	{
	}

	friend inline void 
	_setLength(Segment &, typename Size<THost>::Type)
	{
	}

//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Suffix:
..summary:Suffix sequence type.
..signature:Suffix<T>::Type
..param.T:A sequence type.
...type:Class.String
..returns.param.Type:The suffix type.
..see:Spec.SuffixSegment
..see:Metafunction.Infix
..see:Metafunction.Prefix
*/

struct PrefixSegment;
struct InfixSegment;

template <typename THost>
struct Suffix
{
	typedef Segment<THost, SuffixSegment> Type;
};

template <typename THost>
struct Suffix< Segment<THost, InfixSegment> >
{
	typedef Segment<THost, InfixSegment> Type;
};
template <typename THost>
struct Suffix< Segment<THost, SuffixSegment> >
{
	typedef Segment<THost, SuffixSegment> Type;
};
template <typename THost>
struct Suffix< Segment<THost, PrefixSegment> >
{
	typedef Segment<THost, InfixSegment> Type;
};

template <typename THost, typename TSpec>
struct Suffix< Segment<THost, TSpec> const >:
	Suffix< Segment<THost, TSpec> > {};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TPosition>
inline void
set(Segment<THost, SuffixSegment> & me,
	THost & host_,
	TPosition begin_)
{
SEQAN_CHECKPOINT
	setHost(me, host_);
	setBegin(me, begin_);
}
//____________________________________________________________________________

template <typename THost>
inline void
set(Segment<THost, SuffixSegment> & me,
	THost & host_)
{
SEQAN_CHECKPOINT
	setHost(me, host_);
	setBegin(me, begin(host_, Standard()));
}

//____________________________________________________________________________

template <typename THost, typename TSpec>
inline void
set(Segment<THost, SuffixSegment> & me,
	Segment<THost, TSpec> & source)
{
SEQAN_CHECKPOINT
	setHost(me, host(source));
	setBeginPosition(me, beginPosition(source));
}

template <typename THost, typename TSpec>
inline void
set(Segment<THost, SuffixSegment> & me,
	Segment<THost, TSpec> const & source)
{
SEQAN_CHECKPOINT
	setHost(me, host(source));
	setBeginPosition(me, beginPosition(source));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atBegin(Segment<THost, SuffixSegment> const & segment)
{
SEQAN_CHECKPOINT
	return (beginPosition(segment) == 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atEnd(Segment<THost, SuffixSegment> const & segment)
{
SEQAN_CHECKPOINT
	return (beginPosition(segment) == length(host(segment)));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline void
goBegin(Segment<THost, SuffixSegment> & segment,
		THost &)
{
SEQAN_CHECKPOINT
	goBegin(segment);
}

template <typename THost>
inline void
goBegin(Segment<THost, SuffixSegment> & segment)
{
	setBegin(segment);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline void
goEnd(Segment<THost, SuffixSegment> & segment,
	  THost &)
{
SEQAN_CHECKPOINT
	goEnd(segment);
}

template <typename THost>
inline void
goEnd(Segment<THost, SuffixSegment> & segment)
{
	setBegin(segment, length(host(segment))-1);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, SuffixSegment> &
operator ++(Segment<THost, SuffixSegment> & segment)
{
	setBegin(segment, beginPosition(segment) + 1);
	return segment;
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, SuffixSegment> &
operator --(Segment<THost, SuffixSegment> & segment)
{
	setBegin(segment, beginPosition(segment) - 1);
	return segment;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.suffix:
..cat:Containers
..summary:Creates suffix object.
..signature:suffix(host, begin)
..param.host:The complete sequence.
...type:Class.String
...type:Adaption.char array
..param.begin:Position or iterator of the first element of the segment.
...type:Metafunction.Position
...type:Metafunction.Iterator
..returns:The suffix of $host that begins at $begin$.
...remarks:The type of the suffix is given by @Metafunction.Suffix@.
..remarks:Notational sugar.
..see:Spec.SuffixSegment
..see:Function.infix
*/

template <typename T, typename TPosBegin>
inline typename Suffix<T>::Type
suffix(T & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
	return typename Suffix<T>::Type(t, pos_begin);
}

template <typename T, typename TPosBegin>
inline typename Suffix<T *>::Type
suffix(T * t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
	return typename Suffix<T *>::Type (t, pos_begin);
}

//////////////////////////////////////////////////////////////////////////////
// A suffix of a prefix -> is an infix
template <typename T, typename TPosBegin>
inline typename Suffix<Segment<T, PrefixSegment> >::Type
suffix(Segment<T, PrefixSegment> & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
	return typename Suffix<Segment<T, PrefixSegment> >::Type (
		host(t), 
		beginPosition(t) + pos_begin, 
		endPosition(t));
}
template <typename T, typename TPosBegin>
inline typename Suffix<Segment<T, PrefixSegment> const>::Type
suffix(Segment<T, PrefixSegment> const & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
	return typename Suffix<Segment<T, PrefixSegment> const>::Type (
		host(t), 
		beginPosition(t) + pos_begin, 
		endPosition(t));
}

//////////////////////////////////////////////////////////////////////////////
// A suffix of a infix -> is an infix
template <typename T, typename TPosBegin>
inline typename Suffix<Segment<T, InfixSegment> >::Type
suffix(Segment<T, InfixSegment> & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
	return typename Suffix<Segment<T, InfixSegment> >::Type (
		host(t), 
		beginPosition(t) + pos_begin, 
		endPosition(t));
}
template <typename T, typename TPosBegin>
inline typename Suffix<Segment<T, InfixSegment> const>::Type
suffix(Segment<T, InfixSegment> const & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
	return typename Suffix<Segment<T, InfixSegment> const>::Type (
		host(t), 
		beginPosition(t) + pos_begin, 
		endPosition(t));
}


//////////////////////////////////////////////////////////////////////////////
// A suffix of a suffix -> is a suffix
template <typename T, typename TPosBegin>
inline typename Suffix<Segment<T, SuffixSegment> >::Type
suffix(Segment<T, SuffixSegment> & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
	return typename Suffix<Segment<T, SuffixSegment> >::Type (
		host(t), 
		beginPosition(t) + pos_begin);
}
template <typename T, typename TPosBegin>
inline typename Suffix<Segment<T, SuffixSegment> const>::Type
suffix(Segment<T, SuffixSegment> const & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
	return typename Suffix<Segment<T, SuffixSegment> const>::Type (
		host(t), 
		beginPosition(t) + pos_begin);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
