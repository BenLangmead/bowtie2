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
  $Id: segment_prefix.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEGMENT_PREFIX_H
#define SEQAN_HEADER_SEGMENT_PREFIX_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// PrefixSegment
//////////////////////////////////////////////////////////////////////////////


/**
.Spec.PrefixSegment:
..cat:Segments
..summary:First part of a sequence.
..general:Class.Segment
..signature:Segment<THost, PrefixSegment>
..param.THost:Type of the whole sequence.
...text:Instances of $Segment<THost, PrefixSegment>$ are prefixes of $THost$ objects.
...remarks:Use @Metafunction.Host@ to get the host type for a given class.
..remarks.note:Since the appropriate segment type depends on the host sequence type, 
	it is recommended to use the metafunction @Metafunction.Prefix@ instead of explicitely 
	choose a specialization of @Class.Segment@.
..see:Spec.InfixSegment
..see:Spec.SuffixSegment
..see:Metafunction.Prefix
*/

struct PrefixSegment;

template <typename THost_>
class Segment<THost_, PrefixSegment>
{
protected:
	typedef typename Host<Segment>::Type THost;

	typename _Pointer<THost>::Type data_host;
	typename Position<THost>::Type data_end_position;

//____________________________________________________________________________

public:

/**
.Memfunc.PrefixSegment#Segment:
..class:Spec.PrefixSegment
..summary:Constructor
..signature:Segment<THost, PrefixSegment> ()
..signature:Segment<THost, PrefixSegment> (prefix)
..signature:Segment<THost, PrefixSegment> (host [, end])
..param.prefix:Other prefix object. (copy constructor)
..param.host:The whole sequence.
..param.end:Position in $host$ behind the last item in segment. (optional)
...default:$length(host)$
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
If $begin$ is omitted, the prefix segment corresponding to
the whole sequence $host$ is constructed.
This is the same segment that is returned by @Function.goBegin@.
*/
	Segment():
		data_end_position(0)
	{
SEQAN_CHECKPOINT
	}

	Segment(THost & _host):
		data_host(& _host),
		data_end_position(length(_host))
	{
SEQAN_CHECKPOINT
	}

	Segment(typename _Parameter<THost>::Type _host, typename Position<THost>::Type _end_index):
		data_host(_toPointer(_host)),
		data_end_position(_end_index)
	{
SEQAN_CHECKPOINT
	}
/*
	Segment(typename _Parameter<THost>::Type _host, typename Iterator<THost, Rooted>::Type _end):
		data_host(_toPointer(_host)),
		data_end_position(position(_end))
	{
SEQAN_CHECKPOINT
	}
*/
	Segment(typename _Parameter<THost>::Type _host, typename Iterator<THost, Standard>::Type _end):
		data_host(_toPointer(_host)),
		data_end_position(position(_end, _host))
	{
SEQAN_CHECKPOINT
	}

/*
	Segment(Segment const & _other):
		data_host(_other.data_host),
		data_end_position(_other.data_end_position)
	{
SEQAN_CHECKPOINT
	}
*/
	template <typename THost2, typename TSpec2>
	Segment(Segment<THost2, TSpec2> const & _other):
		data_host(_toPointer(host(_other))),
		data_end_position(endPosition(_other))
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
		return begin(host(me), Standard());
	}
	friend inline typename Iterator<Segment const, Standard>::Type 
	begin(Segment const & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return begin(host(me), Standard());
	}

//____________________________________________________________________________

	friend inline typename Position<Segment const>::Type 
	beginPosition(Segment const & /*me*/)
	{
SEQAN_CHECKPOINT
		return 0;
	}
	friend inline typename Position<Segment>::Type 
	beginPosition(Segment & /*me*/)
	{
SEQAN_CHECKPOINT
		return 0;
	}

//____________________________________________________________________________

	template <typename TIterator>
	friend inline void
	setBegin(Segment &, TIterator)
	{
	}

//____________________________________________________________________________

	friend inline typename Iterator<Segment, Standard>::Type 
	end(Segment & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return begin(host(me), Standard()) + me.data_end_position;
	}
	friend inline typename Iterator<Segment const, Standard>::Type 
	end(Segment const & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return begin(host(me), Standard()) + me.data_end_position;
	}

//____________________________________________________________________________

/* //unnoetig
	friend inline void 
	setEnd(Segment & me)
	{
SEQAN_CHECKPOINT
		me.data_end_position = length(host(me));
	}
*/
	friend inline void 
	setEndPosition(Segment & me, typename Position<Segment>::Type new_end)
	{
SEQAN_CHECKPOINT
		me.data_end_position = new_end;
	}

	friend inline void 
	setEnd(Segment & me, typename Iterator<Segment, Standard>::Type new_end)
	{
SEQAN_CHECKPOINT
		me.data_end_position = new_end - begin(host(me));//, Standard());
	}

	friend inline void 
	setEnd(typename Iterator<Segment, Rooted>::Type new_end)
	{
SEQAN_CHECKPOINT
		container(new_end).data_end_position = hostIterator(new_end) - begin(host(container(new_end)));//, Standard());
	}

//____________________________________________________________________________

	friend inline void 
	_setLength(
		Segment & me, 
		typename Size<THost>::Type new_length)
	{
SEQAN_CHECKPOINT
		me.data_end_position = new_length;
	}

//____________________________________________________________________________

	friend inline typename Position<Segment>::Type 
	endPosition(Segment & me)
	{
SEQAN_CHECKPOINT
		return me.data_end_position;
	}
	friend inline typename Position<Segment const>::Type 
	endPosition(Segment const & me)
	{
SEQAN_CHECKPOINT
		return me.data_end_position;
	}

//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Prefix:
..summary:Prefix sequence type.
..signature:Prefix<T>::Type
..param.T:A sequence type.
...type:Class.String
..returns.param.Type:The prefix type.
..see:Spec.PrefixSegment
..see:Metafunction.Infix
*/

struct InfixSegment;
struct SuffixSegment;

template <typename THost>
struct Prefix
{
	typedef Segment<THost, PrefixSegment> Type;
};

template <typename THost>
struct Prefix< Segment<THost, InfixSegment> >
{
	typedef Segment<THost, InfixSegment> Type;
};
template <typename THost>
struct Prefix< Segment<THost, SuffixSegment> >
{
	typedef Segment<THost, InfixSegment> Type;
};
template <typename THost>
struct Prefix< Segment<THost, PrefixSegment> >
{
	typedef Segment<THost, PrefixSegment> Type;
};

template <typename THost, typename TSpec>
struct Prefix< Segment<THost, TSpec> const >:
	Prefix< Segment<THost, TSpec> > {};


//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TPosition>
inline void
set(Segment<THost, PrefixSegment> & me,
	THost & host_,
	TPosition end_)
{
SEQAN_CHECKPOINT
	setHost(me, host_);
	setEnd(me, end_);
}
//____________________________________________________________________________

template <typename THost>
inline void
set(Segment<THost, PrefixSegment> & me,
	THost & host_)
{
SEQAN_CHECKPOINT
	setHost(me, host_);
	setEnd(me, end(host_));
}

//____________________________________________________________________________

template <typename THost, typename TSpec>
inline void
set(Segment<THost, PrefixSegment> & me,
	Segment<THost, TSpec> & source)
{
SEQAN_CHECKPOINT
	setHost(me, host(source));
	setEndPosition(me, endPosition(source));
}

template <typename THost, typename TSpec>
inline void
set(Segment<THost, PrefixSegment> & me,
	Segment<THost, TSpec> const & source)
{
SEQAN_CHECKPOINT
	setHost(me, host(source));
	setEndPosition(me, endPosition(source));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atBegin(Segment<THost, PrefixSegment> const & segment)
{
SEQAN_CHECKPOINT
	return (endPosition(segment) == length(host(segment)));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atEnd(Segment<THost, PrefixSegment> const & segment)
{
SEQAN_CHECKPOINT
	return (endPosition(segment) == 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline void
goBegin(Segment<THost, PrefixSegment> & segment,
		THost &)
{
SEQAN_CHECKPOINT
	goBegin(segment);
}

template <typename THost>
inline void
goBegin(Segment<THost, PrefixSegment> & segment)
{
	setEnd(segment);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline void
goEnd(Segment<THost, PrefixSegment> & segment,
	  THost &)
{
SEQAN_CHECKPOINT
	goEnd(segment);
}

template <typename THost>
inline void
goEnd(Segment<THost, PrefixSegment> & segment)
{
	setEnd(segment, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, PrefixSegment> &
operator ++(Segment<THost, PrefixSegment> & segment)
{
	setEnd(segment, endPosition(segment) - 1);
	return segment;
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, PrefixSegment> &
operator --(Segment<THost, PrefixSegment> & segment)
{
	setEnd(segment, endPosition(segment) + 1);
	return segment;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.prefix:
..cat:Containers
..summary:Creates prefix object.
..signature:prefix(host, end)
..param.host:The complete sequence.
...type:Class.String
...type:Adaption.char array
..param.end:Position or iterator behind the last element of the segment.
...type:Metafunction.Position
...type:Metafunction.Iterator
..returns:The prefix of $host that begins at $begin$.
...remarks:The type of the prefix is given by @Metafunction.Prefix@.
..remarks:Notational sugar.
..see:Spec.PrefixSegment
..see:Function.suffix
..see:Function.infix
*/

template <typename T, typename TPosEnd>
inline typename Prefix<T>::Type
prefix(T & t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
	return typename Prefix<T>::Type(t, pos_end);
}

template <typename T, typename TPosEnd>
inline typename Prefix<T *>::Type
prefix(T * t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
	return typename Prefix<T *>::Type (t, pos_end);
}

//////////////////////////////////////////////////////////////////////////////
// A prefix of a prefix -> is a prefix
template <typename T, typename TPosEnd>
inline typename Prefix<Segment<T, PrefixSegment> >::Type
prefix(Segment<T, PrefixSegment> & t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
	return typename Prefix<Segment<T, PrefixSegment> >::Type (
		host(t), 
		beginPosition(t) + pos_end);
}
template <typename T, typename TPosEnd>
inline typename Prefix<Segment<T, PrefixSegment> const>::Type
prefix(Segment<T, PrefixSegment> const & t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
	return typename Prefix<Segment<T, PrefixSegment> const>::Type (
		host(t), 
		beginPosition(t) + pos_end);
}

//////////////////////////////////////////////////////////////////////////////
// A prefix of an infix -> is an infix
template <typename T, typename TPosEnd>
inline typename Prefix<Segment<T, InfixSegment> >::Type
prefix(Segment<T, InfixSegment> & t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
	return typename Prefix<Segment<T, InfixSegment> >::Type (
		host(t), 
		beginPosition(t),
		beginPosition(t) + pos_end);
}
template <typename T, typename TPosEnd>
inline typename Prefix<Segment<T, InfixSegment> const>::Type
prefix(Segment<T, InfixSegment> const & t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
	return typename Prefix<Segment<T, InfixSegment> const>::Type (
		host(t), 
		beginPosition(t),
		beginPosition(t) + pos_end);
}


//////////////////////////////////////////////////////////////////////////////
// A prefix of an suffix -> is an infix
template <typename T, typename TPosEnd>
inline typename Prefix<Segment<T, SuffixSegment> >::Type
prefix(Segment<T, SuffixSegment> & t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
	return typename Prefix<Segment<T, SuffixSegment> >::Type (
		host(t), 
		beginPosition(t),
		beginPosition(t) + pos_end);
}
template <typename T, typename TPosEnd>
inline typename Prefix<Segment<T, SuffixSegment> const>::Type
prefix(Segment<T, SuffixSegment> const & t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
	return typename Prefix<Segment<T, SuffixSegment> const>::Type (
		host(t), 
		beginPosition(t),
		beginPosition(t) + pos_end);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
