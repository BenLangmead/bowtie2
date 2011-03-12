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
  $Id: segment_infix.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEGMENT_INFIX_H
#define SEQAN_HEADER_SEGMENT_INFIX_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// InfixSegment
//////////////////////////////////////////////////////////////////////////////


/**
.Spec.InfixSegment:
..cat:Segments
..summary:An arbitrary segment.
..general:Class.Segment
..signature:Segment<THost, InfixSegment>
..param.THost:Type of the whole sequence.
...text:Instances of $Segment<THost, InfixSegment>$ are infixes of $THost$ objects.
...remarks:Use @Metafunction.Host@ to get the host type for a given class.
..remarks.note:Since the appropriate segment type depends on the host sequence type, 
	it is recommended to use the metafunction @Metafunction.Infix@ instead of explicitely 
	choose a specialization of @Class.Segment@.
..see:Metafunction.Infix
*/

template <typename THost_>
class Segment<THost_, InfixSegment>
{
protected:
	typedef typename Host<Segment>::Type THost;

	typename _Pointer<THost>::Type data_host;
	typename Position<THost>::Type data_begin_position;
	typename Position<THost>::Type data_end_position;


//____________________________________________________________________________

public:

/**
.Memfunc.InfixSegment#Segment:
..class:Spec.InfixSegment
..summary:Constructor
..signature:Segment<THost, InfixSegment> ()
..signature:Segment<THost, InfixSegment> (infix)
..signature:Segment<THost, InfixSegment> (host [, begin, end])
..param.infix:Other infix object. (copy constructor)
..param.host:The whole sequence.
..param.begin:Position/iterator in $host$ of the first item in segment.
...type:Metafunction.Position.$Position<THost>::Type$
...type:Metafunction.Iterator.$Iterator<THost>::Type$
..param.end:Position/iterator behind the end of the segment.
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
...text:$begin$ and $end$ must be valid positions/iterators in $host$.
If $begin$ und $end$ are omitted, the infix segment corresponding to
the first character of $host$ is constructed.
This is the same segment that is returned by @Function.goBegin@.
*/
	Segment():
		data_begin_position(0),
		data_end_position(0)
	{
SEQAN_CHECKPOINT
	}

	Segment(typename _Parameter<THost>::Type _host):
		data_host(_toPointer(_host)),
		data_begin_position(0),
		data_end_position(1)
	{
SEQAN_CHECKPOINT
	}

	Segment(typename _Parameter<THost>::Type _host, typename Position<THost>::Type _begin_index, typename Position<THost>::Type _end_index):
		data_host(_toPointer(_host)),
		data_begin_position(_begin_index),
		data_end_position(_end_index)
	{
SEQAN_CHECKPOINT
	}
/*
	Segment(typename _Parameter<THost>::Type _host, typename Iterator<THost, Rooted>::Type _begin, typename Iterator<THost, Rooted>::Type _end):
		data_host(_toPointer(_host)),
		data_begin_position(position(_begin)),
		data_end_position(position(_end))
	{
SEQAN_CHECKPOINT
	}
*/
	Segment(typename _Parameter<THost>::Type _host, typename Iterator<THost, Standard>::Type _begin, typename Iterator<THost, Standard>::Type _end):
		data_host(_toPointer(_host)),
		data_begin_position(position(_begin, _host)),
		data_end_position(position(_end, _host))
	{
SEQAN_CHECKPOINT
	}
	template <typename THost2, typename TSpec2>
	Segment(Segment<THost2, TSpec2> const & _other):
		data_host(_toPointer(host(_other))),
		data_begin_position(beginPosition(_other)),
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

///Function.host.param.object.type:Class.Segment

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

/**
.Function.setHost:
..summary:Sets the host of an object.
..cat:Dependent Objects
..signature:setHost(object, host)
..param.object:The object that will get a new host.
...type:Class.Segment
..param.host:The new host.
..remarks:After this operation, $object$ depends on $host$.
...text:Note that setting the host can invalidate $object$.
For example, if one changes the host of a @Class.Segment@ object, it is possible
that begin- and end-position of the segment does not fit into the new host sequence.
..see:Function.host
*/
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

///.Function.begin.param.object.type:Class.Segment

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

///.Function.beginPosition.param.object.type:Class.Segment

	friend inline typename Position<Segment>::Type 
	beginPosition(Segment & me)
	{
SEQAN_CHECKPOINT
		return me.data_begin_position;
	}
	friend inline typename Position<Segment const>::Type 
	beginPosition(Segment const & me)
	{
SEQAN_CHECKPOINT
		return me.data_begin_position;
	}

//____________________________________________________________________________

/**
.Function.setBegin:
..summary:Sets begin of object in host.
..cat:Dependent Objects
..signature:setBegin(object, new_begin)
..param.object:An object.
...type:Spec.InfixSegment
...type:Spec.SuffixSegment
..param.new_begin:iterator to the new first item in $host(object)$ that belongs of $object$.
...type:Metafunction.Iterator
..see:Function.begin
..see:Function.beginPosition
*/
	template <typename TIterator>
	friend inline void 
	setBegin(Segment & me, TIterator new_begin)
	{
SEQAN_CHECKPOINT
		me.data_begin_position = new_begin - begin(host(me));//, Standard());
	}


//____________________________________________________________________________

/**
.Function.setBeginPosition:
..summary:Sets begin position of object in host.
..cat:Dependent Objects
..signature:setBeginPosition(object, new_begin)
..param.object:An object.
...type:Spec.InfixSegment
...type:Spec.SuffixSegment
..param.new_begin:position of the new first item in $host(object)$ that belongs of $object$.
...type:Metafunction.Position
..see:Function.begin
..see:Function.beginPosition
..see:Function.setBegin
*/

	template <typename TPosition>
	friend inline void 
	setBeginPosition(Segment & me, TPosition new_begin)
	{
SEQAN_CHECKPOINT
		me.data_begin_position = new_begin;
	}

//____________________________________________________________________________

///.Function.begin.param.object.type:Class.Segment

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

///.Function.endPosition.param.object.type:Class.Segment

	friend inline typename Position<Segment>::Type 
	endPosition(Segment & me)
	{
SEQAN_CHECKPOINT
		return me.data_end_position;
	}
	friend inline typename Position<Segment>::Type 
	endPosition(Segment const & me)
	{
SEQAN_CHECKPOINT
		return me.data_end_position;
	}

//____________________________________________________________________________

/**
.Function.setEnd:
..summary:Sets end of object in host.
..cat:Dependent Objects
..signature:setEnd(object, new_end)
..param.object:An object.
...type:Spec.InfixSegment
...type:Spec.PrefixSegment
..param.new_end:Iterator behind the last item in $host(object)$ belongs of $object$.
...type:Metafunction.Iterator
..see:Function.end
..see:Function.endPosition
..see:Function.setBegin
*/

	template <typename TIterator>
	friend inline void 
	setEnd(Segment & me, TIterator new_end)
	{
SEQAN_CHECKPOINT
		me.data_end_position = new_end - begin(host(me));//, Standard());
	}

/* //unnoetig
	friend inline void 
	setEnd(Segment & me)
	{
SEQAN_CHECKPOINT
		setEnd(me, end(host(me)));
	}
*/
//____________________________________________________________________________


/**
.Function.setEndPosition:
..summary:Sets begin position of object in host.
..cat:Dependent Objects
..signature:setEndPosition(object, new_end)
..param.object:An object.
...type:Spec.InfixSegment
...type:Spec.PrefixSegment
..param.new_end:position behind the last item in $host(object)$ that belongs of $object$.
...type:Metafunction.Position
..see:Function.end
..see:Function.endPosition
..see:Function.setBeginPosition
..see:Function.setEnd
*/

	template <typename TPosition>
	friend inline void 
	setEndPosition(Segment & me, TPosition new_end)
	{
SEQAN_CHECKPOINT
		me.data_end_position = new_end;
	}

//____________________________________________________________________________

	friend inline void 
	_setLength(
		Segment & me, 
		typename Size<THost>::Type new_length)
	{
SEQAN_CHECKPOINT
		me.data_end_position = me.data_begin_position + new_length;
	}

//____________________________________________________________________________

};

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Infix:
..summary:Infix sequence type.
..signature:Infix<T>::Type
..param.T:A sequence type.
...type:Class.String
..returns.param.Type:The infix type.
..see:Spec.InfixSegment
*/

template <typename THost>
struct Infix
{
	typedef Segment<THost, InfixSegment> Type;
};

template <typename THost, typename TSpec>
struct Infix< Segment<THost, TSpec> >
{
	typedef Segment<THost, InfixSegment> Type;
};

template <typename THost, typename TSpec>
struct Infix< Segment<THost, TSpec> const >:
	Infix< Segment<THost, TSpec> > {};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TPosition1, typename TPosition2>
inline void
set(Segment<THost, InfixSegment> & me,
	THost & host_,
	TPosition1 begin_,
	TPosition2 end_)
{
SEQAN_CHECKPOINT
	setHost(me, host_);
	setBegin(me, begin_);
	setEnd(me, end_);
}
//____________________________________________________________________________

template <typename THost>
inline void
set(Segment<THost, InfixSegment> & me,
	THost & host_)
{
SEQAN_CHECKPOINT
	setHost(me, host_);
	setBegin(me, begin(host_, Standard()));
	setEnd(me, end(host_, Standard()));
}
template <typename THost>
inline void
set(Segment<THost, InfixSegment> & me,
	THost const & host_)
{
SEQAN_CHECKPOINT
	setHost(me, host_);
	setBegin(me, begin(host_, Standard()));
	setEnd(me, end(host_, Standard()));
}

//____________________________________________________________________________

template <typename THost, typename TSpec>
inline void
set(Segment<THost, InfixSegment> & me,
	Segment<THost, TSpec> & source)
{
SEQAN_CHECKPOINT
	setHost(me, host(source));
	setBeginPosition(me, beginPosition(source));
	setEndPosition(me, endPosition(source));
}
template <typename THost, typename TSpec>
inline void
set(Segment<THost, InfixSegment> & me,
	Segment<THost, TSpec> const & source)
{
SEQAN_CHECKPOINT
	setHost(me, host(source));
	setBeginPosition(me, beginPosition(source));
	setEndPosition(me, endPosition(source));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atBegin(Segment<THost, InfixSegment> & segment)
{
SEQAN_CHECKPOINT
	return (beginPosition(segment) == endPosition(segment));
}
template <typename THost>
inline bool
atBegin(Segment<THost, InfixSegment> const & segment)
{
SEQAN_CHECKPOINT
	return (beginPosition(segment) == endPosition(segment));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atEnd(Segment<THost, InfixSegment> & segment)
{
SEQAN_CHECKPOINT
	return (endPosition(segment) - beginPosition(segment)) > length(host(segment));
}
template <typename THost>
inline bool
atEnd(Segment<THost, InfixSegment> const & segment)
{
SEQAN_CHECKPOINT
	return (endPosition(segment) - beginPosition(segment)) > length(host(segment));
}


//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline void
goBegin(Segment<THost, InfixSegment> & segment)
{
SEQAN_CHECKPOINT
	setBeginPosition(segment, 0);
	setEndPosition(segment, 1);
}
template <typename THost, typename THost2>
inline void
goBegin(Segment<THost, InfixSegment> & segment,
		THost2 &)
{
	goBegin(segment);
}
template <typename THost, typename THost2>
inline void
goBegin(Segment<THost, InfixSegment> & segment,
		THost2 const &)
{
	goBegin(segment);
}

//////////////////////////////////////////////////////////////////////////////


template <typename THost>
inline void
goEnd(Segment<THost, InfixSegment> & segment)
{
SEQAN_CHECKPOINT
	setBeginPosition(segment, 0);
	setEndPosition(segment, length(host(segment)));
}
template <typename THost, typename THost2>
inline void
goEnd(Segment<THost, InfixSegment> & segment,
	  THost2 &)
{
	goEnd(segment);
}
template <typename THost, typename THost2>
inline void
goEnd(Segment<THost, InfixSegment> & segment,
	  THost2 const &)
{
	goEnd(segment);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, InfixSegment> &
operator ++(Segment<THost, InfixSegment> & segment)
{
SEQAN_CHECKPOINT
	if (endPosition(segment) == length(host(segment)))
	{
		setEndPosition(segment, endPosition(segment) - beginPosition(segment) + 1);
		setBeginPosition(segment, 0);
	}
	else
	{
		setBeginPosition(segment, beginPosition(segment) + 1);
		setEndPosition(segment, endPosition(segment) + 1);
	}
	return segment;
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, InfixSegment> &
operator --(Segment<THost, InfixSegment> & segment)
{
SEQAN_CHECKPOINT
	if (!beginPosition(segment))
	{
		typename Size<THost>::Type host_length = length(host(segment));

		setBeginPosition(segment, host_length - endPosition(segment) + beginPosition(segment) + 1);
		setEndPosition(segment, host_length);
	}
	else
	{
		setBeginPosition(segment, beginPosition(segment) - 1);
		setEndPosition(segment, endPosition(segment) - 1);
	}
	return segment;
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec, typename TPos>
inline typename Reference< Segment<THost, TSpec> >::Type 
value(Segment<THost, TSpec> & me, 
	  TPos pos)
{
SEQAN_CHECKPOINT
	return *(begin(me, Standard()) + pos);
}

template <typename THost, typename TSpec, typename TPos>
inline typename Reference< Segment<THost, TSpec> const >::Type 
value(Segment<THost, TSpec> const & me, 
	  TPos pos)
{
SEQAN_CHECKPOINT
	return *(begin(me, Standard()) + pos);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.infix:
..cat:Containers
..summary:Creates infix object.
..signature:infix(host, begin, end)
..param.host:The complete sequence.
...type:Class.String
...type:Adaption.char array
..param.begin:Position or iterator of the first element of the segment.
...type:Metafunction.Position
...type:Metafunction.Iterator
..param.end:Position or iterator behind the last element of the segment.
...remarks:$end$ must have the same type as $begin$.
..returns:The infix of $host$ between $begin$ and $end-1$.
...remarks:The type of the infix is given by @Metafunction.Infix@.
..remarks:Notational sugar.
..see:Spec.InfixSegment
*/

template <typename T, typename TPosBegin, typename TPosEnd>
inline typename Infix<T>::Type
infix(T & t, TPosBegin pos_begin, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
	return typename Infix<T>::Type(t, pos_begin, pos_end);
}

template <typename T, typename TPosBegin, typename TPosEnd>
inline typename Infix<T *>::Type
infix(T * t, TPosBegin pos_begin, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
	return typename Infix<T *>::Type (t, pos_begin, pos_end);
}

template <typename T, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<Segment<T, TSpec> >::Type
infix(Segment<T, TSpec> & t, TPosBegin pos_begin, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
	return typename Infix<Segment<T, TSpec> >::Type (
		host(t), 
		beginPosition(t) + pos_begin, 
		beginPosition(t) + pos_end);
}

template <typename T, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<Segment<T, TSpec> const>::Type
infix(Segment<T, TSpec> const & t, TPosBegin pos_begin, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
	return typename Infix<Segment<T, TSpec> const>::Type (
		host(t), 
		beginPosition(t) + pos_begin, 
		beginPosition(t) + pos_end);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.infixWithLength:
..cat:Containers
..summary:Creates infix object.
..signature:infixWithLength(host, begin, length)
..param.host:The complete sequence.
...type:Class.String
...type:Adaption.char array
..param.begin:Position or iterator of the first element of the segment.
...type:Metafunction.Position
...type:Metafunction.Iterator
..param.length:Length of the returned infix.
..returns:The infix of $host$ between $begin$ and $begin+length-1$.
...remarks:The type of the infix is given by @Metafunction.Infix@.
..remarks:Notational sugar.
..see:Spec.InfixSegment
*/

template <typename T, typename TPosBegin, typename TSize>
inline typename Infix<T>::Type
infixWithLength(T & t, TPosBegin pos_begin, TSize length)
{
SEQAN_CHECKPOINT
	return typename Infix<T>::Type(t, pos_begin, pos_begin + length);
}

template <typename T, typename TPosBegin, typename TSize>
inline typename Infix<T *>::Type
infixWithLength(T * t, TPosBegin pos_begin, TSize length)
{
SEQAN_CHECKPOINT
	return typename Infix<T *>::Type (t, pos_begin, pos_begin + length);
}

template <typename T, typename TSpec, typename TPosBegin, typename TSize>
inline typename Infix<Segment<T, TSpec> >::Type
infixWithLength(Segment<T, TSpec> & t, TPosBegin pos_begin, TSize length)
{
SEQAN_CHECKPOINT
	return typename Infix<Segment<T, TSpec> >::Type (
		host(t), 
		beginPosition(t) + pos_begin, 
		beginPosition(t) + pos_begin + length);
}

template <typename T, typename TSpec, typename TPosBegin, typename TSize>
inline typename Infix<Segment<T, TSpec> const>::Type
infixWithLength(Segment<T, TSpec> const & t, TPosBegin pos_begin, TSize length)
{
SEQAN_CHECKPOINT
	return typename Infix<Segment<T, TSpec> const>::Type (
		host(t), 
		beginPosition(t) + pos_begin, 
		beginPosition(t) + pos_begin + length);
}

//////////////////////////////////////////////////////////////////////////////
//setBegin


template <typename TIterator>
inline void 
setBegin(TIterator new_begin)
{
SEQAN_CHECKPOINT
	setBegin(container(new_begin), hostIterator(new_begin));
}


//////////////////////////////////////////////////////////////////////////////
//setEnd

template <typename TIterator>
inline void 
setEnd(TIterator new_end)
{
SEQAN_CHECKPOINT
	setEnd(container(new_end), new_end);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
