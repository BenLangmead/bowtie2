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
  $Id: string_packed.h,v 1.2 2009/02/19 01:51:23 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEQUENCE_PACKED_H
#define SEQAN_HEADER_SEQUENCE_PACKED_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////

template <typename THostspec = Alloc<> >
struct Packed;


//////////////////////////////////////////////////////////////////////////////
/**
.Spec.Packed String:
..cat:Strings
..general:Class.String
..summary:A string that stores as many values in one machine word as possible.
..signature:String<TValue, Packed<THostspec> >
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.THostspec:The specializing type.
...remarks:This is the specialization of the host string that is used for storing the packed values.
...default:@Spec.Alloc String.Alloc<>@
*/

/*???TODO Optimierungsmöglichkeiten:
- _clearSpace kopiert Zeichenweise im Packed-String, und nicht im Host-String
- _clearSpace verwendet resize, um den Host zu vergrößern, d.h. der Inhalt wird eventuell doppelt kopiert.
*/

//////////////////////////////////////////////////////////////////////////////
//Rooted expandable string
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec>
class String<TValue, Packed<THostspec> >
{
protected:
	typedef typename Host<String>::Type THost;
	typedef typename Size<String>::Type TSize;

	THost data_host;
	TSize data_length;

//____________________________________________________________________________

public:
	String():
		data_length(0)
	{
SEQAN_CHECKPOINT
	}

	template <typename TSource>
	String(TSource & source):
		data_length(0)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
	}
	template <typename TSource>
	String(TSource const & source):
		data_length(0)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
	}
	String(String const & source):
		data_length(0)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
	}

	template <typename TSource>
	String & operator =(TSource const & source)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
		return *this;
	}
	String & operator =(String const & source)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
		return *this;
	}

	~String()
	{
SEQAN_CHECKPOINT
	}


//____________________________________________________________________________

	template <typename TPos>
	inline typename Reference<String>::Type
	operator [] (TPos pos)
	{
SEQAN_CHECKPOINT
		return value(*this, pos);
	}

	template <typename TPos>
	inline typename Reference<String const>::Type
	operator [] (TPos pos) const
	{
SEQAN_CHECKPOINT
		return value(*this, pos);
	}

//____________________________________________________________________________

///.Function.host.param.object.type:Spec.Packed String

	friend inline THost &
	host(String & me)
	{
SEQAN_CHECKPOINT
		return me.data_host;
	}

	friend inline THost const &
	host(String const & me)
	{
SEQAN_CHECKPOINT
		return me.data_host;
	}

//____________________________________________________________________________

	friend inline TSize
	length(String & me)
	{
SEQAN_CHECKPOINT
		return me.data_length;
	}

	friend inline TSize
	length(String const & me)
	{
SEQAN_CHECKPOINT
		return me.data_length;
	}

//____________________________________________________________________________

	friend inline void
	_setLength(
		String & me,
		TSize new_length)
	{
SEQAN_CHECKPOINT
		me.data_length = new_length;
		_setLength(host(me), _PackedConsts<String>::toHostLength(new_length));
	}

//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec>
struct DefaultOverflowImplicit<String<TValue, Packed<THostspec> > >:
	DefaultOverflowImplicit< typename Host<String<TValue, Packed<THostspec> > >::Type >
{
};
template <typename TValue, typename THostspec>
struct DefaultOverflowImplicit<String<TValue, Packed<THostspec> > const>:
	DefaultOverflowImplicit< typename Host<String<TValue, Packed<THostspec> > const>::Type >
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec>
struct DefaultOverflowExplicit<String<TValue, Packed<THostspec> > >:
	DefaultOverflowExplicit< typename Host<String<TValue, Packed<THostspec> > >::Type >
{
};
template <typename TValue, typename THostspec>
struct DefaultOverflowExplicit<String<TValue, Packed<THostspec> > const>:
	DefaultOverflowExplicit< typename Host<String<TValue, Packed<THostspec> > const>::Type >
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec>
struct IsContiguous<String<TValue, Packed<THostspec> > >
{
    typedef False Type;
	enum { VALUE = false };
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Host.param.T.type:Spec.Packed String

template <typename TValue, typename THostspec>
struct Host<String<TValue, Packed<THostspec> > >
{
	typedef String<unsigned int, THostspec> Type;
};
template <typename TValue, typename THostspec>
struct Host<String<TValue, Packed<THostspec> > const>
{
	typedef String<unsigned int, THostspec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec>
struct GetValue<String<TValue, Packed<THostspec> > >:
	Value<String<TValue, Packed<THostspec> > >
{
};
template <typename TValue, typename THostspec>
struct GetValue<String<TValue, Packed<THostspec> > const>:
	Value<String<TValue, Packed<THostspec> > const>
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec>
struct Reference<String<TValue, Packed<THostspec> > >
{
	typedef typename Iterator<String<TValue, Packed<THostspec> >, Standard>::Type TIterator;
	typedef Proxy<IteratorProxy<TIterator> > Type;
};
template <typename TValue, typename THostspec>
struct Reference<String<TValue, Packed<THostspec> > const>
{
	typedef typename Iterator<String<TValue, Packed<THostspec> > const, Standard>::Type TIterator;
	typedef Proxy<IteratorProxy<TIterator> > Type;
};

//////////////////////////////////////////////////////////////////////////////
/*
template <typename TValue, typename THostspec>
struct Size<String<TValue, Packed<THostspec> > >
{
	typedef __int64 Type;
};
template <typename TValue, typename THostspec>
struct Size<String<TValue, Packed<THostspec> > const>
{
	typedef __int64 Type;
};
//*/

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Compute Bit Values
//////////////////////////////////////////////////////////////////////////////

template <typename TPackedContainer>
struct _PackedConsts
{
	typedef typename Value<TPackedContainer>::Type TValue;
	typedef typename Host<TPackedContainer>::Type THost;
	typedef typename Value<THost>::Type THostValue;

	enum
	{
		BITS_PER_VALUE = BitsPerValue<TValue>::VALUE,
		BITS_PER_HOST_VALUE = BitsPerValue<THostValue>::VALUE,
		VALUES_PER_WORD = (BITS_PER_VALUE > BITS_PER_HOST_VALUE) ? 1 : (BITS_PER_HOST_VALUE / BITS_PER_VALUE),
		VALUE_MASK = (1 << BITS_PER_VALUE) - 1,
		MAX_BIT_POS = (VALUES_PER_WORD - 1) * BITS_PER_VALUE
	};

	static typename Size<THost>::Type
	toHostLength(typename Size<TPackedContainer>::Type len)
	{
		return (len + VALUES_PER_WORD - 1) / VALUES_PER_WORD;
	}
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Temporary Copy
//////////////////////////////////////////////////////////////////////////////
//note: this works only, if the copy assignment is done without using _TempCopy

template <typename TValue, typename THostspec>
struct _TempCopy<String<TValue, Packed<THostspec> > >
{
	typedef String<TValue, Packed<THostspec> > Type;
};

//////////////////////////////////////////////////////////////////////////////
//optimized variant for copy assignment. The host sequence is copied instead of
//copying the packed string value by value

template <typename TTarget, typename TSource, typename TTag>
inline void
_assign_copy_packed_string(TTarget & target,
						   TSource & source,
						   Tag<TTag> const tag)
{
	typedef typename Size<TTarget>::Type TSize2;

	assign(host(target), host(source), tag);
	TSize2 new_length_limit = length(host(target)) * _PackedConsts<TTarget>::VALUES_PER_WORD;
	TSize2 new_length = length(source);
	if (new_length > new_length_limit)
	{
		new_length = new_length_limit;
	}
	_setLength(target, new_length);
}
template <typename TTarget, typename TSource, typename TSize, typename TTag>
inline void
_assign_copy_packed_string(TTarget & target,
						   TSource & source,
						   TSize limit,
						   Tag<TTag> const tag)
{
	typedef typename Size<TTarget>::Type TSize2;

	TSize2 host_limit = _PackedConsts<TTarget>::toHostLength(limit);
	assign(host(target), host(source), host_limit, tag);
	TSize2 new_length_limit = length(host(target)) * _PackedConsts<TTarget>::VALUES_PER_WORD;
	TSize2 new_length = length(source);
	if (new_length > new_length_limit)
	{
		new_length = new_length_limit;
	}
	if (new_length > limit)
	{
		new_length = limit;
	}
	_setLength(target, new_length);
}
//____________________________________________________________________________

template <typename TValue, typename THostspec, typename TTag>
inline void
assign(String<TValue, Packed<THostspec> > & target,
	   String<TValue, Packed<THostspec> > & source,
	   Tag<TTag> const tag)
{
	_assign_copy_packed_string(target, source, tag);
}
template <typename TValue, typename THostspec, typename TTag>
inline void
assign(String<TValue, Packed<THostspec> > & target,
	   String<TValue, Packed<THostspec> > const & source,
	   Tag<TTag> const tag)
{
	_assign_copy_packed_string(target, source, tag);
}

template <typename TValue, typename THostspec, typename TSize, typename TTag>
void assign(String<TValue, Packed<THostspec> > & target,
			String<TValue, Packed<THostspec> > & source,
			TSize limit,
			Tag<TTag> const tag)
{
	_assign_copy_packed_string(target, source, limit, tag);
}
template <typename TValue, typename THostspec, typename TSize, typename TTag>
void assign(String<TValue, Packed<THostspec> > & target,
			String<TValue, Packed<THostspec> > const & source,
			TSize limit,
			Tag<TTag> const tag)
{
	_assign_copy_packed_string(target, source, limit, tag);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Function
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec>
inline void const *
id(String<TValue, Packed<THostspec> > const & me)
{
SEQAN_CHECKPOINT
	return id(host(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec, typename TPos, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type
iter(String<TValue, Packed<THostspec> > & me,
	 TPos pos_,
	 Tag<TTag> const)
{
SEQAN_CHECKPOINT
	typedef typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type TIterator;
	return TIterator(me, pos_);
}
template <typename TValue, typename THostspec, typename TPos, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type
iter(String<TValue, Packed<THostspec> > const & me,
	 TPos pos_,
	 Tag<TTag> const)
{
SEQAN_CHECKPOINT
	typedef typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type TIterator;
	return TIterator(me, pos_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type
begin(String<TValue, Packed<THostspec> > & me,
	  Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return iter(me, 0, tag_);
}
template <typename TValue, typename THostspec, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type
begin(String<TValue, Packed<THostspec> > const & me,
	  Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return iter(me, 0, tag_);
}

//////////////////////////////////////////////////////////////////////////////
// end
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type
end(String<TValue, Packed<THostspec> > & me,
	Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return iter(me, length(me), tag_);
}
template <typename TValue, typename THostspec, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type
end(String<TValue, Packed<THostspec> > const & me,
	Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return iter(me, length(me), tag_);
}

//////////////////////////////////////////////////////////////////////////////
// value
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec, typename TPos>
inline typename Reference<String<TValue, Packed<THostspec> > >::Type
value(String<TValue, Packed<THostspec> > & me,
	  TPos pos)
{
SEQAN_CHECKPOINT

	return *iter(me, pos, Standard());
}
template <typename TValue, typename THostspec, typename TPos>
inline typename Reference<String<TValue, Packed<THostspec> > const>::Type
value(String<TValue, Packed<THostspec> > const & me,
	  TPos pos)
{
SEQAN_CHECKPOINT

	return *iter(me, pos, Standard());
}


//////////////////////////////////////////////////////////////////////////////
// capacity
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec>
inline typename Size<String<TValue, Packed<THostspec> > const>::Type
capacity(String<TValue, Packed<THostspec> > const & me)
{
SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Packed<THostspec> > const>::Type TSize;
	TSize len = capacity(host(me));
	len *= _PackedConsts<String<TValue, Packed<THostspec> > >::VALUES_PER_WORD;
	return len;
}

//////////////////////////////////////////////////////////////////////////////
// clear
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec>
inline void
clear(String<TValue, Packed<THostspec> > & me)
{
SEQAN_CHECKPOINT
	clear(host(me));
	_setLength(me, 0);
}

//////////////////////////////////////////////////////////////////////////////
// _clearSpace
//////////////////////////////////////////////////////////////////////////////

//implementation for all expand tags other than "limit"
template <typename TExpand>
struct _ClearSpace_String_Packed_
{
	template <typename T>
	static inline typename Size<T>::Type
	_clearSpace_(
		T & seq,
		typename Size<T>::Type size)
	{
SEQAN_CHECKPOINT
		typedef typename Size<T>::Type TSize;
		TSize wanted_host_length = _PackedConsts<T>::toHostLength(size);
		TSize new_host_length = resize(host(seq), wanted_host_length, TExpand());
		if (new_host_length < wanted_host_length)
		{
			size = new_host_length * _PackedConsts<T>::VALUES_PER_WORD;
		}
		_setLength(seq, size);
		return size;
	}

	template <typename T>
	static inline typename Size<T>::Type
	_clearSpace_(
		T & seq,
		typename Size<T>::Type size,
		typename Size<T>::Type limit)
	{
		if (limit < size)
		{
SEQAN_CHECKPOINT
			size = limit;
		}
		return _clearSpace_(seq, limit);
	}

	template <typename T>
	static inline typename Size<T>::Type
	_clearSpace_(
		T & seq,
		typename Size<T>::Type size,
		typename Size<T>::Type start,
		typename Size<T>::Type end)
	{
SEQAN_CHECKPOINT
		return _clearSpace_(seq, size, start, end, supremumValue<typename Size<T>::Type >());
	}

	template <typename T>
	static typename Size<T>::Type
	_clearSpace_(
		T & seq,
		typename Size<T>::Type size,
		typename Size<T>::Type start,
		typename Size<T>::Type end,
		typename Size<T>::Type limit)
	{
SEQAN_CHECKPOINT
//??? TODO: This function can be accelerated this way:
//				- move values in host
//				- avoid double moving of the rest-part if "resize" allocates a new block

		typedef typename Size<T>::Type TSize;
		typedef typename Iterator<T, Standard>::Type TIterator;

		TSize old_length = length(seq);
		TSize old_size = end - start;
		TSize wanted_new_length = old_length + size - old_size;

		if (wanted_new_length > limit)
		{
			wanted_new_length = limit;
		}

		TSize wanted_host_length = _PackedConsts<T>::toHostLength(wanted_new_length);
		TSize new_host_length = resize(host(seq), wanted_host_length, TExpand());

		TSize new_length;
		if (new_host_length < wanted_host_length)
		{
			new_length = new_host_length * _PackedConsts<T>::VALUES_PER_WORD;
			if (new_length <= start + size)
			{
				goto FINISH;
			}
			old_length = new_length - size + old_size;
		}
		else
		{
			new_length = wanted_new_length;
		}

		//move [end:right_end] to [start + size:..]
		if (old_size > size)
		{//move rest to left
			::std::copy_backward(iter(seq, end, Standard()), iter(seq, old_length, Standard()), iter(seq,  new_length, Standard()));
		}
		else
		{//move rest to right
			::std::copy(iter(seq, end, Standard()), iter(seq, old_length, Standard()), iter(seq, end + size - old_size, Standard()));
		}
FINISH:
		_setLength(seq, new_length);
		return size;
	}
/*
	template <typename T>
	static inline typename Size<T>::Type
	_clearSpace_(
		T & seq,
		typename Size<T>::Type size,
		typename Iterator<T>::Type start,
		typename Iterator<T>::Type end)
	{
SEQAN_CHECKPOINT
		typename Iterator<T>::Type seq_begin = begin(seq);
		return _clearSpace(seq, size, start - seq_begin, end - seq_begin, Insist());
	}

	template <typename T>
	static inline typename Size<T>::Type
	_clearSpace_(
		T & seq,
		typename Size<T>::Type size,
		typename Iterator<T>::Type start,
		typename Iterator<T>::Type end,
		typename Size<T>::Type limit)
	{
SEQAN_CHECKPOINT
		typename Iterator<T>::Type seq_begin = begin(seq);
		return _clearSpace(seq, size, start - seq_begin, end - seq_begin, limit, Insist());
	}
*/
};

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename THostspec, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type
_clearSpace(String<TValue, Packed<THostspec> > & me,
		typename Size< String<TValue, Packed<THostspec> > >::Type size,
		Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	return _ClearSpace_String_Packed_<Tag<TExpand> const>::_clearSpace_(me, size);
}

template<typename TValue, typename THostspec, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type
_clearSpace(String<TValue, Packed<THostspec> > & me,
		typename Size< String<TValue, Packed<THostspec> > >::Type size,
		typename Size< String<TValue, Packed<THostspec> > >::Type limit,
		Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	return _ClearSpace_String_Packed_<Tag<TExpand> const>::_clearSpace_(me, size, limit);
}

template<typename TValue, typename THostspec, typename TPosition, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type
_clearSpace(String<TValue, Packed<THostspec> > & me,
			typename Size< String<TValue, Packed<THostspec> > >::Type size,
			TPosition pos_begin,
			TPosition pos_end,
			Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	return _ClearSpace_String_Packed_<Tag<TExpand> const>::_clearSpace_(me, size, pos_begin, pos_end);
}

template<typename TValue, typename THostspec, typename TPosition, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type
_clearSpace(String<TValue, Packed<THostspec> > & me,
			typename Size< String<TValue, Packed<THostspec> > >::Type size,
			TPosition pos_begin,
			TPosition pos_end,
			typename Size< String<TValue, Packed<THostspec> > >::Type limit,
			Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	return _ClearSpace_String_Packed_<Tag<TExpand> const>::_clearSpace_(me, size, pos_begin, pos_end, limit);
}


//////////////////////////////////////////////////////////////////////////////


///.Function.reserve.param.object.type:Spec.Packed String

template <typename TValue, typename TSpec, typename _TSize, typename TExpand>
inline typename Size< String<TValue, Packed<TSpec> > >::Type
reserve(
	String<TValue, Packed<TSpec> > & seq,
	_TSize new_capacity,
	Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT

	typedef String<TValue, Packed<TSpec> > TString;
	typedef typename Size<TString>::Type TSize;
	TSize ret_value = reserve(host(seq), _PackedConsts<TString>::toHostLength(new_capacity), tag);
	return ret_value * _PackedConsts<TString>::VALUES_PER_WORD;
}

template <typename TValue, typename TSpec, typename _TSize>
inline typename Size< String<TValue, Alloc<TSpec> > >::Type
reserve(
	String<TValue, Packed<TSpec> > & me,
	_TSize new_capacity,
	Limit)
{
SEQAN_CHECKPOINT
	typedef typename Size< String<TValue, Alloc<TSpec> > >::Type TSize;

	TSize me_capacity = capacity(me);
	if (me_capacity < (TSize)new_capacity) return me_capacity;
	return new_capacity;
}

template <typename TValue, typename TSpec, typename _TSize>
inline typename Size< String<TValue, Alloc<TSpec> > >::Type
reserve(
	String<TValue, Packed<TSpec> > & me,
	_TSize new_capacity,
	Insist)
{
SEQAN_CHECKPOINT
	typedef typename Size< String<TValue, Alloc<TSpec> > >::Type TSize;

	return new_capacity;
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Iteration
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THostspec, typename TSpec>
struct Iterator<String<TValue, Packed<THostspec> >, TSpec>
{
	typedef Iter<String<TValue, Packed<THostspec> >, Packed<THostspec> > Type;
};
template <typename TValue, typename THostspec, typename TSpec>
struct Iterator<String<TValue, Packed<THostspec> > const, TSpec>
{
	typedef Iter<String<TValue, Packed<THostspec> > const, Packed<THostspec> > Type;
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Iterator for packed strings
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec>
class Iter<TContainer, Packed<THostspec> >
{
private:
	typedef typename Host<TContainer>::Type THost;
	typedef typename Iterator<THost, Standard>::Type THostIterator;
	typedef typename Position<TContainer>::Type TPosition;

	typename _Pointer<TContainer>::Type data_container;
	THostIterator data_iterator;
	unsigned char data_bitpos;

//____________________________________________________________________________

public:
	Iter()
	{
SEQAN_CHECKPOINT
	}
	Iter(typename _Parameter<TContainer>::Type container_):
		data_container(_toPointer(container_)),
		data_iterator(begin(host(container_))),
		data_bitpos(0)
	{
SEQAN_CHECKPOINT
	}
	Iter(typename _Parameter<TContainer>::Type container_, TPosition pos_):
		data_container(_toPointer(container_))
	{
SEQAN_CHECKPOINT
		setPosition(*this, pos_);
	}
	Iter(Iter const & other_):
		data_container(other_.data_container),
		data_iterator(other_.data_iterator),
		data_bitpos(other_.data_bitpos)
	{
SEQAN_CHECKPOINT
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
		data_iterator = other_.data_iterator;
		data_bitpos = other_.data_bitpos;
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

	friend inline THostIterator &
	hostIterator(Iter & me)
	{
SEQAN_CHECKPOINT
		return me.data_iterator;
	}
	friend inline THostIterator const &
	hostIterator(Iter const & me)
	{
SEQAN_CHECKPOINT
		return me.data_iterator;
	}

//____________________________________________________________________________

	friend inline unsigned char &
	_bitpos(Iter & me)
	{
SEQAN_CHECKPOINT
		return me.data_bitpos;
	}
	friend inline unsigned char
	_bitpos(Iter const & me)
	{
SEQAN_CHECKPOINT
		return me.data_bitpos;
	}

//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// position
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec>
inline typename Position<Iter<TContainer, Packed<THostspec> > const>::Type
position(Iter<TContainer, Packed<THostspec> > const & me)
{
SEQAN_CHECKPOINT
	typedef typename Host<TContainer>::Type THost;
	THost const & host_ = host(container(me));
	return (hostIterator(me) - begin(host_)) * _PackedConsts<TContainer>::VALUES_PER_WORD + _bitpos(me) / _PackedConsts<TContainer>::BITS_PER_VALUE;
}

//////////////////////////////////////////////////////////////////////////////
// setPosition
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec, typename TPosition>
inline void
setPosition(Iter<TContainer, Packed<THostspec> > & me,
			TPosition pos_)
{
SEQAN_CHECKPOINT
	hostIterator(me) = begin(host(container(me))) + pos_ / _PackedConsts<TContainer>::VALUES_PER_WORD;
	_bitpos(me) = (pos_ % _PackedConsts<TContainer>::VALUES_PER_WORD) * _PackedConsts<TContainer>::BITS_PER_VALUE;
}

//////////////////////////////////////////////////////////////////////////////
// value
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec>
inline typename Reference<Iter<TContainer, Packed<THostspec> > >::Type
value(Iter<TContainer, Packed<THostspec> > & me)
{
SEQAN_CHECKPOINT
	return typename Reference<Iter<TContainer, Packed<THostspec> > >::Type(me);
}
template <typename TContainer, typename THostspec>
inline typename Reference<Iter<TContainer, Packed<THostspec> > const>::Type
value(Iter<TContainer, Packed<THostspec> > const & me)
{
SEQAN_CHECKPOINT
	return typename Reference<Iter<TContainer, Packed<THostspec> > const>::Type(me);
}

//////////////////////////////////////////////////////////////////////////////
// getValue
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec>
inline typename GetValue<Iter<TContainer, Packed<THostspec> > >::Type
getValue(Iter<TContainer, Packed<THostspec> > & me)
{
SEQAN_CHECKPOINT
	return (value(hostIterator(me)) >> _bitpos(me)) & _PackedConsts<TContainer>::VALUE_MASK;
}
template <typename TContainer, typename THostspec>
inline typename GetValue<Iter<TContainer, Packed<THostspec> > const>::Type
getValue(Iter<TContainer, Packed<THostspec> > const & me)
{
SEQAN_CHECKPOINT
	return (value(hostIterator(me)) >> _bitpos(me)) & _PackedConsts<TContainer>::VALUE_MASK;
}

/////////////////////////////////////////////////////////////////////////////
// assignValue
//////////////////////////////////////////////////////////////////////////////

template <typename TIter, typename TValue>
inline void
_assignValue_packed_string_iterator(TIter & me,
									TValue & _value)
{
	typedef typename Container<TIter>::Type TContainer;
	typedef typename Host<TContainer>::Type THost;
	typedef typename Value<THost>::Type THostValue;
	THostValue mask_ = _PackedConsts<TContainer>::VALUE_MASK << _bitpos(me);
	THostValue val_ = _value;
	val_ <<= _bitpos(me);

	assignValue(hostIterator(me), (getValue(hostIterator(me)) & ~(mask_)) | val_);
}


template <typename TContainer, typename THostspec, typename TValue>
inline void
assignValue(Iter<TContainer, Packed<THostspec> > & me,
			TValue const & _value)
{
SEQAN_CHECKPOINT
	typedef Iter<TContainer, Packed<THostspec> > TIterator;
	typename Value<TIterator>::Type _temp_value = _value; //conversion
	_assignValue_packed_string_iterator(me, _temp_value);
}
template <typename TContainer, typename THostspec, typename TValue>
inline void
assignValue(Iter<TContainer, Packed<THostspec> > const & me,
			TValue const & _value)
{
SEQAN_CHECKPOINT
	typedef Iter<TContainer, Packed<THostspec> > const TIterator;
	typename Value<TIterator>::Type _temp_value = _value; //conversion
	_assignValue_packed_string_iterator(me, _temp_value);
}

/////////////////////////////////////////////////////////////////////////////
// moveValue
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec, typename TValue>
inline void
moveValue(Iter<TContainer, Packed<THostspec> > & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	assignValue(me, _value);
}
template <typename TContainer, typename THostspec, typename TValue>
inline void
moveValue(Iter<TContainer, Packed<THostspec> > const & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	assignValue(me, _value);
}

//////////////////////////////////////////////////////////////////////////////
// valueConstruct
//////////////////////////////////////////////////////////////////////////////
//emulate construction and destruction

template <typename TContainer, typename THostspec>
inline void
valueConstruct(Iter<TContainer, Packed<THostspec> > const & /*it*/)
{
}
template <typename TContainer, typename THostspec, typename TParam>
inline void
valueConstruct(Iter<TContainer, Packed<THostspec> > const & it,
			   TParam const & param_)
{
	assignValue(it, param_);
}
template <typename TContainer, typename THostspec, typename TParam>
inline void
valueConstruct(Iter<TContainer, Packed<THostspec> > const & it,
			   TParam const & param_,
			   Move tag)
{
	moveValue(it, param_);
}

//////////////////////////////////////////////////////////////////////////////
// valueDestruct
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec>
inline void
valueDestruct(Iter<TContainer, Packed<THostspec> > const & /*it*/)
{
}

//////////////////////////////////////////////////////////////////////////////
// operator ==
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec>
inline bool
operator == (Iter<TContainer, Packed<THostspec> > const & left,
			 Iter<TContainer, Packed<THostspec> > const & right)
{
SEQAN_CHECKPOINT
	return (hostIterator(left) == hostIterator(right)) && (_bitpos(left) == _bitpos(right));
}

//////////////////////////////////////////////////////////////////////////////
// operator !=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec>
inline bool
operator != (Iter<TContainer, Packed<THostspec> > const & left,
			 Iter<TContainer, Packed<THostspec> > const & right)
{
SEQAN_CHECKPOINT
	return (hostIterator(left) != hostIterator(right)) || (_bitpos(left) != _bitpos(right));
}

//////////////////////////////////////////////////////////////////////////////
// operator >
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec>
inline bool
operator > (Iter<TContainer, Packed<THostspec> > const & left,
			Iter<TContainer, Packed<THostspec> > const & right)
{
SEQAN_CHECKPOINT
	return (hostIterator(left) > hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (_bitpos(left) > _bitpos(right)));
}

//////////////////////////////////////////////////////////////////////////////
// operator >=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec>
inline bool
operator >= (Iter<TContainer, Packed<THostspec> > const & left,
			 Iter<TContainer, Packed<THostspec> > const & right)
{
SEQAN_CHECKPOINT
	return (hostIterator(left) > hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (_bitpos(left) >= _bitpos(right)));
}

//////////////////////////////////////////////////////////////////////////////
// operator <
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec>
inline bool
operator < (Iter<TContainer, Packed<THostspec> > const & left,
			Iter<TContainer, Packed<THostspec> > const & right)
{
SEQAN_CHECKPOINT
	return (hostIterator(left) < hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (_bitpos(left) < _bitpos(right)));
}

//////////////////////////////////////////////////////////////////////////////
// operator <=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec>
inline bool
operator <= (Iter<TContainer, Packed<THostspec> > const & left,
			 Iter<TContainer, Packed<THostspec> > const & right)
{
SEQAN_CHECKPOINT
	return (hostIterator(left) < hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (_bitpos(left) <= _bitpos(right)));
}

//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec>
inline void
goNext(Iter<TContainer, Packed<THostspec> > & me)
{
SEQAN_CHECKPOINT
	int new_bitpos = _bitpos(me) + _PackedConsts<TContainer>::BITS_PER_VALUE;
	if (new_bitpos <= _PackedConsts<TContainer>::MAX_BIT_POS)
	{
		_bitpos(me) = (unsigned char) new_bitpos;
	}
	else
	{
		_bitpos(me) = 0;
		goNext(hostIterator(me));
	}
}

//////////////////////////////////////////////////////////////////////////////
// goPrevious
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec>
inline void
goPrevious(Iter<TContainer, Packed<THostspec> > & me)
{
SEQAN_CHECKPOINT
	int new_bitpos = _bitpos(me) - _PackedConsts<TContainer>::BITS_PER_VALUE;
	if (new_bitpos >= 0)
	{
		_bitpos(me) = (unsigned char) new_bitpos;
	}
	else
	{
		_bitpos(me) = _PackedConsts<TContainer>::MAX_BIT_POS;
		goPrevious(hostIterator(me));
	}

	goPrevious(hostIterator(me));
}

//////////////////////////////////////////////////////////////////////////////
// operator +
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec, typename TIntegral>
inline Iter<TContainer, Packed<THostspec> >
operator + (Iter<TContainer, Packed<THostspec> > const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, Packed<THostspec> >(container(left), position(left) + right);
}
template <typename TContainer, typename THostspec, typename TIntegral>
inline Iter<TContainer, Packed<THostspec> >
operator + (TIntegral left,
			Iter<TContainer, Packed<THostspec> > const & right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, Packed<THostspec> >(container(right), position(right) + left);
}

//////////////////////////////////////////////////////////////////////////////
// operator +=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec, typename TIntegral>
inline Iter<TContainer, Packed<THostspec> > &
operator += (Iter<TContainer, Packed<THostspec> > & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	setPosition(left, position(left) + right);
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// operator -
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec, typename TIntegral>
inline Iter<TContainer, Packed<THostspec> >
operator - (Iter<TContainer, Packed<THostspec> > const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, AdaptorIterator<THostspec> >(container(left), position(left) - right);
}

//____________________________________________________________________________

template <typename TContainer, typename THostspec>
inline typename Difference<Iter<TContainer, Packed<THostspec> > >::Type
operator - (Iter<TContainer, Packed<THostspec> > const & left,
			Iter<TContainer, Packed<THostspec> > const & right)
{
SEQAN_CHECKPOINT
	return position(left) - position(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator -=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename THostspec, typename TIntegral>
inline Iter<TContainer, Packed<THostspec> > &
operator -= (Iter<TContainer, Packed<THostspec> > & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	setPosition(left, position(left) - right);
	return left;
}

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
