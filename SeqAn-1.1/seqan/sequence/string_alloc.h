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
  $Id: string_alloc.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEQUENCE_STRING_ALLOC_H
#define SEQAN_HEADER_SEQUENCE_STRING_ALLOC_H


namespace SEQAN_NAMESPACE_MAIN
{

/**
.Spec.Alloc String:
..cat:Strings
..general:Class.String
..summary:Expandable string that is stored on heap.
..signature:String<TValue, Alloc<TSpec> >
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.TSpec:The specializing type.
...default:$void$
*/
//////////////////////////////////////////////////////////////////////////////
//expandable string
//////////////////////////////////////////////////////////////////////////////
//Default: TSpec == void

template <typename TValue>
class String<TValue, Alloc<void> >
{
public:
	typename Value<String>::Type * data_begin;
	typename Value<String>::Type * data_end;
	size_t data_capacity;

//____________________________________________________________________________

public:
	String():
		data_begin(0),
		data_end(0),
		data_capacity(0)
	{
SEQAN_CHECKPOINT
	}

	template <typename TSource>
	String(TSource & source):
		data_begin(0),
		data_end(0),
		data_capacity(0)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
	}
	template <typename TSource>
	String(TSource const & source):
		data_begin(0),
		data_end(0),
		data_capacity(0)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
	}
	String(String const & source):
		data_begin(0),
		data_end(0),
		data_capacity(0)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
	}
	String(String const & source, Move):
		data_begin(0),
		data_end(0),
		data_capacity(0)
	{
SEQAN_CHECKPOINT
		move(*this, source);
	}
	template <typename TSource, typename TSize>
	String(TSource & source, TSize limit):
		data_begin(0),
		data_end(0),
		data_capacity(0)
	{
SEQAN_CHECKPOINT
		assign(*this, source, limit);
	}
	template <typename TSource, typename TSize>
	String(TSource const & source, TSize limit):
		data_begin(0),
		data_end(0),
		data_capacity(0)
	{
SEQAN_CHECKPOINT
		assign(*this, source, limit);
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
		arrayDestruct(this->data_begin, this->data_end);
		_deallocateStorage(*this, this->data_begin, data_capacity);
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

	friend inline typename Iterator<String, Standard>::Type
	begin(String & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return me.data_begin;
	}
	friend inline typename Iterator<String const, Standard>::Type
	begin(String const & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return me.data_begin;
	}

//____________________________________________________________________________

	friend inline typename Iterator<String, Standard>::Type
	end(String & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return me.data_end;
	}
	friend inline typename Iterator<String const, Standard>::Type
	end(String const & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return me.data_end;
	}

//____________________________________________________________________________

	friend inline size_t
	capacity(String & me) 
	{
SEQAN_CHECKPOINT
		return me.data_capacity;
	}

	friend inline size_t
	capacity(String const & me) 
	{
SEQAN_CHECKPOINT
		return me.data_capacity;
	}

//____________________________________________________________________________
/* Entwicklungsschrott?
	inline void 
	move(String & target, 
		 String & source)
	{
		clear(target);
		target.data_begin = source.data_begin;
		target.data_end = source.data_end;
		target.data_capacity = source.data_capacity;

		source.data_begin = 0;
		source.data_end = 0;
		source.data_capacity = 0;
	}
*/
//____________________________________________________________________________

/**
.Internal._setBegin:
*/
	friend inline void 
	_setBegin(
		String & me, 
		typename Value<String>::Type * new_begin)
	{
SEQAN_CHECKPOINT
		me.data_begin = new_begin;
	}

//____________________________________________________________________________

/**
.Internal._setLength:
..cat:Functions
..summary:Set the length of container.
..signature:_setLength(object, new_length)
..param.object:A container.
..param.object.type:Spec.Alloc String
..param.new_length:The new length.
*/
	friend inline void 
	_setLength(
		String & me, 
		size_t new_length)
	{
SEQAN_CHECKPOINT
		me.data_end = me.data_begin + new_length;
	}

//____________________________________________________________________________

/**
.Internal._setCapacity:
*/
	friend inline void 
	_setCapacity(
		String & me, 
		size_t new_capacity)
	{
SEQAN_CHECKPOINT
		me.data_capacity = new_capacity;
	}

//____________________________________________________________________________

/**
.Internal._allocateStorage:
..cat:Functions
..summary:Allocates a new buffer for a container.
..signature:_allocateStorage(object, new_capacity)
..param.object:A container.
..param.object.type:Spec.Alloc String
..param.new_capacity:The capacity of the new allocated buffer.
..returns:The old butter $object$, that is replaced by the new allocated buffer.
..remarks:The returned buffer must be deallocated by @Internal._deallocateStorage@.
..remarks:This function does not construct objects in the allocated buffer.
..see:Internal._reallocateStorage
*/
	friend inline typename Value<String>::Type * 
	_allocateStorage(
		String & me, 
		size_t new_capacity)
	{
SEQAN_CHECKPOINT
		size_t size = _computeSize4Capacity(me, new_capacity);
		typename Value<String>::Type * return_value = me.data_begin;
		allocate(me, me.data_begin, size, TagAllocateStorage());
		me.data_capacity = new_capacity;
		return return_value;
	}

	//____________________________________________________________________________

/**
.Internal._deallocateStorage:
..cat:Functions
..summary:Deallocates a buffer of a container.
..signature:_deallocateStorage(object, buffer, capacity)
..param.object:A container.
..param.object.type:Spec.Alloc String
..param.buffer:The buffer that will be deallocated.
..param.capacity:The capacity of $buffer$.
..remarks:All objects in the buffer must be destroyed before calling $_deallocateStorage$.
..see:Internal._allocateStorage
..see:Internal._reallocateStorage
*/
	friend inline void 
	_deallocateStorage(
		String & me, 
		typename Value<String>::Type * ptr, 
		size_t capacity)
	{
SEQAN_CHECKPOINT
		size_t size = _computeSize4Capacity(me, capacity);
		deallocate(me, ptr, size, TagAllocateStorage());
	}

//____________________________________________________________________________

};


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct DefaultOverflowImplicit<String<TValue, Alloc<TSpec> > >
{
	typedef Generous Type;
};

template <typename TValue, typename TSpec>
struct DefaultOverflowImplicit<String<TValue, Alloc<TSpec> > const >
{
	typedef Generous Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct IsContiguous< String<TValue, Alloc<TSpec> > >
{
    typedef True Type;
	enum { VALUE = true };
};


//////////////////////////////////////////////////////////////////////////////

template <typename TTargetValue, typename TSourceValue, typename TSpec>
inline void 
move(String<TTargetValue, Alloc<TSpec> > & target, 
	 String<TSourceValue, Alloc<TSpec> > & source)
{
	_moveContiguous(target, source);
}
template <typename TTargetValue, typename TSourceValue, typename TSpec>
inline void 
move(String<TTargetValue, Alloc<TSpec> > & target, 
	 String<TSourceValue, Alloc<TSpec> > const & source)
{
	_moveContiguous(target, source);
}

template <typename TValue, typename TSpec>
inline void 
move(String<TValue, Alloc<TSpec> > & target, 
	 String<TValue, Alloc<TSpec> > & source)
{
	clear(target);
	target.data_begin = source.data_begin;
	target.data_end = source.data_end;
	target.data_capacity = source.data_capacity;

	source.data_begin = 0;
	source.data_end = 0;
	source.data_capacity = 0;
}
template <typename TValue, typename TSpec>
inline void 
move(String<TValue, Alloc<TSpec> > & target, 
	 String<TValue, Alloc<TSpec> > const & source)
{
	clear(target);
	target.data_begin = source.data_begin;
	target.data_end = source.data_end;
	target.data_capacity = source.data_capacity;

	source.data_begin = 0;
	source.data_end = 0;
	source.data_capacity = 0;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.reserve.param.object.type:Spec.Alloc String

template <typename TValue, typename TSpec, typename _TSize, typename TExpand>
inline typename Size< String<TValue, Alloc<TSpec> > >::Type
reserve(
	String<TValue, Alloc<TSpec> > & seq, 
	_TSize new_capacity,
	Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	typedef typename Size< String<TValue, Alloc<TSpec> > >::Type TSize;

	TSize old_capacity = capacity(seq);
	if (old_capacity >= (TSize)new_capacity) return new_capacity;

	TSize seq_length = length(seq);
	typename Value< String<TValue, Alloc<TSpec> > >::Type * old_array = _reallocateStorage(seq, new_capacity, tag);
	if (old_array)
	{//buffer was replaced, destruct old buffer
		arrayConstructCopy(old_array, old_array + seq_length, begin(seq, Standard()));
		arrayDestruct(old_array, old_array + seq_length);
		_deallocateStorage(seq, old_array, old_capacity);
		_setLength(seq, seq_length);
	}
	else if (!old_capacity)
	{//new buffer created and the string had no buffer yet
		_setLength(seq, seq_length);
	}
	return new_capacity;
}

template <typename TValue, typename TSpec, typename _TSize>
inline typename Size< String<TValue, Alloc<TSpec> > >::Type
reserve(
	String<TValue, Alloc<TSpec> > & me, 
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
	String<TValue, Alloc<TSpec> > & /*me*/, 
	_TSize new_capacity,
	Insist)
{
SEQAN_CHECKPOINT
	typedef typename Size< String<TValue, Alloc<TSpec> > >::Type TSize;

	return new_capacity;
}

//////////////////////////////////////////////////////////////////////////////



} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
