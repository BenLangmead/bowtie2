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
  $Id: basic_definition.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_DEFINITION_H
#define SEQAN_HEADER_BASIC_DEFINITION_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Tag
{
};

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Default:
..summary:Tag that specifies default behavior.
..tag.Default:Use default behavior. 
*/
struct Default_;
typedef Tag<Default_> const Default;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Move Switch:
..summary:Switch to force move.
..tag.Move:Move instead of assign. 
..remarks.text:The difference between move constructor and copy constructor
is that the source object is not copied but moved into the target object.
The source object can lose its content and will be empty after
this operation in this case.
A move constructor can sigificantly faster than a copy constructor.
..example.code:String source("hello");
String target(source, Move()); // source is moved to target
std::cout << source; //nothing printed since source lost content
std::cout << target; //"hello"
..see:Function.move
*/

struct Move_;
typedef Tag<Move_> const Move;

//////////////////////////////////////////////////////////////////////////////

//Pass to c'tor of iterator to move it to the end
struct GoEnd_;
typedef Tag<GoEnd_> const GoEnd;


//////////////////////////////////////////////////////////////////////////////

//construct without initializing
struct MinimalCtor_;
typedef Tag<MinimalCtor_> const MinimalCtor;

//construct with initializing
struct NonMinimalCtor_;
typedef Tag<NonMinimalCtor_> const NonMinimalCtor;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Logical Values:
..summary:Tag that represents true and false.
..tag.True:The logical value "true".
..tag.False:The logical value "false".
*/
struct True { enum { VALUE = true }; };
struct False { enum { VALUE = false }; };


//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Nothing:
..summary:Tag that represents an absent parameter or an absent type.
..tag.Nothing:Omit parameter.
*/
///Empty Data Class.
struct Nothing {};



//////////////////////////////////////////////////////////////////////////////
// returns TTo const, if TFrom is const, TTo otherwise

template <typename TFrom, typename TTo>
struct _CopyConst
{
	typedef TTo Type;
};
template <typename TFrom, typename TTo>
struct _CopyConst<TFrom const, TTo>
{
	typedef TTo const Type;
};

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._RemoveConst:
..signature:_RemoveConst<T>
..returns:$t$ if $T$ is $t const$, otherwise $T$.
*/
template <typename T>
struct _RemoveConst
{
	typedef T Type;
};
template <typename T>
struct _RemoveConst<T const>:
	public _RemoveConst<T> {};

template <typename T>
struct _RemoveConst<T &>
{
	typedef typename _RemoveConst<T>::Type & Type;
};
template <typename T>
struct _RemoveConst<T *>
{
	typedef typename _RemoveConst<T>::Type * Type;
};
template <typename T, size_t I>
struct _RemoveConst<T const [I]>
{
	typedef T * Type;
};

//////////////////////////////////////////////////////////////////////////////
/**
.Internal._MakeUnsigned:
..signature:_MakeUnsigned<T>
..returns:$unsigned t$ if $T$ is not $unsigned t$, otherwise $T$.
*/
template <typename T>
struct _MakeUnsigned
{
	typedef T Type;
};

template <typename T>
struct _MakeUnsigned<T const> {
	typedef typename _MakeUnsigned<T>::Type const Type;
};

template <>
struct _MakeUnsigned<char>
{
	typedef unsigned char Type;
};

template <>
struct _MakeUnsigned<signed char>
{
	typedef unsigned char Type;
};

template <>
struct _MakeUnsigned<int>
{
	typedef unsigned int Type;
};

template <>
struct _MakeUnsigned<short>
{
	typedef unsigned short Type;
};

template <>
struct _MakeUnsigned<long>
{
	typedef unsigned long Type;
};

/*
template <>
struct _MakeUnsigned<long long>
{
	typedef unsigned long long Type;
};
*/

//////////////////////////////////////////////////////////////////////////////
/**
.Internal._MakeSigned:
..signature:_MakeSigned<T>
..returns:$signed t$ if $T$ is not $signed t$, otherwise $T$.
*/
template <typename T>
struct _MakeSigned
{
	typedef T Type;
};

template <typename T>
struct _MakeSigned<T const> {
	typedef typename _MakeSigned<T>::Type const Type;
};

template <>
struct _MakeSigned<char>
{
	typedef signed char Type;
};

template <>
struct _MakeSigned<unsigned char>
{
	typedef signed char Type;
};

template <>
struct _MakeSigned<unsigned int>
{
	typedef signed int Type;
};

template <>
struct _MakeSigned<unsigned short>
{
	typedef signed short Type;
};

template <>
struct _MakeSigned<unsigned long>
{
	typedef signed long Type;
};

/*
template <>
struct _MakeSigned<unsigned long long>
{
	typedef signed long long Type;
};
*/

//////////////////////////////////////////////////////////////////////////////
/**
.Internal._ClassIdentifier:
..signature:void * _ClassIdentifier<T>::getID()
..returns:A void * that identifies $T$.
...text:The returned values of two calls of $getID$ are equal if and only if
the used type $T$ was the same.
*/
template <typename T>
struct _ClassIdentifier
{
	static inline void *
	getID()
	{
SEQAN_CHECKPOINT
		static bool _id_dummy;
		return &_id_dummy;
	}
};

//////////////////////////////////////////////////////////////////////////////
/**
.Function.log2:
..cat:Miscellaneous
..summary:Computes logarithm of base 2 for integer types
..signature:unsigned int log2(i)
..param.i:An integer type.
..returns:The largest integer smaller or equal than
the logarithm of $i$.
*/

#if 0
template <int BITS_MAX>
struct _Log2_Impl
{
	template <typename T>
	static inline unsigned int
	log2(T val, unsigned int offset)
	{
		unsigned int val2 = val >> (BITS_MAX / 2);
		if (val2)
		{
			val = val2;
			offset += BITS_MAX / 2;
		}
		return _Log2_Impl<BITS_MAX / 2>::log2(val, offset);
	}
};

template <>
struct _Log2_Impl<1>
{
	template <typename T>
	static inline unsigned int
	log2(T /*val*/, unsigned int offset)
	{
		return offset;
	}
};


template <typename T>
inline unsigned int
log2(T val)
{
	enum
	{
//		BITS_PER_VALUE = BitsPerValue<T>::VALUE //TODO???
		BITS_PER_VALUE = sizeof(T) * 8
	};

	return _Log2_Impl<BITS_PER_VALUE>::log2(val, 0);
}
#endif

template <typename TValue, typename TExponent>
inline TValue _intPow(TValue a, TExponent b)
{
SEQAN_CHECKPOINT
	TValue ret = 1;
	while (b != 0)
	{
		if (b & 1) ret *= a;
		a *= a;
		b >>= 1;
	}	
	return ret;
}

//////////////////////////////////////////////////////////////////////////////
// to avoid conflicts with non-standard macros and namespaces
// we define our own Min/Max functions

template<typename _Tx> inline
const _Tx& _min(const _Tx& _Left, const _Tx& _Right)
{	// return smaller of _Left and _Right
	if (_Left < _Right)
		return _Left;
	else
		return _Right;
}

template<typename _Tx, typename _Ty> inline
_Tx _min(const _Tx& _Left, const _Ty& _Right)
{	// return smaller of _Left and _Right
    return (_Right < _Left ? _Right : _Left);
}

template<typename _Ty> inline
const _Ty& _max(const _Ty& _Left, const _Ty& _Right)
{	// return larger of _Left and _Right
	if (_Left < _Right)
		return _Right;
	else
		return _Left;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T1, typename T2>
struct _IsSameType
{
	enum {VALUE = false};
	typedef False Type;
};

template <typename T>
struct _IsSameType<T, T>
{
	enum {VALUE = true};
	typedef True Type;
};

template <typename T1, typename T2>
inline bool 
_isSameType()
{
	return _IsSameType<T1, T2>::VALUE;
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...


