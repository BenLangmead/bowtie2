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
  $Id: string_pointer.h,v 1.3 2009/10/12 16:00:57 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEQUENCE_POINTER_H
#define SEQAN_HEADER_SEQUENCE_POINTER_H

#include <cstring>

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////

/**
.Adaption.char array:
..summary:Zero terminated $char[]$ or $wchar_t[]$.
..remarks:Char arrays only support the Insist @Tag.Overflow Strategy.overflow strategy@.
*/


//////////////////////////////////////////////////////////////////////////////

/**
.Adaption.char array.remarks:The default overflow strategy
(both @Metafunction.DefaultOverflowImplicit@ and @Metafunction.DefaultOverflowExplicit@)
for all operations on char arrays is @Tag.Overflow Strategy.insist@.
*/

template <typename TValue>
struct DefaultOverflowImplicit;

template <typename TValue>
struct DefaultOverflowImplicit< TValue * >
{
	typedef Insist Type;
};
template <typename TValue>
struct DefaultOverflowImplicit< TValue * const>
{
	typedef Insist Type;
};
template <typename TValue, size_t SIZE>
struct DefaultOverflowImplicit< TValue [SIZE] >
{
	typedef Insist Type;
};
template <typename TValue, size_t SIZE>
struct DefaultOverflowImplicit< TValue const [SIZE] >
{
	typedef Insist Type;
};
//____________________________________________________________________________

template <typename TValue>
struct DefaultOverflowExplicit;

template <typename TValue>
struct DefaultOverflowExplicit< TValue * >
{
	typedef Insist Type;
};
template <typename TValue>
struct DefaultOverflowExplicit< TValue * const>
{
	typedef Insist Type;
};
template <typename TValue, size_t SIZE>
struct DefaultOverflowExplicit< TValue [SIZE] >
{
	typedef Insist Type;
};
template <typename TValue, size_t SIZE>
struct DefaultOverflowExplicit< TValue const [SIZE] >
{
	typedef Insist Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IsContiguous.param.T.type:Adaption.char array

template <typename TValue>
struct IsContiguous;

template <typename TValue>
struct IsContiguous< TValue * >
{
    typedef True Type;
	enum { VALUE = true };
};
template <typename TValue, size_t SIZE>
struct IsContiguous< TValue [SIZE] >
{
    typedef True Type;
	enum { VALUE = true };
};
template <typename TValue, size_t SIZE>
struct IsContiguous< TValue const [SIZE] >
{
    typedef True Type;
	enum { VALUE = true };
};

/*DISABLED
.Metafunction.IsString.param.T.type:Adaption.char array
*/

template <typename TValue>
struct IsSequence< TValue * >
{
    typedef True Type;
	enum { VALUE = true };
};
template <typename TValue, size_t SIZE>
struct IsSequence< TValue [SIZE] >
{
    typedef True Type;
	enum { VALUE = true };
};
template <typename TValue, size_t SIZE>
struct IsSequence< TValue const [SIZE] >
{
    typedef True Type;
	enum { VALUE = true };
};

//////////////////////////////////////////////////////////////////////////////


template <typename T>
inline typename Iterator<T *, typename DefaultGetIteratorSpec<T>::Type>::Type
begin(T * me)
{
SEQAN_CHECKPOINT
	return begin(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}

///.Function.begin.param.object.type:Adaption.char array

template <typename TValue>
inline typename Iterator<TValue *, Standard>::Type
begin(TValue * me,
	  Standard)
{
SEQAN_CHECKPOINT
	return me;
}

//folgende Versionen wurde wegen seltsamer Phaenomene bei VC++ 2003 hinzugenommen
template <typename TValue>
inline typename Iterator<TValue const *, Standard>::Type
begin(TValue const * me,
	  Standard)
{
SEQAN_CHECKPOINT
	return me;
}

template <typename TValue, typename TSpec>
inline typename Iterator<TValue *, Tag<TSpec> const>::Type
begin(TValue * me,
	  Tag<TSpec> const)
{
SEQAN_CHECKPOINT
	typedef typename Iterator<TValue *, Tag<TSpec> const>::Type TIterator;
	return TIterator(me, begin(me, Standard()));
}

template <typename TValue, typename TSpec>
inline typename Iterator<TValue const *, Tag<TSpec> const>::Type
begin(TValue const * me,
	  Tag<TSpec> const)
{
SEQAN_CHECKPOINT
	typedef typename Iterator<TValue const *, Tag<TSpec> const>::Type TIterator;
	return TIterator(me, begin(me, Standard()));
}

//////////////////////////////////////////////////////////////////////////////

///.Function.end.param.object.type:Adaption.char array

template <typename TValue>
inline typename Iterator<TValue *, Standard>::Type
end(TValue * me,
	Standard)
{
SEQAN_CHECKPOINT
	return begin(me, Standard()) + length(me);
}

//folgende Version wurde wegen eines seltsamen Phaenomens bei VC++ hinzugenommen
template <typename TValue>
inline typename Iterator<TValue const *, Standard>::Type
end(TValue const * me,
	Standard)
{
SEQAN_CHECKPOINT
	return begin(me, Standard()) + length(me);
}

template <typename TValue, typename TSpec>
inline typename Iterator<TValue *, Tag<TSpec> const>::Type
end(TValue * me,
	  Tag<TSpec> const tag_)
{
SEQAN_CHECKPOINT
	return begin(me, tag_) + length(me);
}

template <typename TValue, typename TSpec>
inline typename Iterator<TValue const *, Tag<TSpec> const>::Type
end(TValue const * me,
	  Tag<TSpec> const tag_)
{
SEQAN_CHECKPOINT
	return begin(me, tag_) + length(me);
}


//////////////////////////////////////////////////////////////////////////////

///.Function.value.param.container.type:Adaption.char array

template <typename TValue, typename TPos>
inline TValue &
value(TValue * me,
	  TPos pos)
{
SEQAN_CHECKPOINT
	return me[pos];
}

template <typename TValue, typename TPos>
inline TValue const &
value(TValue const * me,
	  TPos pos)
{
SEQAN_CHECKPOINT
	return me[pos];
}

//////////////////////////////////////////////////////////////////////////////
// assignValue
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TPos>
inline void
assignValue(TValue * me,
			TPos pos,
			TValue const & _value)
{
SEQAN_CHECKPOINT
	assign(value(me, pos), _value);
}

//////////////////////////////////////////////////////////////////////////////
// moveValue
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TPos>
inline void
moveValue(TValue * me,
		  TPos pos,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	move(value(me, pos), _value);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline bool
atEnd(TValue * pos)
{
SEQAN_CHECKPOINT
	return *pos == 0;
}

//____________________________________________________________________________

template <typename TValue>
inline bool
atEnd(TValue * pos,
	  TValue const * container)
{
SEQAN_CHECKPOINT
	return *pos == 0;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.length.param.object.type:Adaption.char array

template <typename TValue>
inline size_t
length(TValue * me)
{
SEQAN_CHECKPOINT
	TValue * it = me;
	TValue zero = TValue();
	while ( *it != zero) ++it;
	return it - me;
}

template <typename TValue>
inline size_t
length(TValue const * me)
{
SEQAN_CHECKPOINT
	TValue const * it = me;
	TValue const zero = TValue();
	while ( *it != zero) ++it;
	return it - me;
}

inline size_t
length(char * me)
{
SEQAN_CHECKPOINT
	return strlen(me);
}

inline size_t
length(char const * me)
{
SEQAN_CHECKPOINT
	return strlen(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void
_setLength(TValue * me,
		   size_t new_length)
{
SEQAN_CHECKPOINT
	me[new_length] = 0;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.clear.param.object.type:Adaption.char array

template <typename TValue>
inline void
clear(TValue * me)
{
SEQAN_CHECKPOINT
	//arrayDestruct(begin(me), length(me)); //??? Die Laengenbestimmung ist meistens nutzlos, braucht man sowieso nur fuer non-pod
	_setLength(me, 0);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.empty.param.object.type:Adaption.char array

template <typename TValue>
inline bool
empty(TValue * me)
{
SEQAN_CHECKPOINT
	return !me || (*me == TValue());
}

//////////////////////////////////////////////////////////////////////////////


template<typename TValue, typename TExpand>
inline size_t
_clearSpace(TValue * me,
		   size_t size,
		   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	return _ClearSpace_String_Base_<Tag<TExpand> const>::_clearSpace_(me, size);
}

template<typename TValue, typename TExpand>
inline size_t
_clearSpace(TValue * me,
		   size_t size,
		   size_t limit,
		   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	return _ClearSpace_String_Base_<Tag<TExpand> const>::_clearSpace_(me, size, limit);
}

template<typename TValue, typename TPosition, typename TExpand>
inline size_t
_clearSpace(TValue * me,
		   size_t size,
		   TPosition pos_begin,
		   TPosition pos_end,
		   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	return _ClearSpace_String_Base_<Tag<TExpand> const>::_clearSpace_(me, size, pos_begin, pos_end);
}

template<typename TValue, typename TPosition, typename TExpand>
inline size_t
_clearSpace(TValue * me,
		   size_t size,
		   TPosition pos_begin,
		   TPosition pos_end,
		   size_t limit,
		   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	return _ClearSpace_String_Base_<Tag<TExpand> const>::_clearSpace_(me, size, pos_begin, pos_end, limit);
}

//////////////////////////////////////////////////////////////////////////////
// assign
//////////////////////////////////////////////////////////////////////////////

///.Function.assign.param.target.type:.Adaption.char array
///.Function.assign.param.source.type:.Adaption.char array

//overload of binary version for strings:

template<typename TTargetValue, typename TSource>
inline void
assign(TTargetValue * target,
	   TSource & source)
{
SEQAN_CHECKPOINT
	typedef TTargetValue * TTarget;
	assign(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTargetValue, typename TSource>
inline void
assign(TTargetValue * target,
	   TSource const & source)
{
SEQAN_CHECKPOINT
	typedef TTargetValue * TTarget;
	assign(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

//____________________________________________________________________________

template<typename TTargetValue, typename TSource, typename TExpand>
inline void
assign(TTargetValue * target,
	   TSource const & source,
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Assign_String<Tag<TExpand> const>::assign_(target, source);
}

template<typename TTargetValue, typename TSource, typename TExpand>
inline void
assign(TTargetValue * target,
	   TSource const & source,
	   size_t limit,
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Assign_String<Tag<TExpand> const>::assign_(target, source, limit);
}

//____________________________________________________________________________
//this variant is a workaround for the "const array"-bug of VC++

template<typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
assign(TTargetValue * target,
	   TSourceValue const * source,
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Assign_String<Tag<TExpand> const>::assign_(target, source);
}

template<typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
assign(TTargetValue * target,
	   TSourceValue const * source,
	   size_t limit,
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Assign_String<Tag<TExpand> const>::assign_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
// move
//////////////////////////////////////////////////////////////////////////////

//overload of binary version for strings:

template<typename TTargetValue, typename TSource>
inline void
move(TTargetValue * & target,
	 TSource & source)
{
SEQAN_CHECKPOINT
	target = source;
}
template<typename TTargetValue, typename TSource>
inline void
move(TTargetValue * & target,
	 TSource const & source)
{
SEQAN_CHECKPOINT
	target = source;
}


//////////////////////////////////////////////////////////////////////////////
// append
//////////////////////////////////////////////////////////////////////////////

///.Function.append.param.target.type:.Adaption.char array
///.Function.append.param.source.type:.Adaption.char array

template<typename TTargetValue, typename TSource, typename TExpand>
inline void
append(TTargetValue * target,
	   TSource const & source,
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_String<Tag<TExpand> const>::append_(target, source);
}

template<typename TTargetValue, typename TSource, typename TExpand>
inline void
append(TTargetValue * target,
	   TSource const & source,
	   size_t limit,
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_String<Tag<TExpand> const>::append_(target, source, limit);
}

//____________________________________________________________________________
//this variant is a workaround for the "const array"-bug of VC++

template<typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
append(TTargetValue * target,
	   TSourceValue const * source,
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_String<Tag<TExpand> const>::append_(target, source);
}

template<typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
append(TTargetValue * target,
	   TSourceValue const * source,
	   size_t limit,
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_String<Tag<TExpand> const>::append_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
// replace
//////////////////////////////////////////////////////////////////////////////

///.Function.replace.param.target.type:.Adaption.char array
///.Function.replace.param.source.type:.Adaption.char array

template<typename TTargetValue, typename TSource, typename TExpand>
inline void
replace(TTargetValue * target,
		size_t pos_begin,
		size_t pos_end,
		TSource const & source,
		Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Replace_String<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}

template<typename TTargetValue, typename TSource, typename TExpand>
inline void
replace(TTargetValue * target,
		size_t pos_begin,
		size_t pos_end,
		TSource const & source,
		size_t limit,
		Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Replace_String<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}
//____________________________________________________________________________
//this variant is a workaround for the "const array"-bug of VC++

template<typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
replace(TTargetValue * target,
		size_t pos_begin,
		size_t pos_end,
		TSourceValue const * source,
		Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Replace_String<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}

template<typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
replace(TTargetValue * target,
		size_t pos_begin,
		size_t pos_end,
		TSourceValue const * source,
		size_t limit,
		Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Replace_String<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
// handling of iterators as begin and and
/*
template<typename TTargetValue, typename TSource, typename TExpand>
inline void
replace(TTargetValue * target,
		typename Iterator<TTargetValue *, Rooted>::Type pos_begin,
		typename Iterator<TTargetValue *, Rooted>::Type pos_end,
		TSource const & source,
		Tag<TExpand> const tag)
{
	replace(target, position(pos_begin), position(pos_end), source, tag);
}
template<typename TTargetValue, typename TSource, typename TExpand>
inline void
replace(TTargetValue * target,
		typename Iterator<TTargetValue *, Rooted>::Type pos_begin,
		typename Iterator<TTargetValue *, Rooted>::Type pos_end,
		TSource const & source,
		size_t limit,
		Tag<TExpand> const tag)
{
	replace(target, position(pos_begin), position(pos_end), source, limit, tag);
}
*/
//////////////////////////////////////////////////////////////////////////////
///.Function.resize.param.object.type:Adaption.char array

template <typename TValue, typename TExpand>
inline size_t
resize(
	TValue * me,
	size_t new_length,
	Tag<TExpand> const &)
{
SEQAN_CHECKPOINT
	return _Resize_String<Tag<TExpand> const>::resize_(me, new_length);
}

//////////////////////////////////////////////////////////////////////////////
///.Function.fill.param.object.type:Adaption.char array

template <typename TValue, typename TExpand>
inline size_t
fill(
	TValue * me,
	size_t new_length,
	TValue const & val,
	Tag<TExpand> const &)
{
SEQAN_CHECKPOINT
	return _Fill_String<Tag<TExpand> const>::fill_(me, new_length, val);
}

//////////////////////////////////////////////////////////////////////////////
//PROBLEM: ambiguitiy "pointer/iterator" and "c-style string"
//workaround: disable all operators
/*
template <typename TLeftValue, typename TRight >
TLeftValue const *
operator += (TLeftValue * left,
			 TRight const & right)
{
SEQAN_CHECKPOINT
	append(left, right);
	return left;
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TRight >
inline bool
isEqual(TLeftValue * left,
		TRight const & right)
{
SEQAN_CHECKPOINT
	typename Comparator<TLeftValue *>::Type _lex(left, right);
    return isEqual(_lex);
}
/*
template <typename TLeftValue, typename TRight >
inline bool
operator == (TLeftValue * left,
			TRight const & right)
{
SEQAN_CHECKPOINT
	typename Comparator<TLeftValue *>::Type _lex(left, right);
    return isEqual(_lex);
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TRight >
inline bool
isNotEqual(TLeftValue * left,
		   TRight const & right)
{
SEQAN_CHECKPOINT
	typename Comparator<TLeftValue *>::Type _lex(left, right);
    return isNotEqual(_lex);
}
/*
template <typename TLeftValue, typename TRight >
inline bool
operator != (TLeftValue * left,
			 TRight const & right)
{
SEQAN_CHECKPOINT
	typename Comparator<TLeftValue *>::Type _lex(left, right);
    return isNotEqual(_lex);
}
*/

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TRight>
inline bool
isLess(TLeftValue * left,
	   TRight const & right)
{
SEQAN_CHECKPOINT
	return isLess(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
/*
template <typename TLeftValue, typename TRight>
inline bool
operator < (TLeftValue * left,
			TRight const & right)
{
SEQAN_CHECKPOINT
	return isLess(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TRight>
inline bool
isLessOrEqual(TLeftValue * left,
			 TRight const & right)
{
SEQAN_CHECKPOINT
	return isLessOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
/*
template <typename TLeftValue, typename TRight>
inline bool
operator <= (TLeftValue * left,
			 TRight const & right)
{
SEQAN_CHECKPOINT
	return isLessOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TRight>
inline bool
isGreater(TLeftValue * left,
		TRight const & right)
{
SEQAN_CHECKPOINT
	return isGreater(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
/*
template <typename TLeftValue, typename TRight>
inline bool
operator > (TLeftValue * left,
		TRight const & right)
{
SEQAN_CHECKPOINT
	return isGreater(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TRight>
inline bool
isGreaterOrEqual(TLeftValue * left,
		TRight const & right)
{
SEQAN_CHECKPOINT
	return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
/*
template <typename TLeftValue, typename TRight>
inline bool
operator >= (TLeftValue * left,
		TRight const & right)
{
SEQAN_CHECKPOINT
	return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
*/
//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

//____________________________________________________________________________

#endif //#ifndef SEQAN_HEADER_...
