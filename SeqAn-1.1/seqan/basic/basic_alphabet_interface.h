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
  $Id: basic_alphabet_interface.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_ALPHABET_INTERFACE_H
#define SEQAN_HEADER_BASIC_ALPHABET_INTERFACE_H

#include <new>

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//IsSimple
//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.IsSimple:
..summary:Tests type to be simple.
..signature:IsSimple<T>::Type
..param.T:Type that is tested.
..returns.param.Type:@Tag.Logical Values.True@, if $T$ is a simple type, @Tag.Logical Values.False@ otherwise.
...default:@Tag.Logical Values.False@
..remarks:A simple type is a type that does not need constructors to be created,
a destructor to be destroyed, and copy assignment operators or copy constructors
to be copied. All POD ("plain old data") types are simple, but some
non-POD types could be simple too, e.g. some specializations of @Class.SimpleType@.
..see:Class.SimpleType
*/

template <typename T>
struct _IsSimple {
	typedef False Type;
};

template <typename T>
struct IsSimple:
	public _IsSimple<T> {};
template <typename T>
struct IsSimple<T const>:
	public IsSimple<T> {};

//////////////////////////////////////////////////////////////////////////////
//very basic Alphabets

typedef char Ascii;
typedef unsigned char Byte;
typedef wchar_t Unicode;

//////////////////////////////////////////////////////////////////////////////

/**
.Function.valueConstruct:
..cat:Content Manipulation
..summary:Constructs an object at specified position.
..signature:valueConstruct(iterator [, param [, move_tag] ])
..param.iterator:Pointer or iterator to position where the object should be constructed.
..param.param:Parameter that is forwarded to constructor. (optional)
..param.move_tag:Instance of the @Tag.Move Switch.move switch tag@. (optional)
...remarks:If the @Tag.Move Switch.move switch tag@ is specified, it is forwarded to the constructor,
so the constructed object must support move construction.
..remarks:The type of the destructed object is the @Metafunction.Value.value type@ of $iterator$.
*/

struct _ValueConstructor 
{
	template <typename TIterator>
	static inline void
	construct(TIterator it)
	{
		typedef typename Value<TIterator>::Type TValue;
		new( & value(it) ) TValue;
	}

	template <typename TIterator, typename TParam>
	static inline void
	construct(TIterator it,
			  TParam const & param_)
	{
		typedef typename Value<TIterator>::Type TValue;
		new( & value(it) ) TValue(param_);
	}

	template <typename TIterator, typename TParam>
	static inline void
	construct(TIterator it,
			  TParam const & param_,
			  Move tag)
	{
		typedef typename Value<TIterator>::Type TValue;
		new( & value(it) ) TValue(param_, tag);
	}
};

struct _ValueConstructorProxy 
{
	template <typename TIterator>
	static inline void construct(TIterator) {}

	template <typename TIterator, typename TParam>
	static inline void construct(TIterator, TParam const &) {}

	template <typename TIterator, typename TParam>
	static inline void construct(TIterator, TParam const &, Move) {}
};

//____________________________________________________________________________

struct _ValueDestructor 
{
	template <typename TIterator>
	static inline void
	destruct(TIterator it)
	{
		typedef typename Value<TIterator>::Type TValue;
		value(it).~TValue();
	}
};
struct _ValueDestructorProxy 
{
	template <typename TIterator>
	static inline void destruct(TIterator) {}
};

//____________________________________________________________________________

template <typename TIterator>
inline void
valueConstruct(TIterator it)
{
SEQAN_CHECKPOINT
	typedef typename IF<
		TYPECMP<
			typename Value<TIterator>::Type &,
			typename Reference<TIterator>::Type
		>::VALUE,
	// THEN
		_ValueConstructor,			// true,  types are equal
	// ELSE
		_ValueConstructorProxy		// false, types differ -> value() returns a proxy
	>::Type TConstructor;

	TConstructor::construct(it);
}

template <typename TIterator, typename TParam>
inline void
valueConstruct(TIterator it,
			   TParam const & param_)
{
SEQAN_CHECKPOINT
	typedef typename IF<
		TYPECMP<
			typename Value<TIterator>::Type &,
			typename Reference<TIterator>::Type
		>::VALUE,
	// THEN
		_ValueConstructor,			// true,  types are equal
	// ELSE
		_ValueConstructorProxy		// false, types differ -> value() returns a proxy
	>::Type TConstructor;

	TConstructor::construct(it, param_);
}

template <typename TIterator, typename TParam>
inline void
valueConstruct(TIterator it,
			   TParam const & param_,
			   Move tag)
{
SEQAN_CHECKPOINT
	typedef typename IF<
		TYPECMP<
			typename Value<TIterator>::Type &,
			typename Reference<TIterator>::Type
		>::VALUE,
	// THEN
		_ValueConstructor,			// true,  types are equal
	// ELSE
		_ValueConstructorProxy		// false, types differ -> value() returns a proxy
	>::Type TConstructor;

	TConstructor::construct(it, param_, tag);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.valueDestruct:
..cat:Content Manipulation
..summary:Destoys an object at specified position.
..signature:valueDestruct(iterator)
..param.iterator:Pointer or iterator to position where the object should be constructed.
..remarks:The type of the constructed object is the @Metafunction.Value.value type@ of $iterator$.
..see:Function.valueConstruct
*/
template <typename TIterator>
inline void
valueDestruct(TIterator it)
{
SEQAN_CHECKPOINT
	typedef typename IF<
		TYPECMP<
			typename Value<TIterator>::Type &,
			typename Reference<TIterator>::Type
		>::VALUE,
	// THEN
		_ValueDestructor,			// true,  types are equal
	// ELSE
		_ValueDestructorProxy		// false, types differ -> value() returns a proxy
	>::Type TDestructor;

	TDestructor::destruct(it);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.valueConstructMove:
..cat:Content Manipulation
..summary:Move constructs an object at specified position.
..signature:valueConstructMove(iterator, param)
..param.iterator:Pointer or iterator to position where the object should be constructed.
..param.param:Parameter that is moved to the new constructed object.
..remarks:The type of the destructed object is the @Metafunction.Value.value type@ of $iterator$.
..remarks:The default implementation just calls @Function.valueConstruct@.
*/
template <typename TIterator, typename TValue>
inline void
valueConstructMove(TIterator it, TValue const & value)
{
	valueConstruct(it, value);
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//arrayConstruct
//////////////////////////////////////////////////////////////////////////////
/**
.Function.arrayConstruct:
..cat:Array Handling
..summary:Construct objects in a given memory buffer.
..signature:arrayConstruct(begin, end [, value])
..param.begin:Iterator to the begin of the range that is to be constructed.
..param.end:Iterator behind the end of the range.
..param.value:Argument that is forwarded to the constructor. (optional)
...text:An appropriate constructor is required. 
If $value$ is not specified, the default constructor is used. 
..remarks:The type of the constructed Objects is the @Metafunction.Value.value type@
of $begin$ and $end$.
..see:Function.arrayDestruct
..see:Function.arrayConstructCopy
..see:Function.arrayFill
..see:Class.SimpleType
..see:Function.valueConstruct
*/
template<typename TIterator1, typename TIterator2>
inline void 
_arrayConstruct_Default(TIterator1 begin_, 
						TIterator2 end_)
{
SEQAN_CHECKPOINT
	while (begin_ != end_)
	{
		valueConstruct(begin_);
		++begin_;
	}
}
template<typename TIterator1, typename TIterator2>
inline void 
arrayConstruct(TIterator1 begin_, 
			   TIterator2 end_)
{
SEQAN_CHECKPOINT
	_arrayConstruct_Default(begin_, end_);
}

//____________________________________________________________________________

template<typename TIterator1, typename TIterator2, typename TParam>
inline void 
_arrayConstruct_Default(TIterator1 begin_, 
						TIterator2 end_, 
						TParam const & param_)
{
SEQAN_CHECKPOINT
	while (begin_ != end_)
	{
		valueConstruct(begin_, param_);
		++begin_;
	}
}
template<typename TIterator1, typename TIterator2, typename TParam>
inline void 
arrayConstruct(TIterator1 begin_, 
			   TIterator2 end_, 
			   TParam const & param_)
{
SEQAN_CHECKPOINT
	_arrayConstruct_Default(begin_, end_, param_);
}

//////////////////////////////////////////////////////////////////////////////
//arrayConstructCopy
//////////////////////////////////////////////////////////////////////////////

/**
.Function.arrayConstructCopy:
..cat:Array Handling
..summary:Copy constructs an array of objects into in a given memory buffer.
..signature:arrayConstructCopy(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source range.
..param.source_end:Iterator behind the last element of the source range.
...text:$source_end$ should have the same type as $source_begin$.
..param.target:Pointer to the memory block the new objects will be constructed in.
...text:The type of $target$ specifies the type of the constructed objects:
If $T*$ is the type of $target$, then the function constructs objects of type $T$. 
...text:The memory buffer should be large enough to store $source_end$ - $source_begin$ objects.
An appropriate (copy-) constructor that constructs an target objects given a source object is required.
..see:Function.arrayDestruct
..see:Function.arrayCopyForward
..see:Function.arrayCopy
..see:Function.valueConstruct
*/
template<typename TTarget, typename TSource1, typename TSource2>
inline void 
_arrayConstructCopy_Default(TSource1 source_begin, 
							TSource2 source_end, 
							TTarget target_begin)
{
SEQAN_CHECKPOINT
	while (source_begin != source_end)
	{
		valueConstruct(target_begin, *source_begin);
		++source_begin;
		++target_begin;
	}
}

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayConstructCopy(TSource1 source_begin, 
				   TSource2 source_end, 
				   TTarget target_begin)
{
SEQAN_CHECKPOINT
	_arrayConstructCopy_Default(source_begin, source_end, target_begin);
}

//////////////////////////////////////////////////////////////////////////////
//arrayConstructMove
//////////////////////////////////////////////////////////////////////////////

/**
.Function.arrayConstructMove:
..cat:Array Handling
..summary:Move constructs an array of objects into in a given memory buffer.
..signature:arrayConstructMove(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source range.
..param.source_end:Iterator behind the last element of the source range.
...text:$source_end$ should have the same type as $source_begin$.
..param.target:Pointer to the memory block the new objects will be constructed in.
...text:The type of $target$ specifies the type of the constructed objects:
If $T*$ is the type of $target$, then the function constructs objects of type $T$. 
...text:The memory buffer should be large enough to store $source_end$ - $source_begin$ objects.
An appropriate move constructor that constructs an target objects given a source object is required.
..see:Function.arrayDestruct
..see:Function.arrayConstructCopy
..see:Function.arrayMoveForward
..see:Function.arrayMove
..see:Function.valueConstruct
*/
template<typename TTarget, typename TSource1, typename TSource2>
inline void 
_arrayConstructMove_Default(TSource1 source_begin, 
							TSource2 source_end, 
							TTarget target_begin)
{
SEQAN_CHECKPOINT
	while (source_begin < source_end)
	{
		valueConstructMove(target_begin, *source_begin);
		++source_begin;
		++target_begin;
	}
}

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayConstructMove(TSource1 source_begin, 
				   TSource2 source_end, 
				   TTarget target_begin)
{
SEQAN_CHECKPOINT
	_arrayMoveConstruct_Default(source_begin, source_end, target_begin);
}

//////////////////////////////////////////////////////////////////////////////
//arrayDestruct
//////////////////////////////////////////////////////////////////////////////

/**
.Function.arrayDestruct:
..cat:Array Handling
..summary:Destroys an array of objects.
..signature:arrayDestruct(begin, end)
..param.begin:Iterator to the begin of the range that is to be destructed.
..param.end:Iterator behind the end of the range.
..remarks:This function does not deallocates the memory.
..see:Class.SimpleType
..see:Function.valueDestruct
*/
template<typename TIterator1, typename TIterator2>
inline void 
_arrayDestruct_Default(TIterator1 begin_, 
					   TIterator2 end_)
{
SEQAN_CHECKPOINT
	while (begin_ != end_)
	{
		valueDestruct(begin_);
		++begin_;
	}
}
template<typename TIterator1, typename TIterator2>
inline void 
arrayDestruct(TIterator1 begin_, 
			  TIterator2 end_)
{
SEQAN_CHECKPOINT
	_arrayDestruct_Default(begin_, end_);
}

//////////////////////////////////////////////////////////////////////////////
//arrayFill
//////////////////////////////////////////////////////////////////////////////

/**
.Function.arrayFill:
..cat:Array Handling
..summary:Assigns one object to each element of a range.
..signature:arrayFill(begin, end, value)
..param.begin:Iterator to the begin of the range that is to be filled.
..param.end:Iterator behind the end of the range.
..param.value:Argument that is assigned to all $count$ objects in $array$.
..remarks:All objects $target_begin[0]$ to $target_begin[count-1]$ are set to $value$.
..see:Function.arrayCopy
..see:Function.arrayCopyForward
*/
template<typename TIterator1, typename TIterator2, typename TValue>
inline void 
arrayFill(TIterator1 begin_,
		  TIterator2 end_, 
		  TValue const & value)
{
SEQAN_CHECKPOINT
	::std::fill_n(begin_, end_ - begin_, value);
}

//////////////////////////////////////////////////////////////////////////////
//arrayCopyForward
//////////////////////////////////////////////////////////////////////////////

/**
.Function.arrayCopyForward:
..cat:Array Handling
..summary:Copies a range of objects into another range of objects starting from the first element.
..signature:arrayCopyForward(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source array.
..param.source_end:Iterator behind the last element of the source array.
...text:$source_end$ must have the same type as $source_begin$.
..param.target:Iterator to the first element of the target array.
...text:The target capacity should be at least as long as the source range.
..remarks.note:Be careful if source and target range overlap, because in this case
	some source elements could be accidently overwritten before they are moved.
..remarks:If there is no need for the source elements to persist, consider to use 
@Function.arrayMoveForward@ instead to improve performance.
..see:Class.SimpleType
*/
template<typename TTarget, typename TSource1, typename TSource2>
inline void 
_arrayCopyForward_Default(TSource1 source_begin, 
						  TSource2 source_end, 
						  TTarget target_begin)
{
SEQAN_CHECKPOINT
	::std::copy(source_begin, source_end, target_begin);
}
template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayCopyForward(TSource1 source_begin, 
				 TSource2 source_end, 
				 TTarget target_begin)
{
SEQAN_CHECKPOINT
	_arrayCopyForward_Default(source_begin, source_end, target_begin);	
}

//////////////////////////////////////////////////////////////////////////////
//arrayCopyBackward
//////////////////////////////////////////////////////////////////////////////

/**
.Function.arrayCopyBackward:
..cat:Array Handling
..summary:Copies a range of objects into another range of objects starting from the last element.
..signature:arrayCopyBackward(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source array.
..param.source_end:Iterator behind the last element of the source array.
...text:$source_end$ must have the same type as $source_begin$.
..param.target:Iterator to the first element of the target array.
...text:The target capacity should be at least as long as the source range.
..remarks.note:Be careful if source and target range overlap, because in this case
	some source elements could be accidently overwritten before they are moved.
..remarks.text:If source and target do not overlap, consider to use the function
@Function.arrayCopyForward@ instead that is faster in some cases.
..remarks:If there is no need for the source elements to persist, consider to use 
@Function.arrayMoveBackward@ instead to improve performance.
..remarks.note:The semantic of this function's argument $target$ differ from the arguments of $::std::copy_backward$.
..see:Function.arrayCopyForward
..see:Class.SimpleType
*/
template<typename TTarget, typename TSource1, typename TSource2>
inline void 
_arrayCopyBackward_Default(TSource1 source_begin, 
						   TSource2 source_end, 
						   TTarget target_begin)
{
SEQAN_CHECKPOINT
	::std::copy_backward(source_begin, source_end, target_begin + (source_end - source_begin));
}
template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayCopyBackward(TSource1 source_begin, 
				  TSource2 source_end, 
				  TTarget target_begin)
{
SEQAN_CHECKPOINT
	_arrayCopyBackward_Default(source_begin, source_end, target_begin);
}


//////////////////////////////////////////////////////////////////////////////
//arrayCopy
//////////////////////////////////////////////////////////////////////////////

/**
.Function.arrayCopy:
..cat:Array Handling
..summary:Copies a range of objects into another range of objects.
..signature:arrayCopy(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source range.
..param.source_end:Iterator behind the last element of the source range.
...text:$source_end$ must have the same type as $source_begin$.
..param.target:Iterator to the first element of the target range.
...text:The target capacity should be at least as long as the source range.
..remarks.text:If source and target range do not overlap, consider to use
	@Function.arrayCopyForward@ instead to improve performance.
..remarks:If there is no need for the source elements to persist, consider to use 
	@Function.arrayMoveForward@ instead to improve performance.
..DISABLED.remarks.note:Be careful if source and target range overlap and the size of the
	source elements differ from the size of target elements, because in this case
	some source elements could be accidently overwritten before they are moved.
..see:Function.arrayCopyForward
..see:Function.arrayCopyBackward
..see:Class.SimpleType
*/
template<typename TTarget, typename TSource1, typename TSource2>
inline void arrayCopy(TSource1 source_begin, 
					  TSource2 source_end, 
					  TTarget target_begin)
{
	if ((void *) source_begin >= (void *) target_begin)
	{
SEQAN_CHECKPOINT
		arrayCopyForward(source_begin, source_end, target_begin);
	}
	else
	{
SEQAN_CHECKPOINT
		arrayCopyBackward(source_begin, source_end, target_begin);
	}
}

//////////////////////////////////////////////////////////////////////////////
//arrayMoveForward
//////////////////////////////////////////////////////////////////////////////

/**
.Function.arrayMoveForward:
..cat:Array Handling
..summary:Moves a range of objects into another range of objects starting from the first element.
..signature:arrayMoveForward(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source array.
..param.source_end:Iterator behind the last element of the source array.
...text:$source_end$ must have the same type as $source_begin$.
..param.target:Iterator to the first element of the target array.
...text:The target capacity should be at least as long as the source range.
..remarks:The function possibly clears (but does not destroy) the source elements.
	If source elements must persist, consider to use @Function.arrayCopyForward@ instead.
..remarks.note:Be careful if source and target range overlap, because in this case
	some source elements could be accidently overwritten before they are moved.
..see:Function.arrayCopyForward
..see:Class.SimpleType
*/
template<typename TTarget, typename TSource1, typename TSource2>
inline void 
_arrayMoveForward_Default(TSource1 source_begin, 
						  TSource2 source_end, 
						  TTarget target_begin)
{
SEQAN_CHECKPOINT
	while (source_begin != source_end)
	{
		move(*target_begin, *source_begin);
		++source_begin;
		++target_begin;
	}
}
template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayMoveForward(TSource1 source_begin, 
				 TSource2 source_end, 
				 TTarget target_begin)
{
SEQAN_CHECKPOINT
	_arrayMoveForward_Default(source_begin, source_end, target_begin);	
}

//////////////////////////////////////////////////////////////////////////////
//arrayMoveBackward
//////////////////////////////////////////////////////////////////////////////

/**
.Function.arrayMoveBackward:
..cat:Array Handling
..summary:Moves a range of objects into another range of objects starting from the last element.
..signature:arrayMoveBackward(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source array.
..param.source_end:Iterator behind the last element of the source array.
...text:$source_end$ must have the same type as $source_begin$.
..param.target:Iterator to the first element of the target array.
...text:The target capacity should be at least as long as the source range.
..remarks:The function possibly clears (but does not destroy) the source elements.
	If source elements must persist, consider to use @Function.arrayCopyBackward@ instead.
..remarks.note:Be careful if source and target range overlap, because in this case
	some source elements could be accidently overwritten before they are moved.
..remarks.text:If source and target do not overlap, consider to use the function
@Function.arrayMoveForward@ instead that is faster in some cases.
..remarks.note:The semantic of this function's argument $target$ differ from the arguments of $::std::copy_backward$.
..see:Function.arrayMoveForward
..see:Function.arrayCopyBackward
..see:Class.SimpleType
*/
template<typename TTarget, typename TSource1, typename TSource2>
inline void 
_arrayMoveBackward_Default(TSource1 source_begin, 
						   TSource2 source_end, 
						   TTarget target_begin)
{
SEQAN_CHECKPOINT
	target_begin += (source_end - source_begin);
	while (source_end != source_begin)
	{
		--source_end;
		--target_begin;
		move(*target_begin, *source_end);
	}
}
template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayMoveBackward(TSource1 source_begin, 
				  TSource2 source_end, 
				  TTarget target_begin)
{
SEQAN_CHECKPOINT
	_arrayMoveBackward_Default(source_begin, source_end, target_begin);
}

//////////////////////////////////////////////////////////////////////////////
//arrayMove
//////////////////////////////////////////////////////////////////////////////

/**
.Function.arrayMove:
..cat:Array Handling
..summary:Moves a range of objects into another range of objects.
..signature:arrayMove(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source range.
..param.source_end:Iterator behind the last element of the source range.
...text:$source_end$ must have the same type as $source_begin$.
..param.target:Iterator to the first element of the target range.
...text:The target capacity should be at least as long as the source range.
..remarks:The function possibly clears (but does not destroy) the source elements.
	If source elements must persist, consider to use @Function.arrayCopy@ instead.
..remarks.text:If source and target range do not overlap, consider to use
	@Function.arrayMoveForward@ instead to improve performance.
..DISABLED.remarks.note:Be careful if source and target range overlap and the size of the
	source elements differ from the size of target elements, because in this case
	some source elements could be accidently overwritten before they are moved.
..remarks.note:Don't confuse this function with the standard $move$ function that
resembles @Function.arrayCopy@.
..see:Function.arrayMoveForward
..see:Function.arrayMoveBackward
..see:Function.arrayCopy
..see:Class.SimpleType
*/
template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayMove(TSource1 source_begin, 
		  TSource2 source_end,
		  TTarget target_begin)
{
	if ((void *) source_begin >= (void *) target_begin)
	{
SEQAN_CHECKPOINT
		arrayMoveForward(source_begin, source_end, target_begin);
	}
	else
	{
SEQAN_CHECKPOINT
		arrayMoveBackward(source_begin, source_end, target_begin);
	}
}

//////////////////////////////////////////////////////////////////////////////
//arrayClearSpace
//////////////////////////////////////////////////////////////////////////////

/**
.Function.arrayClearSpace:
..cat:Array Handling
..summary:Destroys the begin of an array and keeps the rest.
..signature:arrayClearSpace(arr_begin, arr_length, keep_from, move_to)
..param.arr_begin:Pointer to the first element of the array.
..param.arr_length:Length of the array.
..param.keep_from:Offset of the first object that will be kept.
..param.move_to:Offset the first kept object will get at the end of the function. 
..remarks.text:The objects $arr[keep_from]$ to $arr[arr_length-1]$
are moved to the area beginning at positions $move_to$. 
All objects in $arr[0]$ to $arr[keep_from-1]$ are destroyed.
After this function, the first $move_to$ positions of the array
are free and dont contain objects. 
..remarks.text:The array must have at least enough space to store $arr_length + move_to - keep_from$ objects.
..see:Function.arrayCopy
..see:Function.arrayDestruct
..see:Function.arrayCopyForward
..see:Class.SimpleType
*/
template <typename TIterator>
void _arrayClearSpace_Default(TIterator array_begin, 
							  size_t array_length, 
							  size_t keep_from, 
							  size_t move_to)
{
	if (keep_from == array_length)
	{
		arrayDestruct(array_begin, array_begin + array_length);
		return;
	}

	SEQAN_ASSERT(keep_from < array_length)

	if (keep_from == move_to)
	{
		arrayDestruct(array_begin, array_begin + move_to);
	}
	else if (keep_from < move_to) 
	{
		if (array_length > move_to)
		{
SEQAN_CHECKPOINT
			size_t middle = array_length - (move_to - keep_from);
			arrayConstructMove(array_begin + middle, array_begin + array_length, array_begin + array_length);
			arrayMove(array_begin + keep_from, array_begin + middle, array_begin + move_to);
			arrayDestruct(array_begin, array_begin + move_to);
		}
		else
		{
SEQAN_CHECKPOINT
			arrayConstructMove(array_begin + keep_from, array_begin + array_length, array_begin + move_to);
			arrayDestruct(array_begin, array_begin + array_length);
		}
	}
	else
	{
SEQAN_CHECKPOINT
		arrayMove(array_begin + keep_from, array_begin + array_length, array_begin + move_to);
		arrayDestruct(array_begin, array_begin + move_to);
		arrayDestruct(array_begin + array_length - (keep_from - move_to), array_begin + array_length);
	}
}

template <typename TIterator>
void arrayClearSpace(TIterator array_begin, 
					 size_t array_length, 
					 size_t keep_from, 
					 size_t move_to)
{
	_arrayClearSpace_Default(array_begin, array_length, keep_from, move_to);
}


//////////////////////////////////////////////////////////////////////////////
//BitsPerValue
//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.BitsPerValue:
..summary:Number of bits needed to store a value.
..signature:BitsPerValue<T>::VALUE
..param.T:A class.
..returns.param.VALUE:Number of bits needed to store $T$.
...default:$sizeof<T> * 8$
..see:Metafunction.ValueSize
*/
template <typename TValue>
struct BitsPerValue
{
	enum { VALUE = sizeof(TValue) * 8 };
};
template <typename TValue>
struct BitsPerValue<TValue const>:
	public BitsPerValue<TValue> {};

//////////////////////////////////////////////////////////////////////////////
//ValueSize
//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.ValueSize:
..summary:Number of different values a value type object can have.
..signature:ValueSize<T>::VALUE
..param.T:A class.
..returns.param.VALUE:Value size of $T$.
..remarks
...text:This function is only defined for integral types like $unsigned int$, $double$ or @Spec.Dna@.
..see:Metafunction.Value
*/
template <typename T>
struct ValueSize
{
	enum { VALUE = 1 << BitsPerValue<T>::VALUE };
};
template <typename TValue>
struct ValueSize<TValue const>:
	public ValueSize<TValue> {};



template < typename T >
struct _SupremumValueUnsigned {	static const T VALUE; };
template < typename T >
struct _SupremumValueSigned {	static const T VALUE; };

template < typename T >
struct _InfimumValueUnsigned {	static const T VALUE; };
template < typename T >
struct _InfimumValueSigned {	static const T VALUE; };



template < typename T >
const T _SupremumValueUnsigned<T>::VALUE = ~(T)0;
template < typename T >
const T _SupremumValueSigned<T>::VALUE = ( (((T)1 << (BitsPerValue<T>::VALUE - 2)) - 1) << 1) + 1;

template < typename T >
const T _InfimumValueUnsigned<T>::VALUE = 0;
template < typename T >
const T _InfimumValueSigned<T>::VALUE = ~(T)_SupremumValueSigned<T>::VALUE;



template < 
	typename T, 
	typename TParent = typename IF<
		TYPECMP< typename _MakeSigned<T>::Type, T >::VALUE,
		_SupremumValueSigned<T>,
		_SupremumValueUnsigned<T> >::Type >
struct SupremumValue:
	public TParent 
{
};

template < 
	typename T, 
	typename TParent = typename IF<
		TYPECMP< typename _MakeSigned<T>::Type, T >::VALUE,
		_InfimumValueSigned<T>,
		_InfimumValueUnsigned<T> >::Type >
struct InfimumValue:
	public TParent 
{
};

//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
