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
  $Id: sequence_interface.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEQUENCE_INTERFACE_H
#define SEQAN_HEADER_SEQUENCE_INTERFACE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Expand
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Overflow Strategy:
..summary:The strategy for resizing containers.
..tag.Insist:No capacity check. 
...remarks:The user has to ensure that the container's capacity is large enough. 
..tag.Limit:Limit the contents to current capacity. 
...remarks: All entries that exceed the capacity are lost. 
..tag.Exact:Expand as far as needed.
...remarks: The capacity is only changed if the current capacity is not large enough.
 If the capacity can only be expanded up to a certain ammount, it will be increased as far as possible
 and the contents are limited to the new capacity.
...remarks:Note that the capacity will never be shrinked. 
 Use @Function.shrinkToFit@ to resize the capacity down to the current length.
..tag.Generous:Expand if needed, get precautionary extra space. 
...remarks:Whenever the capacity has to be increased, the new capacity is choosen somewhat large than actually needed.
 This strategy limits the number of capacity changes, so that resizing takes armotized constant time.
 Use this strategy if the total amount of storage is unkown at first.  
...remarks:The new capacity is computed by @Function.computeGenerousCapacity@.
By default, it is guaranteed not to exceed about 
 tree halfs of the space that is used to store the data. 
 The user can overload @Function.computeGenerousCapacity@ in order to change this behavior.
..remarks:Changing the capacity of a container can invalidate the iterators of this container.
..remarks:If no overflow tag is specified, most operations use the default overflow strategy given by @Metafunction.DefaultOverflowImplicit@
or @Metafunction.DefaultOverflowExplicit@, depending on the kind of operation.
*/
struct TagInsist_;
typedef Tag<TagInsist_> const Insist;
typedef Tag<TagInsist_> const Tight;
//Insist INSIST;

struct TagLimit_;
typedef Tag<TagLimit_> const Limit;
//Limit LIMIT;

struct TagGenerous_;
typedef Tag<TagGenerous_> const Generous;
//Generous GENEROUS;

struct TagExact_;
typedef Tag<TagExact_> const Exact;
//Exact EXACT;

//____________________________________________________________________________

/**
.Metafunction.DefaultOverflowImplicit:
..hidefromindex
..summary:The default overflow strategy for implicit resize.
..signature:DefaultOverflowImplicit<T>::Type
..param.T:Type for which the overflow strategy is determined.
...type:Class.String
..returns.param.Type:Expansion tag for type of $T$.
..remarks:This function is used for functions that cause an implicit change of a container's size, like
e.g. @Function.assign@, @Function.append@, and @Function.replace@.
..see:Tag.Overflow Strategy
*/
template <typename T>
struct DefaultOverflowImplicit
{
	typedef Insist Type;
};

//____________________________________________________________________________

/**
.Metafunction.DefaultOverflowExplicit:
..hidefromindex
..summary:The default overflow strategy for explicit resize.
..signature:DefaultOverflowExplicit<T>::Type
..param.T:Type for which the overflow strategy is determined.
...type:Class.String
..returns.param.Type:Expansion tag for type of $T$.
..remarks:This function is used for functions that change a container's size explicit, like e.g. @Function.resize@.
..see:Tag.Overflow Strategy
..see:Metafunction.DefaultOverflowImplicit
*/
template <typename T>
struct DefaultOverflowExplicit
{
	typedef Exact Type;
};

//////////////////////////////////////////////////////////////////////////////
// IsContiguous
//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.IsContiguous:
..summary:Determines whether a container stores its elements in a contiguous array.
..signature:IsContiguous<T>::VALUE
..param.T:Type that is tested for being a string.
..returns.param.VALUE:$true$ if $T$ is a string, $false$ otherwise.
..remarks:Definition: A sequence container is "contiguous", if its elements
	are stored in a single contiguous array.
	Examples for contiguous sequences are @Spec.Alloc String@ or @Adaption.char array@.
..remarks:If an object $obj$ is a contiguous sequence, then $begin(obj)$ can be
	converted to a pointer to the first element of the content array.
*/
template <typename T>
struct IsContiguous
{
    typedef False Type;
	enum { VALUE = false };
};
template <typename T>
struct IsContiguous<T const>:
	public IsContiguous<T> {};

//////////////////////////////////////////////////////////////////////////////
// IsSequence
//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.IsSequence:
..summary:Determines whether a container stores its elements in sequential order.
..signature:IsSequence<T>::VALUE
..param.T:Type that is tested for being a sequence.
..returns.param.VALUE:$true$ if $T$ is a sequence, $false$ otherwise.
..remarks:For example @Class.String@ and @Class.Segment@ return $true$.
*/
template <typename T>
struct IsSequence
{
    typedef False Type;
	enum { VALUE = false };
};
template <typename T>
struct IsSequence<T const>:
	public IsSequence<T> {};

//////////////////////////////////////////////////////////////////////////////
// AllowsFastRandomAccess
//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.AllowsFastRandomAccess:
..summary:Determines whether a sequence efficiently supports random access.
..signature:AllowsFastRandomAccess<T>::VALUE
..param.T:Type that is tested for fast random access.
..returns.param.VALUE:$true$ if $T$ supports fast random access, $false$ otherwise.
..remarks:For example @Spec.Alloc String@, @Class.Segment@, and @Spec.Block String@ return $true$.
*/
template <typename T>
struct AllowsFastRandomAccess
{
    typedef True Type;
	enum { VALUE = true };
};
template <typename T>
struct AllowsFastRandomAccess<T const>:
	public AllowsFastRandomAccess<T> {};

//////////////////////////////////////////////////////////////////////////////
// identification
//////////////////////////////////////////////////////////////////////////////
/**
.Function.id:
..cat:Miscellaneous
..summary:A value that identifies the underlying sequence. 
..signature:void const * id(object)
..param.object:The object for which the id will be determined.
..returns:The id of $sequence$. 
..remarks.text:Two sequences should have the same id, if they share the same resource, e.g. the same memory buffer.
..remarks.text:The exact semantic of the returned id can vary for different classes.
Typically, the id of a string is a $void const *$ to the end of the string.
..remarks.note:The id of a single character need not to be the id of its container.
..example.code:String<char> str = "hallo seqan";
bool b1 = (id(str) == id(infix(str, 3, 7));   //true
bool b2 = (id(str) == id(String<char>(str))); //false
bool b3 = (id(str) == id(toCString(str))); 
..example.text:In this example, $b1$ is $true$, since the segment object returned by $infix()$
is just a filter and uses the buffer of it's host object $str$.
..example.text:$String<char>(str)$ constructs a temporary copy of $str$, so these two
strings have different id values.
..example.text:The result of the last comparison depends on the implementation of $toCString$
and cannot be predicted at compile time.
*/
template <typename T>
inline void const * 
id(T const & me)
{
SEQAN_CHECKPOINT
	return end(me, Standard());
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.shareResources:
..cat:Miscellaneous
..summary:Determines whether two sequences share the same resource.
..signature:bool shareResources(sequence1, sequence2)
..param.sequence1, sequence2:Two sequences.
..returns:$false$ if it can be guaranteed that $sequence1$ and $sequence2$ can be modified without changing each other, $true$ otherwise.
..remarks:Non-sequences are interpreted as sequences of size 1. 
..remarks:Note that this function may not work properly for argument types that are not listed here.
*/

template <typename T1, typename T2>
inline bool 
shareResources(T1 const & obj1,
			   T2 const & obj2)
{
SEQAN_CHECKPOINT
	return id(obj1) == id(obj2);
}

//////////////////////////////////////////////////////////////////////////////
// begin
//////////////////////////////////////////////////////////////////////////////

/**
.Function.begin:
..cat:Iteration
..cat:Containers
..summary:The begin of a container. 
..signature:Iterator begin(object [, tag])
..param.object:A container.
...type:Class.String
...concept:Concept.Container
..param.tag:An @Tag.Iterator Spec.iterator spec@ tag that specifies the kind of the iterator returned. (optional)
...default:Given by @Metafunction.DefaultGetIteratorSpec@.
..returns:An iterator to the first item in $object$. 
...metafunction:Metafunction.Iterator
..remarks.text:If the container does not contain any items at all, the function may return 0.
..see:Function.end
..see:Metafunction.Iterator
*/
template <typename T>
inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type 
begin(T & me)
{
SEQAN_CHECKPOINT
	return begin(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}
template <typename T>
inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type
begin(T const & me)
{
SEQAN_CHECKPOINT
	return begin(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}

//____________________________________________________________________________

//* ???Anti Default Sequences
template <typename T>
inline typename Iterator<T, Standard>::Type 
_begin_default(T & me,
			   Standard)
{
SEQAN_CHECKPOINT
	return & me;
}
template <typename T>
inline typename Iterator<T const, Standard>::Type 
_begin_default(T const & me,
			   Standard)
{
SEQAN_CHECKPOINT
	return & me;
}
//*/

//____________________________________________________________________________

template <typename T>
inline typename Iterator<T, Rooted>::Type 
_begin_default(T & me,
			   Rooted)
{
SEQAN_CHECKPOINT
	typedef typename Iterator<T, Rooted>::Type TIterator;
	return TIterator(me, begin(me, Standard()));
}
template <typename T>
inline typename Iterator<T const, Rooted>::Type 
_begin_default(T const & me,
			   Rooted)
{
SEQAN_CHECKPOINT
	typedef typename Iterator<T const, Rooted>::Type TIterator;
	return TIterator(me, begin(me, Standard()));
}
//____________________________________________________________________________

//folgende forward Deklaration wurde wegen Phaenomene bei VC++ 2003 hinzugenommen
//implemented in string_pointer.h
template <typename TValue>
inline typename Iterator<TValue const *, Standard>::Type  
begin(TValue const * me, 
	  Standard);

//____________________________________________________________________________

template <typename T, typename TSpec>
inline typename Iterator<T, Tag<TSpec> const>::Type 
begin(T & me,
	  Tag<TSpec> const tag_)
{
SEQAN_CHECKPOINT
	return _begin_default(me, tag_);
}
template <typename T, typename TSpec>
inline typename Iterator<T const, Tag<TSpec> const>::Type 
begin(T const & me,
	  Tag<TSpec> const tag_)
{
SEQAN_CHECKPOINT
	return _begin_default(me, tag_);
}


/*
template <typename TValue>
inline typename Iterator<TValue *, Standard>::Type  
begin(TValue * me, 
	  Standard)
{
SEQAN_CHECKPOINT
	return me;
}

//folgende Version wurde wegen eines seltsamen Phaenomens bei VC++ hinzugenommen
template <typename TValue>
inline typename Iterator<TValue const *, Standard>::Type  
begin(TValue const * me, 
	  Standard)
{
SEQAN_CHECKPOINT
	return me;
}

template <typename TValue, typename TSpec>
inline typename Iterator<TValue *, Standard>::Type  
begin(TValue * me, 
	  Tag<TSpec> const tag_)
//	  Standard)
{
SEQAN_CHECKPOINT
	return me;
}

//folgende Version wurde wegen eines seltsamen Phaenomens bei VC++ hinzugenommen
template <typename TValue, typename TSpec>
inline typename Iterator<TValue const *, Standard>::Type  
begin(TValue const * me, 
	  Tag<TSpec> const tag_)
//	  Standard)
{
SEQAN_CHECKPOINT
	return me;
}

*/



//////////////////////////////////////////////////////////////////////////////
// beginPosition
//////////////////////////////////////////////////////////////////////////////
/**
.Function.beginPosition:
..cat:Containers
..summary:Begin position of object in host.
..signature:Position beginPosition(object)
..param.object:An object.
...type:Class.String
...concept:Concept.Container
..returns:The position of the first item in $host(object)$ that belongs of $object$.
...metafunction:Metafunction.Position
..remarks
...text:For most classes $beginPosition$ always returns 0. Exceptions are e.g. @Spec.InfixSegment@ and @Spec.SuffixSegment@.
..see:Function.begin
*/
template <typename T>
inline typename Position<T>::Type 
beginPosition(T &)
{
SEQAN_CHECKPOINT
	return 0;
}
template <typename T>
inline typename Position<T>::Type 
beginPosition(T const &)
{
SEQAN_CHECKPOINT
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// end
//////////////////////////////////////////////////////////////////////////////

/**
.Function.end:
..cat:Iteration
..cat:Containers
..summary:The end of a container. 
..signature:Iterator end(object [, tag])
..param.object:A container.
...type:Class.String
...concept:Concept.Container
..param.tag:An @Tag.Iterator Spec.iterator spec@ tag that specifies the kind of the iterator returned. (optional)
...default:Given by @Metafunction.DefaultGetIteratorSpec@.
..returns:An iterator that points behind the last item in $object$.
...metafunction:Metafunction.Iterator
..remarks.text:If the container does not contain any items at all, the function may return 0.
..see:Function.begin
..see:Metafunction.Iterator
*/
template <typename T>
inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type 
end(T & me)
{
SEQAN_CHECKPOINT
	return end(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}
template <typename T>
inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type 
end(T const & me)
{
SEQAN_CHECKPOINT
	return end(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}

//____________________________________________________________________________

//* ???Anti Default Sequences
template <typename T>
inline typename Iterator<T, Standard>::Type 
_end_default(T & me,
			 Standard)
{
SEQAN_CHECKPOINT
	return (& me) + 1;
}
template <typename T>
inline typename Iterator<T const, Standard>::Type 
_end_default(T const & me,
			 Standard)
{
SEQAN_CHECKPOINT
	return (& me) + 1;
}
//*/

//____________________________________________________________________________

template <typename T>
inline typename Iterator<T, Rooted>::Type 
_end_default(T & me,
			 Rooted)
{
SEQAN_CHECKPOINT
	typedef typename Iterator<T, Rooted>::Type TIterator;
	return TIterator(me, end(me, Standard()));
}
template <typename T>
inline typename Iterator<T const, Rooted>::Type 
_end_default(T const & me,
			 Rooted)
{
SEQAN_CHECKPOINT
	typedef typename Iterator<T const, Rooted>::Type TIterator;
	return TIterator(me, end(me, Standard()));
}

//____________________________________________________________________________

template <typename T, typename TSpec>
inline typename Iterator<T, Tag<TSpec> const>::Type 
end(T & me,
	Tag<TSpec> const tag_)
{
SEQAN_CHECKPOINT
	return _end_default(me, tag_);
}
template <typename T, typename TSpec>
inline typename Iterator<T const, Tag<TSpec> const>::Type 
end(T const & me,
	Tag<TSpec> const tag_)
{
SEQAN_CHECKPOINT
	return _end_default(me, tag_);
}

//////////////////////////////////////////////////////////////////////////////
// endPosition
//////////////////////////////////////////////////////////////////////////////

/**
.Function.endPosition:
..cat:Containers
..summary:End position of object in host.
..signature:Position endPosition(object)
..param.object:An object.
...concept:Concept.Container
...type:Class.String
..returns:The position behind the last item in $host(object)$ that belongs of $object$.
...metafunction:Metafunction.Position
..see:Function.end
..see:Function.beginPosition
*/
template <typename T>
inline typename Position<T>::Type 
endPosition(T & me)
{
SEQAN_CHECKPOINT
	return length(me);
}
template <typename T>
inline typename Position<T>::Type 
endPosition(T const & me)
{
SEQAN_CHECKPOINT
	return length(me);
}

//////////////////////////////////////////////////////////////////////////////
// value
//////////////////////////////////////////////////////////////////////////////

/**
.Function.value:
..cat:Iteration
..cat:Containers
..summary:Reference to the value.
..signature:Reference value(container, position)
..param.container:A container of values.
..param.position:A position in $container$ on which the value should be accessed.
..returns:A reference or proxy to the value.
...metafunction:Metafunction.Reference
*/

//* ???Anti Default Sequences
template <typename T, typename TPos>
inline typename Reference<T>::Type
value(T & me, 
	  TPos /*pos*/)
{
SEQAN_CHECKPOINT
	return me;
} 
template <typename T, typename TPos>
inline typename Reference<T const>::Type
value(T const & me, 
	  TPos /*pos*/)
{
SEQAN_CHECKPOINT
	return me;
} 
//*/

//////////////////////////////////////////////////////////////////////////////
// getValue
//////////////////////////////////////////////////////////////////////////////

/**
.Function.getValue:
..summary:Access to the value.
..cat:Containers
..cat:Content Manipulation
..signature:GetValue getValue(container, pos)
..param.container:A container.
...concept:Concept.Container
..param.pos:The position of an item in $object$.
...remarks:$pos$ should be convertible to $Position<T>::Type$ for $container$-type $T$.
..returns:The item at position $pos$ in $container$.
This can either be a reference to the item or a temporary copy of the item.
...metafunction:Metafunction.GetValue
..remarks:
...text:If $pos$ is out of range, then the behavior of the function is undefined.
..see:Metafunction.GetValue
..see:Metafunction.Position
..see:Function.value
*/

template <typename T, typename TPos>
inline typename GetValue<T>::Type
getValue(T & me, 
		 TPos pos)
{
SEQAN_CHECKPOINT
	return (typename GetValue<T>::Type) value(me, pos);
} 
template <typename T, typename TPos>
inline typename GetValue<T const>::Type
getValue(T const & me, 
		 TPos pos)
{
SEQAN_CHECKPOINT
	return value(me, pos);
} 

//////////////////////////////////////////////////////////////////////////////
// front
//////////////////////////////////////////////////////////////////////////////

/**
.Function.Container#front:
..cat:Containers
..summary:The first item in container. 
..signature:Iterator front(container)
..param.container:A container.
...concept:Concept.Container
..returns:A @Metafunction.Reference.reference@ of the first item in $container$.
...metafunction:Metafunction.Reference
..remarks:This function is equivalent to $value(me, beginPosition(me))$.
..see:Function.value
..see:Function.begin
*/


template <typename T>
inline typename Reference<T>::Type
front(T & me)
{
SEQAN_CHECKPOINT
    return value(me, beginPosition(me));
}
template <typename T>
inline typename Reference<T const>::Type
front(T const & me)
{
SEQAN_CHECKPOINT
    return value(me, beginPosition(me));
}

//////////////////////////////////////////////////////////////////////////////
// back
//////////////////////////////////////////////////////////////////////////////

/**
.Function.back:
..cat:Containers
..summary:The last item in container. 
..signature:Iterator back(container)
..param.container:A container.
...concept:Concept.Container
..returns:A @Metafunction.Reference.reference@ of the last item in $container$.
...metafunction:Metafunction.Reference
..remarks:This function is equivalent to $value(me, endPosition(me) - 1)$.
..see:Function.value
..see:Function.end
..see:Function.Container#front
*/

template <typename T>
inline typename Reference<T const>::Type
back(T const & me)
{
SEQAN_CHECKPOINT
    return value(me, endPosition(me) - 1); 
}

template <typename T>
inline typename Reference<T>::Type
back(T & me)
{
SEQAN_CHECKPOINT
    return value(me, endPosition(me) - 1);
}

//////////////////////////////////////////////////////////////////////////////
// iter
//////////////////////////////////////////////////////////////////////////////

/**
.Function.iter:
..cat:Containers
..summary:Iterator to item at given position. 
..signature:Iterator iter(object, pos [, tag])
..param.object:A container.
...type:Class.String
..param.pos:The position of an item in $object$.
...metafunction:Metafunction.Position
..param.tag:An @Tag.Iterator Spec.iterator spec@ tag that specifies the kind of the iterator returned. (optional)
...default:Given by @Metafunction.DefaultGetIteratorSpec@.
..returns:An iterator to the item at position $pos$ in $object$.
...metafunction:Metafunction.Iterator
..remarks:
...text:If $pos$ is out of range, then the behavior of the function is undefined.
..see:Function.value
..see:Metafunction.Iterator
..see:Metafunction.Position
*/

template <typename T, typename TPos>
inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type
iter(T & me, 
	 TPos pos)
{
SEQAN_CHECKPOINT
	return iter(me, pos, typename DefaultGetIteratorSpec<T>::Type());
} 
template <typename T, typename TPos>
inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type
iter(T const & me, 
	 TPos pos)
{
SEQAN_CHECKPOINT
	return iter(me, pos, typename DefaultGetIteratorSpec<T const>::Type());
} 


template <typename T, typename TPos, typename TTag>
inline typename Iterator<T, Tag<TTag> const>::Type
iter(T & me, 
	 TPos pos,
	 Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return begin(me, tag_) + pos;
} 
template <typename T, typename TPos, typename TTag>
inline typename Iterator<T const, Tag<TTag> const>::Type
iter(T const & me, 
	 TPos pos,
	 Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return begin(me, tag_) + pos;
} 

//////////////////////////////////////////////////////////////////////////////
// assignValue (3)
//////////////////////////////////////////////////////////////////////////////
/**
.Function.assignValue:
..cat:Content Manipulation
..signature:assignValue(container, pos, value)
..param.container:A container.
...concept:Concept.Container
..param.pos:Position of the item in $container$ to that $value$ is assigned.
..remarks:If $object$ is a container (that is $pos$ is not specified), 
	the whole content of $object$ is replaced by $value$.
..remarks.text:
	If $value$ is not used again after calling this function, 
	then consider to use @Function.moveValue@ that could be faster in some cases instead.
*/

template <typename T, typename TValue, typename TPos>
inline void
assignValue(T & me,
			TPos pos, 
			TValue const & _value)
{
SEQAN_CHECKPOINT
	assign(value(me, pos), _value);
} 

//////////////////////////////////////////////////////////////////////////////
// moveValue (3)
//////////////////////////////////////////////////////////////////////////////
/**
.Function.moveValue:
..cat:Content Manipulation
..signature:moveValue(container, pos, value)
..param.object:
...concept:Concept.Container
..param.container:A container.
...concept:Concept.Container
..param.pos:Position of the item in $container$ to that $value$ is moved to.
..remarks:If $object$ is a container (that is $pos$ is not specified), 
the whole content of $object$ is replaced by $value$.
..remarks.text:
	This function possibly clears $value$.
	If $value$ should be used further, consider to use @Function.assignValue@ instead.
*/

template <typename T, typename TValue, typename TPos>
inline void
moveValue(T & me,
		  TPos pos, 
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	move(value(me, pos), _value);
} 

//////////////////////////////////////////////////////////////////////////////
// length
//////////////////////////////////////////////////////////////////////////////
/**
.Function.length:
..cat:Containers
..summary:The number of items/characters. 
..signature:Size length(object)
..param.object:A container.
...concept:Concept.Container
..returns:The number of items/characters in $object$.
...metafunction:Metafunction.Size
..remarks.text:The length of a sequence can never exceed it's capacity.
..see:Function.capacity
*/

//* ???Anti Default Sequences
template <typename T> 
inline typename Size<T const>::Type
length(T const & /*me*/)
{
SEQAN_CHECKPOINT
	return 1;
}
//*/

//////////////////////////////////////////////////////////////////////////////
// capacity
//////////////////////////////////////////////////////////////////////////////

/**
.Function.capacity:
..cat:Containers
..summary:The maximal length.
..signature:Size capacity(object)
..param.object:A container.
...remarks: If $object$ cannot be converted to one of these types, the function returns 1.
..returns:The maximal number of items/characters that can be stored in $object$.
...metafunction:Metafunction.Size
..remarks.text:The size of a sequence can never exceed it's capacity, but some containers support 
resizing of the capacity. 
Some functions do that implicitely if they are called with a suitable @Tag.Overflow Strategy.overflow strategy@.
The function @Function.reserve@ can be used to change the capacity explicitely.
*/
template <typename T> 
inline typename Size<T const>::Type
capacity(T const & me)
{
SEQAN_CHECKPOINT
	return length(me);
}

//////////////////////////////////////////////////////////////////////////////
// empty
//////////////////////////////////////////////////////////////////////////////

/**
.Function.empty:
..cat:Containers
..summary:Test a container for being empty.
..signature:bool empty(object)
..param.object:A container.
..returns:$true$ if $object$ contains no elements, otherwise $false$.
..remarks.text:$empty(x)$ is guaranteed to be at least as fast as $length(me) == 0$, 
but can be significantly faster in some cases.
..see:Function.length
*/
template <typename T>
inline bool
empty(T const & me)
{
SEQAN_CHECKPOINT
	return (length(me) == 0);
}

//////////////////////////////////////////////////////////////////////////////
// _computeSize4Capacity
//////////////////////////////////////////////////////////////////////////////

// note: for value types of size 1 or 2,
// an extra position for the termination character is allocated.
// This speeds up a conversion to a c style string (see Spec.CStyle String)
// note that this extra position is necessary not only for char and wchar_t,
// but also for all other value types of size 1 and 2 to make the application
// of the funciton move for in-place alphabet conversion.


template <typename T, typename TSize>
inline TSize 
_computeSize4Capacity(T const & /*me*/, 
					  TSize capacity)
{
SEQAN_CHECKPOINT
	if (sizeof(T) <= 2) return capacity + 1;
	else return capacity;
}


//////////////////////////////////////////////////////////////////////////////
// computeGenerousCapacity
//////////////////////////////////////////////////////////////////////////////
/**
.Function.computeGenerousCapacity:
..hidefromindex
..cat:Containers
..summary:Capacity for generous expansion.
..signature:Size computeGenerousCapacity(container, capacity)
..param.container:A container that should be expanded.
..param.capacity:Minimal capacity needed.
..returns:A value larger than $capacity$ that should be used as new capacity for $container$
when it is expanded using the @Tag.Overflow Strategy."Generous" overflow strategy@.
...metafunction:Metafunction.Size
..see:Tag.Overflow Strategy
*/
template <typename T, typename TSize>
inline TSize 
computeGenerousCapacity(T const & /*me*/, 
						 TSize capacity)
{
SEQAN_CHECKPOINT
	if (capacity <= 32) return 32;
	return capacity + (capacity >> 1);
}

//////////////////////////////////////////////////////////////////////////////
// _storageUpdated
//////////////////////////////////////////////////////////////////////////////

/*
template <typename T>
inline void 
_storageUpdated(T & me,
				void const *)
{
}

template <typename T>
inline void 
_storageUpdated(T & me)
{
	_storageUpdated_(me, (T *) 0);
}

template <typename T>
inline void 
_storageUpdated(T const & me)
{
	_storageUpdated_(me, (T const *) 0);
}
*/


//////////////////////////////////////////////////////////////////////////////
// assign
//////////////////////////////////////////////////////////////////////////////

template<typename TTarget, typename TSource>
inline void 
assign(TTarget & target, 
	   TSource & source,
	   typename Size<TTarget>::Type limit)
{
SEQAN_CHECKPOINT
	assign(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTarget, typename TSource>
inline void 
assign(TTarget const & target, 
	   TSource & source,
	   typename Size<TTarget>::Type limit)
{
SEQAN_CHECKPOINT
	assign(target, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}
template<typename TTarget, typename TSource>
inline void 
assign(TTarget & target, 
	   TSource const & source,
	   typename Size<TTarget>::Type limit)
{
SEQAN_CHECKPOINT
	assign(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTarget, typename TSource>
inline void 
assign(TTarget const & target, 
	   TSource const & source,
	   typename Size<TTarget>::Type limit)
{
SEQAN_CHECKPOINT
	assign(target, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}

//////////////////////////////////////////////////////////////////////////////
// append
//////////////////////////////////////////////////////////////////////////////

/**
.Function.append:
..summary:Concatenate two containers.
..cat:Content Manipulation
..signature:append(target, source [, limit] [,resize_tag])
..param.target: A container $source$ is append to.
..param.source: A container that is append to $target$.
...remarks:The function does not modify this container.
..param.limit: The maximal length of $target$ after the operation. (optional)
..param.resize_tag: Specifies the strategy that is applied if $target$ has not enough capacity to store the complete content. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowImplicit@ of the $target$ type.
..remarks:The result of this operation is stored in $target$.
..see:Function.assign
*/
template<typename TTarget, typename TSource>
inline void 
append(TTarget & target, 
	   TSource & source)
{
SEQAN_CHECKPOINT
	append(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTarget, typename TSource>
inline void 
append(TTarget const & target, 
	   TSource & source)
{
SEQAN_CHECKPOINT
	append(target, source, typename DefaultOverflowImplicit<TTarget const>::Type());
}
template<typename TTarget, typename TSource>
inline void 
append(TTarget & target, 
	   TSource const & source)
{
SEQAN_CHECKPOINT
	append(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTarget, typename TSource>
inline void 
append(TTarget const & target, 
	   TSource const & source)
{
SEQAN_CHECKPOINT
	append(target, source, typename DefaultOverflowImplicit<TTarget const>::Type());
}

//____________________________________________________________________________

template<typename TTarget, typename TSource>
inline void 
append(TTarget & target, 
	   TSource & source,
	   typename Size<TTarget>::Type limit)
{
SEQAN_CHECKPOINT
	append(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTarget, typename TSource>
inline void 
append(TTarget const & target, 
	   TSource & source,
	   typename Size<TTarget>::Type limit)
{
SEQAN_CHECKPOINT
	append(target, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}
template<typename TTarget, typename TSource>
inline void 
append(TTarget & target, 
	   TSource const & source,
	   typename Size<TTarget>::Type limit)
{
SEQAN_CHECKPOINT
	append(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTarget, typename TSource>
inline void 
append(TTarget const & target, 
	   TSource const & source,
	   typename Size<TTarget>::Type limit)
{
SEQAN_CHECKPOINT
	append(target, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}


//////////////////////////////////////////////////////////////////////////////
// appendValue
//////////////////////////////////////////////////////////////////////////////
/**
.Function.appendValue:
..signature:appendValue(target, value [, resize_tag])
..cat:Content Manipulation
..summary:Appends a value to a container.
..param.target:A container.
..param.value:Value that is appended to $target$.
..param.resize_tag:
..param.resize_tag: Specifies the strategy that is applied if $target$ has not enough capacity to store the complete content. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowImplicit@ of the $target$ type.
*/

template <typename T, typename TValue>
inline void
appendValue(T & me, 
			TValue const & _value)
{
SEQAN_CHECKPOINT
	appendValue(me, _value, typename DefaultOverflowImplicit<T>::Type());
} 
template <typename T, typename TValue>
inline void
appendValue(T const & me, 
			TValue const & _value)
{
SEQAN_CHECKPOINT
	appendValue(me, _value, typename DefaultOverflowImplicit<T const>::Type());
} 

//////////////////////////////////////////////////////////////////////////////
// insertValue
//////////////////////////////////////////////////////////////////////////////
/**
.Function.insertValue:
..cat:Content Manipulation
..summary:Inserts a single value into a container.
..signature:insertValue(target, pos, value [, resize_tag])
..param.target:The container
..param.pos:Position within $target$ at which $value$ is to be inserted.
..param.value:Value that will be inserted into $target$.
..param.resize_tag:Strategy that is applied if $target$ has not enough capacity to store the complete content.
...type:Tag.Overflow Strategy
..see:Function.assignValue
..see:Function.appendValue
*/

template <typename T, typename TPosition, typename TValue>
inline void
insertValue(T & me, 
			TPosition pos, 
			TValue const & _value)
{
SEQAN_CHECKPOINT
	insertValue(me, pos, _value, typename DefaultOverflowImplicit<T>::Type());
} 
template <typename T, typename TPosition, typename TValue>
inline void
insertValue(T const & me, 
			TPosition pos, 
			TValue const & _value)
{
SEQAN_CHECKPOINT
	insertValue(me, pos, _value, typename DefaultOverflowImplicit<T const>::Type());
} 

//////////////////////////////////////////////////////////////////////////////
// replace
//////////////////////////////////////////////////////////////////////////////

/**
.Function.replace:
..summary:Replaces a part of a container with another container.
..cat:Content Manipulation
..signature:replace(target, pos_begin, pos_end, source [, limit] [,resize_tag])
..param.target: A container that is modified.
..param.pos_begin: Begin of replaced area.
...text:The first position in $target$ of the area that is replaced by $source$.
..param.pos_end: End of replaced area.
...text:The position behind the last position in $target$ of the area that is replaced by $source$.
..param.source: A container that is inserted into $target$.
...remarks:The function does not modify this container.
..param.limit: The maximal length of $target$ after the operation. (optional)
..param.resize_tag: Specifies the strategy that is applied if $target$ has not enough capacity to store the complete content. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowImplicit@ of the $target$ type.
..see:Function.assign
..see:Function.append
..remarks.text:Some compilers have difficulties if $pos_begin$ and $pos_end$ are both 0, since 0 can be
both a position or an iterator. The workaround is to convert at least one of these arguments
explicite to the position or to the interator type.
*/
template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void 
replace(TTarget & target,
		TPositionBegin pos_begin, 
		TPositionEnd pos_end,
		TSource & source)
{
	replace(target, pos_begin, pos_end, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void 
replace(TTarget const & target,
		TPositionBegin pos_begin, 
		TPositionEnd pos_end,
		TSource & source)
{
	replace(target, pos_begin, pos_end, source, typename DefaultOverflowImplicit<TTarget const>::Type());
}
template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void 
replace(TTarget & target,
		TPositionBegin pos_begin, 
		TPositionEnd pos_end,
		TSource const & source)
{
	replace(target, pos_begin, pos_end, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void 
replace(TTarget const & target,
		TPositionBegin pos_begin, 
		TPositionEnd pos_end,
		TSource const & source)
{
	replace(target, pos_begin, pos_end, source, typename DefaultOverflowImplicit<TTarget const>::Type());
}

//____________________________________________________________________________

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void 
replace(TTarget & target,
		TPositionBegin pos_begin, 
		TPositionEnd pos_end,
		TSource & source,
		typename Size<TTarget>::Type limit)
{
	replace(target, pos_begin, pos_end, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void 
replace(TTarget const & target,
		TPositionBegin pos_begin, 
		TPositionEnd pos_end,
		TSource & source,
		typename Size<TTarget>::Type limit)
{
	replace(target, pos_begin, pos_end, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}
template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void 
replace(TTarget & target,
		TPositionBegin pos_begin, 
		TPositionEnd pos_end,
		TSource const & source,
		typename Size<TTarget>::Type limit)
{
	replace(target, pos_begin, pos_end, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void 
replace(TTarget const & target,
		TPositionBegin pos_begin, 
		TPositionEnd pos_end,
		TSource const & source,
		typename Size<TTarget>::Type limit)
{
	replace(target, pos_begin, pos_end, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}


//////////////////////////////////////////////////////////////////////////////
// reserve
//////////////////////////////////////////////////////////////////////////////

/**
.Function.reserve:
..cat:Containers
..summary:Increases the capacity.
..signature:Size reserve(object, new_capacity [, resize_tag])
..param.object: A container.
..param.new_capacity: The new capacity $object$ will get.
..param.resize_tag: Specifies the strategy that is applied for changing the capacity. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowExplicit@.
..returns:The ammout of the requested capacity that was available.
That is the function returns the minimum of $new_capacity$ and $capacity(me)$.
...metafunction:Metafunction.Size
..remarks:At the end of the operation, $capacity(me)$ can be larger than $new_capacity$.
If $new_capacity$ is smaller than $capacity(me)$ at the beginning of the operation, 
the operation need not to change the capacity at all.
..remarks:This operation does not changes the content of $object$.
...note:This operation may invalidate iterators of $object$.
..see:Function.capacity
*/
template <typename T, typename TSize>
inline typename Size<T>::Type 
reserve(
	T & /*me*/, 
	TSize new_capacity,
	Insist)
{
SEQAN_CHECKPOINT
	return new_capacity;
}

template <typename T, typename TSize>
inline typename Size<T>::Type 
reserve(
	T & me, 
	TSize new_capacity,
	Limit)
{
SEQAN_CHECKPOINT
	typename Size<T>::Type me_capacity = capacity(me);
	if (me_capacity < (typename Size<T>::Type) new_capacity) return me_capacity;
	return new_capacity;
}

template <typename T, typename TSize>
inline typename Size<T>::Type 
reserve(
	T & me, 
	TSize new_capacity)
{
SEQAN_CHECKPOINT
	return reserve(me, new_capacity, typename DefaultOverflowExplicit<T>::Type());
}

//////////////////////////////////////////////////////////////////////////////
// resize
//////////////////////////////////////////////////////////////////////////////

/**
.Function.resize:
..cat:Containers
..summary:Changes the length.
..signature:Size resize(object, new_length [, resize_tag])
..param.object: A container.
...type:Class.String
..param.new_length: The new length $object$ will get.
..param.resize_tag: Specifies the strategy that is applied if the capacity of $object$ is less than $new_length$. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowExplicit@.
..returns:The new length $length(object)$.
...metafunction:Metafunction.Size
..remarks:This function can be used both for expanding and for shrinking $object$.
...text:
If $new_length$ is too large for $object$ (i.e. the @Function.capacity@ of $object$ is too small), then $object$ is
expanded as far as possible. The resulting length
of $object$ could be less than $new_length$, depending on the type of $object$ and the available storage.
..see:Function.length
..see:Function.reserve
*/
template <typename T, typename TSize>
inline typename Size<T>::Type  
resize(
	T & me, 
	TSize new_length)
{
SEQAN_CHECKPOINT
	return resize(me, new_length, typename DefaultOverflowExplicit<T>::Type());
}

//////////////////////////////////////////////////////////////////////////////
// resizeSpace
//////////////////////////////////////////////////////////////////////////////

/**
.Function.resizeSpace:
..cat:Containers
..summary:Makes free space in container
..signature:Size resizeSpace(object, size, pos_begin, pos_end [, limit] [, resize_tag])
..param.object:The container.
...type:Class.String
..param.size:Number of characters that should be freed.
..param.pos_begin:Position of the first item in $object$ that is to be destroyed.
..param.pos_end:Position behind the last item in $object$ that is to be destroyed.
...remarks:If $pos_end == pos_begin$, no item in $object$ will be destroyed.
..param.limit:Maximal length $object$ can get after this operation. (optional)
..param.resize_tag:Strategy that is applied if $object$ has not enough capacity to store the complete content. (optional)
...metafunction:Metafunction.DefaultOverflowExplicit
..returns:The number of free characters.
...metafunction:Metafunction.Size
...remarks:Depeding on the @Tag.Overflow Strategy.overflow strategy@ specified by $resize_tag$,
this could be $size$ or less than $size$ if $object$ has not enough @Function.capacity@.
*/

template<typename T, typename TSize, typename TPosition>
inline TSize 
resizeSpace(T & me, 
			TSize size, 
			TPosition pos_begin, 
			TPosition pos_end)
{
SEQAN_CHECKPOINT
	return resizeSpace(me, size, pos_begin, pos_end, typename DefaultOverflowExplicit<T>::Type());
}

template<typename T, typename TSize, typename TPosition>
inline TSize 
resizeSpace(T & me, 
			TSize size, 
			TPosition pos_begin, 
			TPosition pos_end,
			TSize limit)
{
SEQAN_CHECKPOINT
	return resizeSpace(me, size, pos_begin, pos_end, limit, typename DefaultOverflowExplicit<T>::Type());
}

//////////////////////////////////////////////////////////////////////////////
// erase
//////////////////////////////////////////////////////////////////////////////

/**
.Function.erase:
..summary:Erases a part of a container
..cat:Containers
..signature:erase(object, pos [, pos_end])
..param.object:The container.
...type:Class.String
..param.pos:Position of the first item in $object$ that is to be destroyed.
..param.pos_end:Position behind the last item in $object$ that is to be destroyed. (optional)
...default:$pos + 1$
...remarks:If $pos_end$ is omitted, only one element in $object$ at position $pos$ is destroyed.
..remarks:$erase(object, pos, pos_end)$ is semantically the same as @Function.resizeSpace.resizeSpace(object, 0, pos, pos_end)@.
*/

template<typename T, typename TPosition>
inline void 
erase(T & me, 
	  TPosition pos, 
	  TPosition pos_end)
{
SEQAN_CHECKPOINT
	resizeSpace(me, 0, pos, pos_end);
}

template<typename T, typename TPosition>
inline void 
erase(T & me, 
	  TPosition pos)
{
SEQAN_CHECKPOINT
	resizeSpace(me, 0, pos, pos + 1);
}

//////////////////////////////////////////////////////////////////////////////
// fill
//////////////////////////////////////////////////////////////////////////////

/**
.Function.fill:
..cat:Containers
..summary:Resizes and fills a container.
..signature:Size fill(object, new_length, value [, resize_tag])
..param.object: A container.
...type:Class.String
..param.new_length: The new length $object$ will get.
..param.value: Value that is copied if new items are created in $object$.
...remarks:If the $value$ argument is omitted, the default constructor is used to create
new items in $object$.
..param.resize_tag: Specifies the strategy that is applied if the capacity of $object$ is less than $new_length$. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowExplicit@.
..returns:The new length $length(object)$.
...metafunction:Metafunction.Size
..remarks:This function can be used both for expanding and for shrinking $object$.
...text:
	If $new_length$ is too large for $object$ (i.e. the @Function.capacity@ of $object$ is too small), then $object$ is
	expanded as far as possible. The resulting length
	of $object$ could be less than $new_length$, depending on the type of $object$ and the available storage.
..see:Function.length
..see:Function.reserve
..see:Function.resize
*/

template <typename T, typename TSize, typename TValue>
inline typename Size<T>::Type  
fill(
	T & me, 
	TSize new_length,
	TValue const & val)
{
SEQAN_CHECKPOINT
	return fill(me, new_length, val, typename DefaultOverflowExplicit<T>::Type());
}

//////////////////////////////////////////////////////////////////////////////
// shrinkToFit
//////////////////////////////////////////////////////////////////////////////

/**
.Function.shrinkToFit:
..cat:Containers
..summary:Resizes container to minimum capacity
..signature:shrinkToFit(object) 
..param.object: A container.
..remarks
...text:$shrinkToFit(object)$ is equivalent to $reserve(object, length(object), Exact())$.
..see:Function.capacity
..see:Function.length
..see:Function.reserve
*/

template <typename T, typename TSize, typename TValue>
inline void  
shrinkToFit(T & me)
{
SEQAN_CHECKPOINT
	reserve(me, length(me), Exact());
}
//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
