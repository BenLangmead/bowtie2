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
  $Id: basic_allocator_interface.h,v 1.1 2008/08/25 16:20:02 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_ALLOCATOR_INTERFACE_H
#define SEQAN_HEADER_BASIC_ALLOCATOR_INTERFACE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//Allocator
//////////////////////////////////////////////////////////////////////////////


/**
.Class.Allocator:
..cat:Basic
..summary:Manager for allocated memory.
..signature:Allocator<TSpec>
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
..implements:Concept.Allocator
..include:basic.h
..remarks:There are two reasons for using non-trivial allocators:
...text:1. Allocators support the function @Function.Allocator#clear@ for a fast deallocation of all 
allocated memory blocks. 
...text:2. Some allocators are faster in allocating an deallocating memory.
Pool allocators like e.g. @Spec.Single Pool Allocator@ or @Spec.Multi Pool Allocator@
speed up @Function.allocate@, @Function.deallocate@, and @Function.Allocator#clear@ for
pooled memory blocks.
*/

template <typename TSpec>
struct Allocator;

///.Function.allocate.param.object.type:Class.Allocator
///.Function.deallocate.param.object.type:Class.Allocator


//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////

//.Metafunction.Spec.param.T.type:Class.Allocator

template <typename TSpec>
struct Spec<Allocator<TSpec> >
{
	typedef TSpec Type;
};


//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Allocator Usage:
..summary:The purpose of an allocated memory block.
..tag.TagAllocateTemp:Temporary memory. 
..tag.TagAllocateStorage:Memory for storing container content. 
..see:Function.allocate
..see:Function.deallocate
*/
struct TagAllocateUnspecified_; //< usage not specified
typedef Tag<TagAllocateUnspecified_> const TagAllocateUnspecified;

struct TagAllocateTemp_; //< allocate temporary memory
typedef Tag<TagAllocateTemp_> const TagAllocateTemp;

struct TagAllocateStorage_; //< allocate memory for storing member data
typedef Tag<TagAllocateStorage_> const TagAllocateStorage;


//////////////////////////////////////////////////////////////////////////////
//allocates memory on heap. No c'tors are called.
//////////////////////////////////////////////////////////////////////////////

/**
.Function.allocate:
..cat:Memory
..summary:Allocates memory from heap.
..signature:allocate(object, data, count [, usage_tag])
..param.object:Allocator object.
...remarks:$object$ is conceptually the "owner" of the allocated memory.
 Objects of all types can be used as allocators. If no special behavior is implemented,
 default functions allocation/deallocation are applied that uses standard
 $new$ and $delete$ operators.
..param.count:Number of items that could be stored in the allocated memory.
...text:The type of the allocated items is given by the type of $data$.
..param.usage_tag:A tag the specifies the purpose for the allocated memory.
...value:@Tag.Allocator Usage@
..returns.param.data:Pointer to allocated memory.
...remarks:The value of this pointer is overwritten by the function.
..remarks:
...text:The function allocates at least $count*sizeof(data)$ bytes. 
 The allocated memory is large enough 
 to hold $count$ objects of type $T$, where $T *$ is type of $data$.
...note:These objects are not constructed by $allocate$.
...text:Use e.g. one of the functions @Function.valueConstruct@, @Function.arrayConstruct@, @Function.arrayConstructCopy@ or @Function.arrayFill@
to construct the objects.
A $new$ operator which is part of the C++ standard (defined in $<new>$)
 can also be used to construct objects at a given memory address.
..note:All allocated memory blocks should be deallocated by the corresponding function @Function.deallocate@.
..see:Function.deallocate
..see:Function.valueConstruct
..see:Function.arrayFill
..see:Function.arrayConstruct
..see:Function.arrayConstructCopy
*/
template <typename T, typename TValue, typename TSize>
inline void
allocate(T const & me,
		 TValue * & data,
		 TSize count)
{
	allocate(me, data, count, TagAllocateUnspecified());
}
template <typename T, typename TValue, typename TSize>
inline void
allocate(T & me,
		 TValue * & data,
		 TSize count)
{
	allocate(me, data, count, TagAllocateUnspecified());
}

template <typename T, typename TValue, typename TSize, typename TUsage>
inline void
allocate(T const &, 
		 TValue * & data,
		 TSize count,
		 Tag<TUsage> const)
{
	data = (TValue *) operator new(count * sizeof(TValue));
	if (data)
	    SEQAN_PROADD(SEQAN_PROMEMORY, count * sizeof(TValue));
}
template <typename T, typename TValue, typename TSize, typename TUsage>
inline void
allocate(T &, 
		 TValue * & data,
		 TSize count,
		 Tag<TUsage> const)
{
	data = (TValue *) operator new(count * sizeof(TValue));
	if (data)
	    SEQAN_PROADD(SEQAN_PROMEMORY, count * sizeof(TValue));
}


//////////////////////////////////////////////////////////////////////////////
//deallocates memory that was allocates using allocate(.)

/**
.Function.deallocate:
..cat:Memory
..summary:Deallocates memory.
..signature:deallocate(object, data, count [, usage_tag])
..param.object:Allocator object.
...remarks:$object$ is conceptually the "owner" of the allocated memory.
 Objects of all types can be used as allocators. If no special behavior is implemented,
 default functions allocation/deallocation are applied that uses standard
 $new$ and $delete$ operators.
..param.data:Pointer to allocated memory that was allocated by $allocate$.
..param.count:Number of items that could be stored in the allocated memory.
..param.usage_tag:A tag the specifies the purpose for the allocated memory.
...value:@Tag.Allocator Usage@
..remarks:
...text:The values for $object$, $count$ and $usage_tag$ should be the same that was 
used when $allocate$ was called. The value of $data$ should be the same that was
returned by $allocate$.
...note:$deallocate$ does not destruct objects.
...text:Use e.g. one of the functions @Function.valueDestruct@ or @Function.arrayDestruct@ to destruct the objects.
$delete$ and $delete []$ operators which are part of the C++ standard (defined in $<new>$)
 can also be used to destruct objects at a given memory address.
..see:Function.valueDestruct
..see:Function.arrayDestruct
*/
template <typename T, typename TValue, typename TSize>
inline void 
deallocate(T const & me, 
		   TValue * data, 
		   TSize const count)
{
	deallocate(me, data, count, TagAllocateUnspecified());
}
template <typename T, typename TValue, typename TSize>
inline void 
deallocate(T & me, 
		   TValue * data, 
		   TSize const count)
{
	deallocate(me, data, count, TagAllocateUnspecified());
}

template <typename T, typename TValue, typename TSize, typename TUsage>
inline void 
deallocate(T const & /*me*/,
		   TValue * data, 
		   TSize count,
		   Tag<TUsage> const)
{
	if (data && count)	// .. to use count if SEQAN_PROFILE is not defined
	    SEQAN_PROSUB(SEQAN_PROMEMORY, count * sizeof(TValue));
	operator delete ((void *) data);
}
template <typename T, typename TValue, typename TSize, typename TUsage>
inline void 
deallocate(T & /*me*/,
		   TValue * data, 
		   TSize count,
		   Tag<TUsage> const)
{
	if (data && count)	// .. to use count if SEQAN_PROFILE is not defined
	    SEQAN_PROSUB(SEQAN_PROMEMORY, count * sizeof(TValue));
	operator delete ((void *) data);
}
//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
