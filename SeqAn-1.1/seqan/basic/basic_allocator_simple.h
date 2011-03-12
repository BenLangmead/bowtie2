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
  $Id: basic_allocator_simple.h,v 1.1 2008/08/25 16:20:02 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_ALLOCATOR_SIMPLE_H
#define SEQAN_HEADER_BASIC_ALLOCATOR_SIMPLE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// SimpleAlloc Allocator
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Simple Allocator:
..cat:Allocators
..general:Class.Allocator
..summary:General purpose allocator.
..signature:Allocator< SimpleAlloc<ParentAllocator> >
..param.ParentAllocator:An allocator that is by the simple allocator used to allocate memory.
...default:@Tag.Default@
...remarks:@Tag.Default@ used as allocator means that the default implementations
of @Function.allocate@ and @Function.deallocate@ are used.
*/

template <typename TParentAllocator = Default>
struct SimpleAlloc;

//////////////////////////////////////////////////////////////////////////////


typedef Allocator<SimpleAlloc<Default> > SimpleAllocator;

template <typename TParentAllocator>
struct Allocator<SimpleAlloc<TParentAllocator> >
{
	struct Header
	{
		Header * left;
		Header * right;
		size_t size;
	};

	Header * data_storages;
	Holder<TParentAllocator> data_parent_allocator;

	Allocator():
		data_storages(0)
	{
SEQAN_CHECKPOINT
	}

	Allocator(TParentAllocator & parent_alloc):
		data_storages(0)
	{
SEQAN_CHECKPOINT
		setValue(data_parent_allocator, parent_alloc);
	}

	//Dummy copy
	Allocator(Allocator const &):
		data_storages(0)
	{
	}
	inline Allocator &
	operator = (Allocator const &)
	{
		clear(*this);
		return *this;
	}

	~Allocator()
	{
SEQAN_CHECKPOINT
		clear(*this);
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator>
inline TParentAllocator &
parentAllocator(Allocator<SimpleAlloc<TParentAllocator> > & me)
{
SEQAN_CHECKPOINT
	return value(me.data_parent_allocator);
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.Allocator#clear:
..cat:Memory
..summary:Deallocates all memory blocks.
..signature:clear(allocator)
..param.allocator:Allocator object.
...type:Class.Allocator
...concept:Concept.Allocator
..remarks:This function deallocates all memory blocks 
that was allocated using @Function.allocate@ for $allocator$.
The memory is not pooled but directly passed back to the heap manager.
..see:Function.allocate
..see:Function.deallocate
*/
template <typename TParentAllocator>
void
clear(Allocator<SimpleAlloc<TParentAllocator> > & me)
{
SEQAN_CHECKPOINT
	typedef Allocator<SimpleAlloc<TParentAllocator> > TAllocator;

	while (me.data_storages)
	{
		typename TAllocator::Header * next_storage = me.data_storages->right;
		deallocate(parentAllocator(me), reinterpret_cast<char *>(me.data_storages), me.data_storages->size);
		me.data_storages = next_storage;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void
allocate(Allocator<SimpleAlloc<TParentAllocator> > & me, 
		 TValue * & data,
		 TSize count,
		 Tag<TUsage> const)
{
SEQAN_CHECKPOINT
	typedef Allocator<SimpleAlloc<TParentAllocator> > TAllocator;
	typedef typename TAllocator::Header THeader;

	//compute needed bytes
	size_t bytes_needed = count * sizeof(TValue) + sizeof(THeader);

	//allocate storage from parent
	char * ptr;
	allocate(parentAllocator(me), ptr, bytes_needed, TagAllocateStorage());

	THeader * new_block = reinterpret_cast<THeader *>(ptr);
	new_block->left = 0;
	new_block->right = me.data_storages;
	new_block->size = bytes_needed;

	if (me.data_storages)
	{
		me.data_storages->left = new_block;
	}
	me.data_storages = new_block;

	//return data
	data = reinterpret_cast<TValue *>(ptr + sizeof(THeader));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void 
deallocate(Allocator<SimpleAlloc<TParentAllocator> > & me,
		   TValue * data, 
		   TSize,
		   Tag<TUsage> const)
{
SEQAN_CHECKPOINT
	typedef Allocator<SimpleAlloc<TParentAllocator> > TAllocator;
	typedef typename TAllocator::Header THeader;

	//update links
	THeader & header = *(reinterpret_cast<THeader *>(data) - 1);
	if (header.left)
	{
		header.left->right = header.right;
	}
	else
	{
		me.data_storages = header.right;
	}
	if (header.right)
	{
		header.right->left = header.left;
	}

	//deallocate storage using parent
	char * ptr = reinterpret_cast<char *>(& header);
	deallocate(parentAllocator(me), ptr, header.size);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
