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
  $Id: basic_allocator_singlepool.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_ALLOCATOR_SINGLE_POOL_H
#define SEQAN_HEADER_BASIC_ALLOCATOR_SINGLE_POOL_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// SinglePool Allocator
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Single Pool Allocator:
..cat:Allocators
..general:Class.Allocator
..summary:Allocator that pools memory blocks of specific size.
..signature:Allocator< SinglePool<SIZE, ParentAllocator> >
..param.SIZE:Size of memory blocks that are pooled.
...value:An unsigned integer with $SIZE >= sizeof(void *)$.
..param.ParentAllocator:An allocator that is by the pool allocator used to allocate memory.
...default:@Spec.Simple Allocator@
...note:The single pool allocator only supports @Function.clear@ if this function is also implemented for $ParentAllocator$.
..remarks:A pool allocator allocates several memory blocks at once. 
Freed blocks are not immediately deallocated but recycled in subsequential allocations.
This way, the number of calls to the heap manager is reduced, and that speeds up memory management.
...text:The single pool allocator only pools memory blocks of size $SIZE$.
Blocks of other sizes are allocated and deallocated using an allocator of type $ParentAllocator$.
...text:Using the single pool allocator for blocksizes larger than some KB is not advised.
*/

template <size_t SIZE, typename TParentAllocator = SimpleAllocator>
struct SinglePool;

//////////////////////////////////////////////////////////////////////////////

template <size_t SIZE, typename TParentAllocator>
struct Allocator<SinglePool<SIZE, TParentAllocator> >
{
	enum
	{
		SIZE_PER_ITEM = SIZE,
		ITEMS_PER_BLOCK = (SIZE_PER_ITEM < 0x0100) ? 0x01000 / SIZE_PER_ITEM : 16,
		STORAGE_SIZE = SIZE * ITEMS_PER_BLOCK,

		STORAGE_SIZE_MIN = SIZE
	};

	char * data_recycled_blocks;
	char * data_current_begin;
	char * data_current_end;
	char * data_current_free;
	Holder<TParentAllocator> data_parent_allocator;

	Allocator()
	{
SEQAN_CHECKPOINT
		data_recycled_blocks = data_current_end = data_current_free = 0;
		//dont need to initialize data_current_begin
	}

	Allocator(size_t reserve_item_count)
	{
SEQAN_CHECKPOINT
		data_recycled_blocks = 0;

		size_t storage_size = (reserve_item_count * SIZE > STORAGE_SIZE_MIN) ? reserve_item_count * SIZE : STORAGE_SIZE_MIN;
		allocate( parentAllocator( *this ), data_current_begin, storage_size );
		data_current_end = data_current_begin + storage_size;
		data_current_free = data_current_begin;
	}

	Allocator(TParentAllocator & parent_alloc)
	{
SEQAN_CHECKPOINT
		setValue(data_parent_allocator, parent_alloc);

		data_recycled_blocks = data_current_end = data_current_free = 0;
		//dont need to initialize data_current_begin
	}

	Allocator(size_t reserve_item_count, TParentAllocator & parent_alloc)
	{
SEQAN_CHECKPOINT
		data_recycled_blocks = 0;

		setValue(data_parent_allocator, parent_alloc);

		size_t storage_size = (reserve_item_count * SIZE > STORAGE_SIZE_MIN) ? reserve_item_count * SIZE : STORAGE_SIZE_MIN;
		allocate( parentAllocator( *this ), data_current_begin, storage_size );
		data_current_end = data_current_begin + storage_size;
		data_current_free = data_current_begin;
	}

	//Dummy copy
	Allocator(Allocator const &)
	{
		data_recycled_blocks = data_current_end = data_current_free = 0;
		//dont need to initialize data_current_begin
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

template <size_t SIZE, typename TParentAllocator>
inline TParentAllocator &
parentAllocator(Allocator<SinglePool<SIZE, TParentAllocator> > & me)
{
SEQAN_CHECKPOINT
	return value(me.data_parent_allocator);
}

//////////////////////////////////////////////////////////////////////////////

template <size_t SIZE, typename TParentAllocator>
void
clear(Allocator<SinglePool<SIZE, TParentAllocator> > & me)
{
SEQAN_CHECKPOINT

	me.data_recycled_blocks = me.data_current_end = me.data_current_free = 0;

	clear(parentAllocator(me));
}

//////////////////////////////////////////////////////////////////////////////

template <size_t SIZE, typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void
allocate(Allocator<SinglePool<SIZE, TParentAllocator> > & me, 
		 TValue * & data,
		 TSize count,
		 Tag<TUsage> const tag_)
{
SEQAN_CHECKPOINT
	typedef Allocator<SinglePool<SIZE, TParentAllocator> > TAllocator;
	size_t bytes_needed = count * sizeof(TValue);

	if (bytes_needed != TAllocator::SIZE_PER_ITEM)
	{//no blocking
		allocate(parentAllocator(me), data, count, tag_);
		return;
	}

	char * ptr;
	if (me.data_recycled_blocks)
	{//use recycled
		ptr = me.data_recycled_blocks;
		me.data_recycled_blocks = * reinterpret_cast<char **>(ptr);
	}
	else
	{//use new
		ptr = me.data_current_free;
		if (ptr + bytes_needed > me.data_current_end)
		{//not enough free space in current storage: allocate new
			allocate(parentAllocator(me), ptr, (size_t) TAllocator::STORAGE_SIZE, tag_);
			me.data_current_begin = ptr;
			me.data_current_end = ptr + TAllocator::STORAGE_SIZE;
		}
		me.data_current_free = ptr + bytes_needed;
	}

	data = reinterpret_cast<TValue *>(ptr);
}

//////////////////////////////////////////////////////////////////////////////

template <size_t SIZE, typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void 
deallocate(Allocator<SinglePool<SIZE, TParentAllocator> > & me,
		   TValue * data, 
		   TSize count,
		   Tag<TUsage> const tag_)
{
SEQAN_CHECKPOINT
	typedef Allocator<SinglePool<SIZE, TParentAllocator> > TAllocator;

	size_t bytes_needed = count * sizeof(TValue);

	if (bytes_needed != TAllocator::SIZE_PER_ITEM)
	{//no blocking
		deallocate(parentAllocator(me), data, count, tag_);
		return;
	}

	//link in recycling list
	*reinterpret_cast<char **>(data) = me.data_recycled_blocks;
	me.data_recycled_blocks = reinterpret_cast<char *>(data);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// alternative Interface that takes a Type instead of a SIZE
//////////////////////////////////////////////////////////////////////////////


template <typename TValue, typename TParentAllocator = SimpleAllocator>
struct SinglePool2;

template <typename TValue, typename TParentAllocator>
struct Allocator<SinglePool2<TValue, TParentAllocator> >
{
	Allocator<SinglePool<sizeof(TValue), TParentAllocator> > data_alloc;


	Allocator(size_t reserve_item_count)
		: data_alloc(reserve_item_count)
	{
	}

	Allocator(TParentAllocator & parent_alloc)
		: data_alloc(parent_alloc)
	{
	}

	Allocator(size_t reserve_item_count, TParentAllocator & parent_alloc)
		: data_alloc(reserve_item_count, parent_alloc)

	{
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TParentAllocator>
inline TParentAllocator &
parentAllocator(Allocator<SinglePool2<TValue, TParentAllocator> > & me)
{
SEQAN_CHECKPOINT
	return parentAllocator(me.data_alloc);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TParentAllocator>
void
clear(Allocator<SinglePool2<TValue, TParentAllocator> > & me)
{
SEQAN_CHECKPOINT
	clear(me.data_alloc);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TParentAllocator, typename TValue2, typename TSize, typename TUsage>
inline void
allocate(Allocator<SinglePool2<TValue, TParentAllocator> > & me, 
		 TValue2 * & data,
		 TSize count,
		 Tag<TUsage> const tag_)
{
SEQAN_CHECKPOINT
	allocate(me.data_alloc, data, count, tag_);
}

template <typename TValue, typename TParentAllocator, typename TValue2, typename TSize, typename TUsage>
inline void 
deallocate(Allocator<SinglePool2<TValue, TParentAllocator> > & me,
		   TValue2 * data, 
		   TSize count,
		   Tag<TUsage> const tag_)
{
SEQAN_CHECKPOINT
	deallocate(me.data_alloc, data, count, tag_);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
