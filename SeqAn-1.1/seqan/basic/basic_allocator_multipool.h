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
  $Id: basic_allocator_multipool.h,v 1.2 2009/02/19 01:51:23 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_ALLOCATOR_MULTIPOOL_H
#define SEQAN_HEADER_BASIC_ALLOCATOR_MULTIPOOL_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// MultiPool Allocator
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Multi Pool Allocator:
..cat:Allocators
..general:Class.Allocator
..summary:Allocator that pools memory blocks.
..signature:Allocator< MultiPool<ParentAllocator, BLOCKING_LIMIT> >
..param.ParentAllocator:An allocator that is by the pool allocator used to allocate memory.
...default:@Spec.Simple Allocator@
...note:The multi pool allocator only supports @Function.clear@ if this function is also implemented for $ParentAllocator$.
..remarks:A pool allocator allocates several memory blocks at once.
..param.BLOCKING_LIMIT:The maximum size for memory blocks to be pooled.
...default:256
Freed blocks are not immediately deallocated but recycled in subsequential allocations.
This way, the number of calls to the heap manager is reduced, and that speeds up memory management.
...text:Note that memory blocks larger than $BLOCKING_LIMIT$ are not pooled
but immediately allocated and deallocated using $ParentAllocator$.
*/


template <typename TParentAllocator = Allocator<SimpleAlloc<Default> >, unsigned int BLOCKING_LIMIT = 0x100>
struct MultiPool;

//////////////////////////////////////////////////////////////////////////////

typedef Allocator<MultiPool<Allocator<SimpleAlloc<Default> >, 0x100> > PoolAllocator;

template <typename TParentAllocator, unsigned int BLOCKING_LIMIT_>
struct Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT_> >
{
	enum
	{
		BLOCKING_LIMIT = BLOCKING_LIMIT_,
		GRANULARITY_BITS = 2,
		BLOCKING_COUNT = BLOCKING_LIMIT >> GRANULARITY_BITS,
		STORAGE_SIZE = 0xf80
	};

	char * data_recycled_blocks [BLOCKING_COUNT];
	char * data_current_begin [BLOCKING_COUNT];
	char * data_current_free [BLOCKING_COUNT];
	Holder<TParentAllocator> data_parent_allocator;

	Allocator()
	{
SEQAN_CHECKPOINT
		memset(data_recycled_blocks, 0, sizeof(data_recycled_blocks));
		memset(data_current_begin, 0, sizeof(data_current_begin));
		memset(data_current_free, 0, sizeof(data_current_free));
	}

	Allocator(TParentAllocator & parent_alloc)
	{
SEQAN_CHECKPOINT
		memset(data_recycled_blocks, 0, sizeof(data_recycled_blocks));
		memset(data_current_begin, 0, sizeof(data_current_begin));
		memset(data_current_free, 0, sizeof(data_current_free));

		setValue(data_parent_allocator, parent_alloc);
	}

	//Dummy copy
	Allocator(Allocator const &)
	{
		memset(data_recycled_blocks, 0, sizeof(data_recycled_blocks));
		memset(data_current_begin, 0, sizeof(data_current_begin));
		memset(data_current_free, 0, sizeof(data_current_free));
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

template <typename TParentAllocator, unsigned int BLOCKING_LIMIT>
inline TParentAllocator &
parentAllocator(Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > & me)
{
SEQAN_CHECKPOINT
	return value(me.data_parent_allocator);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator, unsigned int BLOCKING_LIMIT>
void
clear(Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > & me)
{
SEQAN_CHECKPOINT
	memset(me.data_recycled_blocks, 0, sizeof(me.data_recycled_blocks));
	memset(me.data_current_begin, 0, sizeof(me.data_current_begin));
	memset(me.data_current_free, 0, sizeof(me.data_current_free));

	clear(parentAllocator(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator, unsigned int BLOCKING_LIMIT>
inline unsigned int
_allocatorBlockNumber(Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > &,
					  size_t size_)
{
SEQAN_CHECKPOINT
	typedef Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > TAllocator;

	SEQAN_ASSERT(size_)

	if (size_ < BLOCKING_LIMIT)
	{//blocks
		return size_ >> TAllocator::GRANULARITY_BITS;
	}
	else
	{//no blocking
		return TAllocator::BLOCKING_COUNT;
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator, unsigned int BLOCKING_LIMIT, typename TValue, typename TSize, typename TUsage>
inline void
allocate(Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > & me,
		 TValue * & data,
		 TSize count,
		 Tag<TUsage> const tag_)
{
SEQAN_CHECKPOINT
	typedef Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > TAllocator;

	size_t bytes_needed = count * sizeof(TValue);
	char * ptr;

	unsigned int block_number =  _allocatorBlockNumber(me, bytes_needed);
	if (block_number == TAllocator::BLOCKING_COUNT)
	{//no blocking
		return allocate(parentAllocator(me), data, count, tag_);
	}

	bytes_needed = (block_number + 1) << TAllocator::GRANULARITY_BITS;

	if (me.data_recycled_blocks[block_number])
	{//use recycled
		ptr = me.data_recycled_blocks[block_number];
		me.data_recycled_blocks[block_number] = * reinterpret_cast<char **>(ptr);
	}
	else
	{//use new
		ptr = me.data_current_free[block_number];
		if (!ptr || (ptr + bytes_needed > me.data_current_begin[block_number] + TAllocator::STORAGE_SIZE))
		{//not enough free space in current storage: allocate new
			allocate(parentAllocator(me), ptr, (size_t) TAllocator::STORAGE_SIZE, tag_);
			me.data_current_begin[block_number] = ptr;
		}
		me.data_current_free[block_number] = ptr + bytes_needed;
	}

	data = reinterpret_cast<TValue *>(ptr);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator, unsigned int BLOCKING_LIMIT, typename TValue, typename TSize, typename TUsage>
inline void
deallocate(Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > & me,
		   TValue * data,
		   TSize count,
		   Tag<TUsage> const tag_)
{
SEQAN_CHECKPOINT
	typedef Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > TAllocator;

	size_t bytes_needed = count * sizeof(TValue);

	unsigned int block_number = _allocatorBlockNumber(me, bytes_needed);
	if (block_number == TAllocator::BLOCKING_COUNT)
	{//no blocking
		return deallocate(parentAllocator(me), data, count, tag_);
	}

	bytes_needed = (block_number + 1) << TAllocator::GRANULARITY_BITS;

	//link in recycling list
	*reinterpret_cast<char **>(data) = me.data_recycled_blocks[block_number];
	me.data_recycled_blocks[block_number] = reinterpret_cast<char *>(data);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
