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
  $Id: string_stack.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_STACK_H
#define SEQAN_HEADER_GRAPH_STACK_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// TAGS
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Block String:
..cat:Strings
..general:Class.String
..summary:String optimized for push_back, top, and pop (Stack behaviour).
..signature:String<TValue, Block<size> >
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.size:A positive integer that specifies the number of values in each allocated block.
...remarks: Size should be a power of 2, e.g., 1024.
*/

	template<unsigned int SPACE = 4096>
	struct Block;

//////////////////////////////////////////////////////////////////////////////


	template<typename TValue, unsigned int SPACE>
	class String<TValue, Block<SPACE> > 
	{
		typedef String<TValue, Array<SPACE> >				TBlock;
		typedef TBlock*										PBlock;
		typedef Allocator< SinglePool<sizeof(TBlock)> >		TAllocator;

	public:
		typedef typename Iterator<TBlock, Standard>::Type	TBlockIter;
		typedef String<PBlock>								TBlockTable;

		TBlockTable		blocks;
		TBlockIter		blockFirst, blockLast;	// current block boundaries
		TBlockIter		lastValue;				// pointer to top value
		TAllocator		alloc;
	    
	//____________________________________________________________________________
	      
		public:
			String():
				blockFirst(TBlockIter()),
				blockLast(TBlockIter()),
				lastValue(TBlockIter()) {}

			template<typename TSource>
			String(TSource const& source):
				blockFirst(TBlockIter()),
				blockLast(TBlockIter()),
				lastValue(TBlockIter())
			{
			SEQAN_CHECKPOINT
				assign(*this, source);
			} 

			String(String const & source):
				blockFirst(TBlockIter()),
				blockLast(TBlockIter()),
				lastValue(TBlockIter())
			{
			SEQAN_CHECKPOINT
				assign(*this, source);
			}

			template<typename TSource>
			String & operator =(TSource const& source) 
			{
			SEQAN_CHECKPOINT
				assign(*this, source);
				return *this;
			}

			String & operator =(String const& _other)	
			{
			SEQAN_CHECKPOINT
				if (this == &_other) return *this;
				assign(*this, _other);
				return *this;
			}

			~String() 
			{
				clear(*this);
			}

	//____________________________________________________________________________

		public:
			template<typename TPos>
			inline typename Reference<String>::Type 
				operator[] (TPos pos) 
			{
			SEQAN_CHECKPOINT
				return value(*this, pos);
			}

			template<typename TPos>
			inline typename Reference<String const>::Type 
				operator[] (TPos pos) const 
			{
			SEQAN_CHECKPOINT
				return value(*this, pos);
			}
	};


	template<typename TValue, unsigned int SPACE>
	struct DefaultOverflowImplicit< String<TValue, Block<SPACE> > >
	{
		typedef Generous Type;
	};


//////////////////////////////////////////////////////////////////////////////
// Block metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Iterators
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Spec.Block String

	template<typename TValue, unsigned int SPACE>
	struct Iterator<String<TValue, Block<SPACE> >, Standard> 
	{
		typedef Iter<String<TValue, Block<SPACE> >, PositionIterator> Type;
	};

	template<typename TValue, unsigned int SPACE>
	struct Iterator<String<TValue, Block<SPACE> > const, Standard> 
	{
		typedef Iter<String<TValue, Block<SPACE> > const, PositionIterator> Type;
	};

	template<typename TValue, unsigned int SPACE>
	struct Iterator<String<TValue, Block<SPACE> >, Rooted> 
	{
		typedef Iter<String<TValue, Block<SPACE> >, PositionIterator> Type;
	};

	template<typename TValue, unsigned int SPACE>
	struct Iterator<String<TValue, Block<SPACE> > const, Rooted> 
	{
		typedef Iter<String<TValue, Block<SPACE> > const, PositionIterator> Type;
	};


///////////////////////////////////////////////////////////////
// Block interface
///////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// begin
//////////////////////////////////////////////////////////////////////////////

	template<typename TValue, unsigned int SPACE, typename TSpec>
	inline typename Iterator<String<TValue, Block<SPACE> >, Tag<TSpec> const >::Type 
	begin(String<TValue, Block<SPACE> > &me, Tag<TSpec> const)
	{
	SEQAN_CHECKPOINT
		return Iter<String<TValue, Block<SPACE> >, PositionIterator>(me, 0);
	}

	template<typename TValue, unsigned int SPACE, typename TSpec>
	inline typename Iterator<String<TValue, Block<SPACE> > const, Tag<TSpec> const>::Type 
	begin(String<TValue, Block<SPACE> > const &me, Tag<TSpec> const)
	{
	SEQAN_CHECKPOINT
		return Iter<String<TValue, Block<SPACE> > const, PositionIterator>(me, 0);
	}


	template<typename TValue, unsigned int SPACE, typename TSpec>
	inline typename Iterator<String<TValue, Block<SPACE> >, Tag<TSpec> const >::Type 
	end(String<TValue, Block<SPACE> > &me, Tag<TSpec> const)
	{
	SEQAN_CHECKPOINT
		return Iter<String<TValue, Block<SPACE> >, PositionIterator>(me, length(me));
	}

	template<typename TValue, unsigned int SPACE, typename TSpec>
	inline typename Iterator<String<TValue, Block<SPACE> > const, Tag<TSpec> const>::Type 
	end(String<TValue, Block<SPACE> > const &me, Tag<TSpec> const)
	{
	SEQAN_CHECKPOINT
		return Iter<String<TValue, Block<SPACE> > const, PositionIterator>(me, length(me));
	}

	template<typename TValue, unsigned int SPACE, typename TSource>
	inline void 
	assign(
		String<TValue, Block<SPACE> >& target, 
		TSource const& source) 
	{
	SEQAN_CHECKPOINT
		clear(target);
		typedef typename Iterator<TSource const, Standard>::Type TIter;
		for(TIter it = begin(source, Standard()); !atEnd(it, source); goNext(it))
			push(target, *it);
	}


	template<typename TValue, unsigned int SPACE, typename TPos>
	inline typename Reference<String<TValue, Block<SPACE> > >::Type 
	value(
		String<TValue, Block<SPACE> >& stack, 
		TPos const pos) 
	{
	SEQAN_CHECKPOINT
		return value(*(stack.blocks[pos / SPACE]), pos % SPACE);
	}

	template<typename TValue, unsigned int SPACE, typename TPos>
	inline typename Reference<String<TValue, Block<SPACE> > >::Type 
	value(
		String<TValue, Block<SPACE> > const& stack, 
		TPos const pos) 
	{
	SEQAN_CHECKPOINT
		return value(*(stack.blocks[pos / SPACE]), pos % SPACE);
	}

	template<typename TValue, unsigned int SPACE, typename TIteratorSpec>
	inline bool 
	atEnd(
		Iter<String<TValue, Block<SPACE> >, TIteratorSpec>& it, 
		String<TValue, Block<SPACE> >& container) 
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<String<TValue, Block<SPACE> >, Standard>::Type TIter;
		TIter endIt = end(container, Standard());
		return (it == endIt);
	}


	template<typename TValue, unsigned int SPACE>
	inline void 
	clear(String<TValue, Block<SPACE> >& me)
	{
	SEQAN_CHECKPOINT
		typedef String<TValue, Block<SPACE>	>			TBlockString;
		typedef typename TBlockString::TBlockTable		TBlockTable;
		typedef typename Iterator<TBlockTable, Standard>::Type	TIter;
		
		TIter it = begin(me.blocks), itEnd = end(me.blocks);
		while (it != itEnd) {
			deallocate(me.alloc, *it, 1);
			++it;
		}
		clear(me.blocks);
		me.lastValue = me.blockLast = typename TBlockString::TBlockIter();
	}


//////////////////////////////////////////////////////////////////////////////
///.Function.reserve.param.object.type:Spec.Block String
/*
	template <typename TValue, unsigned int SPACE, typename TSize, typename TExpand>
	inline typename Size< String<TValue, Block<SPACE> > >::Type
	reserve(
		String<TValue, Block<SPACE> >& me, 
		TSize new_capacity,
		Tag<TExpand> const tag)
	{
	SEQAN_CHECKPOINT
		reserve(me.blocks, (new_capacity + SPACE - 1) / SPACE, tag);
		return capacity(me.blocks) * SPACE;
	}
*/

	template<typename TValue, unsigned int SPACE, typename TSize2, typename TExpand>
	inline typename Size< String<TValue, Block<SPACE> > >::Type
	resize(String<TValue, Block<SPACE> > & me,
		TSize2 new_length,
		Tag<TExpand> const)
	{
	SEQAN_CHECKPOINT
		typedef String<TValue, Block<SPACE>	>			TBlockString;
		typedef typename Size<TBlockString>::Type		TSize;
		TSize len = length(me);

		if (new_length > len)
		{
			for (; len < new_length; ++len) push(me);
		}
		else if (new_length < len)
		{
			for (; len > new_length; --len) pop(me);
		}
		return new_length;
	}
	template<typename TValue, unsigned int SPACE, typename TSize2>
	inline typename Size< String<TValue, Block<SPACE> > >::Type
	resize(String<TValue, Block<SPACE> > & me,
		TSize2 new_length,
		Limit)
	{
	SEQAN_CHECKPOINT
		typedef String<TValue, Block<SPACE>	>			TBlockString;
		typedef typename Size<TBlockString>::Type		TSize;
		TSize len = length(me);

		if (new_length > capacity(me)) new_length = capacity(me);

		if (new_length > len)
		{
			TValue val;
			for (; len < new_length; ++len) push(me, val);
		}
		else if (new_length < len)
		{
			for (; len > new_length; --len) pop(me);
		}
		return new_length;
	}

	
	//dummy implementation
	template<typename TValue, unsigned int SPACE, typename TSize, typename TExpand>
	inline typename Size< String<TValue, Block<SPACE> > >::Type
	reserve(String<TValue, Block<SPACE> > & me,
		TSize new_capacity,
		Tag<TExpand> const)
	{
	SEQAN_CHECKPOINT
		return new_capacity;
	}


	template<typename TValue, unsigned int SPACE, typename TSource, typename TExpand>
	inline void 
	append(
		String<TValue, Block<SPACE> >& me,
		TSource const& source,
		Tag<TExpand> const tag)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TSource const, Standard>::Type TIter;
		for(TIter it = begin(source, Standard()); !atEnd(it, source); goNext(it))
			appendValue(me, *it);
	}

//////////////////////////////////////////////////////////////////////////////
///.Function.appendValue.param.target.type:Spec.Block String

	template<typename TValue, unsigned int SPACE, typename TVal, typename TExpand>
	inline void 
	appendValue(
		String<TValue, Block<SPACE> >& me, 
		TVal const& source,
		Tag<TExpand> const tag)
	{
	SEQAN_CHECKPOINT
		if (me.lastValue == me.blockLast) {
			typename Size< String<TValue, Block<SPACE> > >::Type last = length(me.blocks);

			resize(me.blocks, last + 1, tag);
			allocate(me.alloc, me.blocks[last], 1);
			me.lastValue = me.blockFirst = begin(*me.blocks[last]);
			me.blockLast = (me.blockFirst + (SPACE - 1));
		} else
			++me.lastValue;
		valueConstruct(me.lastValue, source);
	}
	 
	template<typename TValue, unsigned int SPACE, typename TVal>
	inline void 
	push(
		String<TValue, Block<SPACE> >& me, 
		TVal const& source)
	{
		appendValue(me, source);
	}

	template<typename TValue, unsigned int SPACE>
	inline void 
	push(String<TValue, Block<SPACE> >& me)
	{
	SEQAN_CHECKPOINT
		if (me.lastValue == me.blockLast) {
			typename Size< String<TValue, Block<SPACE> > >::Type last = length(me.blocks);

			resize(me.blocks, last + 1, typename DefaultOverflowImplicit<String<TValue, Block<SPACE> > >::Type());
			allocate(me.alloc, me.blocks[last], 1);
			me.lastValue = me.blockFirst = begin(*me.blocks[last]);
			me.blockLast = (me.blockFirst + (SPACE - 1));
		} else
			++me.lastValue;
		valueConstruct(me.lastValue);
	}
	 
	template<typename TValue, unsigned int SPACE, typename TVal>
	inline void 
	push_back(
		String<TValue, Block<SPACE> >& me, 
		TVal const& source)
	{
		appendValue(me, source);
	}

	template<typename TValue, unsigned int SPACE>
	inline TValue &
	top(String<TValue, Block<SPACE> > & me) 
	{
	SEQAN_CHECKPOINT
		return *me.lastValue;
	}

	template<typename TValue, unsigned int SPACE>
	inline TValue const &
	top(String<TValue, Block<SPACE> > const& me) 
	{
	SEQAN_CHECKPOINT
		return *me.lastValue;
	}

	template<typename TValue, unsigned int SPACE>
	inline TValue &
	topPrev(String<TValue, Block<SPACE> > & me) 
	{
	SEQAN_CHECKPOINT
		if (me.lastValue != me.blockFirst)
			return *(me.lastValue - 1);
		else
			return *(begin(*me.blocks[length(me.blocks) - 1]) + (SPACE - 1));
	}

	template<typename TValue, unsigned int SPACE>
	inline TValue const &
	topPrev(String<TValue, Block<SPACE> > const& me) 
	{
	SEQAN_CHECKPOINT
		if (me.lastValue != me.blockFirst)
			return *(me.lastValue - 1);
		else
			return *(begin(*me.blocks[length(me.blocks) - 1]) + (SPACE - 1));
	}

	template<typename TValue, unsigned int SPACE>
	inline void 
	pop(String<TValue, Block<SPACE> >& me) 
	{
	SEQAN_CHECKPOINT
		if (me.lastValue == me.blockFirst) {
			typename Size< String<TValue, Block<SPACE> > >::Type last = length(me.blocks);

			if (last) {
				valueDestruct(me.lastValue);
				deallocate(me.alloc, me.blocks[--last], 1);
				resize(me.blocks, last);
				if (last) {
					me.blockFirst = begin(*me.blocks[--last]);
					me.lastValue = me.blockLast = (me.blockFirst + (SPACE - 1));
				}
			}
		} else {
			valueDestruct(me.lastValue);
			--me.lastValue;
		}
	}

	template<typename TValue, unsigned int SPACE>
	inline void 
	pop_back(String<TValue, Block<SPACE> >& me) {
		pop(me);
	}

	template<typename TValue, unsigned int SPACE>
	inline bool 
	empty(String<TValue, Block<SPACE> > const& me) 
	{
	SEQAN_CHECKPOINT
		return length(me.blocks) == 0;
	}

	template<typename TValue, unsigned int SPACE>
	inline typename Size<String<TValue, Block<SPACE> > >::Type
	length(String<TValue, Block<SPACE> > const & me) 
	{
	SEQAN_CHECKPOINT
		if (length(me.blocks))
			return (length(me.blocks) - 1) * SPACE + (me.lastValue - me.blockFirst) + 1;
		else
			return 0;
	}

	template<typename TValue, unsigned int SPACE>
	inline typename Size<String<TValue, Block<SPACE> > >::Type
	capacity(String<TValue, Block<SPACE> > const & me) 
	{
	SEQAN_CHECKPOINT
		if (length(me.blocks))
			return length(me.blocks) * SPACE;
		else
			return 0;
	}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
