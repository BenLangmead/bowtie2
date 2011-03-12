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
  $Id: basic_allocator_to_std.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_ALLOCATOR_TO_STD_H
#define SEQAN_HEADER_BASIC_ALLOCATOR_TO_STD_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//helper caller for calling functions that have same name as member functions

template <typename TMe, typename TValue, typename TSize>
inline void call_allocate(TMe & me, TValue * & data, TSize const count)
{
	allocate(me, data, count);
}
template <typename TMe, typename TValue, typename TSize>
inline void call_deallocate(TMe & me, TValue * data, TSize const count)
{
	deallocate(me, data, count);
}

//////////////////////////////////////////////////////////////////////////////
//Filter that adapts seqan allocator zu std allocator
/**
.Class.ToStdAllocator:
..summary:Emulates standard conform allocator.
..signature:ToStdAllocator<THost, TValue>
..param.THost:Type of the host allocator object.
...text:This object is used to call @Function.allocate@ and @Function.deallocate@.
..param.TValue:Type of allocated items.
..remarks:The member functions $allocate$ and $deallocate$ of $ToStdAllocator$ call
the (globale) functions @Function.allocate@ and @Function.deallocate@, respectively. The globale functions
get an allocator object as their first arguments. This allocator object is not the $ToStdAllocator$ object itself,
but the host object that was given to the constructor. 
..remarks:
..see:Function.allocate
..see:Function.deallocate
*/
template <typename THost, typename TValue>
struct ToStdAllocator
{
	typedef TValue value_type;
	typedef value_type * pointer;
	typedef value_type & reference;
	typedef value_type const * const_pointer;
	typedef value_type const & const_reference;

//	typedef typename THost::Size size_type;
//	typedef typename THost::Difference difference_type;
	typedef size_t size_type;
	typedef ptrdiff_t difference_type;

/**
.Memfunc.ToStdAllocator:
..summary:Constructor
..signature:ToStdAllocator(host)
..class:Class.ToStdAllocator
..param.host:The host object that is used as allocator for @Function.allocate@ and @Function.deallocate@.
*/
	ToStdAllocator(THost & host): m_host(& host)
	{
	}
	ToStdAllocator(ToStdAllocator const & alloc): m_host(alloc.m_host)
	{
	}
	ToStdAllocator & operator= (ToStdAllocator const & alloc)
	{
		m_host = alloc.m_host;
	}
	~ToStdAllocator()
	{
	}

/**
.Function.host:
..summary:The object a given object depends on.
..cat:Dependent Objects
..signature:host(object)
..param.object:An object.
...type:Class.ToStdAllocator
..returns:The host object.
*/
    friend THost & host(ToStdAllocator & me)
    {
        return *me.m_host;
    }

	pointer allocate(size_type count)
	{
		value_type * ptr;
		call_allocate(*m_host, ptr, count);
		return pointer(ptr);
	}
	pointer allocate(size_type count, const void *)
	{
		value_type * ptr;
		call_allocate(*m_host, ptr, count);
		return pointer(ptr);
	}

	void deallocate(pointer data, size_type count)
	{
		call_deallocate(*m_host, data, count);
	}

	void construct(pointer ptr, const_reference data)
	{
		new(ptr) TValue(data);
	}

	void destroy(pointer ptr)
	{
		ptr->~TValue();
	}

	pointer address(reference value) const
	{
		return (&value);
	}
	const_pointer address(const_reference value) const
	{
		return (&value);
	}

	size_type max_size() const
	{
		return ~0UL / sizeof(value_type);
	}

	template<class TValue2>
	struct rebind
	{
		typedef ToStdAllocator<THost, TValue2> other;
	};

	private:
		THost * m_host;
};
//////////////////////////////////////////////////////////////////////////////



//returns std-allocator type (for allocators)
template <typename T, typename TData>
struct StdAllocator
{
	typedef ToStdAllocator<T, TData> Type;
};


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN


//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
