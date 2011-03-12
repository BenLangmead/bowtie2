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
  $Id: chunk_collector.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_CHUNK_COLLECTOR_H
#define SEQAN_HEADER_CHUNK_COLLECTOR_H


namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////

/**
.Internal._ChunkCollector:
..cat:Classes
..summary:Reads piecewise from stream, collects pieces (chunks) in a vector.
..signature:_ChunkCollector<Host>
..param.Host:Type of host object that is used as allocator.
*/

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct ChunkLength
{
	enum { VALUE = 1024 };
};

//////////////////////////////////////////////////////////////////////////////
// _StreamChunkCollector class: collects content of a stream in chunks
//////////////////////////////////////////////////////////////////////////////

template <typename THost>
class _ChunkCollector
{
protected:
	THost * data_host;
	typename Size<THost>::Type data_length;

	typedef ::std::vector<typename Value<THost>::Type *, ToStdAllocator<THost, typename Value<THost>::Type *> > Chunk_Holder;
	Chunk_Holder data_chunks; 

public:
	static int const CHUNK_LENGTH = ChunkLength<_ChunkCollector>::VALUE;

public:
	_ChunkCollector(THost & _host):
		data_host(& _host),
		data_length(0),
		data_chunks(typename Chunk_Holder::allocator_type(_host))
	{
	}

	~_ChunkCollector()
	{
		clear(*this);
	}

public:

	friend inline void
	clear(_ChunkCollector & me)
	{
		typename Chunk_Holder::iterator it = me.data_chunks.begin();
		typename Chunk_Holder::iterator it_end = me.data_chunks.end();

		for (; it != it_end; ++it)
		{
			deallocate(me.data_host, *it, CHUNK_LENGTH);
		}

		me.data_chunks.clear();
		me.data_length = 0;
	}

	friend inline typename Size<THost>::Type
	length(_ChunkCollector const & me)
	{
		return me.data_length;
	}

	friend inline void
	_setLength(_ChunkCollector & me, typename Size<THost>::Type new_length)
	{
		me.data_length = new_length;
	}

	friend inline int
	chunkCount(_ChunkCollector const & me)
	{
		return me.data_chunks.size();
	}

	friend inline typename Value<THost>::Type *
	getChunk(_ChunkCollector const & me, int chunk_number)
	{
		return me.data_chunks[chunk_number];
	}

	friend inline typename Value<THost>::Type *
	createChunk(_ChunkCollector & me)
	{
		typename Value<THost>::Type * new_chunk;
		allocate(me.data_host, new_chunk, CHUNK_LENGTH);
		me.data_chunks.push_back(new_chunk);
		return new_chunk;
	}
};


//////////////////////////////////////////////////////////////////////////////

template <typename THost>
struct Host<_ChunkCollector<THost> >
{
	typedef THost Type;
};

template <typename THost>
struct Host<_ChunkCollector<THost> const >
{
	typedef THost Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
struct Value<_ChunkCollector<THost> >
{
	typedef typename Value<THost>::Type Type;
};

template <typename THost>
struct Value<_ChunkCollector<THost> const >
{
	typedef typename Value<THost>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
struct GetValue<_ChunkCollector<THost> >
{
	typedef typename GetValue<THost>::Type Type;
};

template <typename THost>
struct GetValue<_ChunkCollector<THost> const >
{
	typedef typename GetValue<THost>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
struct Size<_ChunkCollector<THost> >
{
	typedef typename Size<THost>::Type Type;
};

template <typename THost>
struct Size<_ChunkCollector<THost> const >
{
	typedef typename Size<THost>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

struct _Assign_Stream_2_ChunkCollector
{
	template <typename THost, typename TSource>
	static inline void 
	assign_(_ChunkCollector<THost> & target,
		TSource & source)
	{
		clear(target);

		while (!_streamEOF(source))
		{
			typename Value<THost>::Type * chunk = createChunk(target);
			typename Size<THost>::Type count = _streamRead(chunk, source, ChunkLength< _ChunkCollector<THost> >::VALUE);
			_setLength(target, length(target) + count);
		}
	}

	template <typename THost, typename TSource>
	static inline void 
	assign_(_ChunkCollector<THost> & target,
		TSource & source,
		typename Size< _ChunkCollector<THost> >::Type limit)
	{
		clear(target);

		while (!_streamEOF(source))
		{
			typename Value<THost>::Type * chunk = createChunk(target);
			typename Size<THost>::Type count = _streamRead(chunk, source, ChunkLength< _ChunkCollector<THost> >::VALUE);
			_setLength(target, length(target) + count);

			if (length(target) >= limit)
			{
				_setLength(target, limit);
				break;
			}
		}
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSource>
inline void 
assign(_ChunkCollector<THost> & target,
	   TSource & source)
{
	_Assign_Stream_2_ChunkCollector::assign_(target, source);
}
template <typename THost, typename TSource>
inline void 
assign(_ChunkCollector<THost> & target,
	   TSource const & source)
{
	_Assign_Stream_2_ChunkCollector::assign_(target, source);
}

template <typename THost, typename TSource, typename TSize>
inline void 
assign(_ChunkCollector<THost> & target,
	   TSource & source,
	   TSize limit)
{
	_Assign_Stream_2_ChunkCollector::assign_(target, source, limit);
}
template <typename THost, typename TSource, typename TSize>
inline void 
assign(_ChunkCollector<THost> & target,
	   TSource const & source,
	   TSize limit)
{
	_Assign_Stream_2_ChunkCollector::assign_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct _Assign_ChunkCollector_2_String
{
	template <typename TTarget, typename TSource>
	static void assign_(
		TTarget & target, 
		TSource & source)
	{
		typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), TExpand());

		int i_end = chunkCount(source);
		typename Value<TTarget>::Type * pos = begin(target);
		for (int i = 0; i < i_end; ++i)
		{
			bool is_last_chunk = ( part_length <= ChunkLength<TSource>::VALUE);
			typename Size<TTarget>::Type chunk_length = (is_last_chunk) ? part_length : ChunkLength<TSource>::VALUE;
			typename Value<TSource>::Type * chunk = getChunk(source, i);
			
			arrayConstructCopy(chunk, chunk + chunk_length, pos);
			if (is_last_chunk) break;
			pos += chunk_length;
		}
	}

	template <typename TTarget, typename TSource>
	static void assign_(
		TTarget & target, 
		TSource & source,
		typename Size<TTarget>::Type limit)
	{
		typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), limit, TExpand());

		int i_end = chunkCount(source);
		typename Value<TTarget>::Type * pos = begin(target);
		for (int i = 0; i < i_end; ++i)
		{
			bool is_last_chunk = ( part_length <= ChunkLength<TSource>::VALUE);
			typename Size<TTarget>::Type chunk_length = (is_last_chunk) ? part_length : ChunkLength<TSource>::VALUE;
			typename Value<TSource>::Type * chunk = getChunk(source, i);
			
			arrayConstructCopy(chunk, chunk + chunk_length, pos);
			if (is_last_chunk) break;
			pos += chunk_length;
		}
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
assign(String<TTargetValue, TTargetSpec> & target,
	   _ChunkCollector<TSourceHost> const & source,
	   Tag<TExpand> const tag)
{
	_Assign_ChunkCollector_2_String<Tag<TExpand> const>::assign_(target, source);
}
template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
assign(String<TTargetValue, TTargetSpec> & target,
	   _ChunkCollector<TSourceHost> const & source,
	   typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
	   Tag<TExpand> const tag)
{
	_Assign_ChunkCollector_2_String<Tag<TExpand> const>::assign_(target, source, limit);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct _Append_ChunkCollector_2_String
{
	template <typename TTarget, typename TSource>
	static void append_(
		TTarget & target, 
		TSource & source)
	{
		typedef typename Size<TTarget>::Type TSize;
		TSize target_length_old = length(target);
		TSize part_length = _clearSpace(target, length(source), target_length_old, target_length_old, TExpand());

		int i_end = chunkCount(source);
		typename Value<TTarget>::Type * pos = begin(target) + target_length_old; //begin(target) was possibly changed by _clearSpace
		for (int i = 0; i < i_end; ++i)
		{
			bool is_last_chunk = ( part_length <= ChunkLength<TSource>::VALUE);
			typename Size<TTarget>::Type chunk_length = (is_last_chunk) ? part_length : (TSize) ChunkLength<TSource>::VALUE;
			typename Value<TSource>::Type * chunk = getChunk(source, i);
			
			arrayConstructCopy(chunk, chunk + chunk_length, pos);
			if (is_last_chunk) break;
			pos += chunk_length;
		}
	}

	template <typename TTarget, typename TSource>
	static void append_(
		TTarget & target, 
		TSource & source,
		typename Size<TTarget>::Type limit)
	{
		typedef typename Size<TTarget>::Type TSize;
		TSize target_length_old = length(target);
		TSize part_length = _clearSpace(target, length(source), target_length_old, target_length_old, limit, TExpand());

		int i_end = chunkCount(source);
		typename Value<TTarget>::Type * pos = begin(target) + target_length_old; //begin(target) was possibly changed by _clearSpace
		for (int i = 0; i < i_end; ++i)
		{
			bool is_last_chunk = ( part_length <= ChunkLength<TSource>::VALUE);
			typename Size<TTarget>::Type chunk_length = (is_last_chunk) ? part_length : (TSize) ChunkLength<TSource>::VALUE;
			typename Value<TSource>::Type * chunk = getChunk(source, i);
			
			arrayConstructCopy(chunk, chunk + chunk_length, pos);
			if (is_last_chunk) break;
			pos += chunk_length;
		}
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
append(String<TTargetValue, TTargetSpec> & target,
	   _ChunkCollector<TSourceHost> const & source,
	   Tag<TExpand> const )
{
	_Append_ChunkCollector_2_String<Tag<TExpand> const>::append_(target, source);
}
template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
append(String<TTargetValue, TTargetSpec> & target,
	   _ChunkCollector<TSourceHost> const & source,
	   typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
	   Tag<TExpand> const )
{
	_Append_ChunkCollector_2_String<Tag<TExpand> const>::append_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct _Replace_ChunkCollector_2_String
{
	template <typename TTarget, typename TSource>
	static void replace_(
		TTarget & target,
		typename Size<TTarget>::Type pos_begin,
		typename Size<TTarget>::Type pos_end,
		TSource & source)
	{
		typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), pos_begin, pos_end, TExpand());

		int i_end = chunkCount(source);
		typename Value<TTarget>::Type * pos = begin(target) + pos_begin;
		for (int i = 0; i < i_end; ++i)
		{
			bool is_last_chunk = ( part_length <= ChunkLength<TSource>::VALUE);
			typename Size<TTarget>::Type chunk_length = (is_last_chunk) ? part_length : ChunkLength<TSource>::VALUE;
			typename Value<TSource>::Type * chunk = getChunk(source, i);
			
			arrayConstructCopy(chunk, chunk + chunk_length, pos);
			if (is_last_chunk) break;
			pos += chunk_length;
		}
	}

	template <typename TTarget, typename TSource>
	static void replace_(
		TTarget & target, 
		typename Size<TTarget>::Type pos_begin,
		typename Size<TTarget>::Type pos_end,
		TSource & source,
		typename Size<TTarget>::Type limit)
	{
		typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), pos_begin, pos_end, limit, TExpand());

		int i_end = chunkCount(source);
		typename Value<TTarget>::Type * pos = begin(target) + pos_begin;
		for (int i = 0; i < i_end; ++i)
		{
			bool is_last_chunk = ( part_length <= ChunkLength<TSource>::VALUE);
			typename Size<TTarget>::Type chunk_length = (is_last_chunk) ? part_length : ChunkLength<TSource>::VALUE;
			typename Value<TSource>::Type * chunk = getChunk(source, i);
			
			arrayConstructCopy(chunk, chunk + chunk_length, pos);
			if (is_last_chunk) break;
			pos += chunk_length;
		}
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
replace(String<TTargetValue, TTargetSpec> & target,
		typename Size< String<TTargetValue, TTargetSpec> >::Type pos_begin,
		typename Size< String<TTargetValue, TTargetSpec> >::Type pos_end,
	   _ChunkCollector<TSourceHost> const & source,
	   Tag<TExpand> const tag)
{
	_Replace_ChunkCollector_2_String<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}

template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
replace(String<TTargetValue, TTargetSpec> & target,
		typename Size< String<TTargetValue, TTargetSpec> >::Type pos_begin,
		typename Size< String<TTargetValue, TTargetSpec> >::Type pos_end,
	   _ChunkCollector<TSourceHost> const & source,
	   typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
	   Tag<TExpand> const tag)
{
	_Replace_ChunkCollector_2_String<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTargetValue, typename TSourceHost, typename TExpand>
inline void 
replace(TTargetValue * target,
		size_t pos_begin,
		size_t pos_end,
		_ChunkCollector<TSourceHost> const & source,
		Tag<TExpand> const tag)
{
	_Replace_ChunkCollector_2_String<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}

template <typename TTargetValue, typename TSourceHost, typename TExpand>
inline void 
replace(TTargetValue * target,
		size_t pos_begin,
		size_t pos_end,
		_ChunkCollector<TSourceHost> const & source,
		size_t limit,
		Tag<TExpand> const tag)
{
	_Replace_ChunkCollector_2_String<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}
//____________________________________________________________________________
/*
template <typename TTargetValue, typename TSourceHost, typename TExpand>
inline void 
replace(TTargetValue * target,
		size_t pos_begin,
		size_t pos_end,
	   _ChunkCollector<TSourceHost> const & source,
	   Tag<TExpand> const tag)
{
	_Replace_ChunkCollector_2_String<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}

template <typename TTargetValue, typename TTargetSpec, typename TSourceHost, typename TExpand>
inline void 
replace(String<TTargetValue, TTargetSpec> & target,
		typename Size< String<TTargetValue, TTargetSpec> >::Type pos_begin,
		typename Size< String<TTargetValue, TTargetSpec> >::Type pos_end,
	   _ChunkCollector<TSourceHost> const & source,
	   typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
	   Tag<TExpand> const tag)
{
	_Replace_ChunkCollector_2_String<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}
*/
//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
