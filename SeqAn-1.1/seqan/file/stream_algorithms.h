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
  $Id: stream_algorithms.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_STREAM_ALGORITHMS_H
#define SEQAN_HEADER_STREAM_ALGORITHMS_H


namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////


/**
.Internal._streamPutInt:
..summary:Converts an integer to a character and writes it to stream.
..cat:Streams
..signature:_streamPutInt(stream, number [, format_string])
..param.target:An output stream.
...type:Adaption."std::iostream"
..param.number:A number that is written to $stream$.
*/
template <typename TStream>
inline void
_streamPutInt(TStream & target,
			  int number, 
			  char const * format_string)
{
SEQAN_CHECKPOINT
	char str[BitsPerValue<int>::VALUE];
	sprintf(str, format_string, number);
	_streamWrite(target, str);
}
template <typename TStream>
inline void
_streamPutInt(TStream & target,
			  int number)
{
SEQAN_CHECKPOINT
	_streamPutInt(target, number, "%d");
}

/**
.Internal._streamPutFloat:
..summary:Converts a float to a character and writes it to stream.
..cat:Streams
..signature:_streamPutFloat(stream, number [, format_string])
..param.target:An output stream.
...type:Adaption."std::iostream"
..param.number:A number that is written to $stream$.
*/
template <typename TStream>
inline void
_streamPutFloat(TStream & target,
			  float number, 
			  char const * format_string)
{
SEQAN_CHECKPOINT
	char str[BitsPerValue<float>::VALUE];
	sprintf(str, format_string, number);
	_streamWrite(target, str);
}
template <typename TStream>
inline void
_streamPutFloat(TStream & target,
				float number)
{
SEQAN_CHECKPOINT
	_streamPutFloat(target, number, "%f");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename T1, typename T2, typename TCompression>
inline void
_streamWrite(TTarget & target,
			 Pair<T1, T2, TCompression> const & source)
{
SEQAN_CHECKPOINT
	_streamWrite(target, getValueI1(source));
	_streamWrite(target, getValueI2(source));
}

template <typename TTarget, typename T1, typename T2, typename T3, typename TCompression>
inline void
_streamWrite(TTarget & target,
			 Triple<T1, T2, T3, TCompression> const & source)
{
SEQAN_CHECKPOINT
	_streamWrite(target, getValueI1(source));
	_streamWrite(target, getValueI2(source));
	_streamWrite(target, getValueI3(source));
}

//////////////////////////////////////////////////////////////////////////////



/**
.Internal._streamWrite:
..summary:Writes a sequence to stream.
..cat:Streams
..signature:_streamWrite(stream, sequence)
..param.stream:An input stream.
..param.sequence:A sequence that is written to $stream$.
*/

template <typename TTarget, typename TSource>
inline void
_streamWrite(TTarget & target,
			 TSource const & source)
{
SEQAN_CHECKPOINT
	_streamWriteSeq(target, source, typename IsSequence<TSource const>::Type());
}

//____________________________________________________________________________

template <typename TTarget, typename TSource>
inline void
_streamWriteSeq(TTarget & target,
				TSource const & source,
				False const)
{
	_streamPut(target, source);
}

//____________________________________________________________________________

template <typename TTarget, typename TSource>
inline void
_streamWriteSeq(TTarget & target,
				TSource const & source,
				True const)
{
SEQAN_CHECKPOINT
	typename Iterator<TSource const, Standard>::Type it = begin(source, Standard());
	typename Iterator<TSource const, Standard>::Type it_end = end(source, Standard());

	for (; it < it_end; ++it)
	{
		typename GetValue<TSource const>::Type val_ = getValue(it);
		_streamWrite(target, val_);
	}
}

template <typename TTarget, typename TSourceValue>
inline void
_streamWriteSeq(TTarget & target,
			    TSourceValue const * source,
				True const)
{
SEQAN_CHECKPOINT

	for (; !atEnd(source); ++source)
	{
		_streamWrite(target, *source);
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamWriteRange:
..summary:Writes a range to stream.
..cat:Streams
..signature:_streamWriteRange(stream, begin_iterator, end_iterator)
..param.stream:An input stream.
..param.sequence:A sequence that is written to $stream$.
*/

template <typename TTarget, typename TIterator>
inline void
_streamWriteRange(TTarget & target,
				  TIterator begin_,
				  TIterator end_)
{
SEQAN_CHECKPOINT
	for (; begin_ != end_; ++begin_)
	{
		_streamPut(target, *begin_);
	}
}


	

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
