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
  $Id: cstream.h,v 1.2 2009/03/03 18:47:37 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_CSTREAM_H
#define SEQAN_HEADER_CSTREAM_H

#include <cstdio>

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
/**
.Adaption."std::FILE *":
..summary:Standard library C style streams.
*/

//////////////////////////////////////////////////////////////////////////////
// Position is now defined in file/file_cstyle.h
/*
template <>
struct Position<FILE *>
{
	typedef long Type;
};
*/
//////////////////////////////////////////////////////////////////////////////

template <>
struct Value<FILE *>
{
	typedef char Type;
};

//////////////////////////////////////////////////////////////////////////////
/*
template <>
struct Position<FILE *>
{
	typedef ::std::fpos_t Type;
};
*/

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct _IsTellSeekStream;

template <>
struct _IsTellSeekStream<FILE *>
{
	typedef True Type;
};

//////////////////////////////////////////////////////////////////////////////

inline bool
_streamOpen(::std::FILE * & me, String<char> path, bool for_read = true)
{
SEQAN_CHECKPOINT
	size_t plen = length(path);
	char *s = new char[plen + 1];
	for(size_t i = 0; i < plen; i++) {
		s[i] = path[i];
	}
	s[plen] = '\0';
	if (for_read)
	{
		me = fopen(s, "rb");
	}
	else
	{
		me = fopen(s, "wb");
	}
	delete[] s;
	return (me != 0);
}


//////////////////////////////////////////////////////////////////////////////

inline void
_streamClose(::std::FILE * & me)
{
SEQAN_CHECKPOINT
	if (me)
	{
		fclose(me);
		me = 0;
	}
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamEOF.param.stream.type:Adaption."std::FILE *"

inline bool
_streamEOF(::std::FILE * me)
{
SEQAN_CHECKPOINT
	return feof(me) || ferror(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamRead.param.stream.type:Adaption."std::FILE *"

template <typename TValue>
inline size_t
_streamRead(TValue * target,
			::std::FILE * source,
			size_t limit)
{
SEQAN_CHECKPOINT
	return ::std::fread(target, sizeof(TValue), limit, source);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamGet.param.stream.type:Adaption."std::FILE *"

inline char
_streamGet(::std::FILE * source)
{
SEQAN_CHECKPOINT
	return getc(source);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamPut.param.stream.type:Adaption."std::FILE *"

inline void
_streamPut(::std::FILE * target,
		   char character)
{
SEQAN_CHECKPOINT
	putc(character, target);
}


//////////////////////////////////////////////////////////////////////////////

///.Internal._streamPut.param.stream.type:Adaption."std::FILE *"

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamTellG.param.stream.type:Adaption."std::FILE *"

inline Position<FILE *>::Type
_streamTellG(FILE * me)
{
SEQAN_CHECKPOINT
	return ::std::ftell(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamTellP.param.stream.type:Adaption."std::FILE *"

inline Position<FILE *>::Type
_streamTellP(FILE * me)
{
SEQAN_CHECKPOINT
	return ::std::ftell(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamSeekG.param.stream.type:Adaption."std::FILE *"

inline void
_streamSeekG(FILE * me,
			 Position<FILE *>::Type pos)
{
SEQAN_CHECKPOINT
	::std::fseek(me, pos, SEEK_SET);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamSeekP.param.stream.type:Adaption."std::FILE *"

inline void
_streamSeekP(FILE * me,
			 Position<FILE *>::Type pos)
{
SEQAN_CHECKPOINT
	::std::fseek(me, pos, SEEK_SET);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamSeek2G.param.stream.type:Adaption."std::FILE *"

inline void
_streamSeek2G(FILE * me,
	 int off)
{
SEQAN_CHECKPOINT
	::std::fseek(me, off, SEEK_CUR);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamUnget.param.stream.type:Adaption."std::FILE *"

inline void
_streamUnget(::std::FILE * stream)
{
SEQAN_CHECKPOINT
	_streamSeek2G(stream, -1);
}

//////////////////////////////////////////////////////////////////////////////
// Stream operators for FILE *
//////////////////////////////////////////////////////////////////////////////

// ISO C++ operators are only allowed for classes, not for pointers

/*
template <typename TSource>
inline FILE *
operator << (FILE * target,
			 TSource & source)
{
SEQAN_CHECKPOINT
	write(target, source);
	return target;
}
template <typename TSource>
inline FILE *
operator << (FILE * target,
			 TSource const & source)
{
SEQAN_CHECKPOINT
	write(target, source);
	return target;
}

//____________________________________________________________________________

template <typename TTarget>
inline FILE *
operator >> (FILE * source,
			 TTarget & target)
{
SEQAN_CHECKPOINT
	read(source, target);
	return source;
}
template <typename TTarget>
inline FILE *
operator >> (FILE * source,
			 TTarget const & target)
{
SEQAN_CHECKPOINT
	read(source, target);
	return source;
}
*/

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
