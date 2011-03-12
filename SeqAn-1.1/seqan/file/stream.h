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
  $Id: stream.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_STREAM_H
#define SEQAN_HEADER_STREAM_H

#include <iosfwd>

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
/**
.Adaption."std::iostream":
..summary:Standard library stream classes.
*/

//////////////////////////////////////////////////////////////////////////////
	
template <typename TValue, typename TTraits>
struct Position< ::std::basic_ios<TValue, TTraits> >
{
	typedef typename ::std::basic_ios<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_streambuf<TValue, TTraits> >
{
	typedef typename ::std::basic_streambuf<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_istream<TValue, TTraits> >
{
	typedef typename ::std::basic_istream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_ostream<TValue, TTraits> >
{
	typedef typename ::std::basic_ostream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_iostream<TValue, TTraits> >
{
	typedef typename ::std::basic_iostream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_stringbuf<TValue, TTraits> >
{
	typedef typename ::std::basic_stringbuf<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_istringstream<TValue, TTraits> >
{
	typedef typename ::std::basic_istringstream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_ostringstream<TValue, TTraits> >
{
	typedef typename ::std::basic_ostringstream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_stringstream<TValue, TTraits> >
{
	typedef typename ::std::basic_stringstream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_filebuf<TValue, TTraits> >
{
	typedef typename ::std::basic_filebuf<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_ifstream<TValue, TTraits> >
{
	typedef typename ::std::basic_ifstream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_ofstream<TValue, TTraits> >
{
	typedef typename ::std::basic_ofstream<TValue, TTraits>::pos_type Type;
};
template <typename TValue, typename TTraits>
struct Position< ::std::basic_fstream<TValue, TTraits> >
{
	typedef typename ::std::basic_fstream<TValue, TTraits>::pos_type Type;
};

//////////////////////////////////////////////////////////////////////////////
	
template <typename TValue, typename TTraits>
struct Value< ::std::basic_ios<TValue, TTraits> >
{
	typedef typename ::std::basic_ios<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_streambuf<TValue, TTraits> >
{
	typedef typename ::std::basic_streambuf<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_istream<TValue, TTraits> >
{
	typedef typename ::std::basic_istream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_ostream<TValue, TTraits> >
{
	typedef typename ::std::basic_ostream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_iostream<TValue, TTraits> >
{
	typedef typename ::std::basic_iostream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_stringbuf<TValue, TTraits> >
{
	typedef typename ::std::basic_stringbuf<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_istringstream<TValue, TTraits> >
{
	typedef typename ::std::basic_istringstream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_ostringstream<TValue, TTraits> >
{
	typedef typename ::std::basic_ostringstream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_stringstream<TValue, TTraits> >
{
	typedef typename ::std::basic_stringstream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_filebuf<TValue, TTraits> >
{
	typedef typename ::std::basic_filebuf<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_ifstream<TValue, TTraits> >
{
	typedef typename ::std::basic_ifstream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_ofstream<TValue, TTraits> >
{
	typedef typename ::std::basic_ofstream<TValue, TTraits>::char_type Type;
};
template <typename TValue, typename TTraits>
struct Value< ::std::basic_fstream<TValue, TTraits> >
{
	typedef typename ::std::basic_fstream<TValue, TTraits>::char_type Type;
};

//////////////////////////////////////////////////////////////////////////////

/**.interal._IsTellSeekStream:
..summary:Determines whether stream supports tell and seek functions.
..cat:Metafunction
*/

template <typename T>
struct _IsTellSeekStream
{
	typedef False Type;
};


template <typename TValue, typename TTraits>
struct _IsTellSeekStream< ::std::basic_ifstream<TValue, TTraits> >
{
	typedef True Type;
};
template <typename TValue, typename TTraits>
struct _IsTellSeekStream< ::std::basic_fstream<TValue, TTraits> >
{
	typedef True Type;
};

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamEOF:
..summary:Test stream for being in eof or error state.
..cat:Streams
..signature:_streamEOF(stream)
..param.stream:A stream object.
...type:Adaption."std::iostream"
..returns:$true$, if stream is at end of file or was set to error state, $false$ otherwise.
*/
template <typename TValue, typename TTraits>
inline bool 
_streamEOF(::std::basic_ios<TValue, TTraits> const & me)
{
SEQAN_CHECKPOINT
	return me.eof() || me.fail();
}

//////////////////////////////////////////////////////////////////////////////
 
/**
.Internal._streamRead:
..summary:Read some characters from stream into a buffer.
..cat:Streams
..signature:_streamRead(target, stream, limit)
..param.target:A buffer that is filled.
..param.stream:An input stream.
...type:Adaption."std::iostream"
..param.limit:The maximal number of characters that is read from $stream$.
..returns:The number of characters read from $stream$.
*/
template <typename TValue, typename TTraits>
inline ::std::streamsize 
_streamRead(TValue * target,
			::std::basic_istream<TValue, TTraits> & source,
			::std::streamsize limit)
{
SEQAN_CHECKPOINT
	source.read(target, limit);
	return source.gcount();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamGet:
..summary:Read one character from stream.
..cat:Streams
..signature:_streamGet(stream)
..param.stream:An input stream.
...type:Adaption."std::iostream"
..returns:The character read.
*/

template <typename TValue, typename TTraits>
inline TValue 
_streamGet(::std::basic_istream<TValue, TTraits> & source)
{
SEQAN_CHECKPOINT
	return source.get();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamPeek:
..summary:Return the next character to be read from stream.
..cat:Streams
..signature:_streamPeek(stream)
..param.stream:An input stream.
...type:Adaption."std::iostream"
..returns:The character to be read.
*/

template <typename TValue, typename TTraits>
inline TValue 
_streamPeek(::std::basic_istream<TValue, TTraits> & source)
{
SEQAN_CHECKPOINT
	return source.peek();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamUnget:
..summary:Put the last read character back into stream.
..cat:Streams
..signature:_streamUnget(stream)
..param.stream:An input stream.
...type:Adaption."std::iostream"
*/

template <typename TValue, typename TTraits>
inline void
_streamUnget(::std::basic_istream<TValue, TTraits> & source)
{
SEQAN_CHECKPOINT
	source.unget();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamPut:
..summary:Writes one character to stream.
..cat:Streams
..signature:_streamPut(stream, character)
..param.stream:An input stream.
...type:Adaption."std::iostream"
..param.character:A character that is written to $stream$.
*/

template <typename TValue, typename TTraits, typename TChar>
inline void
_streamPut(::std::basic_ostream<TValue, TTraits> & target,
		   TChar character)
{
SEQAN_CHECKPOINT
	target.put(convert<TValue>(character));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamTellG:
..cat:Streams
..summary:Gets current position of input stream.
..signature:_streamTellG(stream)
..param.stream:An input stream.
...type:Adaption."std::iostream"
..returns:The current position in $stream$.
*/
template <typename TValue, typename TTraits>
inline typename Position< ::std::basic_istream<TValue, TTraits> >::Type
_streamTellG(::std::basic_istream<TValue, TTraits> & me)
{
SEQAN_CHECKPOINT
	return me.tellg();
}

//////////////////////////////////////////////////////////////////////////////
/**
.Internal._streamTellP:
..cat:Streams
..summary:Gets current position of output stream.
..signature:_streamTellP(stream)
..param.stream:An ouput stream.
...type:Adaption."std::iostream"
..returns:The current position in $stream$.
..see:Internal._streamTellG
*/
template <typename TValue, typename TTraits>
inline typename Position< ::std::basic_ostream<TValue, TTraits> >::Type
_streamTellP(::std::basic_ostream<TValue, TTraits> & me)
{
SEQAN_CHECKPOINT
	return me.tellp();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamSeekG:
..summary:Moves input stream to a position.
..cat:Streams
..signature:_streamSeekG(stream, position)
..param.stream:An input stream.
...type:Adaption."std::iostream"
..param.position:A position within the stream.
...remarks:Use @Function._streamTellG@ to get valid stream positions.
..see:Internal._streamTellG
*/
template <typename TValue, typename TTraits>
inline void
_streamSeekG(::std::basic_istream<TValue, TTraits> & me,
	 typename Position< ::std::basic_istream<TValue, TTraits> >::Type pos)
{
SEQAN_CHECKPOINT
	me.clear();
	me.seekg(pos);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamSeekP:
..summary:Moves output stream to a position.
..cat:Streams
..signature:_streamSeekP(stream, position)
..param.stream:An output stream.
...type:Adaption."std::iostream"
..param.position:A position within the stream.
...remarks:Use @Function._streamTellP@ to get valid stream positions.
..see:Internal._streamTellP
..see:Internal._streamSeekG
*/
template <typename TValue, typename TTraits>
inline void
_streamSeekP(::std::basic_ostream<TValue, TTraits> & me,
	 typename Position< ::std::basic_ostream<TValue, TTraits> >::Type pos)
{
SEQAN_CHECKPOINT
	me.clear();
	me.seekp(pos);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamSeek2G:
..summary:Moves input stream position relative to current position.
..cat:Streams
..signature:_streamSeek2G(stream, offset)
..param.stream:An input stream.
...type:Adaption."std::iostream"
..param.offset:The amout the position is changed.
...remarks:If this value is negative.
..see:Internal._streamSeekG
*/
template <typename TValue, typename TTraits>
inline void
_streamSeek2G(::std::basic_istream<TValue, TTraits> & me,
	 int off)
{
SEQAN_CHECKPOINT
	me.seekg(off, ::std::ios_base::cur);
}


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
