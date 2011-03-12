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
  $Id: file_format.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_FORMAT_H
#define SEQAN_HEADER_FILE_FORMAT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format:
..summary:A file format.
*/


//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//Base Class for all FileFormat classes
//////////////////////////////////////////////////////////////////////////////

/**
.Class.FileFormat:
..cat:Input/Output
..summary:Object that stores a file format.
..signature:FileFormat<File, Data [, Format [, Meta] ]>
..see:Tag.File Format
*/

template <
	typename TFile, 
	typename TData,
	typename TMeta,
	typename TFormat = void >
struct FileFormat:
	public FileFormat<TFile, TData, TMeta, void>
{
public:
	typedef typename Size<TData>::Type TSize;

	FileFormat() {}
	FileFormat(FileFormat const &) {}
	~FileFormat() {}
	FileFormat const & operator =(FileFormat const &) {}

	inline void * 
	formatID_() const
	{
SEQAN_CHECKPOINT
		return _ClassIdentifier<TFormat>::getID();
	}

	virtual void
	read_(TFile & file, TData & data) const
	{
SEQAN_CHECKPOINT
		read(file, data, TFormat());
	}
	virtual void
	read_(TFile & file, TData & data, TSize limit) const
	{
SEQAN_CHECKPOINT
		read(file, data, limit, TFormat());
	}

	virtual void
	readMeta_(TFile & file, TMeta & meta) const
	{
SEQAN_CHECKPOINT
		readMeta(file, meta, TFormat());
	}

	virtual void
	goNext_(TFile & file) const
	{
SEQAN_CHECKPOINT
		goNext(file, TFormat());
	}

	virtual TSize
	length_(TFile & file) const
	{
SEQAN_CHECKPOINT
		length(file, TFormat());
	}

	virtual void
	write_(TFile & file, TData & data) const
	{
SEQAN_CHECKPOINT
		write(file, data, TFormat());
	}
	virtual void
	write_(TFile & file, TData & data, TMeta & meta) const
	{
SEQAN_CHECKPOINT
		write(file, data, meta, TFormat());
	}
};

//____________________________________________________________________________

//base class for all file format classes 

template <typename TFile, typename TData, typename TMeta>
struct FileFormat<TFile, TData, TMeta, void>
{
public:
	typedef typename Size<TData>::Type TSize;

	FileFormat() {}
	FileFormat(FileFormat const &) {}
	~FileFormat() {};
	FileFormat const & operator =(FileFormat const &) {}

	virtual void *
	formatID_() const = 0;

	virtual void
	read_(TFile & file, TData & data) const = 0;
	virtual void
	read_(TFile & file, TData & data, TSize limit) const = 0;

	virtual void
	readMeta_(TFile & file, TMeta & meta) const = 0;

	virtual void
	goNext_(TFile & file) const = 0;

	virtual TSize
	length_(TFile & file) const = 0;

	virtual void
	write_(TFile & file, TData & data) const = 0;
	virtual void
	write_(TFile & file, TData & data, TMeta & meta) const = 0;

};

//////////////////////////////////////////////////////////////////////////////
// Wrapper for functions to virtuals
//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void *
formatID(FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	return file_format.formatID_();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.Fileformat#read:
..cat:Input/Output
..summary:Loads a record from file.
..signature:read(file, data [, meta], format)
..signature:read(file, data [, meta], tag)
..param.file:An input file.
..param.data:A container that gets the data read from $file$.
..param.meta:A container that gets meta data from $file$. (optional)
..param.format:A file format object.
...type:Class.FileFormat.File Format object
..param.tag:A file format tag.
...type:Tag.File Format.File Format tag
..remarks:The result of this operation is stored in $data$.
..remarks:The function leaves $file$ at the position for reading the next record.
..see:Function.assign
*/
template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
read(TFile & file,
	 TData & data,
	 FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.read_(file, data);
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TSize>
inline void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.read_(file, data, limit);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.readMeta:
..cat:Input/Output
..summary:Read meta information from file.
..signature:readMeta(file, meta, file_format)
..param.file:A file that contains data in the format specified by $file_format$.
..param.meta:A data structure that is able to store meta informations stored in $file$.
..param.file_format:A file format.
..returns.param.meta:The meta data read from $file$.
...type:Tag.File Format
*/

template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
readMeta(TFile & file,
		 TMeta & meta,
		 FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.readMeta_(file, meta);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.goNext:
..cat:Input/Output
*/

template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
goNext(TFile & file,
	   FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.goNext_(file);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.length:
..cat:Input/Output
*/

template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
length(TFile & file,
	   FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.length_(file);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.Fileformat#write:
..cat:Input/Output
..summary:Writes to stream.
..signature:write(stream, source)
..signature:write(stream, begin, end)
..param.stream: A stream object.
...type:Adaption."std::iostream"
..param.source: Container that is written to $stream$.
..param.begin: Iterator to the first character of the range.
..param.end: Iterator behind the last character of the range.
..remarks:The content of $source$ is written 'as-is' to $stream$.
*/

template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
write(TFile & file,
	  TData & data,
	  FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.write_(file, data);
}
template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
write(TFile & file,
	  TData & data,
	  TMeta & meta,
	  FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.write_(file, data, meta);
}




//////////////////////////////////////////////////////////////////////////////
// Comparison of two FileFormat objects
//////////////////////////////////////////////////////////////////////////////

template <typename TFileLeft, typename TDataLeft, typename TMetaLeft, typename TFormatLeft, typename TFileRight, typename TDataRight, typename TMetaRight, typename TFormatRight>
inline bool
operator == (FileFormat<TFileLeft, TDataLeft, TMetaLeft, TFormatLeft> const & left, 
			 FileFormat<TFileRight, TDataRight, TMetaRight, TFormatRight> const & right)
{
SEQAN_CHECKPOINT
	return formatID(left) == formatID(right);
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TFormat2>
inline bool
operator == (FileFormat<TFile, TData, TMeta, TFormat> const & left, 
			 Tag<TFormat2> const)
{
SEQAN_CHECKPOINT
	return formatID(left) == _ClassIdentifier<Tag<TFormat2> const>::getID();
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TFormat2>
inline bool
operator == (Tag<TFormat2> const,
			 FileFormat<TFile, TData, TMeta, TFormat> const & right)
{
SEQAN_CHECKPOINT
	return _ClassIdentifier<Tag<TFormat2> const>::getID() == formatID(right);
}

//____________________________________________________________________________

template <typename TFileLeft, typename TDataLeft, typename TMetaLeft, typename TFormatLeft, typename TFileRight, typename TDataRight, typename TMetaRight, typename TFormatRight>
inline bool
operator != (FileFormat<TFileLeft, TDataLeft, TMetaLeft, TFormatLeft> const & left, 
			 FileFormat<TFileRight, TDataRight, TMetaRight, TFormatRight> const & right)
{
SEQAN_CHECKPOINT
	return formatID(left) != formatID(right);
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TFormat2>
inline bool
operator != (FileFormat<TFile, TData, TMeta, TFormat> const & left, 
			 Tag<TFormat2> const)
{
SEQAN_CHECKPOINT
	return formatID(left) != _ClassIdentifier<Tag<TFormat2> const>::getID();
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TFormat2>
inline bool
operator != (Tag<TFormat2> const,
			 FileFormat<TFile, TData, TMeta, TFormat> const & right)
{
SEQAN_CHECKPOINT
	return _ClassIdentifier<Tag<TFormat2> const>::getID() != formatID(right);
}

//////////////////////////////////////////////////////////////////////////////
// allgemeine Funktionen fuer Streams
//////////////////////////////////////////////////////////////////////////////
//TODO??? Das muss in eine extra Datei


/*
template <typename TStream, typename TIterator>
inline void
write(TStream & target,
	  TIterator begin_,
	  TIterator end_)
{
	while (begin_ != end_)
	{
		_streamPut(target, convert<char>(*begin_));
		++begin_;
	}
}

//____________________________________________________________________________

template <typename TStream, typename TSource>
inline void
write(TStream & target,
	  TSource const & source)
{
	write(target, begin(source), end(source));
}
//TODO???: Spezialisierungen zum blockweise schreiben bei contiguous strings von char
//Anmerkungen: write wird nach dem zweiten Argument (source) spezialisiert!

//____________________________________________________________________________

template <typename TStream, typename TSource>
inline void
write(TStream & target,
	  TSource const & source,
	  typename Size<TSource>::Type limit_)
{
	if (length(source) > limit_)
	{
		write(target, begin(source), begin(source) + limit_);
	}
	else
	{
		write(target, begin(source), end(source));
	}
}

*/
//////////////////////////////////////////////////////////////////////////////

// Helper function for scanning a stream
// c = next character, pass it to the next call of the function

template <typename TFile, typename TString, typename TChar>
inline void
_stream_appendLine(TFile & file,
				   TString & str,
				   TChar & c)
{
	while (true)
	{
		if (_streamEOF(file)) break;

		if (c == '\r')
		{
			c = _streamGet(file);
			if (c == '\n') 
			{
				c = _streamGet(file);
			}
			break;
		}
		if (c == '\n')
		{
			c = _streamGet(file);
			break;
		}

		appendValue(str, c);

		c = _streamGet(file);
	}
}
//____________________________________________________________________________

template <typename TFile, typename TChar>
inline void
_stream_countLine(TFile & file,
				  TChar & c)

{
	while (true)
	{
		if (_streamEOF(file)) break;

		if (c == '\r')
		{
			c = _streamGet(file);
			if (c == '\n') 
			{
				c = _streamGet(file);
			}
			break;
		}
		if (c == '\n')
		{
			c = _streamGet(file);
			break;
		}

		c = _streamGet(file);
	}
}

//____________________________________________________________________________

template <typename TFile, typename TChar>
inline typename Size<TFile>::Type
_stream_skipLine(TFile & file,
				 TChar & c)

{
	typename Size<TFile>::Type count = 0;
	while (true)
	{
		if (_streamEOF(file)) break;

		if (c == '\r')
		{
			c = _streamGet(file);
			if (c == '\n') 
			{
				c = _streamGet(file);
			}
			break;
		}
		if (c == '\n')
		{
			c = _streamGet(file);
			break;
		}

		++count;

		c = _streamGet(file);
	}

	return count;
}




////////////////////////////////////////////////////////////////////////////
//new ones

//new ones for streams
template<typename TFile, typename TChar>
inline void 
_stream_skipWhitespace(TFile& file, TChar& c)
{
	if ((c!=' ') && (c != '\t')) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((c!=' ') && (c != '\t')) break;
	}
}

////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TChar>
inline String<char>
_stream_readWord(TFile & file, TChar& c)
{
	// Read word
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_stream_isLetter(c)) break;
		append(str, c);
	}
	return str;
}

////////////////////////////////////////////////////////////////////////////

template<typename TChar>
inline bool
_stream_isLetter(TChar const c)
{
	return ((c == 'a') || (c == 'b') || (c == 'c') || (c == 'd') || (c == 'e') || 
			(c == 'f') || (c == 'g') || (c == 'h') || (c == 'i') || (c == 'j') ||
			(c == 'k') || (c == 'l') || (c == 'm') || (c == 'n') || (c == 'o') || 
			(c == 'p') || (c == 'q') || (c == 'r') || (c == 's') || (c == 't') ||
			(c == 'u') || (c == 'v') || (c == 'w') || (c == 'x') || (c == 'y') || 
			(c == 'z') || (c == 'A') || (c == 'B') || (c == 'C') || (c == 'D') ||
			(c == 'E') || (c == 'F') || (c == 'G') || (c == 'H') || (c == 'I') || 
			(c == 'J') || (c == 'K') || (c == 'L') || (c == 'M') || (c == 'N') ||
			(c == 'O') || (c == 'P') || (c == 'Q') || (c == 'R') || (c == 'S') || 
			(c == 'T') || (c == 'U') || (c == 'V') || (c == 'W') || (c == 'X') ||
			(c == 'Y') || (c == 'Z'));
}


//////////////////////////////////////////////////////////////////////////////

//new ones for strings

template <typename TString, typename TIter>
inline typename Size<TString>::Type
_string_skipLine(TString & str,
				 TIter & it)

{
	typename Size<TString>::Type count = 0;
	typename Iterator<TString,Standard>::Type end_it = end(str,Standard());
	while (true)
	{
		if (it == end_it) break;

		if (*it == '\r')
		{
			++it;
			if (*it == '\n') 
			{
				++it;
			}
			break;
		}
		if (*it == '\n')
		{
			++it;
			break;
		}

		++count;
		++it;
	}

	return count;
}

/////////////////////////////////////////////////////////////////////////

template <typename TString1, typename TString2, typename TIter>
inline void
_string_appendLine(TString1 & str,
				   TString2 & a_str,
				   TIter & it)
{
	typename Iterator<TString1,Standard>::Type end_it = end(str,Standard());
	while (true)
	{
		if (it == end_it) break;

		if (*it == '\r')
		{
			++it; 
			if (*it == '\n') 
			{
				++it;
			}
			break;
		}
		if (*it == '\n')
		{
			++it;
			break;
		}

		appendValue(a_str, getValue(it));
		++it;
	}
}

////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TIter>
inline void 
_string_skipWhitespace(TString& str, TIter& it)
{
	typename Iterator<TString,Standard>::Type end_it = end(str,Standard())-1;
	while (it != end_it) {
		if ((*it!=' ') && (*it != '\t')) break;
		++it;
	}
}

////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TIter>
inline int
_string_readNumber(TString & str, TIter& it)
{
	// Read number
	typename Iterator<TString,Standard>::Type end_it = end(str,Standard())-1;
	String<char> numstr(getValue(it));
	while (it != end_it) {
		++it;
		if (!_parse_isDigit(*it)) break;
		append(numstr, getValue(it));
	}
 	return atoi(toCString(numstr));
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
