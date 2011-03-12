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
  $Id: file_format_fasta_align.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_FASTA_ALIGN_H
#define SEQAN_HEADER_FILE_FASTA_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - Fasta alignment format
//////////////////////////////////////////////////////////////////////////////

//forward declarations
template <typename T>
struct Row;

template <typename T>
struct Rows;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Fasta alignment:
	FASTA alignment file format for sequences.
*/
struct TagFastaAlign_;
typedef Tag<TagFastaAlign_> const FastaAlign;


/////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TSize>
void _fasta_align_scan_line(TFile & file, TSize & count) {

	SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))

	while (true) {
		typename Value<TFile>::Type c = _streamGet(file);

		if (_streamEOF(file)) return;
		if (c == '\n') return;

		if ((c != '\r') && (c!='-')) {
			++count;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// read
//////////////////////////////////////////////////////////////////////////////
template <typename TFile, typename TSource, typename TSpec>
void read(TFile & file, Align<TSource, TSpec> & align, FastaAlign) {
SEQAN_CHECKPOINT

	SEQAN_ASSERT(!_streamEOF(file))
	
	typedef typename Value<TSource>::Type TSourceValue;
	typedef typename Size<TSourceValue>::Type TSize;
	TSize limit = supremumValue<TSize>();

	//Determine begin position, end position and length of each sequence
	String<TSize> beg_end_length;
	
	typename Position<TFile>::Type begin_pos;
	typename Position<TFile>::Type end_pos;
	typename Value<TFile>::Type c;
	TSize count;

	while (!_streamEOF(file)) {
		begin_pos = _streamTellG(file);
		count = 0;
		SEQAN_ASSERT(!_streamEOF(file))

	
		c = _streamGet(file);
		
		// Skip id
		if (c == '>') {
			_fasta_align_scan_line(file, count);
			begin_pos = _streamTellG(file);
			count = 0;
		} else {  //If no id first letter belongs to sequence
			count = 1;
		}

		// Count letters
		while (true) {
			_fasta_align_scan_line(file, count);

			typename Value<TFile>::Type c = _streamGet(file);
			if (c == '>') {
				_streamSeek2G(file, -1);
				end_pos = _streamTellG(file);
				break;
			}
			if (_streamEOF(file)) {
				end_pos = _streamTellG(file);
				break;
			}
			if ((c != '\n') && (c != '\r') && (c!='-'))	{
				++count;
			}
		}
		if (count > limit) {
			count = limit;
		}

		append(beg_end_length, begin_pos);
		append(beg_end_length, end_pos);
		append(beg_end_length, count);
	}

	// Resize alignment data structure
	TSize numRows=length(beg_end_length) / 3;
	resize(rows(align), numRows);	//rows
		
	typedef Align<TSource, TSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	
	for(TSize i=0;i<numRows;++i) {
		TSize begin = beg_end_length[i*3];
//		TSize end = beg_end_length[i*3+1];
		count = beg_end_length[i*3+2];
		
		//Reserve space
		clear(row(align,i));
		createSource(row(align,i));
		resize(source(row(align,i)),count);
		if (length(source(row(align,i))) < count) {
			count = length(source(row(align,i)));
		}
		setSourceEndPosition(row(align,i),count);
		
		//Read sequence
		_streamSeekG(file, begin);

		typename Position<TSource>::Type pos;
		for (pos = 0; pos < count; ) {
			c = _streamGet(file);
			if ((c != '\n') && (c != '\r') && (c != '-'))	{
				source(row(align,i))[pos] = c;
				++pos;
			}
			if (c=='-') {
				insertGap(row(align,i), toViewPosition(row(align,i), pos));
			}
		}
	}

	_streamSeekG(file, 0);
}

//////////////////////////////////////////////////////////////////////////////
// readIDs
//////////////////////////////////////////////////////////////////////////////
 
template <typename TFile, typename TStringContainer>
void readIDs(TFile& file, TStringContainer& ids, FastaAlign) {
	
	SEQAN_CHECKPOINT
	
	SEQAN_ASSERT(!_streamEOF(file))

	typedef typename Value<TStringContainer>::Type TString;
	typename Position<TFile>::Type start_pos;
	typename Value<TFile>::Type c;


	TString id;
	while(true) {
		c = _streamGet(file);
		while ((!_streamEOF(file)) && (c != '>')) c = _streamGet(file);
		if (!_streamEOF(file)) {
			start_pos = _streamTellG(file);
			typename Size<TString>::Type count = 0;
			_fasta_align_scan_line(file, count);
			if (! count) clear(id);
			else {
				resize(id, count);
				if (length(id) < count)	{
					count = length(id);
				}
				_streamSeekG(file, start_pos);
				for (typename Position<TString>::Type pos = 0; pos<count; ++pos) {
					id[pos] = _streamGet(file);
				}
			}
			appendValue(ids, id);
		} else {
			break;
		}
	}
	_streamSeekG(file, 0);
}

//////////////////////////////////////////////////////////////////////////////
// readMeta
//////////////////////////////////////////////////////////////////////////////

//Fasta file records have no meta data

template <typename TFile, typename TMeta>
void readMeta(TFile & file, TMeta & meta, FastaAlign) {
	SEQAN_CHECKPOINT
	clear(meta);
}


//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////
template <typename TFile>
void goNext(TFile & file, FastaAlign) {
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))
	
	return;
}


//////////////////////////////////////////////////////////////////////////////
// write
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
void _write_impl(TFile& file, Align<TSource, TSpec>& align, TStringContainer& ids, FastaAlign) {
	SEQAN_CHECKPOINT

	typedef Align<TSource, TSpec> const TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
	typedef typename Position<TAlign>::Type TPosition;
	TRowsPosition row_count = length(rows(align));

	for(TRowsPosition i=0;i<row_count;++i) {
		TRow & row_ = row(align, i);	
	
		typedef typename Iterator<typename Row<TAlign>::Type const, Standard>::Type TIter;
		TIter begin_ = iter(row_, beginPosition(cols(align)));
		TIter end_ = iter(row_, endPosition(cols(align)));
	
		_streamPut(file, '>');
		_streamWrite(file, getValue(ids,i));
		_streamPut(file, '\n');

		int chars=0;
		while(begin_ != end_) {
			if (chars == 60) {
				_streamPut(file, '\n');
				chars = 0;
			}
			if (isGap(begin_)) _streamPut(file, gapValue<char>());
			else _streamPut(file, getValue(source(begin_)));
			chars++;
			++begin_;
		}
		_streamPut(file, '\n');
	}
}

//____________________________________________________________________________

template <typename TFile, typename TSource, typename TSpec>
void write(TFile & file, Align<TSource, TSpec>& align, FastaAlign) {
	SEQAN_CHECKPOINT
	_write_impl(file, align, String<String<char> >(), FastaAlign());
}

//____________________________________________________________________________

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
void write(TFile & file, Align<TSource, TSpec> & align, TStringContainer& ids, FastaAlign) {
	SEQAN_CHECKPOINT
	_write_impl(file, align, ids, FastaAlign());
}


//VisualC++ const array bug workaround
template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
void write(TFile & file, Align<TSource, TSpec>* align, TStringContainer & ids, FastaAlign) {
	SEQAN_CHECKPOINT
	_write_impl(file, align, ids, FastaAlign());
}

//____________________________________________________________________________

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec, typename TMeta>
void write(TFile & file, Align<TSource, TSpec> & align, TStringContainer& ids, TMeta &, FastaAlign) {
	SEQAN_CHECKPOINT
	_write_impl(file, align, ids, FastaAlign());
}



//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
