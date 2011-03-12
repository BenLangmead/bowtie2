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
  $Id: file_format_cgviz.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_CGVIZ_H
#define SEQAN_HEADER_FILE_CGVIZ_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - CGViz
//////////////////////////////////////////////////////////////////////////////


/**
.Tag.File Format.tag.CGViz:
	CGViz file format for sequences. Only output.
*/
struct TagCGViz_;
typedef Tag<TagCGViz_> const CGViz;

/////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
void goNext(TFile & file, CGViz) {
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))
	
	return;
}



//////////////////////////////////////////////////////////////////////////////
// write
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
void _write_impl(TFile& target, Align<TSource, TSpec>& align, TStringContainer& ids, CGViz) {
	SEQAN_CHECKPOINT

	typedef Align<TSource, TSpec> const TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
	typedef typename Position<TAlign>::Type TPosition;
	TRowsPosition row_count = length(rows(align));
	if (row_count < 2) return;

	unsigned int pair=1;
	unsigned int count=0;
	for(TRowsPosition i=0;i<row_count-1;++i) {
		for(TRowsPosition j=i+1;j<row_count;++j) {
			
			// Print header
			_streamWrite(target, "{DATA dat"); 
			_streamPutInt(target, pair);
			_streamPut(target, '\n');
			_streamWrite(target, "[__GLOBAL__] dimension=2:\n"); 
						
			TPosition begin_ = beginPosition(cols(align));
			TPosition end_ = endPosition(cols(align));
		
			bool match = false;
			while(begin_ < end_) {
				if ((row(align, i)[begin_]==row(align, j)[begin_]) && (row(align, i)[begin_]!='-')) {
					if (!match) {
						match=true;
						_streamPutInt(target, toSourcePosition(row(align,i),begin_+1));
						_streamPut(target, ' ');
						_streamPutInt(target, toSourcePosition(row(align,j),begin_+1));
						_streamPut(target, ' ');
					}
				}
				if ((row(align, i)[begin_]!=row(align, j)[begin_]) || (row(align, i)[begin_]=='-') || (row(align, j)[begin_]=='-')) {
					if (match) {
						_streamPutInt(target, toSourcePosition(row(align,i),begin_));
						_streamPut(target, ' ');
						_streamPutInt(target, toSourcePosition(row(align,j),begin_));
						_streamPut(target, '\n');
						match=false;
					}
				}
				begin_++;
			}
			if (match) {
				_streamPutInt(target, toSourcePosition(row(align,i),begin_));
				_streamPut(target, ' ');
				_streamPutInt(target, toSourcePosition(row(align,j),begin_));
				_streamPut(target, '\n');
				match=false;
			}
			_streamPut(target, '}');
			_streamPut(target, '\n');

			// Write footer
			_streamWrite(target, "{GLYPH Glyph");
			_streamPutInt(target, pair);
			_streamPut(target, '\n');
			_streamWrite(target, "drawerName=Lines\n");
			_streamWrite(target, "lineWidth=3\n");
			_streamPut(target, '}');
			_streamPut(target, '\n');
			_streamWrite(target, "{PANE Pane");
			_streamPutInt(target, pair);
			_streamPut(target, '\n');
			_streamWrite(target, "uLabel=");
			_streamWrite(target, getValue(ids,i));
			_streamPut(target, '\n');
			_streamWrite(target, "uStop=");
			_streamPutInt(target, length(source(row(align,i))));
			_streamPut(target, '\n');
			_streamWrite(target, "vLabel=");
			_streamWrite(target, getValue(ids,j));
			_streamPut(target, '\n');
			_streamWrite(target, "vStop=");
			_streamPutInt(target, length(source(row(align,j))));
			_streamPut(target, '\n');
			_streamPut(target, '}');
			_streamPut(target, '\n');
			_streamWrite(target, "{WINDOW Window");
			_streamPutInt(target, pair);
			_streamPut(target, '\n');
			_streamPut(target, '}');
			_streamPut(target, '\n');
			_streamWrite(target, "{FEEDER Feeder<");
			_streamPutInt(target, pair);
			_streamPut(target, '>');
			_streamPut(target, ' ');
			_streamPutInt(target, count);
			_streamPut(target, ' ');
			_streamPutInt(target, count+1);
			_streamPut(target, '\n');
			_streamPut(target, '}');
			_streamPut(target, '\n');
			++count;
			_streamWrite(target, "{THREADER Threader<");
			_streamPutInt(target, pair);
			_streamPut(target, '>');
			_streamPut(target, ' ');
			_streamPutInt(target, count);
			_streamPut(target, ' ');
			_streamPutInt(target, count+1);
			_streamPut(target, '\n');
			_streamPut(target, '}');
			_streamPut(target, '\n');
			++count;
			_streamWrite(target, "{ANCHOR Anchor<");
			_streamPutInt(target, pair);
			_streamPut(target, '>');
			_streamPut(target, ' ');
			_streamPutInt(target, count);
			_streamPut(target, ' ');
			_streamPutInt(target, count+1);
			_streamPut(target, '\n');
			_streamPut(target, '}');
			_streamPut(target, '\n');
			count+=2;
			++pair;
		}
	}
}


//____________________________________________________________________________

template <typename TFile, typename TSource, typename TSpec>
void write(TFile & file, Align<TSource, TSpec>& align, CGViz) {
	SEQAN_CHECKPOINT
	_write_impl(file, align, String<String<char> >(), CGViz());
}

//____________________________________________________________________________

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
void write(TFile & file, Align<TSource, TSpec> & align, TStringContainer& ids, CGViz) {
	SEQAN_CHECKPOINT
	_write_impl(file, align, ids, CGViz());
}


//VisualC++ const array bug workaround
template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
void write(TFile & file, Align<TSource, TSpec>* align, TStringContainer & ids, CGViz) {
	SEQAN_CHECKPOINT
	_write_impl(file, align, ids, CGViz());
}

//____________________________________________________________________________

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec, typename TMeta>
void write(TFile & file, Align<TSource, TSpec> & align, TStringContainer& ids, TMeta &, CGViz) {
	SEQAN_CHECKPOINT
	_write_impl(file, align, ids, CGViz());
}



//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
