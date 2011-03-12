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
  $Id: file_format_guess.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_GUESS_H
#define SEQAN_HEADER_FILE_GUESS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// guessFileFormat
//////////////////////////////////////////////////////////////////////////////

//guessFileFormat braucht auch data, weil die FileFormat-Klasse von TData
//abhaengig, und das ist so, weil sonst die Kombination von Templates mit
//virtuellen Funktionen nicht funktionieren wuerde.
/**
.Function.guessFileFormat:
..cat:Input/Output
..summary:Tries to determine the format of a file.
..signature:guessFileFormat(file, data)
..param.file: An input file.
..param.data: The target container.
...remarks:This container is not modified by this function.
..returns:A file format object instance that represents the determined file format.
...type:Class.FileFormat
..remarks:The $data$-argument is used here as a tag to determine the type of the target.
..see:Function.Fileformat#read
..see:Tag.File Format
*/
template <typename TFile, typename TData, typename TMeta>
inline FileFormat<TFile, TData, TMeta, void> &
guessFileFormat(TFile & file,
				TData & data)
{
SEQAN_CHECKPOINT
	typename Position<TFile>::Type old_pos = _streamTellG(file);
	typename Value<TFile>::Type c;

	_streamSeekG(file, 0); /// move to beginning of file
	c = _streamGet(file);
		
	if (c=='>') 
	{
		_streamSeekG(file, old_pos);
		return getFileFormatInstance<TFile, TData, Fasta, TMeta>();
	}
	
	if (c=='L')
	{
		_streamSeekG(file, old_pos);
		return getFileFormatInstance<TFile, TData, Genbank, TMeta>();
	}

	if (c=='I')
	{
		_streamSeekG(file, old_pos);
		return getFileFormatInstance<TFile, TData, Embl, TMeta>();
	}

	else
	{
		_streamSeekG(file, old_pos);
		return getFileFormatInstance<TFile, TData, Raw, TMeta>();
	}
}

//////////////////////////////////////////////////////////////////////////////

/* DOCH NICHT:
template <typename TTarget, typename TSource>
inline void
read(TTarget & target,
	 TSource & source)
{
SEQAN_CHECKPOINT
	read(target, source, guessFileFormat(target, source));
}
*/

//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
