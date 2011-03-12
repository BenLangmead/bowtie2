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
  $Id: sequence_shortcuts.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEQUENCE_SHORTCUTS_H
#define SEQAN_HEADER_SEQUENCE_SHORTCUTS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
/**
.Shortcut.CharString:
..cat:Strings
..summary:A string of $char$.
..signature:CharString
..shortcutfor:Spec.Alloc String
...signature:String<char, Alloc<> >
*/

typedef String<char, Alloc<void> > CharString;

//____________________________________________________________________________

/**
.Shortcut.CharIterator:
..cat:Iterators
..summary:Iterator for @Shortcut.CharString@.
..signature:CharIterator
..shortcutfor:Concept.Rooted Iterator
...signature:Iterator<CharString, Rooted>::Type
..see:Shortcut.CharString
*/

typedef Iterator<CharString, Rooted>::Type CharIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Shortcut.UnicodeString:
..cat:Strings
..summary:A string of $wchar_t$.
..signature:UnicodeString
..shortcutfor:Spec.Alloc String
...signature:String<wchar_t, Alloc<> >
*/

typedef String<wchar_t, Alloc<void> > UnicodeString;

//____________________________________________________________________________

/**
.Shortcut.UnicodeIterator:
..cat:Iterators
..summary:Iterator for @Shortcut.UnicodeString@.
..signature:UnicodeIterator
..shortcutfor:Concept.Rooted Iterator
...signature:Iterator<UnicodeString, Rooted>::Type
..see:Shortcut.UnicodeString
*/

typedef Iterator<UnicodeString, Rooted>::Type UnicodeIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Shortcut.DnaString:
..cat:Strings
..summary:A string of @Spec.Dna@.
..signature:DnaString
..shortcutfor:Spec.Alloc String
...signature:String<Dna, Alloc<> >
..see:Spec.Dna
*/

typedef String<Dna, Alloc<void> > DnaString;

//____________________________________________________________________________

/**
.Shortcut.DnaIterator:
..cat:Iterators
..summary:Iterator for @Shortcut.DnaString@.
..signature:DnaIterator
..shortcutfor:Concept.Rooted Iterator
...signature:Iterator<DnaString, Rooted>::Type
..see:Spec.Dna
..see:Shortcut.DnaString
*/

typedef Iterator<DnaString, Rooted>::Type DnaIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Shortcut.Dna5String:
..cat:Strings
..summary:A string of @Spec.Dna5@.
..signature:Dna5String
..shortcutfor:Spec.Alloc String
...signature:String<Dna5, Alloc<> >
..see:Spec.Dna5
..see:Shortcut.DnaString
*/

typedef String<Dna5, Alloc<void> > Dna5String;

//____________________________________________________________________________

/**
.Shortcut.Dna5Iterator:
..cat:Iterators
..summary:Iterator for @Shortcut.Dna5String@.
..signature:Dna5Iterator
..shortcutfor:Concept.Rooted Iterator
...signature:Iterator<Dna5String, Rooted>::Type
..see:Spec.Dna5
..see:Shortcut.Dna5String
..see:Shortcut.DnaIterator
*/

typedef Iterator<Dna5String, Rooted>::Type Dna5Iterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Shortcut.IupacString:
..cat:Strings
..summary:A string of @Spec.Iupac@.
..signature:IupacString
..shortcutfor:Spec.Alloc String
...signature:String<Iupac, Alloc<> >
..see:Spec.Iupac
*/

typedef String<Iupac, Alloc<void> > IupacString;

//____________________________________________________________________________

/**
.Shortcut.IupacIterator:
..cat:Iterators
..summary:Iterator for @Shortcut.IupacString@.
..signature:IupacIterator
..shortcutfor:Concept.Rooted Iterator
...signature:Iterator<IupacString, Rooted>::Type
..see:Spec.Iupac
..see:Shortcut.IupacString
*/

typedef Iterator<IupacString, Rooted>::Type IupacIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Shortcut.Peptide:
..cat:Strings
..summary:A string of @Spec.AminoAcid@.
..signature:IupacString
..shortcutfor:Spec.Alloc String
...signature:String<AminoAcid, Alloc<> >
..see:Spec.AminoAcid
*/

typedef String<AminoAcid, Alloc<void> > Peptide;

//____________________________________________________________________________

/**
.Shortcut.PeptideIterator:
..cat:Iterators
..summary:Iterator for @Shortcut.Peptide@.
..signature:PeptideIterator
..shortcutfor:Concept.Rooted Iterator
...signature:Iterator<Peptide, Rooted>::Type
..see:Spec.AminoAcid
..see:Shortcut.Peptide
*/

typedef Iterator<Peptide, Rooted>::Type PeptideIterator;

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
