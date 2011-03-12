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
  $Id: find_multiple_shiftand.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_MULTIPLESHIFTAND_H
#define SEQAN_HEADER_FIND_MULTIPLESHIFTAND_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Multiple ShiftAnd
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.MultipleShiftAnd:
..summary: Multiple exact string matching using bit parallelism. The total size of the patterns should fit into a computer word.
..general:Class.Pattern
..cat:Searching
..signature:Pattern<TNeedle, MultipleShiftAnd>
..param.TNeedle:The needle type, a string of keywords.
...type:Class.String
..remarks.text:The types of all keywords in the needle and the haystack have to match.
*/

///.Class.Pattern.param.TSpec.type:Spec.MultipleShiftAnd

struct _MultipleShiftAnd;
typedef Tag<_MultipleShiftAnd> MultipleShiftAnd;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, MultipleShiftAnd> {
//____________________________________________________________________________
private:
	Pattern(Pattern const& other);
	Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
	typedef unsigned int TWord;
	typedef typename Size<TNeedle>::Type TSize;
	Holder<TNeedle> data_needle;
	TWord* table;			// Look up table for each character in the alphabet (called B in "Navarro")
	TWord* prefSufMatch;	// Set of all the prefixes of needle that match a suffix of haystack (called D in "Navarro")
	TWord* di;				// Initialization word
	TWord* df;				// Final test word
	TWord alphabetSize;		// e.g., char --> 256
	TWord totalLength;		// Lenght of concatenated keywords
	TWord blockCount;		// #unsigned ints required to store needle	
	std::deque<Pair<TSize, TSize> > data_keyword;  // All keywords that produced a hit here
	TSize data_keywordIndex;  // Last keyword index
	TSize data_needleLength;  // Last needle length

//____________________________________________________________________________

	Pattern() {
		table = 0;
		prefSufMatch=0;
		di = 0;
		df = 0;
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
		table = 0;
		prefSufMatch=0;
		di = 0;
		df = 0;
		setHost(*this, ndl);
	}

	~Pattern() {
		SEQAN_CHECKPOINT
		if (table != 0) {
			deallocate(this, table, alphabetSize * blockCount);
			deallocate(this, prefSufMatch, blockCount);
			deallocate(this, di, blockCount);
			deallocate(this, df, blockCount);
		}
	}		
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Host Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct Host< Pattern<TNeedle, MultipleShiftAnd> >
{
	typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, MultipleShiftAnd> const>
{
	typedef TNeedle const Type;
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, MultipleShiftAnd> & me, TNeedle2 const & needle) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef typename Value<TNeedle>::Type TKeyword;
	typedef typename Value<TKeyword>::Type TAlphabet;
	if (me.table != 0) {
		deallocate(me, me.table, me.alphabetSize * me.blockCount);
		deallocate(me, me.di, me.blockCount);
		deallocate(me, me.df, me.blockCount);
	}

	typename Iterator<TNeedle2 const, Rooted>::Type it = begin(needle);
	me.totalLength = 0;
	for(;!atEnd(it);goNext(it)) {
		me.totalLength += length(*it);
	}
	me.alphabetSize = ValueSize<TAlphabet>::VALUE;
	if (me.totalLength<1) me.blockCount=1;
	else me.blockCount=((me.totalLength-1) / BitsPerValue<TWord>::VALUE)+1;
			
	allocate (me, me.table, me.blockCount * me.alphabetSize);
	arrayFill (me.table, me.table + me.blockCount * me.alphabetSize, 0);

	allocate (me, me.di, me.blockCount);
	arrayFill (me.di, me.di + me.blockCount, 0);

	allocate (me, me.df, me.blockCount);
	arrayFill (me.df, me.df + me.blockCount, 0);

	goBegin(it);
	TWord j = 0;
	for(;!atEnd(it);goNext(it)) {
		me.di[j / BitsPerValue<TWord>::VALUE] |= (1<<(j%BitsPerValue<TWord>::VALUE));
		for (TWord posInKeyword = 0; posInKeyword < length(*it); ++posInKeyword) {
			// Determine character position in array table
			TWord index = convert<TWord>(getValue(*it,posInKeyword));
			me.table[me.blockCount*index + j / BitsPerValue<TWord>::VALUE] |= (1<<(j%BitsPerValue<TWord>::VALUE));
			++j;
		}
		me.df[(j - 1) / BitsPerValue<TWord>::VALUE] |= (1<<((j-1)%BitsPerValue<TWord>::VALUE));
	}
	setValue(me.data_needle, needle);

	/*
	// Debug code
	std::cout << "Alphabet size: " << me.alphabetSize << ::std::endl;
	std::cout << "Needle length: " << me.totalLength << ::std::endl;
	std::cout << "Block count: " << me.blockCount << ::std::endl;
	std::cout << "K: ";
	goBegin(it);
	for(;!atEnd(it);goNext(it)) {
		std::cout << *it;
	}
	std::cout << ::std::endl;
	std::cout << "Table: " << ::std::endl;
	for(unsigned int i=0;i<me.alphabetSize;++i) {
		//if ((i<97) || (i>122)) continue;
		std::cout << TAlphabet(i) << ": ";
		for(int j=0;j<me.blockCount;++j) {
			for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
				std::cout << ((me.table[me.blockCount*i+j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
	}
	std::cout << "DI and DF: " << ::std::endl;
	std::cout << "I: ";
	for(int j=0;j<me.blockCount;++j) {
		for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
			std::cout << ((me.di[j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
		}
	}
	std::cout << ::std::endl;
	std::cout << "F: ";
	for(int j=0;j<me.blockCount;++j) {
		for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
			std::cout << ((me.df[j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
		}
	}
	std::cout << ::std::endl;
	*/
}

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, MultipleShiftAnd> & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _patternInit (Pattern<TNeedle, MultipleShiftAnd> & me) 
{
SEQAN_CHECKPOINT
	typedef unsigned int TWord;

	if (me.prefSufMatch != 0) {
		deallocate(me, me.prefSufMatch, me.blockCount);
	}
	allocate (me, me.prefSufMatch, me.blockCount);
	arrayFill (me.prefSufMatch, me.prefSufMatch + me.blockCount, 0);
	me.data_keyword.clear();
	me.data_keywordIndex = 0;
}

//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, MultipleShiftAnd>const>::Type & 
host(Pattern<TNeedle, MultipleShiftAnd> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, MultipleShiftAnd>const>::Type & 
host(Pattern<TNeedle, MultipleShiftAnd> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________


template <typename TNeedle>
inline typename Size<TNeedle>::Type
position(Pattern<TNeedle, MultipleShiftAnd> & me)
{
	return me.data_keywordIndex;
}


template <typename TFinder, typename TNeedle>
bool _findShiftAnd_SmallNeedle(TFinder & finder, Pattern<TNeedle, MultipleShiftAnd> & me) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef typename Size<TNeedle>::Type TSize;
	while (!atEnd(finder)) {
		TWord pos = convert<TWord>(*finder);
		me.prefSufMatch[0] = ((me.prefSufMatch[0] << 1) | me.di[0]) & me.table[me.blockCount*pos];

		/*
		// Debug code
		std::cout << "   ";
		for(int j=0;j<me.blockCount;++j) {
			for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
				std::cout << ((me.prefSufMatch[j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
		*/

		if ((me.prefSufMatch[0] & me.df[0]) != 0) {
			// Check which pattern has matched
			typename Iterator<TNeedle, Rooted>::Type it = begin(value(me.data_needle));
			TWord j = 0;
			for(;!atEnd(it);goNext(it)) {
				j += length(*it);
				TWord test = (1<<((j-1)%BitsPerValue<TWord>::VALUE));
				/*
				// Debug code
				std::cout << "Tes";
				for(int j=0;j<me.blockCount;++j) {
					for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
						std::cout << ((test & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
					}
				}
				std::cout << ::std::endl;
				*/
				if ((me.prefSufMatch[0] & test) != 0) {
					me.data_keyword.push_back(Pair<TSize,TSize>(position(it),length(*it)));
				}
			}
			me.data_keywordIndex = (me.data_keyword.front()).i1;
			me.data_needleLength = (me.data_keyword.front()).i2;
			me.data_keyword.pop_front();
			finder -= (me.data_needleLength - 1);
			return true;
		}
		goNext(finder);
	}
	return false;
}

template <typename TFinder, typename TNeedle>
bool _findShiftAnd_LargeNeedle(TFinder & finder, Pattern<TNeedle, MultipleShiftAnd> & me) {
	SEQAN_CHECKPOINT
	typedef typename Size<TNeedle>::Type TSize;
	typedef unsigned int TWord;
	while (!atEnd(finder)) {
		TWord pos = convert<TWord>(*finder);
		TWord carry = 1;
		for(TWord block=0;block<me.blockCount;++block) {
			bool newCarry = ((me.prefSufMatch[block] & (1<< (BitsPerValue<TWord>::VALUE - 1)))!=0); 
			me.prefSufMatch[block]<<=1;
			me.prefSufMatch[block]|=carry;
			carry = newCarry;
		}
		for(TWord block=0;block<me.blockCount;++block) me.prefSufMatch[block] |= me.di[block];
		for(TWord block=0;block<me.blockCount;++block) me.prefSufMatch[block] &= me.table[me.blockCount*pos+block];

		/*
		// Debug code
		std::cout << "   ";
		for(int j=0;j<me.blockCount;++j) {
			for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
				std::cout << ((me.prefSufMatch[j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
		*/

		bool match = false;
		for(TWord block=0;block<me.blockCount;++block) {
			if ((me.prefSufMatch[block] & me.df[block]) != 0) {
				match = true;
				break;
			}
		}
		if (match) {
			// Check which pattern has matched
			typename Iterator<TNeedle, Rooted>::Type it = begin(value(me.data_needle));
			TWord j = 0;
			for(;!atEnd(it);goNext(it)) {
				j += length(*it);
				TWord* test;
				allocate (me, test, me.blockCount);
				arrayFill (test, test + me.blockCount, 0);
				test[(j - 1) / BitsPerValue<TWord>::VALUE] |= (1<<((j-1)%BitsPerValue<TWord>::VALUE));

				/*
				// Debug code
				std::cout << "Tes";
				for(int i=0;i<me.blockCount;++i) {
					for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
						std::cout << ((test[i] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
					}
				}
				std::cout << ::std::endl;
				*/

				if ((me.prefSufMatch[(j - 1) / BitsPerValue<TWord>::VALUE] & test[(j - 1) / BitsPerValue<TWord>::VALUE]) != 0) {
					me.data_keyword.push_back(Pair<TSize,TSize>(position(it),length(*it)));			
				}
				deallocate(me, test, me.blockCount);
			}
			me.data_keywordIndex = (me.data_keyword.front()).i1;
			me.data_needleLength = (me.data_keyword.front()).i2;
			me.data_keyword.pop_front();
			finder -= (me.data_needleLength - 1);
			return true;
		}
		goNext(finder);
	}
	return false;
}

template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, Pattern<TNeedle, MultipleShiftAnd> & me) {
	SEQAN_CHECKPOINT

	// Check for left-over keywords
	if ((!empty(finder)) &&
		(!me.data_keyword.empty())) {
		finder += me.data_needleLength - 1;
		me.data_keywordIndex = (me.data_keyword.front()).i1;
		me.data_needleLength = (me.data_keyword.front()).i2;
		me.data_keyword.pop_front();
		finder -= (me.data_needleLength - 1);
		return true;
	}


	if (empty(finder)) {
		_patternInit(me);
		_finderSetNonEmpty(finder);
	} else
		finder += me.data_needleLength;

	// Fast algorithm for needles < machine word?
	if (me.blockCount == 1) {
		return _findShiftAnd_SmallNeedle(finder, me);
	} else {
		return _findShiftAnd_LargeNeedle(finder, me);
	}
	return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_MULTIPLESHIFTAND_H
