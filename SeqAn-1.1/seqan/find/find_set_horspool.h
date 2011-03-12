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
  $Id: find_set_horspool.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_SETHORSPOOL_H
#define SEQAN_HEADER_FIND_SETHORSPOOL_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Set Horspool Algorithm
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.SetHorspool:
..summary: Multiple exact string matching using set horspool algorithm.
..general:Class.Pattern
..cat:Searching
..signature:Pattern<TNeedle, SetHorspool>
..param.TNeedle:The needle type, a string of keywords.
...type:Class.String
..remarks.text:The types of all keywords in the needle and the haystack have to match.
*/

///.Class.Pattern.param.TSpec.type:Spec.SetHorspool

struct _SetHorspool;
typedef Tag<_SetHorspool> SetHorspool;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, SetHorspool> {
//____________________________________________________________________________
private:
	Pattern(Pattern const& other);
	Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
	typedef typename Size<TNeedle>::Type TSize;
	typedef typename Value<TNeedle>::Type TValue;
	typedef typename Value<TValue>::Type TAlphabet;
	typedef Graph<Automaton<TAlphabet> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	
	Holder<TNeedle> data_needle;
	Graph<Automaton<TAlphabet> > data_reverseTrie;  // Search trie
	String<String<TSize> > data_terminalStateMap;
	String<TSize> data_dMap;	// Jump table
	TSize data_lmin;
	String<TSize> data_endPositions;	// All remaining keyword indices
	TSize data_keywordIndex;			// Current keyword that produced a hit
	TSize data_needleLength;			// Last length of needle to reposition finder
	TVertexDescriptor data_lastState;   // Last state in the trie

//____________________________________________________________________________

	Pattern() {
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
		setHost(*this, ndl);
	}

	~Pattern() {
		SEQAN_CHECKPOINT
	}		
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Host Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct Host< Pattern<TNeedle, SetHorspool> >
{
	typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, SetHorspool> const>
{
	typedef TNeedle const Type;
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, SetHorspool> & me, TNeedle2 const & needle) {
	SEQAN_CHECKPOINT
	typedef typename Value<TNeedle>::Type TKeyword;
	typedef typename Size<TKeyword>::Type TSize;
	typedef typename Value<TKeyword>::Type TAlphabet;

	// clean-up
	clear(me.data_reverseTrie);
	clear(me.data_terminalStateMap);
	clear(me.data_endPositions);
	clear(me.data_dMap);
	me.data_lmin=0;

	// Create Trie
	createTrieOnReverse(me.data_reverseTrie,me.data_terminalStateMap,needle);
	assignRoot(me.data_reverseTrie,0);
	setValue(me.data_needle, needle);

	// Create jump map
	TSize alphabet_size = ValueSize<TAlphabet>::VALUE;
	resize(me.data_dMap, alphabet_size);
	me.data_lmin = _getInfinity<TSize>();
	typename Iterator<TNeedle2 const, Rooted>::Type it = begin(needle);
	for(;!atEnd(it);goNext(it)) {
		TSize tmp = length(*it);
		if (tmp<me.data_lmin) me.data_lmin = tmp;
	}
	for(TSize i=0;i<alphabet_size;++i) {
		me.data_dMap[i]=me.data_lmin;
	}
	goBegin(it);
	for(;!atEnd(it);goNext(it)) {
		for(TSize pos = 0;pos < length(*it) - 1; ++pos) {
			TSize ind = ordValue((TAlphabet)(*it)[pos]);	
			if ((length(*it)- 1 - pos) < me.data_dMap[ind]) {
				me.data_dMap[ind] = (length(*it) - 1 - pos);
			}
		}
	}

	/*
	fstream strm;
	strm.open(TEST_PATH "my_trie.dot", ios_base::out | ios_base::trunc);
	String<String<char> > nodeMap;
	_createTrieNodeNames(me.data_reverseTrie, me.data_terminalStateMap, nodeMap);
	String<String<char> > edgeMap;
	_createEdgeNames(me.data_reverseTrie,edgeMap);
	write(strm,me.data_reverseTrie,nodeMap,edgeMap,DotDrawing());
	strm.close();
	*/
}

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, SetHorspool> & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _patternInit (Pattern<TNeedle, SetHorspool> & me) 
{
SEQAN_CHECKPOINT
	clear(me.data_endPositions);
	me.data_keywordIndex = 0;
	me.data_lastState = getRoot(me.data_reverseTrie);
}


//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, SetHorspool>const>::Type & 
host(Pattern<TNeedle, SetHorspool> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, SetHorspool>const>::Type & 
host(Pattern<TNeedle, SetHorspool> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________


template <typename TNeedle>
inline typename Size<TNeedle>::Type
position(Pattern<TNeedle, SetHorspool> & me)
{
	return me.data_keywordIndex;
}


template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, Pattern<TNeedle, SetHorspool> & me) {
	SEQAN_CHECKPOINT
	typedef typename Value<TNeedle>::Type TKeyword;
	typedef typename Size<TKeyword>::Type TSize;
	typedef typename Value<TKeyword>::Type TAlphabet;
	typedef Graph<Automaton<TAlphabet> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	TVertexDescriptor current = getRoot(me.data_reverseTrie); 

	// Process left-over hits
	if ((!empty(finder)) &&
		(!empty(me.data_endPositions))) {
		finder += me.data_needleLength;
		current = me.data_lastState;
		me.data_keywordIndex = me.data_endPositions[length(me.data_endPositions)-1];
		me.data_needleLength = length(getValue(host(me), me.data_keywordIndex))-1;
		if (length(me.data_endPositions) > 1) resize(me.data_endPositions, (length(me.data_endPositions)-1));
		else clear(me.data_endPositions);
		me.data_lastState = current;
		finder -= me.data_needleLength;
		return true;
	}

	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TSize j = 0;
	if (empty(finder)) {
		_patternInit(me);
		_finderSetNonEmpty(finder);
		finder += me.data_lmin - 1;
	} else {
		finder += me.data_needleLength;
		j = me.data_needleLength + 1;
		current = me.data_lastState;
	}

	TSize haystackLength = length(container(finder));
	bool oldMatch = true;
	// Do not change to !atEnd(finder) because of jump map!
	while(position(finder) < haystackLength) {
		while ((position(finder)>=j) && 
				(getSuccessor(me.data_reverseTrie, current, *(finder-j))!= nilVal))
		{
			me.data_endPositions = getProperty(me.data_terminalStateMap,current);
			if ((!oldMatch) && (!empty(me.data_endPositions))) break;
			current = getSuccessor(me.data_reverseTrie, current, *(finder-j));
			if (current == nilVal) break;
			++j;
			oldMatch = false;
		}
		me.data_endPositions = getProperty(me.data_terminalStateMap,current);
		if ((!oldMatch) &&
			(!empty(me.data_endPositions)))
		{
			me.data_keywordIndex = me.data_endPositions[length(me.data_endPositions)-1];
			me.data_needleLength = length(getValue(host(me), me.data_keywordIndex))-1;
			if (length(me.data_endPositions) > 1) resize(me.data_endPositions, length(me.data_endPositions)-1);
			else clear(me.data_endPositions);
			me.data_lastState = current;
			finder -= me.data_needleLength;
			return true;
		}
		oldMatch = false;
		TSize ind = ordValue(*finder);
		setPosition(finder, position(finder) + getValue(me.data_dMap, ind));
		j = 0;
		current = getRoot(me.data_reverseTrie);
	}
	return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SETHORSPOOL_H
