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
  $Id: find_score.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_SCORE_H
#define SEQAN_HEADER_FIND_SCORE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// DPSearch
//////////////////////////////////////////////////////////////////////////////

template <typename TScore>
struct DPSearch;


/**
.Spec.DPSearch:
..cat:Searching
..general:Class.Pattern
..summary:A dynamic programming algorithm for approximate string-matching with a user-definable scoring function.
..signature:Pattern<TNeedle, DPSearch<TScore> >
..param.TNeedle:The needle type.
...type:Class.String
..param.TScore:The scoring function.
...type:Class.Score
..remarks.text:The algorithm is based on the Sellers/Needleman-Wunsch dynamic progamming algorithm. 
The $Pattern$ object only contains the right-most column of the dynamic programming matrix.
...note:At the moment, the algorithm only works on linear gap costs.
*/

///.Class.Pattern.param.TSpec.type:Class.Score


template <typename TNeedle, typename TScore>
class Pattern<TNeedle, DPSearch<TScore> >
{
public:
	typedef typename Value<TScore>::Type TScoreValue;

	Holder<TNeedle>		data_needle;
	TScore				data_score;
	TScoreValue			data_limit;
	String<TScoreValue>	data_tab;

public: 
	Pattern(): 
		data_limit(0)
	{ 
SEQAN_CHECKPOINT
	}

	Pattern(TNeedle & _needle, 
			TScore & _score_func, 
			TScoreValue _limit = 0): 
		data_score(_score_func),
		data_limit(_limit)
	{ 
SEQAN_CHECKPOINT
		setHost(*this, _needle);
	}

	Pattern(TNeedle & _needle,
			TScoreValue _limit = 0): 
		data_limit(_limit)
	{ 
SEQAN_CHECKPOINT
		setHost(*this, _needle);
	}

	Pattern(TScoreValue _limit): 
		data_limit(_limit)
	{ 
SEQAN_CHECKPOINT
		create(data_score);
	}

	Pattern(Pattern const & other): 
		data_needle( other.data_needle ),
		data_score( other.data_score ), 
		data_limit( other.data_limit ),
		data_tab( other.data_tab )
	{
SEQAN_CHECKPOINT
	}

	inline Pattern & 
	operator = (Pattern const & other) 
	{ 
SEQAN_CHECKPOINT
		this->data_needle = other.data_needle;
		this->data_score = other.data_score;
		this->data_limit = other.data_limit;
		this->data_tab = other.data_tab;

		return *this;
	}
};


//////////////////////////////////////////////////////////////////////////////
// Host
//////////////////////////////////////////////////////////////////////////////

/* //see find_base.h

template <typename TNeedle, typename TScore>
struct Host< Pattern<TNeedle, DPSearch<TScore> > >
{
	typedef TNeedle Type;
};

template <typename TNeedle, typename TScore>
struct Host< Pattern<TNeedle, DPSearch<TScore> > const>
{
	typedef TNeedle const Type;
};
*/

//////////////////////////////////////////////////////////////////////////////


template <typename TNeedle, typename TScore>
inline typename Host<Pattern<TNeedle, DPSearch<TScore> > >::Type & 
host(Pattern<TNeedle, DPSearch<TScore> > & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle, typename TScore>
inline typename Host<Pattern<TNeedle, DPSearch<TScore> > const>::Type & 
host(Pattern<TNeedle, DPSearch<TScore> >  const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}


//____________________________________________________________________________

template <typename TNeedle, typename TScore, typename TNeedle2>
void 
setHost(Pattern<TNeedle, DPSearch<TScore> > & me, 
		TNeedle2 & ndl)
{
	me.data_needle = ndl;
	clear(me.data_tab);
}
template <typename TNeedle, typename TScore, typename TNeedle2>
void 
setHost(Pattern<TNeedle, DPSearch<TScore> > & me, 
		TNeedle2 const & ndl)
{
	me.data_needle = ndl;
	clear(me.data_tab);
}


//____________________________________________________________________________

/**.Function.scoringScheme
..cat:Searching
..summary:The @glos:scoring scheme@ used for finding or aligning.
..signature:scoringScheme(obj)
..param.obj:Object that holds a @glos:scoring scheme@
...type:Spec.DPSearch
..returns:The @glos:scoring scheme@ used in $obj$
..see:glos:scoring scheme
*/

template <typename TNeedle, typename TScore>
inline TScore const & 
scoringScheme(Pattern<TNeedle, DPSearch<TScore> > & me)
{
SEQAN_CHECKPOINT
	return me.data_score;
}

//____________________________________________________________________________

/**.Function.setScoringScheme
..cat:Searching
..summary:Sets the @glos:scoring scheme@ used for finding or aligning.
..signature:setScoringScheme(obj, score)
..param.obj:Object that holds a @glos:scoring scheme@.
...type:Spec.DPSearch
..param.score:The new @glos:scoring scheme@ used by $obj$.
..see:glos:scoring scheme
..see:Function.scoringScheme
*/

template <typename TNeedle, typename TScore, typename TScore2>
inline void
setScoringScheme(Pattern<TNeedle, DPSearch<TScore> > & me, 
				 TScore2 & score)
{
SEQAN_CHECKPOINT
	me.data_score = score;
	clear(me.data_tab);
}
template <typename TNeedle, typename TScore, typename TScore2>
inline void
setScoringScheme(Pattern<TNeedle, DPSearch<TScore> > & me, 
				 TScore2 const & score)
{
SEQAN_CHECKPOINT
	me.data_score = score;
	clear(me.data_tab);
}

//____________________________________________________________________________


/**.Function.scoreLimit
..cat:Searching
..summary:The minimal score a match must reach in approximate searching.
..signature:scoreLimit(pattern)
..param.pattern:A @Concept.Pattern|pattern@ that can be used for approximate searching.
...type:Spec.DPSearch
..returns:The current score limit of $pattern$.
*/

template <typename TNeedle, typename TScore>
inline typename Value<TScore>::Type 
scoreLimit(Pattern<TNeedle, DPSearch<TScore> > const & me)
{
SEQAN_CHECKPOINT
	return me.data_limit;
}

//____________________________________________________________________________

/**.Function.setScoreLimit
..cat:Searching
..summary:Sets the minimal score a match must reach in approximate searching.
..signature:setScoreLimit(pattern, limit)
..param.pattern:A @Concept.Pattern|pattern@ that can be used for approximate searching.
...type:Spec.DPSearch
..param.limit:The new score limit.
..see:Function.scoreLimit
*/

template <typename TNeedle, typename TScore, typename TScoreValue>
inline void 
setScoreLimit(Pattern<TNeedle, DPSearch<TScore> > & me, 
			  TScoreValue _limit)
{
SEQAN_CHECKPOINT
	me.data_limit = _limit;
}

//____________________________________________________________________________
// returns the score of the last hit position found (note:position = end of occurrence in haystack)

/**.Function.getScore
..cat:Searching
..summary:Score of the last found match in approximate searching.
..signature:getScore(pattern)
..param.pattern:A @Concept.Pattern|pattern@ that can be used for approximate searching.
...type:Spec.DPSearch
..returns:The score of the last match found using $pattern$.
...remarks:If no match was found, the value is undefined.
..see:Function.scoreLimit
..see:Function.setScoreLimit
..see:Function.find
*/

template <typename TNeedle, typename TScore>
inline typename Value<TScore>::Type
getScore(Pattern<TNeedle, DPSearch<TScore> > & me)
{
	return front(me.data_tab);
}

//////////////////////////////////////////////////////////////////////////////


template <typename TNeedle, typename TScore>
inline void _patternInit (Pattern<TNeedle, DPSearch<TScore> > & me) 
{
	typedef Pattern<TNeedle, DPSearch<TScore> > TPattern;
	typedef typename Value<TScore>::Type TScoreValue;

	typedef typename Size<TPattern>::Type TSize;

	typedef String<TScoreValue> TTab;
	typedef typename Iterator<TTab, Standard>::Type TIterator;

	TScoreValue score_gap = scoreGapExtend(scoringScheme(me));

	TTab & string_tab = me.data_tab;

	//allocate enough memory for one column of DP matrix
	TSize need_length = length(needle(me));
	SEQAN_ASSERT(need_length);

	TSize got_length = resize(string_tab, need_length);
	SEQAN_ASSERT(got_length >= need_length);

//	if (length(_dataNeedle(me)) < got_length) throw(0); //???TODO: Throw "not enough memory" exception

	//init matrix
	//note: The column is stored in reverse order
	TIterator tab_end = begin(string_tab, Standard());
	TIterator tab = end(string_tab, Standard());
	
	TScoreValue x = score_gap;

	while (tab > tab_end)
	{
		--tab;
		*tab = x;
		x += score_gap;
	}
}



//////////////////////////////////////////////////////////////////////////////
// find, findNext
//////////////////////////////////////////////////////////////////////////////

//proportional gap cost: Needleman-Wunsch 

//???TODO: Ukkonen trick?
//???TODO: Finder for affine gap costs?

template <typename TFinder, typename TNeedle, typename TScore>
bool 
_find_score_simple_proportional(TFinder & finder, Pattern<TNeedle, DPSearch<TScore> > & me)
{
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TScoreValue> TTab;
	typedef typename Iterator<TTab, Standard>::Type TTabIterator;
	typedef typename Iterator<TNeedle const, Standard>::Type TNeedleIterator;
	typedef typename Value<typename Haystack<TFinder>::Type>::Type THaystackValue;

	String<TScoreValue> & string_tab = me.data_tab;

	TScoreValue score_gap = scoreGapExtend(scoringScheme(me));
	TScoreValue score_match = scoreMatch(scoringScheme(me));
	TScoreValue score_mismatch = scoreMismatch(scoringScheme(me));

	//init table

	if (empty(finder))
	{
		clear(me.data_tab);
		_finderSetNonEmpty(finder);
	}
	else
	{
		goNext(finder);
	}

	if (! length(me.data_tab))
	{
		_patternInit(me);
	}
	//start searching

	TTabIterator tab_begin = end(string_tab, Standard());

	TNeedleIterator it_begin = begin(host(me), Standard());
	TNeedleIterator it_end = end(host(me), Standard());

	//for each character in haystack, do...
	for (; !atEnd(finder); ++finder)
	{
		//get character
		THaystackValue c = *finder;

		//init some variables
		TNeedleIterator it = it_begin;
		TScoreValue * tab = tab_begin;
		TScoreValue h = 0;
		TScoreValue v = 0;

		//fill the column
		while (it < it_end)
		{
			--tab; //note: the column is stored in "reverse order"

			TScoreValue m2 = (c == *it) ? h + score_match : h + score_mismatch;
			h = *tab;
			TScoreValue m1 = (h > v) ? h + score_gap : v + score_gap;

			v = (m1 > m2) ? m1 : m2;
			*tab = v;

			++it;
		}

		if (*tab >= scoreLimit(me) )
		{//found a hit
			return true;
		}

	}

	//found nothing
	return false;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFinder, typename TNeedle, typename TScore>
inline bool 
find(TFinder & finder, 
	 Pattern<TNeedle, DPSearch<TScore> > & me)
{
	SEQAN_ASSERT(scoreGapOpen(scoringScheme(me)) == scoreGapExtend(scoringScheme(me))) //this finder is only defined for linear gap costs
	return _find_score_simple_proportional(finder, me);
}

template <typename TFinder, typename TNeedle, typename TScore>
inline bool 
find(TFinder & finder, 
	 Pattern<TNeedle, DPSearch<TScore> > & me,
	 int const limit_)
{
	SEQAN_ASSERT(scoreGapOpen(scoringScheme(me)) == scoreGapExtend(scoringScheme(me))) //this finder is only defined for linear gap costs
	setScoreLimit(me, limit_);
	return _find_score_simple_proportional(finder, me);
}


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
