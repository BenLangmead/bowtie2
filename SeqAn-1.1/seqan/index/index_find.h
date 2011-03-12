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
  $Id: index_find.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_FIND_H
#define SEQAN_HEADER_INDEX_FIND_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

	template < typename TText, typename TSpec, typename TSpecFinder >
	struct Position< Finder< Index<TText, TSpec>, TSpecFinder > >:
		SAValue< Index<TText, TSpec> > {};


//////////////////////////////////////////////////////////////////////////////
// generic Finder class for all indices containing a suffix array or 
// similar table, where a query result is an interval in this table
//
// your index must specialize the function _findFirstIndex and set
// finder.range to the interval containing the query hits. See:
//
//	template < typename TText, typename TSpec, typename TSpecFinder, typename TPattern >
//	inline void _findFirstIndex(
//		Finder< Index<TText, TSpec>, TSpecFinder > &finder,
//		TPattern const &pattern,
//		ESA_FIND_MLR const)
//	{
//		Index<TText, TSpec> &index = haystack(finder);
//		indexRequire(index, ESA_SA());
//		finder.range = equalRangeSAIterator(indexText(index), indexSA(index), pattern);
//	}

	template < typename TText, typename TSpec, typename TSpecFinder >
	class Finder< Index<TText, TSpec>, TSpecFinder >
	{
    protected:
		typedef Index<TText, TSpec>								TIndex;
		typedef typename Fibre<TIndex, Fibre_SA>::Type			TSA;
		typedef typename Iterator<TSA const, Standard>::Type	TIterator;

	public:
		Holder<TIndex>	index;
		Pair<TIterator>	range;
		TIterator		data_iterator;

		Finder() 
		{
			clear(*this);
		}
		Finder(TIndex &_index): index(_index) 
		{
			clear(*this);
		}
		Finder(TIndex const &_index): index(_index)
		{
			clear(*this);
		}
	};

//____________________________________________________________________________

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline typename _Parameter< Index<TText, TSpec> >::Type 
	host(Finder< Index<TText, TSpec>, TSpecFinder > & me)
	{
SEQAN_CHECKPOINT
		return value(me.index);
	}

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline typename _Parameter< Index<TText, TSpec> >::Type 
	host(Finder< Index<TText, TSpec>, TSpecFinder > const & me)
	{
SEQAN_CHECKPOINT
		return value(me.index);
	}

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline typename _Parameter< Index<TText, TSpec> >::Type 
	container(Finder< Index<TText, TSpec>, TSpecFinder > & me)
	{
SEQAN_CHECKPOINT
		return value(me.index);
	}

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline typename _Parameter< Index<TText, TSpec> >::Type 
	container(Finder< Index<TText, TSpec>, TSpecFinder > const & me)
	{
SEQAN_CHECKPOINT
		return value(me.index);
	}

//____________________________________________________________________________

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline void
	setHost(
		Finder< Index<TText, TSpec>, TSpecFinder > & me, 
		typename _Parameter<Index<TText, TSpec> >::Type container_)
	{
SEQAN_CHECKPOINT
		me.index = container;
	}

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline void
	setContainer(
		Finder< Index<TText, TSpec>, TSpecFinder > & me, 
		typename _Parameter<Index<TText, TSpec> >::Type container_)
	{
SEQAN_CHECKPOINT
		me.index = container;
	}

//____________________________________________________________________________

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline typename Iterator< typename Fibre<Index<TText, TSpec>, Fibre_SA>::Type, Standard>::Type &
	hostIterator(Finder< Index<TText, TSpec>, TSpecFinder > & me)
	{
SEQAN_CHECKPOINT
		return me.data_iterator;
	}

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline typename Iterator< typename Fibre<Index<TText, TSpec>, Fibre_SA>::Type, Standard>::Type const &
	hostIterator(Finder< Index<TText, TSpec>, TSpecFinder > const & me)
	{
SEQAN_CHECKPOINT
		return me.data_iterator;
	}


//____________________________________________________________________________

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline bool
	empty(Finder< Index<TText, TSpec>, TSpecFinder > & me)
	{
SEQAN_CHECKPOINT
		return me.range.i1 == me.range.i2;
	}

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline void
	clear(Finder< Index<TText, TSpec>, TSpecFinder > & me)
	{
SEQAN_CHECKPOINT
		typedef Index<TText, TSpec>						TIndex;
		typedef typename Fibre<TIndex, Fibre_SA>::Type	TSA;
		typedef typename Iterator<TSA, Standard>::Type	TIterator;
		me.range.i1 = me.range.i2 = TIterator();
	}

//____________________________________________________________________________

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline bool
	atBegin(Finder< Index<TText, TSpec>, TSpecFinder > & me)
	{
SEQAN_CHECKPOINT
		return (empty(me) || hostIterator(me) == me.range.i1);
	}

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline bool
	atEnd(Finder< Index<TText, TSpec>, TSpecFinder > & me)
	{
SEQAN_CHECKPOINT
		return (empty(me) || hostIterator(me) == me.range.i2);
	}

//____________________________________________________________________________

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline void
	goBegin(Finder< Index<TText, TSpec>, TSpecFinder > & me)
	{
SEQAN_CHECKPOINT
		hostIterator(me) = me.range.i1;
	}

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline void
	goEnd(Finder< Index<TText, TSpec>, TSpecFinder > & me)
	{
SEQAN_CHECKPOINT
		hostIterator(me) = me.range.i2;
	}

//____________________________________________________________________________
/*
	template < typename TText, typename TSpec, typename TSpecFinder, typename TPosition >
	inline void 
	setPosition(Finder< Index<TText, TSpec>, TSpecFinder > & me, TPosition pos_)
	{
SEQAN_CHECKPOINT
		hostIterator(me) = me.range.i1 + pos_;
	}
*/
//____________________________________________________________________________

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline typename Position< Finder< Index<TText, TSpec>, TSpecFinder > >::Type
	position(Finder< Index<TText, TSpec>, TSpecFinder > & me)
	{
SEQAN_CHECKPOINT
		SEQAN_ASSERT(!empty(me))
		return *me.data_iterator;
	}

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline typename Position< Finder< Index<TText, TSpec>, TSpecFinder > >::Type
	position(Finder< Index<TText, TSpec>, TSpecFinder > const & me)
	{
SEQAN_CHECKPOINT
		SEQAN_ASSERT(!empty(me))
		return hostIterator(me) - begin(container(me), Rooted());
	}


//////////////////////////////////////////////////////////////////////////////
// find

	template < typename TText, typename TSpec, typename TSpecFinder, typename TPattern >
	inline bool find(
		Finder<Index<TText, TSpec>, TSpecFinder> &finder,
		TPattern const &pattern)
	{
		if (empty(finder)) {
			_findFirstIndex(finder, needle(pattern), TSpecFinder());
			hostIterator(finder) = finder.range.i1;
		} else
			++hostIterator(finder);
		return !atEnd(finder);
	}

	template < typename TText, typename TSpec, typename TSpecFinder >
	inline bool find(Finder<Index<TText, TSpec>, TSpecFinder> &finder)
	{
		if (empty(finder)) return false;
		++hostIterator(finder);
		return !atEnd(finder);
	}

}

#endif
