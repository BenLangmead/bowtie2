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
  $Id: find_multi.h,v 1.1 2008/08/25 16:20:06 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_MULTI_H
#define SEQAN_HEADER_FIND_MULTI_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

struct _MultipatternFinder;
typedef Tag<_MultipatternFinder> MultipatternFinder;
	
//____________________________________________________________________________

template <typename THaystack>
class Finder<THaystack, MultipatternFinder>
{
//____________________________________________________________________________
private:
	unsigned int data_pattern;

public:
	Finder():
		data_pattern(0)
	{
SEQAN_CHECKPOINT
	}

	Finder(Finder const & other_):
		data_pattern(other_.data_pattern)
	{
SEQAN_CHECKPOINT
	}

	~Finder()
	{
SEQAN_CHECKPOINT
	}
//____________________________________________________________________________

	Finder & 
	operator = (Finder const & other_)
	{
SEQAN_CHECKPOINT
		data_pattern = other_.data_pattern;
		return *this;
	}
//____________________________________________________________________________

	friend inline unsigned int &
	needle(Finder & me)
	{
SEQAN_CHECKPOINT
		return me.data_pattern;
	}
	friend inline unsigned int const &
	needle(Finder const & me)
	{
SEQAN_CHECKPOINT
		return me.data_pattern;
	}
//____________________________________________________________________________

	friend inline void
	setNeedle(Finder & me, unsigned int const needleIndex_)
	{
SEQAN_CHECKPOINT
		me.data_pattern = needleIndex_;
	}

//____________________________________________________________________________

	friend inline void
	init(Finder & me)
	{
SEQAN_CHECKPOINT
		me.data_pattern = 0;
	}
//____________________________________________________________________________


//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TNeedle>
friend inline bool
find(Finder & me,
	 THaystack & hstk,
	 TNeedle const & ndl)
{
SEQAN_CHECKPOINT
	while ( needle(me) < length(ndl) )
	{
		Finder<THaystack, Horspool> horspool(ndl[needle(me)]);
		bool found = find(horspool, hstk, ndl[needle(me)]);
		if (found)
		{
			return true;
		}
		setPosition(hstk, 0);
		++needle(me);
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////
/*
template <typename THaystack, typename TNeedle>
bool
findNext(Finder & me,
		 THaystack & hstk,
		 TNeedle const & ndl)
{
SEQAN_CHECKPOINT
	++hstk;
	return find(me, hstk, ndl);
}*/

//////////////////////////////////////////////////////////////////////////////

};

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
