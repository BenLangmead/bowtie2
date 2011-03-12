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
  $Id: basic_host.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_HOST_H
#define SEQAN_HEADER_BASIC_HOST_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Host Functions
//////////////////////////////////////////////////////////////////////////////
//these functions assume that the hosted object exports a function "_dataHost"
//that returns a reference to a holder type of Host<T>::Type & 

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline bool
emptyHost(T const & me)
{
SEQAN_CHECKPOINT
	return empty(_dataHost(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline bool
dependentHost(T const & me)
{
SEQAN_CHECKPOINT
	return dependent(_dataHost(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline void
clearHost(T & me)
{
SEQAN_CHECKPOINT
	clear(_dataHost(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline void
createHost(T & me)
{
SEQAN_CHECKPOINT
	create(_dataHost(me));
}

//____________________________________________________________________________

template <typename T, typename THost>
inline void
createHost(T & me,
		   THost & host_)
{
SEQAN_CHECKPOINT
	create(_dataHost(me), host_);
}
template <typename T, typename THost>
inline void
createHost(T & me,
		   THost const & host_)
{
SEQAN_CHECKPOINT
	create(_dataHost(me), host_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, typename THost>
inline void
setHost(T & me,
		THost & host_)
{
SEQAN_CHECKPOINT
	setValue(_dataHost(me), host_);
}
template <typename T, typename THost>
inline void
setHost(T & me,
		THost const & host_)
{
SEQAN_CHECKPOINT
	setValue(_dataHost(me), host_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline typename Host<T>::Type &
host(T & me)
{
SEQAN_CHECKPOINT
	return value(_dataHost(me));
}
template <typename T>
inline typename Host<T const>::Type &
host(T const & me)
{
SEQAN_CHECKPOINT
	return value(_dataHost(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, typename THost>
inline void
assignHost(T & me,
		   THost & host_)
{
SEQAN_CHECKPOINT
	assignValue(_dataHost(me), host_);
}
template <typename T, typename THost>
inline void
assignHost(T & me,
		   THost const & host_)
{
SEQAN_CHECKPOINT
	assignValue(_dataHost(me), host_);
}
//////////////////////////////////////////////////////////////////////////////

template <typename T, typename THost>
inline void
moveHost(T & me,
		 THost & host_)
{
SEQAN_CHECKPOINT
	moveValue(_dataHost(me), host_);
}
template <typename T, typename THost>
inline void
moveHost(T & me,
		 THost const & host_)
{
SEQAN_CHECKPOINT
	moveValue(_dataHost(me), host_);
}

//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...


