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
  $Id: basic_iterator_base.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_ITERATOR_BASE_H
#define SEQAN_HEADER_BASIC_ITERATOR_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Iter
//////////////////////////////////////////////////////////////////////////////
/**
.Class.Iter:
..cat:Basic
..summary:Iterator that is used to traverse containers.
..signature:Iter<TContainer, TSpec>
..param.TContainer:Type of the container that can be iterated by $Iter$.
...metafunction:Metafunction.Container
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
..implements:Concept.Iterator
*/
template <typename TContainer, typename TSpec>
class Iter;

//////////////////////////////////////////////////////////////////////////////
///.Metafunction.Spec.param.T.type:Class.Iter

template <typename TContainer, typename TSpec>
struct Spec<Iter<TContainer, TSpec> >
{
	typedef TSpec Type;
};
template <typename TContainer, typename TSpec>
struct Spec<Iter<TContainer, TSpec> const>
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Iter

template <typename TContainer, typename TSpec>
struct Value<Iter<TContainer, TSpec> >:
	Value<TContainer>
{
};
template <typename TContainer, typename TSpec>
struct Value<Iter<TContainer, TSpec> const>:
	Value<TContainer>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.GetValue.param.T.type:Class.Iter

template <typename TContainer, typename TSpec>
struct GetValue<Iter<TContainer, TSpec> >:
	GetValue<TContainer>
{
};
template <typename TContainer, typename TSpec>
struct GetValue<Iter<TContainer, TSpec> const>:
	GetValue<TContainer>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Reference.param.T.type:Class.Iter

template <typename TContainer, typename TSpec>
struct Reference<Iter<TContainer, TSpec> >:
	Reference<TContainer>
{
};
template <typename TContainer, typename TSpec>
struct Reference<Iter<TContainer, TSpec> const>:
	Reference<TContainer>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Container.param.T.type:Class.Iter

template <typename T> struct Container;

template <typename TContainer, typename TSpec>
struct Container<Iter<TContainer, TSpec> >
{
	typedef TContainer Type;
};
template <typename TContainer, typename TSpec>
struct Container<Iter<TContainer, TSpec> const>
{
	typedef TContainer Type;
};

//////////////////////////////////////////////////////////////////////////////

/*
///.Metafunction.Host.param.T.type:Class.Iter

template <typename TContainer, typename TSpec>
struct Host<Iter<TContainer, TSpec> >:
	Container<Iter<TContainer, TSpec> >
{
};
template <typename TContainer, typename TSpec>
struct Host<Iter<TContainer, TSpec> const>:
	Container<Iter<TContainer, TSpec> const>
{
};
*/

//////////////////////////////////////////////////////////////////////////////
// operator *
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline typename Reference<Iter<TContainer, TSpec> >::Type 
operator * (Iter<TContainer, TSpec> & me)
{
SEQAN_CHECKPOINT
	return value(me);
}
template <typename TContainer, typename TSpec>
inline typename Reference<Iter<TContainer, TSpec> const>::Type 
operator * (Iter<TContainer, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return value(me);
}

//////////////////////////////////////////////////////////////////////////////
// operator ++
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline Iter<TContainer, TSpec> const &
operator ++ (Iter<TContainer, TSpec> & me)
{
SEQAN_CHECKPOINT
	goNext(me);
	return me;
}

template <typename TContainer, typename TSpec>
inline Iter<TContainer, TSpec> const
operator ++ (Iter<TContainer, TSpec> & me, int)
{
SEQAN_CHECKPOINT
	Iter<TContainer, TSpec> temp_(me);
	goNext(me);
	return temp_;
}

//////////////////////////////////////////////////////////////////////////////
// operator --
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline Iter<TContainer, TSpec> const &
operator -- (Iter<TContainer, TSpec> & me)
{
SEQAN_CHECKPOINT
	goPrevious(me);
	return me;
}

template <typename TContainer, typename TSpec>
inline Iter<TContainer, TSpec> const
operator -- (Iter<TContainer, TSpec> & me, int)
{
SEQAN_CHECKPOINT
	Iter<TContainer, TSpec> temp_(me);
	goPrevious(me);
	return temp_;
}

//////////////////////////////////////////////////////////////////////////////
// position
//////////////////////////////////////////////////////////////////////////////

//most Iter classes are rooted strings

template <typename TContainer, typename TSpec, typename TContainer2>
inline typename Position<Iter<TContainer, TSpec> const>::Type 
position(Iter<TContainer, TSpec> const & me,
		 TContainer2 const &)
{
SEQAN_CHECKPOINT
	return position(me);
}


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
