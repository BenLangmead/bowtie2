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
  $Id: segment_base.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEGMENT_BASE_H
#define SEQAN_HEADER_SEGMENT_BASE_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Segment
//////////////////////////////////////////////////////////////////////////////

/**
.Class.Segment:
..cat:Sequences
..summary:A contiguous part of a sequence.
..signature:Segment<THost, TSpec>
..param.THost:Type of the whole sequence.
...metafunction:Metafunction.Host
...text:Instances of $Segment<THost, TSpec>$ are subsequences of $THost$ objects.
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:@Spec.InfixSegment@.
*/

struct InfixSegment;

template <typename THost, typename TSpec = InfixSegment>
class Segment
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Host.param.T.type:Class.Segment

template <typename THost, typename TSpec>
struct Host<Segment<THost, TSpec> >
{
	typedef THost Type;
};

template <typename THost, typename TSpec>
struct Host<Segment<THost, TSpec> const >
{
	typedef THost Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.Segment

template <typename THost, typename TSpec>
struct Spec<Segment<THost, TSpec> >
{
	typedef TSpec Type;
};
template <typename THost, typename TSpec>
struct Spec<Segment<THost, TSpec> const>
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Segment

template <typename THost, typename TSpec>
struct Value<Segment<THost, TSpec> >
{
	typedef typename Value<THost>::Type Type;
};

template <typename THost, typename TSpec>
struct Value<Segment<THost, TSpec> const >
{
	typedef typename Value<THost const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.GetValue.param.T.type:Class.Segment

template <typename THost, typename TSpec>
struct GetValue<Segment<THost, TSpec> >
{
	typedef typename GetValue<THost>::Type Type;
};

template <typename THost, typename TSpec>
struct GetValue<Segment<THost, TSpec> const >
{
	typedef typename GetValue<THost const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Class.Segment

template <typename THost, typename TSpec>
struct Iterator<Segment<THost, TSpec>, Rooted>
{
	typedef Segment<THost, TSpec> TSequence;
	typedef typename Iterator<THost, Standard>::Type TIterator;
	typedef Iter<TSequence, AdaptorIterator<TIterator> > Type;
};
template <typename THost, typename TSpec>
struct Iterator<Segment<THost, TSpec> const, Rooted>
{
	typedef Segment<THost, TSpec> const TSequence;
	typedef typename Iterator<THost const, Standard>::Type TIterator;
	typedef Iter<TSequence, AdaptorIterator<TIterator> > Type;
};

template <typename THost, typename TSpec>
struct Iterator<Segment<THost, TSpec>, Standard>:
	Iterator<THost, Standard>
{
};
template <typename THost, typename TSpec>
struct Iterator<Segment<THost, TSpec> const, Standard>:
	Iterator<THost, Standard>
{
};



//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Size.param.T.type:Class.Segment

template <typename THost, typename TSpec>
struct Size<Segment<THost, TSpec> >
{
	typedef typename Size<THost>::Type Type;
};

template <typename THost, typename TSpec>
struct Size<Segment<THost, TSpec> const >
{
	typedef typename Size<THost>::Type Type;
};

template <typename THost, typename TSpec>
struct Position<Segment<THost, TSpec> >
{
	typedef typename Position<THost>::Type Type;
};

template <typename THost, typename TSpec>
struct Position<Segment<THost, TSpec> const >
{
	typedef typename Position<THost>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.DefaultOverflowImplicit.param.T.type:Class.Segment

template <typename THost, typename TSpec>
struct DefaultOverflowImplicit<Segment<THost, TSpec > >:
	DefaultOverflowImplicit<THost>
{
};

template <typename THost, typename TSpec>
struct DefaultOverflowImplicit<Segment<THost, TSpec > const >:
	DefaultOverflowImplicit<THost>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.DefaultOverflowExplicit.param.T.type:Class.Segment

template <typename THost, typename TSpec>
struct DefaultOverflowExplicit<Segment<THost, TSpec > >:
	DefaultOverflowExplicit<THost>
{
};

template <typename THost, typename TSpec>
struct DefaultOverflowExplicit<Segment<THost, TSpec > const >:
	DefaultOverflowExplicit<THost>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IsContiguous.param.T.type:Class.Segment

template <typename THost, typename TSpec>
struct IsContiguous< Segment<THost, TSpec> >:
	public IsContiguous<THost> {};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IsSequence.param.T.type:Class.Segment

template <typename THost, typename TSpec>
struct IsSequence< Segment<THost, TSpec> > {
    typedef True Type;
	enum { VALUE = true };
};

//////////////////////////////////////////////////////////////////////////////

///.Function.atBegin.param.iterator.type:Class.Segment
///.Function.atEnd.param.iterator.type:Class.Segment
///.Function.goBegin.param.iterator.type:Class.Segment
///.Function.goEnd.param.iterator.type:Class.Segment
///.Function.goNext.param.iterator.type:Class.Segment
///.Function.goPrevious.param.iterator.type:Class.Segment
///.Function.value.param.container.type:Class.Segment

///.Function.shareResources.param.sequence1, sequence2.type:Class.Segment

//////////////////////////////////////////////////////////////////////////////
// functions for all Segment classes
//////////////////////////////////////////////////////////////////////////////

///.Function.id.param.object.type:Class.Segment

template <typename THost, typename TSpec>
inline void const * 
id(Segment<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return id(host(me));
}

//////////////////////////////////////////////////////////////////////////////

///.Function.length.param.object.type:Class.Segment

template <typename THost, typename TSpec>
inline typename Size<Segment<THost, TSpec> const>::Type 
length(Segment<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return endPosition(me) - beginPosition(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.capacity.param.object.type:Class.Segment

template <typename THost, typename TSpec>
inline typename Size< Segment<THost, TSpec> const>::Type 
capacity(Segment<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return capacity(host(me)) + length(me) - length(host(me));
}

//////////////////////////////////////////////////////////////////////////////
// assign
//////////////////////////////////////////////////////////////////////////////

/**
.Function.assign:
..remarks:If $target$ is a @Class.Segment@ object, then
$limit$ denotes the maximal length of @Function.host.$host(target)$@ after the operation.
..param.target.type:Class.Segment
..param.source.type:Class.Segment
*/

//overload of binary version for strings: 

template<typename THost, typename TSpec, typename TSource>
inline void 
assign(Segment<THost, TSpec> & target, 
	  TSource & source)
{
SEQAN_CHECKPOINT
	typedef Segment<THost, TSpec> TTarget;
	assign(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename THost, typename TSpec, typename TSource>
inline void 
assign(Segment<THost, TSpec> & target, 
	  TSource const & source)
{
SEQAN_CHECKPOINT
	typedef Segment<THost, TSpec> TTarget;
	assign(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

//(for temporary targets)

template<typename THost, typename TSpec, typename TSource>
inline void 
assign(Segment<THost, TSpec> const & target, 
	  TSource & source)
{
SEQAN_CHECKPOINT
	typedef Segment<THost, TSpec> const TTarget;
	assign(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename THost, typename TSpec, typename TSource>
inline void 
assign(Segment<THost, TSpec> const & target, 
	  TSource const & source)
{
SEQAN_CHECKPOINT
	typedef Segment<THost, TSpec> const TTarget;
	assign(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct _Assign_Segment
{
	template <typename THost, typename TSpec, typename TSource>
	static inline void 
	assign_(
		Segment<THost, TSpec> & target, 
		TSource & source)
	{
SEQAN_CHECKPOINT
		if ((void *) &target == (void *) &source) return;

		typedef Segment<THost, TSpec> Target;

		replace(host(target), beginPosition(target), endPosition(target), source, TExpand());

		typename Iterator<Target, Standard>::Type new_end = begin(target, Standard()) + length(source);
		typename Iterator<THost, Standard>::Type host_end = end(host(target), Standard());
		if (new_end > host_end) new_end = host_end;
		setEnd(target, new_end);
	}

	template <typename THost, typename TSpec, typename TSource>
	static inline void 
	assign_(
		Segment<THost, TSpec> & target, 
		TSource & source, 
		typename Size< Segment<THost, TSpec> >::Type limit)
	{
SEQAN_CHECKPOINT
		if ((void *) &target == (void *) &source) return;

		typedef Segment<THost, TSpec> Target;

		replace(host(target), beginPosition(target), endPosition(target), source, limit, TExpand());

		typename Iterator<Target, Standard>::Type new_end = begin(target, Standard()) + length(source);
		typename Iterator<THost, Standard>::Type host_end = end(host(target), Standard());
		if (begin(target, Standard()) > host_end) setBegin(target, host_end);
		if (new_end > host_end) new_end = host_end;
		setEnd(target, new_end);
	}

	template <typename THost, typename TSpec, typename TSource>
	static inline void 
	assign_(
		Segment<THost, TSpec> const & target, 
		TSource & source)
	{
SEQAN_CHECKPOINT
		replace(host(target), beginPosition(target), endPosition(target), source, TExpand());
	}

	template <typename THost, typename TSpec, typename TSource>
	static inline void 
	assign_(
		Segment<THost, TSpec> const & target, 
		TSource & source, 
		typename Size< Segment<THost, TSpec> >::Type limit)
	{
SEQAN_CHECKPOINT
		replace(host(target), beginPosition(target), endPosition(target), source, limit, TExpand());
	}
};

//____________________________________________________________________________

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
assign(Segment<THost, TSpec> & target, 
	   TSource & source, 
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Assign_Segment<Tag<TExpand> const>::assign_(target, source);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
assign(Segment<THost, TSpec> & target, 
	   TSource const & source, 
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Assign_Segment<Tag<TExpand> const>::assign_(target, source);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
assign(Segment<THost, TSpec> & target, 
	   TSource & source, 
	   typename Size< Segment<THost, TSpec> >::Type limit, 
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Assign_Segment<Tag<TExpand> const>::assign_(target, source, limit);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
assign(Segment<THost, TSpec> & target, 
	   TSource const & source, 
	   typename Size< Segment<THost, TSpec> >::Type limit, 
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Assign_Segment<Tag<TExpand> const>::assign_(target, source, limit);
}

//(for temporary targets)

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
assign(Segment<THost, TSpec> const & target, 
	   TSource & source, 
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Assign_Segment<Tag<TExpand> const>::assign_(target, source);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
assign(Segment<THost, TSpec> const & target, 
	   TSource const & source, 
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Assign_Segment<Tag<TExpand> const>::assign_(target, source);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
assign(Segment<THost, TSpec> const & target, 
	   TSource & source, 
	   typename Size< Segment<THost, TSpec> >::Type limit, 
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Assign_Segment<Tag<TExpand> const>::assign_(target, source, limit);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
assign(Segment<THost, TSpec> const & target, 
	   TSource const & source, 
	   typename Size< Segment<THost, TSpec> >::Type limit, 
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Assign_Segment<Tag<TExpand> const>::assign_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
// move
//////////////////////////////////////////////////////////////////////////////

//overload of binary version: 

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
move(Segment<THost, TSpec> & target, 
	 TSource & source)
{
SEQAN_CHECKPOINT
	typedef Segment<THost, TSpec> TTarget;
	move(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
move(Segment<THost, TSpec> & target, 
	 TSource const & source)
{
SEQAN_CHECKPOINT
	typedef Segment<THost, TSpec> TTarget;
	move(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

//(for temporary targets)

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
move(Segment<THost, TSpec> const & target, 
	 TSource & source)
{
SEQAN_CHECKPOINT
	typedef Segment<THost, TSpec> const TTarget;
	move(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
move(Segment<THost, TSpec> const & target, 
	 TSource const & source)
{
SEQAN_CHECKPOINT
	typedef Segment<THost, TSpec> const TTarget;
	move(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}



//////////////////////////////////////////////////////////////////////////////
// append
//////////////////////////////////////////////////////////////////////////////

/**
.Function.append:
..remarks:If $target$ is a @Class.Segment@ object, then
$limit$ denotes the maximal length of @Function.host.$host(target)$@ after the operation.
..param.target.type:Class.Segment
..param.source.type:Class.Segment
*/


template <typename TExpand>
struct _Append_Sequence_2_Segment
{
	template <typename THost, typename TSpec, typename TSource>
	static inline void 
	append_(
		Segment<THost, TSpec> & target, 
		TSource & source)
	{
SEQAN_CHECKPOINT
		typedef Segment<THost, TSpec> Target;

		replace(host(target), endPosition(target), endPosition(target), source, TExpand());

		typename Iterator<Target, Standard>::Type new_end = end(target, Standard()) + length(source);
		typename Iterator<THost, Standard>::Type host_end = end(host(target), Standard());
		if (new_end > host_end) new_end = host_end;
		setEnd(target, new_end);
	}

	template <typename THost, typename TSpec, typename TSource>
	static inline void 
	append_(
		Segment<THost, TSpec> & target, 
		TSource const & source, 
		typename Size< Segment<THost, TSpec> >::Type limit)
	{
SEQAN_CHECKPOINT
		typedef Segment<THost, TSpec> Target;

		replace(host(target), endPosition(target), endPosition(target), source, limit, TExpand());
		typename Iterator<Target, Standard>::Type new_end = end(target, Standard()) + length(source);
		typename Iterator<THost, Standard>::Type host_end = end(host(target), Standard());
		if (begin(target) > host_end) setBegin(target, host_end);
		if (new_end > host_end) new_end = host_end;
		setEnd(target, new_end);
	}

	template <typename THost, typename TSpec, typename TSource>
	static inline void 
	append_(
		Segment<THost, TSpec> const & target, 
		TSource & source)
	{
SEQAN_CHECKPOINT
		replace(host(target), endPosition(target), endPosition(target), source, TExpand());
	}

	template <typename THost, typename TSpec, typename TSource>
	static inline void 
	append_(
		Segment<THost, TSpec> const & target, 
		TSource const & source, 
		typename Size< Segment<THost, TSpec> >::Type limit)
	{
SEQAN_CHECKPOINT
		replace(host(target), endPosition(target), endPosition(target), source, limit, TExpand()); //??? INSERT
	}
};
//____________________________________________________________________________


template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
append(
	Segment<THost, TSpec> & target, 
	TSource & source, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_Sequence_2_Segment<Tag<TExpand> const>::append_(target, source);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
append(
	Segment<THost, TSpec> & target, 
	TSource const & source, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_Sequence_2_Segment<Tag<TExpand> const>::append_(target, source);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
append(
	Segment<THost, TSpec> & target, 
	TSource & source, 
	typename Size< Segment<THost, TSpec> >::Type limit, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_Sequence_2_Segment<Tag<TExpand> const>::append_(target, source, limit);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
append(
	Segment<THost, TSpec> & target, 
	TSource const & source, 
	typename Size< Segment<THost, TSpec> >::Type limit, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_Sequence_2_Segment<Tag<TExpand> const>::append_(target, source, limit);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
append(
	Segment<THost, TSpec> const & target, 
	TSource & source, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_Sequence_2_Segment<Tag<TExpand> const>::append_(target, source);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
append(
	Segment<THost, TSpec> const & target, 
	TSource const & source, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_Sequence_2_Segment<Tag<TExpand> const>::append_(target, source);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
append(
	Segment<THost, TSpec> const & target, 
	TSource & source, 
	typename Size< Segment<THost, TSpec> >::Type limit, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_Sequence_2_Segment<Tag<TExpand> const>::append_(target, source, limit);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
append(
	Segment<THost, TSpec> const & target, 
	TSource const & source, 
	typename Size< Segment<THost, TSpec> >::Type limit, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_Sequence_2_Segment<Tag<TExpand> const>::append_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
// appendValue
//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct _Append_Value_2_Segment
{
	template <typename T, typename TValue>
	static inline void 
	appendValue_(T & me,
				TValue & _value)
	{
SEQAN_CHECKPOINT
		insertValue(host(me), endPosition(me), TExpand());
		if (endPosition(me) < length(host(me)) ) //this could be false for some TExpand
		{
			setEndPosition(me, endPosition(me) + 1);
		}
	}
};

//____________________________________________________________________________

template <typename THost, typename TSpec, typename TValue, typename TExpand>
inline void
appendValue(Segment<THost, TSpec> & me, 
			TValue const & _value,
			Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_Value_2_Segment<Tag<TExpand> const>::appendValue_(me, _value);
}
template <typename THost, typename TSpec, typename TValue, typename TExpand>
inline void
appendValue(Segment<THost, TSpec> const & me, 
			TValue const & _value,
			Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Append_Value_2_Segment<Tag<TExpand> const>::appendValue_(me, _value);
}

//////////////////////////////////////////////////////////////////////////////
// insertValue
//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct _Insert_Value_2_Segment
{
	template <typename T, typename TPosition, typename TValue>
	static inline void 
	insertValue_(T & me,
				TPosition pos,
				TValue & _value)
	{
SEQAN_CHECKPOINT
		insertValue(host(me), beginPosition(me) + pos, TExpand());
		if (endPosition(me) < length(host(me)) ) //this could be false for some TExpand
		{
			setEndPosition(me, endPosition(me) + 1);
		}
	}
};

//____________________________________________________________________________

template <typename THost, typename TSpec, typename TPosition, typename TValue, typename TExpand>
inline void
insertValue(Segment<THost, TSpec> & me, 
			TPosition pos,
			TValue const & _value,
			Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Insert_Value_2_Segment<Tag<TExpand> const>::insertValue_(me, pos, _value);
}
template <typename THost, typename TSpec, typename TPosition, typename TValue, typename TExpand>
inline void
insertValue(Segment<THost, TSpec> const & me, 
			TPosition pos,
			TValue const & _value,
			Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Insert_Value_2_Segment<Tag<TExpand> const>::insertValue_(me, pos, _value);
}

//////////////////////////////////////////////////////////////////////////////
// replace
//////////////////////////////////////////////////////////////////////////////

/**
.Function.replace:
..remarks:If $target$ is a @Class.Segment@ object, then
$limit$ denotes the maximal length of @Function.host.$host(target)$@ after the operation.
..param.target.type:Class.Segment
..param.source.type:Class.Segment
*/

template <typename TExpand>
struct _Replace_Sequence_2_Segment
{
	template <typename THost, typename TSpec, typename TSource>
	static inline void 
	replace_(
		Segment<THost, TSpec> & target, 
		typename Position< Segment<THost, TSpec> >::Type pos_begin,
		typename Position< Segment<THost, TSpec> >::Type pos_end,
		TSource & source)
	{
SEQAN_CHECKPOINT
		typedef Segment<THost, TSpec> Target;

		replace(host(target), beginPosition(target) + pos_begin, beginPosition(target) + pos_end, source, TExpand());

		typename Iterator<Target, Standard>::Type new_end = begin(target, Standard()) + length(target) - pos_end + pos_begin + length(source);
		typename Iterator<THost, Standard>::Type host_end = end(host(target), Standard());
		if (new_end > host_end) new_end = host_end;
		setEnd(target, new_end);
	}

	template <typename THost, typename TSpec, typename TSource>
	static inline void 
	replace_(
		Segment<THost, TSpec> & target, 
		typename Position< Segment<THost, TSpec> >::Type pos_begin,
		typename Position< Segment<THost, TSpec> >::Type pos_end,
		TSource & source, 
		typename Size< Segment<THost, TSpec> >::Type limit)
	{
SEQAN_CHECKPOINT
		typedef Segment<THost, TSpec> Target;

		replace(host(target), beginPosition(target) + pos_begin, beginPosition(target) + pos_end, source, limit, TExpand());

		typename Iterator<Target, Standard>::Type new_end = begin(target, Standard()) + length(target) - pos_end + pos_begin + length(source);
		typename Iterator<THost, Standard>::Type host_end = end(host(target), Standard());
		if (begin(target, Standard()) > host_end) setBegin(target, host_end);
		if (new_end > host_end) new_end = host_end;
		setEnd(target, new_end);
	}

	template <typename THost, typename TSpec, typename TSource>
	static inline void 
	replace_(
		Segment<THost, TSpec> const & target, 
		typename Position< Segment<THost, TSpec> const>::Type pos_begin,
		typename Position< Segment<THost, TSpec> const>::Type pos_end,
		TSource & source)
	{
SEQAN_CHECKPOINT
		replace(host(target), beginPosition(target) + pos_begin, beginPosition(target) + pos_end, source, TExpand());
	}

	template <typename THost, typename TSpec, typename TSource>
	static inline void 
	replace_(
		Segment<THost, TSpec> const & target, 
		typename Position< Segment<THost, TSpec> const>::Type pos_begin,
		typename Position< Segment<THost, TSpec> const>::Type pos_end,
		TSource & source, 
		typename Size< Segment<THost, TSpec> >::Type limit)
	{
SEQAN_CHECKPOINT
		replace(host(target), beginPosition(target) + pos_begin, beginPosition(target) + pos_end, source, limit, TExpand()); //??? INSERT
	}
};
//____________________________________________________________________________


template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
replace(
	Segment<THost, TSpec> & target, 
	typename Position< Segment<THost, TSpec> >::Type pos_begin,
	typename Position< Segment<THost, TSpec> >::Type pos_end,
	TSource & source, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Replace_Sequence_2_Segment<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
replace(
	Segment<THost, TSpec> & target, 
	typename Position< Segment<THost, TSpec> >::Type pos_begin,
	typename Position< Segment<THost, TSpec> >::Type pos_end,
	TSource const & source, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Replace_Sequence_2_Segment<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
replace(
	Segment<THost, TSpec> & target, 
	typename Position< Segment<THost, TSpec> >::Type pos_begin,
	typename Position< Segment<THost, TSpec> >::Type pos_end,
	TSource & source, 
	typename Size< Segment<THost, TSpec> >::Type limit, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Replace_Sequence_2_Segment<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
replace(
	Segment<THost, TSpec> & target, 
	typename Position< Segment<THost, TSpec> >::Type pos_begin,
	typename Position< Segment<THost, TSpec> >::Type pos_end,
	TSource const & source, 
	typename Size< Segment<THost, TSpec> >::Type limit, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Replace_Sequence_2_Segment<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
replace(
	Segment<THost, TSpec> const & target, 
	typename Position< Segment<THost, TSpec> const>::Type pos_begin,
	typename Position< Segment<THost, TSpec> const>::Type pos_end,
	TSource & source, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Replace_Sequence_2_Segment<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
replace(
	Segment<THost, TSpec> const & target, 
	typename Position< Segment<THost, TSpec> const>::Type pos_begin,
	typename Position< Segment<THost, TSpec> const>::Type pos_end,
	TSource const & source, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Replace_Sequence_2_Segment<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
replace(
	Segment<THost, TSpec> const & target, 
	typename Position< Segment<THost, TSpec> const>::Type pos_begin,
	typename Position< Segment<THost, TSpec> const>::Type pos_end,
	TSource & source, 
	typename Size< Segment<THost, TSpec> >::Type limit, 
	Tag<TExpand> const )
{
SEQAN_CHECKPOINT
	_Replace_Sequence_2_Segment<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void 
replace(
	Segment<THost, TSpec> const & target, 
	typename Position< Segment<THost, TSpec> const>::Type pos_begin,
	typename Position< Segment<THost, TSpec> const>::Type pos_end,
	TSource const & source, 
	typename Size< Segment<THost, TSpec> >::Type limit, 
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Replace_Sequence_2_Segment<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
// handling of iterators as begin and end

/*
template<typename THost, typename TTargetSpec, typename TSource, typename TExpand>
inline void 
replace(Segment<THost, TTargetSpec> & target,
		typename Iterator< Segment<THost, TTargetSpec>, Rooted>::Type pos_begin,
		typename Iterator< Segment<THost, TTargetSpec>, Rooted>::Type pos_end,
		TSource & source,
		Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	replace(target, position(pos_begin), position(pos_end), source, tag);
}

template<typename THost, typename TTargetSpec, typename TSource, typename TExpand>
inline void 
replace(Segment<THost, TTargetSpec> & target,
		typename Iterator< Segment<THost, TTargetSpec>, Rooted>::Type pos_begin,
		typename Iterator< Segment<THost, TTargetSpec>, Rooted>::Type pos_end,
		TSource & source,
		typename Size< Segment<THost, TTargetSpec> >::Type limit,
		Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	replace(target, position(pos_begin), position(pos_end), source, limit, tag);
}

template<typename THost, typename TTargetSpec, typename TSource, typename TExpand>
inline void 
replace(Segment<THost, TTargetSpec> const & target,
		typename Iterator< Segment<THost, TTargetSpec> const, Rooted>::Type pos_begin,
		typename Iterator< Segment<THost, TTargetSpec> const, Rooted>::Type pos_end,
		TSource & source,
		Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	replace(target, position(pos_begin), position(pos_end), source, tag);
}

template<typename THost, typename TTargetSpec, typename TSource, typename TExpand>
inline void 
replace(Segment<THost, TTargetSpec> const & target,
		typename Iterator< Segment<THost, TTargetSpec> const, Rooted>::Type pos_begin,
		typename Iterator< Segment<THost, TTargetSpec> const, Rooted>::Type pos_end,
		TSource & source,
		typename Size< Segment<THost, TTargetSpec> >::Type limit,
		Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	replace(target, position(pos_begin), position(pos_end), source, limit, tag);
}
*/

//////////////////////////////////////////////////////////////////////////////
///.Function.resize.param.object.type:Class.Segment

template <typename THost, typename TSpec, typename TExpand>
inline typename Size< Segment<THost, TSpec> >::Type 
resize(
	Segment<THost, TSpec> & me,
	typename Size< Segment<THost, TSpec> >::Type new_length,
	Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT

	typename Size<Segment<THost, TSpec> >::Type me_length = length(me);
	typename Position<THost>::Type me_end_pos = endPosition(me);
	if (new_length > me_length)
	{
		new_length = me_length + resizeSpace(host(me), new_length - me_length, me_end_pos, me_end_pos, tag);
	}
	else if (new_length < me_length)
	{
		new_length = resizeSpace(host(me), 0, me_end_pos - (me_length - new_length), me_end_pos, tag);
	}
	_setLength(me, new_length);
	return new_length;
}

//////////////////////////////////////////////////////////////////////////////
//??? TODO: fill (kopie von resize anpassen)

//////////////////////////////////////////////////////////////////////////////
///.Function.clear.param.object.type:Class.Segment

template <typename THost, typename TSpec>
inline void 
clear(Segment<THost, TSpec> & target)
{
SEQAN_CHECKPOINT
	assign(target, "");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TLeftSpec, typename TRight>
Segment<TLeftValue, TLeftSpec> const & 
operator += (Segment<TLeftValue, TLeftSpec> & left,
			 TRight const & right)
{
SEQAN_CHECKPOINT
	append(left, right);
	return left;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight >
inline bool
operator == (Segment<TLeftHost, TLeftSpec> const & left, 
			TRight const & right)
{
SEQAN_CHECKPOINT
	typename Comparator<Segment<TLeftHost, TLeftSpec> >::Type _lex(left, right);
    return isEqual(_lex);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight >
inline bool
operator != (Segment<TLeftHost, TLeftSpec> const & left, 
			TRight const & right)
{
SEQAN_CHECKPOINT
	typename Comparator<Segment<TLeftHost, TLeftSpec> >::Type _lex(left, right);
    return isNotEqual(_lex);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight>
inline bool
operator < (Segment<TLeftHost, TLeftSpec> const & left, 
			TRight const & right)
{
SEQAN_CHECKPOINT
	return isLess(left, right, typename DefaultPrefixOrder<Segment<TLeftHost, TLeftSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight>
inline bool
operator <= (Segment<TLeftHost, TLeftSpec> const & left, 
			 TRight const & right)
{
SEQAN_CHECKPOINT
	return isLessOrEqual(left, right, typename DefaultPrefixOrder<Segment<TLeftHost, TLeftSpec> >::Type());
}
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight>
inline bool
operator > (Segment<TLeftHost, TLeftSpec> const & left, 
			TRight const & right)
{
SEQAN_CHECKPOINT
	return isGreater(left, right, typename DefaultPrefixOrder<Segment<TLeftHost, TLeftSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight>
inline bool
operator >= (Segment<TLeftHost, TLeftSpec> const & left, 
		TRight const & right)
{
SEQAN_CHECKPOINT
	return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<Segment<TLeftHost, TLeftSpec> >::Type());
}


//////////////////////////////////////////////////////////////////////////////
// stream operators
//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename THost, typename TSpec>
inline TStream &
operator << (TStream & target, 
			 Segment<THost, TSpec> const & source)
{
SEQAN_CHECKPOINT
	write(target, source);
	return target;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename THost, typename TSpec>
inline TStream &
operator >> (TStream & source, 
			 Segment<THost, TSpec> & target)
{
SEQAN_CHECKPOINT
	read(source, target);
	return source;
}
template <typename TStream, typename THost, typename TSpec>
inline TStream &
operator >> (TStream & source, 
			 Segment<THost, TSpec> const & target)
{
SEQAN_CHECKPOINT
	read(source, target);
	return source;
}

//////////////////////////////////////////////////////////////////////////////


} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
