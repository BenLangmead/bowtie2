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
  $Id: find_base.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_BASE_H
#define SEQAN_HEADER_FIND_BASE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.DefaultFinder:
..cat:Searching
..summary:Default @Class.Finder@ specialization type.
..signature:DefaultFinder<THaystack>::Type
..param.THaystack:The given haystack type.
..returns:Is $void$ by default and @Tag.Index Find Algorithm.ESA_FIND_MLR@ if $THaystack$ is an @Class.Index@.
*/
	template < typename TObject >
	struct DefaultFinder {
		typedef void Type;
	};

/**
.Metafunction.DefaultPattern:
..cat:Searching
..summary:Default @Class.Pattern@ specialization type.
..signature:DefaultPattern<TNeedle>::Type
..param.TNeedle:The given needle type.
..returns:Is $void$ by default.
*/
	template < typename TObject >
	struct DefaultPattern {
		typedef void Type;
	};

/**
.Metafunction.Haystack:
..summary:Returns the haystack type of a @Class.Finder@ type.
..cat:Searching
..signature:Haystack<TFinder>::Type
..param.TFinder:A @Class.Finder@ type.
...type:Class.Finder
..returns:The haystack type of $TFinder$, i.e. $THaystack$ for $Finder<THaystack, TSpec>$.
*/

	template <typename TFinder>
	struct Haystack {
		typedef typename Container<TFinder>::Type Type;
	};

/**
.Metafunction.Needle:
..summary:Returns the needle type of a @Class.Pattern@ type.
..cat:Searching
..signature:Needle<TPattern>::Type
..param.TPattern:A @Class.Pattern@ type.
...type:Class.Pattern
..returns:The needle type of $TPattern$, i.e. $TNeedle$ for $Pattern<TNeedle, TSpec>$.
*/

	template <typename TPattern>
	struct Needle {
		typedef typename Host<TPattern>::Type Type;
	};

	template <typename THost, typename TSpec>
	struct Needle<Segment<THost, TSpec> > {
		typedef Segment<THost, TSpec> Type;
	};

	template <typename THost, typename TSpec>
	struct Needle<Segment<THost, TSpec> const> {
		typedef Segment<THost, TSpec> const Type;
	};

//////////////////////////////////////////////////////////////////////////////

/**
.Class.Pattern:
..summary:Holds the needle and preprocessing data (depends on algorithm).
..cat:Searching
..signature:Pattern<TNeedle[, TSpec]>
..param.TNeedle:The needle type.
...type:Class.String
..param.TSpec:The online-algorithm to search with.
...remarks:Leave empty for index-based pattern matching (see @Class.Index@).
...default:The result of @Metafunction.DefaultPattern@
..remarks:If $TNeedle$ is a set of strings, then $position(pattern)$ returns the index of the currently matching needle.
*/

template < typename TNeedle, typename TSpec = typename DefaultPattern<TNeedle>::Type >
class Pattern;

//default implementation
template < typename TNeedle >
class Pattern<TNeedle, void>
{
public:
	typedef typename Position<TNeedle>::Type TNeedlePosition;

	Holder<TNeedle> data_host;
	TNeedlePosition data_begin_position;
	TNeedlePosition data_end_position;

	Pattern() {}

	template <typename _TNeedle>
	Pattern(_TNeedle & ndl):
		data_host(ndl) {}

	template <typename _TNeedle>
	Pattern(_TNeedle const & ndl):
		data_host(ndl) {}

};
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
inline Holder<TNeedle> & 
_dataHost(Pattern<TNeedle, TSpec> & me) 
{ 
	return me.data_host;
}
template <typename TNeedle, typename TSpec>
inline Holder<TNeedle> & 
_dataHost(Pattern<TNeedle, TSpec> const & me) 
{
	return const_cast<Holder<TNeedle> &>(me.data_host);
}

//host access: see basic_host.h


//???TODO: Diese Funktion entfernen! (sobald setHost bei anderen pattern nicht mehr eine Art "assignHost" ist)
template <typename TNeedle, typename TSpec, typename TNeedle2>
inline void 
setHost(Pattern<TNeedle, TSpec> & me,
		TNeedle2 const & ndl) 
{
	 me.data_host = ndl; //assign => Pattern haelt eine Kopie => doof!
}
template <typename TNeedle, typename TSpec, typename TNeedle2>
inline void 
setHost(Pattern<TNeedle, TSpec> & me,
		TNeedle2 & ndl) 
{ 
	 me.data_host = ndl; //assign => Pattern haelt eine Kopie => doof!
}
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
inline typename Position<Pattern<TNeedle, TSpec> >::Type & 
beginPosition(Pattern<TNeedle, TSpec> & me) 
{
	return me.data_begin_position;
}
template <typename TNeedle, typename TSpec>
inline typename Position<Pattern<TNeedle, TSpec> const >::Type & 
beginPosition(Pattern<TNeedle, TSpec> const & me) 
{
	return me.data_begin_position;
}


template <typename TNeedle, typename TSpec, typename TPosition>
inline void
setBeginPosition(Pattern<TNeedle, TSpec> & me, 
				 TPosition _pos) 
{
	me.data_begin_position = _pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
inline typename Position<Pattern<TNeedle, TSpec> >::Type & 
endPosition(Pattern<TNeedle, TSpec> & me) 
{
	return me.data_end_position;
}
template <typename TNeedle, typename TSpec>
inline typename Position<Pattern<TNeedle, TSpec> const >::Type & 
endPosition(Pattern<TNeedle, TSpec> const & me) 
{
	return me.data_end_position;
}

template <typename TNeedle, typename TSpec, typename TPosition>
inline void
setEndPosition(Pattern<TNeedle, TSpec> & me, 
			   TPosition _pos) 
{
	me.data_end_position = _pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
inline typename Infix<TNeedle>::Type 
segment(Pattern<TNeedle, TSpec> & me) 
{
	typedef typename Infix<TNeedle>::Type TInfix;
	return TInfix(host(me), me.data_begin_position, me.data_end_position);
}
template <typename TNeedle, typename TSpec>
inline typename Infix<TNeedle>::Type 
segment(Pattern<TNeedle, TSpec> const & me) 
{
	typedef typename Infix<TNeedle>::Type TInfix;
	return TInfix(host(me), me.data_begin_position, me.data_end_position);
}

//////////////////////////////////////////////////////////////////////////////


/**
.Function.needle:
..summary:Returns the needle of a @Class.Pattern@ object (not implemented for some online-algorithms).
..cat:Searching
..signature:needle(pattern)
..param.pattern:The @Class.Pattern@ object to search with.
...type:Class.Pattern
..returns:The needle object to search for.
..remarks:The result type is @Metafunction.Needle@$<TPattern>::Type$ for pattern of type $TPattern$.
*/

template < typename TObject >
inline typename Needle<TObject>::Type &
needle(TObject &obj) 
{
	return obj;
}

template < typename TObject >
inline typename Needle<TObject const>::Type &
needle(TObject const &obj) 
{
	return obj;
}


///.Function.position.param.iterator.type:Class.Pattern

template < typename TNeedle, typename TSpec >
inline typename Needle< Pattern<TNeedle, TSpec> >::Type &
needle(Pattern<TNeedle, TSpec> & obj) 
{
	return host(obj);
}

template < typename TNeedle, typename TSpec >
inline typename Needle< Pattern<TNeedle, TSpec> const>::Type &
needle(Pattern<TNeedle, TSpec> const & obj) 
{
	return host(obj);
}

/**
.Function.setNeedle:
..summary:Sets the needle of a @Class.Pattern@ object and optionally induces preprocessing.
..cat:Searching
..signature:setNeedle(pattern, needle)
..param.pattern:The @Class.Pattern@ object to search with.
...type:Class.Pattern
..param.needle:The needle object to search for.
...type:Class.String
*/

	template < typename TNeedle, typename TSpec >
	inline void
	setNeedle(Pattern<TNeedle, TSpec> &obj, TNeedle const &ndl) {
		setHost(obj, ndl);
	}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.find:
..summary:Search for a @Class.Pattern@ in a @Class.Finder@ object.
..cat:Searching
..signature:find(finder, pattern)
..signature:find(finder, pattern, k)
..param.finder:The @Class.Finder@ object to search through.
...remarks:For online-algorithm $patterns$, finder can also be an arbitrary @Concept.Rooted Iterator@.
...type:Class.Finder
...type:Concept.Rooted Iterator
..param.pattern:The @Class.Pattern@ object to search for.
...remarks:For index $finders$, pattern can also be a Sequence.
...type:Class.Pattern
..param.k:Desired minimal score (for approximate matching).
...remarks:$k$ has to be a number <= 0.
...remarks:Differences are deletions, insertions and substitutions.
..returns:$boolean$ that indicates whether an occurence of $pattern$ was found or not.
..remarks:Repeated calls of this function iterate through all occurences of $pattern$.
*/

/**
.Class.Finder:
..summary:Holds the haystack and a current search context.
..cat:Searching
..signature:Finder<THaystack[, TSpec]>
..param.THaystack:The haystack type.
...type:Class.String
...type:Class.Index
..param.TSpec:The index-algorithm to search with (Optional).
...default:The result of @Metafunction.DefaultFinder@
...remarks:Leave empty for online pattern matching (see @Class.Pattern@).
...remarks:If $THaystack$ is an @Class.Index@, then $TSpec$ specifies the index search algorithm.
..remarks:$position(finder)$ returns the position of the current hit in the haystack.
If $THaystack$ is a set of strings or an index of a set of strings, then $position(finder)$ returns a @Class.Pair@ $(hayNo, pos)$,
in which $hayNo$ is the haystack index and $pos$ the local position of the hit.
..remarks:Use $clear(finder)$ to reset a finder object and search from the beginning.
*/

///.Function.clear.param.object.type:Class.Finder
///.Function.position.param.iterator.type:Class.Finder

	template < typename THaystack, typename TSpec = typename DefaultFinder<THaystack>::Type >
	class Finder
	{
		typedef typename Iterator<THaystack, Rooted>::Type TIterator;

	public:
		TIterator data_iterator;
		bool _needReinit;					// if true, the Pattern needs to be reinitialized

		Finder():
			_needReinit(true) {}

		Finder(THaystack &haystack):
			data_iterator(begin(haystack, Rooted())),
			_needReinit(true) {}

		Finder(TIterator &iter):
			data_iterator(iter),
			_needReinit(true) {}

		Finder(TIterator const &iter):
			data_iterator(iter),
			_needReinit(true) {}

		Finder(Finder const &orig):
			data_iterator(orig.data_iterator),
			_needReinit(orig._needReinit) {};

//____________________________________________________________________________

		inline typename Reference<TIterator>::Type 
		operator* () 
		{
SEQAN_CHECKPOINT
			return value(hostIterator(*this));
		}

		inline typename Reference<TIterator const>::Type 
		operator* () const
		{
SEQAN_CHECKPOINT
			return value(hostIterator(*this));
		}

//____________________________________________________________________________

		operator TIterator () const
		{
SEQAN_CHECKPOINT
			return data_iterator;
		}

//____________________________________________________________________________

	};



//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline typename _Parameter<THaystack>::Type 
	host(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		return container(hostIterator(me));
	}

	template <typename THaystack, typename TSpec>
	inline typename _Parameter<THaystack>::Type 
	host(Finder<THaystack, TSpec> const & me)
	{
SEQAN_CHECKPOINT
		return container(hostIterator(me));
	}

	template <typename THaystack, typename TSpec>
	inline typename _Parameter<THaystack>::Type 
	container(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		return container(hostIterator(me));
	}

	template <typename THaystack, typename TSpec>
	inline typename _Parameter<THaystack>::Type 
	container(Finder<THaystack, TSpec> const & me)
	{
SEQAN_CHECKPOINT
		return container(hostIterator(me));
	}

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline void
	setHost(Finder<THaystack, TSpec> & me, typename _Parameter<THaystack>::Type container_)
	{
SEQAN_CHECKPOINT
		setContainer(hostIterator(me), container_);
		goBegin(me);
	}

	template <typename THaystack, typename TSpec>
	inline void
	setContainer(Finder<THaystack, TSpec> & me, typename _Parameter<THaystack>::Type container_)
	{
SEQAN_CHECKPOINT
		setContainer(hostIterator(me), container_);
		goBegin(me);
	}

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline typename Iterator<THaystack, Rooted>::Type &
	hostIterator(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		return me.data_iterator;
	}

	template <typename THaystack, typename TSpec>
	inline typename Iterator<THaystack, Rooted>::Type const &
	hostIterator(Finder<THaystack, TSpec> const & me)
	{
SEQAN_CHECKPOINT
		return me.data_iterator;
	}

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline bool
	empty(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		return me._needReinit;
	}

	template <typename THaystack, typename TSpec>
	inline void
	clear(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		me._needReinit = true;
	}

//____________________________________________________________________________

	template <typename T>
	inline void
	_finderSetNonEmpty(T & me)
	{
SEQAN_CHECKPOINT
		goBegin(me);
	}


	template <typename THaystack, typename TSpec>
	inline void
	_finderSetNonEmpty(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		me._needReinit = false;
	}

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline bool
	atBegin(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		return (!empty(me) && atBegin(hostIterator(me)));
	}

	template <typename THaystack, typename TSpec>
	inline bool
	atEnd(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		return (!empty(me) && atEnd(hostIterator(me)));
	}

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline void
	goBegin(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		//_finderSetNonEmpty(me);
		goBegin(hostIterator(me));
	}

	template <typename THaystack, typename TSpec>
	inline void
	goEnd(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		//_finderSetNonEmpty(me);
		goEnd(hostIterator(me));
	}

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline typename Position<Finder<THaystack, TSpec> >::Type
	position(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		if (empty(me)) return 0;
		return position(hostIterator(me));
	}

	template <typename THaystack, typename TSpec>
	inline typename Position<Finder<THaystack, TSpec> >::Type
	position(Finder<THaystack, TSpec> const & me)
	{
SEQAN_CHECKPOINT
		if (empty(me)) return 0;
		return position(hostIterator(me));
	}

//____________________________________________________________________________
/**
.Function.setPosition:
..cat:Searching
..summary:Sets the position of a finder.
..signature:setPosition(finder, pos)
..param.finder:A finder.
...class:Class.Finder
..param.pos:A position.
...metafunction:Metafunction.Position
..see:Function.position
*/

	template <typename THaystack, typename TSpec, typename TPosition>
	inline void 
	setPosition(Finder<THaystack, TSpec> & me, TPosition pos_)
	{
SEQAN_CHECKPOINT
		setPosition(hostIterator(me), pos_);
	}

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline Finder<THaystack, TSpec> &
	operator--(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		--hostIterator(me);
		return me;
	}

	template <typename THaystack, typename TSpec>
	inline Finder<THaystack, TSpec> &
	operator++(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
/*			if (beforeBegin()) {
			goBegin(hostIterator(me));
		} else*/
			++hostIterator(me);
		return me;
	}

//////////////////////////////////////////////////////////////////////////////
// operator +
//////////////////////////////////////////////////////////////////////////////

	template <typename THaystack, typename TSpec, typename TIntegral>
	inline Finder<THaystack, TSpec> const
	operator + (Finder<THaystack, TSpec> const & left, TIntegral right)
	{
SEQAN_CHECKPOINT
		return Finder<THaystack, TSpec>(hostIterator(left) + right);
	}

//////////////////////////////////////////////////////////////////////////////
// operator +=
//////////////////////////////////////////////////////////////////////////////

	template <typename THaystack, typename TSpec, typename TIntegral>
	inline Finder<THaystack, TSpec> &
	operator += (Finder<THaystack, TSpec> & left,
					TIntegral right)
	{
SEQAN_CHECKPOINT
		hostIterator(left) += right;
		return left;
	}

//////////////////////////////////////////////////////////////////////////////
// operator -
//////////////////////////////////////////////////////////////////////////////

	template <typename THaystack, typename TSpec, typename TIntegral>
	inline Finder<THaystack, TSpec> const
	operator - (Finder<THaystack, TSpec> const & left, TIntegral right)
	{
SEQAN_CHECKPOINT
		return Finder<THaystack, TSpec>(hostIterator(left) - right);
	}

	template <typename THaystack, typename TSpec, typename TIntegral>
	inline typename Difference<Finder<THaystack, TSpec> const>::Type
	operator - (Finder<THaystack, TSpec> const & left, Finder<THaystack, TSpec> const & right)
	{
SEQAN_CHECKPOINT
		return hostIterator(left) - hostIterator(right);
	}

//////////////////////////////////////////////////////////////////////////////
// operator -=
//////////////////////////////////////////////////////////////////////////////

	template <typename THaystack, typename TSpec, typename TIntegral>
	inline Finder<THaystack, TSpec> &
	operator -= (Finder<THaystack, TSpec> & left,
					TIntegral right)
	{
SEQAN_CHECKPOINT
		hostIterator(left) -= right;
		return left;
	}

//____________________________________________________________________________


/**
.Function.setHaystack:
..summary:Sets the haystack of a @Class.Finder@ object.
..cat:Searching
..signature:setHaystack(finder, haystack)
..param.finder:The @Class.Finder@ object to search with.
...type:Class.Finder
..param.haystack:The haystack object the finder searches through.
...type:Class.String
*/

	template < typename THaystack, typename TSpec >
	inline void
	setHaystack(Finder<THaystack, TSpec> &obj, THaystack const &hstk) {
		setHost(obj, hstk);
	}

/**
.Function.haystack:
..summary:Returns the haystack of a @Class.Finder@ object.
..cat:Searching
..signature:haystack(finder)
..param.finder:The @Class.Finder@ object to search through.
...type:Class.Finder
..returns:The haystack object.
..remarks:The result type is @Metafunction.Haystack@$<TFinder>::Type$ for finder of type $TFinder$.
*/

	template < typename TObject >
	inline typename Haystack<TObject>::Type &
	haystack(TObject &obj) {
		return container(obj);
	}

	template < typename TObject >
	inline typename Haystack<TObject const>::Type &
	haystack(TObject const &obj) {
		return container(obj);
	}


//////////////////////////////////////////////////////////////////////////////


	template <typename THaystack, typename TSpec>
	struct Container< Finder<THaystack, TSpec> > {
		typedef THaystack Type;
	};

	template <typename THaystack, typename TSpec>
	struct Container< Finder<THaystack, TSpec> const> {
		typedef THaystack const Type;
	};

	template <typename THaystack, typename TSpec>
	struct Host< Finder<THaystack, TSpec> > {
		typedef THaystack Type;
	};

	template <typename THaystack, typename TSpec>
	struct Host< Finder<THaystack, TSpec> const> {
		typedef THaystack const Type;
	};


	template <typename THaystack, typename TSpec>
	struct Value< Finder<THaystack, TSpec> > {
		typedef typename Value<THaystack>::Type Type;
	};

	template <typename THaystack, typename TSpec>
	struct Position< Finder<THaystack, TSpec> >:
		Position<THaystack> {};

	template <typename THaystack, typename TSpec>
	struct Difference< Finder<THaystack, TSpec> > {
		typedef typename Difference<THaystack>::Type Type;
	};

	template <typename THaystack, typename TSpec>
	struct Size< Finder<THaystack, TSpec> > {
		typedef typename Size<THaystack>::Type Type;
	};

//____________________________________________________________________________


	template <typename TNeedle, typename TSpec>
	struct Container< Pattern<TNeedle, TSpec> > {
		typedef TNeedle Type;
	};

	template <typename TNeedle, typename TSpec>
	struct Container< Pattern<TNeedle, TSpec> const > {
		typedef TNeedle const Type;
	};

	template <typename TNeedle, typename TSpec>
	struct Host< Pattern<TNeedle, TSpec> > {
		typedef TNeedle Type;
	};

	template <typename TNeedle, typename TSpec>
	struct Host< Pattern<TNeedle, TSpec> const > {
		typedef TNeedle const Type;
	};


	template <typename TPattern, typename TSpec>
	struct Value< Pattern<TPattern, TSpec> > {
		typedef typename Value<TPattern>::Type Type;
	};

	template <typename TPattern, typename TSpec>
	struct Position< Pattern<TPattern, TSpec> > {
		typedef typename Position<TPattern>::Type Type;
	};

	template <typename TPattern, typename TSpec>
	struct Difference< Pattern<TPattern, TSpec> > {
		typedef typename Difference<TPattern>::Type Type;
	};

	template <typename TPattern, typename TSpec>
	struct Size< Pattern<TPattern, TSpec> > {
		typedef typename Size<TPattern>::Type Type;
	};


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
