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
  $Id: sequence_multiple.h,v 1.2 2009/02/19 01:51:23 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_STRING_SET_H
#define SEQAN_HEADER_STRING_SET_H

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// StringSet specs
	//////////////////////////////////////////////////////////////////////////////

	//struct Generous;						// large string-string, large id-string
	//struct Tight;							// small string-string, small id-string
	//struct Optimal??						// small string-string, large id-string

/**
.Spec.Dependent:
..summary:A string set storing references of the strings.
..cat:Sequences
..general:Class.StringSet
..signature:StringSet<TString, Dependent<TSpec> >
..param.TString:The string type.
...type:Class.String
..param.TSpec:The specializing type for the dependent string set.
...remarks:Possible values are $Tight$ or $Generous$
...remarks:$Tight$ is very space efficient whereas $Generous$ provides fast access to the strings in the container via ids.
*/

	// Default id holder string set
	template<typename TSpec = Generous >
	struct Dependent;						// holds references of its elements


/**
.Spec.ConcatDirect:
..summary:A string set storing the concatenation of all strings within one string.
..cat:Sequences
..general:Spec.Owner
..signature:StringSet<TString, Owner<ConcatDirect<> > >
..param.TString:The string type.
...type:Class.String
..remarks:The strings are internally stored in a $TString$ object and the character position type is a
a single integer value between 0 and the sum of string lengths minus 1.
..remarks:The position type can be returned or modified by the meta-function @Metafunction.SAValue@ called with the @Class.StringSet@ type.
*/

	template < typename TDelimiter = void >
	struct ConcatDirect;					// contains 1 string (the concatenation of n strings)
	//struct Default;							// contains n strings in a string

/**
.Spec.Owner:
..summary:A string set storing the strings as members.
..cat:Sequences
..general:Class.StringSet
..signature:StringSet<TString, Owner<> >
..signature:StringSet<TString, Owner<Default> >
..param.TString:The string type.
...type:Class.String
..remarks:The strings are internally stored in a $String<TString>$ object and the character position type is a
@Class.Pair@ $(seqNo,seqOfs)$ where seqNo identifies the string within the stringset and seqOfs identifies the position within this string.
..remarks:The position type can be returned or modified by the meta-function @Metafunction.SAValue@ called with the @Class.StringSet@ type.
*/
	template < typename TSpec = Default >
	struct Owner;							// owns its elements


    //////////////////////////////////////////////////////////////////////////////
	// Forwards
	//////////////////////////////////////////////////////////////////////////////

	template < typename TString, typename TSpec = Owner<> >
	class StringSet;

    template <typename TObject>
	struct Concatenator {
		typedef TObject Type;
	};

    template <typename TObject>
	struct Concatenator<TObject const> {
		typedef typename Concatenator<TObject>::Type const Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// StringSet limits
	//////////////////////////////////////////////////////////////////////////////

	template <typename TString>
	struct StringSetLimits {
		typedef Nothing Type;
	};

	template <typename TString>
	struct StringSetLimits<TString const> {
		typedef typename StringSetLimits<TString>::Type const Type;
	};

	template <typename TString>
	struct StringSetPosition {
		typedef typename Size<TString>::Type	Type;
	};

	template <typename TString, typename TSpec>
	struct StringSetLimits< StringSet<TString, TSpec> > {
		typedef typename Size<TString>::Type	TSize;
		typedef String<TSize>					Type;
	};

	template <typename TString, typename TSpec>
	struct StringSetPosition< StringSet<TString, TSpec> > {
		typedef typename Size<TString>::Type	TSize;
		typedef Pair<TSize>						Type;
	};



	//////////////////////////////////////////////////////////////////////////////
	// get StringSet limits
	//////////////////////////////////////////////////////////////////////////////

/**
.Function.stringSetLimits:
..cat:Sequences
..summary:Retrieves a string of delimiter positions of a @Class.StringSet@ which is needed for local<->global position conversions.
..signature:stringSetLimits(me)
..param.me:A string or string set.
...type:Class.String
...type:Class.StringSet
..returns:A reference to a string.
...remarks:If $me$ is a @Class.StringSet@ then the returned string is of size $length(me)+1$ and contains the ascending (virtual) delimiter positions of the concatenation of all strings in the string set.
...remarks:If $me$ is a @Class.String@, @Tag.Nothing@ is returned.
*/

    template <typename TStringSet>
	inline typename StringSetLimits<TStringSet>::Type
	stringSetLimits(TStringSet &) {
		return typename StringSetLimits<TStringSet>::Type();
	}

    template <typename TString, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSpec> >::Type &
	stringSetLimits(StringSet<TString, TSpec> &stringSet) {
		if (!_validStringSetLimits(stringSet))
			_refreshStringSetLimits(stringSet);
		return stringSet.limits;
	}

    template <typename TString, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSpec> const>::Type &
	stringSetLimits(StringSet<TString, TSpec> const &stringSet) {
		if (!_validStringSetLimits(stringSet))
			_refreshStringSetLimits(const_cast< StringSet<TString, TSpec>& >(stringSet));
		return stringSet.limits;
	}

	//////////////////////////////////////////////////////////////////////////////
	// accessing local positions
	//////////////////////////////////////////////////////////////////////////////

/**
.Function.getSeqNo:
..cat:Sequences
..summary:Returns the sequence number of a position.
..signature:getSeqNo(pos[, limits])
..param.pos:A position.
...type:Class.Pair
..param.limits:The limits string returned by @Function.stringSetLimits@.
..returns:A single integer value that identifies the string within the stringset $pos$ points at.
...remarks:If $limits$ is omitted or @Tag.Nothing@ $getSeqNo$ returns 0.
...remarks:If $pos$ is a local position (of class @Class.Pair@) then $i1$ is returned.
...remarks:If $pos$ is a global position (integer type and $limits$ is a @Class.String@) then $pos$ is converted to a local position and $i1$ is returned.
*/

/**
.Function.getSeqOffset:
..cat:Sequences
..summary:Returns the local sequence offset of a position.
..signature:getSeqOffset(pos[, limits])
..param.pos:A position.
...type:Class.Pair
..param.limits:The limits string returned by @Function.stringSetLimits@.
..returns:A single integer value that identifies the position within the string $pos$ points at.
...remarks:If $limits$ is omitted or @Tag.Nothing@ $getSeqNo$ returns $pos$.
...remarks:If $pos$ is a local position (of class @Class.Pair@) then $i2$ is returned.
...remarks:If $pos$ is a global position (integer type and $limits$ is a @Class.String@) then $pos$ is converted to a local position and $i2$ is returned.
*/

	//////////////////////////////////////////////////////////////////////////////
	// 1 sequence
	//////////////////////////////////////////////////////////////////////////////

	template <typename TPosition>
	inline TPosition getSeqNo(TPosition const &, Nothing const &) {
		return 0;
	}

	template <typename TPosition>
	inline TPosition getSeqOffset(TPosition const &pos, Nothing const &) {
		return pos;
	}
//____________________________________________________________________________

	template <typename TPosition>
	inline TPosition getSeqNo(TPosition const &) {
		return 0;
	}

	template <typename TPosition>
	inline TPosition getSeqOffset(TPosition const &pos) {
		return pos;
	}

	//////////////////////////////////////////////////////////////////////////////
	// n sequences (position type is Pair)
	//////////////////////////////////////////////////////////////////////////////

	template <typename T1, typename T2, typename TCompression, typename TLimitsString>
	inline T1 getSeqNo(Pair<T1, T2, TCompression> const &pos, TLimitsString const &) {
		return getValueI1(pos);
	}

	template <typename T1, typename T2, typename TCompression>
	inline T1 getSeqNo(Pair<T1, T2, TCompression> const &pos) {
		return getValueI1(pos);
	}
//____________________________________________________________________________

	template <typename T1, typename T2, typename TCompression, typename TLimitsString>
	inline T2 getSeqOffset(Pair<T1, T2, TCompression> const &pos, TLimitsString const &) {
		return getValueI2(pos);
	}

	template <typename T1, typename T2, typename TCompression>
	inline T1 getSeqOffset(Pair<T1, T2, TCompression> const &pos) {
		return getValueI2(pos);
	}

	//////////////////////////////////////////////////////////////////////////////
	// n sequences (position type is an integral type)
	//////////////////////////////////////////////////////////////////////////////

	template <typename TPos, typename TLimitsString>
	inline TPos getSeqNo(TPos const &pos, TLimitsString const &limits) {
		typedef typename Iterator<TLimitsString const, Standard>::Type TIter;
		typedef typename Value<TLimitsString>::Type TSize;
        TIter _begin = begin(limits, Standard());
        TIter _upper = ::std::upper_bound(_begin, end(limits, Standard()), (TSize)pos) - 1;
        return difference(_begin, _upper);
	}

	template <typename TPos, typename TLimitsString>
	inline TPos getSeqOffset(TPos const &pos, TLimitsString const &limits) {
		typedef typename Iterator<TLimitsString const, Standard>::Type TIter;
		typedef typename Value<TLimitsString>::Type TSize;
        TIter _begin = begin(limits, Standard());
        TIter _upper = ::std::upper_bound(_begin, end(limits, Standard()), (TSize)pos) - 1;
        return pos - *_upper;
	}

	//////////////////////////////////////////////////////////////////////////////
	// local -> global conversions
	//////////////////////////////////////////////////////////////////////////////

/**
.Function.posGlobalize:
..cat:Sequences
..summary:Converts a local/global to a global position.
..signature:posGlobalize(pos, limits)
..param.pos:A local or global position (pair or integer value).
...type:Class.Pair
..param.limits:The limits string returned by @Function.stringSetLimits@.
..returns:The corresponding global position of $pos$.
...remarks:If $pos$ is an integral type $pos$ is returned.
...remarks:If not, $limits[getSeqNo(pos, limits)] + getSeqOffset(pos, limits)$ is returned.
*/

	// any_position and no limits_string -> any_position
	template <typename TPosition>
	inline TPosition posGlobalize(TPosition const &pos, Nothing const &) {
		return pos;
	}

	// local_position (0,x) and no limits_string -> global_position x
	template <typename T1, typename T2, typename TCompression>
	inline T2 posGlobalize(Pair<T1, T2, TCompression> const &pos, Nothing const &) {
		return getSeqOffset(pos);
	}

	// any_position and no limits_string -> any_position
	template <typename TLimitsString, typename TPosition>
	inline TPosition posGlobalize(TPosition const &pos, TLimitsString const &) {
		return pos;
	}

	// local_position and limits_string -> global_position
	template <typename TLimitsString, typename T1, typename T2, typename TCompression>
	inline typename Value<TLimitsString>::Type
	posGlobalize(Pair<T1, T2, TCompression> const &pos, TLimitsString const &limits) {
		return limits[getSeqNo(pos, limits)] + getSeqOffset(pos, limits);
	}

	//////////////////////////////////////////////////////////////////////////////
	// global -> local position
	//////////////////////////////////////////////////////////////////////////////

/**
.Function.posLocalize:
..cat:Sequences
..summary:Converts a local/global to a local position.
..signature:posLocalize(result, pos, limits)
..param.pos:A local or global position (pair or integer value).
...type:Class.Pair
..param.limits:The limits string returned by @Function.stringSetLimits@.
..param.result:Reference to the resulting corresponding local position of $pos$.
...remarks:If $pos$ is an integral type and $limits$ is omitted or @Tag.Nothing@, $pos$ is returned.
...remarks:If $pos$ is a local position (of class @Class.Pair@) then $pos$ is returned.
...remarks:If $pos$ is a global position (integer type and $limits$ is a @Class.String@) then $pos$ is converted to a local position.
*/

	// any_position and no limits_string -> any_position
	template <typename TResult, typename TPosition>
	inline void posLocalize(TResult &result, TPosition const &pos, Nothing const &) {
		result = pos;
	}

	template <typename T1, typename T2, typename TCompression, typename TPosition>
	inline void posLocalize(Pair<T1, T2, TCompression> &result, TPosition const &pos, Nothing const &) {
		result.i1 = 0;
		result.i2 = pos;
	}

	// global_position and limits_string -> local_position
	template <typename TResult, typename TSize, typename TSpec, typename TPosition>
	inline void posLocalize(TResult &result, TPosition const &pos, String<TSize, TSpec> const &limits) {
		typedef typename Iterator<String<TSize> const, Standard>::Type TIter;
        TIter _begin = begin(limits, Standard());
		TIter _upper = ::std::upper_bound(_begin, end(limits, Standard()), (TSize)pos) - 1;
        result.i1 = difference(_begin, _upper);
        result.i2 = pos - *_upper;
	}

	// local_position -> local_position
	template <typename TResult, typename TSize, typename TSpec, typename T1, typename T2, typename TCompression>
	inline void posLocalize(TResult &result, Pair<T1, T2, TCompression> const &pos, String<TSize, TSpec> const &/*limits*/) {
		result = pos;
	}

	//////////////////////////////////////////////////////////////////////////////
	// prefix
	//////////////////////////////////////////////////////////////////////////////

	template < typename TString, typename TSpec, typename TPosition >
	inline typename Prefix<TString>::Type
	prefix(StringSet< TString, TSpec > &me, TPosition pos)
	{
		typedef StringSet<TString, TSpec>				TStringSet;
		typedef typename Size<TStringSet>::Type			TSetSize;
		typedef typename Size<TString>::Type			TStringSize;
		typedef Pair<TSetSize, TStringSize, Compressed>	TPair;

		TPair lPos;
		posLocalize(lPos, pos, stringSetLimits(me));
		return prefix(me[getSeqNo(lPos)], getSeqOffset(lPos));
	}

	template < typename TString, typename TSpec, typename TPosition >
	inline typename Prefix<TString const>::Type
	prefix(StringSet< TString, TSpec > const &me, TPosition pos)
	{
		typedef StringSet<TString, TSpec>				TStringSet;
		typedef typename Size<TStringSet>::Type			TSetSize;
		typedef typename Size<TString>::Type			TStringSize;
		typedef Pair<TSetSize, TStringSize, Compressed>	TPair;

		TPair lPos;
		posLocalize(lPos, pos, stringSetLimits(me));
		return prefix(me[getSeqNo(lPos)], getSeqOffset(lPos));
	}

	template < typename TString, typename TDelimiter, typename TPosition >
	inline typename Infix<TString>::Type
	prefix(StringSet< TString, Owner<ConcatDirect<TDelimiter> > > &me, TPosition pos)
	{
		return infix(me.concat, stringSetLimits(me)[getSeqNo(pos, stringSetLimits(me))], posGlobalize(pos, stringSetLimits(me)));
	}

	template < typename TString, typename TDelimiter, typename TPosition >
	inline typename Infix<TString const>::Type
	prefix(StringSet< TString, Owner<ConcatDirect<TDelimiter> > > const &me, TPosition pos)
	{
		return infix(me.concat, stringSetLimits(me)[getSeqNo(pos, stringSetLimits(me))], posGlobalize(pos, stringSetLimits(me)));
	}

	//////////////////////////////////////////////////////////////////////////////
	// suffix
	//////////////////////////////////////////////////////////////////////////////

	template < typename TString, typename TSpec, typename TPosition >
	inline typename Suffix<TString>::Type
	suffix(StringSet< TString, TSpec > &me, TPosition pos)
	{
		typedef StringSet<TString, TSpec>				TStringSet;
		typedef typename Size<TStringSet>::Type			TSetSize;
		typedef typename Size<TString>::Type			TStringSize;
		typedef Pair<TSetSize, TStringSize, Compressed>	TPair;

		TPair lPos;
		posLocalize(lPos, pos, stringSetLimits(me));
		return suffix(me[getSeqNo(lPos)], getSeqOffset(lPos));
	}

	template < typename TString, typename TSpec, typename TPosition >
	inline typename Suffix<TString const>::Type
	suffix(StringSet< TString, TSpec > const &me, TPosition pos)
	{
		typedef StringSet<TString, TSpec>				TStringSet;
		typedef typename Size<TStringSet>::Type			TSetSize;
		typedef typename Size<TString>::Type			TStringSize;
		typedef Pair<TSetSize, TStringSize, Compressed>	TPair;

		TPair lPos;
		posLocalize(lPos, pos, stringSetLimits(me));
		return suffix(me[getSeqNo(lPos)], getSeqOffset(lPos));
	}

	template < typename TString, typename TDelimiter, typename TPosition >
	inline typename Infix<TString>::Type
	suffix(StringSet< TString, Owner<ConcatDirect<TDelimiter> > > &me, TPosition pos)
	{
		return infix(me.concat, posGlobalize(pos, stringSetLimits(me)), stringSetLimits(me)[getSeqNo(pos, stringSetLimits(me)) + 1]);
	}

	template < typename TString, typename TDelimiter, typename TPosition >
	inline typename Infix<TString const>::Type
	suffix(StringSet< TString, Owner<ConcatDirect<TDelimiter> > > const &me, TPosition pos)
	{
		return infix(me.concat, posGlobalize(pos, stringSetLimits(me)), stringSetLimits(me)[getSeqNo(pos, stringSetLimits(me)) + 1]);
	}

	//////////////////////////////////////////////////////////////////////////////
	// infixWithLength
	//////////////////////////////////////////////////////////////////////////////

	template < typename TString, typename TSpec, typename TPosition, typename TSize >
	inline typename Infix<TString>::Type
	infixWithLength(StringSet< TString, TSpec > &me, TPosition pos, TSize length)
	{
		typedef StringSet<TString, TSpec>				TStringSet;
		typedef typename Size<TStringSet>::Type			TSetSize;
		typedef typename Size<TString>::Type			TStringSize;
		typedef Pair<TSetSize, TStringSize, Compressed>	TPair;

		TPair lPos;
		posLocalize(lPos, pos, stringSetLimits(me));
		return infixWithLength(me[getSeqNo(lPos)], getSeqOffset(lPos), length);
	}

	template < typename TString, typename TSpec, typename TPosition, typename TSize >
	inline typename Infix<TString const>::Type
	infixWithLength(StringSet< TString, TSpec > const &me, TPosition pos, TSize length)
	{
		typedef StringSet<TString, TSpec>				TStringSet;
		typedef typename Size<TStringSet>::Type			TSetSize;
		typedef typename Size<TString>::Type			TStringSize;
		typedef Pair<TSetSize, TStringSize, Compressed>	TPair;

		TPair lPos;
		posLocalize(lPos, pos, stringSetLimits(me));
		return infixWithLength(me[getSeqNo(lPos)], getSeqOffset(lPos), length);
	}

	template < typename TString, typename TDelimiter, typename TPosition, typename TSize >
	inline typename Infix<TString>::Type
	infixWithLength(StringSet< TString, Owner<ConcatDirect<TDelimiter> > > &me, TPosition pos, TSize length)
	{
		return infixWithLength(me.concat, posGlobalize(pos, stringSetLimits(me)), length);
	}

	template < typename TString, typename TDelimiter, typename TPosition, typename TSize >
	inline typename Infix<TString const>::Type
	infixWithLength(StringSet< TString, Owner<ConcatDirect<TDelimiter> > > const &me, TPosition pos, TSize length)
	{
		return infixWithLength(me.concat, posGlobalize(pos, stringSetLimits(me)), length);
	}

	//////////////////////////////////////////////////////////////////////////////
	// position arithmetics
	//////////////////////////////////////////////////////////////////////////////

	// posAtFirstLocal
	template <typename TPos, typename TLimitsString>
	inline bool posAtFirstLocal(TPos pos, TLimitsString const &limits) {
		return getSeqOffset(pos, limits) == 0;
	}
	template <typename TPos>
	inline bool posAtFirstLocal(TPos pos) {
		return getSeqOffset(pos) == 0;
	}


	// posPrev
	template <typename TPos>
	inline TPos posPrev(TPos pos) {
		return pos - 1;
	}

	template <typename T1, typename T2, typename TCompression>
	inline Pair<T1, T2, TCompression> posPrev(Pair<T1, T2, TCompression> const &pos) {
		return Pair<T1, T2, TCompression>(getValueI1(pos), getValueI2(pos) - 1);
	}

	// posNext
	template <typename TPos>
	inline TPos posNext(TPos pos) {
		return pos + 1;
	}

	template <typename T1, typename T2, typename TCompression>
	inline Pair<T1, T2, TCompression>
	posNext(Pair<T1, T2, TCompression> const &pos) {
		return Pair<T1, T2, TCompression>(getValueI1(pos), getValueI2(pos) + 1);
	}

	// posAdd
	template <typename TPos, typename TDelta>
	inline TPos posAdd(TPos pos, TDelta delta) {
		return pos + delta;
	}

	template <typename T1, typename T2, typename TCompression, typename TDelta>
	inline Pair<T1, T2, TCompression>
	posAdd(Pair<T1, T2, TCompression> const &pos, TDelta delta) {
		return Pair<T1, T2, TCompression>(getValueI1(pos), getValueI2(pos) + delta);
	}

	// posSub
	template <typename TA, typename TB>
	inline TA posSub(TA a, TB b) {
		return a - b;
	}

	template <
		typename TA1, typename TA2, typename TACompression,
		typename TB1, typename TB2, typename TBCompression
	>
	inline TA2
	posSub(Pair<TA1, TA2, TACompression> const &a, Pair<TB1, TB2, TBCompression> const &b) {
		return getValueI2(a) - getValueI2(b);
	}

	//////////////////////////////////////////////////////////////////////////////
	// position relations
	//////////////////////////////////////////////////////////////////////////////

	template <typename TPos>
	inline bool posLess(TPos const &a, TPos const &b) {
		return a < b;
	}

	template <typename T1, typename T2, typename TCompression>
	inline bool posLess(Pair<T1, T2, TCompression> const &a, Pair<T1, T2, TCompression> const &b) {
		return
			 (getValueI1(a) <  getValueI1(b)) ||
			((getValueI1(a) == getValueI1(b)) && (getValueI2(a) < getValueI2(b)));
	}

	template <typename TPos>
	inline int posCompare(TPos const &a, TPos const &b) {
		if (a < b) return -1;
		if (a > b) return 1;
		return 0;
	}

	template <typename T1, typename T2, typename TCompression>
	inline int posCompare(Pair<T1, T2, TCompression> const &a, Pair<T1, T2, TCompression> const &b) {
		if (getValueI1(a) < getValueI1(b)) return -1;
		if (getValueI1(a) > getValueI1(b)) return 1;
		return posCompare(getValueI2(a), getValueI2(b));
	}

	//////////////////////////////////////////////////////////////////////////////

	template <typename TPos, typename TString>
	inline typename Size<TString>::Type
	suffixLength(TPos pos, TString const &string) {
		return length(string) - pos;
	}

	template <typename TPos, typename TString, typename TSpec>
	inline typename Size<TString>::Type
	suffixLength(TPos pos, StringSet<TString, TSpec> const &stringSet) {
		return length(stringSet[getSeqNo(pos, stringSetLimits(stringSet))]) - getSeqOffset(pos, stringSetLimits(stringSet));
	}

	//////////////////////////////////////////////////////////////////////////////

	template <typename TString>
	inline unsigned
	countSequences(TString const &) {
		return 1;
	}

	template <typename TString, typename TSpec>
	inline typename Size<StringSet<TString, TSpec> >::Type
	countSequences(StringSet<TString, TSpec> const &stringSet) {
		return length(stringSet);
	}

	//////////////////////////////////////////////////////////////////////////////

	template <typename TSeqNo, typename TString>
	inline typename Size<TString>::Type
	sequenceLength(TSeqNo /*seqNo*/, TString const &string) {
		return length(string);
	}

	template <typename TSeqNo, typename TString, typename TSpec>
	inline typename Size<StringSet<TString, TSpec> >::Type
	sequenceLength(TSeqNo seqNo, StringSet<TString, TSpec> const &stringSet) {
		return length(stringSet[seqNo]);
	}

	//////////////////////////////////////////////////////////////////////////////
    // StringSet Container
    //////////////////////////////////////////////////////////////////////////////

/**
.Class.StringSet:
..cat:Sequences
..summary:A container class for a set of strings.
..signature:StringSet<TString, TSpec>
..param.TString:The string type.
...type:Class.String
..param.TSpec:The specializing type for the StringSet.
...metafunction:Metafunction.Spec
...default:$Generous$.
..include:sequence.h
*/

	//////////////////////////////////////////////////////////////////////////////
    // StringSet with individual sequences in a tight string of string pointers and corr. IDs
	template <typename TString>
	class StringSet<TString, Dependent<Tight> >
	{
		public:
			typedef String<TString*>							TStrings;
			typedef typename Id<StringSet>::Type				TIdType;
			typedef String<TIdType>								TIds;
			typedef typename StringSetLimits<StringSet>::Type	TLimits;
			typedef typename Concatenator<StringSet>::Type		TConcatenator;
	//____________________________________________________________________________

			TStrings		strings;
			TIds			ids;
			TLimits			limits;
			bool			limitsValid;		// is true if limits contains the cumulative sum of the sequence lengths
			TConcatenator	concat;
	//____________________________________________________________________________

			StringSet():
				limitsValid(true)
			{
			SEQAN_CHECKPOINT
				appendValue(limits, 0);
				concat.set = this;
			}

			template <typename TDefault>
			StringSet(StringSet<TString, Owner<TDefault> > const& _other) :
				limitsValid(true)
			{
				SEQAN_CHECKPOINT
				appendValue(limits, 0);
				concat.set = this;
				for(unsigned int i = 0; i<length(_other); ++i) appendValue(*this, _other[i]);
			}

			template <typename TPos>
			inline typename Reference<StringSet>::Type
			operator [] (TPos pos)
			{
		SEQAN_CHECKPOINT
				return value(*this, pos);
			}

			template <typename TPos>
			inline typename Reference<StringSet const>::Type
			operator [] (TPos pos) const
			{
				return value(*this, pos);
			}
	};

	//////////////////////////////////////////////////////////////////////////////
    // StringSet with individual sequences in a string of string pointers
	template <typename TString>
	class StringSet<TString, Dependent<Generous> >
	{
		public:
			typedef String<TString*>							TStrings;
			typedef typename Size<StringSet>::Type				TSize;
			typedef typename StringSetLimits<StringSet>::Type	TLimits;
			typedef typename Concatenator<StringSet>::Type		TConcatenator;
	//____________________________________________________________________________

			TStrings		strings;
			TLimits			limits;
			bool			limitsValid;		// is true if limits contains the cumulative sum of the sequence lengths
			TConcatenator	concat;
	//____________________________________________________________________________

			StringSet():
				limitsValid(true)
			{
			SEQAN_CHECKPOINT
				appendValue(limits, 0);
				concat.set = this;
			}

			template <typename TDefault>
			StringSet(StringSet<TString, Owner<TDefault> > const& _other) :
				limitsValid(true)
			{
				SEQAN_CHECKPOINT
				appendValue(limits, 0);
				concat.set = this;
				for(unsigned int i = 0; i<length(_other); ++i) appendValue(*this, _other[i]);
			}

			template <typename TPos>
			inline typename Reference<StringSet>::Type
			operator [] (TPos pos)
			{
				return value(*this, pos);
			}

			template <typename TPos>
			inline typename Reference<StringSet const>::Type
			operator [] (TPos pos) const
			{
				return value(*this, pos);
			}
	};

	//////////////////////////////////////////////////////////////////////////////
    // StringSet with individual sequences in a string of strings
    template < typename TString >
    class StringSet< TString, Owner<Default> >
    {
	public:

        typedef String<TString>								TStrings;
		typedef typename StringSetLimits<StringSet>::Type	TLimits;
		typedef typename Concatenator<StringSet>::Type		TConcatenator;
//____________________________________________________________________________

		TStrings		strings;
        TLimits			limits;
		bool			limitsValid;		// is true if limits contains the cumulative sum of the sequence lengths
		TConcatenator	concat;
//____________________________________________________________________________

		StringSet():
			limitsValid(true)
		{
			appendValue(limits, 0);
			concat.set = this;
		};
//____________________________________________________________________________

		template <typename TPos>
		inline typename Reference<StringSet>::Type
		operator [] (TPos pos)
		{
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<StringSet const>::Type
		operator [] (TPos pos) const
		{
			return value(*this, pos);
		}
	};


	//////////////////////////////////////////////////////////////////////////////
    // StringSet with directly concatenated sequences
    template < typename TString, typename TDelimiter >
    class StringSet< TString, Owner<ConcatDirect<TDelimiter> > >
    {
	public:
		typedef typename StringSetLimits<StringSet>::Type	TLimits;
		typedef typename Concatenator<StringSet>::Type		TConcatenator;
//____________________________________________________________________________

		TLimits			limits;
		TConcatenator	concat;
//____________________________________________________________________________

		StringSet()	{
			appendValue(limits, 0);
		}
//____________________________________________________________________________

		template <typename TPos>
		inline typename Reference<StringSet>::Type
		operator [] (TPos pos)
		{
		SEQAN_CHECKPOINT
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<StringSet const>::Type
		operator [] (TPos pos) const
		{
		SEQAN_CHECKPOINT
			return value(*this, pos);
		}
	};


//////////////////////////////////////////////////////////////////////////////
// meta functions

	template < typename TString, typename TSpec >
    struct Value< StringSet< TString, TSpec > > {
        typedef TString Type;
    };

	template < typename TString, typename TSpec >
    struct Value< StringSet< TString, TSpec > const> {
        typedef TString Type;
    };

	template < typename TString, typename TSpec >
	struct Size< StringSet< TString, TSpec > >:
		Size< typename StringSetLimits< StringSet<TString, TSpec> >::Type > {};

	template < typename TString, typename TSpec >
	struct Prefix< StringSet< TString, TSpec > >:
		Prefix< TString > {};

	template < typename TString, typename TSpec >
	struct Prefix< StringSet< TString, TSpec > const >:
		Prefix< TString const > {};

	template < typename TString, typename TSpec >
	struct Suffix< StringSet< TString, TSpec > >:
		Suffix< TString > {};

	template < typename TString, typename TSpec >
	struct Suffix< StringSet< TString, TSpec > const >:
		Suffix< TString const > {};

	template < typename TString, typename TSpec >
	struct Infix< StringSet< TString, TSpec > >:
		Infix< TString > {};

	template < typename TString, typename TSpec >
	struct Infix< StringSet< TString, TSpec > const >:
		Infix< TString const > {};


	// direct concatenation
	template < typename TString, typename TSpec >
	struct Value< StringSet< TString, Owner<ConcatDirect<TSpec> > > >:
		Infix<TString> {};

    template < typename TString, typename TSpec >
	struct GetValue< StringSet< TString, Owner<ConcatDirect<TSpec> > > >:
		Infix<TString> {};

    template < typename TString, typename TSpec >
	struct GetValue< StringSet< TString, Owner<ConcatDirect<TSpec> > > const >:
		Infix<TString const> {};

    template < typename TString, typename TSpec >
	struct Reference< StringSet< TString, Owner<ConcatDirect<TSpec> > > >:
		Infix<TString> {};

    template < typename TString, typename TSpec >
	struct Reference< StringSet< TString, Owner<ConcatDirect<TSpec> > > const >:
		Infix<TString const> {};

    template <typename TString, typename TSpec>
	struct AllowsFastRandomAccess< StringSet< TString, TSpec > >:
		AllowsFastRandomAccess<TString> {};

	template < typename TString, typename TSpec >
	struct Prefix< StringSet< TString, Owner<ConcatDirect<TSpec> > > >:
		Infix< TString > {};

	template < typename TString, typename TSpec >
	struct Prefix< StringSet< TString, Owner<ConcatDirect<TSpec> > > const >:
		Infix< TString const > {};

	template < typename TString, typename TSpec >
	struct Suffix< StringSet< TString, Owner<ConcatDirect<TSpec> > > >:
		Infix< TString > {};

	template < typename TString, typename TSpec >
	struct Suffix< StringSet< TString, Owner<ConcatDirect<TSpec> > > const >:
		Infix< TString const > {};

	template < typename TString, typename TSpec >
	struct Infix< StringSet< TString, Owner<ConcatDirect<TSpec> > > >:
		Infix< TString > {};

	template < typename TString, typename TSpec >
	struct Infix< StringSet< TString, Owner<ConcatDirect<TSpec> > > const >:
		Infix< TString const > {};


//////////////////////////////////////////////////////////////////////////////
// validStringSetLimits

	template < typename T >
    inline bool _validStringSetLimits(T const &) {
        return true;
    }

	template < typename TString, typename TSpec >
    inline bool _validStringSetLimits(StringSet< TString, TSpec > const &me) {
        return me.limitsValid;
    }

	template < typename TString, typename TSpec >
    inline bool _validStringSetLimits(StringSet< TString, Owner<ConcatDirect<TSpec> > > const &) {
        return true;
    }

//////////////////////////////////////////////////////////////////////////////
// _refreshStringSetLimits

	template < typename T >
	inline void _refreshStringSetLimits(T &) {}

	template < typename TString, typename TSpec >
    inline void _refreshStringSetLimits(StringSet< TString, Owner<ConcatDirect<TSpec> > > &) {}

	template < typename TString, typename TSpec >
    inline void _refreshStringSetLimits(StringSet< TString, TSpec > &me)
	{
		typedef StringSet< TString, TSpec >					TStringSet;
		typedef typename StringSetLimits<TStringSet>::Type	TLimits;

		typename Value<TLimits>::Type	sum = 0;
		typename Size<TStringSet>::Type	len = length(me);
		typename Size<TStringSet>::Type	i = 0;

//		SEQAN_ASSERT(length(me.limits) == len + 1);
//		resize(me.limits, len + 1);
		for(; i < len; ++i) {
			me.limits[i] = sum;
			sum += length(me[i]);
		}
		me.limits[i] = sum;
		me.limitsValid = true;
    }

//////////////////////////////////////////////////////////////////////////////
// find the i-th non-zero value of a string me

	template < typename TValue, typename TSpec, typename TPos >
	inline typename Size< String<TValue, TSpec> >::Type
	_findIthNonZeroValue(String<TValue, TSpec> const &me, TPos i)
	{
		typename Iterator< String<TValue, TSpec> const, Standard >::Type it = begin(me, Standard());
		typename Iterator< String<TValue, TSpec> const, Standard >::Type itEnd = end(me, Standard());

		for(; it != itEnd; ++it)
			if (*it) {
				if (i) {
					--i;
				}
				else {
					return position(it, me);
				}
			}
		return length(me);
	}

//////////////////////////////////////////////////////////////////////////////
// count non-zero values before position i

	template < typename TValue, typename TSpec, typename TPos >
	inline typename Size< String<TValue, TSpec> >::Type
	_countNonZeroValues(String<TValue, TSpec> const &me, TPos i)
	{
		typename Iterator< String<TValue, TSpec> const, Standard >::Type it = begin(me, Standard());
		typename Iterator< String<TValue, TSpec> const, Standard >::Type itEnd = begin(me, Standard()) + i;
		typename Size< String<TValue, TSpec> >::Type counter = 0;

		for(; it != itEnd; ++it)
			if (*it) ++counter;
		return counter;
	}

//////////////////////////////////////////////////////////////////////////////
// lengthSum

	template < typename TString >
    inline typename Size<TString>::Type lengthSum(TString const &me) {
        return length(me);
    }

	template < typename TString, typename TSpec >
    inline typename Size<TString>::Type lengthSum(StringSet< TString, TSpec > const &me) {
        return back(stringSetLimits(me));
    }

///.Function.appendValue.param.target.type:Class.StringSet
//////////////////////////////////////////////////////////////////////////////
// appendValue

	// Default
	template < typename TString, typename TString2, typename TExpand >
    inline void appendValue(
		StringSet< TString, Owner<Default> > &me,
		TString2 const &obj,
		Tag<TExpand> const)
	{
        appendValue(me.strings, obj);
        appendValue(me.limits, lengthSum(me) + length(obj));
    }

	// ConcatDirect
	template < typename TString, typename TString2, typename TExpand >
    inline void appendValue(
		StringSet< TString, Owner<ConcatDirect<void> > > &me,
		TString2 const &obj,
		Tag<TExpand> const)
	{
        append(me.concat, obj);
        appendValue(me.limits, lengthSum(me) + length(obj));
    }

    template < typename TString, typename TDelimiter, typename TString2, typename TExpand >
    inline void appendValue(
		StringSet< TString, Owner<ConcatDirect<TDelimiter> > > &me,
		TString2 const &obj,
		Tag<TExpand> const)
	{
        append(me.concat, obj);
        appendValue(me.concat, TDelimiter());
        appendValue(me.limits, lengthSum(me) + length(obj) + 1);
    }

	// Generous
	template < typename TString, typename TExpand >
	inline void appendValue(
		StringSet<TString, Dependent<Generous> > &me,
		TString const &obj,
		Tag<TExpand> const)
	{
		SEQAN_CHECKPOINT
		appendValue(me.strings, const_cast<TString*>(&obj));
        appendValue(me.limits, lengthSum(me) + length(obj));
	}

	// Tight
	template < typename TString, typename TExpand >
	inline void appendValue(
		StringSet<TString, Dependent<Tight> > &me,
		TString const &obj,
		Tag<TExpand> const)
	{
		SEQAN_CHECKPOINT
		appendValue(me.strings, const_cast<TString*>(&obj));
		appendValue(me.ids, length(me.strings) - 1);
        appendValue(me.limits, lengthSum(me) + length(obj));
	}

/*
    inline void append(TString *_objs[], unsigned count) {
        for(unsigned i = 0; i < count; ++i)
            add(_objs[i]);
    }
*/

///.Function.clear.param.object.type:Class.StringSet
//////////////////////////////////////////////////////////////////////////////
// clear

	template < typename TString >
    inline void	clear(StringSet< TString, Owner<Default> > &me)
	{
	SEQAN_CHECKPOINT
		clear(me.strings);
		resize(me.limits, 1);
		me.limitsValid = true;
    }

    template < typename TString, typename TDelimiter >
    inline void	clear(StringSet< TString, Owner<ConcatDirect<TDelimiter> > > &me)
	{
	SEQAN_CHECKPOINT
		clear(me.concat);
		resize(me.limits, 1);
    }

    template < typename TString >
	inline void	clear(StringSet< TString, Dependent<Generous> > & me)
	{
	SEQAN_CHECKPOINT
		clear(me.strings);
		resize(me.limits, 1);
		me.limitsValid = true;
	}

    template < typename TString >
	inline void	clear(StringSet<TString, Dependent<Tight> >& me)
	{
	SEQAN_CHECKPOINT
		clear(me.strings);
		resize(me.limits, 1);
		me.limitsValid = true;

		clear(me.ids);
	}

///.Function.length.param.object.type:Class.StringSet
//////////////////////////////////////////////////////////////////////////////
// length

    template < typename TString, typename TSpec >
    inline typename Size< StringSet< TString, TSpec > >::Type
	length(StringSet< TString, TSpec > const &me) {
        return length(me.limits) - 1;
    }

	template <typename TString>
	inline typename Size<StringSet<TString, Dependent<Tight> > >::Type
	length(StringSet<TString, Dependent<Tight> > const &me)
	{
		return length(me.strings);
	}

///.Function.resize.param.object.type:Class.StringSet
//////////////////////////////////////////////////////////////////////////////
// resize

	template < typename TString, typename TSpec, typename TSize >
    inline typename Size< StringSet< TString, TSpec > >::Type
	resize(StringSet< TString, TSpec > &me, TSize new_size) {
		resize(me.limits, new_size + 1);
		me.limitsValid = (new_size == 0);
		return resize(me.strings, new_size);
    }

	template < typename TString, typename TSpec, typename TSize >
    inline typename Size< StringSet< TString, Owner<ConcatDirect<TSpec> > > >::Type
	resize(StringSet< TString, Owner<ConcatDirect<TSpec> > > &me, TSize new_size) {
		return resize(me.limits, new_size + 1) - 1;
    }

///.Function.value.param.object.type:Class.StringSet
//////////////////////////////////////////////////////////////////////////////
// value

	// Default
	template < typename TString, typename TPos >
	inline typename Reference< StringSet< TString, Owner<Default> > >::Type
	value(StringSet< TString, Owner<Default> > & me, TPos pos)
	{
		return me.strings[pos];
	}

	template < typename TString, typename TPos >
	inline typename Reference< StringSet< TString, Owner<Default> > const >::Type
	value(StringSet< TString, Owner<Default> > const & me, TPos pos)
	{
		return me.strings[pos];
	}

	// ConcatDirect
	template < typename TString, typename TSpec, typename TPos >
	inline typename Infix<TString>::Type
	value(StringSet< TString, Owner<ConcatDirect<TSpec> > > & me, TPos pos)
	{
		return infix(me.concat, me.limits[pos], me.limits[pos + 1]);
	}

	template < typename TString, typename TSpec, typename TPos >
	inline typename Infix<TString const>::Type
	value(StringSet< TString, Owner<ConcatDirect<TSpec> > > const & me, TPos pos)
	{
		return infix(me.concat, me.limits[pos], me.limits[pos + 1]);
	}

	// Tight
	template < typename TString, typename TPos >
	inline typename Reference<StringSet< TString, Dependent<Tight> > >::Type
	value(StringSet< TString, Dependent<Tight> >& me, TPos pos)
	{
	SEQAN_CHECKPOINT
		if (me.strings[pos])
			return *me.strings[pos];
		static TString tmp = "";
		return tmp;
	}

	template < typename TString, typename TPos >
	inline typename Reference<StringSet< TString, Dependent<Tight> > const >::Type
	value(StringSet< TString, Dependent<Tight> >const & me, TPos pos)
	{
		if (me.strings[pos])
			return *me.strings[pos];
		static TString tmp = "";
		return tmp;
	}

	// Generous
	template < typename TString, typename TPos >
	inline typename Reference<StringSet< TString, Dependent<Generous> > >::Type
	value(StringSet< TString, Dependent<Generous> >& me, TPos pos)
	{
	SEQAN_CHECKPOINT
		unsigned i = _findIthNonZeroValue(me.strings, pos);
		if (i < length(me.strings))
			return *me.strings[i];
		static TString tmp = "";
		return tmp;
	}

	template < typename TString, typename TPos >
	inline typename Reference< StringSet< TString, Dependent<Generous> > const >::Type
	value(StringSet< TString, Dependent<Generous> > const & me, TPos pos)
	{
	SEQAN_CHECKPOINT
		unsigned i = _findIthNonZeroValue(me.strings, pos);
		if (i < length(me.strings))
			return *me.strings[i];
		static TString tmp = "";
		return tmp;
	}

//////////////////////////////////////////////////////////////////////////////
// getValueById

/**
.Function.getValueById:
..cat:Sequences
..summary:Retrieves a string from the StringSet given an id.
..signature:getValueById(me, id)
..param.me:A StringSet.
...type:Class.StringSet
..param.id:An id.
...type:Metafunction.Id
..returns:A reference to a string.
..see:Function.assignValueById
..see:Function.valueById
*/

	template <typename TString, typename TSpec, typename TId>
	inline typename Reference<StringSet<TString, Owner<TSpec> > >::Type
	getValueById(StringSet<TString, Owner<TSpec> >& me,
				TId const id)
	{
	SEQAN_CHECKPOINT
		if (id < (TId) length(me)) return value(me, id);
		static TString tmp = "";
		return tmp;
	}

	template <typename TString, typename TId>
	inline typename Reference<StringSet<TString, Dependent<Generous> > >::Type
	getValueById(StringSet<TString, Dependent<Generous> >& me,
				TId const id)
	{
	SEQAN_CHECKPOINT
		if (me.strings[id])
			return *me.strings[id];
		static TString tmp = "";
		return tmp;
	}

	template <typename TString, typename TId>
	inline typename Reference<StringSet<TString, Dependent<Tight> > >::Type
	getValueById(StringSet<TString, Dependent<Tight> >&me,
				TId const id)
	{
	SEQAN_CHECKPOINT
		for(unsigned i = 0; i < length(me.strings); ++i)
			if ((TId) me.ids[i] == id)
				return value(me, i);
		static TString tmp = "";
		return tmp;
	}


//////////////////////////////////////////////////////////////////////////////
// valueById

/**
.Function.valueById:
..cat:Sequences
..summary:Retrieves a string from the StringSet given an id.
..signature:valueById(me, id)
..param.me:A StringSet.
...type:Class.StringSet
..param.id:An id.
...type:Metafunction.Id
..returns:A reference to a string.
..see:Function.assignValueById
..see:Function.getValueById
*/

	template<typename TString, typename TSpec, typename TId>
	inline typename Reference<StringSet<TString, TSpec> >::Type
	valueById(StringSet<TString, TSpec>& me,
			TId const id)
	{
	SEQAN_CHECKPOINT
		return getValueById(me, id);
	}


//////////////////////////////////////////////////////////////////////////////
// assignValueById

/**
.Function.assignValueById:
..cat:Sequences
..summary:Adds a new string to the StringSet and returns an id.
..signature:assignValueById(dest, str, [id])
..signature:assignValueById(dest, source, id)
..param.dest:A StringSet.
...type:Class.StringSet
..param.source:A StringSet.
...type:Class.StringSet
..param.str:A new string.
...type:Metafunction.Value
..param.id:An associated id.
...type:Metafunction.Id
..returns:A new id
...type:Metafunction.Id
..see:Function.getValueById
..see:Function.valueById
*/

	template<typename TString, typename TSpec, typename TString2>
	inline typename Id<StringSet<TString, TSpec> >::Type
	assignValueById(StringSet<TString, TSpec>& me,
					TString2& obj)
	{
	SEQAN_CHECKPOINT
		appendValue(me, obj);
		SEQAN_ASSERT(length(me.limits) == length(me) + 1);
		return length(me.strings) - 1;
	}

	template <typename TString, typename TSpec, typename TId>
	inline typename Id<StringSet<TString, Owner<TSpec> > >::Type
	assignValueById(StringSet<TString, Owner<TSpec> >& me,
					TString& obj,
					TId id)
	{
	SEQAN_CHECKPOINT
		if (id >= (TId) length(me.strings)) {
			fill(me.strings, id+1, TString());
			resize(me.limits, length(me.limits) + 1);
		}
		assignValue(me, id, obj);
		me.limitsValid = false;
		return id;
	}

	template<typename TString, typename TId>
	inline typename Id<StringSet<TString, Dependent<Generous> > >::Type
	assignValueById(StringSet<TString, Dependent<Generous> >& me,
					TString& obj,
					TId id)
	{
	SEQAN_CHECKPOINT
		SEQAN_ASSERT(length(me.limits) == length(me) + 1);
		if (id >= (TId) length(me.strings)) fill(me.strings, id+1, (TString*) 0);
		if ((TString*) me.strings[id] == (TString*) 0)
			resize(me.limits, length(me.limits) + 1);
		me.strings[id] = &obj;
		me.limitsValid = false;
		SEQAN_ASSERT(length(me.limits) == length(me) + 1);
		return id;
	}

	//////////////////////////////////////////////////////////////////////////////

	template<typename TString, typename TId>
	inline typename Id<StringSet<TString, Dependent<Tight> > >::Type
	assignValueById(StringSet<TString, Dependent<Tight> >& me,
					TString& obj,
					TId id)
	{
	SEQAN_CHECKPOINT
		typedef StringSet<TString, Dependent<Tight> > TStringSet;
		typedef typename Size<TStringSet>::Type TSize;

		for(TSize i = 0; i < length(me.ids); ++i)
			if ((TId) me.ids[i] == id) {
				me.strings[i] = &obj;
				me.limitsValid = false;
				return id;
			}
		appendValue(me.strings, &obj);
		appendValue(me.ids, id);
		return id;
	}

	template<typename TString, typename TSpec1, typename TSpec2, typename TId>
	inline typename Id<StringSet<TString, TSpec1> >::Type
	assignValueById(StringSet<TString, TSpec1>& dest,
					StringSet<TString, TSpec2>& source,
					TId id)
	{
	SEQAN_CHECKPOINT
		return assignValueById(dest, getValueById(source, id), id);
	}

//////////////////////////////////////////////////////////////////////////////
// removeValueById

/**
.Function.removeValueById:
..cat:Sequences
..summary:Removes a string from the StringSet given an id.
..signature:removeValueById(me, id)
..param.me:A StringSet.
...type:Class.StringSet
..param.id:An id.
...type:Metafunction.Id
..returns:void
..see:Function.assignValueById
*/

	template<typename TString, typename TSpec, typename TId>
	inline void
	removeValueById(StringSet<TString, Owner<TSpec> >& me, TId const id)
	{
	SEQAN_CHECKPOINT
		erase(me.strings, id);
		resize(me.limits, length(me.limits) - 1);
		me.limitsValid = empty(me);
	}

	template<typename TString, typename TId>
	inline void
	removeValueById(StringSet<TString, Dependent<Generous> >& me, TId const id)
	{
	SEQAN_CHECKPOINT
		if (me.strings[id] != (TString*) 0) {
			resize(me.limits, length(me.limits) - 1);
			me.limitsValid = empty(me);
		}
		me.strings[id] = 0;
		while (!empty(me.strings) && !me.strings[length(me.strings) - 1])
			resize(me.strings, length(me.strings) - 1);
	}

	template<typename TString, typename TId>
	inline void
	removeValueById(StringSet<TString, Dependent<Tight> >& me, TId const id)
	{
	SEQAN_CHECKPOINT
		typedef StringSet<TString, Dependent<Tight> > TStringSet;
		typedef typename Size<TStringSet>::Type TSize;

		SEQAN_ASSERT(length(me.limits) == length(me) + 1);
		for(TSize i = 0; i < length(me.strings); ++i)
			if (me.ids[i] == id) {
				erase(me.strings, i);
				erase(me.ids, i);
				resize(me.limits, length(me.limits) - 1);
				me.limitsValid = empty(me);
			}
		SEQAN_ASSERT(length(me.limits) == length(me) + 1);
	}

//////////////////////////////////////////////////////////////////////////////
//

/**
.Function.positionToId:
..cat:Sequences
..summary:Retrieves the id of a string in the StringSet given a position.
..signature:positionToId(string_set, pos)
..param.string_set:A StringSet.
...type:Class.StringSet
..param.pos:A position that is transfored into an id.
..returns:An id that corresponds to $pos$ within $string_set$
..see:Function.assignValueById
..see:Function.valueById
*/

	template <typename TString, typename TSpec, typename TPos>
	inline typename Id<StringSet<TString, Owner<TSpec> > >::Type
	positionToId(StringSet<TString, Owner<TSpec> >&,
				TPos const pos)
	{
	SEQAN_CHECKPOINT
		return pos;
	}

	template <typename TString, typename TPos>
	inline typename Id<StringSet<TString, Dependent<Generous> > >::Type
	positionToId(StringSet<TString, Dependent<Generous> >& me,
			TPos const pos)
	{
	SEQAN_CHECKPOINT
		return _findIthNonZeroValue(me.strings,pos);
	}

	template <typename TString, typename TPos>
	inline typename Id<StringSet<TString, Dependent<Generous> > >::Type
	positionToId(StringSet<TString, Dependent<Generous> > const& me,
				TPos const pos)
	{
	SEQAN_CHECKPOINT
		return _findIthNonZeroValue(me.strings,pos);
	}

	template <typename TString, typename TPos>
	inline typename Id<StringSet<TString, Dependent<Tight> > >::Type
	positionToId(StringSet<TString, Dependent<Tight> >&me,
				TPos const pos)
	{
	SEQAN_CHECKPOINT
		return me.ids[pos];
	}

	template <typename TString, typename TPos>
	inline typename Id<StringSet<TString, Dependent<Tight> > >::Type
	positionToId(StringSet<TString, Dependent<Tight> > const&me,
				TPos const pos)
	{
	SEQAN_CHECKPOINT
		return me.ids[pos];
	}


/**
.Function.idToPosition:
..cat:Sequences
..summary:Retrieves the position of a string in the StringSet given an id.
..signature:idToPosition(me, id)
..param.me:A StringSet.
...type:Class.StringSet
..param.id:An id.
...type:Metafunction.Id
..returns:A reference to a string.
..see:Function.assignValueById
..see:Function.valueById
*/

	template <typename TString, typename TSpec, typename TId>
	inline typename Id<StringSet<TString, Owner<TSpec> > >::Type
	idToPosition(StringSet<TString, Owner<TSpec> >&,
				TId const id)
	{
	SEQAN_CHECKPOINT
		return id;
	}

	template <typename TString, typename TId>
	inline typename Id<StringSet<TString, Dependent<Generous> > >::Type
	idToPosition(StringSet<TString, Dependent<Generous> >& me,
				TId const id)
	{
	SEQAN_CHECKPOINT
		return _countNonZeroValues(me.strings,id);
	}

	template <typename TString, typename TId>
	inline typename Id<StringSet<TString, Dependent<Tight> > >::Type
	idToPosition(StringSet<TString, Dependent<Tight> >&me,
				TId const id)
	{
	SEQAN_CHECKPOINT
		for(unsigned i = 0; i < length(me.ids); ++i)
			if ((TId) me.ids[i] == id)
				return i;
		return 0;
	}




//////////////////////////////////////////////////////////////////////////////
// subset

/**
.Function.subset:
..cat:Sequences
..summary:Creates a subset of a given StringSet.
..signature:subset(source, dest, id_array [, len])
..param.source:In-parameter:The source StringSet.
...type:Class.StringSet
..param.dest:Out-parameter:The destination StringSet (the subset).
...type:Class.StringSet
..param.id_array:In-parameter:An array of ids. Each id corresponds to a sequence that is supposed to be in the subset.
..param.len:In-parameter:Optional length of the id array.
...remarks:If len is not defined the length function must be valid for this array or string type.
..returns:void
*/

	template <typename TString, typename TSpec, typename TDestSpec, typename TIds, typename TLength>
	inline void
	subset(StringSet<TString, Owner<TSpec> >& source,
		StringSet<TString, TDestSpec>& dest,
		TIds ids,
		TLength len)
	{
	SEQAN_CHECKPOINT
	}

	template <typename TString, typename TIds, typename TLength>
	inline void
	subset(StringSet<TString, Dependent<Generous> >& source,
		StringSet<TString, Dependent<Generous> >& dest,
		TIds ids,
		TLength len)
	{
	SEQAN_CHECKPOINT
		typedef StringSet<TString, Dependent<Generous> > TStringSet;
		typedef typename Id<TStringSet>::Type TId;
		typedef typename Size<TStringSet>::Type TSize;

		clear(dest);
		resize(dest.limits, len + 1);
		dest.limitsValid = (len == 0);
		fill(dest.strings, length(source.strings), (TString*) 0);
		for(TSize i = 0; i < len; ++i)
			dest.strings[ids[i]] = source.strings[ids[i]];
	}

	template <typename TString, typename TIds, typename TLength>
	inline void
	subset(StringSet<TString, Dependent<Tight> >& source,
		StringSet<TString, Dependent<Tight> >& dest,
		TIds ids,
		TLength len)
	{
	SEQAN_CHECKPOINT
		typedef StringSet<TString, Dependent<Tight> > TStringSet;
		typedef typename Id<TStringSet>::Type TId;
		typedef typename Size<TStringSet>::Type TSize;

		clear(dest);
		resize(dest.limits, len + 1);
		dest.limitsValid = (len == 0);
		TLength upperBound = length(source.ids);
		for(TSize i=0;i<len;++i) {
			TId id = ids[i];
			if ((upperBound > id) &&
				(source.ids[id] == id)) {
					appendValue(dest.strings, source.strings[id]);
					appendValue(dest.ids, id);
			} else {
				typedef String<TId> TIdString;
				typedef typename Iterator<TIdString, Rooted>::Type TIter;
				TIter it = begin(source.ids);
				for(;!atEnd(it);goNext(it)) {
					if (*it == id) {
						appendValue(dest.strings, source.strings[position(it)]);
						appendValue(dest.ids, id);
					}
				}
			}
		}
	}

	template <typename TString, typename TSpec, typename TIds>
	inline void
	subset(StringSet<TString, TSpec>& source,
		StringSet<TString, TSpec>& dest,
		TIds ids)
	{
	SEQAN_CHECKPOINT
		subset(source, dest, ids, length(ids));
	}

//////////////////////////////////////////////////////////////////////////////




	//////////////////////////////////////////////////////////////////////////////
	// ConcatenatorNto1 - a StringSet to String converter
	//////////////////////////////////////////////////////////////////////////////

    template <typename TStringSet>
	struct ConcatenatorNto1 {
		TStringSet *set;
		ConcatenatorNto1 () {}
		ConcatenatorNto1 (TStringSet &_set): set(&_set) {}

//____________________________________________________________________________
// WARNING:
// operator[] conducts a binary search and should be avoided
// you better use StringSet<.., Owner<ConcatDirect<..> > > for random access
// or ConcatenatorNto1's iterators for sequential access

		template <typename TPos>
		inline typename Reference<ConcatenatorNto1>::Type
		operator [] (TPos pos)
		{
	SEQAN_CHECKPOINT
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<ConcatenatorNto1 const>::Type
		operator [] (TPos pos) const
		{
	SEQAN_CHECKPOINT
			return value(*this, pos);
		}
	};
//____________________________________________________________________________

    template <typename TStringSet>
	struct Value< ConcatenatorNto1<TStringSet> > {
		typedef typename Value< typename Value<TStringSet>::Type >::Type Type;
	};

    template <typename TStringSet>
	struct Value< ConcatenatorNto1<TStringSet> const >:
		Value< ConcatenatorNto1<TStringSet> > {};
//____________________________________________________________________________

    template <typename TStringSet>
    struct Size< ConcatenatorNto1<TStringSet> > {
		typedef typename Size< typename Value<TStringSet>::Type >::Type Type;
    };
//____________________________________________________________________________

    template <typename TStringSet>
	struct AllowsFastRandomAccess< ConcatenatorNto1<TStringSet> >
	{
		typedef False Type;
		enum { VALUE = false };
	};

//////////////////////////////////////////////////////////////////////////////
// value

	template < typename TStringSet, typename TPos >
	inline typename Reference< ConcatenatorNto1<TStringSet> >::Type
	value(ConcatenatorNto1<TStringSet> &me, TPos globalPos)
	{
		Pair<unsigned, typename Size< typename Value<TStringSet>::Type >::Type> localPos;
		posLocalize(localPos, globalPos, stringSetLimits(*me.set));
        return value(value(*me.set, getValueI1(localPos)), getValueI2(localPos));
	}

	template < typename TStringSet, typename TPos >
	inline typename Reference< ConcatenatorNto1<TStringSet> const >::Type
	value(ConcatenatorNto1<TStringSet> const &me, TPos globalPos)
	{
		typedef typename Value<TStringSet>::Type TString;
		Pair<unsigned, typename Size<TString>::Type> localPos;
		posLocalize(localPos, globalPos, stringSetLimits(*me.set));
        return value(value(*(TStringSet const*)me.set, getValueI1(localPos)), getValueI2(localPos));
	}

//////////////////////////////////////////////////////////////////////////////
// length

	template < typename TStringSet >
    inline typename Size< ConcatenatorNto1<TStringSet> >::Type
	length(ConcatenatorNto1<TStringSet> const &me) {
        return lengthSum(*me.set);
    }

//////////////////////////////////////////////////////////////////////////////
// begin

	template < typename TStringSet, typename TSpec >
	inline typename Iterator< ConcatenatorNto1<TStringSet const>, Tag<TSpec> const >::Type
	begin(ConcatenatorNto1<TStringSet const> concat, Tag<TSpec> const)
	{
		return typename Iterator< ConcatenatorNto1<TStringSet const>, Tag<TSpec> const >::Type (*concat.set);
	}

	template < typename TStringSet, typename TSpec >
	inline typename Iterator< ConcatenatorNto1<TStringSet>, Tag<TSpec> const >::Type
	begin(ConcatenatorNto1<TStringSet> concat, Tag<TSpec> const)
	{
		return typename Iterator< ConcatenatorNto1<TStringSet>, Tag<TSpec> const >::Type (*concat.set);
	}

//////////////////////////////////////////////////////////////////////////////
// end

	template < typename TStringSet, typename TSpec >
	inline typename Iterator< ConcatenatorNto1<TStringSet const>, Tag<TSpec> const >::Type
	end(ConcatenatorNto1<TStringSet const> concat, Tag<TSpec> const)
	{
		return typename Iterator< ConcatenatorNto1<TStringSet>, Tag<TSpec> const >::Type
			(*concat.set, length(*concat.set), 0);
	}

	template < typename TStringSet, typename TSpec >
	inline typename Iterator< ConcatenatorNto1<TStringSet>, Tag<TSpec> const >::Type
	end(ConcatenatorNto1<TStringSet> concat, Tag<TSpec> const)
	{
		return typename Iterator< ConcatenatorNto1<TStringSet>, Tag<TSpec> const >::Type
			(*concat.set, length(*concat.set), 0);
	}

//////////////////////////////////////////////////////////////////////////////
// Concatenator metafunction

	template < typename TString, typename TSpec >
	struct Concatenator< StringSet<TString, TSpec> > {
		typedef ConcatenatorNto1< StringSet<TString, TSpec> > Type;
	};

	template < typename TString, typename TSpec >
	struct Concatenator< StringSet<TString, Owner<ConcatDirect<TSpec> > > > {
		typedef TString Type;
	};

//////////////////////////////////////////////////////////////////////////////
// concat

	template <typename TString>
	inline typename Concatenator<TString>::Type &
	concat(TString &string) {
		return string;
	}

	template <typename TString, typename TSpec>
	inline typename Concatenator< StringSet<TString, TSpec> >::Type &
	concat(StringSet<TString, TSpec> &set) {
		return set.concat;
	}

	template <typename TString, typename TSpec>
	inline typename Concatenator< StringSet<TString, TSpec> const>::Type &
	concat(StringSet<TString, TSpec> const &set) {
 		return set.concat;
	}


	//////////////////////////////////////////////////////////////////////////////
	// This iterator sequentially iterates through the elements of TStringSet
	// as if they were directly concatenated (compare StringSet<.., Owner<ConcatDirect<> > >
    //////////////////////////////////////////////////////////////////////////////

	template < typename TDelimiter = void >
	struct ConcatVirtual;

    template < typename TStringSet, typename TSpec >
	class Iter< TStringSet, ConcatVirtual<TSpec> >
	{
	public:
		typedef typename Value<TStringSet>::Type		TString;
		typedef typename Value<TString>::Type			TValue;
        typedef typename Size<TString>::Type			TSize;
//____________________________________________________________________________

	public:
        typedef typename Iterator<TString, Standard>::Type			obj_iterator;
        typedef typename Iterator<TString const, Standard>::Type	const_obj_iterator;

		//////////////////////////////////////////////////////////////////////////////
		// STL compatible public iterator interface

        typedef Iter								iterator;
		typedef ::std::bidirectional_iterator_tag	iterator_category;
        typedef TValue								value_type;
		typedef TValue &							reference;
		typedef TValue const &						const_reference;
		typedef TValue*								pointer;
		typedef TSize								size_type;
		typedef typename Difference<TString>::Type	difference_type;
//____________________________________________________________________________

		TStringSet		*host;
        unsigned		objNo;
        obj_iterator	_begin, _cur, _end;
//____________________________________________________________________________

		inline Iter() {}

        inline Iter(TStringSet &_host):
            host(&_host)
        {
            objNo = 0;
            _begin = _cur = begin(_host[objNo]);
            _end = end(_host[objNo]);
            _testEnd();
        }

        inline Iter(TStringSet &_host, unsigned _objNo, difference_type _offset):
            host(&_host)
        {
            if (_objNo < length(_host)) {
	            objNo = _objNo;
				_begin = _cur = begin(_host[objNo]);
				_end = end(_host[objNo]);
                goFurther(_cur, _offset);
                _testEnd();
            } else {
				objNo = length(_host) - 1;
				_begin = _cur = _end = end(_host[objNo]);
            }
        }

        inline operator obj_iterator() {
            return _cur;
        }
//____________________________________________________________________________

		inline bool _atEndOfSequence() {
			if (_cur == _begin && objNo > 0) return true;
			if (_cur == _end) return true;
			return false;
		}

		inline void _testBegin() {
            while (_cur == _begin && objNo > 0) {
                --objNo;
				_begin = host->_begin(objNo);
				_end = _cur = host->_end(objNo);
            }
        }

        inline void _testEnd() {
			while (_cur == _end && objNo < (length(*host) - 1)) {
				++objNo;
				_begin = _cur = begin((*host)[objNo]);
				_end = end((*host)[objNo]);
            };
        }

        inline TSize _tell() const {
			typedef Pair<unsigned, TSize> TPair;
			return posGlobalize(TPair(objNo, difference(_begin, _cur)), stringSetLimits(*host));
        }
    };

	//////////////////////////////////////////////////////////////////////////////
	// ConcatenatorNto1 meta functions
	//////////////////////////////////////////////////////////////////////////////
//____________________________________________________________________________
// default concatenator iterators

    template <typename TString, typename TSpec >
    struct Iterator< ConcatenatorNto1< StringSet<TString, TSpec> >, Standard > {
        typedef Iter<StringSet<TString, TSpec>, ConcatVirtual<> > Type;
    };

    template <typename TString, typename TSpec >
    struct Iterator< ConcatenatorNto1< StringSet<TString, TSpec> const >, Standard > {
        typedef Iter<StringSet<TString, TSpec> const, ConcatVirtual<> > Type;
    };

    template <typename TString, typename TSpec >
    struct Iterator< ConcatenatorNto1< StringSet<TString, TSpec> >, Rooted > {
        typedef Iter<StringSet<TString, TSpec>, ConcatVirtual<> > Type;
    };

    template <typename TString, typename TSpec >
    struct Iterator< ConcatenatorNto1< StringSet<TString, TSpec> const >, Rooted > {
        typedef Iter<StringSet<TString, TSpec> const, ConcatVirtual<> > Type;
    };
//____________________________________________________________________________

	template <typename TStringSet >
    struct Iterator< ConcatenatorNto1<TStringSet> const, Standard > {
		typedef typename Iterator< ConcatenatorNto1<TStringSet>, Standard >::Type Type;
    };

    template <typename TStringSet >
    struct Iterator< ConcatenatorNto1<TStringSet> const, Rooted > {
		typedef typename Iterator< ConcatenatorNto1<TStringSet>, Rooted >::Type Type;
    };

	//////////////////////////////////////////////////////////////////////////////
	// meta functions
	//////////////////////////////////////////////////////////////////////////////

	template <typename TStringSet, typename TSpec>
	struct Value< Iter< TStringSet, ConcatVirtual<TSpec> > >:
		Value< typename Value<TStringSet>::Type > {};

	template <typename TStringSet, typename TSpec>
	struct Value< Iter< TStringSet, ConcatVirtual<TSpec> > const >:
		Value< typename Value<TStringSet>::Type > {};

	template <typename TStringSet, typename TSpec>
	struct GetValue< Iter< TStringSet, ConcatVirtual<TSpec> > >:
		GetValue< typename Value<TStringSet>::Type > {};

	template <typename TStringSet, typename TSpec>
	struct GetValue< Iter< TStringSet, ConcatVirtual<TSpec> > const >:
		GetValue< typename Value<TStringSet>::Type > {};

	template <typename TStringSet, typename TSpec>
	struct Size< Iter< TStringSet, ConcatVirtual<TSpec> > >:
		Size< typename Value<TStringSet>::Type > {};

	template <typename TStringSet, typename TSpec>
	struct Reference< Iter< TStringSet, ConcatVirtual<TSpec> > >:
		Reference< typename Value<TStringSet>::Type > {};

	template <typename TStringSet, typename TSpec>
	struct Reference< Iter< TStringSet, ConcatVirtual<TSpec> > const >:
		Reference< typename Value<TStringSet>::Type > {};


	//////////////////////////////////////////////////////////////////////////////
	// operator *
	//////////////////////////////////////////////////////////////////////////////

	template <typename TStringSet, typename TSpec>
	inline typename Reference< Iter< TStringSet, ConcatVirtual<TSpec> > const>::Type
	value(Iter<TStringSet, ConcatVirtual<TSpec> > const & me) {
        return *me._cur;
    }

	template <typename TStringSet, typename TSpec>
	inline typename Reference< Iter< TStringSet, ConcatVirtual<TSpec> > >::Type
	value(Iter<TStringSet, ConcatVirtual<TSpec> > & me) {
        return *me._cur;
    }

	template <typename TStringSet, typename TSpec>
	inline typename Reference< Iter< TStringSet, ConcatVirtual<TSpec> > const>::Type
	operator * (Iter<TStringSet, ConcatVirtual<TSpec> > const & me) {
        return *me._cur;
    }

	template <typename TStringSet, typename TSpec>
	inline typename Reference< Iter< TStringSet, ConcatVirtual<TSpec> > >::Type
	operator * (Iter<TStringSet, ConcatVirtual<TSpec> > & me) {
        return *me._cur;
    }

	//////////////////////////////////////////////////////////////////////////////
	// operator ++
	//////////////////////////////////////////////////////////////////////////////

	template <typename TStringSet, typename TSpec>
	inline void
	goNext(Iter<TStringSet, ConcatVirtual<TSpec> > & me) {
        ++me._cur;
        me._testEnd();
    }

	template <typename TStringSet, typename TSpec>
	inline Iter<TStringSet, ConcatVirtual<TSpec> > const &
	operator ++ (Iter<TStringSet, ConcatVirtual<TSpec> > & me) {
        goNext(me);
		return me;
    }

	template <typename TStringSet, typename TSpec>
	inline Iter<TStringSet, ConcatVirtual<TSpec> > const &
	operator ++ (Iter<TStringSet, ConcatVirtual<TSpec> > & me, int) {
		Iter<TStringSet, ConcatVirtual<TSpec> > before = me;
        goNext(me);
		return before;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator --
	//////////////////////////////////////////////////////////////////////////////

	template <typename TStringSet, typename TSpec>
	inline void
	goPrevious(Iter<TStringSet, ConcatVirtual<TSpec> > & me) {
        me._testBegin();
        --me._cur;
    }

	template <typename TStringSet, typename TSpec>
	inline Iter<TStringSet, ConcatVirtual<TSpec> > const &
	operator -- (Iter<TStringSet, ConcatVirtual<TSpec> > & me) {
        goPrevious(me);
		return me;
    }

	template <typename TStringSet, typename TSpec>
	inline Iter<TStringSet, ConcatVirtual<TSpec> > const &
	operator -- (Iter<TStringSet, ConcatVirtual<TSpec> > & me, int) {
		Iter<TStringSet, ConcatVirtual<TSpec> > before = me;
        goPrevious(me);
		return before;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator +
	//////////////////////////////////////////////////////////////////////////////

	template <typename TStringSet, typename TSpec, typename TDelta>
	inline Iter<TStringSet, ConcatVirtual<TSpec> >
	operator + (Iter<TStringSet, ConcatVirtual<TSpec> > const & me, TDelta delta) {
		Pair<unsigned, typename Size< typename Value<TStringSet>::Type >::Type> pos;
		posLocalize(pos, me._tell() + delta, stringSetLimits(*me.host));
        return Iter<TStringSet, ConcatVirtual<TSpec> > (*me.host, getValueI1(pos), getValueI2(pos));
    }

	template <typename TStringSet, typename TSpec, typename T1, typename T2, typename TCompression>
	inline Iter<TStringSet, ConcatVirtual<TSpec> >
	operator + (Iter<TStringSet, ConcatVirtual<TSpec> > const & me, Pair<T1, T2, TCompression> delta) {
		Pair<unsigned, typename Size< typename Value<TStringSet>::Type >::Type> pos;
		posLocalize(pos, me._tell() + delta, stringSetLimits(*me.host));
        return Iter<TStringSet, ConcatVirtual<TSpec> > (*me.host, getValueI1(pos), getValueI2(pos));
    }

	//////////////////////////////////////////////////////////////////////////////
	// operator -
	//////////////////////////////////////////////////////////////////////////////

	template <typename TSSetL, typename TSpecL, typename TSSetR, typename TSpecR>
	typename Difference<Iter<TSSetL, ConcatVirtual<TSpecL> > >::Type
	operator - (
		Iter<TSSetL, ConcatVirtual<TSpecL> > const &L,
		Iter<TSSetR, ConcatVirtual<TSpecR> > const &R)
	{
        return L._tell() - R._tell();
    }

	template <typename TStringSet, typename TSpec, typename TDelta>
	inline Iter<TStringSet, ConcatVirtual<TSpec> >
	operator - (Iter<TStringSet, ConcatVirtual<TSpec> > const & me, TDelta delta) {
		Pair<unsigned, typename Size< typename Value<TStringSet>::Type >::Type> pos;
		posLocalize(pos, me._tell() - delta, stringSetLimits(*me.host));
        return Iter<TStringSet, ConcatVirtual<TSpec> > (*me.host, getValueI1(pos), getValueI2(pos));
    }

	//////////////////////////////////////////////////////////////////////////////
	// operator ==
	//////////////////////////////////////////////////////////////////////////////

	template <typename TSSetL, typename TSpecL, typename TSSetR, typename TSpecR>
	inline bool
	operator == (
		Iter<TSSetL, ConcatVirtual<TSpecL> > const &L,
		Iter<TSSetR, ConcatVirtual<TSpecR> > const &R)
	{
		SEQAN_ASSERT(L.host == R.host);
		return L.objNo == R.objNo && L._cur == R._cur;
	}

	template <typename TSSetL, typename TSpecL, typename TSSetR, typename TSpecR>
	inline bool
	operator != (
		Iter<TSSetL, ConcatVirtual<TSpecL> > const &L,
		Iter<TSSetR, ConcatVirtual<TSpecR> > const &R)
	{
		SEQAN_ASSERT(L.host == R.host);
		return L.objNo != R.objNo || L._cur != R._cur;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator <
	//////////////////////////////////////////////////////////////////////////////

	template <typename TSSetL, typename TSpecL, typename TSSetR, typename TSpecR>
	inline bool
	operator < (
		Iter<TSSetL, ConcatVirtual<TSpecL> > const &L,
		Iter<TSSetR, ConcatVirtual<TSpecR> > const &R)
	{
		SEQAN_ASSERT(L.host == R.host);
		return L.objNo < R.objNo || (L.objNo == R.objNo && L._cur < R._cur);
	}

	template <typename TSSetL, typename TSpecL, typename TSSetR, typename TSpecR>
	inline bool
	operator > (
		Iter<TSSetL, ConcatVirtual<TSpecL> > const &L,
		Iter<TSSetR, ConcatVirtual<TSpecR> > const &R)
	{
		SEQAN_ASSERT(L.host == R.host);
		return L.objNo > R.objNo || (L.objNo == R.objNo && L._cur > R._cur);
	}

	//////////////////////////////////////////////////////////////////////////////
	// container
	//////////////////////////////////////////////////////////////////////////////

	template <typename TSSet, typename TSpec>
	inline typename Concatenator<TSSet>::Type
	container(Iter<TSSet, ConcatVirtual<TSpec> > &me)
	{
		return concat(*me.host);
	}

	template <typename TSSet, typename TSpec>
	inline typename Concatenator<TSSet>::Type
	container(Iter<TSSet, ConcatVirtual<TSpec> > const &me)
	{
		return concat(*me.host);
	}

	//////////////////////////////////////////////////////////////////////////////
	// atBegin
	//////////////////////////////////////////////////////////////////////////////

	template <typename TSSet, typename TSpec>
	inline bool
	atBegin(Iter<TSSet, ConcatVirtual<TSpec> > &me)
	{
		return me._cur == me._begin && me.objNo == 0;
	}

	template <typename TSSet, typename TSpec>
	inline bool
	atBegin(Iter<TSSet, ConcatVirtual<TSpec> > const &me)
	{
		return me._cur == me._begin && me.objNo == 0;
	}

	//////////////////////////////////////////////////////////////////////////////
	// atEnd
	//////////////////////////////////////////////////////////////////////////////

	template <typename TSSet, typename TSpec>
	inline bool
	atEnd(Iter<TSSet, ConcatVirtual<TSpec> > &me)
	{
		return me._cur == me._end && me.objNo == (length(*me.host) - 1);
	}

	template <typename TSSet, typename TSpec>
	inline bool
	atEnd(Iter<TSSet, ConcatVirtual<TSpec> > const &me)
	{
		return me._cur == me._end && me.objNo == (length(*me.host) - 1);
	}

	//////////////////////////////////////////////////////////////////////////////
	// atEndOfSequence
	//////////////////////////////////////////////////////////////////////////////

	template <typename TIterator>
	inline bool
	atEndOfSequence(TIterator const &me)
	{
		return atEnd(me);
	}

	template <typename TSSet, typename TSpec>
	inline bool
	atEndOfSequence(Iter<TSSet, ConcatVirtual<TSpec> > const &me)
	{
		return me._atEndOfSequence();
	}

	template <typename TIterator>
	inline bool
	atEndOfSequence(TIterator &me)
	{
		return atEndOfSequence(reinterpret_cast<TIterator const &>(me));
	}

}

#endif
