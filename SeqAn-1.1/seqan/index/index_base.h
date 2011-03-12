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
  $Id: index_base.h,v 1.3 2009/03/13 14:34:32 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_BASE_H
#define SEQAN_HEADER_INDEX_BASE_H

//#define SEQAN_TEST_INDEX

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// needful forward declarations

	// suffix array construction specs
	struct Skew3;
	struct Skew7;
	struct LarssonSadakane;
	struct ManberMyers;
	struct SAQSort;
	struct QGram_Alg;

	// lcp table construction algorithms
	struct Kasai;
	struct KasaiOriginal;	// original, but more space-consuming algorithm

	// enhanced suffix array construction algorithms
	struct ChildTab;
	struct BWT;

	template <typename TSpec = void>
	struct Index_ESA;


//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.DefaultIndexSpec:
..cat:Index
..summary:Default @Class.Index@ specialization type.
..signature:DefaultIndexSpec<TText>::Type
..param.TText:The given text type.
..returns:Can be @Spec.Index_ESA@ or $Index_QGram$, etc.
..remarks:Currently @Spec.Index_ESA@ is default if $TText$ is a @Class.String@.
*/
    template < typename TObject >
    struct DefaultIndexSpec {
        typedef Index_ESA<> Type;
    };

//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.DefaultIndexStringSpec:
..cat:Index
..summary:Default @Class.String@ specialization type of the @Metafunction.Fibre@ of an @Class.Index@.
..signature:DefaultIndexStringSpec<TIndex>::Type
..param.TIndex:An @Class.Index@ Type.
..returns:If the underlying text is a @Class.String@ or a set of Strings (see @Class.StringSet@) the String's spec. type is returned.
..remarks:Most of the @Class.Index@ fibres are strings. The @Class.String@ specialization type is chosen by this meta-function.
*/
    template < typename TIndex >
    struct DefaultIndexStringSpec {
        typedef Alloc<> Type;
    };

//    template < typename TValue, typename TSpec >
//    struct DefaultIndexStringSpec< String<TValue, External<TSpec> > > {
//        typedef External<TSpec> Type;
//    };

	template < typename TString, typename TSpec >
	struct DefaultIndexStringSpec< StringSet<TString, TSpec> >:
		DefaultIndexStringSpec<TString> {};


//////////////////////////////////////////////////////////////////////////////
/**
.Class.Index:
..summary:Contains preprocessing data of a fixed text. Allows fast dictionary look-up and advanced computations.
..cat:Index
..signature:Index<TText[, TSpec]>
..param.TText:The text type.
...type:Class.String
...metafunction:Metafunction.Host
..param.TSpec:The index type.
...default:The result of @Metafunction.DefaultIndexSpec@
...metafunction:Metafunction.Spec
..remarks:An index contains various arrays or objects, also called fibres (see @Metafunction.Fibre@).
..remarks:These fibres are created on demand depending on the requirements of an algorithm.
*/

///.Function.setHaystack.param.haystack.type:Class.Index

	// index as a haystack
	template <
        typename TObject,
        typename TSpec = typename DefaultIndexSpec<TObject>::Type >
	class Index;

	template <typename TObject, typename TSpec>
	struct Host< Index<TObject, TSpec> > {
		typedef TObject Type;
	};

	template <typename TObject, typename TSpec>
	struct Spec< Index<TObject, TSpec> > {
		typedef TSpec Type;
	};


//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.Fibre:
..summary:Type of a specific bundle member (fibre).
..signature:Fibre<TIndex, TSpec>::Type
..cat:Index
..param.TIndex:The fibre container type.
..param.TSpec:Type to specify the fibre.
..returns:Fibre type.
..remarks:An @Class.Index@ can be seen as a bundle consisting of various fibres. In most cases this type is $String<Size<TIndex>::Type>$.
..remarks:A @Metafunction.Fibre@ need not to be a real container. It can also be view (see @Tag.ESA Index Fibres.ESA_RawText@).
*/
	// meta function to get the type of a bundle fibre
	template < typename TIndex, typename TSpec >
	struct Fibre {
		typedef String< typename Size<TIndex>::Type > Type;
	};

	template < typename TIndex, typename TSpec >
	struct Fibre<TIndex const, TSpec> {
		typedef typename Fibre<TIndex, TSpec>::Type const Type;
	};

	struct FibreRecord {
		unsigned	id;
		void*		ptr;
		bool		owner;
	};

	// less function to search in sorted list for fibre id
	struct FibreLess: public ::std::binary_function<FibreRecord, unsigned, bool>
	{	// functor for operator>
		inline bool operator()(FibreRecord const & _Left, unsigned const _Right) const
		{	// apply operator> to operands
			return (_Left.id < _Right);
		}
	};

//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.DefaultIndexCreator:
..cat:Index
..summary:Default algorithm to create a demanded and not yet existing @Metafunction.Fibre@.
..signature:DefaultIndexCreator<TIndex, TFibre>::Type
..param.TIndex:An @Class.Index@ Type.
..param.TFibre:A tag specifying the fibre (e.g. @Tag.ESA Index Fibres.ESA_SA@).
..returns:A tag specifying the default algorithm to create the fibre with.
*/
    // standard algorithm for indices creation
    template < typename TIndex, typename TFibre >
	struct DefaultIndexCreator {
		typedef Default Type;
	};

//////////////////////////////////////////////////////////////////////////////
/**
	.Class.Bundle:
	..summary:General purpose container of various members.
	..signature:Bundle<TValue, TSize>
	..param.TValue:The value type, that is the type of the items/characters stored in the string.
	...remarks:Use @Metafunction.Value@ to get the value type for a given class.
	..param.TSpec:The specializing type.
	...default:$Alloc<>$, see @Spec.Alloc String@.
*/
/*
	template < typename TSpec = void >
	struct Bundle {
		typedef ::std::vector<FibreRecord>	TFibreRecords;
		TFibreRecords						fibres;
	};

	template < typename TBundleSpec, typename TFibreSpec >
	inline FibreRecord& getRecord(Bundle<TBundleSpec> &bundle, TFibreSpec const) {
		unsigned id = (unsigned)_ClassIdentifier<TFibreSpec>::getID();

		typename Bundle<TBundleSpec>::TFibreRecords::iterator first = lower_bound(bundle.fibres.begin(), bundle.fibres.end(), id, FibreLess());
		if (!first->id != id) {
			FibreRecord rec;
			rec.id = id;
			rec.ptr = NULL;
			rec.owner = true;
			bundle.fibres.insert(first, rec);
		} else
			return *first;
	}

	template < typename TBundleSpec, typename TFibreSpec >
	inline typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type & getFibre(Bundle<TBundleSpec> &bundle, TFibreSpec const) {
		typedef typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type Type;
		unsigned id = (unsigned)_ClassIdentifier<TFibreSpec>::getID();

		FibreRecord &rec = getRecord(bundle, TFibreSpec());
		if (!rec.ptr)
			rec.ptr = new Type();
		return *reinterpret_cast<Type*>(rec.ptr);
	}

	template < typename TBundleSpec, typename TFibreSpec >
	inline typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type const & getFibre(Bundle<TBundleSpec> const &bundle, TFibreSpec const) {
		typedef typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type Type;
		unsigned id = (unsigned)_ClassIdentifier<TFibreSpec>::getID();

		FibreRecord &rec = getRecord(bundle, TFibreSpec());
		return *reinterpret_cast<Type*>(rec.ptr);
	}
*/

//////////////////////////////////////////////////////////////////////////////
// various fibre specs for enhanced suffix arrays

	struct _Fibre_Text;		// Original text. Can be a String or a StringSet
	struct _Fibre_RawText;	// Concatenation of the strings above
	struct _Fibre_SA;		// suffix array (of raw text with virtual $-delimiters) with Pair entries
	struct _Fibre_RawSA;	// suffix array with integer entries
	struct _Fibre_SAE;		// suffix array reordered in a b-tree
	struct _Fibre_LCP;		// lcp table of raw text
	struct _Fibre_LCPE;		// lcp interval tree
	struct _Fibre_ChildTab;	// childtab (Kurtz et al.) of raw text
	struct _Fibre_BWT;		// burrows wheeler table of raw text

	typedef Tag<_Fibre_Text> const		Fibre_Text;
	typedef Tag<_Fibre_RawText> const	Fibre_RawText;
	typedef Tag<_Fibre_SA> const		Fibre_SA;
	typedef Tag<_Fibre_RawSA> const		Fibre_RawSA;
	typedef Tag<_Fibre_SAE> const		Fibre_SAE;
	typedef Tag<_Fibre_LCP> const		Fibre_LCP;
	typedef Tag<_Fibre_LCPE> const		Fibre_LCPE;
	typedef Tag<_Fibre_ChildTab> const	Fibre_ChildTab;
	typedef Tag<_Fibre_BWT> const		Fibre_BWT;

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.SAValue:
..cat:Index
..summary:The default alphabet type of a suffix array, i.e. the type to store a position of a string or string set.
..signature:SAValue<TObject>::Type
..param.TObject:A string, string set, or index type.
...type:Class.String
...type:Class.StringSet
...type:Class.Index
..returns:A type to store a position.
...text:If $TObject$ is a @Class.String@, it is a single integer value. By default this is the @Metafunction.Size@ type of $TObject$.
...text:If $TObject$ is a @Class.StringSet@, it could be a single integer too (called global position, see @Spec.ConcatDirect@) or a @Class.Pair@ (called local position, see @Spec.Owner@).
Currently SeqAn defaults to a local position for @Class.StringSet@ classes (index_base.h):
...code:template < typename TString, typename TSpec >
struct SAValue< StringSet<TString, TSpec> > {
	typedef Pair<
		typename Size< StringSet<TString, TSpec> >::Type,
		typename SAValue<TString>::Type,
		Compressed
	> Type;
};
..note:SAValue is the return type of various function, e.g. @Function.position@ for the @Class.Index@ @Class.Finder@ class, @Function.getOccurrence@, @Function.getOccurrences@ etc.
You should always use the type of this meta-function to store the return values.
If you want to write algorithms for both variants (local and global positions) you
should use the functions @Function.posLocalize@, @Function.posGlobalize@, @Function.getSeqNo@ and @Function.getSeqOffset@.
..note:If $TObject$ is an @Class.Index@, @Metafunction.Position@ returns the same value as $SAValue$. You can change the position type of an index by overloading $SAValue$, not @Metafunction.Position@.
*/
	template <typename TObject>
	struct SAValue:
		Size<TObject> {};

	template <typename TObject>
	struct SAValue<TObject const>:
		SAValue<TObject> {};

	// to speed up sequence number computation
	// we use a pair of seqNo and localPosition
	template < typename TString, typename TSpec >
	struct SAValue< StringSet<TString, TSpec> > {
		typedef Pair<
			typename Size< StringSet<TString, TSpec> >::Type,
			typename SAValue<TString>::Type,
			Compressed
		> Type;
	};

/*
	template < typename TString, typename TSpec >
	struct SAValue< StringSet<TString, TSpec> > {
		typedef Pair<
			typename Size< StringSet<TString, TSpec> >::Type,
			typename SAValue<TString>::Type,
			CutCompressed<4>						// max. 4 sequences
		> Type;										// max. 2^30 characters each
	};
*/
	template < typename TText, typename TSpec >
	struct SAValue< Index<TText, TSpec> >:
		SAValue<TText> {};

	template < typename TObject, typename TSpec >
	struct DefaultIndexStringSpec< Index<TObject, TSpec> >:
		DefaultIndexStringSpec<TObject> {};

//////////////////////////////////////////////////////////////////////////////
// value and size type of an index

	template < typename TText, typename TSpec >
    struct Value< Index<TText, TSpec> > {
		typedef typename Value<
			typename Fibre< Index<TText, TSpec>, Fibre_RawText>::Type
		>::Type Type;
    };

	template < typename TText, typename TSpec >
    struct Size< Index<TText, TSpec> > {
		typedef typename Size<
			typename Fibre< Index<TText, TSpec>, Fibre_RawText>::Type
		>::Type Type;
    };

	template < typename TText, typename TSpec >
	struct Position< Index<TText, TSpec> >:
		SAValue< Index<TText, TSpec> > {};

//////////////////////////////////////////////////////////////////////////////
// default table type

	template < typename TObject, typename TSpec, typename TFibre >
	struct Fibre< Index<TObject, TSpec>, Tag<TFibre> const > {
		typedef String<
			typename Size< Index<TObject, TSpec> >::Type,
			typename DefaultIndexStringSpec< Index<TObject, TSpec> >::Type
		> Type;
	};

//////////////////////////////////////////////////////////////////////////////
// original text

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, Fibre_Text> {
		typedef TText Type;
	};

//////////////////////////////////////////////////////////////////////////////
// concatenated text

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, Fibre_RawText> {
		typedef typename Concatenator<TText>::Type Type;
	};

//////////////////////////////////////////////////////////////////////////////
// suffix array type

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, Fibre_SA> {
		typedef String<
			typename SAValue< Index<TText, TSpec> >::Type,
			typename DefaultIndexStringSpec< Index<TText, TSpec> >::Type
		> Type;
	};

//////////////////////////////////////////////////////////////////////////////
// globalize functor

	template <typename InType, typename TLimitsString, typename Result = typename Value<TLimitsString>::Type>
	struct FunctorGlobalize : public ::std::unary_function<InType,Result> {
		TLimitsString const *limits;

		FunctorGlobalize() {}
		FunctorGlobalize(TLimitsString const &_limits) : limits(&_limits) {}
        inline Result operator()(const InType& x) const
        {
			return posGlobalize(x, *limits);
		}
    };

//////////////////////////////////////////////////////////////////////////////
// raw suffix array contains integer offsets relative to raw text

	template < typename TString, typename TSSetSpec, typename TSpec >
	struct Fibre< Index<StringSet<TString, TSSetSpec>, TSpec>, Fibre_RawSA>
	{
		typedef Index< StringSet<TString, TSSetSpec>, TSpec> TIndex;
//		typedef ModifiedString<
//			typename Fibre<TIndex, Fibre_SA>::Type,
//			ModView< FunctorGlobalize<
//				typename Value< typename Fibre<TIndex, Fibre_SA>::Type >::Type,
//				typename StringSetLimits<StringSet<TString, TSSetSpec> >::Type >
//			>
//		> Type;
	};

//////////////////////////////////////////////////////////////////////////////
// default burrows-wheeler table

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, Fibre_BWT> {
		typedef String <
			typename Value< Index<TText, TSpec> >::Type,
			typename DefaultIndexStringSpec< Index<TText, TSpec> >::Type
		> Type;
	};


//////////////////////////////////////////////////////////////////////////////
// default fibre creators

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, Fibre_SA> {
        typedef Skew7 Type;							// standard suffix array creator is skew7
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, Fibre_LCP> {
        typedef Kasai Type;
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, Fibre_BWT> {
        typedef BWT Type;
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, Fibre_ChildTab> {
        typedef ChildTab Type;
    };


//////////////////////////////////////////////////////////////////////////////
// fibre interface to access the enhanced suffix array tables

/**
.Function.getFibre:
..summary:Returns a specific @Metafunction.Fibre@ of an @Class.Index@ object.
..cat:Index
..signature:getFibre(index, fibre_tag)
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..param.fibre_tag:A tag that identifies the @Metafunction.Fibre@ (e.g. @Tag.ESA Index Fibres.ESA_SA@).
..returns:A reference to the @Metafunction.Fibre@ object.
*/

	template <typename TText, typename TSpec>
	inline Holder<TText> & _dataHost(Index<TText, TSpec> &index) {
		return index.text;
	}
	template <typename TText, typename TSpec>
	inline Holder<TText> const & _dataHost(Index<TText, TSpec> const &index) {
		return index.text;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_Text>::Type &
	getFibre(Index<TText, TSpec> &index, Fibre_Text) {
		return value(index.text);
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_Text>::Type &
	getFibre(Index<TText, TSpec> const &index, Fibre_Text) {
		return value(index.text);
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_RawText>::Type &
	getFibre(Index<TText, TSpec> &index, Fibre_RawText) {
		return concat(value(index.text));
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_RawText>::Type &
	getFibre(Index<TText, TSpec> const &index, Fibre_RawText) {
		return concat(value(index.text));
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_SA>::Type &
	getFibre(Index<TText, TSpec> &index, Fibre_SA) {
		return index.sa;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_SA>::Type &
	getFibre(Index<TText, TSpec> const &index, Fibre_SA) {
		return index.sa;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_SA>::Type &
	getFibre(Index<TText, TSpec> &index, Fibre_RawSA) {
		return indexSA(index);
	}
/*
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_SA>::Type &
	getFibre(Index<TText, TSpec> const &index, Fibre_RawSA) {
		return indexSA(index);
	}
*/
	template <typename TString, typename TSSetSpec, typename TSpec>
	inline typename Fibre<Index<StringSet<TString, TSSetSpec>, TSpec>, Fibre_RawSA>::Type
	getFibre(Index<StringSet<TString, TSSetSpec>, TSpec> &index, Fibre_RawSA)
	{
		typedef Index< StringSet<TString, TSSetSpec>, TSpec> TIndex;

		typedef FunctorGlobalize<
			typename Value< typename Fibre<TIndex, Fibre_SA>::Type >::Type,
			typename StringSetLimits<StringSet<TString, TSSetSpec> >::Type
		> TFunctor;

//		typedef ModifiedString<
//			typename Fibre<Index<StringSet<TString, TSSetSpec>, TSpec>, Fibre_SA>::Type,
//			ModView< TFunctor >
//		> ModString;

		return ModString(indexSA(index), TFunctor(stringSetLimits(indexText(index))));
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_LCP>::Type &
	getFibre(Index<TText, TSpec> &index, Fibre_LCP) {
		return index.lcp;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_LCP>::Type &
	getFibre(Index<TText, TSpec> const &index, Fibre_LCP) {
		return index.lcp;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_LCPE>::Type &
	getFibre(Index<TText, TSpec> &index, Fibre_LCPE) {
		return index.lcpe;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_LCPE>::Type &
	getFibre(Index<TText, TSpec> const &index, Fibre_LCPE) {
		return index.lcpe;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_ChildTab>::Type &
	getFibre(Index<TText, TSpec> &index, Fibre_ChildTab) {
		return index.childtab;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_ChildTab>::Type &
	getFibre(Index<TText, TSpec> const &index, Fibre_ChildTab) {
		return index.childtab;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_BWT>::Type &
	getFibre(Index<TText, TSpec> &index, Fibre_BWT) {
		return index.bwt;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_BWT>::Type &
	getFibre(Index<TText, TSpec> const &index, Fibre_BWT) {
		return index.bwt;
	}

//////////////////////////////////////////////////////////////////////////////
///.Function.length.param.object.type:Class.Index

	template <typename TText, typename TSpec>
	inline typename Size<Index<TText, TSpec> >::Type
	length(Index<TText, TSpec> const &index) {
		return length(indexRawText(index));
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Size<TText>::Type
	countSequences(Index<TText, TSpec> const &index) {
		return countSequences(indexText(index));
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TSeqNo, typename TText, typename TSpec>
	inline typename Size<Index<TText, TSpec> >::Type
	sequenceLength(TSeqNo seqNo, Index<TText, TSpec> const &index) {
		return sequenceLength(seqNo, indexText(index));
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TPos, typename TText, typename TSpec>
	inline typename Size<Index<TText, TSpec> >::Type
	suffixLength(TPos pos, Index<TText, TSpec> const &index) {
		return sequenceLength(getSeqNo(pos, stringSetLimits(index)), index) - getSeqOffset(pos, stringSetLimits(index));
	}


//////////////////////////////////////////////////////////////////////////////
// unified textAt interface

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, Fibre_RawText>::Type>::Type
	textAt(TPos i, TIndex &index) {
		return value(getFibre(index, Fibre_RawText()), i);
	}
	template <typename TPos, typename TString, typename TSSetSpec, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, TSSetSpec>, TSpec>, Fibre_RawText>::Type>::Type
	textAt(TPos i, Index< StringSet<TString, TSSetSpec>, TSpec> &index) {
		return value(getFibre(index, Fibre_RawText()), posGlobalize(i, stringSetLimits(index)));
	}
	template <typename TPos, typename TString, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, Owner<Default> >, TSpec>, Fibre_RawText>::Type>::Type
	textAt(TPos i, Index< StringSet<TString, Owner<Default> >, TSpec> &index) {
		Pair <
			typename Size< StringSet<TString, Owner<Default> > >::Type,
			typename Size< TString >::Type > locPos;
		posLocalize(locPos, i, stringSetLimits(index));
		return value(value(getFibre(index, Fibre_Text()), getValueI1(locPos)), getValueI2(locPos));
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.rawtextAt:
..summary:Shortcut for $value(indexRawText(..), ..)$.
..cat:Index
..signature:rawtextAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference or proxy to the value.
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, Fibre_RawText>::Type>::Type rawtextAt(TPos i, TIndex &index) {
		return value(getFibre(index, Fibre_RawText()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, Fibre_RawText>::Type>::Type rawtextAt(TPos i, TIndex const &index) {
		return value(getFibre(index, Fibre_RawText()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.saAt:
..summary:Shortcut for $value(indexSA(..), ..)$.
..cat:Index
..signature:saAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference or proxy to the value.
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, Fibre_SA>::Type>::Type saAt(TPos i, TIndex &index) {
		return value(getFibre(index, Fibre_SA()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, Fibre_SA>::Type>::Type saAt(TPos i, TIndex const &index) {
		return value(getFibre(index, Fibre_SA()), i);
	}

	template <typename TPos, typename TIndex>
	inline typename Value<typename Fibre<TIndex const, Fibre_RawSA>::Type>::Type rawsaAt(TPos i, TIndex const &index) {
		return posGlobalize(saAt(i, index), stringSetLimits(indexText(index)));
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.lcpAt:
..summary:Shortcut for $value(indexLCP(..), ..)$.
..cat:Index
..signature:lcpAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference or proxy to the value.
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, Fibre_LCP>::Type>::Type lcpAt(TPos i, TIndex &index) {
		return value(getFibre(index, Fibre_LCP()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, Fibre_LCP>::Type>::Type lcpAt(TPos i, TIndex const &index) {
		return value(getFibre(index, Fibre_LCP()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.lcpeAt:
..summary:Shortcut for $value(indexLCPE(..), ..)$.
..cat:Index
..signature:lcpeAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference or proxy to the value.
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, Fibre_LCPE>::Type>::Type lcpeAt(TPos i, TIndex &index) {
		return value(getFibre(index, Fibre_LCPE()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, Fibre_LCPE>::Type>::Type lcpeAt(TPos i, TIndex const &index) {
		return value(getFibre(index, Fibre_LCPE()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.childAt:
..summary:Shortcut for $value(indexChildTab(..), ..)$.
..cat:Index
..signature:childAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference or proxy to the value.
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, Fibre_ChildTab>::Type>::Type childAt(TPos i, TIndex &index) {
		return value(getFibre(index, Fibre_ChildTab()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, Fibre_ChildTab>::Type>::Type childAt(TPos i, TIndex const &index) {
		return value(getFibre(index, Fibre_ChildTab()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.bwtAt:
..summary:Shortcut for $value(indexBWT(..), ..)$.
..cat:Index
..signature:bwtAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference or proxy to the value.
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, Fibre_BWT>::Type>::Type bwtAt(TPos i, TIndex &index) {
		return value(getFibre(index, Fibre_BWT()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, Fibre_BWT>::Type>::Type bwtAt(TPos i, TIndex const &index) {
		return value(getFibre(index, Fibre_BWT()), i);
	}

//////////////////////////////////////////////////////////////////////////////
// interface for infinity/invalid values

	template <typename TValue>
	inline void _setSizeInval(TValue &v) {
		v = SupremumValue<TValue>::VALUE;
	}

	template <typename TValue>
	inline bool _isSizeInval(TValue const &v) {
		return v == SupremumValue<TValue>::VALUE;
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexText:
..summary:Shortcut for $getFibre(.., ESA_Text)$.
..cat:Index
..signature:indexText(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA Index Fibres.ESA_Text@ fibre (original text).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_Text>::Type & indexText(Index<TText, TSpec> &index) { return getFibre(index, Fibre_Text()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_Text>::Type & indexText(Index<TText, TSpec> const &index) { return getFibre(index, Fibre_Text()); }

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename StringSetLimits<TText const>::Type
	stringSetLimits(Index<TText, TSpec> &) {
		return Nothing();
	}

	template <typename TText, typename TSpec>
	inline typename StringSetLimits<TText const>::Type
	stringSetLimits(Index<TText, TSpec> const &) {
		return Nothing();
	}

	template <typename TString, typename TSSetSpec, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSSetSpec> const >::Type &
	stringSetLimits(Index<StringSet<TString, TSSetSpec>, TSpec> &index) {
		return stringSetLimits(indexText(index));
	}

	template <typename TString, typename TSSetSpec, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSSetSpec> const >::Type &
	stringSetLimits(Index<StringSet<TString, TSSetSpec>, TSpec> const &index) {
		return stringSetLimits(indexText(index));
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexRawText:
..summary:Shortcut for $getFibre(.., ESA_RawText)$.
..cat:Index
..signature:indexRawText(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA Index Fibres.ESA_RawText@ fibre (concatenated input text).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_RawText>::Type & indexRawText(Index<TText, TSpec> &index) { return getFibre(index, Fibre_RawText()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_RawText>::Type & indexRawText(Index<TText, TSpec> const &index) { return getFibre(index, Fibre_RawText()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexSA:
..summary:Shortcut for $getFibre(.., ESA_SA)$.
..cat:Index
..signature:indexSA(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA Index Fibres.ESA_SA@ fibre (suffix array).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_SA>::Type & indexSA(Index<TText, TSpec> &index) { return getFibre(index, Fibre_SA()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_SA>::Type & indexSA(Index<TText, TSpec> const &index) { return getFibre(index, Fibre_SA()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexRawSA:
..summary:Shortcut for $getFibre(.., ESA_RawSA)$.
..cat:Index
..signature:indexRawSA(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA Index Fibres.ESA_RawSA@ fibre (suffix array).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_RawSA>::Type indexRawSA(Index<TText, TSpec> &index) { return getFibre(index, Fibre_RawSA()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_RawSA>::Type indexRawSA(Index<TText, TSpec> const &index) { return getFibre(index, Fibre_RawSA()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexLCP:
..summary:Shortcut for $getFibre(.., ESA_LCP)$.
..cat:Index
..signature:indexLCP(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA Index Fibres.ESA_LCP@ fibre (lcp table).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_LCP>::Type & indexLCP(Index<TText, TSpec> &index) { return getFibre(index, Fibre_LCP()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_LCP>::Type & indexLCP(Index<TText, TSpec> const &index) { return getFibre(index, Fibre_LCP()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexLCPE:
..summary:Shortcut for $getFibre(.., ESA_LCPE)$.
..cat:Index
..signature:indexLCPE(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA Index Fibres.ESA_LCPE@ fibre (enhanced lcp table).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_LCPE>::Type & indexLCPE(Index<TText, TSpec> &index) { return getFibre(index, Fibre_LCPE()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_LCPE>::Type & indexLCPE(Index<TText, TSpec> const &index) { return getFibre(index, Fibre_LCPE()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexBWT:
..summary:Shortcut for $getFibre(.., ESA_BWT)$.
..cat:Index
..signature:indexBWT(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA Index Fibres.ESA_BWT@ fibre (Burrows-Wheeler table).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_BWT>::Type & indexBWT(Index<TText, TSpec> &index) { return getFibre(index, Fibre_BWT()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_BWT>::Type & indexBWT(Index<TText, TSpec> const &index) { return getFibre(index, Fibre_BWT()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexChildTab:
..summary:Shortcut for $getFibre(.., ESA_ChildTab)$.
..cat:Index
..signature:indexChildTab(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA Index Fibres.ESA_ChildTab@ fibre (child table).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Fibre_ChildTab>::Type & indexChildTab(Index<TText, TSpec> &index) { return getFibre(index, Fibre_ChildTab()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Fibre_ChildTab>::Type & indexChildTab(Index<TText, TSpec> const &index) { return getFibre(index, Fibre_ChildTab()); }

}

#endif
