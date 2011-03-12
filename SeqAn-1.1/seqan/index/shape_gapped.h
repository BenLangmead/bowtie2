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
  $Id: shape_gapped.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SHAPE_GAPPED_H
#define SEQAN_HEADER_SHAPE_GAPPED_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// HardwiredShape allows compiler-time defined gapped shape

/**
.Class.HardwiredShape:
..cat:Index
..summary:A structure to define a fixed gapped shape.
..signature:HardwiredShape<P1, P2, ..., Pn>
..param.P1, P2, ..., Pn:Px is the distance of the x'th '1' to the next '1' in the shape.
...remarks:At most 20 parameters are allowed, so the maximal shape weight is 21.
..remarks:You can use this structure to define your one gapped shapes in conjunction with @Spec.FixedGappedShape@.
..note:The shape $1100101$ corresponds to $HardwiredShape<1,3,2>$.
..note:The following predefined shapes are already available in $seqan/index/shape_predefined.h$:
..file:../projects/library/seqan/index/shape_predefined.h
*/

	// Pxx = spaces between '1's
	template <
		int P00 = 0, int P01 = 0, int P02 = 0, int P03 = 0, int P04 = 0, 
		int P05 = 0, int P06 = 0, int P07 = 0, int P08 = 0, int P09 = 0,
		int P10 = 0, int P11 = 0, int P12 = 0, int P13 = 0, int P14 = 0,
		int P15 = 0, int P16 = 0, int P17 = 0, int P18 = 0, int P19 = 0	
	>
	struct HardwiredShape {
		static const int DIFFS[];
	};

	template <
		int P00, int P01, int P02, int P03, int P04, 
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19	
	>
	const int HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	>::DIFFS[] = {
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19, 0 
	};


//////////////////////////////////////////////////////////////////////////////
// LENGTH meta-function for fixed gapped shapes

	template <>
	struct LENGTH< HardwiredShape<
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0> >
	{
		enum { VALUE = 1 };
	};

	template <
		int P00, int P01, int P02, int P03, int P04, 
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19	
	>
	struct LENGTH< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19> >
	{
		enum { VALUE = LENGTH< HardwiredShape<
			P01,P02,P03,P04,P05,
			P06,P07,P08,P09,P10,
			P11,P12,P13,P14,P15,
			P16,P17,P18,P19, 0 > >::VALUE + P00 };
	};

	template <
	    typename TValue,
		int P00, int P01, int P02, int P03, int P04, 
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19	
	>
	struct LENGTH< Shape<TValue, FixedGappedShape< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	> > > >:
	LENGTH< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	> > {};


//////////////////////////////////////////////////////////////////////////////
// WEIGHT meta-function for fixed gapped shapes

	template <>
	struct WEIGHT< HardwiredShape<
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0> >
	{
		enum { VALUE = 1 };
	};

	template <
		int P00, int P01, int P02, int P03, int P04, 
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19	
	>
	struct WEIGHT< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19> >
	{
		enum { VALUE = WEIGHT< HardwiredShape<
			P01,P02,P03,P04,P05,
			P06,P07,P08,P09,P10,
			P11,P12,P13,P14,P15,
			P16,P17,P18,P19, 0 > >::VALUE + 1 };
	};

	template <
	    typename TValue,
		int P00, int P01, int P02, int P03, int P04, 
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19	
	>
	struct WEIGHT< Shape<TValue, FixedGappedShape< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	> > > >:
	WEIGHT< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	> > {};


//////////////////////////////////////////////////////////////////////////////

/**
.Spec.GappedShape:
..cat:Index
..summary:A variable gapped shape.
..general:Class.Shape
..signature:Shape<TValue, GappedShape>
..param.TValue:The @Metafunction.Value@ type of the string the shape is applied to (e.g. $Dna$).
..remarks:A GappedShape must be initialized first with a valid shape. To do so, call @Function.stringToShape@.
..see:Spec.FixedGappedShape
*/

	//////////////////////////////////////////////////////////////////////////////
	// variable gapped shape
	//////////////////////////////////////////////////////////////////////////////

	template <typename TValue>
	class Shape<TValue, GappedShape>
	{
	public:
//____________________________________________________________________________

		unsigned span;
		unsigned weight;
		String<int> diffs;
	
		typename Value<Shape>::Type	hValue;		// current hash value
//____________________________________________________________________________

/**
.Memfunc.GappedShape#Shape:
..class:Spec.GappedShape
..summary:Constructor
..signature:Shape<TValue, GappedShape> ()
..signature:Shape<TValue, GappedShape> (q)
..signature:Shape<TValue, GappedShape> (shape)
..signature:Shape<TValue, GappedShape> (predefined)
..param.q:Creates an ungapped q-gram.
..param.shape:Any other gapped/ungapped shape.
..param.predefined:Any instance of a predefined shape spec (e.g. $ShapePatternHunter$).
..see:Class.HardwiredShape
*/
		Shape()	{}

		// c'tor for ungapped shapes
		Shape(unsigned _span):
			span(_span),
			weight(_span)
		{
		SEQAN_CHECKPOINT
			resize(diffs, _span);
			for(unsigned i = 0; i < _span; ++i)
				diffs[i] = 1;
		}

		Shape(Shape const &other):
			span(other.span),
			weight(other.weight),
			diffs(other.diffs),
			hValue(other.hValue) {}	

		template <typename TSpec>
		Shape(Shape<TValue, TSpec> const &other)
		{
			*this = other;
		}	

		template <typename TSpec>
		Shape(FixedGappedShape<TSpec> const &other)
		{
			*this = other;
		}

		template <typename TStringValue, typename TSpec>
		Shape(String<TStringValue, TSpec> const &bitmap)
		{
			*this = bitmap;
		}	

//____________________________________________________________________________

		template <unsigned q>
		inline Shape &
		operator=(Shape<TValue, FixedShape<q> > const &other)
		{
			span = length(other);
			weight = weight(other);
			resize(diffs, weight);
			for(unsigned i = 1; i < weight; ++i)
				diffs[i] = 1;
			hValue = other.hValue;
			return *this;
		}

		template <typename TSpec>
		inline Shape &
		operator=(Shape<TValue, FixedGappedShape<TSpec> > const &other)
		{
			span = other.span;
			weight = other.weight;
			diffs = other.diffs;
			hValue = other.hValue;
			return *this;
		}

		template <typename TSpec>
		inline Shape &
		operator=(FixedGappedShape<TSpec> const)
		{
			typedef Shape<TValue, FixedGappedShape<TSpec> > TShape;
			return *this = TShape();
		}

		template <typename TStringValue, typename TSpec>
		inline Shape &
		operator=(String<TStringValue, TSpec> const &bitmap)
		{
			stringToShape(*this, bitmap);
			return *this;
		}
	};


//////////////////////////////////////////////////////////////////////////////

/**
.Spec.FixedGappedShape:
..cat:Index
..summary:A fixed gapped shape.
..general:Class.Shape
..signature:Shape<TValue, FixedGappedShape<TSpec> >
..param.TValue:The @Metafunction.Value@ type of the string the shape is applied to (e.g. $Dna$).
..param.TSpec:A structure to store the shape at compile-time.
...type:Class.HardwiredShape
..remarks:There are predefined shapes in $index/shape_predefined.h$.
You can simply use them with $Shape<TValue, ShapePatternHunter>$ for example.
..see:Class.HardwiredShape
*/

	//////////////////////////////////////////////////////////////////////////////
	// fixed gapped shape
	//////////////////////////////////////////////////////////////////////////////

	template <typename TValue, typename TSpec>
	class Shape<TValue, FixedGappedShape<TSpec> >
	{
	public:
//____________________________________________________________________________

		typedef FixedGappedShape<TSpec>	TShapeSpec;

		enum { span = LENGTH<Shape>::VALUE };
		enum { weight = WEIGHT<Shape>::VALUE };
		const int *diffs;
	
		typename Value<Shape>::Type	hValue;		// current hash value
//____________________________________________________________________________

		Shape():
			diffs(TSpec::DIFFS) {}

		Shape(Shape const &other):
			diffs(other.diffs),	
			hValue(other.hValue) {}
//____________________________________________________________________________

		inline Shape &
		operator=(Shape const &other)
		{
			hValue = other.hValue;
		}
	};

//////////////////////////////////////////////////////////////////////////////

	template <typename TValue, typename TSpec>
	inline typename Size< Shape<TValue, FixedGappedShape<TSpec> > >::Type
	weight(Shape<TValue, FixedGappedShape<TSpec> > const & me)
	{
	SEQAN_CHECKPOINT
		return me.weight;
	}

//____________________________________________________________________________

	template <typename TValue, typename TIter>
	inline typename Value< Shape<TValue, GappedShape> >::Type
	hash(Shape<TValue, GappedShape> &me, TIter it)	
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, GappedShape> >::Type	THValue;
		typedef typename Size< Shape<TValue, GappedShape> >::Type	TSize;

		me.hValue = ordValue((TValue)*it);
		TSize iEnd = me.weight - 1;
		for(TSize i = 0; i < iEnd; ++i) {
			goFurther(it, me.diffs[i]);
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
		}
		return me.hValue;
	}

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, FixedGappedShape<TSpec> > >::Type
	hash(Shape<TValue, FixedGappedShape<TSpec> > &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, FixedGappedShape<TSpec> > >::Type	THValue;

		TSize iEnd = me.weight;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = ordValue((TValue)*it);
			--iEnd;
			for(; i < iEnd; ++i) {
				goFurther(it, me.diffs[i]);
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
			}
			++i;
		} else
			return me.hValue = 0;

		// fill shape with zeros
		for(; i < (TSize)me.weight; ++i)
			me.hValue *= ValueSize<TValue>::VALUE;
		return me.hValue;
	}

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, FixedGappedShape<TSpec> > >::Type
	hashUpper(Shape<TValue, FixedGappedShape<TSpec> > &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, FixedGappedShape<TSpec> > >::Type	THValue;

		TSize iEnd = me.weight;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = ordValue((TValue)*it);
			--iEnd;
			for(; i < iEnd; ++i) {
				goFurther(it, me.diffs[i]);
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
			}
			++i;
			++me.hValue;
		} else
			return me.hValue = 1;

		// fill shape with zeros
		for(; i < (TSize)me.weight; ++i)
			me.hValue *= ValueSize<TValue>::VALUE;
		return me.hValue;
	}

//____________________________________________________________________________

	template <typename THValue, typename TValue, typename TIter>
	inline THValue
	_hashHardwiredShape(THValue hash, TIter &, TValue const, HardwiredShape<
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0 > const)
	{
		return hash;
	}

	template <
		         int P01, int P02, int P03, int P04,
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19,
		typename THValue, typename TValue, typename TIter
	>
	inline THValue
	_hashHardwiredShape(THValue hash, TIter &it, TValue const, HardwiredShape<
		 1 ,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 > const)
	{
		++it;
		return _hashHardwiredShape(hash * ValueSize<TValue>::VALUE + ordValue((TValue)*it),
			it, TValue(), HardwiredShape<
				P01,P02,P03,P04,P05,
				P06,P07,P08,P09,P10,
				P11,P12,P13,P14,P15,
				P16,P17,P18,P19, 0 >());
	}

	template <
		int P00, int P01, int P02, int P03, int P04,
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19,
		typename THValue, typename TValue, typename TIter
	>
	inline THValue
	_hashHardwiredShape(THValue hash, TIter &it, TValue const, HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 > const)
	{
		it += P00;
		return _hashHardwiredShape(hash * ValueSize<TValue>::VALUE + ordValue((TValue)*it),
			it, TValue(), HardwiredShape<
				P01,P02,P03,P04,P05,
				P06,P07,P08,P09,P10,
				P11,P12,P13,P14,P15,
				P16,P17,P18,P19, 0 >());
	}

	template <
		int P00, int P01, int P02, int P03, int P04,
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19,
		typename TValue, typename TIter
	>
	inline typename Value< Shape<TValue, FixedGappedShape< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	> > > >::Type
	hash(Shape<TValue, FixedGappedShape< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	> > > &me, TIter it)
	{
	SEQAN_CHECKPOINT
		typedef HardwiredShape<
			P00,P01,P02,P03,P04,
			P05,P06,P07,P08,P09,
			P10,P11,P12,P13,P14,
			P15,P16,P17,P18,P19 >								TSpec;
		typedef FixedGappedShape<TSpec>							TShape;
		typedef typename Value< Shape<TValue, TShape> >::Type	THValue;

		me.hValue = (THValue)ordValue((TValue)*it);
		return me.hValue = _hashHardwiredShape(me.hValue, it, TValue(), TSpec());
	}

//____________________________________________________________________________

	template <typename TValue, typename TSpec, typename TIter>
	inline typename Value< Shape<TValue, FixedGappedShape<TSpec> > >::Type
	hashNext(Shape<TValue, FixedGappedShape<TSpec> > &me, TIter it)	
	{
	SEQAN_CHECKPOINT
		return hash(me, it);
	}


//____________________________________________________________________________

/**.Function.stringToShape:
..cat:Index
..summary:Takes a shape given as a string of '1' (relevant position) and '0' 
(irrelevant position) and converts it into a Shape object.
..signature:stringToShape(shape, bitmap)
..param.shape:Shape object that is manipulated.
...type:Spec.GappedShape
..param.bitmap:A character string of '1' and '0' representing relevant and irrelevant positions (blanks) respectively.
...remarks:This string must begin with a '1'.
...type:Class.String
*/

	template <typename TValue, typename TSpec, typename TShapeString>
	inline void
	stringToShape(
		Shape<TValue, FixedGappedShape<TSpec> > &me, 
		TShapeString const &bitmap)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TShapeString const>::Type		TIter;
		typedef typename Iterator<String<int> >::Type			TShapeIter;

		me.span = length(bitmap);

		unsigned oneCount = 0;
		TIter it = begin(bitmap, Standard());
		TIter itEnd = end(bitmap, Standard());
		for(; it != itEnd; ++it)
			if (*it == '1')
				++oneCount;

		me.weight = oneCount;
		resize(me.diffs, oneCount);

		unsigned diff = 0;
		it = begin(bitmap, Standard());
		TShapeIter itS = begin(me.diffs, Standard());
		for(; it != itEnd; ++it) {
			if (*it == '1') {
				*itS = diff;
				++itS;
				diff = 0;
			}
			++diff;
		}
	}

}	// namespace seqan

#endif
