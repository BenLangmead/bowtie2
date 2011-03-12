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
  $Id: shape_base.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SHAPE_BASE_H
#define SEQAN_HEADER_SHAPE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	template <unsigned q>
	struct FixedShape {};
	typedef FixedShape<0> SimpleShape;

	template <typename TSpec>
	struct FixedGappedShape {};
	typedef FixedGappedShape<Default> GappedShape;


/**
.Class.Shape:
..cat:Index
..summary:Stores hash value and shape for an ungapped or gapped q-gram.
..signature:Shape<TValue, TSpec>
..param.TValue:The @Metafunction.Value@ type of the string the shape is applied to (e.g. $Dna$).
..param.TSpec:The specializing type.
...default:@Spec.SimpleShape@, for ungapped q-grams.
..remarks:The @Metafunction.ValueSize@ of Shape is the ValueSize of TValue which is the alphabet size.
..remarks:To get the span or the weight of a shape call @Function.length@ or @Function.weight@.
.Memfunc.Shape#Shape:
..class:Class.Shape
..summary:Constructor
..signature:Shape<TValue, TSpec> ()
..signature:Shape<TValue, TSpec> (shape)
..param.shape:Other Shape object. (copy constructor)
*/
	template <typename TValue = Dna, typename TSpec = SimpleShape>
	class Shape;

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Shape
	template <typename TValue, typename TSpec>
	struct Value<Shape<TValue,TSpec> >
	{
		typedef unsigned Type;
	};

///.Metafunction.Size.param.T.type:Class.Shape
	template <typename TValue, typename TSpec>
	struct Size<Shape<TValue,TSpec> >
	{
		typedef unsigned Type;
	};

///.Metafunction.LENGTH.param.T.type:Class.Shape
    template <typename TValue, unsigned q>
	struct LENGTH< Shape<TValue, FixedShape<q> > >
	{
		enum { VALUE = q };
	};

///.Metafunction.WEIGHT.param.T.type:Class.Shape
    template <typename TValue, unsigned q>
	struct WEIGHT< Shape<TValue, FixedShape<q> > >
	{
		enum { VALUE = q };
	};

///.Metafunction.ValueSize.param.T.type:Class.Shape
	template <typename TValue, typename TSpec>
	struct ValueSize< Shape<TValue, TSpec> > {
		enum { VALUE = Power<
						ValueSize<TValue>::VALUE, 
						WEIGHT< Shape<TValue, TSpec> >::VALUE >::VALUE };
	};


//////////////////////////////////////////////////////////////////////////////

/**
.Spec.SimpleShape:
..cat:Index
..summary:A variable length ungapped shape (also called q-gram or k-mer).
..general:Class.Shape
..signature:Shape<TValue, SimpleShape>
..param.TValue:The @Metafunction.Value@ type of the string the shape is applied to (e.g. $Dna$).
..remarks:A SimpleShape must be resized first to a valid length. To do so, call @Function.resize@.
..see:Spec.FixedShape
*/

	//////////////////////////////////////////////////////////////////////////////
	// ungapped shape with variable length
	//////////////////////////////////////////////////////////////////////////////

	template <typename TValue>
	class Shape<TValue, SimpleShape>
	{
	public:
//____________________________________________________________________________

		unsigned					span;
		typename Value<Shape>::Type	hValue;
		typename Value<Shape>::Type	XValue;
		typename Value<Shape>::Type	leftFactor;
		typename Value<Shape>::Type	leftFactor2;
		TValue						leftChar;
//____________________________________________________________________________
		
/**
.Memfunc.SimpleShape#Shape:
..class:Spec.SimpleShape
..summary:Constructor
..signature:Shape<TValue, SimpleShape> ()
..signature:Shape<TValue, SimpleShape> (shape)
..signature:Shape<TValue, SimpleShape> (q)
..param.shape:Other Shape object. (copy constructor)
..param.q:Length of the ungapped q-gram.
*/
		Shape() {}

		Shape(unsigned _span)
		{
			resize(*this, _span);
		}

		template <unsigned q>
		Shape(Shape<TValue, FixedShape<q> > const &other)
		{
			*this = other;
		}	

//____________________________________________________________________________

		template <unsigned q>
		inline Shape &
		operator=(Shape<TValue, FixedShape<q> > const &other)
		{
			span = other.span;
			hValue = other.hValue;
			XValue = other.XValue;
			leftFactor = other.leftFactor;
			leftFactor2 = other.leftFactor2;
			leftChar = other.leftChar;
			return *this;
		}
	};

	//////////////////////////////////////////////////////////////////////////////
	// ungapped shape with fixed length q
	//////////////////////////////////////////////////////////////////////////////

/**
.Spec.FixedShape:
..cat:Index
..summary:A fixed length ungapped shape (also called q-gram or k-mer).
..general:Class.Shape
..signature:Shape<TValue, FixedShape<q> >
..param.TValue:The @Metafunction.Value@ type of the sequence the shape is applied to (e.g. $Dna$).
..param.q:The length of the shape.
*/

	template <typename TValue, unsigned q>
	class Shape<TValue, FixedShape<q> >
	{
	public:
//____________________________________________________________________________

		enum { span = q };
		enum { leftFactor = Power<ValueSize<TValue>::VALUE, q - 1>::VALUE };
		enum { leftFactor2 = (Power<ValueSize<TValue>::VALUE, q>::VALUE - 1) / (ValueSize<TValue>::VALUE - 1) };
		// Sigma^(q-1) + Sigma^(q-2) + ... + Sigma + 1

		typename Value<Shape>::Type	hValue;		// current hash value
		typename Value<Shape>::Type	XValue;		// Sum_{i=0..q-1} (x_i + 1)
		TValue						leftChar;	// left-most character
//____________________________________________________________________________
		
	};



//////////////////////////////////////////////////////////////////////////////

///.Function.value.param.object.type:Class.Shape
	template <typename TValue, typename TSpec>
	inline typename Value< Shape<TValue, TSpec> >::Type
	value(Shape<TValue, TSpec> &me)
	{
		return me.hValue;
	}

//____________________________________________________________________________

///.Function.length.param.object.type:Class.Shape
	template <typename TValue, typename TSpec>
	inline typename Size< Shape<TValue, TSpec> >::Type
	length(Shape<TValue, TSpec> const &me)
	{
	SEQAN_CHECKPOINT
		return me.span;
	}

//____________________________________________________________________________

/**.Function.weight:
..cat:Index
..summary:Number of relevant positions in a shape.
..signature:weight(shape)
..param.shape:Shape object for which the number of relevant positions is determined.
...type:Class.Shape
..returns:Number of relevant positions.
..remarks.text:For ungapped shapes the return value is the result of the @Function.length@ function.
For gapped shapes this is the number of '1's.
*/
	template <typename TValue, typename TSpec>
	inline typename Size< Shape<TValue, TSpec> >::Type
	weight(Shape<TValue, TSpec> const &me)
	{
	SEQAN_CHECKPOINT
		return length(me);
	}

//____________________________________________________________________________

///.Function.resize.param.object.type:Spec.SimpleShape
	template <typename TValue, typename TSize>
	inline typename Size< Shape<TValue, SimpleShape> >::Type
	resize(Shape<TValue, SimpleShape> & me, TSize new_length)
	{
	SEQAN_CHECKPOINT
		me.leftFactor = _intPow((unsigned)ValueSize<TValue>::VALUE, new_length - 1);
		me.leftFactor2 = (_intPow((unsigned)ValueSize<TValue>::VALUE, new_length) - 1) / (ValueSize<TValue>::VALUE - 1);
		return me.span = new_length;
	}

//____________________________________________________________________________

/**.Function.hash:
..cat:Index
..summary:Computes a (lower) hash value for a shape applied to a sequence.
..signature:hash(shape, it)
..signature:hash(shape, it, charsLeft)
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the shape.
..param.charsLeft:The distance of $it$ to the string end. 
If $charsLeft$ is smaller than the shape's span, the hash value corresponds to the smallest shape beginning with $charsLeft$ characters.
..returns:Hash value of the shape.
*/

	template <typename TValue, typename TIter>
	typename Value< Shape<TValue, SimpleShape> >::Type
	hash(Shape<TValue, SimpleShape> &me, TIter it)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, SimpleShape> >::Type	THValue;
		typedef typename Size< Shape<TValue, SimpleShape> >::Type	TSize;

		me.hValue = ordValue(me.leftChar = *it);
		for(TSize i = 1; i < me.span; ++i) {
			++it;
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
		}
		return me.hValue;
	}

//____________________________________________________________________________
// fixed ungapped shapes

	// loop unrolling ...
	template <typename THValue, typename TValue, typename TIter>
	inline THValue
	_hashFixedShape(THValue hash, TIter &, TValue const, FixedShape<1> const) {
		return hash;
	}

	template <typename THValue, typename TValue, typename TIter, unsigned q>
	inline THValue
	_hashFixedShape(THValue hash, TIter &it, TValue const, FixedShape<q> const) {
		++it;
		return _hashFixedShape(
			hash * ValueSize<TValue>::VALUE + ordValue((TValue)*it),
			it, TValue(), FixedShape<q - 1>());
	}

	// ... for fixed ungapped shapes
	template <typename TValue, unsigned q, typename TIter>
	inline typename Value< Shape<TValue, FixedShape<q> > >::Type
	hash(Shape<TValue, FixedShape<q> > &me, TIter it)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, FixedShape<q> > >::Type	THValue;
		typedef typename Size< Shape<TValue, FixedShape<q> > >::Type	TSize;

		me.hValue = ordValue(me.leftChar = *it);
		return me.hValue = _hashFixedShape(me.hValue, it, TValue(), FixedShape<q>());
	}


	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = ordValue(me.leftChar = *it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
			}
		} else
			return me.hValue = 0;

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			me.hValue *= ValueSize<TValue>::VALUE;
		return me.hValue;
	}

//____________________________________________________________________________
// Tuple -> fixed ungapped shapes

	template <typename THValue, typename TValue, typename TTValue, unsigned SIZE, typename TCompressed>
	inline THValue
	_hashTuple2FixedShape(
		THValue const, 
		Tuple<TTValue, SIZE, TCompressed> const &tuple,
		TValue const,
		FixedShape<1> const) 
	{
		return ordValue(tuple[0]);
	}

	template <typename THValue, typename TValue, typename TTValue, unsigned SIZE, typename TCompressed, unsigned q>
	inline THValue
	_hashTuple2FixedShape(
		THValue const, 
		Tuple<TTValue, SIZE, TCompressed> const &tuple,
		TValue const,
		FixedShape<q> const) 
	{
		return _hashTuple2FixedShape(THValue(), tuple, TValue(), FixedShape<q - 1>()) 
			* ValueSize<TValue>::VALUE + ordValue(tuple[q-1]);
	}

	// ... for fixed ungapped shapes
	template <
		typename TValue,
		typename TTValue, 
		unsigned SIZE, 
		unsigned q>
	typename Value< Shape<TValue, FixedShape<q> > >::Type
	hash(
		Shape<TValue, FixedShape<q> > &me, 
		Tuple<TTValue, SIZE, Compressed> const &tuple)
	{
	SEQAN_CHECKPOINT
		if (ValueSize<TValue>::VALUE == (1 << BitsPerValue<TTValue>::VALUE))
			if (q == SIZE)
				return tuple.i;
			else
				return tuple >> (q - SIZE);
		else
			return me.hValue = _hashTuple2FixedShape(me.hValue, tuple, TValue(), FixedShape<q>());
	}

	template <
		typename TValue,
		typename TTValue, 
		unsigned SIZE, 
		typename TCompressed, 
		unsigned q>
	typename Value< Shape<TValue, FixedShape<q> > >::Type
	hash(
		Shape<TValue, FixedShape<q> > &me, 
		Tuple<TTValue, SIZE, TCompressed> const &tuple)
	{
	SEQAN_CHECKPOINT
		return me.hValue = _hashTuple2FixedShape(me.hValue, tuple, TValue(), FixedShape<q>());
	}

//____________________________________________________________________________

/**.Function.hashUpper:
..cat:Index
..summary:Computes an upper hash value for a shape applied to a sequence.
..signature:hashUpper(shape, it, charsLeft)
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the shape.
..param.charsLeft:The distance of $it$ to the string end. 
If $charsLeft$ is smaller than the shape's span, the hash value corresponds to the biggest shape beginning with $charsLeft$ characters + 1.
..returns:Upper hash value of the shape.
*/

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hashUpper(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = ordValue(me.leftChar = *it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
			}
			++me.hValue;
		} else
			me.hValue = 1;

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			me.hValue *= ValueSize<TValue>::VALUE;
		return me.hValue;
	}

//____________________________________________________________________________

/**
.Function.hashNext:
..cat:Index
..summary:Computes the hash value for the adjacent shape.
..signature:hashNext(shape, it)
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the adjacent shape.
..returns:Hash value of the q-gram.
..remarks:@Function.hash@ has to be called before.
*/

	template <typename TValue, typename TSpec, typename TIter>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hashNext(Shape<TValue, TSpec> &me, TIter &it)
	{
	SEQAN_CHECKPOINT
		// remove first, shift left, and add next character
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;
		me.hValue = 
			(me.hValue - ordValue(me.leftChar) * (THValue)me.leftFactor) * ValueSize<TValue>::VALUE
			+ ordValue((TValue)*(it + (THValue)me.span - 1));
		me.leftChar = *it;
		return me.hValue;
	}

//____________________________________________________________________________

/**.Function.hash2:
..cat:Index
..summary:Computes a unique hash value of a shape, even if it is shorter than its span.
..signature:hash2(shape, it, charsLeft)
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the shape.
..param.charsLeft:The distance of $it$ to the string end. 
..returns:Hash value of the shape.
*/

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash2(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = me.XValue = ordValue(me.leftChar = *it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				// update sum of x_i
				me.XValue += ordValue((TValue)*it);
				// shift hash
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + me.XValue;
			}
		} else
			return me.hValue = me.XValue = 0;

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + me.XValue;
		return me.hValue += iEnd;
	}

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash2Upper(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		THValue hValue, XValue;
		TSize i = 0;
		if (iEnd > 0) {
			hValue = XValue = ordValue((TValue)*it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				// update sum of x_i
				XValue += ordValue((TValue)*it);
				// shift hash
				hValue = hValue * ValueSize<TValue>::VALUE + XValue;
			}
		} else
			hValue = XValue = 0;

		if (charsLeft <= me.span) {
			++XValue;
			++hValue;
		}

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			hValue = hValue * ValueSize<TValue>::VALUE + XValue;
		return hValue += iEnd;
	}

//____________________________________________________________________________

/**
.Function.hash2Next:
..cat:Index
..summary:Computes a unique hash value for the adjacent shape, even if it is shorter than q.
..signature:hash2Next(shape, it)
..param.shape:Shape to be used for hashing the q-gram.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the adjacent shape.
..returns:Hash value of the shape.
..remarks:@Function.hash@ has to be called before with $shape$ on the left adjacent q-gram.
*/

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash2Next(Shape<TValue, TSpec> &me, TIter &it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		// remove first, shift left, and add next character
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		if (charsLeft >= me.span) {
			// update sum of x_i
			me.XValue = me.XValue + ordValue((TValue)*(it + me.span - 1)) - ordValue(me.leftChar);
			// shift hash
			me.hValue = (me.hValue - ordValue(me.leftChar) * (THValue)me.leftFactor2) * ValueSize<TValue>::VALUE + me.XValue
						- me.span * (ValueSize<TValue>::VALUE - 1);
		} else {
			// update sum of x_i
			me.XValue -= ordValue(me.leftChar);
			// shift hash
			me.hValue = (me.hValue - ordValue(me.leftChar) * (THValue)me.leftFactor2) * ValueSize<TValue>::VALUE + me.XValue
				        - charsLeft * (ValueSize<TValue>::VALUE - 1) - ValueSize<TValue>::VALUE;
		}

		me.leftChar = *it;
		return me.hValue;
	}

	template <typename TString, typename THash>
	inline void unhash(TString &result, THash hash, unsigned q)
	{
	SEQAN_CHECKPOINT
		typedef typename Value<TString>::Type	TValue;

		resize(result, q);
		for (unsigned i = q; i > 0; ) 
		{
			result[--i] = (TValue)(hash % ValueSize<TValue>::VALUE);
			hash /= ValueSize<TValue>::VALUE;
		}
	}

}	// namespace seqan

#endif
