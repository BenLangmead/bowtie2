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
  $Id: basic_alphabet_simple.h,v 1.1 2008/08/25 16:20:02 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_ALPHABET_SIMPLE_H
#define SEQAN_HEADER_BASIC_ALPHABET_SIMPLE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//Class that is used for various simple value types
//////////////////////////////////////////////////////////////////////////////

/**
.Class.SimpleType:
..cat:Basic
..summary:Implementation for "simple" types.
..signature:SimpleType<TValue, TSpec>
..param.TValue:Type that stores the values of an instance.
...remarks:TValue must be a simple type.
...metafunction:Metafunction.Value
..param.TSpec:Specialization tag.
...metafunction:Metafunction.Spec
..remarks:
...text:A "simple type" is a C++ type that can be constructed without constructor,
destructed without destructor and copied without copy constructor or assignment operator.
All basic types (like $char$, $int$ or $float$) are simple. Pointers, references and arrays of
simple types are simple.
POD types ("plain old data types"), that are - simplified spoken - C++-types that already existed in C,
are simple too. 
...text:Arrays of simple types can be copied very fast by memory manipulation routines, 
but the default implementation of functions like @Function.arrayCopyForward@ and @Function.arrayCopy@
are not optimized for simple types this way.
But for classes derived from $SimpleType$, optimized variants of array manipulation functions are applied. 
...text:Note that simple types need not to be derived or specialized from $SimpleType$, but
it could be convenient to do so.
..implements:Concept.Simple Type
*/
template <typename TValue, typename TSpec>
struct SimpleType
{
//____________________________________________________________________________

	TValue value;

//____________________________________________________________________________

	SimpleType() 
	{
SEQAN_CHECKPOINT
	}

//____________________________________________________________________________

	SimpleType(SimpleType const & other)
	{
SEQAN_CHECKPOINT
		assign(*this, other);
	}

	template <typename T> 
	SimpleType(T const & other) 
	{
SEQAN_CHECKPOINT
		assign(*this, other);
	}


//____________________________________________________________________________

	SimpleType & operator=(SimpleType const & other) 
	{ 
SEQAN_CHECKPOINT
		assign(*this, other);
		return *this;
	}
	template <typename T>
	SimpleType & operator=(T const & other) 
	{ 
SEQAN_CHECKPOINT
		assign(*this, other);
		return *this;
	}
//____________________________________________________________________________

	~SimpleType()
	{
SEQAN_CHECKPOINT
	}
//____________________________________________________________________________

	//this cannot be a template since a template would be in conflict to
	//the template c'tor


	operator long() const
	{
SEQAN_CHECKPOINT
		long c;
		assign(c, *this);
		return c;
	}
	operator unsigned long() const
	{
SEQAN_CHECKPOINT
		unsigned long c;
		assign(c, *this);
		return c;
	}
	operator int() const
	{
SEQAN_CHECKPOINT
		int c;
		assign(c, *this);
		return c;
	}
	operator unsigned int() const
	{
SEQAN_CHECKPOINT
		unsigned int c;
		assign(c, *this);
		return c;
	}
	operator short() const
	{
SEQAN_CHECKPOINT
		short c;
		assign(c, *this);
		return c;
	}
	operator unsigned short() const
	{
SEQAN_CHECKPOINT
		unsigned short c;
		assign(c, *this);
		return c;
	}
	operator char() const
	{
SEQAN_CHECKPOINT
		char c;
		assign(c, *this);
		return c;
	}
	operator signed char() const
	{
SEQAN_CHECKPOINT
		signed char c;
		assign(c, *this);
		return c;
	}
	operator unsigned char() const
	{
SEQAN_CHECKPOINT
		unsigned char c;
		assign(c, *this);
		return c;
	}

//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// METAFUNCTIONS
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IsSimple.param.T.type:Class.SimpleType

template <typename TValue, typename TSpec>
struct IsSimple<SimpleType<TValue, TSpec> >
{
	typedef True Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.SimpleType
template <typename TValue, typename TSpec>
struct Value<SimpleType<TValue, TSpec> >
{
	typedef TValue Type;
};

template <typename TValue, typename TSpec>
struct Value<SimpleType<TValue, TSpec> const >
{
	typedef TValue const Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.SimpleType
template <typename TValue, typename TSpec>
struct Spec<SimpleType<TValue, TSpec> >
{
	typedef TSpec Type;
};

template <typename TValue, typename TSpec>
struct Spec<SimpleType<TValue, TSpec> const >
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Iterator<SimpleType<TValue, TSpec>, Standard>
{
	typedef SimpleType<TValue, TSpec> * Type;
//	typedef Iter<SimpleType<TValue, TSpec>, SimpleIterator> * Type;
};

template <typename TValue, typename TSpec>
struct Iterator<SimpleType<TValue, TSpec> const, Standard>
{
	typedef SimpleType<TValue, TSpec> const * Type;
//	typedef Iter<SimpleType<TValue, TSpec> const, SimpleIterator> * Type;
};


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename T, typename TSourceValue, typename TSourceSpec>
inline typename _RemoveConst<TTarget>::Type
convertImpl(Convert<TTarget, T> const,
			SimpleType<TSourceValue, TSourceSpec> const & source_)
{
SEQAN_CHECKPOINT
	typename _RemoveConst<TTarget>::Type target_;
	assign(target_, source_);
	return target_;
}



//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TValue, typename TSpec>
inline TStream &
operator << (TStream & stream, 
			 SimpleType<TValue, TSpec> const & data)
{
SEQAN_CHECKPOINT
	stream << convert<char>(data);
	return stream;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TValue, typename TSpec>
inline TStream &
operator >> (TStream & stream, 
			 SimpleType<TValue, TSpec> & data)
{
SEQAN_CHECKPOINT
	char c;
	stream >> c;
	assign(data, c);
	return stream;
}

//////////////////////////////////////////////////////////////////////////////
// assign
//////////////////////////////////////////////////////////////////////////////

///.Function.assign.param.target.type:Class.SimpleType
///.Function.assign.param.source.type:Class.SimpleType


template <typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   SimpleType<TSourceValue, TSourceSpec> & source)
{
SEQAN_CHECKPOINT
	target.value = source.value;
}
template <typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   SimpleType<TSourceValue, TSourceSpec> const & source)
{
SEQAN_CHECKPOINT
	target.value = source.value;
}

//____________________________________________________________________________

template <typename TTargetValue, typename TTargetSpec, typename TSource>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   TSource & source)
{
SEQAN_CHECKPOINT
	target.value = source;
}
template <typename TTargetValue, typename TTargetSpec, typename TSource>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   TSource const & source)
{
SEQAN_CHECKPOINT
	target.value = source;
}

//____________________________________________________________________________
// Assign Proxy to SimpleType 
//??? Diese Funktionen wurden noetig wegen eines seltsamen VC++-Verhaltens

template <typename TTargetValue, typename TTargetSpec, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   Proxy<TSourceSpec> & source)
{
SEQAN_CHECKPOINT
	target.value = getValue(source);
}

template <typename TTargetValue, typename TTargetSpec, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   Proxy<TSourceSpec> const & source)
{
SEQAN_CHECKPOINT
	target.value = getValue(source);
}

//____________________________________________________________________________
//INTEGRAL TYPES
//note: it is not possible to write a single function here since "assign"
//must be specialized for the first argument at the first place

//int
template <typename TValue, typename TSpec>
inline void 
assign(int & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(int & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//unsigned int
template <typename TValue, typename TSpec>
inline void 
assign(unsigned int & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(unsigned int & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//short
template <typename TValue, typename TSpec>
inline void 
assign(short & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(short & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//unsigned short
template <typename TValue, typename TSpec>
inline void 
assign(unsigned short & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(unsigned short & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//char
template <typename TValue, typename TSpec>
inline void 
assign(char & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(char & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//signed char
template <typename TValue, typename TSpec>
inline void 
assign(signed char & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(signed char & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//unsigned char
template <typename TValue, typename TSpec>
inline void 
assign(unsigned char & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(unsigned char & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// CompareType
//////////////////////////////////////////////////////////////////////////////

/**.Metafunction.CompareType:
..summary:Type to convert other types for comparisons.
..signature:CompareType<TLeft, TRight>::Type
..param.TLeft:Type of the left operand of a comparison.
..param.TRight:Type of the right operand of a comparison.
..return.Type:The Type in which the arguments are converted in order to compare them.
..remarks:Comparisons are for example operators like $==$ or $<$.
..remarks.text:Note that there is no rule that guarantees that $CompareType<T1, T2>::Type$
is the same as $CompareType<T2, T1>::Type$. It is also possible, that only one of these
two types is defined.
..remarks.text:This metafunction is used for the implementation of
comparisons that involve @Class.SimpleType@.
*/
//???TODO: muss geprueft werden, ob diese Metafunktion noch ausgeweitet oder aber versteckt wird.

template <typename TLeft, typename TRight>
struct CompareType;

template <typename T>
struct CompareType<T, T>
{
	typedef T Type;
};

//____________________________________________________________________________

template <typename TValue, typename TSpec, typename TRight>
struct CompareType<SimpleType<TValue, TSpec>, TRight>
{
	typedef TRight Type;
};

//////////////////////////////////////////////////////////////////////////////
// operator ==

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator == (SimpleType<TValue, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator == (TLeft const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator == (SimpleType<TLeftValue, TLeftSpec> const & left_, 
			 SimpleType<TRightValue, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
	typedef SimpleType<TRightValue, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator == (SimpleType<TValue, TSpec> const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return convert<TValue>(left_) == convert<TValue>(right_);
}


template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator == (Proxy<TSpec> const & left_, 
			 SimpleType<TValue, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator == (SimpleType<TValue, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}


//____________________________________________________________________________
// operator !=

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator != (SimpleType<TValue, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator != (TLeft const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator != (SimpleType<TLeftValue, TLeftSpec> const & left_, 
			 SimpleType<TRightValue, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
	typedef SimpleType<TRightValue, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator != (SimpleType<TValue, TSpec> const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return convert<TValue>(left_) != convert<TValue>(right_);
}


template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator != (Proxy<TSpec> const & left_, 
			 SimpleType<TValue, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator != (SimpleType<TValue, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}


//____________________________________________________________________________
// operator <

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator < (SimpleType<TValue, TSpec> const & left_, 
			TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator < (TLeft const & left_, 
			SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator < (SimpleType<TLeftValue, TLeftSpec> const & left_, 
			SimpleType<TRightValue, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
	typedef SimpleType<TRightValue, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator < (SimpleType<TValue, TSpec> const & left_, 
			SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return convert<TValue>(left_) < convert<TValue>(right_);
}


template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator < (Proxy<TSpec> const & left_, 
			 SimpleType<TValue, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator < (SimpleType<TValue, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}


//____________________________________________________________________________
// operator <=

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator <= (SimpleType<TValue, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator <= (TLeft const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator <= (SimpleType<TLeftValue, TLeftSpec> const & left_, 
			 SimpleType<TRightValue, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
	typedef SimpleType<TRightValue, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator <= (SimpleType<TValue, TSpec> const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return convert<TValue>(left_) <= convert<TValue>(right_);
}


template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator <= (Proxy<TSpec> const & left_, 
			 SimpleType<TValue, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator <= (SimpleType<TValue, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}



//____________________________________________________________________________
// operator >

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator > (SimpleType<TValue, TSpec> const & left_, 
			TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator > (TLeft const & left_, 
			SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator > (SimpleType<TLeftValue, TLeftSpec> const & left_, 
			SimpleType<TRightValue, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
	typedef SimpleType<TRightValue, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator > (SimpleType<TValue, TSpec> const & left_, 
			SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return convert<TValue>(left_) > convert<TValue>(right_);
}


template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator > (Proxy<TSpec> const & left_, 
			 SimpleType<TValue, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator > (SimpleType<TValue, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}


//____________________________________________________________________________
// operator >=

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator >= (SimpleType<TValue, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator >= (TLeft const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator >= (SimpleType<TLeftValue, TLeftSpec> const & left_, 
			 SimpleType<TRightValue, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
	typedef SimpleType<TRightValue, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator >= (SimpleType<TValue, TSpec> const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return convert<TValue>(left_) >= convert<TValue>(right_);
}


template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator >= (Proxy<TSpec> const & left_, 
			 SimpleType<TValue, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator >= (SimpleType<TValue, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}


//////////////////////////////////////////////////////////////////////////////

template<typename _T, typename TSpec> 
inline
bool lexLess(SimpleType<_T, TSpec> const &_Left, SimpleType<_T, TSpec> const &_Right)
{	// return lexicographical _Left < _Right
	typedef typename _MakeUnsigned<_T>::Type TUnsigned;
    return (TUnsigned)(_Left.value) < (TUnsigned)(_Right.value);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec> &
operator ++ (SimpleType<TValue, TSpec> & me)
{
	++me.value;
	return me;
}
template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec>
operator ++ (SimpleType<TValue, TSpec> & me,
			 int)
{
	SimpleType<TValue, TSpec> dummy = me;
	++me.value;
	return dummy;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec> &
operator -- (SimpleType<TValue, TSpec> & me)
{
	--me.value;
	return me;
}
template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec>
operator -- (SimpleType<TValue, TSpec> & me,
			 int)
{
	SimpleType<TValue, TSpec> dummy = me;
	--me.value;
	return dummy;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Dna:
..cat:Alphabets
..summary:Alphabet for DNA.
..general:Class.SimpleType
..signature:Dna
..remarks:
...text:The @Metafunction.ValueSize@ of $Dna$ is 4. 
The nucleotides are enumerated this way: $'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3$.
...text:Objects of type $Dna$ can be converted to various other types and vice versa. 
An object that has a value not in ${'A', 'C', 'G', 'T'}$ is converted to $'A'$.
...text:$Dna$ is typedef for $SimpleType<char,_Dna>$, while $_Dna$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
..see:Spec.Dna5
*/
struct _Dna {};
typedef SimpleType<unsigned char,_Dna> Dna;

template <> struct ValueSize< Dna > { enum { VALUE = 4 }; };
template <> struct BitsPerValue< Dna > { enum { VALUE = 2 }; };

//____________________________________________________________________________

/**
.Spec.Dna5:
..cat:Alphabets
..summary:Alphabet for DNA including 'N' character.
..general:Class.SimpleType
..signature:Dna5
..remarks:
...text:The @Metafunction.ValueSize@ of $Dna5$ is 5. 
The nucleotides are enumerated this way: $'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3$. 
The 'N' character ("unkown nucleotide") is encoded by 4.
...text:Objects of type $Dna5$ can be converted to various other types and vice versa. 
An object that has a value not in ${'A', 'C', 'G', 'T'}$ is converted to $'N'$.
...text:$Dna5$ is typedef for $SimpleType<char,_Dna5>$, while $_Dna5$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
*/
struct _Dna5 {};
typedef SimpleType<unsigned char, _Dna5> Dna5;

template <> struct ValueSize< Dna5 > { enum { VALUE = 5 }; };
template <> struct BitsPerValue< Dna5 > { enum { VALUE = 3 }; };

//____________________________________________________________________________

/**
.Spec.Iupac:
..cat:Alphabets
..summary:Iupac code for DNA.
..general:Class.SimpleType
..signature:Iupac
..remarks:
...text:The @Metafunction.ValueSize@ of $Iupac$ is 16. 
The nucleotides are enumerated from 0 to 15 in this order: 
'U'=0, 'T', 'A', 'W', 'C', 'Y', 'M', 'H', 'G', 'K', 'R', 'D', 'S', 'B', 'V', 'N'=15. 
...text:Objects of type $Iupac$ can be converted to various other types and vice versa. 
Unkown values are converted to $'N'$.
...text:$Iupac$ is typedef for $SimpleType<char,_Iupac>$, while $_Iupac$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
*/
struct _Iupac {};
typedef SimpleType<unsigned char, _Iupac> Iupac;

template <> struct ValueSize< Iupac > { enum { VALUE = 16 }; };
template <> struct BitsPerValue< Iupac > { enum { VALUE = 4 }; };


//____________________________________________________________________________

/**
.Spec.AminoAcid:
..cat:Alphabets
..summary:Iupac code for amino acids.
..general:Class.SimpleType
..signature:AminoAcid
..remarks:
...text:The @Metafunction.ValueSize@ of $AminoAcid$ is 24. 
...text:The amino acids are enumerated from 0 to 15 in this order: 
...text:'A'=0, 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'=19.
...text:The remaining 4 symbols are:
...text: 'B'=20 (Aspartic Acid, Asparagine), 'Z'=21 (Glutamic Acid, Glutamine), 'X'=22 (unknown), '*'=23 (terminator)
...text:Objects of type $AminoAcid$ can be converted to $char$ and vice versa. 
Unkown values are converted to $'X'$.
...text:$AminoAcid$ is typedef for $SimpleType<char,_AminoAcid>$, while $_AminoAcid$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
*/
struct _AminoAcid {};
typedef SimpleType<unsigned char, _AminoAcid> AminoAcid;

template <> struct ValueSize< AminoAcid > { enum { VALUE = 24 }; };
template <> struct BitsPerValue< AminoAcid > { enum { VALUE = 5 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign(Ascii & c_target, 
				   Dna const & source)
{
SEQAN_CHECKPOINT
	c_target = _Translate_Table_Dna5_2_Ascii<>::VALUE[source.value];
}
//____________________________________________________________________________

inline void assign(Ascii & c_target, 
				   Dna5 const & source)
{
SEQAN_CHECKPOINT
	c_target = _Translate_Table_Dna5_2_Ascii<>::VALUE[source.value];
}
//____________________________________________________________________________

inline void assign(Ascii & c_target, Iupac const & source)
{
SEQAN_CHECKPOINT
	c_target = _Translate_Table_Iupac_2_Ascii<>::VALUE[source.value];
}
//____________________________________________________________________________

inline void assign(Ascii & c_target, AminoAcid const & source)
{
SEQAN_CHECKPOINT
	c_target = _Translate_Table_AA_2_Ascii<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////
//DNA (4 letters)

template <>
struct CompareType<Dna, Byte> { typedef Dna Type; };
inline void assign(Dna & target, Byte c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Byte_2_Dna<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Dna, Ascii> { typedef Dna Type; };
inline void assign(Dna & target, Ascii c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_Dna<>::VALUE[(unsigned char)c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Dna, Unicode> { typedef Dna Type; };
inline void assign(Dna & target, Unicode c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_Dna<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Dna, Dna5> { typedef Dna Type; };
inline void assign(Dna & target, Dna5 const & c_source)
{
SEQAN_CHECKPOINT
	target.value = c_source.value & 0x03;
}
//____________________________________________________________________________

template <>
struct CompareType<Dna, Iupac> { typedef Dna Type; };
inline void assign(Dna & target, Iupac const & source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Iupac_2_Dna<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////
//DNA (5 letters)

template <>
struct CompareType<Dna5, Byte> { typedef Dna5 Type; };
inline void assign(Dna5 & target, Byte c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Byte_2_Dna5<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Dna5, Ascii> { typedef Dna5 Type; };
inline void assign(Dna5 & target, Ascii c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_Dna5<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Dna5, Unicode> { typedef Dna5 Type; };
inline void assign(Dna5 & target, Unicode c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_Dna5<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Dna5, Iupac> { typedef Dna5 Type; };
inline void assign(Dna5 & target, Iupac const & source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Iupac_2_Dna5<>::VALUE[source.value];
}

//____________________________________________________________________________

template <>
struct CompareType<Dna5, Dna> { typedef Dna Type; };
inline void assign(Dna5 & target, Dna const & c_source)
{
SEQAN_CHECKPOINT
	target.value = c_source.value;
}

//////////////////////////////////////////////////////////////////////////////
//IUPAC (4 bits)

template <>
struct CompareType<Iupac, Byte> { typedef Iupac Type; };
inline void assign(Iupac & target, Byte c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Byte_2_Iupac<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Iupac, Ascii> { typedef Iupac Type; };
inline void assign(Iupac & target, Ascii c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_Iupac<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Iupac, Unicode> { typedef Iupac Type; };
inline void assign(Iupac & target, Unicode c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_Iupac<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

inline void assign(Iupac & target, Dna const & source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Dna5_2_Iupac<>::VALUE[source.value];
}
//____________________________________________________________________________

inline void assign(Iupac & target, Dna5 const & source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Dna5_2_Iupac<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////
//Amino Acid (5 bits)

template <>
struct CompareType<AminoAcid, Byte> { typedef AminoAcid Type; };
inline void assign(AminoAcid & target, Byte c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Byte_2_AA<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<AminoAcid, Ascii> { typedef AminoAcid Type; };
inline void assign(AminoAcid & target, Ascii c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<AminoAcid, Unicode> { typedef AminoAcid Type; };
inline void assign(AminoAcid & target, Unicode c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
