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
  $Id: lexical.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_LEXICAL_H
#define SEQAN_HEADER_LEXICAL_H


namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Switches for prefix ordering mode
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Prefix Order:
..summary:Specify whether a prefix is smaller or greater.
..tag.TagPrefixLess:A prefix is smaller.
...text:For example: $"abc" < "abcde"$.
..tag.TagPrefixGreater:A prefix is greater.
...text:For example: $"abc" > "abcde"$.
..remarks:The default for all comparison functions is $TagPrefixLess$.
*/
struct TagPrefixLess_ {};
typedef Tag<TagPrefixLess_> const TagPrefixLess;

struct TagPrefixGreater_ {};
typedef Tag<TagPrefixGreater_> const TagPrefixGreater;


/**
.Metafunction.DefaultPrefixOrder:
..hidefromindex
..summary:The default prefix order.
..signature:DefaultPrefixOrder<T>::Type
..param.T:Type for which the prefix order is determined.
..returns.param.Type:Prefix order tag for type of $T$.
..see:Tag.Prefix Order
*/
template <typename T>
struct DefaultPrefixOrder
{
	typedef TagPrefixLess Type;
};

//////////////////////////////////////////////////////////////////////////////
// Lexical
//////////////////////////////////////////////////////////////////////////////

/**
.Class.Lexical:
..cat:Basic
..summary:Comparator for lexical comparison.
..signature:Lexical<TSpec>
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...text:This type can be used for specializations of $Lexical$.
...remarks:$TSpec$ is by default interpreted as size-type.
...default:$size_t$
..remarks:
...text:This class implement comparator objects that perform (lexical) comparisons between two sequences.
The result of the comparison is stored in the data members of the instance an can be
accessed by some functions, for example @Function.isLess@ or @Function.isEqual@.
...text:In most cases, there is no need for an explicite use of comparators,
but sometimes this concept provide the opportunity to speed up the code.
..example:
...text:This program compares the strings $str1$ and $str2$:
...code:if (isLess(str1, str2)) //first comparison
{
	//str1 < str2
}
else if (isGreater(str1, str2)) //second comparison
{
	//str1 > str2
}
else
{
	//str == str2
}
...text:Using a comparator, the same program only needs one comparison instead of two:
...code:Lexical <> comparator(str1, str2); //comparison is executed here
if (isLess(comparator))
{
	//str1 < str2
}
else if (lexGreater(comparator))
{
	//str1 > str2
}
else
{
	//str == str2
}
...text:The state of a default constructed $Lexical$ instance is undefined until
it is set by a call of @Function.compare@.
..see:Metafunction.Comparator
*/

template <typename TSpec = size_t>
struct Lexical
{
public:
	typename Size<Lexical>::Type data_lcp;
	char data_compare;

public:
	Lexical()
	{
SEQAN_CHECKPOINT
	}

	template <typename TLeft, typename TRight>
	Lexical(TLeft const & left, TRight const & right)
	{
SEQAN_CHECKPOINT
		compare(*this, left, right);
	}

	Lexical(Lexical const & other):
		data_lcp(other.data_lcp),
		data_compare(other.data_compare)
	{
SEQAN_CHECKPOINT
	};

	Lexical & operator=(Lexical const & other)
	{
SEQAN_CHECKPOINT
		data_compare = other.data_compare;
		data_lcp = other.data_lcp;
		return *this;
	}

	~Lexical() {}
//____________________________________________________________________________

	enum
	{
		EQUAL = 1,
		LESS = 2,
		GREATER = 4,
		LEFT_IS_PREFIX = 8,
		RIGHT_IS_PREFIX = 16
	};
};


//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////
// Comparator: returns object that can compare objects of type T

/**
.Metafunction.Comparator:
..summary:Type of comparator object
..signature:Comparator<T>::Type
..param.T:Type for which the comparator type is to be determined.
..returns.param.Type:Comparator type
..remarks:Comparators are objects that can be used to compare other objects and store the
result of comparisons.
*/
template <typename T>
struct Comparator
{
	typedef Lexical<typename Size<T>::Type> Type;
};

//////////////////////////////////////////////////////////////////////////////
// Size

template <typename TSpec>
struct Size<Lexical<TSpec> >
{
	typedef TSpec Type;
};

template <typename TSpec>
struct Size<Lexical<TSpec> const>
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////
// Spec

template <typename TSpec>
struct Spec<Lexical<TSpec> >
{
	typedef TSpec Type;
};

template <typename TSpec>
struct Spec<Lexical<TSpec> const>
{
	typedef TSpec Type;
};


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// compare
//////////////////////////////////////////////////////////////////////////////

/**
.Function.compare:
..cat:Comparisons
..summary:Compares two objects.
..signature:compare(comparator, left, right)
..param.left:The first objects.
..param.right:The second objects that is compared to $left$.
..param.comparator:Object that stores the results.
...type:Class.Lexical
..see:Metafunction.Comparator
*/

template <typename TSpec, typename TLeft, typename TRight>
inline void
compare_(Lexical<TSpec> & lexical, 
		 TLeft & left, 
		 TRight & right)
{
SEQAN_CHECKPOINT
	typedef typename Value<TLeft>::Type TLeftValue;

	typename Iterator<TLeft, Standard>::Type left_it = begin(left, Standard());
	typename Size<TLeft>::Type left_length = length(left);
	typename Iterator<TRight, Standard>::Type right_it = begin(right, Standard());
	typename Size<TRight>::Type right_length = length(right);

	if (left_length == right_length) lexical.data_compare = Lexical<TSpec>::EQUAL;
	else if (left_length < right_length) lexical.data_compare = Lexical<TSpec>::LEFT_IS_PREFIX;
	else
	{
		lexical.data_compare = Lexical<TSpec>::RIGHT_IS_PREFIX;
		left_length = right_length;
	}

	lexical.data_lcp = 0;
	for (lexical.data_lcp = 0; lexical.data_lcp < left_length; ++lexical.data_lcp)
	{
		if (*left_it < *right_it)
		{
			lexical.data_compare = Lexical<TSpec>::LESS;
			break;
		}
		if (*left_it > *right_it)
		{
			lexical.data_compare = Lexical<TSpec>::GREATER;
			break;
		}
		++left_it;
		++right_it;
	}
}
//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TLeft, typename TRight>
inline void
compare(Lexical<TSpec> & lexical, 
		TLeft const & left, 
		TRight const & right)
{
	compare_(lexical, left, right);
}

//workaround for VC++ "const arrays" bug
template <typename TSpec, typename TLeftValue, typename TRight>
inline void
compare(Lexical<TSpec> & lexical, 
		TLeftValue const * left, 
		TRight const & right)
{
	compare_(lexical, left, right);
}
template <typename TSpec, typename TLeftValue, typename TRightValue>
inline void
compare(Lexical<TSpec> & lexical, 
		TLeftValue const * left, 
		TRightValue const * right)
{
	compare_(lexical, left, right);
}
template <typename TSpec, typename TLeft, typename TRightValue>
inline void
compare(Lexical<TSpec> & lexical, 
		TLeft const & left, 
		TRightValue const * right)
{
	compare_(lexical, left, right);
}

//////////////////////////////////////////////////////////////////////////////
// isEqual
//////////////////////////////////////////////////////////////////////////////

/**
.Function.isEqual:
..cat:Comparisons
..summary:Operator "==".
..signature:isEqual(left, right)
..signature:isEqual(comparator)
..param.left:The first parameter.
..param.right:The second parameter that is compared to $left$.
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ equals $right$, $false$ otherwise.
..see:Metafunction.Comparator
*/
template <typename TLeft, typename TRight >
inline bool
isEqual(TLeft const & left, 
		TRight const & right)
{
SEQAN_CHECKPOINT
	return left == right;
}

template <typename TSpec>
inline bool
isEqual(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
	return (_lex.data_compare & Lexical<TSpec>::EQUAL);
}

//////////////////////////////////////////////////////////////////////////////
// isNotEqual
//////////////////////////////////////////////////////////////////////////////

/**
.Function.isNotEqual:
..cat:Comparisons
..summary:Operator "!=".
..signature:isNotEqual(left, right)
..signature:isNotEqual(comparator)
..param.left:The first parameter.
..param.right:The second parameter that is compared to $left$.
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ is not equal to $right$, $false$ otherwise.
..see:Metafunction.Comparator
*/
template <typename TLeft, typename TRight >
inline bool
isNotEqual(TLeft const & left, 
		 TRight const & right)
{
SEQAN_CHECKPOINT
	return left != right;
}

template <typename TSpec>
inline bool
isNotEqual(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
	return !(_lex.data_compare & Lexical<TSpec>::EQUAL);
}

//////////////////////////////////////////////////////////////////////////////
// isLess
//////////////////////////////////////////////////////////////////////////////

/**
.Function.isLess:
..cat:Comparisons
..summary:Operator "<".
..signature:isLess(left, right [, prefix_order_tag])
..signature:isLess(comparator)
..param.left:The first parameter.
..param.right:The second parameter that is compared to $left$.
..param.prefix_order_tag:Tag that specify whether prefixes are less or greater. (optional)
...text:If omitted, the default tag is determined by @Metafunction.DefaultPrefixOrder@ for the type of $left$.
...see:Tag.Prefix Order
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ is less than $right$, $false$ otherwise.
..see:Metafunction.Comparator
..remarks:
...text:Sequences are compared in lexicographical order.
..see:Tag.Prefix Order
..see:Metafunction.DefaultPrefixOrder
*/
template <typename TLeft, typename TRight, typename TPrefixOrder >
inline bool
isLess(TLeft const & left, 
	   TRight const & right,
	   Tag<TPrefixOrder> const tag)
{
SEQAN_CHECKPOINT
	typename Comparator<TLeft>::Type _lex(left, right);
    return isLess(_lex, tag);
}
template <typename TLeft, typename TRight>
inline bool
isLess(TLeft const & left, 
	   TRight const & right)
{
SEQAN_CHECKPOINT
	return left < right;
}

template <typename TSpec>
inline bool
isLess(Lexical<TSpec> const & _lex,
	   TagPrefixLess)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::LESS | Lexical<TSpec>::LEFT_IS_PREFIX));
}
template <typename TSpec>
inline bool
isLess(Lexical<TSpec> const & _lex,
	   TagPrefixGreater)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::LESS | Lexical<TSpec>::RIGHT_IS_PREFIX));
}
template <typename TSpec>
inline bool
isLess(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
	return isLess(_lex, typename DefaultPrefixOrder< Lexical<TSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////
// isLessOrEqual
//////////////////////////////////////////////////////////////////////////////

/**
.Function.isLessOrEqual:
..cat:Comparisons
..summary:Operator "<=".
..signature:isLessOrEqual(left, right [, prefix_order_tag])
..signature:isLessOrEqual(comparator)
..param.left:The first parameter.
..param.right:The second parameter that is compared to $left$.
..param.prefix_order_tag:Tag that specify whether prefixes are less or greater. (optional)
...text:If omitted, the default tag is determined by @Metafunction.DefaultPrefixOrder@ for the type of $left$.
...see:Tag.Prefix Order
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ is less than or equal to $right$, $false$ otherwise.
..see:Metafunction.Comparator
..remarks:
...text:Sequences are compared in lexicographical order.
..see:Tag.Prefix Order
..see:Metafunction.DefaultPrefixOrder
*/

template <typename TLeft, typename TRight, typename TPrefixOrder >
inline bool
isLessOrEqual(TLeft const & left, 
		TRight const & right,
		Tag<TPrefixOrder> const tag)
{
SEQAN_CHECKPOINT
	typename Comparator<TLeft>::Type _lex(left, right);
    return isLessOrEqual(_lex, tag);
}
template <typename TLeft, typename TRight>
inline bool
isLessOrEqual(TLeft const & left, 
		TRight const & right)
{
SEQAN_CHECKPOINT
	return left <= right;
}

template <typename TSpec>
inline bool
isLessOrEqual(Lexical<TSpec> const & _lex,
		TagPrefixLess)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::LESS | Lexical<TSpec>::EQUAL | Lexical<TSpec>::LEFT_IS_PREFIX));
}
template <typename TSpec>
inline bool
isLessOrEqual(Lexical<TSpec> const & _lex,
		TagPrefixGreater)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::LESS | Lexical<TSpec>::EQUAL | Lexical<TSpec>::RIGHT_IS_PREFIX));
}
template <typename TSpec>
inline bool
isLessOrEqual(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
	return isLessOrEqual(_lex, typename DefaultPrefixOrder< Lexical<TSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////
// isGreater
//////////////////////////////////////////////////////////////////////////////

/**
.Function.isGreater:
..cat:Comparisons
..summary:Operator ">".
..signature:isGreater(left, right [, prefix_order_tag])
..signature:isGreater(comparator)
..param.left:The first parameter.
..param.right:The second parameter that is compared to $left$.
..param.prefix_order_tag:Tag that specify whether prefixes are less or greater. (optional)
...text:If omitted, the default tag is determined by @Metafunction.DefaultPrefixOrder@ for the type of $left$.
...see:Tag.Prefix Order
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ is greater than $right$, $false$ otherwise.
..see:Metafunction.Comparator
..remarks:
...text:Sequences are compared in lexicographical order.
..see:Tag.Prefix Order
..see:Metafunction.DefaultPrefixOrder
*/
template <typename TLeft, typename TRight, typename TPrefixOrder >
inline bool
isGreater(TLeft const & left, 
		TRight const & right,
		Tag<TPrefixOrder> const tag)
{
SEQAN_CHECKPOINT
	typename Comparator<TLeft>::Type _lex(left, right);
    return isGreater(_lex, tag);
}
template <typename TLeft, typename TRight>
inline bool
isGreater(TLeft const & left, 
		TRight const & right)
{
SEQAN_CHECKPOINT
	return left > right;
}

template <typename TSpec>
inline bool
isGreater(Lexical<TSpec> const & _lex,
		TagPrefixLess)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::GREATER | Lexical<TSpec>::RIGHT_IS_PREFIX));
}
template <typename TSpec>
inline bool
isGreater(Lexical<TSpec> const & _lex,
		TagPrefixGreater)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::GREATER | Lexical<TSpec>::LEFT_IS_PREFIX));
}
template <typename TSpec>
inline bool
isGreater(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
	return isGreater(_lex, typename DefaultPrefixOrder< Lexical<TSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////
// isGreaterOrEqual
//////////////////////////////////////////////////////////////////////////////

/**
.Function.isGreaterOrEqual:
..cat:Comparisons
..summary:Operator ">=".
..signature:isGreaterOrEqual(left, right [, prefix_order_tag])
..signature:isGreaterOrEqual(comparator)
..param.left:The first parameter.
..param.right:The second parameter that is compared to $left$.
..param.prefix_order_tag:Tag that specify whether prefixes are less or greater. (optional)
...text:If omitted, the default tag is determined by @Metafunction.DefaultPrefixOrder@ for the type of $left$.
...see:Tag.Prefix Order
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ is greater than or equal to $right$, $false$ otherwise.
..see:Metafunction.Comparator
..remarks:
...text:Sequences are compared in lexicographical order.
..see:Tag.Prefix Order
..see:Metafunction.DefaultPrefixOrder
*/

template <typename TLeft, typename TRight, typename TPrefixOrder >
inline bool
isGreaterOrEqual(TLeft const & left, 
		TRight const & right,
		Tag<TPrefixOrder> const tag)
{
SEQAN_CHECKPOINT
	typename Comparator<TLeft>::Type _lex(left, right);
    return isGreaterOrEqual(_lex, tag);
}
template <typename TLeft, typename TRight>
inline bool
isGreaterOrEqual(TLeft const & left, 
		TRight const & right)
{
SEQAN_CHECKPOINT
	return left >= right;
}

template <typename TSpec>
inline bool
isGreaterOrEqual(Lexical<TSpec> const & _lex,
		TagPrefixLess)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::GREATER | Lexical<TSpec>::EQUAL | Lexical<TSpec>::RIGHT_IS_PREFIX));
}
template <typename TSpec>
inline bool
isGreaterOrEqual(Lexical<TSpec> const & _lex,
		TagPrefixGreater)
{
SEQAN_CHECKPOINT
   return (_lex.data_compare & (Lexical<TSpec>::GREATER | Lexical<TSpec>::EQUAL | Lexical<TSpec>::LEFT_IS_PREFIX));
}
template <typename TSpec>
inline bool
isGreaterOrEqual(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
	return isGreaterOrEqual(_lex, typename DefaultPrefixOrder< Lexical<TSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////
// isPrefix
//////////////////////////////////////////////////////////////////////////////

/**
.Function.isPrefix:
..cat:Comparisons
..summary:Test whether a sequence is prefix of another sequence.
..signature:isPrefix(left, right)
..signature:isPrefix(comparator)
..param.left:The first sequence, the putative prefix.
..param.right:The second sequence.
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $left$ is a prefix of $right$, $false$ otherwise.
..see:Metafunction.Comparator
..remarks:By definition, the whole sequence is a prefix of itself too: $isPrefix("abc", "abc") == true$.
*/

template <typename TLeft, typename TRight >
inline bool
isPrefix(TLeft const & left, 
		TRight const & right)
{
SEQAN_CHECKPOINT
	typename Comparator<TLeft>::Type _lex(left, right);
    return isPrefix(_lex);
}
template <typename TSpec>
inline bool
isPrefix(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
    return (_lex.data_compare & (Lexical<TSpec>::LEFT_IS_PREFIX | Lexical<TSpec>::EQUAL));
}


//////////////////////////////////////////////////////////////////////////////
// hasPrefix
//////////////////////////////////////////////////////////////////////////////

/**
.Function.hasPrefix:
..cat:Comparisons
..summary:Test whether a sequence is prefix of another sequence.
..signature:hasPrefix(left, right)
..signature:hasPrefix(comparator)
..param.left:The first sequence.
..param.right:The second sequence, the putative prefix.
..param.comparator:A comparator.
...type:Class.Lexical
..returns:$true$ if $right$ is a prefix of $left$, $false$ otherwise.
..see:Metafunction.Comparator
..see:Function.isPrefix
..remarks:By definition, the whole sequence is a prefix of itself too: $hasPrefix("abc", "abc") == true$.
*/

template <typename TLeft, typename TRight >
inline bool
hasPrefix(TLeft const & left, 
		TRight const & right)
{
SEQAN_CHECKPOINT
	typename Comparator<TLeft>::Type _lex(left, right);
    return hasPrefix(_lex);
}
template <typename TSpec>
inline bool
hasPrefix(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
    return (_lex.data_compare & (Lexical<TSpec>::RIGHT_IS_PREFIX | Lexical<TSpec>::EQUAL));
}

//////////////////////////////////////////////////////////////////////////////
// lcpLength
//////////////////////////////////////////////////////////////////////////////
/**
.Function.lcpLength:
..summary:Length of longest common prefix.
..cat:Comparisons
..signature:lcpLength(left, right)
..signature:lcpLength(comparator)
..param.left:The first sequence.
..param.right:The second sequence that is compared to $left$.
..param.comparator:A comparator.
...type:Class.Lexical
..returns:The length of the longest common prefix of $left$ and $right$.
..see:Metafunction.Comparator
*/
template <typename TLeft, typename TRight >
inline typename Size<TLeft>::Type
lcpLength(TLeft const & left, TRight const & right)
{
SEQAN_CHECKPOINT
	typename Comparator<TLeft>::Type _lex(left, right);
    return lcpLength(_lex);
}

template <typename TSpec>
inline typename Size< Lexical<TSpec> >::Type
lcpLength(Lexical<TSpec> const & _lex)
{
SEQAN_CHECKPOINT
    return _lex.data_lcp;
}

//////////////////////////////////////////////////////////////////////////////
// lcpLength
//////////////////////////////////////////////////////////////////////////////
/**
.Function.ordValue:
..summary:Maps an alphabet 1-to-1 to the interval [0..ValueSize).
..cat:Alphabets
..signature:ordValue(value)
..param.value:Arbitrary character value.
...type:Class.SimpleType
..returns:An $unsigned int$ between 0 and @Metafunction.ValueSize@ of the type of value.
..note:This function first converts value to its unsigned value type and after that to an $unsigned int$.
You can't use $(unsigned int)c$ for a character $c$ as on some systems $char$ is signed and a $-1$ would be mapped to $0xffffffff$ instead of $0x000000ff$.
*/

template <typename TValue>
inline unsigned ordValue(TValue const &c) 
{
	return (typename _MakeUnsigned<TValue>::Type const &)c;
}

template <typename TValue, typename TSpec>
inline unsigned ordValue(SimpleType<TValue,TSpec> const &c) 
{
	return c;
}



//////////////////////////////////////////////////////////////////////////////


} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
