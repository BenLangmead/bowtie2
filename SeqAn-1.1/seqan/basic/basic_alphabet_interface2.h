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
  $Id: basic_alphabet_interface2.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_ALPHABET_INTERFACE2_H
#define SEQAN_HEADER_BASIC_ALPHABET_INTERFACE2_H

#include <new>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// gapValue, gapValueImpl
//////////////////////////////////////////////////////////////////////////////
/**
.Function.gapValueImpl:
..hidefromindex
..cat:Alphabets
..summary:Implements @Function.gapValue@.
..signature:gapValueImpl(value_pointer_tag)
..param.value_pointer_tag:A pointer that is used as a tag to specify the value type.
...remarks:The pointer needs not to point to a valid object, so it is possible to use a null pointer here.
..returns:A gap character.
..remarks.text:This function implements @Function.getValue@. 
It is recommended to use @Function.gapValue@ rather than $gapValueImpl$.
*/

template <typename T>
inline T const &
gapValueImpl(T *)
{
SEQAN_CHECKPOINT
	static T const _gap = T();
	return _gap;
}


/**
.Function.gapValue:
..cat:Alphabets
..cat:Alignments
..summary:Returns reference to a value that is used as gap character.
..signature:gapValue<TValue>()
..param.TValue:Value type.
..returns:A gap character.
..remarks.text:The function is implemented in @Function.gapValueImpl@. 
Do not specialize $gapValue$, specialize @Function.gapValueImpl@ instead!
..see:Function.gapValueImpl
*/

template <typename T>
inline T const &
gapValue()
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return gapValueImpl(_tag);
}


//////////////////////////////////////////////////////////////////////////////
// supremumValue, supremumValueImpl
//////////////////////////////////////////////////////////////////////////////

/**
.Function.supremumValueImpl:
..hidefromindex
..cat:Alphabets
..summary:Implements @Function.supremumValue@.
..signature:supremumValueImpl(value_pointer_tag)
..param.value_pointer_tag:A pointer that is used as a tag to specify the value type.
...remarks:The pointer needs not to point to a valid object, so it is possible to use a null pointer here.
..returns:A value $inf$ that holds: $inf >= i$ for all values $i$.
..remarks.text:This function implements @Function.supremumValue@. 
It is recommended to use @Function.supremumValue@ rather than $supremumValueImpl$.
*/

/*
template <typename T>
inline T const &
supremumValueImpl(T *)
{
	static T const _value = -1;
	return _value;
}
*/

/**
.Function.supremumValue:
..cat:Alphabets
..summary:Supremum for a given type.
..signature:supremumValue<T>()
..param.T:An ordered type.
..returns:A value $inf$ that holds: $inf >= i$ for all values $i$ of type $T$.
..remarks.text:The function is implemented in @Function.supremumValueImpl@. 
Do not specialize $supremumValue$, specialize @Function.supremumValueImpl@ instead!
..see:Function.supremumValueImpl
*/

template <typename T>
inline T const &
supremumValue()
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return supremumValueImpl(_tag);
}

//////////////////////////////////////////////////////////////////////////////
// infimumValue, infimumValueImpl
//////////////////////////////////////////////////////////////////////////////

/**
.Function.infimumValueImpl:
..hidefromindex
..cat:Alphabets
..summary:Implements @Function.infimumValue@.
..signature:infimumValueImpl(value_pointer_tag)
..param.value_pointer_tag:A pointer that is used as a tag to specify the value type.
...remarks:The pointer needs not to point to a valid object, so it is possible to use a null pointer here.
..returns:A value $inf$ that holds: $inf <= i$ for all values $i$.
..remarks.text:This function implements @Function.infimumValue@. 
It is recommended to use @Function.infimumValue@ rather than $infimumValueImpl$.
*/

/*
template <typename T>
inline T const &
infimumValueImpl(T *)
{
	static T const _value = -1;
	return _value;
}
*/

/**
.Function.infimumValue:
..cat:Alphabets
..summary:Infimum for a given type.
..signature:infimumValue<T>()
..param.T:An ordered type.
..returns:A value $inf$ that holds: $inf <= i$ for all values $i$ of type $T$.
..remarks.text:The function is implemented in @Function.infimumValueImpl@. 
Do not specialize $infimumValue$, specialize @Function.infimumValueImpl@ instead!
..see:Function.infimumValueImpl
..see:Function.supremumValue
*/

template <typename T>
inline T const &
infimumValue()
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return infimumValueImpl(_tag);
}

//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
