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
  $Id: basic_iterator.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_ITERATOR_H
#define SEQAN_HEADER_BASIC_ITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// TAGS
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Iterator Spec:
..summary:Specifies the kind of an iterator.
..tag.Rooted:Rooted iterator. 
...remarks
....text:This iterator implements some more advanced functions like
@Function.container@ and @Function.position@.
....concept:Concept.Rooted Iterator
..tag.Standard:Standard conform iterator. 
...remarks
....text:Note that standard iterators need not to implement all functions
that are available for rooted iterators.
....concept:Concept.Iterator
..remarks.text:The default iterator spec is given by @Metafunction.DefaultIteratorSpec@.
..see:Metafunction.DefaultIteratorSpec
..see:Concept.Iterator
*/

struct TagRooted_;
typedef Tag<TagRooted_> const Rooted;

struct TagStandard_;
typedef Tag<TagStandard_> const Standard;


//////////////////////////////////////////////////////////////////////////////
// METAFUNCTIONS
//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.DefaultIteratorSpec:
..hidefromindex
..summary:Specifies default kind of iterator.
..signature:DefaultIteratorSpec<T>::Type
..param.T:Container type for which the default iterator spec is determined.
...concept:Concept.Container
..returns.param.Type:Iterator spec of $T$.
..see:Metafunction.Iterator
*/

template <typename T>
struct DefaultIteratorSpec
{
	typedef Standard Type;
};

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.DefaultGetIteratorSpec:
..hidefromindex
..summary:Specifies default kind of iterator returned by functions.
..signature:DefaultGetIteratorSpec<T>::Type
..param.T:Container type for which the spec is determined.
...concept:Concept.Container
..returns.param.Type:Iterator spec of $T$.
..remarks:This metafunction returns the iterator spec of iterators that are returned by functions like 
@Function.begin@, @Function.end@, or @Function.iter@.
..see:Metafunction.Iterator
..see:Metafunction.DefaultIteratorSpec
*/

template <typename T>
struct DefaultGetIteratorSpec
{
	typedef Rooted Type;
};

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Iterator:
..summary:Type of iterator objects that are used to traverse the container.
..signature:Iterator<T, TSpec>::Type
..param.T:Type for which the iterator type is determined.
...concept:Concept.Container
...type:Class.Iter
..param.TSpec:Specifies an @Tag.Iterator Spec.iterator spec@.
...default:The default iterator spec is given by @Metafunction.DefaultIteratorSpec@.
..returns.param.Type:Iterator type of $T$.
..remarks.text:Iterators behave like pointers in some respects. 
 For example, you can use $*it$ to access the value object the iterator $it$ points to.
 But note that $Iterator<T>::Type$ can differ from $T *$, depending on $T$.
..see:Metafunction.Position
*/

//____________________________________________________________________________

template <typename T, typename TSpec>
struct Iterator_Default_Imp;

//Iterator_Default_Imp<T, Standard> is implemented in basic_iterator_simple.h
//Iterator_Default_Imp<T, Rooted> is implemented in basic_iterator_adaptor.h 

//____________________________________________________________________________

template <typename T, typename TSpec = typename DefaultIteratorSpec<T>::Type>
struct Iterator:
	Iterator_Default_Imp<T, TSpec>
{
};


//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Container:
..summary:Type of the container given an iterator.
..signature:Container<T>::Type
..param.T:Iterator type.
...type:Class.Iter
...concept:Concept.Iterator
..returns.param.Type:The container type to $T$.
*/

template <typename T>
struct Container
{
	typedef T Type;
};


//////////////////////////////////////////////////////////////////////////////
// GENERAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// value
//////////////////////////////////////////////////////////////////////////////

/**
.Function.value:
..signature:Reference value(object)
..param.object:An object that holds a value or an iterator that points to a value.
...type:Class.Iter
...concept:Concept.Iterator
*/

template <typename T>
inline typename Reference<T>::Type
value(T & me)
{
SEQAN_CHECKPOINT
	return *me;
} 
template <typename T>
inline typename Reference<T const>::Type
value(T const & me)
{
SEQAN_CHECKPOINT
	return *me;
} 


template <typename T>
inline T &
value(T * me)
{
SEQAN_CHECKPOINT
	return *me;
} 

//////////////////////////////////////////////////////////////////////////////
// getValue
//////////////////////////////////////////////////////////////////////////////

//unary getValue
/**
.Function.getValue:
..cat:Iteration
..signature:GetValue getValue(object)
..param.object:An object that holds a value or points to a value.
...type:Class.Iter
...concept:Concept.Iterator
..see:Metafunction.GetValue
*/

template <typename T>
inline typename GetValue<T>::Type
getValue(T & me)
{
SEQAN_CHECKPOINT
	return value(me);
} 
template <typename T>
inline typename GetValue<T const>::Type
getValue(T const & me)
{
SEQAN_CHECKPOINT
	return value(me);
} 

template <typename T>
inline T &
getValue(T * me)
{
SEQAN_CHECKPOINT
	return value(me);
} 

//////////////////////////////////////////////////////////////////////////////
//toGetValue
//////////////////////////////////////////////////////////////////////////////
//Nimmt eine Reference und macht daraus einen GetValue
//???TODO toGetValue()

//////////////////////////////////////////////////////////////////////////////
// assignValue
//////////////////////////////////////////////////////////////////////////////
/**
.Function.assignValue:
..cat:Iteration
..summary:Assigns value to item.
..signature:assignValue(object, value)
..param.object:An object that holds a value or points to a value.
...type:Class.Iter
...concept:Concept.Iterator
..param.value:A value that is assigned to the item $object$ holds or points to.
..remarks.text:This function is similar to @Function.assign@.
The difference is, that $assignValue$ just changes a value stored in $object$ or the value $object$ points to, 
while @Function.assign@ changes the whole object.
..see:Function.assign
*/

template <typename T, typename TValue>
inline void
assignValue(T & me,
			TValue const & _value)
{
SEQAN_CHECKPOINT
	assign(value(me), _value);
} 

//const version for iterators as targets
template <typename T, typename TValue>
inline void
assignValue(T const & me,
			TValue const & _value)
{
SEQAN_CHECKPOINT
	assign(value(me), _value);
} 

//////////////////////////////////////////////////////////////////////////////
// moveValue
//////////////////////////////////////////////////////////////////////////////

/**
.Function.moveValue:
..cat:Iteration
..summary:Assigns value to item.
..signature:moveValue(object, value)
..param.object:An object that holds a value or points to a value.
...type:Class.Iter
...concept:Concept.Iterator
..param.value:A value that is handed over to the item $object$ holds or points to.
..remarks.text:This function is similar to @Function.move@.
The difference is, that $moveValue$ just changes a value stored in $object$ or the value $object$ points to, 
while @Function.move@ changes the whole object.
..see:Function.move
..see:Function.assignValue
*/

template <typename T, typename TValue>
inline void
moveValue(T & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	move(value(me), _value);
}
//const version for iterators as targets
template <typename T, typename TValue>
inline void
moveValue(T const & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	move(value(me), _value);
} 

//////////////////////////////////////////////////////////////////////////////

/**
.Function.container:
..cat:Iteration
..summary:Container of an iterator.
..signature:Container container(iterator)
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.Rooted Iterator
..returns:The container that $iterator$ traverses.
...metafunction:Metafunction.Container
*/

template <typename T>
inline typename Container<T>::Type 
container(T me)
{
SEQAN_CHECKPOINT
	return me;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.position:
..summary:Position of an iterator.
..cat:Iteration
..signature:Position position(iterator [, container])
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.Iterator
..param.container:A container.
...concept:Concept.Container
...remarks:If $iterator$ implements @Concept.Rooted Iterator@, then $container$ is optional.
...remarks:If $container$ is specified, $iterator$ must be a container of $container$.
..returns:The position of the value in the container $iterator$ points to.
...metafunction:Metafunction.Position
*/

template <typename T>
inline typename Position<T>::Type 
position(T * me)
{
SEQAN_CHECKPOINT
	return 0;
}

template <typename TContainer, typename TIterator>
inline typename Position<TContainer>::Type 
position(TIterator const & it,
		 TContainer const & me)
{
SEQAN_CHECKPOINT
	return it - begin(me, Standard());
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// atBegin
//////////////////////////////////////////////////////////////////////////////

/**
.Function.atBegin:
..cat:Iteration
..summary:Determines whether an iterator is at the beginning position.
..signature:bool atBegin(iterator [, container])
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.Iterator
..param.container:Container of $iterator$. (optional)
...remarks.text:If $iterator$ implements @Concept.Rooted Iterator@ then $container$ is optional otherwise $container$ is required.
..returns:$true$ if $iterator$ points to the fist item of the container, otherwise $false$.
..see:Function.begin
*/

//TODO???: Was, wenn der Container leer ist?

template <typename T, typename TContainer>
inline bool
atBegin(T const & it, TContainer const & cont)
{
SEQAN_CHECKPOINT
	return it == begin(cont, Standard());	
}

//____________________________________________________________________________

template <typename T>
inline bool
atBegin(T const & it)
{
SEQAN_CHECKPOINT
	return atBegin(it, container(it));	
}



//////////////////////////////////////////////////////////////////////////////
// atEnd
//////////////////////////////////////////////////////////////////////////////
/**
.Function.atEnd:
..cat:Iteration
..summary:Determines whether an iterator is at the end position. 
..signature:bool atEnd(iterator [, container])
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.Iterator
..param.container:Container of $iterator$.
...remarks.text:If $iterator$ implements @Concept.Rooted Iterator@ then $container$ is optional.
....text:$container$ is also optional for iterators to @Adaption.char array.char arrays@.
....text:Otherwise, $container$ is required.
..returns:$true$ if $iterator$ points behind the last item of the container, otherwise $false$.
..see:Function.atBegin
..see:Function.end
*/

template <typename T, typename TContainer>
inline bool
atEnd(T & it, 
	  TContainer const & cont)
{
SEQAN_CHECKPOINT
	return it == end(cont, Standard());	
}
template <typename T, typename TContainer>
inline bool
atEnd(T const & it, 
	  TContainer const & cont)
{
SEQAN_CHECKPOINT
	return it == end(cont, Standard());	
}
//____________________________________________________________________________

template <typename T>
inline bool
atEnd(T & it)
{
SEQAN_CHECKPOINT
	return atEnd(it, container(it));	
}
template <typename T>
inline bool
atEnd(T const & it)
{
SEQAN_CHECKPOINT
	return atEnd(it, container(it));	
}

//////////////////////////////////////////////////////////////////////////////
// goBegin
//////////////////////////////////////////////////////////////////////////////
/**
.Function.goBegin:
..cat:Iteration
..summary:Iterates to the first position of a container. 
..signature:goBegin(iterator [, container])
..param.iterator:Object that iterates through $container$.
...type:Class.Iter
...concept:Concept.Iterator
...text:$iterator$ is set to the position of the first item in $container$.
..param.container:Container of $iterator$.
...remarks.text:If $iterator$ implements @Concept.Rooted Iterator@ then $container$ is optional,
otherwise $container$ is required.
..remarks:This function is equivalent to $iterator = begin(container)$.
..see:Function.begin
..see:Function.atBegin
..see:Function.goEnd
*/
template <typename TIterator, typename TContainer>
inline void
goBegin(TIterator & it,
		TContainer & container)
{
SEQAN_CHECKPOINT
	it = begin(container);
}
/*
template <typename TIterator, typename TContainer>
inline void
goBegin(TIterator & it,
		TContainer const & container)
{
SEQAN_CHECKPOINT
	it = begin(container);
}
*/

template <typename TIterator>
inline void
goBegin(TIterator & it)
{
SEQAN_CHECKPOINT
	goBegin(it, container(it));
}

//////////////////////////////////////////////////////////////////////////////
// goEnd
//////////////////////////////////////////////////////////////////////////////
/**
.Function.goEnd:
..cat:Iteration
..summary:Iterates to the last position of a container. 
..signature:goEnd(iterator [, container])
..param.iterator:Object that iterates through $container$.
...type:Class.Iter
...concept:Concept.Iterator
...text:$iterator$ is set to the position behin the last item in $container$.
..param.container:Container of $iterator$.
...remarks.text:If $iterator$ implements @Concept.Rooted Iterator@ then $container$ is optional,
otherwise $container$ is required.
..remarks:This function is equivalent to $iterator = end(container)$.
..see:Function.end
..see:Function.atEnd
..see:Function.goBegin
..see:Function.goEnd
*/
template <typename TIterator, typename TContainer>
inline void
goEnd(TIterator & it,
	  TContainer & container)
{
SEQAN_CHECKPOINT
	it = end(container);
}
template <typename TIterator, typename TContainer>
inline void
goEnd(TIterator & it,
	  TContainer const & container)
{
SEQAN_CHECKPOINT
	it = end(container);
}

template <typename TIterator>
inline void
goEnd(TIterator & it)
{
SEQAN_CHECKPOINT
	goEnd(it, container(it));
}

//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////
/**
.Function.goNext:
..cat:Iteration
..summary:Iterates to next position. 
..signature:goNext(iterator)
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.Iterator
...text:$iterator$ is set to the next position of an iteration through its container.
..remarks:This function is equivalent to $++iterator$.
..see:Function.goBegin
..see:Function.goEnd
*/
template <typename TIterator>
inline void
goNext(TIterator & it)
{
SEQAN_CHECKPOINT
	++it;
}

//////////////////////////////////////////////////////////////////////////////
// goFurther
//////////////////////////////////////////////////////////////////////////////
/**
.Function.goFurther:
..cat:Iteration
..summary:Iterates some steps further. 
..signature:goFurther(iterator, steps)
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.Iterator
...text:$iterator$ is set $steps$ positions further in the iteration through the container.
..param.steps:Number of steps $iterator$ should be moved further.
...remarks:If $iterator$ supports bidirectional iteration, $steps$ could also be negativ.
..remarks:This function is equivalent to $iterator += steps$ for random access iterators.
..see:Function.goNext
..see:Function.goPrevious
*/

template <typename TIterator, typename TDiff> inline
void goFurther(TIterator & it, TDiff steps)
{	// return distance type from arbitrary argument
    it += steps;
}

//////////////////////////////////////////////////////////////////////////////
// goPrevious
//////////////////////////////////////////////////////////////////////////////
/**
.Function.goPrevious:
..cat:Iteration
..summary:Iterates to pevious position. 
..signature:goPrevious(iterator)
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.Iterator
...text:$iterator$ is set to the pevious position of an iteration through its container.
..remarks:This function is equivalent to $--iterator$.
..see:Function.goBegin
..see:Function.goEnd
..see:Function.goNext
*/
template <typename TIterator>
inline void
goPrevious(TIterator & it)
{
SEQAN_CHECKPOINT
	--it;
}

//////////////////////////////////////////////////////////////////////////////
// difference
//////////////////////////////////////////////////////////////////////////////
/**
.Function.difference:
..cat:Iteration
..summary:The difference between two iterators. 
..signature:difference(begin, end)
..param.begin:Iterator to the first position of a range.
...type:Class.Iter
...Concept.Iterator
..param.end:Iterator behind the last position of a range.
...type:Class.Iter
...Concept.Iterator
..returns:Length of the range between $begin$ and $end$.
..remarks:This function is equivalent to $begin - end$.
...text:Usually, $begin$ and $end$ have the same type.
..see:Function.begin
..see:Function.end
..see:Function.length
*/

template <typename TIterator> inline
typename Difference<TIterator>::Type difference(
	TIterator const & begin, 
	TIterator const & end)
{	// return distance type from arbitrary argument
SEQAN_CHECKPOINT
    return end - begin;
}

//////////////////////////////////////////////////////////////////////////////
// goNil
//////////////////////////////////////////////////////////////////////////////


/**
.Function.goNil:
..cat:Iteration
..summary:Moves iterator to nil position.
..signature:goNil(iterator)
..param.iterator:The iterator that will be moved.
...type:Class.String
..remarks:$iterator$ is set to an invalid position, e.g. $NULL$ for pointer types.
..see:Function.clear
*/

template <typename TIterator>
inline void
goNil(TIterator & me)
{
SEQAN_CHECKPOINT
	me = TIterator();
}

template <typename TIterator>
inline void
goNil(TIterator * & me)
{
SEQAN_CHECKPOINT
	me = 0;
}

//////////////////////////////////////////////////////////////////////////////
// atNil
//////////////////////////////////////////////////////////////////////////////


/**
.Function.atNil:
..cat:Iteration
..summary:Tests whether iterator is at nil position.
..signature:bool atNil(iterator)
..param.iterator:An iterator.
...type:Class.String
..returns:$true$ if $iterator$ points to an ivalid position, e.g. $iterator$ is a $NULL$ pointer.
$false$ otherwise.
..see:Function.goNil
*/

template <typename TIterator>
inline bool
atNil(TIterator & me)
{
SEQAN_CHECKPOINT
	return me == TIterator();
}


template <typename TIterator>
inline bool
atNil(TIterator * me)
{
SEQAN_CHECKPOINT
	return me == 0;
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
