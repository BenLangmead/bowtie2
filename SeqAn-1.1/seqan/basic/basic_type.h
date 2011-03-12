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
  $Id: basic_type.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_TYPE_H
#define SEQAN_HEADER_BASIC_TYPE_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Value:
..summary:Type of the items in the container. 
..signature:Value<T>::Type
..param.T:Type for which the value type is determined.
..returns.param.Type:Value type of $T$.
..remarks.text:The value type of a container $T$ is the type of the elements in $T$.
    For example, the value type of a sequence of $int$ is $int$.
..example.code:Value<String<char> >::Type c; //c has type char
*/

template <typename T, const int i = 0>
struct Value
{
	typedef T Type;
};
/*
template <typename T>
struct Value<T const>
{
	typedef T Type;
};
*/
//____________________________________________________________________________

/**
.Metafunction.GetValue:
..summary:Type for reading values. 
..signature:GetValue<T>::Type
..param.T:Type of container that holds a value.
..returns.param.Type:GetValue type of $T$.
..remarks.text:Depending on $T$, the $GetValue$-type can either be $Value<T>::Type &$ or $Value<T>::Type$.
..text:$GetValue$ is the return type of @Function.getValue@ that allows a (read-only) access to objects.
Do not confuse it with @Function.value@ that returns a @Metafunction.Reference.reference@ to the value.
..see:Metafunction.Value
..see:Function.getValue
*/
template <typename T>
struct GetValue
{
	typedef typename Value<T>::Type const & Type;
};
template <typename T>
struct GetValue<T const>:
	public GetValue<T>
{
};

//____________________________________________________________________________

/**
.Metafunction.Reference:
..summary:Reference type. 
..signature:Reference<T>::Type
..param.T:A Type.
..returns.param.Type:Either $T &$ or a proxy object @Class.Proxy@ for $T$.
..see:Metafunction.Value
..see:Metafunction.GetValue
*/
template <typename T>
struct Reference
{
	typedef typename Value<T>::Type & Type;
};
template <typename T>
struct Reference<T const>
{
	typedef typename Value<T>::Type const & Type;
};

//____________________________________________________________________________


/**
.Metafunction.Size:
..summary:Type of an object that is suitable to hold size information.
..signature:Size<T>::Type
..param.T:Type for which the size type is determined.
..returns.param.Type:Size type of $T$.
..remarks.text:In most cases this type is $size_t$.
*/
template <typename T>
struct Size
{
	typedef size_t Type;
};
template <typename T>
struct Size<T const>:
	Size<T>
{
};

//____________________________________________________________________________


/**
.Metafunction.Difference:
..summary:Type of an object that stores the difference between two iterators.
..signature:Difference<T>::Type
..param.T:Type for which the difference type is determined.
...type:Class.Iter
..returns.param.Type:Difference type of $T$.
..remarks.text:In most cases this type is $ptrdiff_t$.
..see:Metafunction.Size
*/
template <typename T>
struct Difference
{
	typedef ptrdiff_t Type;
};
template <typename T>
struct Difference<T const>:
	Difference<T>
{
};

//____________________________________________________________________________


/**
.Metafunction.Position:
..summary:Type of an object that represents a position in a container.
..signature:Position<T>::Type
..param.T:Type for which the position type is determined.
...type:Class.Iter
...type:Class.String
..returns.param.Type:Position type of $T$.
..see:Metafunction.Iterator
*/
template <typename T>
struct Position
{
	typedef typename Size<T>::Type Type;
};
template <typename T>
struct Position<T const>:
	Position<T>
{
};

//____________________________________________________________________________

/**
.Metafunction.Host:
..summary:Type of the object a given object depends on.
..signature:Host<T>::Type
..param.T:Type for which the host type is determined.
..returns.param.Type:Host type of $T$.
*/
template <typename T>
struct Host
{
	typedef T Type;
};

//____________________________________________________________________________

/**
.Metafunction.Spec:
..summary:The spec of a class. 
..signature:Spec<T>::Type
..param.T:Type for which the spec is determined.
..returns.param.Type:Spec of $T$.
..remarks:The spec of a SeqAn type is the class that is used in template subclassing 
 to specify the specialization. 
 For example, the spec of $String<char, Alloc<> >$ is $Alloc<>$.
*/

// default case
template <typename T>
struct Spec {
	typedef void Type;
};


// one argument case
template <template <typename> class T, typename TSpec>
struct Spec< T<TSpec> > {
	typedef TSpec Type;
};

template <typename T>
struct Spec<T const>:
	public Spec<T> {};

//____________________________________________________________________________

/**
.Metafunction.DeepestSpec:
..summary:The deepest spec of a class with nested template arguments.
..signature:DeepestSpec<T>::Type
..param.T:Type for which the deepest spec is determined.
..returns.param.Type:Deepest spec of $T$.
..remarks:The spec of a SeqAn type is the innermost class that is used in nested subclassing.
 For example, the deepest spec of $Iter<..., VSTree<BottomUp<MUMs> > >$ is $MUMs$.
*/

// default case
template <typename T>
struct DeepestSpec {
	typedef T Type;
};

// recursion for 1 argument
template <
	template <typename> class T, 
	typename T1 >
struct DeepestSpec< T<T1> > {
	typedef typename 
		IF<
			TYPECMP<T1, void>::VALUE,										// is T1 void?
			T<T1>,															// yes, end of recursion
			typename DeepestSpec< typename Spec< T<T1> >::Type >::Type		// no,  recurse
		>::Type Type;
};

// recursion for 2 arguments
template <
	template <typename, typename> class T, 
	typename T1, typename T2 >
struct DeepestSpec< T<T1,T2> >:
	DeepestSpec< typename Spec< T<T1,T2> >::Type > {};

// recursion for 3 arguments
template <
	template <typename, typename, typename> class T, 
	typename T1, typename T2, typename T3 >
struct DeepestSpec< T<T1,T2,T3> >:
	DeepestSpec< typename Spec< T<T1,T2,T3> >::Type > {};

// recursion for 4 arguments
template <
	template <typename, typename, typename, typename> class T, 
	typename T1, typename T2, typename T3, typename T4 >
struct DeepestSpec< T<T1,T2,T3,T4> >:
	DeepestSpec< typename Spec< T<T1,T2,T3,T4> >::Type > {};

// recursion for 5 arguments
template <
	template <typename, typename, typename, typename, typename> class T, 
	typename T1, typename T2, typename T3, typename T4, typename T5 >
struct DeepestSpec< T<T1,T2,T3,T4,T5> >:
	DeepestSpec< typename Spec< T<T1,T2,T3,T4,T5> >::Type > {};

template <typename T>
struct DeepestSpec<T const>:
	public DeepestSpec<T> {};

//____________________________________________________________________________

/**
.Metafunction.Cargo:
..summary:Type of additional data stored in an object. 
..signature:Cargo<T>::Type
..param.T:Type for which the cargo tyoe is determined.
..returns.param.Type:Cargo of $T$.
..remarks:The definition of Cargo allows the addition of user specific data to existing data structures.
*/

template <typename T>
struct Cargo {
	typedef Nothing Type;
};
template <typename T>
struct Cargo<T const> {
	typedef typename Cargo<T>::Type const Type;
};

//____________________________________________________________________________

/**
.Metafunction.VertexDescriptor:
..summary:Type of an object that represents a vertex descriptor.
..signature:VertexDescriptor<T>::Type
..param.T:Type T must be a graph. All graphs currently use ids as vertex descriptors.
..returns.param.Type:VertexDescriptor type.
..remarks.text:The vertex descriptor is a unique handle to a vertex in a graph.
It is used in various graph functions, e.g., to add edges, to create OutEdge Iterators or to remove a vertex.
It is also used to attach properties to vertices.
..example.code:VertexDescriptor<Graph<> >::Type vD; //vD is a vertex descriptor
*/

template <typename T>
struct VertexDescriptor {
	typedef void* Type;
};
template <typename T>
struct VertexDescriptor<T const>:
	public VertexDescriptor<T> {};


//____________________________________________________________________________

	
/**
.Metafunction.Id:
..summary:Type of an object that represents an id.
..signature:Id<T>::Type
..param.T:Type for which a suitable id type is determined.
..returns.param.Type:Id type.
..remarks.text:The id type of a container is the type that is used to uniquely identify its elements.
In most cases this type is unsigned int.
..example.code:Id<Graph<> >::Type id; //id has type unsigned int
*/
template<typename T>
struct Id {
	typedef unsigned int Type;
};

//____________________________________________________________________________

template<typename T>
struct Id<T const> {
	typedef unsigned int Type;
};

//____________________________________________________________________________

/**
.Metafunction.Key:
..summary:Key type of a key to cargo mapping.
..signature:Key<T>::Type
..param.T:Type for which a key type is determined.
..returns.param.Type:Key type.
...default:The type $T$ itself.
*/
template< typename T >
struct Key
{
	typedef T Type;
};

template <typename T>
struct Key<T const>:
	Key<T> {};

//____________________________________________________________________________

/*VERALTET
.Metafunction.Object:
..summary:Object type of a key to object mapping.
..signature:Object<T>::Type
..param.T:Type for which a object type is determined.
..returns.param.Type:Object type.
*/

template<typename T>
struct Object; 

template <typename T>
struct Object<T const>:
	Object<T> {};


//____________________________________________________________________________

/**
.Metafunction.Source
*/

template < typename TSpec = void >
struct Source
{
	typedef TSpec Type;
};

template <typename T>
struct Source<T const>:
	Source<T>
{
};

//____________________________________________________________________________

/**
.Internal._Parameter:
..cat:Metafunctions
..summary:Type for function parameters and return values.
..signature:_Parameter<T>::Type
..param.T:A type.
..returns.param.Type:The parameter type for arguments of type $T$.
...text:If $T$ is a pointer or array type, then $_Parameter<T>::Type$ is $T$, 
otherwise $_Parameter<T>::Type$ is $T &$.
*/
template <typename T>
struct _Parameter
{
	typedef T & Type;
};

template <typename T>
struct _Parameter<T *>
{
	typedef T * Type;
};
template <typename T, size_t I>
struct _Parameter<T [I]>
{
	typedef T * Type;
};


/**
.Internal._toParameter:
..cat:Functions
..summary:Transforms pointers to parameter types.
..signature:_toParameter<T>(pointer)
..param.pointer:A pointer.
..param.T:A Type.
...text:$object$ is transformed into the parameter type of $T$ that is given by @Internal._Parameter@.
...note:This type must be explicitely specified.
..returns:To $TParameter$ transformed $object$.
..see:Internal._Parameter
*/
template <typename T>
typename _Parameter<T>::Type
_toParameter(T * _object)
{
SEQAN_CHECKPOINT
	return * _object;
}
template <typename T>
typename _Parameter<T>::Type
_toParameter(T _object)
{
SEQAN_CHECKPOINT
	return _object;
}

//____________________________________________________________________________

/**
.Internal._ConstParameter:
..cat:Metafunctions
..summary:Type for constant function parameters and return values.
..signature:_ConstParameter<T>::Type
..param.T:A type.
..returns.param.Type:The const parameter type for arguments of type $T$.
...text:If $T$ is a pointer or array type, then $_Parameter<T>::Type$ is a pointer to a const array, 
otherwise $_Parameter<T>::Type$ is $T const &$.
..see:Internal._Parameter
*/
template <typename T>
struct _ConstParameter
{
	typedef T const & Type;
};
template <typename T>
struct _ConstParameter<T const>:
	public _ConstParameter<T> {};

template <typename T>
struct _ConstParameter<T *>
{
	typedef T const * Type;
};
template <typename T>
struct _ConstParameter<T const *>
{
	typedef T const * Type;
};

template <typename T, size_t I>
struct _ConstParameter<T [I]>
{
	typedef T const * Type;
};
template <typename T, size_t I>
struct _ConstParameter<T const [I]>
{
	typedef T const * Type;
};

//____________________________________________________________________________

/**
.Internal._Pointer:
..cat:Metafunctions
..summary:The associated pointer type.
..signature:_Pointer<T>::Type
..param.T:A type.
..returns.param.Type:A pointer type for $T$.
...text:if $T$ is already a pointer type, then $_Pointer<T>::Type$ is $T$,
otherwise $_Pointer<T>::Type$ is $T *$.
..see:Internal._Parameter
..see:Internal._toParameter
*/
template <typename T>
struct _Pointer
{
	typedef T * Type;
};

template <typename T>
struct _Pointer<T *>
{
	typedef T * Type;
};
template <typename T>
struct _Pointer<T * const>
{
	typedef T * const Type;
};

template <typename T, size_t I>
struct _Pointer<T [I]>
{
	typedef T * Type;
};

/**
.Internal._toPointer:
..cat:Functions
..summary:Transforms types into pointers.
..signature:_toPointer(object)
..param.object:An object.
..returns:$object$, transformed to a pointer. 
...text:The type of the returned pointer is given by @Internal._Pointer@.
..see:Internal._Pointer
*/
template <typename T>
typename _Pointer<T>::Type
_toPointer(T & _object)
{
SEQAN_CHECKPOINT
	return & _object;
}
template <typename T>
typename _Pointer<T const>::Type
_toPointer(T const & _object)
{
SEQAN_CHECKPOINT
	return & _object;
}

template <typename T>
typename _Pointer<T *>::Type
_toPointer(T * _object)
{
SEQAN_CHECKPOINT
	return _object;
}

//____________________________________________________________________________


/**
.Metafunction.LENGTH:
..summary:Number of elements in a fixed-size container.
..signature:LENGTH<T>::Type
..param.T:Type for which the number of elements is determined.
..returns.param.VALUE:Number of elements.
..remarks.text:The default return value is 1 for dynamic-size containers.
*/
template <typename T>
struct LENGTH
{
	enum { VALUE = 1 };
};
template <typename T>
struct LENGTH<T const>:
	LENGTH<T>
{
};

/**
.Metafunction.WEIGHT:
..summary:Number of relevant positions in a shape.
..signature:WEIGHT<T>::Type
..param.T:Shape type for which the number of relevant positions is determined.
...type:Class.Shape
..returns.param.VALUE:Number of relevant positions.
..remarks.text:The default return value is the result of the @Metafunction.LENGTH@ function.
For gapped shapes this is the number of '1's.
*/
template <typename T>
struct WEIGHT:
	LENGTH<T>
{
};
template <typename T>
struct WEIGHT<T const>:
	WEIGHT<T>
{
};

//////////////////////////////////////////////////////////////////////////////

//Iterator: see basic_iterator.h

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
