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
  $Id: basic_holder.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_HOLDER_H
#define SEQAN_HEADER_BASIC_HOLDER_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// addRef
//////////////////////////////////////////////////////////////////////////////
/**
.Function.addRef:
..summary:Called when dependency is added. 
..cat:Dependent Objects
..signature:addRef(host)
..param.host:The host object.
..remarks.text:A call of this function denotes that a client object is about to become
dependent on $host$.
..remarks.text:The default behavior is: Do nothing.
..see:Class.Holder
*/

template <typename T>
inline void 
addRef(T & /*me*/)
{// general: do nothing
SEQAN_CHECKPOINT
}
template <typename T>
inline void 
addRef(T const & /*me*/)
{// general: do nothing
SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////
// releaseRef
//////////////////////////////////////////////////////////////////////////////
/**
.Function.releaseRef:
..summary:Called when dependency is released. 
..cat:Dependent Objects
..signature:releaseRef(host)
..param.host:The host object.
..remarks.text:A call of this function denotes that a former dependent client object 
ceases to be dependent on $host$.
..remarks.text:The default behavior is: Do nothing.
..see:Class.Holder
..see:Function.addRef
*/
template <typename T>
inline void 
releaseRef(T & /*me*/)
{// general: do nothing
SEQAN_CHECKPOINT
}
template <typename T>
inline void 
releaseRef(T const & /*me*/)
{// general: do nothing
SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////
// Tags

struct Simple;
struct Tristate;


//////////////////////////////////////////////////////////////////////////////
// Holder
//////////////////////////////////////////////////////////////////////////////

/**
.Class.Holder:
..cat:Basic
..summary:Manages relationship to another object.
..signature:Holder<TValue, TSpec>
..param.TValue:Type of the managed object.
...metafunction:Metafunction.Value
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Tristate$
..remarks.text:The main purpose of this class is to facilitate the handling of
member objects. If we want class $A$ to be dependent on or the owner of another object of class $B$, 
then we add a data member of type $Holder<B>$ to $A$. 
$Holder$ offers some useful access functions, stores the kind of relationship between $A$ and $B$,
and executes all needed @Function.addRef@ and @Function.releaseRef@ calls.
*/

template <typename TValue, typename TSpec = Tristate>
struct Holder;


//////////////////////////////////////////////////////////////////////////////
// METAFUNCTIONS
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Holder

template <typename TValue, typename TSpec>
struct Value< Holder<TValue, TSpec> >
{
	typedef typename _RemoveConst<TValue>::Type Type;
};
template <typename TValue, typename TSpec>
struct Value< Holder<TValue, TSpec> const>
{
	typedef typename _RemoveConst<TValue>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.Holder

template <typename TValue, typename TSpec>
struct Spec< Holder<TValue, TSpec> >
{
	typedef TSpec Type;
};
template <typename TValue, typename TSpec>
struct Spec< Holder<TValue, TSpec> const>
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Reference.param.T.type:Class.Holder

template <typename TValue, typename TSpec>
struct Reference< Holder<TValue, TSpec> >
{
	typedef typename Value< Holder<TValue, TSpec> >::Type & Type;
};
template <typename TValue, typename TSpec>
struct Reference< Holder<TValue, TSpec> const>
{
	typedef typename Value< Holder<TValue, TSpec> const>::Type & Type;
};


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Tristate Holder
//////////////////////////////////////////////////////////////////////////////
/**
.Spec.Tristate Holder
..cat:Holders
..summary:Holder that can be empty, dependent, or owner.
..signature:Holder<TValue, Tristate>
..param.TValue:Type of the managed object.
..general:Class.Holder
..remarks.text:A tristate holder $A$ that holds an object $B$ has one of the following states:
..remarks.text:- owner: $A$ is the owner of $B$. If $A$ is destroyed, $B$ will be destroyed automatically.
..remarks.text:- dependent: $A$ depends on $B$. $B$ should not be destroyed as long as $A$ is used.
..remarks.text:- empty: there is currently no object reference stored in the holder $A$.
..remarks.text:The state of the holder can be determined by @Function.empty@ and @Function.dependent@.
*/

template <typename TValue>
struct Holder<TValue, Tristate>
{
public:
	enum EHolderState
	{
		EMPTY = 0,
		OWNER = 1,
		DEPENDENT = ~0
	};

//	typedef typename _RemoveConst<TValue>::Type TValue_NotConst;
	typedef typename Value<Holder>::Type THostValue;

	typename _Pointer<THostValue>::Type data_value;
	EHolderState data_state;

//____________________________________________________________________________

/**
.Memfunc.Holder:
..class:Class.Holder
..summary:Constructor
..signature:Holder<TValue, TInfix> ()
..signature:Holder<TValue, TInfix> (holder)
..signature:Holder<TValue, TInfix> (value)
..param.holder:Another holder object.
..param.value:An object of type $TValue$.
..remarks.text:
The default constructor creates a holder that is in state 'empty'.
If a $value$ is passed to the constructor, the holder will be in state 'dependent'.
*/

	Holder():
		data_state(EMPTY)
	{
SEQAN_CHECKPOINT
	}
	Holder(Holder const & source_):
		data_state(EMPTY)
	{
SEQAN_CHECKPOINT
		assign(*this, source_);
	}
	Holder(typename _Parameter<THostValue>::Type value_):
		data_state(EMPTY)
	{
SEQAN_CHECKPOINT
		setValue(*this, value_);
	}
	Holder(typename _ConstParameter<THostValue>::Type value_):
		data_state(EMPTY)
	{
SEQAN_CHECKPOINT
		assignValue(*this, value_);
	}

/**
.Memfunc.~Holder:
..class:Class.Holder
..summary:Destructor
..signature:~Holder()
..remarks.text:
If the holder is in state 'owner', the holded object will be destoyed too.
*/
	~Holder()
	{
SEQAN_CHECKPOINT
		clear(*this);
	}

//____________________________________________________________________________

	Holder const &
	operator = (Holder const & source_)
	{
SEQAN_CHECKPOINT
		assign(*this, source_);
		return *this;
	}

	Holder const &
	operator = (typename _ConstParameter<THostValue>::Type value_)
	{
SEQAN_CHECKPOINT
		assignValue(*this, value_);
		return *this;
	}

	operator typename _Parameter<THostValue>::Type()
	{
SEQAN_CHECKPOINT
		return *data_value;
	}
//____________________________________________________________________________

};

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

///.Function.empty.param.object.type:Class.Holder

template <typename TValue>
inline bool
empty(Holder<TValue, Tristate> const & me)
{
SEQAN_CHECKPOINT
	return (me.data_state == Holder<TValue, Tristate>::EMPTY);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.dependent.param.object.type:Class.Holder

template <typename TValue>
inline bool
dependent(Holder<TValue, Tristate> const & me)
{
SEQAN_CHECKPOINT
	return (me.data_state == Holder<TValue, Tristate>::DEPENDENT);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.clear.param.object.type:Class.Holder
///.Function.clear.remarks.text:If $clear$ is applied on a @Class.Holder@ object,
///the state of this object is set to 'empty'.

template <typename TValue>
inline void
clear(Holder<TValue, Tristate> & me)
{
	switch (me.data_state)
	{
	case Holder<TValue, Tristate>::EMPTY:
		break;

	case Holder<TValue, Tristate>::DEPENDENT:
		{
SEQAN_CHECKPOINT
		releaseRef(_toParameter<TValue>(me.data_value));
		me.data_state = Holder<TValue, Tristate>::EMPTY;
		}
		break;

	default: /*Holder<TValue, TSpec>::OWNER*/
		{
SEQAN_CHECKPOINT
		valueDestruct(me.data_value);
		deallocate(me, me.data_value, 1);
		me.data_state = Holder<TValue, Tristate>::EMPTY;
		}
		break;
	}
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.create:
..summary:Makes an object to owner of its content.
..cat:Dependent Objects
..signature:create(holder [, object])
..param.holder:A holder object.
...type:Class.Holder
..param.object:Object from which a copy is made and stored in $holder$. (optional)
...type:Metafunction.Value.Value<Holder>::Type
..remarks.text:After this operation, $holder$ will be in state 'owner'.
If $object$ is specified, $holder$ will hold a copy of $object$ at the end of this function.
If $object$ is not specified, the action depends on the former state of $holder$:
..remarks.text:- If the state of $holder$ was 'empty', a new object is default constructed and stored into $holder$.
..remarks.text:- If the state of $holder$ was 'dependent', a copy of the former object is made and stored into $holder$. 
..remarks.text:- If the state of $holder$ was already 'owner', nothing happens.
..see:Class.Holder
*/

template <typename TValue>
inline void
create(Holder<TValue, Tristate> & me)
{
	typedef Holder<TValue, Tristate> THolder;

	switch (me.data_state)
	{
	case Holder<TValue, Tristate>::EMPTY:
		{
SEQAN_CHECKPOINT
		allocate(me, me.data_value, 1);
		valueConstruct(me.data_value);
		me.data_state = THolder::OWNER;
		}
		break;

	case THolder::DEPENDENT:
		{
SEQAN_CHECKPOINT
		typename _Parameter<TValue>::Type old_value = value(me);
		allocate(me, me.data_value, 1);
		valueConstruct(me.data_value, old_value);
		me.data_state = THolder::OWNER;
		releaseRef(old_value);
		}
		break;
	default:;
	}
}

//____________________________________________________________________________

template <typename TValue>
inline void
create(Holder<TValue, Tristate> & me,
	   typename _Parameter<TValue const>::Type value_)
{
SEQAN_CHECKPOINT

	if (me.data_state == Holder<TValue, Tristate>::OWNER)
	{
		assign(_toParameter(me.data_value), value_);
		return;
	}

	clear(me);
	allocate(me, me.data_value, 1);
	valueConstruct(me.data_value, value_);
	me.data_state = Holder<TValue, Tristate>::OWNER;
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.detach:
..summary:Makes an object independent from other objects.
..cat:Dependent Objects
..signature:detach(object)
..param.object:An object.
...type:Class.Holder
..remarks:
After this function, $object$ does not depends from any other entity outside of $object$,
like a @Function.source@ or a @Function.host@, and @Function.dependent.dependent(object)@ returns $false$ 
..see:Function.source
..see:Function.host
..see:Function.createSource
..see:Function.create
*/

template <typename TValue>
inline void
detach(Holder<TValue, Tristate> & me)
{
SEQAN_CHECKPOINT
	create(me);
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.setValue:
..cat:Content Manipulation
..summary:Makes holder dependent.
..signature:setValue(holder, object)
..param.holder:A holder object.
...type:Class.Holder
..param.object:Object from which $holder$ will be dependent.
...type:Metafunction.Value.Value<Holder>::Type
..remarks.text:After this operation, $holder$ will be dependent in state 'dependent'.
..see:Class.Holder
*/

template <typename TValue>
inline void
setValue(Holder<TValue, Tristate> & me,
		 typename _Parameter<TValue>::Type value_)
{
SEQAN_CHECKPOINT
	typedef typename Value<Holder<TValue, Tristate> >::Type THolderType;

	clear(me);
	me.data_value = _toPointer(value_);
	me.data_state = Holder<TValue, Tristate>::DEPENDENT;
	addRef(_toParameter<THolderType>(me.data_value));
}

template <typename TValue, typename TValue2>
inline void
setValue(Holder<TValue, Tristate> & me,
		 TValue2 const & value_)
{
SEQAN_CHECKPOINT
	set(value(me), value_);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.value.param.object.type:Class.Holder

template <typename TValue>
inline typename Reference<Holder<TValue, Tristate> >::Type
value(Holder<TValue, Tristate> & me)
{
SEQAN_CHECKPOINT
	typedef Holder<TValue, Tristate> THolder;

	if (empty(me))
	{
		allocate(me, me.data_value, 1);
		valueConstruct(me.data_value);
		me.data_state = THolder::OWNER;
	}

	typedef typename Value<Holder<TValue, Tristate> >::Type THolderType;
	return _toParameter<THolderType>(me.data_value);
}
template <typename TValue>
inline typename Reference<Holder<TValue, Tristate> const>::Type
value(Holder<TValue, Tristate> const & me)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!empty(me));

	typedef typename Value<Holder<TValue, Tristate> >::Type THolderType;
	return _toParameter<THolderType>(me.data_value);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.assignValue.param.object.type:Class.Holder

template <typename TValue, typename TSource>
inline void
assignValue(Holder<TValue, Tristate> & me,
			TSource const & value_)
{
SEQAN_CHECKPOINT
	typedef typename Value<Holder<TValue, Tristate> >::Type THostValue;
	if (empty(me))
	{
		create(me, value_);
	}
	else
	{
		assign(_toParameter<THostValue>(me.data_value), value_);
	}
}

//////////////////////////////////////////////////////////////////////////////

///.Function.moveValue.param.object.type:Class.Holder

template <typename TValue, typename TSource>
inline void
moveValue(Holder<TValue, Tristate> & me,
		  TSource const & value_)
{
SEQAN_CHECKPOINT
	if (empty(me))
	{
		create(me, value_);
	}
	else
	{
		move(value(me), value_);
	}
}

//////////////////////////////////////////////////////////////////////////////

///.Function.assign.param.target.type:Class.Holder
///.Function.assign.param.source.type:Class.Holder

template <typename TValue>
inline void
assign(Holder<TValue, Tristate> & target_,
	   Holder<TValue, Tristate> const & source_)
{
SEQAN_CHECKPOINT
	switch(source_.data_state)
	{
	case Holder<TValue, Tristate>::EMPTY:
		{
		clear(target_);
		}
		break;

	case Holder<TValue, Tristate>::OWNER:
		{
		assignValue(target_, value(source_));
		}
		break;

	default: /*case Holder<TValue, Tristate>::DEPENDENT*/
		{
		setValue(target_, value(source_));
		}
		break;
	}
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Simple Holder
//////////////////////////////////////////////////////////////////////////////
//??? TODO: Documentation of Simple Holder

template <typename TValue>
struct Holder<TValue, Simple>
{
	typedef typename Value<Holder>::Type THolderValue;
	typedef typename _Parameter<THolderValue>::Type THolderParameter;

	mutable THolderValue data_value;
//____________________________________________________________________________

	Holder()
	{
SEQAN_CHECKPOINT
	}
	Holder(Holder & source_):
		data_value(source_.data_value)
	{
SEQAN_CHECKPOINT
	}
	Holder(Holder const & source_):
		data_value(source_.data_value)
	{
SEQAN_CHECKPOINT
	}
	template <typename TSource>
	Holder(TSource & value_):
		data_value(value_)
	{
SEQAN_CHECKPOINT
	}
	template <typename TSource>
	Holder(TSource const & value_):
		data_value(value_)
	{
SEQAN_CHECKPOINT
	}
/*
	Holder(TValue const & value_):
		data_value(value_)
	{
SEQAN_CHECKPOINT
	}
*/
	~Holder()
	{
SEQAN_CHECKPOINT
	}

//____________________________________________________________________________

	Holder const &
	operator = (Holder const & source_)
	{
SEQAN_CHECKPOINT
		data_value = source_.data_value;
		return *this;
	}

	Holder const &
	operator = (THolderValue const & value_)
	{
SEQAN_CHECKPOINT
		data_value = value_;
		return *this;
	}

	operator THolderParameter()
	{
SEQAN_CHECKPOINT
		return *data_value;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline bool
empty(Holder<TValue, Simple> const & me)
{
SEQAN_CHECKPOINT
	return false;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline bool
dependent(Holder<TValue, Simple> const & me)
{
SEQAN_CHECKPOINT
	return false;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void
clear(Holder<TValue, Simple> & me)
{
SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void
create(Holder<TValue, Simple> & me)
{
SEQAN_CHECKPOINT
}

//____________________________________________________________________________

template <typename TValue>
inline void
create(Holder<TValue, Simple> & me,
	   TValue const & value_)
{
SEQAN_CHECKPOINT
	me.data_value = value_;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void
detach(Holder<TValue, Simple> & me)
{
SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void
setValue(Holder<TValue, Simple> & me,
		 TValue const & value_)
{
SEQAN_CHECKPOINT
	me.data_value = value_;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline typename Reference<Holder<TValue, Simple> >::Type
value(Holder<TValue, Simple> & me)
{
SEQAN_CHECKPOINT
	return me.data_value;
}
template <typename TValue>
inline typename Reference<Holder<TValue, Simple> const>::Type
value(Holder<TValue, Simple> const & me)
{
SEQAN_CHECKPOINT
	return me.data_value;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSource>
inline void
assignValue(Holder<TValue, Simple> & me,
			TSource const & value_)
{
SEQAN_CHECKPOINT
	assignValue(me.data_value, value_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSource>
inline void
moveValue(Holder<TValue, Simple> & me,
		  TSource const & value_)
{
SEQAN_CHECKPOINT
	move(me.data_value, value_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void
assign(Holder<TValue, Simple> & target_,
	   Holder<TValue, Simple> const & source_)
{
SEQAN_CHECKPOINT
	assignValue(target_, source_);
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
// New Tristate Holder that works also on pointers
//////////////////////////////////////////////////////////////////////////////

struct Tristate2;

template <typename TValue>
struct Holder<TValue, Tristate2>
{
public:
	enum EHolderState
	{
		EMPTY = 0,
		OWNER = 1,
		DEPENDENT = ~0
	};

	typedef typename Value<Holder>::Type THostValue;

	TValue * data_value;
	EHolderState data_state;

//____________________________________________________________________________

	Holder():
		data_state(EMPTY)
	{
SEQAN_CHECKPOINT
	}
	Holder(Holder const & source_):
		data_state(EMPTY)
	{
SEQAN_CHECKPOINT
		assign(*this, source_);
	}
	Holder(typename _Parameter<THostValue>::Type value_):
		data_state(EMPTY)
	{
SEQAN_CHECKPOINT
		setValue(*this, value_);
	}
	Holder(typename _ConstParameter<THostValue>::Type value_):
		data_state(EMPTY)
	{
SEQAN_CHECKPOINT
		assignValue(*this, value_);
	}
	~Holder()
	{
SEQAN_CHECKPOINT
		clear(*this);
	}

//____________________________________________________________________________

	Holder const &
	operator = (Holder const & source_)
	{
SEQAN_CHECKPOINT
		assign(*this, source_);
		return *this;
	}

	Holder const &
	operator = (TValue const & value_)
	{
SEQAN_CHECKPOINT
		assignValue(*this, value_);
		return *this;
	}

	operator TValue &()
	{
SEQAN_CHECKPOINT
		return *data_value;
	}
//____________________________________________________________________________

};

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

///.Function.empty.param.object.type:Class.Holder

template <typename TValue>
inline bool
empty(Holder<TValue, Tristate2> const & me)
{
SEQAN_CHECKPOINT
	return (me.data_state == Holder<TValue, Tristate2>::EMPTY);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.dependent.param.object.type:Class.Holder

template <typename TValue>
inline bool
dependent(Holder<TValue, Tristate2> const & me)
{
SEQAN_CHECKPOINT
	return (me.data_state == Holder<TValue, Tristate2>::DEPENDENT);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.clear.param.object.type:Class.Holder
///.Function.clear.remarks.text:If $clear$ is applied on a @Class.Holder@ object,
///the state of this object is set to 'empty'.

template <typename TValue>
inline void
clear(Holder<TValue, Tristate2> & me)
{
	switch (me.data_state)
	{
	case Holder<TValue, Tristate2>::EMPTY:
		break;

	case Holder<TValue, Tristate2>::DEPENDENT:
		{
SEQAN_CHECKPOINT
		releaseRef(*(me.data_value));
		me.data_state = Holder<TValue, Tristate2>::EMPTY;
		}
		break;

	default: /*Holder<TValue, TSpec>::OWNER*/
		{
SEQAN_CHECKPOINT
		valueDestruct(me.data_value);
		deallocate(me, me.data_value, 1);
		me.data_state = Holder<TValue, Tristate2>::EMPTY;
		}
		break;
	}
}

//////////////////////////////////////////////////////////////////////////////


template <typename TValue>
inline void
create(Holder<TValue, Tristate2> & me)
{
	typedef Holder<TValue, Tristate2> THolder;

	switch (me.data_state)
	{
	case Holder<TValue, Tristate2>::EMPTY:
		{
SEQAN_CHECKPOINT
		allocate(me, me.data_value, 1);
		valueConstruct(me.data_value);
		me.data_state = THolder::OWNER;
		}
		break;

	case THolder::DEPENDENT:
		{
SEQAN_CHECKPOINT
		TValue & old_value = value(me);
		allocate(me, me.data_value, 1);
		valueConstruct(me.data_value, old_value);
		me.data_state = THolder::OWNER;
		releaseRef(old_value);
		}
		break;
	default:;
	}
}

//____________________________________________________________________________

template <typename TValue>
inline void
create(Holder<TValue, Tristate2> & me,
	   TValue const & value_)
{
SEQAN_CHECKPOINT

	if (me.data_state == Holder<TValue, Tristate2>::OWNER)
	{
		assign(*(me.data_value), value_);
		return;
	}

	clear(me);
	allocate(me, me.data_value, 1);
	valueConstruct(me.data_value, value_);
	me.data_state = Holder<TValue, Tristate2>::OWNER;
}

//////////////////////////////////////////////////////////////////////////////


template <typename TValue>
inline void
detach(Holder<TValue, Tristate2> & me)
{
SEQAN_CHECKPOINT
	create(me);
}

//////////////////////////////////////////////////////////////////////////////


template <typename TValue>
inline void
setValue(Holder<TValue, Tristate2> & me,
		 TValue & value_)
{
SEQAN_CHECKPOINT
	typedef typename Value<Holder<TValue, Tristate2> >::Type THolderType;

	clear(me);
	me.data_value = & value_;
	me.data_state = Holder<TValue, Tristate2>::DEPENDENT;
	addRef(value_);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.value.param.object.type:Class.Holder

template <typename TValue>
inline typename Reference<Holder<TValue, Tristate2> >::Type
value(Holder<TValue, Tristate2> & me)
{
SEQAN_CHECKPOINT
	typedef Holder<TValue, Tristate2> THolder;

	if (empty(me))
	{
		allocate(me, me.data_value, 1);
		valueConstruct(me.data_value);
		me.data_state = THolder::OWNER;
	}

	typedef typename Value<Holder<TValue, Tristate2> >::Type THolderType;
	return *(me.data_value);
}
template <typename TValue>
inline typename Reference<Holder<TValue, Tristate2> const>::Type
value(Holder<TValue, Tristate2> const & me)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!empty(me));

	return *(me.data_value);
}


//////////////////////////////////////////////////////////////////////////////

///.Function.assignValue.param.object.type:Class.Holder

template <typename TValue, typename TSource>
inline void
assignValue(Holder<TValue, Tristate2> & me,
			TSource const & value_)
{
SEQAN_CHECKPOINT
	typedef typename Value<Holder<TValue, Tristate2> >::Type THostValue;
	if (empty(me))
	{
		create(me, value_);
	}
	else
	{
		assign(*(me.data_value), value_);
	}
}

//////////////////////////////////////////////////////////////////////////////

///.Function.moveValue.param.object.type:Class.Holder

template <typename TValue, typename TSource>
inline void
moveValue(Holder<TValue, Tristate2> & me,
		  TSource const & value_)
{
SEQAN_CHECKPOINT
	if (empty(me))
	{
		create(me, value_);
	}
	else
	{
		move(value(me), value_);
	}
}

//////////////////////////////////////////////////////////////////////////////

///.Function.assign.param.target.type:Class.Holder
///.Function.assign.param.source.type:Class.Holder

template <typename TValue>
inline void
assign(Holder<TValue, Tristate2> & target_,
	   Holder<TValue, Tristate2> const & source_)
{
SEQAN_CHECKPOINT
	switch(source_.data_state)
	{
	case Holder<TValue, Tristate2>::EMPTY:
		{
		clear(target_);
		}
		break;

	case Holder<TValue, Tristate2>::OWNER:
		{
		assignValue(target_, value(source_));
		}
		break;

	default: /*case Holder<TValue, Tristate2>::DEPENDENT*/
		{
		setValue(target_, value(source_));
		}
		break;
	}
}
//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...


