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
  $Id: string_cstyle.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEQUENCE_CSTYLE_H
#define SEQAN_HEADER_SEQUENCE_CSTYLE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//string that exports an interface to cstyle strings.

/**
.Spec.CStyle String:
..cat:Strings
..general:Class.String
..summary:Allows adaption of strings to C-style strings.
..signature:String<TValue, CStyle>
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
...note:$TValue$ must be a simple type.???(Link)
..remarks:
...text:The purpose of this class is to access to the content of a sequence 
in a "zero terminated string" style. 
This can be useful if SEQAN classes has to be integrated in programs that use $char$ arrays
to store strings.
Instances of $String<TValue, CStyle>$ can implicitely converted to a $TValue *$ that
points to a zero terminated CStyle of $TValue$. 
...text:The stored c-style string object can be set by constructors or assignment. 
The content of a c-style string can eighter be stored in a separate buffer, that is the source string
is copied. Or the buffer of the source string itself is used instead, in this case the c-style string
depends on the source string and gets invalid as soon as the buffer of the source string is destroyed.
...text:Hence, this class is a kind of adaptor from an arbitrary SEQAN string to char arrays.
Of course, the opposite way is possible too: 
Read @Adaption.char array.here@ about adapting char arrays to SEQAN strings.
..example:
...code://Create a string str:
String<char> str = "this is a test string";

//Create a c-style string object for str:
String<char, CStyle> c_style = str;

//Now use c_style as char array:
strcmp(c_style, "compare it to this string");
...text:If the c-style string is needed only temporarily, the function $toCString$ can be used:
...code:String<char> str = "this is a test string";
strcmp(toCString(str), "compare it to this string");
*/

struct CStyle;

template <typename TValue>
class String <TValue, CStyle >
{
protected:
	TValue * data_begin;
	TValue * data_end;
	size_t data_size; //if data_size > 0, then the buffer is owned by me and must be deallocated

public:
	static TValue EMPTY_STRING;
//____________________________________________________________________________

public:
	String():
		data_begin(&EMPTY_STRING),
		data_end(&EMPTY_STRING),
		data_size(0)
	{
SEQAN_CHECKPOINT
	}

//____________________________________________________________________________

	template <typename TString>
	String(TString & str):
		data_size(0)
	{
SEQAN_CHECKPOINT
		assign(*this, str);
	}
	template <typename TString>
	String(TString const & str):
		data_size(0)
	{
SEQAN_CHECKPOINT
		assign(*this, str);
	}

	String(String & str):
		data_size(0)
	{
SEQAN_CHECKPOINT
		assign(*this, str);
	}
	String(String const & str):
		data_size(0)
	{
SEQAN_CHECKPOINT
		assign(*this, str);
	}

	String(TValue * str):
		data_begin(str),
		data_end(end(str)),
		data_size(0)
	{
SEQAN_CHECKPOINT
	}

//____________________________________________________________________________

	template <typename TString>
	String & operator = (TString & str)
	{
SEQAN_CHECKPOINT
		assign(*this, str);
		return *this;
	}
	template <typename TString>
	String & operator = (TString const & str)
	{
SEQAN_CHECKPOINT
		assign(*this, str);
		return *this;
	}
	String & operator = (String & str)
	{
SEQAN_CHECKPOINT
		assign(*this, str);
		return *this;
	}
	String & operator = (String const & str)
	{
SEQAN_CHECKPOINT
		assign(*this, str);
		return *this;
	}

	~String()
	{
SEQAN_CHECKPOINT
		clear(*this);
	}

//____________________________________________________________________________

	operator TValue * ()
	{
SEQAN_CHECKPOINT
		return data_begin;
	}

	operator TValue const * () const
	{
SEQAN_CHECKPOINT
		return data_begin;
	}

//____________________________________________________________________________

	friend inline void 
	move(
		String & target,
		String & source)
	{
SEQAN_CHECKPOINT
		clear(target);

		target.data_begin = source.data_begin;
		target.data_end = source.data_end;
		target.data_size = source.data_size;

		source.data_begin = 0;
		source.data_end = 0;
		source.data_size = 0;
	}

//____________________________________________________________________________

	friend inline typename Iterator<String, Standard>::Type
	begin(String & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return me.data_begin;
	}
	friend inline typename Iterator<String const, Standard>::Type
	begin(String const & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return me.data_begin;
	}

//____________________________________________________________________________

	friend inline void
	_setBegin(String & me, TValue * new_begin)
	{
SEQAN_CHECKPOINT
		me.data_begin = new_begin;
	}

//____________________________________________________________________________

	friend inline typename Iterator<String, Standard>::Type
	end(String & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return me.data_end;
	}
	friend inline typename Iterator<String const, Standard>::Type
	end(String const & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return me.data_end;
	}

//____________________________________________________________________________

	friend inline void
	_setEnd(String & me, TValue * new_end)
	{
SEQAN_CHECKPOINT
		me.data_end = new_end;
		*new_end = TValue(); //??? ist das wirklich sinnvoll fuer typen, die weder char noch wchar_t sind?
	}

//____________________________________________________________________________

	friend inline size_t 
	capacity(String const & me)
	{
SEQAN_CHECKPOINT
		if (me.data_size) return me.data_size -1;
		else return me.data_end - me.data_begin;
	}


//____________________________________________________________________________

	friend inline void
	clear(String & me)
	{
		if (me.data_size)
		{
SEQAN_CHECKPOINT
//			arrayDestruct(me, length(me)); 
			deallocate(me, me.data_begin, me.data_size);
			me.data_size = 0;
		}
		me.data_begin = me.data_end = &EMPTY_STRING;
	}
//____________________________________________________________________________

//??? TODO: reserve

//____________________________________________________________________________

///.Internal._reallocateStorage.param.object.type:Spec.CStyle String
///.Internal._reallocateStorage.param.resize_tag.remarks:@Spec.CStyle String@ only supports @Tag.Overflow Strategy.exact@.
//this function works also for dependent buffers
	friend inline TValue *
	_reallocateStorage(
		String & me, 
		size_t new_capacity,
		Exact)
	{
SEQAN_CHECKPOINT
		TValue * return_value;
		if (me.data_size)
		{//dependent
			return_value = me.data_begin;
		}
		else
		{//not dependent
			return_value = 0;
		}

		me.data_size = new_capacity + 1; //+1 for zero termination
		allocate(me, me.data_begin, me.data_size, TagAllocateStorage());
		return return_value;
	}
//____________________________________________________________________________

///.Internal._deallocateStorage.param.object.type:Spec.CStyle String

	friend inline void 
	_deallocateStorage(
		String & me, 
		TValue * ptr, 
		size_t capacity)
	{
SEQAN_CHECKPOINT
		size_t size = capacity + 1;
		deallocate(me, ptr, size, TagAllocateStorage());
	}

//____________________________________________________________________________

/**
.Function.dependent:
..summary:Test whether object depends on other objects.
..cat:Dependent Objects
..signature:bool dependent(object)
..param.object:An object.
...type:Spec.CStyle String
..returns:$true$ if $object$ depends one some other object, $false$ otherwise.
..remarks:An object "$a$" depends on another object "$b$", if changing "$b$" can invalidate "$a$";
especially the destruction of "$b$" invalidates "$a$".
*/
	friend inline bool
	dependent(String & me)
	{
SEQAN_CHECKPOINT
		return (me.data_size == 0);
	}
//____________________________________________________________________________

//special implementation for char array sources
	friend inline void
	assign(String & target,
		TValue * source)
	{
	SEQAN_CHECKPOINT
		clear(target);
		target.data_begin = source;
		target.data_end = end(source);
	}

//____________________________________________________________________________

};

//////////////////////////////////////////////////////////////////////////////
// Define the static member

template <typename TValue>
TValue String<TValue, CStyle >::EMPTY_STRING = TValue();

//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
struct DefaultOverflowImplicit<String<TValue, CStyle> >
{
	typedef Exact Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
struct IsContiguous< String<TValue, CStyle > >
{
    typedef True Type;
	enum { VALUE = true };
};

//////////////////////////////////////////////////////////////////////////////
// assign
//////////////////////////////////////////////////////////////////////////////

template <typename TTargetValue, typename TSource, typename TExpand>
inline void
assign(String<TTargetValue, CStyle> & target,
	   TSource & source,
	   Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	create(target, source, tag);
}

template <typename TTargetValue, typename TSource, typename TExpand>
inline void
assign(String<TTargetValue, CStyle> & target,
	   TSource const & source,
	   Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	create(target, source, tag);
}

template <typename TTargetValue, typename TSource, typename TSize, typename TExpand>
inline void
assign(String<TTargetValue, CStyle> & target,
	   TSource & source,
	   TSize limit,
	   Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	create(target, source, tag);
}

template <typename TTargetValue, typename TSource, typename TSize, typename TExpand>
inline void
assign(String<TTargetValue, CStyle> & target,
	   TSource const & source,
	   TSize limit,
	   Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	create(target, source, limit, tag);
}


//____________________________________________________________________________
//this variant is a workaround for the "const array"-bug of VC++

template <typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
assign(String<TTargetValue, CStyle> & target,
	   TSourceValue const * source, 
	   Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	create(target, source, tag);
}

template <typename TTargetValue, typename TSourceValue, typename TSize, typename TExpand>
inline void
assign(String<TTargetValue, CStyle> & target,
	   TSourceValue const * source,
	   TSize limit,
	   Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	create(target, source, limit, tag);
}

//////////////////////////////////////////////////////////////////////////////

//If source is non-const String, then there could be the possibility
//to use the source buffer

template <typename TExpand, bool IS_CONTIGUOUS>
struct _Assign_String_2_StringArray;

//____________________________________________________________________________

template <typename TExpand>
struct _Assign_String_2_StringArray<TExpand, true>
{
	template <typename TValue, typename TSourceSpec>
	static inline void
	assign_(String<TValue, CStyle> & target,
		String<TValue, TSourceSpec> & source)
	{
		if (capacity(source) > length(source))
		{//use source's buffer
SEQAN_CHECKPOINT
			clear(target);
			_setBegin(target, begin(source));
			_setEnd(target, end(source));
		}
		else
		{
			create(target, source, TExpand());
		}
	}

//special treatment of char:
//_computeSize4Capacity is specialized for char such that there
//is enough place for the zero termination

	template <typename TSourceSpec>
	static inline void
	assign_(String<char, CStyle> & target,
		String<char, TSourceSpec> & source)
	{
SEQAN_CHECKPOINT
		clear(target);
		_setBegin(target, begin(source));
		_setEnd(target, end(source));
	}
};

//____________________________________________________________________________

template <typename TExpand>
struct _Assign_String_2_StringArray<TExpand, false>
{
	template <typename TValue, typename TSourceSpec>
	static inline void
	assign_(String<TValue, CStyle> & target,
		String<TValue, TSourceSpec> & source)
	{
SEQAN_CHECKPOINT
		create(target, source, TExpand());
	}
};

//____________________________________________________________________________

template <typename TValue, typename TSourceSpec, typename TExpand>
inline void
assign(String<TValue, CStyle> & target,
	String<TValue, TSourceSpec> & source,
	Tag<TExpand> const)
{
	typedef String<TValue, TSourceSpec> TSource;
	_Assign_String_2_StringArray<Tag<TExpand> const, IsContiguous<TSource>::VALUE>::assign_(target, source);
}


//////////////////////////////////////////////////////////////////////////////
// create
//////////////////////////////////////////////////////////////////////////////

//see basic_holder
/**
.Function.create:
..signature:create(target, source [, limit] [,resize_tag])
..param.target: Gets a copy of the content of $source$.
...type:Spec.CStyle String
..param.source: Is copied to $target$.
..param.limit: The maximal length of $target$ after the operation. (optional)
..param.resize_tag: Specifies the strategy that is applied if $target$ has not enough capacity to store the complete content. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowImplicit@ of the $target$ type.
..remarks.text:It is guaranteed, that after calling this function $source$ and $target$ can be used independently.
..see:Spec.CStyle String
*/

template <typename TExpand>
struct _Create_ArrayString_Expand
{
	template <typename TTarget, typename TSource>
	static inline void
	create_(TTarget & target, 
		TSource & source)
	{
		typename Size<TTarget>::Type source_length = length(source);
		if (dependent(target) || (capacity(target) < source_length))
		{
SEQAN_CHECKPOINT
			typename Size<TTarget>::Type old_target_capacity = capacity(target);
			typename Value<TTarget>::Type * buf = _reallocateStorage(target, source_length, TExpand());
			if (buf)
			{
				_deallocateStorage(target, buf, old_target_capacity);
			}
		}
		if (length(source) > 0)
		{
			assignValue(begin(target, Standard()), 0); //set target length to 0
			assign(begin(target, Standard()), source, Insist());
			_setEnd(target, begin(target) + source_length);
		}
	}

	template <typename TTarget, typename TSource, typename TLimit>
	static inline void
	create_(TTarget & target, 
		TSource & source,
		TLimit limit)
	{
		typename Size<TTarget>::Type copy_length = length(source);
		if (limit < copy_length)
		{
			copy_length = limit;
		}
		if (dependent(target) || (capacity(target) < copy_length))
		{
SEQAN_CHECKPOINT
			typename Size<TTarget>::Type old_target_capacity = capacity(target);
			TTarget * buf = _reallocateStorage(target, copy_length, TExpand());
			if (buf)
			{
				_deallocateStorage(target, buf, old_target_capacity);
			}
		}
		assign(begin(target, Standard()), source, copy_length, Insist());
		_setEnd(target, begin(target, Standard()) + copy_length);
	}
};
//____________________________________________________________________________

template <typename TExpand>
struct _Create_ArrayString;
//____________________________________________________________________________

template <>
struct _Create_ArrayString<Insist>
{
	template <typename TTarget, typename TSource>
	static inline void
	create_(TTarget & target, 
		TSource & source)
	{
SEQAN_CHECKPOINT
		typename Size<TTarget>::Type source_length = length(source);
		if (dependent(target))
		{
			TTarget * buf = _reallocateStorage(target, source_length, Exact());
		}
		assign(begin(target, Standard()), source, source_length, Insist());
		_setEnd(target, begin(target, Standard()) + source_length);
	}

	template <typename TTarget, typename TSource, typename TSize>
	static inline void
	create_(TTarget & target, 
		TSource & source,
		TSize limit)
	{
SEQAN_CHECKPOINT
		typename Size<TTarget>::Type copy_size = length(source);
		if (limit < copy_size)
		{
			copy_size = limit;
		}
		if (dependent(target))
		{
			TTarget * buf = _reallocateStorage(target, copy_size, Exact());
		}
		assign(begin(target, Standard()), source, copy_size, Insist());
		_setEnd(target, begin(target, Standard()) + copy_size);
	}
};

//____________________________________________________________________________

template <>
struct _Create_ArrayString<Limit>
{
	template <typename TTarget, typename TSource>
	static inline void
	create_(TTarget & target, 
		TSource & source)
	{
SEQAN_CHECKPOINT
		_Create_ArrayString<Insist>::create_(target, source, capacity(target));
	}

	template <typename TTarget, typename TSource, typename TSize>
	static inline void
	create_(TTarget & target, 
		TSource & source,
		TSize & limit)
	{
SEQAN_CHECKPOINT
		typename Size<TTarget>::Type copy_size = capacity(target);
		if (copy_size > limit)
		{
			copy_size = limit;
		}
		_Create_ArrayString<Insist>::create_(target, source, copy_size);
	}
};
//____________________________________________________________________________

template <>
struct _Create_ArrayString<Exact>:
	_Create_ArrayString_Expand<Exact>
{
};

//____________________________________________________________________________

template <>
struct _Create_ArrayString<Generous>:
	_Create_ArrayString_Expand<Generous>
{
};

//____________________________________________________________________________

template <typename TTargetValue, typename TSource>
inline void
create(String<TTargetValue, CStyle> & target,
	   TSource & source)
{
SEQAN_CHECKPOINT
	typedef String<TTargetValue, CStyle> TTarget;
	create(target, source, typename DefaultOverflowImplicit<TTarget>::Type()); 
}

template <typename TTargetValue, typename TSource, typename TSize>
inline void
create(String<TTargetValue, CStyle> & target,
	   TSource & source,
	   TSize limit)
{
SEQAN_CHECKPOINT
	typedef String<TTargetValue, CStyle> TTarget;
	create(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type()); 
}

//____________________________________________________________________________

template <typename TTargetValue, typename TSource, typename TExpand>
inline void
create(String<TTargetValue, CStyle> & target,
	   TSource & source, 
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Create_ArrayString<Tag<TExpand> const>::create_(target, source);
}

template <typename TTargetValue, typename TSource, typename TSize, typename TExpand>
inline void
create(String<TTargetValue, CStyle> & target,
	   TSource & source,
	   TSize limit,
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Create_ArrayString<Tag<TExpand> const>::create_(target, source, limit);
}

template <typename TTargetValue, typename TSource, typename TExpand>
inline void
create(String<TTargetValue, CStyle> & target,
	   TSource const & source, 
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Create_ArrayString<Tag<TExpand> const>::create_(target, source);
}

template <typename TTargetValue, typename TSource, typename TSize, typename TExpand>
inline void
create(String<TTargetValue, CStyle> & target,
	   TSource const & source,
	   TSize limit,
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Create_ArrayString<Tag<TExpand> const>::create_(target, source, limit);
}

//____________________________________________________________________________
//this variant is a workaround for the "const array"-bug of VC++

template <typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
create(String<TTargetValue, CStyle> & target,
	   TSourceValue const * source, 
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Create_ArrayString<Tag<TExpand> const>::create_(target, source);
}

template <typename TTargetValue, typename TSourceValue, typename TSize, typename TExpand>
inline void
create(String<TTargetValue, CStyle> & target,
	   TSourceValue const * source,
	   TSize limit,
	   Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	_Create_ArrayString<Tag<TExpand> const>::create_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
// Shotcut

/**
.Function.toCString:
..cat:Containers
..summary:Access sequence as c-style string.
..signature:toCString(object)
..param.object:A sequence.
...type:Class.String
...type:Adaption.char array
..returns:A temporary @Spec.CStyle String@ object of that was constructed for $object$.
..remarks:Notational sugar.
*/

template <typename T>
inline typename Value<T>::Type *
toCString(T & me)
{
SEQAN_CHECKPOINT
	return String<typename Value<T>::Type, CStyle>(me);
}


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

//____________________________________________________________________________

#endif //#ifndef SEQAN_HEADER_...
