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
  $Id: string_value_expand.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEQUENCE_STRING_VALUEEXPAND_H
#define SEQAN_HEADER_SEQUENCE_STRING_VALUEEXPAND_H


namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Tags

template <typename THost, typename TMap, typename TSpec = Default>
struct ValueExpand;

struct ValueExpandIter;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename THost, typename TMap, typename TSpec>
class String<TValue, ValueExpand<THost, TMap, TSpec> >
{
//____________________________________________________________________________
private:

	typedef typename Value<String>::Type TLargeValue;
	typedef typename Value<THost>::Type TSmallValue;

	Holder<THost> data_host;
	Holder<TMap> data_map;

//____________________________________________________________________________

public:
	String() 
	{
	}
	String(String const & other_):
		data_host(other_.data_host),
		data_map(other_.data_map)
	{
	}
	~String()
	{
	}
	String const &
	operator = (String const & other_)
	{
		data_host = other_.data_host;
		data_map = other_.data_map;
	}
//____________________________________________________________________________

	template <typename TPos>
	inline typename Reference<String>::Type
	operator [](TPos pos)
	{
		return value(*this, pos);
	}

//____________________________________________________________________________

	friend inline Holder<THost> &
	_dataHost(String & me)
	{
		return me.data_host;
	}
//____________________________________________________________________________

	friend inline Holder<TMap> &
	_dataMap(String & me)
	{
		return me.data_map;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Metafunctions

template <typename TValue, typename THost, typename TMap, typename TSpec>
struct Value<String< TValue, ValueExpand<THost, TMap, TSpec> > >
{
	typedef TValue Type;
};
template <typename TValue, typename THost, typename TMap, typename TSpec>
struct Value<String< TValue, ValueExpand<THost, TMap, TSpec> > const>
{
	typedef TValue const Type;
};

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec>
struct GetValue<String< TValue, ValueExpand<THost, TMap, TSpec> > >:
	Value<String< TValue, ValueExpand<THost, TMap, TSpec> > >
{
};
template <typename TValue, typename THost, typename TMap, typename TSpec>
struct GetValue<String< TValue, ValueExpand<THost, TMap, TSpec> > const>:
	Value<String< TValue, ValueExpand<THost, TMap, TSpec> > const>
{
};

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec>
struct Reference<String< TValue, ValueExpand<THost, TMap, TSpec> > >	
{
	typedef String< TValue, ValueExpand<THost, TMap, TSpec> > TMe;
	typedef typename Iterator<TMe, Standard>::Type TIterator;
	typedef Proxy<IteratorProxy<TIterator> > Type;
};
template <typename TValue, typename THost, typename TMap, typename TSpec>
struct Reference<String< TValue, ValueExpand<THost, TMap, TSpec> > const>	
{
	typedef String< TValue, ValueExpand<THost, TMap, TSpec> > const TMe;
	typedef typename Iterator<TMe, Standard>::Type TIterator;
	typedef Proxy<IteratorProxy<TIterator> > Type;
};

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec, typename TIteratorSpec>
struct Iterator<String< TValue, ValueExpand<THost, TMap, TSpec> >, TIteratorSpec>	
{
	typedef ValueExpand<THost, TMap, TSpec> TValueExpand;
	typedef String< TValue, TValueExpand> TMe;
	typedef typename Iterator<THost, Standard>::Type THostIterator;

	typedef Iter<TMe, AdaptorIterator<THostIterator, ValueExpandIter> > Type;
};
template <typename TValue, typename THost, typename TMap, typename TSpec, typename TIteratorSpec>
struct Iterator<String< TValue, ValueExpand<THost, TMap, TSpec> > const, TIteratorSpec>	
{
	typedef ValueExpand<THost, TMap, TSpec> TValueExpand;
	typedef String< TValue, TValueExpand> const TMe;
	typedef typename Iterator<THost, Standard>::Type THostIterator;

	typedef Iter<TMe, AdaptorIterator<THostIterator, ValueExpandIter> > Type;
};

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec>
struct DefaultOverflowImplicit<String< TValue, ValueExpand<THost, TMap, TSpec> > >:
	DefaultOverflowImplicit< typename Host<String< TValue, ValueExpand<THost, TMap, TSpec> > >::Type >
{
};
template <typename TValue, typename THost, typename TMap, typename TSpec>
struct DefaultOverflowImplicit<String< TValue, ValueExpand<THost, TMap, TSpec> > const>:
	DefaultOverflowImplicit< typename Host<String< TValue, ValueExpand<THost, TMap, TSpec> > const>::Type >
{
};

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec>
struct DefaultOverflowExplicit<String< TValue, ValueExpand<THost, TMap, TSpec> > >:
	DefaultOverflowExplicit< typename Host<String< TValue, ValueExpand<THost, TMap, TSpec> > >::Type >
{
};
template <typename TValue, typename THost, typename TMap, typename TSpec>
struct DefaultOverflowExplicit<String< TValue, ValueExpand<THost, TMap, TSpec> > const>:
	DefaultOverflowExplicit< typename Host<String< TValue, ValueExpand<THost, TMap, TSpec> > const>::Type >
{
};

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec>
struct IsContiguous<String< TValue, ValueExpand<THost, TMap, TSpec> > >
{
    typedef False Type;
	enum { VALUE = false };
};

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec>
struct Host<String< TValue, ValueExpand<THost, TMap, TSpec> > >
{
	typedef THost Type;
};
template <typename TValue, typename THost, typename TMap, typename TSpec>
struct Host<String< TValue, ValueExpand<THost, TMap, TSpec> > const>
{
	typedef THost const Type;
};

//____________________________________________________________________________

template <typename T>
struct MapType;

template <typename TValue, typename THost, typename TMap, typename TSpec>
struct MapType<String< TValue, ValueExpand<THost, TMap, TSpec> > >
{
	typedef TMap Type;
};
template <typename TValue, typename THost, typename TMap, typename TSpec>
struct MapType<String< TValue, ValueExpand<THost, TMap, TSpec> > const>
{
	typedef TMap const Type;
};

//////////////////////////////////////////////////////////////////////////////
// Special Functions

template <typename TValue, typename THost, typename TMap, typename TSpec>
inline TMap &
_getMap(String< TValue, ValueExpand<THost, TMap, TSpec> > & me)
{
	return value(_dataMap(me));
}
template <typename TValue, typename THost, typename TMap, typename TSpec>
inline TMap &
_getMap(String< TValue, ValueExpand<THost, TMap, TSpec> > const & me)
{
	return value(_dataMap(me));
}

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec>
inline typename Value<THost>::Type
_getValueExpandFlagValue(String< TValue, ValueExpand<THost, TMap, TSpec> > const & me)
{
	typedef typename Value<THost>::Type TSmallValue;
	return supremumValue<TSmallValue>();
}

//////////////////////////////////////////////////////////////////////////////
// Public Functions


template <typename TValue, typename THost, typename TMap, typename TSpec>
inline void const * 
id(String< TValue, ValueExpand<THost, TMap, TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return id(host(me));
}

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TPos, typename TSpec, typename TTag>
inline typename Iterator<String< TValue, ValueExpand<THost, TMap, TSpec> >, Tag<TTag> const>::Type 
iter(String< TValue, ValueExpand<THost, TMap, TSpec> > & me,
	 TPos pos_,
	 Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	typedef String< TValue, ValueExpand<THost, TMap, TSpec> > TMe;
	typedef typename Iterator<TMe, Tag<TTag> const>::Type TIterator;
	return TIterator(me, begin(host(me), Standard()) + pos_);
}
template <typename TValue, typename THost, typename TMap, typename TPos, typename TSpec, typename TTag>
inline typename Iterator<String< TValue, ValueExpand<THost, TMap, TSpec> > const, Tag<TTag> const>::Type 
iter(String< TValue, ValueExpand<THost, TMap, TSpec> > const & me,
	 TPos pos_,
	 Tag<TTag> const tag_ )
{
SEQAN_CHECKPOINT
	typedef String< TValue, ValueExpand<THost, TMap, TSpec> > const TMe;
	typedef typename Iterator<TMe, Tag<TTag> const>::Type TIterator;
	return TIterator(me, begin(host(me), Standard()) + pos_);
}

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec, typename TTag>
inline typename Iterator<String< TValue, ValueExpand<THost, TMap, TSpec> >, Tag<TTag> const>::Type 
begin(String< TValue, ValueExpand<THost, TMap, TSpec> > & me,
	  Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return iter(me, 0, tag_);
}
template <typename TValue, typename THost, typename TMap, typename TSpec, typename TTag>
inline typename Iterator<String< TValue, ValueExpand<THost, TMap, TSpec> > const, Tag<TTag> const>::Type 
begin(String< TValue, ValueExpand<THost, TMap, TSpec> > const & me,
	  Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return iter(me, 0, tag_);
}

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec, typename TTag>
inline typename Iterator<String< TValue, ValueExpand<THost, TMap, TSpec> >, Tag<TTag> const>::Type 
end(String< TValue, ValueExpand<THost, TMap, TSpec> > & me,
	Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return iter(me, length(me), tag_);
}
template <typename TValue, typename THost, typename TMap, typename TSpec, typename TTag>
inline typename Iterator<String< TValue, ValueExpand<THost, TMap, TSpec> > const, Tag<TTag> const>::Type 
end(String< TValue, ValueExpand<THost, TMap, TSpec> > const & me,
	Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return iter(me, length(me), tag_);
}

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec, typename TPos>
inline typename Reference<String< TValue, ValueExpand<THost, TMap, TSpec> > >::Type
value(String< TValue, ValueExpand<THost, TMap, TSpec> > & me, 
	  TPos pos)
{
SEQAN_CHECKPOINT
	
	return *iter(me, pos, Standard());
} 
template <typename TValue, typename THost, typename TMap, typename TSpec, typename TPos>
inline typename Reference<String< TValue, ValueExpand<THost, TMap, TSpec> > const>::Type
value(String< TValue, ValueExpand<THost, TMap, TSpec> > const & me, 
	  TPos pos)
{
SEQAN_CHECKPOINT
	
	return *iter(me, pos, Standard());
} 

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec>
inline typename Size<String< TValue, ValueExpand<THost, TMap, TSpec> > const>::Type
length(String< TValue, ValueExpand<THost, TMap, TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return length(host(me));
}

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec>
inline typename Size<String< TValue, ValueExpand<THost, TMap, TSpec> > const>::Type
capacity(String< TValue, ValueExpand<THost, TMap, TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return capacity(host(me));
}

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec, typename TSize>
inline TSize
resize(String< TValue, ValueExpand<THost, TMap, TSpec> > & me,
	   TSize new_length)
{
SEQAN_CHECKPOINT
	return resize(host(me), new_length);
}

//____________________________________________________________________________

template <typename TValue, typename THost, typename TMap, typename TSpec, typename TSize, typename TExpand>
inline TSize 
reserve(String< TValue, ValueExpand<THost, TMap, TSpec> > & me,
		TSize new_length,
		Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	return reserve(host(me), new_length, tag);
}

//////////////////////////////////////////////////////////////////////////////
// Iterator for ValueExpandIterator string: Subclass of AdaptorIterator
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator>
inline typename Value<TContainer>::Type
getValue(Iter<TContainer, AdaptorIterator<TIterator, ValueExpandIter> > & me)
{
SEQAN_CHECKPOINT
	typedef typename Host<TContainer>::Type THost;
	typedef typename Value<THost>::Type TSmallValue;
	TSmallValue c = value(hostIterator(me));
	if (c == _getValueExpandFlagValue(container(me)))
	{//value is large 
		return _getMap(container(me))[position(me)];
	}
	else
	{
		return c;
	}
}
template <typename TContainer, typename TIterator>
inline void
getValue(Iter<TContainer, AdaptorIterator<TIterator, ValueExpandIter> > const & me)
{
SEQAN_CHECKPOINT
	typedef typename Host<TContainer>::Type THost;
	typedef typename Value<THost>::Type TSmallValue;
	TSmallValue c = value(hostIterator(me));
	if (c == _getValueExpandFlagValue(container(me)))
	{//value is large 
		return _getMap(container(me))[position(me)];
	}
	else
	{
		return c;
	}
}

//____________________________________________________________________________

template <typename TContainer, typename TIterator>
inline typename Reference<Iter<TContainer, AdaptorIterator<TIterator, ValueExpandIter> > >::Type 
value(Iter<TContainer, AdaptorIterator<TIterator, ValueExpandIter> > & me)
{
SEQAN_CHECKPOINT
	typedef Iter<TContainer, AdaptorIterator<TIterator, ValueExpandIter> > TMe;
	return typename Reference<TMe>::Type(me);
}
template <typename TContainer, typename TIterator>
inline typename Reference<Iter<TContainer, AdaptorIterator<TIterator, ValueExpandIter> > >::Type 
value(Iter<TContainer, AdaptorIterator<TIterator, ValueExpandIter> > const & me)
{
SEQAN_CHECKPOINT
	typedef Iter<TContainer, AdaptorIterator<TIterator, ValueExpandIter> > const TMe;
	return typename Reference<TMe>::Type(me);
}

//____________________________________________________________________________

template <typename TContainer, typename TIterator, typename TValue>
inline void
assignValue(Iter<TContainer, AdaptorIterator<TIterator, ValueExpandIter> > & me,
			TValue const & _value)
{
SEQAN_CHECKPOINT
	typedef typename Host<TContainer>::Type THost;
	typedef typename Value<THost>::Type TSmallValue;
	
	TSmallValue flagValue = _getValueExpandFlagValue(container(me));
	if (_value >= flagValue)
	{//use map to store LargeValue
		assignValue(hostIterator(me), flagValue);
		_getMap(container(me))[position(me)] = _value;
	}
	else
	{
		assignValue(hostIterator(me), _value);
	}
}
template <typename TContainer, typename TIterator, typename TValue>
inline void
assignValue(Iter<TContainer, AdaptorIterator<TIterator, ValueExpandIter> > const & me,
			TValue const & _value)
{
SEQAN_CHECKPOINT
	typedef typename Host<TContainer>::Type THost;
	typedef typename Value<THost>::Type TSmallValue;
	
	TSmallValue flagValue = _getValueExpandFlagValue(container(me));
	if (_value >= flagValue)
	{//use map to store LargeValue
		assignValue(hostIterator(me), flagValue);
		_getMap(container(me))[position(me)] = _value;
	}
	else
	{
		assignValue(hostIterator(me), _value);
	}
}

//____________________________________________________________________________

template <typename TContainer, typename TIterator, typename TValue>
inline void
moveValue(Iter<TContainer, AdaptorIterator<TIterator, ValueExpandIter> > & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	assignValue(me, _value);
}
template <typename TContainer, typename TIterator, typename TValue>
inline void
moveValue(Iter<TContainer, AdaptorIterator<TIterator, ValueExpandIter> > const & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	assignValue(me, _value);
}


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
