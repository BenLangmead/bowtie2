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
  $Id: basic_aggregates.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_AGGREGATES_H
#define SEQAN_HEADER_BASIC_AGGREGATES_H

namespace SEQAN_NAMESPACE_MAIN
{

//____________________________________________________________________________

    struct _Compressed;
	typedef Tag<_Compressed> Compressed;

	// for Pairs with small i1-values
	// store i1 and i2 in one word of type i2
	// use the upper bits for i1 and the lower bits for i2
	template <unsigned valueSizeI1 = 16>
	struct CutCompressed {
		enum { bitSizeI1 = Log2<valueSizeI1>::VALUE };
	};

/**
.Class.Pair:
..cat:Aggregates
..summary:Stores two arbitrary objects.
..signature:Pair<T1, T2[, Compression]>
..param.T1:The type of the first object.
..param.T2:The type of the second object.
..param.Compression:If $Compressed$, the pair is stored in a more space efficient way (useful for external storage).
...note:When compression is enabled, referring to members is not allowed.
...default:$void$, no compression (faster access).
.Memfunc.Pair#Pair:
..class:Class.Pair
..summary:Constructor
..signature:Pair<T1, T2> ()	
..signature:Pair<T1, T2> (pair)
..signature:Pair<T1, T2> (i1, i2)
..param.pair:Other Pair object. (copy constructor)
..param.i1:T1 object.
..param.i2:T2 object.
.Memvar.Pair#i1:
..class:Class.Pair
..summary:T1 object
.Memvar.Pair#i2:
..class:Class.Pair
..summary:T2 object
*/

	// standard storage 
	template <typename _T1, typename _T2 = _T1, typename TCompression = void>
    struct Pair {
        typedef _T1 T1;
        typedef _T2 T2;
	    _T1 i1;
	    _T2 i2;
		inline Pair() {}
		inline Pair(Pair const &_p): i1(_p.i1), i2(_p.i2) {}
		inline Pair(_T1 const &_i1, _T2 const &_i2): i1(_i1), i2(_i2) {}

		template <typename __T1, typename __T2, typename __TCompression>
		inline Pair(Pair<__T1, __T2, __TCompression> const &_p):
			i1(getValueI1(_p)), i2(getValueI2(_p)) {}
    };



	// unaligned and unpadded storage (space efficient)
#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif
    template <typename _T1, typename _T2>
    struct Pair<_T1, _T2, Compressed> {
        typedef _T1 T1;
        typedef _T2 T2;
        _T1 i1;
        _T2 i2;
		inline Pair() {}
		inline Pair(Pair const &_p): i1(_p.i1), i2(_p.i2) {}
		inline Pair(_T1 const &_i1, _T2 const &_i2): i1(_i1), i2(_i2) {}

		template <typename __T1, typename __T2, typename __TCompression>
		inline Pair(Pair<__T1, __T2, __TCompression> const &_p):
			i1(getValueI1(_p)), i2(getValueI2(_p)) {}
	}
#ifndef PLATFORM_WINDOWS
	__attribute__((packed))
#endif
	;
#ifdef PLATFORM_WINDOWS
    #pragma pack(pop)
#endif



#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif
    template <typename _T1, typename _T2, unsigned valueSizeI1>
    struct Pair<_T1, _T2, CutCompressed<valueSizeI1> > {
        typedef _T1 T1;
        typedef _T2 T2;

		typedef _T2 T12;

        T12 i12;

		enum { bitSizeI1 = CutCompressed<valueSizeI1>::bitSizeI1 };
        enum { bitShiftI1 = BitsPerValue<T12>::VALUE - bitSizeI1 };

		inline Pair() {}
		inline Pair(Pair const &_p): i12(_p.i12) {}
		inline Pair(_T1 const &_i1, _T2 const &_i2):
			i12(((T12)_i1 << bitShiftI1) + (T12)_i2) {}

		template <typename __T1, typename __T2, typename __TCompression>
		inline Pair(Pair<__T1, __T2, __TCompression> const &_p):
			i12(((T12)getValueI1(_p) << bitShiftI1) + (T12)getValueI2(_p)) {}
	}
#ifndef PLATFORM_WINDOWS
	__attribute__((packed))
#endif
	;
#ifdef PLATFORM_WINDOWS
    #pragma pack(pop)
#endif



    template <typename _T1, typename _T2, typename TCompression>
	std::ostream& operator<<(std::ostream &out, Pair<_T1,_T2,TCompression> const &p) {
		out << "< " << getValueI1(p) << " , " << getValueI2(p) << " >";
		return out;
	}

	template <typename T1, typename T2, typename TCompression>
	struct Value< Pair<T1, T2, TCompression>, 1 > {
		typedef T1 Type;
	};

	template <typename T1, typename T2, typename TCompression>
	struct Value< Pair<T1, T2, TCompression>, 2 > {
		typedef T2 Type;
	};

	template <typename T1, typename T2, typename TCompression>
	struct Spec< Pair<T1, T2, TCompression> > {
		typedef TCompression Type;
	};


//____________________________________________________________________________

	template <typename TKey, typename TObject, typename TSpec>
	struct Key< Pair<TKey, TObject, TSpec> > 
	{
		typedef TKey Type;
	};

	template <typename TKey, typename TCargo, typename TSpec>
	struct Cargo< Pair<TKey, TCargo, TSpec> > 
	{
		typedef TCargo Type;
	};
//____________________________________________________________________________

/**
.Class.Triple:
..cat:Aggregates
..summary:Stores three arbitrary objects.
..signature:Triple<T1, T2, T3[, Compression]>
..param.T1:The type of the first object.
..param.T2:The type of the second object.
..param.T3:The type of the third object.
..param.Compression:If $Compressed$, the triple is stored in a more space efficient way (useful for external storage).
...note:When compression is enabled, referring to members is not allowed.
...default:$void$, no compression (faster access).
.Memfunc.Triple#Triple:
..class:Class.Triple
..summary:Constructor
..signature:Triple<T1, T2, T3> ()
..signature:Triple<T1, T2, T3> (triple)
..signature:Triple<T1, T2, T3> (i1, i2, i3)
..param.triple:Other Triple object. (copy constructor)
..param.i1:T1 object.
..param.i2:T2 object.
..param.i3:T3 object.
.Memvar.Triple#i1:
..class:Class.Triple
..summary:T1 object
.Memvar.Triple#i2:
..class:Class.Triple
..summary:T2 object
.Memvar.Triple#i3:
..class:Class.Triple
..summary:T3 object
*/

	// standard storage 
	template <typename _T1, typename _T2 = _T1, typename _T3 = _T1, typename TCompression = void>
    struct Triple {
        typedef _T1 T1;
        typedef _T2 T2;
        typedef _T3 T3;
        _T1 i1;
        _T2 i2;
        _T3 i3;
		inline Triple() {}
		inline Triple(Triple const &_p):
			i1(_p.i1), i2(_p.i2), i3(_p.i3) {}
		inline Triple(_T1 const &_i1, _T2 const &_i2, _T3 const &_i3):
			i1(_i1), i2(_i2), i3(_i3) {}

		template <typename __T1, typename __T2, typename __T3, typename __TCompression>
		inline Triple(Triple<__T1, __T2, __T3, __TCompression> const &_p):
			i1(getValueI1(_p)), i2(getValueI2(_p)), i3(getValueI3(_p)) {}
	};

	// unaligned and unpadded storage (space efficient)
#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif
    template <typename _T1, typename _T2, typename _T3>
    struct Triple<_T1, _T2, _T3, Compressed> {
        typedef _T1 T1;
        typedef _T2 T2;
        typedef _T3 T3;
        _T1 i1;
        _T2 i2;
        _T3 i3;
		inline Triple() {}
		inline Triple(Triple const &_p):
			i1(_p.i1), i2(_p.i2), i3(_p.i3) {}
		inline Triple(_T1 const &_i1, _T2 const &_i2, _T3 const &_i3):
			i1(_i1), i2(_i2), i3(_i3) {}

		template <typename __T1, typename __T2, typename __T3, typename __TCompression>
		inline Triple(Triple<__T1, __T2, __T3, __TCompression> const &_p):
			i1(getValueI1(_p)), i2(getValueI2(_p)), i3(getValueI3(_p)) {}
	}
#ifndef PLATFORM_WINDOWS
	__attribute__((packed))
#endif
	;
#ifdef PLATFORM_WINDOWS
    #pragma pack(pop)
#endif

	template <typename _T1, typename _T2, typename _T3, typename TCompression>
	std::ostream& operator<<(std::ostream &out, Triple<_T1,_T2,_T3,TCompression> const &t) {
		out << "< " << getValueI1(t) << " , " << getValueI2(t) << " , " << getValueI3(t) << " >";
		return out;
	}

	template <typename T1, typename T2, typename T3, typename TCompression>
	struct Value< Triple<T1, T2, T3, TCompression>, 1 > {
		typedef T1 Type;
	};

	template <typename T1, typename T2, typename T3, typename TCompression>
	struct Value< Triple<T1, T2, T3, TCompression>, 2 > {
		typedef T2 Type;
	};

	template <typename T1, typename T2, typename T3, typename TCompression>
	struct Value< Triple<T1, T2, T3, TCompression>, 3 > {
		typedef T3 Type;
	};

	template <typename T1, typename T2, typename T3, typename TCompression>
	struct Spec< Triple<T1, T2, T3, TCompression> > {
		typedef TCompression Type;
	};


//____________________________________________________________________________

/**
.Class.Tuple:
..cat:Aggregates
..summary:A plain fixed-length string.
..signature:Tuple<T, size[, compress]>
..param.T:The value type, that is the type of characters stored in the tuple.
..param.size:The size/length of the tuple.
...remarks:In contrast to @Class.String@ the length of Tuple is fixed.
..param.compress:Enable/Disable compression.
..param.compress:If $void$, no compression is used.
..param.compress:If $Compressed$, the characters are stored as a bit sequence in an ordinal type (char, ..., __int64)
...remarks:Only useful for small alphabets and small tuple sizes (|Sigma|^size <= 2^64) as for DNA or protein m-grams)
...default:void.
..see:Spec.Sampler
*/

	// standard storage 
	template <typename _T, unsigned _size, typename TCompression = void>
    struct Tuple {
        typedef _T T;
        enum { size = _size };
        _T i[_size];

		template <typename TPos>
        inline _T& operator[](TPos k) {
            SEQAN_ASSERT(k >= 0 && k < size);
            return i[k];
        }
		template <typename TPos>
        inline const _T& operator[](TPos k) const {
            SEQAN_ASSERT(k >= 0 && k < size);
            return i[k];
        }
		inline _T* operator&() { return i; }
		inline const _T* operator&() const { return i; }

		// has to be inline because elements (like this tuple) of packed structs can't be arguments
		template <typename TPos, typename SSS>
		inline SSS const assignValueAt(TPos k, SSS const source) {
			return i[k] = source;
		}
    };


    template < unsigned char _size >
	struct _BitVector {
        typedef typename _BitVector<_size + 1>::Type Type;
    };

    template <> struct _BitVector<8> { typedef unsigned char Type; };
    template <> struct _BitVector<16> { typedef unsigned short Type; };
    template <> struct _BitVector<32> { typedef unsigned int Type; };
    template <> struct _BitVector<64> { typedef __int64 Type; };
    template <> struct _BitVector<255> { typedef __int64 Type; };

	// bit-compressed storage (space efficient)
#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif
    template <typename _T, unsigned _size>
    struct Tuple<_T, _size, Compressed> {
        typedef _T T;
        enum { size = _size };
        enum { bitSize = BitsPerValue<_T>::VALUE };
        enum { bitMask = (1 << bitSize) - 1 };
        enum { mask = (1 << (size * bitSize)) - 1 };
        typedef typename _BitVector< bitSize * size >::Type CT;
        
        CT i;
/*
		inline Tuple() {
			SEQAN_ASSERT(bitSize * size <= sizeof(CT) * 8);
		}
*/
		template <typename TPos>
        inline const _T operator[](TPos k) const {
            SEQAN_ASSERT(k >= 0 && k < size);
            return (i >> (size - 1 - k) * bitSize) & bitMask;
        }
		template <unsigned __size>
		inline Tuple operator=(Tuple<_T, __size, Compressed> const &_right) {
			i = _right.i;
			return *this;
		}
		template <typename TShiftSize>
        inline CT operator<<=(TShiftSize shift) {
            return i = (i << (shift * bitSize)) & mask;
        }
		template <typename TShiftSize>
        inline CT operator<<(TShiftSize shift) const {
            return (i << (shift * bitSize)) & mask;
        }
		template <typename TShiftSize>
        inline CT operator>>=(TShiftSize shift) {
            return i = (i >> (shift * bitSize));
        }
		template <typename TShiftSize>
        inline CT operator>>(TShiftSize shift) const {
            return i >> (shift * bitSize);
        }
        template <typename T>
        inline void operator|=(T const &t) {
            i |= t;
        }
        template <typename T, typename TSpec>
        inline void operator|=(SimpleType<T, TSpec> const &t) {
            i |= t.value;
        }
		inline CT* operator&() { return &i; }
		inline const CT* operator&() const { return &i; }

		// has to be inline because elements (like this tuple) of packed structs can't be arguments
		template <typename TPos, typename SSS>
		inline SSS const assignValueAt(TPos k, SSS const source) {
			typedef Tuple<_T, _size, Compressed> Tup;
			typename Tup::CT mask = Tup::bitMask << ((_size - 1 - k) * bitSize);
			i = (i & ~mask) | ((CT)source << ((_size - 1 - k) * bitSize));
			return source;
		}
    }
#ifndef PLATFORM_WINDOWS
	__attribute__((packed))
#endif
	;
#ifdef PLATFORM_WINDOWS
    #pragma pack(pop)
#endif


//////////////////////////////////////////////////////////////////////////////
// length

    template <typename _T, unsigned _size, typename TCompression>
	inline unsigned length(Tuple<_T, _size, TCompression> const &) { return _size; }

	///.Metafunction.LENGTH.param.T.type:Class.Tuple
    template <typename _T, unsigned _size, typename TCompression>
	struct LENGTH< Tuple<_T, _size, TCompression> >
	{
		enum { VALUE = _size };
	};

//////////////////////////////////////////////////////////////////////////////
// assignValueAt

    template <typename TObject, typename TPos, typename TSource>
    inline TSource & 
	assignValueAt(TObject &me, TPos k, TSource &source) {
        assign(value(me, k), source);
		return source;
    }

    template <typename TObject, typename TPos, typename TSource>
    inline TSource const & 
	assignValueAt(TObject &me, TPos k, TSource const &source) {
        assign(value(me, k), source);
		return source;
    }

    template <typename TTT, unsigned _size, typename SSS, typename TPos>
    inline SSS const assignValueAt(Tuple<TTT, _size, void> &me, TPos k, SSS const source) {
        return me.i[k] = source;
    }

    template <typename TTT, unsigned _size, typename SSS, typename TPos>
    inline SSS const assignValueAt(Tuple<TTT, _size, Compressed> &me, TPos k, SSS const source) {
        typedef Tuple<TTT, _size, Compressed> Tup;
        typename Tup::CT mask = Tup::bitMask << ((_size - 1 - k) * me.bitSize);
        me.i = (me.i & ~mask) | source << ((_size - 1 - k) * me.bitSize);
        return source;
    }

    template <typename TTT, typename SSS, typename SSSpec, unsigned _size, typename TPos>
    inline SimpleType<SSS, SSSpec> const & assignValueAt(Tuple<TTT, _size, Compressed> &me, TPos k, SimpleType<SSS, SSSpec> const &source) {
        typedef Tuple<TTT, _size, Compressed> Tup;
        typename Tup::CT mask = Tup::bitMask << ((_size - 1 - k) * me.bitSize);
        me.i = (me.i & ~mask) | source.value << ((_size - 1 - k) * me.bitSize);
        return source;
    }

//////////////////////////////////////////////////////////////////////////////
// clear

	template <typename TTT, unsigned _size, typename TCompression>
	inline void clear(Tuple<TTT, _size, TCompression> &me) {
        memset<sizeof(me.i), 0>(&(me.i));
	}
    template <typename TTT, unsigned _size>
	inline void clear(Tuple<TTT, _size, Compressed> &me) {
		me.i = 0; 
	}

//////////////////////////////////////////////////////////////////////////////
// optimized compares

	template <typename TTT, unsigned _sizeL, unsigned _sizeR>
	inline bool operator<(Tuple<TTT, _sizeL, Compressed> const &_left, Tuple<TTT, _sizeR, Compressed> const &_right) {
		return _left.i < _right.i;
	}
	template <typename TTT, unsigned _sizeL, unsigned _sizeR>
	inline bool operator>(Tuple<TTT, _sizeL, Compressed> const &_left, Tuple<TTT, _sizeR, Compressed> const &_right) {
		return _left.i > _right.i;
	}
	template <typename TTT, unsigned _sizeL, unsigned _sizeR>
	inline bool operator==(Tuple<TTT, _sizeL, Compressed> const &_left, Tuple<TTT, _sizeR, Compressed> const &_right) {
		return _left.i == _right.i;
	}
	template <typename TTT, unsigned _sizeL, unsigned _sizeR>
	inline bool operator!=(Tuple<TTT, _sizeL, Compressed> const &_left, Tuple<TTT, _sizeR, Compressed> const &_right) {
		return _left.i != _right.i;
	}

//////////////////////////////////////////////////////////////////////////////
// optimized shifts

    struct _TupleShiftLeftWorker {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg[I-1] = arg[I];
        }
    };

    struct _TupleShiftRightWorker {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg[I] = arg[I-1];
        }
    };

	template <typename _T, unsigned _size, typename TCompression>
	inline void shiftLeft(Tuple<_T, _size, TCompression> &me) {
		LOOP<_TupleShiftLeftWorker, _size - 1>::run(me);
	}

	template <typename _T, unsigned _size, typename TCompression>
	inline void shiftRight(Tuple<_T, _size, TCompression> &me) {
		LOOP_REVERSE<_TupleShiftRightWorker, _size - 1>::run(me);
	}

	template <typename _T, unsigned _size>
	inline void shiftLeft(Tuple<_T, _size, Compressed> &me) {
		me<<=1;
	}

	template <typename _T, unsigned _size>
	inline void shiftRight(Tuple<_T, _size, Compressed> &me) {
		me>>=1;
	}

//////////////////////////////////////////////////////////////////////////////
// standard output

	template <typename _T, unsigned _size, typename TCompression>
	std::ostream& operator<<(std::ostream& out, Tuple<_T,_size,TCompression> const &a) {
		out << "[";
		if (a.size > 0)
			out << a[0];
		for(unsigned j = 1; j < a.size; ++j)
			out << " " << a[j];
		out << "]";
		return out;
	}

	template <typename _T, unsigned _size, typename TCompression>
	struct Value< Tuple<_T, _size, TCompression> > {
		typedef _T Type;
	};

	template <typename _T, unsigned _size, typename TCompression>
	struct Spec< Tuple<_T, _size, TCompression> > {
		typedef TCompression Type;
	};

//////////////////////////////////////////////////////////////////////////////
// getValueIx

	template <typename T1, typename T2, typename TCompression>
	inline T1 getValueI1(Pair<T1, T2, TCompression> const &pair) {
		return pair.i1;
	}

	template <typename T1, typename T2, typename TCompression>
	inline T2 getValueI2(Pair<T1, T2, TCompression> const &pair) {
		return pair.i2;
	}

	template <typename T1, typename T2, unsigned valueSizeI1>
	inline T1 getValueI1(Pair<T1, T2, CutCompressed<valueSizeI1> > const &pair) {
		typedef Pair<T1, T2, CutCompressed<valueSizeI1> > TPair;
		return pair.i12 >> TPair::bitShiftI1;
	}

	template <typename T1, typename T2, unsigned valueSizeI1>
	inline T2 getValueI2(Pair<T1, T2, CutCompressed<valueSizeI1> > const &pair) {
		typedef Pair<T1, T2, CutCompressed<valueSizeI1> > TPair;		 
		return pair.i12 & (((typename TPair::T12)1 << TPair::bitShiftI1) - 1);
	}
//____________________________________________________________________________

	template <typename T1, typename T2, typename T3, typename TCompression>
	inline T1 getValueI1(Triple<T1, T2, T3, TCompression> const &triple) {
		return triple.i1;
	}

	template <typename T1, typename T2, typename T3, typename TCompression>
	inline T2 getValueI2(Triple<T1, T2, T3, TCompression> const &triple) {
		return triple.i2;
	}

	template <typename T1, typename T2, typename T3, typename TCompression>
	inline T3 getValueI3(Triple<T1, T2, T3, TCompression> const &triple) {
		return triple.i3;
	}

//////////////////////////////////////////////////////////////////////////////
// assignValueIx

	template <typename T1, typename T2, typename TCompression, typename T>
	inline void assignValueI1(Pair<T1, T2, TCompression> &pair, T const &_i) {
		pair.i1 = _i;
	}

	template <typename T1, typename T2, typename TCompression, typename T>
	inline void assignValueI2(Pair<T1, T2, TCompression> &pair, T const &_i) {
		pair.i2 = _i;
	}

	template <typename T1, typename T2, unsigned valueSizeI1, typename T>
	inline void assignValueI1(Pair<T1, T2, CutCompressed<valueSizeI1> > &pair, T const &_i) 
	{
		typedef Pair<T1, T2, CutCompressed<valueSizeI1> > TPair;
		pair.i12 = ((typename TPair::T12)_i << TPair::bitShiftI1) |
		           (pair.i12 & (((typename TPair::T12)1 << TPair::bitShiftI1) - 1));
	}

	template <typename T1, typename T2, unsigned valueSizeI1, typename T>
	inline void assignValueI2(Pair<T1, T2, CutCompressed<valueSizeI1> > &pair, T const &_i) {
		typedef Pair<T1, T2, CutCompressed<valueSizeI1> > TPair;
		pair.i12 = (pair.i12 & ~(((typename TPair::T12)1 << TPair::bitShiftI1) - 1)) | _i;
	}
//____________________________________________________________________________

	template <typename T1, typename T2, typename T3, typename TCompression, typename T>
	inline T const assignValueI1(Triple<T1, T2, T3, TCompression> &triple, T const &_i) {
		return triple.i1 = _i;
	}

	template <typename T1, typename T2, typename T3, typename TCompression, typename T>
	inline T const assignValueI2(Triple<T1, T2, T3, TCompression> &triple, T const &_i) {
		return triple.i2 = _i;
	}

	template <typename T1, typename T2, typename T3, typename TCompression, typename T>
	inline T const assignValueI3(Triple<T1, T2, T3, TCompression> &triple, T const &_i) {
		return triple.i3 = _i;
	}

//////////////////////////////////////////////////////////////////////////////
// operator ==/!= for pairs and triples

	template <typename L1, typename L2, typename LCompression, typename R1, typename R2, typename RCompression>
	inline bool operator==(Pair<L1, L2, LCompression> const &_left, Pair<R1, R2, RCompression> const &_right) {
		return _left.i1 == _right.i1 && _left.i2 == _right.i2;
	}
	template <typename L1, typename L2, typename LCompression, typename R1, typename R2, typename RCompression>
	inline bool operator!=(Pair<L1, L2, LCompression> const &_left, Pair<R1, R2, RCompression> const &_right) {
		return _left.i1 != _right.i1 || _left.i2 != _right.i2;
	}

	template <typename L1, typename L2, unsigned LSizeI1, typename R1, typename R2, unsigned RSizeI1>
	inline bool operator==(Pair<L1, L2, CutCompressed<LSizeI1> > const &_left, Pair<R1, R2, CutCompressed<RSizeI1> > const &_right) {
		return _left.i12 == _right.i12;
	}
	template <typename L1, typename L2, unsigned LSizeI1, typename R1, typename R2, unsigned RSizeI1>
	inline bool operator!=(Pair<L1, L2, CutCompressed<LSizeI1> > const &_left, Pair<R1, R2, CutCompressed<RSizeI1> > const &_right) {
		return _left.i12 != _right.i12;
	}
//____________________________________________________________________________

	template <
		typename L1, typename L2, typename L3, typename LCompression, 
		typename R1, typename R2, typename R3, typename RCompression>
	inline bool operator==(Triple<L1, L2, L3, LCompression> const &_left, Triple<R1, R2, R3, RCompression> const &_right) {
		return _left.i1 == _right.i1 && _left.i2 == _right.i2 && _left.i3 == _right.i3;
	}
	template <
		typename L1, typename L2, typename L3, typename LCompression, 
		typename R1, typename R2, typename R3, typename RCompression>
	inline bool operator!=(Triple<L1, L2, L3, LCompression> const &_left, Triple<R1, R2, R3, RCompression> const &_right) {
		return _left.i1 != _right.i1 || _left.i2 != _right.i2 || _left.i3 != _right.i3;
	}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
