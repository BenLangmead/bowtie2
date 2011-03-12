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
  $Id: basic_metaprogramming.h,v 1.2 2009/02/19 01:51:23 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_BASIC_METAPROGRAMMING_H
#define SEQAN_BASIC_METAPROGRAMMING_H

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// generic "if" (using meta-programming)
	// if Flag is true,  the resulting type is Type1
	// if Flag is false, the resulting type is Type2
	//////////////////////////////////////////////////////////////////////////////

	template <bool Flag,class Type1, class Type2>
	struct IF
	{
		typedef Type1 Type;
	};

	template <class Type1, class Type2>
	struct IF<false,Type1,Type2>
	{
		typedef Type2 Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// generic type comparison (using meta-programming)
	// if Type1 equals Type2,		VALUE is true
	// if Type1 differs from Type2, VALUE is false
	//////////////////////////////////////////////////////////////////////////////

	template <class Type1, class Type2>
	struct TYPECMP
	{
		typedef False Type;
		enum { VALUE = false };
	};

	template <class Type1>
	struct TYPECMP<Type1, Type1>
	{
		typedef True Type;
		enum { VALUE = true };
	};

	//////////////////////////////////////////////////////////////////////////////
	// generic "switch" (using meta-programming)
	//////////////////////////////////////////////////////////////////////////////

	const int DEFAULT = ~(~0u >> 1); // initialize with the smallest int

	struct NilCase {};

	template <int tag_,class Type_,class Next_ = NilCase>
	struct CASE
	{
		enum { tag = tag_ };
		typedef Type_ Type;
		typedef Next_ Next;
	};

	template <int tag,class Case>
	class SWITCH
	{
		typedef typename Case::Next NextCase;
		enum
		{
			caseTag = Case::tag,
			found   = (caseTag == tag || caseTag == DEFAULT)
		};
	public:
		typedef typename
			IF<
				found,
				typename Case::Type,
				typename SWITCH<tag,NextCase>::Type
			>::Type Type;
	};

	template <int tag>
	class SWITCH<tag,NilCase>
	{
	public:
		typedef NilCase Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// generic loops (using meta-programming)
	// corresponds to for(i=1; i<=I; ++i) ...
	//////////////////////////////////////////////////////////////////////////////

	// example of a loop Worker class
	struct WorkerNothing
	{
		template <typename Arg>
		static inline void body(Arg &arg, int I) {}
	};

	template <typename Worker, int I>
	class LOOP {
	public:
		template <typename Arg>
		static inline void run(Arg &arg) {
			LOOP<Worker, I - 1>::run(arg);
			Worker::body(arg, I);
		}
	};

	template <typename Worker>
	class LOOP<Worker, 0> {
	public:
		// end of loop
		template <typename Arg>
		static inline void run(Arg &) {}
	};

	//////////////////////////////////////////////////////////////////////////////
	// generic reverse loops (using meta-programming)
	// corresponds to for(i=I; i>0; --i) ...
	//////////////////////////////////////////////////////////////////////////////

	template <typename Worker, int I>
	class LOOP_REVERSE {
	public:
		template <typename Arg>
		static inline void run(Arg &arg) {
			Worker::body(arg, I);
			LOOP_REVERSE<Worker, I - 1>::run(arg);
		}
	};

	template <typename Worker>
	class LOOP_REVERSE<Worker, 0> {
	public:
		// end of loop
		template <typename Arg>
		static inline void run(Arg &) {}
	};

	//////////////////////////////////////////////////////////////////////////////
	// logarithmus dualis (using meta-programming)
	//////////////////////////////////////////////////////////////////////////////

	template < __int64 numerus >
	struct Log2 {
		enum { VALUE = Log2<(numerus + 1) / 2>::VALUE + 1 };		// ceil(log_2(n))
	};

	template < __int64 numerus >
	struct Log2Floor {
		enum { VALUE = Log2Floor<numerus / 2>::VALUE + 1 };		// floor(log_2(n))
	};

	template <> struct Log2<1> { enum { VALUE = 0 }; };
	template <> struct Log2<0> { enum { VALUE = 0 }; };
	template <> struct Log2Floor<1> { enum { VALUE = 0 }; };
	template <> struct Log2Floor<0> { enum { VALUE = 0 }; };


	//////////////////////////////////////////////////////////////////////////////
	// exponentiation (using meta-programming)
	//////////////////////////////////////////////////////////////////////////////

	template < __int64 base, __int64 exponent >
	struct Power {
		enum {
			VALUE =
				Power<base, exponent / 2>::VALUE *
				Power<base, exponent - (exponent / 2)>::VALUE
		};
	};

	template < __int64 base > struct Power<base, 1> { enum { VALUE = base }; };
	template < __int64 base > struct Power<base, 0> { enum { VALUE = 1 }; };


	//////////////////////////////////////////////////////////////////////////////
	// memset with fill size (using meta-programming)
	//////////////////////////////////////////////////////////////////////////////

	using ::memset;

	template <unsigned SIZE, bool direct>
	struct MemsetWorker {
		finline static void run(unsigned char* ptr, unsigned char c) { memset(ptr, c, SIZE); }
	};

	template <unsigned  SIZE>
	struct MemsetWorker<SIZE, true> {
		finline static void run(unsigned char* ptr, unsigned char c) {
			*((unsigned*)ptr) = ((unsigned)c << 24) + ((unsigned)c << 16) + ((unsigned)c << 8) + (unsigned)c;
			MemsetWorker<SIZE - 4, true>::run(ptr + 4, c);
		}
	};

	template <>
	struct MemsetWorker<0, true> {
		finline static void run(unsigned char*, unsigned char) {}
	};

	template <>
	struct MemsetWorker<1, true> {
		finline static void run(unsigned char* ptr, unsigned char c) { *ptr = c; }
	};

	template <>
	struct MemsetWorker<2, true> {
		finline static void run(unsigned char* ptr, unsigned char c) { *(unsigned short *)ptr = ((unsigned short)c << 8) + (unsigned short)c; }
	};

	template <>
	struct MemsetWorker<3, true> {
		finline static void run(unsigned char* ptr, unsigned char c) {
			MemsetWorker<2, true>::run(ptr, c);
			MemsetWorker<1, true>::run(ptr + 2, c);
		}
	};

	template <unsigned SIZE>
	finline void memset(void* ptr, unsigned char c) {
		MemsetWorker<SIZE, SIZE <= 32>::run((unsigned char*)ptr, c);
	}


	//////////////////////////////////////////////////////////////////////////////
	// memset with fill value (using meta-programming)
	//////////////////////////////////////////////////////////////////////////////

	template <unsigned SIZE, bool direct, unsigned char c>
	struct MemsetConstValueWorker {
		finline static void run(unsigned char* ptr) { memset(ptr, c, SIZE); }
	};

	template <unsigned  SIZE, unsigned char c>
	struct MemsetConstValueWorker<SIZE, true, c> {
		finline static void run(unsigned char* ptr) {
			*((unsigned*)ptr) = ((unsigned)c << 24) + ((unsigned)c << 16) + ((unsigned)c << 8) + (unsigned)c;
			MemsetConstValueWorker<SIZE - 4, true, c>::run(ptr + 4);
		}
	};

	template <unsigned char c>
	struct MemsetConstValueWorker<0, true, c> {
		finline static void run(unsigned char* ptr) {}
	};

	template <unsigned char c>
	struct MemsetConstValueWorker<1, true, c> {
		finline static void run(unsigned char* ptr) { *ptr = c; }
	};

	template <unsigned char c>
	struct MemsetConstValueWorker<2, true, c> {
		finline static void run(unsigned char* ptr) { *(unsigned short *)ptr = ((unsigned short)c << 8) + (unsigned short)c; }
	};

	template <unsigned char c>
	struct MemsetConstValueWorker<3, true, c> {
		finline static void run(unsigned char* ptr) {
			MemsetConstValueWorker<2, true, c>::run(ptr);
			MemsetConstValueWorker<1, true, c>::run(ptr + 2);
		}
	};

	template <unsigned SIZE, unsigned char c>
	finline void memset(void* ptr) {
		MemsetConstValueWorker<SIZE, SIZE <= 32, c>::run((unsigned char*)ptr);
	}

}

#endif
