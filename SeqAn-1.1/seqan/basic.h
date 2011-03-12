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
  $Id: basic.h,v 1.2 2009/05/06 20:32:59 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_H
#define SEQAN_HEADER_BASIC_H

//____________________________________________________________________________
// prerequisites

#include <seqan/platform.h>

//#include <cstring>
#ifdef PLATFORM_WINDOWS
#include <limits>	// limits include file exists only for g++ >= 3.0
#endif

#include <cstddef>	// size_t
#include <cstdio>	// FILE, basic_debug
#include <ctime>
#include <iterator>
#include <algorithm>
#include <memory.h> // memset
#include <string>	// basic_profile

#define SEQAN_NAMESPACE_MAIN seqan

//____________________________________________________________________________

#include <seqan/basic/basic_forwards.h>
#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/basic/basic_generated_forwards.h>
#endif

#include <seqan/basic/basic_debug.h>
#include <seqan/basic/basic_profile.h>
#include <seqan/basic/basic_definition.h>
#include <seqan/basic/basic_metaprogramming.h>
#include <seqan/basic/basic_type.h>
#include <seqan/basic/basic_tag.h>

//____________________________________________________________________________
// allocators

#include <seqan/basic/basic_allocator_interface.h>
#include <seqan/basic/basic_allocator_to_std.h>

#include <seqan/basic/basic_holder.h>

#include <seqan/basic/basic_allocator_simple.h>
#include <seqan/basic/basic_allocator_singlepool.h>
#include <seqan/basic/basic_allocator_multipool.h>
//#include <seqan/basic/basic_allocator_chunkpool.h>

//____________________________________________________________________________

#include <seqan/basic/basic_converter.h>
#include <seqan/basic/basic_compare.h>
#include <seqan/basic/basic_operator.h>

#include <seqan/basic/basic_host.h>

//____________________________________________________________________________
// iterators

#include <seqan/basic/basic_iterator.h>
#include <seqan/basic/basic_iterator_base.h>

#include <seqan/basic/basic_transport.h>

#include <seqan/basic/basic_iterator_simple.h>
#include <seqan/basic/basic_iterator_adaptor.h>
#include <seqan/basic/basic_iterator_position.h>
#include <seqan/basic/basic_iterator_adapt_std.h>
//#include <seqan/basic_identifier.h>

#include <seqan/basic/basic_proxy.h>

#include <seqan/basic/basic_pointer.h>

//____________________________________________________________________________
// alphabets

#include <seqan/basic/basic_alphabet_interface.h>
#include <seqan/basic/basic_alphabet_trait_basic.h>

#include <seqan/basic/basic_alphabet_interface2.h>

#include <seqan/basic/basic_alphabet_simple_tabs.h>
#include <seqan/basic/basic_alphabet_simple.h>

//____________________________________________________________________________

//#include <seqan/basic/basic_counted_ptr>
#include <seqan/basic/basic_volatile_ptr.h>

#include <seqan/basic/basic_aggregates.h>

#endif //#ifndef SEQAN_HEADER_...
