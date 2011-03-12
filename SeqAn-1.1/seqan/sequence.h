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
  $Id: sequence.h,v 1.1 2008/08/25 16:20:06 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEQUENCE_H
#define SEQAN_HEADER_SEQUENCE_H

//____________________________________________________________________________
// prerequisites

#include <seqan/basic.h>

//____________________________________________________________________________

#include <seqan/sequence/sequence_forwards.h>

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/sequence/sequence_generated_forwards.h>
#endif

#include <seqan/sequence/sequence_interface.h>
#include <seqan/sequence/lexical.h>

//____________________________________________________________________________
// segments (suffix, ...)

#include <seqan/sequence/segment_base.h>
#include <seqan/sequence/segment_infix.h>
#include <seqan/sequence/segment_suffix.h>
#include <seqan/sequence/segment_prefix.h>

//____________________________________________________________________________
// strings

#include <seqan/sequence/string_base.h>
#include <seqan/sequence/string_pointer.h>
#include <seqan/sequence/string_alloc.h>
#include <seqan/sequence/string_array.h>
#include <seqan/sequence/string_cstyle.h>
#include <seqan/sequence/string_stack.h>
#include <seqan/sequence/string_packed.h>
#include <seqan/sequence/string_value_expand.h>

#include <seqan/sequence/std_string.h>

#include <seqan/sequence/sequence_multiple.h>
#include <seqan/sequence/sequence_shortcuts.h>

#endif //#ifndef SEQAN_HEADER_...
