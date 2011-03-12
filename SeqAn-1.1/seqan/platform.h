#define PLATFORM "gcc"

#ifndef PLATFORM_GCC
  #define PLATFORM_GCC
#endif

// should be set before including anything
#ifndef _FILE_OFFSET_BITS
  #define _FILE_OFFSET_BITS 64
#endif

#ifndef _LARGEFILE_SOURCE
  #define _LARGEFILE_SOURCE
#endif

//#include <unistd.h>
#include <inttypes.h>

#define finline __inline__

// default 64bit type
#ifndef __int64
typedef int64_t __int64;
#endif

//define SEQAN_SWITCH_USE_FORWARDS to use generated forwards 
#define SEQAN_SWITCH_USE_FORWARDS
