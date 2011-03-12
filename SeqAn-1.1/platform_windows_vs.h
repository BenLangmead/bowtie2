#define PLATFORM "windows"

#ifndef PLATFORM_WINDOWS
  #define PLATFORM_WINDOWS
#endif

#pragma warning( disable : 4675 )
#pragma warning( disable : 4503 )

#define finline __forceinline

//define SEQAN_SWITCH_USE_FORWARDS to use generated forwards 
//#define SEQAN_SWITCH_USE_FORWARDS
