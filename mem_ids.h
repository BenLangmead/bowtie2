/*
 *  mem_ids.h
 *
 * Constants that identify the various consumers of heap memory, for
 * tracking purposes.
 */

// For holding index data
#define EBWT_CAT  ((int) 1)
// For holding index-building data
#define EBWTB_CAT ((int) 2)
// For holding cache data
#define CA_CAT    ((int) 3)
// For holding group-walk-left bookkeeping data
#define GW_CAT    ((int) 4)
// For holding alignment bookkeeping data
#define AL_CAT    ((int) 5)
// For holding dynamic programming bookkeeping data
#define DP_CAT    ((int) 6)
// For holding alignment results and other hit objects
#define RES_CAT   ((int) 7)
#define MISC_CAT  ((int) 9)
#define DEBUG_CAT ((int)10)
