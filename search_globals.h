/*
 * search_globals.h
 *
 *  Created on: Dec 5, 2009
 *      Author: Ben Langmead
 */

#ifndef SEARCH_GLOBALS_H_
#define SEARCH_GLOBALS_H_

#include <stdint.h>

// declared in ebwt_search.cpp
extern bool     gColor;
extern bool     gColorExEnds;
extern bool     gReportOverhangs;
extern bool     gColorSeq;
extern bool     gColorEdit;
extern bool     gColorQual;
extern int      gSnpPhred;
extern bool     gNoMaqRound;
extern bool     gStrandFix;
extern bool     gRangeMode;
extern int      gVerbose;
extern int      gQuiet;
extern bool     gNofw;
extern bool     gNorc;
extern bool     gMate1fw;
extern bool     gMate2fw;
extern int      gMinInsert;
extern int      gMaxInsert;
extern int      gTrim5;
extern int      gTrim3;
extern int      gGapBarrier;
extern int      gAllowRedundant;

#endif /* SEARCH_GLOBALS_H_ */
