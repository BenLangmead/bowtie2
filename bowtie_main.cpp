/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include "tokenize.h"
#include "ds.h"

#ifdef ENABLE_x86_64_v3
#include <unistd.h>
#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <cpuid.h>
#endif
#endif

using namespace std;

extern "C" {
	int bowtie(int argc, const char **argv);
}

#ifdef ENABLE_x86_64_v3

/*
 * Feature bits making up the x86-64-v3 microarchitecture level, as listed in
 * the Intel SDM Vol. 2A under CPUID.  The -v256 binary is compiled with
 * -march=x86-64-v3, so it may emit any instruction in that level and not just
 * AVX2; testing a subset would let it launch on a CPU that then dies with
 * SIGILL on the first BMI2 or FMA instruction.
 *
 * AVX2, BMI1 and BMI2 are reported in leaf 7 sub-leaf 0, and LZCNT in
 * extended leaf 0x80000001.  Only the SSE levels are in leaf 1.
 */

/* CPUID.(EAX=1, ECX=0):ECX */
#define V3_L1_ECX_SSE3    (1u <<  0)
#define V3_L1_ECX_SSSE3   (1u <<  9)
#define V3_L1_ECX_FMA     (1u << 12)
#define V3_L1_ECX_CX16    (1u << 13)
#define V3_L1_ECX_SSE4_1  (1u << 19)
#define V3_L1_ECX_SSE4_2  (1u << 20)
#define V3_L1_ECX_MOVBE   (1u << 22)
#define V3_L1_ECX_POPCNT  (1u << 23)
#define V3_L1_ECX_OSXSAVE (1u << 27)
#define V3_L1_ECX_AVX     (1u << 28)
#define V3_L1_ECX_F16C    (1u << 29)

/* CPUID.(EAX=1, ECX=0):EDX */
#define V3_L1_EDX_SSE     (1u << 25)
#define V3_L1_EDX_SSE2    (1u << 26)

/* CPUID.(EAX=7, ECX=0):EBX */
#define V3_L7_EBX_BMI1    (1u <<  3)
#define V3_L7_EBX_AVX2    (1u <<  5)
#define V3_L7_EBX_BMI2    (1u <<  8)

/* CPUID.(EAX=0x80000001):ECX */
#define V3_E1_ECX_LZCNT   (1u <<  5)

/* XCR0 bits that have to be set for the OS to preserve YMM state */
#define V3_XCR0_YMM       0x6

/*
 * Read one CPUID leaf into regs[] = {EAX, EBX, ECX, EDX}.  Returns 0 if the
 * CPU does not implement the leaf.  Masking the leaf with 0x80000000 picks
 * the basic or the extended range, which is what the maximum-leaf query
 * expects.
 */
static int cpuid_leaf(unsigned int leaf, unsigned int subleaf, unsigned int regs[4]) {
#if defined(_MSC_VER)
	int r[4];
	__cpuid(r, (int)(leaf & 0x80000000u));
	if((unsigned int)r[0] < leaf) return 0;
	__cpuidex(r, (int)leaf, (int)subleaf);
	regs[0] = (unsigned int)r[0]; regs[1] = (unsigned int)r[1];
	regs[2] = (unsigned int)r[2]; regs[3] = (unsigned int)r[3];
	return 1;
#else
	if(__get_cpuid_max(leaf & 0x80000000u, NULL) < leaf) return 0;
	__cpuid_count(leaf, subleaf, regs[0], regs[1], regs[2], regs[3]);
	return 1;
#endif
}

/* Only valid once CPUID has reported OSXSAVE. */
static unsigned long long read_xcr0(void) {
#if defined(_MSC_VER)
	return _xgetbv(0);
#else
	unsigned int lo, hi;
	__asm__ __volatile__("xgetbv" : "=a" (lo), "=d" (hi) : "c" (0));
	return ((unsigned long long)hi << 32) | lo;
#endif
}

static int has_x86_64_v3(void) {
	unsigned int regs[4];
	const unsigned int need_l1_ecx =
		V3_L1_ECX_SSE3 | V3_L1_ECX_SSSE3 | V3_L1_ECX_FMA | V3_L1_ECX_CX16 |
		V3_L1_ECX_SSE4_1 | V3_L1_ECX_SSE4_2 | V3_L1_ECX_MOVBE |
		V3_L1_ECX_POPCNT | V3_L1_ECX_OSXSAVE | V3_L1_ECX_AVX | V3_L1_ECX_F16C;
	const unsigned int need_l1_edx = V3_L1_EDX_SSE | V3_L1_EDX_SSE2;
	const unsigned int need_l7_ebx = V3_L7_EBX_BMI1 | V3_L7_EBX_AVX2 | V3_L7_EBX_BMI2;

	if(!cpuid_leaf(1, 0, regs)) return 0;
	if((regs[2] & need_l1_ecx) != need_l1_ecx) return 0;
	if((regs[3] & need_l1_edx) != need_l1_edx) return 0;

	// The CPU can have AVX and the OS still not preserve the upper halves of
	// the YMM registers across a context switch, in which case any AVX
	// instruction faults.  OSXSAVE above says XGETBV is safe to execute.
	if((read_xcr0() & V3_XCR0_YMM) != V3_XCR0_YMM) return 0;

	if(!cpuid_leaf(7, 0, regs)) return 0;
	if((regs[1] & need_l7_ebx) != need_l7_ebx) return 0;

	if(!cpuid_leaf(0x80000001u, 0, regs)) return 0;
	if(!(regs[2] & V3_E1_ECX_LZCNT)) return 0;

	return 1;
}

void check_x86_64_v3(int argc, const char **argv) {
	if (has_x86_64_v3() && (argc<126)) {
		const char* new_argv[128]; // should always be enough, but above check enforces it, too
		const char * org_path = argv[0];
		// Append -v256 to the original path 
		const int fn_len = strlen(org_path);
		char *new_path = (char*) malloc(fn_len+16);
		memcpy(new_path,org_path,fn_len);
		strncpy(new_path+fn_len,"-v256",15);

		for (int i=1; i<=argc; i++) new_argv[i] = argv[i]; // all but first the same, also copy final NULL
		new_argv[0] = new_path;
		// now replace the executable with new variant
		// assuming the executable exists... execvp will gracefully fail else, which is OK
		execvp(new_argv[0], (char *const *)new_argv);
		// we should never get out of the above call
		fprintf(stderr,"[WARNING] Failed to launch x86-64-v3 version, staying with default\n");
	}
}
#endif

/**
 * Bowtie main function.  It is placed in a separate source file to
 * make it slightly easier to compile Bowtie as a library.
 *
 * If the user specifies -A <file> as the first two arguments, main
 * will interpret that file as having one set of command-line arguments
 * per line, and will dispatch each batch of arguments one at a time to
 * bowtie.
 */
int main(int argc, const char **argv) {
#ifdef ENABLE_x86_64_v3
	check_x86_64_v3(argc, argv);
#endif
	if(argc > 2 && strcmp(argv[1], "-A") == 0) {
		const char *file = argv[2];
		ifstream in;
			in.open(file);
		char buf[4096];
		int lastret = -1;
		while(in.getline(buf, 4095)) {
			EList<string> args;
			args.push_back(string(argv[0]));
			tokenize(buf, " \t", args);
			const char **myargs = (const char**)malloc(sizeof(char*)*args.size());
			for(size_t i = 0; i < args.size(); i++) {
				myargs[i] = args[i].c_str();
			}
			if(args.size() == 1) continue;
			lastret = bowtie((int)args.size(), myargs);
			free(myargs);
		}
		if(lastret == -1) {
			cerr << "Warning: No arg strings parsed from " << file << endl;
			return 0;
		}
		return lastret;
	} else {
		return bowtie(argc, argv);
	}
}
