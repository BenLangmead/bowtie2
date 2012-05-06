/*
 * Copyright 2011, Ben Langmead <blangmea@jhsph.edu>
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

#ifndef OPTS_H_
#define OPTS_H_

enum {
	ARG_ORIG = 256,             // --orig
	ARG_SEED,                   // --seed
	ARG_SOLEXA_QUALS,           // --solexa-quals
	ARG_VERBOSE,                // --verbose
	ARG_STARTVERBOSE,           // --startverbose
	ARG_QUIET,                  // --quiet
	ARG_METRIC_IVAL,            // --met
	ARG_METRIC_FILE,            // --met-file
	ARG_METRIC_STDERR,          // --met-stderr
	ARG_METRIC_PER_READ,        // --met-per-read
	ARG_REFIDX,                 // --refidx
	ARG_SANITY,                 // --sanity
	ARG_PARTITION,              // --partition
	ARG_INTEGER_QUALS,          // --int-quals
	ARG_USE_SPINLOCK,           // --nospin
	ARG_FILEPAR,                // --filepar
	ARG_SHMEM,                  // --shmem
	ARG_MM,                     // --mm
	ARG_MMSWEEP,                // --mmsweep
	ARG_FF,                     // --ff
	ARG_FR,                     // --fr
	ARG_RF,                     // --rf
	ARG_NO_MIXED,               // --no-mixed
	ARG_NO_DISCORDANT,          // --no-discordant
	ARG_CACHE_LIM,              // --
	ARG_CACHE_SZ,               // --
	ARG_NO_FW,                  // --nofw
	ARG_NO_RC,                  // --norc
	ARG_SKIP,                   // --skip
	ARG_ONETWO,                 // --12
	ARG_PHRED64,                // --phred64
	ARG_PHRED33,                // --phred33
	ARG_HADOOPOUT,              // --hadoopout
	ARG_FUZZY,                  // --fuzzy
	ARG_FULLREF,                // --fullref
	ARG_USAGE,                  // --usage
	ARG_SNPPHRED,               // --snpphred
	ARG_SNPFRAC,                // --snpfrac
	ARG_SAM_NO_QNAME_TRUNC,     // --sam-no-qname-trunc
	ARG_SAM_OMIT_SEC_SEQ,       // --sam-omit-sec-seq
	ARG_SAM_NOHEAD,             // --sam-noHD/--sam-nohead
	ARG_SAM_NOSQ,               // --sam-nosq/--sam-noSQ
	ARG_SAM_RG,                 // --sam-rg
	ARG_SAM_RGID,               // --sam-rg-id
	ARG_GAP_BAR,                // --gbar
	ARG_QUALS1,                 // --Q1
	ARG_QUALS2,                 // --Q2
	ARG_QSEQ,                   // --qseq
	ARG_SEED_SUMM,              // --seed-summary
	ARG_OVERHANG,               // --overhang
	ARG_NO_CACHE,               // --no-cache
	ARG_USE_CACHE,              // --cache
	ARG_NOISY_HPOLY,            // --454/--ion-torrent
	ARG_LOCAL,                  // --local
	ARG_END_TO_END,             // --end-to-end
	ARG_SCAN_NARROWED,          // --scan-narrowed
	ARG_QC_FILTER,              // --qc-filter
	ARG_BWA_SW_LIKE,            // --bwa-sw-like
	ARG_MULTISEED_IVAL,         // --multiseed
	ARG_SCORE_MIN,              // --score-min
	ARG_SCORE_MA,               // --ma
	ARG_SCORE_MMP,              // --mm
	ARG_SCORE_NP,               // --nm
	ARG_SCORE_RDG,              // --rdg
	ARG_SCORE_RFG,              // --rfg
	ARG_N_CEIL,                 // --n-ceil
	ARG_DPAD,                   // --dpad
	ARG_SAM_PRINT_YI,           // --mapq-print-inputs
	ARG_ALIGN_POLICY,           // --policy
	ARG_PRESET_VERY_FAST,       // --very-fast
	ARG_PRESET_FAST,            // --fast
	ARG_PRESET_SENSITIVE,       // --sensitive
	ARG_PRESET_VERY_SENSITIVE,  // --very-sensitive
	ARG_PRESET_VERY_FAST_LOCAL,      // --very-fast-local
	ARG_PRESET_FAST_LOCAL,           // --fast-local
	ARG_PRESET_SENSITIVE_LOCAL,      // --sensitive-local
	ARG_PRESET_VERY_SENSITIVE_LOCAL, // --very-sensitive-local
	ARG_NO_SCORE_PRIORITY,      // --no-score-priority
	ARG_IGNORE_QUALS,           // --ignore-quals
	ARG_DESC,                   // --arg-desc
	ARG_TAB5,                   // --tab5
	ARG_TAB6,                   // --tab6
	ARG_WRAPPER,                // --wrapper
	ARG_DOVETAIL,               // --dovetail
	ARG_NO_DOVETAIL,            // --no-dovetail
	ARG_CONTAIN,                // --contain
	ARG_NO_CONTAIN,             // --no-contain
	ARG_OVERLAP,                // --overlap
	ARG_NO_OVERLAP,             // --no-overlap
	ARG_MAPQ_V,                 // --mapq-v
	ARG_SSE8,                   // --sse8
	ARG_SSE8_NO,                // --no-sse8
	ARG_UNGAPPED,               // --ungapped
	ARG_UNGAPPED_NO,            // --no-ungapped
	ARG_TIGHTEN,                // --tighten
	ARG_UNGAP_THRESH,           // --ungap-thresh
	ARG_EXACT_UPFRONT,          // --exact-upfront
	ARG_1MM_UPFRONT,            // --1mm-upfront
	ARG_EXACT_UPFRONT_NO,       // --no-exact-upfront
	ARG_1MM_UPFRONT_NO,         // --no-1mm-upfront
	ARG_1MM_MINLEN,             // --1mm-minlen
	ARG_VERSION,                // --version
	ARG_SEED_OFF,               // --seed-off
	ARG_SEED_BOOST_THRESH,      // --seed-boost
	ARG_READ_TIMES,             // --read-times
	ARG_EXTEND_ITERS,           // --extends
	ARG_DP_MATE_STREAK_THRESH,  // --db-mate-streak
	ARG_DP_FAIL_STREAK_THRESH,  // --dp-fail-streak
	ARG_UG_FAIL_STREAK_THRESH,  // --ug-fail-streak
	ARG_EE_FAIL_STREAK_THRESH,  // --ee-fail-streak
	ARG_DP_FAIL_THRESH,         // --dp-fails
	ARG_UG_FAIL_THRESH,         // --ug-fails
	ARG_MAPQ_EX,                // --mapq-extra
	ARG_NO_EXTEND,              // --no-extend
	ARG_REORDER,                // --reorder
	ARG_SHOW_RAND_SEED,         // --show-rand-seed
	ARG_READ_PASSTHRU,          // --passthrough
	ARG_SAMPLE,                 // --sample
	ARG_CP_MIN,                 // --cp-min
	ARG_CP_IVAL,                // --cp-ival
	ARG_TRI,                    // --tri
	ARG_LOCAL_SEED_CACHE_SZ,    // --local-seed-cache-sz
	ARG_CURRENT_SEED_CACHE_SZ,  // --seed-cache-sz
	ARG_SAM_NO_UNAL             // --no-unal
};

#endif

