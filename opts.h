/**
 * opts.h
 */

enum {
	ARG_ORIG = 256,
	ARG_SEED,
	ARG_DUMP_PATS,
	ARG_SOLEXA_QUALS,
	ARG_VERBOSE,
	ARG_STARTVERBOSE,
	ARG_QUIET,
	ARG_FAST,
	ARG_METRIC_IVAL,
	ARG_METRIC_FILE,
	ARG_METRIC_STDERR,
	ARG_REFIDX,
	ARG_SANITY,
	ARG_PARTITION,
	ARG_INTEGER_QUALS,
	ARG_NOMAQROUND,
	ARG_USE_SPINLOCK,
	ARG_FILEPAR,
	ARG_SHMEM,
	ARG_MM,
	ARG_MMSWEEP,
	ARG_FF,
	ARG_FR,
	ARG_RF,
	ARG_NO_MIXED,
	ARG_NO_DISCORDANT,
	ARG_CACHE_LIM,
	ARG_CACHE_SZ,
	ARG_NO_FW,
	ARG_NO_RC,
	ARG_SKIP,
	ARG_ONETWO,
	ARG_PHRED64,
	ARG_PHRED33,
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
	ARG_SAM_RG,                 // --sam-RG
	ARG_SUPPRESS_FIELDS,        // --suppress
	ARG_COLOR_SEQ,              // --col-cseq
	ARG_COLOR_EDIT,             // --col-cedit
	ARG_COLOR_QUAL,             // --col-cqual
	ARG_PRINT_PLACEHOLDERS,     // --print-placeholders
	ARG_PRINT_FLAGS,            // --print-flags
	ARG_PRINT_PARAMS,           // --print-params
	ARG_COLOR_KEEP_ENDS,        // --col-keepends
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
	ARG_SCAN_NARROWED,          // --scan-narrowed
	ARG_NO_SSE,                 // --no-sse
	ARG_QC_FILTER,              // --qc-filter
	ARG_BWA_SW_LIKE,            // --bwa-sw-like
	ARG_MULTISEED_IVAL,         // --multiseed
	ARG_MULTISEED_IVAL_CONST,   // --multiseed-const
	ARG_MULTISEED_IVAL_LINEAR,  // --multiseed-linear
	ARG_MULTISEED_IVAL_SQRT,    // --multiseed-sqrt
	ARG_MULTISEED_IVAL_LOG,     // --multiseed-log
	ARG_SCORE_MIN_LINEAR,       // --score-min
	ARG_SCORE_MIN_CONST,        // --score-min-const
	ARG_SCORE_MIN_SQRT,         // --score-min-sqrt
	ARG_SCORE_MIN_LOG,          // --score-min-log
	ARG_SCORES,                 // --scoring
	ARG_N_CEIL,                 // --n-ceil
	ARG_DPAD,                   // --dpad
	ARG_MAPQ_TOP_COEFF,         // --mapq-top-coeff
	ARG_MAPQ_BOT_COEFF,         // --mapq-bot-coeff
	ARG_MAPQ_MAX,               // --mapq-max
	ARG_SAM_PRINT_YI,           // --mapq-print-inputs
	ARG_ALIGN_POLICY,           // --policy
	ARG_PRESET_VERY_FAST,       // --very-fast
	ARG_PRESET_FAST,            // --fast
	ARG_PRESET_SENSITIVE,       // --sensitive
	ARG_PRESET_VERY_SENSITIVE,  // --very-sensitive
	ARG_NO_SCORE_PRIORITY,      // --no-score-priority
	ARG_IGNORE_QUALS            // --ignore-quals
};
