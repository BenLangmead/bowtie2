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

#ifndef SAM_H_
#define SAM_H_

#include <string>
#include <vector>
#include "ds.h"
#include "read.h"
#include "util.h"
#include "aligner_result.h"
#include "scoring.h"

enum {
	// Comments use language from v1.4-r962 spec
	SAM_FLAG_PAIRED         = 1,   // templ. having mult. frag.s in sequencing
	SAM_FLAG_MAPPED_PAIRED  = 2,   // each frag properly aligned
	SAM_FLAG_UNMAPPED       = 4,   // fragment unmapped
	SAM_FLAG_MATE_UNMAPPED  = 8,   // next fragment in template unmapped
	SAM_FLAG_QUERY_STRAND   = 16,  // SEQ is reverse comp'ed from original
	SAM_FLAG_MATE_STRAND    = 32,  // next fragment SEQ reverse comp'ed
	SAM_FLAG_FIRST_IN_PAIR  = 64,  // first fragment in template
	SAM_FLAG_SECOND_IN_PAIR = 128, // last fragment in template
	SAM_FLAG_NOT_PRIMARY    = 256, // secondary alignment
	SAM_FLAG_FAILS_CHECKS   = 512, // not passing quality controls
	SAM_FLAG_DUPLICATE      = 1024 // PCR or optical duplicate
};

class AlnRes;
class AlnFlags;
class AlnSetSumm;

/**
 * Encapsulates all the various ways that a user may wish to customize SAM
 * output.
 */
class SamConfig {

	typedef EList<std::string> StrList;
	typedef EList<size_t> LenList;

public:

	SamConfig(
		const StrList& refnames,  // reference sequence names
		const LenList& reflens,   // reference sequence lengths
		bool truncQname,          // truncate read name to 255?
		bool appendComment,	  // append FASTA/Q comment to sam record
		bool omitsec,             // omit secondary SEQ/QUAL
		bool noUnal,              // omit unaligned reads
		const std::string& pg_id, // id
		const std::string& pg_pn, // name
		const std::string& pg_vn, // version
		const std::string& pg_cl, // command-line
		const std::string& rgs,   // read groups string
		bool print_as,
		bool print_xs,
		bool print_xss,
		bool print_yn,
		bool print_xn,
		bool print_x0,
		bool print_x1,
		bool print_xm,
		bool print_xo,
		bool print_xg,
		bool print_nm,
		bool print_md,
		bool print_yf,
		bool print_yi,
		bool print_ym,
		bool print_yp,
		bool print_yt,
		bool print_ys,
		bool print_zs,
		bool print_xr,
		bool print_xt,
		bool print_xd,
		bool print_xu,
		bool print_ye, // streak of failed DPs at end
		bool print_yl, // longest streak of failed DPs
		bool print_yu, // index of last succeeded DP
		bool print_xp, // print seed hit information
		bool print_yr, // # redundant seed hits
		bool print_zb, // # Ftab lookups
		bool print_zr, // # redundant path checks
		bool print_zf, // # FM Index ops
		bool print_zm, // FM Index op string for best-first search
		bool print_zi, // # seed extend loop iters
		bool print_zp,
		bool print_zu,
		bool print_zt) :
		truncQname_(truncQname),
		appendComment_(appendComment),
		omitsec_(omitsec),
		noUnal_(noUnal),
		pg_id_(pg_id),
		pg_pn_(pg_pn),
		pg_vn_(pg_vn),
		pg_cl_(pg_cl),
		rgs_(rgs),
		refnames_(refnames),
		reflens_(reflens),
		print_as_(print_as), // alignment score of best alignment
		print_xs_(print_xs), // alignment score of second-best alignment
		print_xss_(print_xss),
		print_yn_(print_yn), // minimum valid score and perfect score
		print_xn_(print_xn),
		print_x0_(print_x0),
		print_x1_(print_x1),
		print_xm_(print_xm),
		print_xo_(print_xo),
		print_xg_(print_xg),
		print_nm_(print_nm),
		print_md_(print_md),
		print_yf_(print_yf),
		print_yi_(print_yi),
		print_ym_(print_ym),
		print_yp_(print_yp),
		print_yt_(print_yt),
		print_ys_(print_ys),
		print_zs_(print_zs),
		print_xr_(print_xr),
		print_xt_(print_xt), // time elapsed in microseconds
		print_xd_(print_xd), // DP extend attempts
		print_xu_(print_xu), // ungapped extend attempts
		print_ye_(print_ye), // streak of failed DPs at end
		print_yl_(print_yl), // longest streak of failed DPs
		print_yu_(print_yu), // index of last succeeded DP
		print_xp_(print_xp), // print seed hit information
		print_yr_(print_yr), // index of last succeeded DP
		print_zb_(print_zb), // # Ftab lookups
		print_zr_(print_zr), // # redundant path checks
		print_zf_(print_zf), // # FM Index ops
		print_zm_(print_zm), // FM Index op string for best-first search
		print_zi_(print_zi), // # seed extend loop iters
		print_zp_(print_zp), // # seed extend loop iters
		print_zu_(print_zu), // # seed extend loop iters
		print_zt_(print_zt)  // extra features for MAPQ estimation
	{
		assert_eq(refnames_.size(), reflens_.size());
	}

	/**
	 * Print a reference name in a way that doesn't violate SAM's character
	 * constraints. \*|[!-()+-<>-~][!-~]*
	 */
	void printRefName(
		BTString& o,
		const std::string& name)
		const;

	/**
	 * Print a :Z optional field where certain characters (whitespace, colon
	 * and percent) are escaped using % escapes.
	 */
	template<typename T>
	void printOptFieldEscapedZ(BTString& o, const T& s) const {
		size_t len = s.length();
		for(size_t i = 0; i < len; i++) {
			if(s[i] < 33 || s[i] > 126 || s[i] == ':' || s[i] == '%') {
				// percent-encode it
				o.append('%');
				int ms = s[i] >> 4;
				int ls = s[i] & 15;
				assert_range(0, 15, ms);
				assert_range(0, 15, ls);
				o.append("0123456789ABCDEF"[ms]);
				o.append("0123456789ABCDEF"[ls]);
			} else {
				o.append(s[i]);
			}
		}
	}

	/**
	 * Print a :Z optional field where newline characters are escaped using %
	 * escapes.
	 */
	template<typename T>
	void printOptFieldNewlineEscapedZ(BTString& o, const T& s) const {
		size_t len = s.length();
		for(size_t i = 0; i < len; i++) {
			if(s[i] == 10 || s[i] == 13 || s[i] == '%') {
				// percent-encode it
				o.append('%');
				int ms = s[i] >> 4;
				int ls = s[i] & 15;
				assert_range(0, 15, ms);
				assert_range(0, 15, ls);
				o.append("0123456789ABCDEF"[ms]);
				o.append("0123456789ABCDEF"[ls]);
			} else {
				o.append(s[i]);
			}
		}
	}

	/**
	 * Print a read name in a way that doesn't violate SAM's character
	 * constraints. [!-?A-~]{1,255} (i.e. [33, 63], [65, 126])
	 */
	template<typename TStr>
	void printReadName(
		BTString& o,
		const TStr& name,
		bool omitSlashMate)
		const
	{
		size_t namelen = name.length();
		if(truncQname_ && namelen > 255) {
			namelen = 255;
		}
		for(size_t i = 0; i < namelen; i++) {
			if(truncQname_ && isspace(name[i])) {
				break;
			}
			o.append(name[i]);
		}
		size_t olen = o.length();
		if(omitSlashMate &&
		   olen >= 2 &&
		   o[olen-2] == '/' &&
		   (o[olen-1] == '1' || o[olen-1] == '2' || o[olen-1] == '3'))
		{
			o.resize(olen-2);
		}
	}

	/**
	 * Print a reference name given a reference index.
	 */
	void printRefNameFromIndex(
		BTString& o,
		size_t i)
		const;

	/**
	 * Print SAM header to given output buffer.
	 */
	void printHeader(
		BTString& o,
		const std::string& rgid,
		const std::string& rgs,
		bool printHd,
		bool printSq,
		bool printPg)
		const;

	/**
	 * Print the @HD header line to the given string.
	 */
	void printHdLine(BTString& o, const char *samver) const;

	/**
	 * Print the @SQ header lines to the given string.
	 */
	void printSqLines(BTString& o) const;

	/**
	 * Print the @PG header line to the given string.
	 */
	void printPgLine(BTString& o) const;

	/**
	 * Print the optional flags to the given string.
	 */
	void printAlignedOptFlags(
		BTString& o,               // output buffer
		bool first,                // first opt flag printed is first overall?
		const Read& rd,            // the read
		const Read* rdo,           // the opposite read
		AlnRes& res,               // individual alignment result
		StackedAln& staln,         // stacked alignment
		const AlnFlags& flags,     // alignment flags
		const AlnSetSumm& summ,    // summary of alignments for this read
		const SeedAlSumm& ssm,     // seed alignment summary
		const PerReadMetrics& prm, // per-read metics
		const Scoring& sc,         // scoring scheme
		const char *mapqInp)       // inputs to MAPQ calculation
		const;

	/**
	 * Print the optional flags to the given string.
	 */
	void printEmptyOptFlags(
		BTString& o,               // output buffer
		bool first,                // first opt flag printed is first overall?
		const Read& rd,            // the read
		const AlnFlags& flags,     // alignment flags
		const AlnSetSumm& summ,    // summary of alignments for this read
		const SeedAlSumm& ssm,     // seed alignment summary
		const PerReadMetrics& prm, // per-read metrics
		const Scoring& sc)         // scoring scheme
		const;

	/**
	 * Print optional flags that that have been preserve from BAM input
	 */
	void printPreservedOptFlags(BTString& o, const Read& rd) const;

	/**
	 * Append FASTA/Q comment to SAM record
	 */
        template <typename TStr>
        void printComment(BTString &o, TStr &name) const {
		if (appendComment_) {
			size_t i;

			for (i = 0; i < name.length() && !isspace(name[i]); i++) ;
			o.append('\t');
			if (i < name.length()) {
				const char *comment = name.toZBuf() + i + 1;
				if (isIllumina(comment))
					o.append("BC:Z:");
				o.append(comment);
			}
		}
	}

	bool isIllumina(const char *str) const {
		const char *start, *end;
		int field = 0;

		for (start = str, end = str; *end != '\0' && *end != ' '; end++) {
			if (*end == ':') {
				switch (field) {
                                case 0: {
					char *endptr;
					long read = strtol(start, &endptr, 10);
					if (endptr != end || (read != 1 && read != 2))
						return false;
					break;
                                }
                                case 1:
					if (*start != 'N' && *start != 'Y')
						return false;
					break;
				case 2: {
					char *endptr;
					long controlNumber = strtol(start, &endptr, 10);
					if (endptr != end || controlNumber % 2 != 0)
						return false;
					break;

				}
				default:
					return false;
				}
				start = end + 1, field++;
			}

		}
		return true;
	}

	template<typename T>
	static void readTagVal(BTString& o, const char *data, size_t &offset, size_t count) {
		std::vector<T> val(count);
		size_t i = 0;

		memcpy(val.data(), data + offset, sizeof(T) * count);
		do {
			std::string str = std::to_string(val[i]);
			o.append(str.c_str(), str.length());
			if (i < (count - 1))
				o.append(",");
		} while (++i < count);
		offset += sizeof(T) * count;
	}


        /**
	 * Return true iff we should try to obey the SAM spec's recommendations
	 * that:
	 *
	 * SEQ and QUAL of secondary alignments should be set to ‘*’ to reduce the
	 * file size.
	 */
	bool omitSecondarySeqQual() const {
		return omitsec_;
	}

	bool omitUnalignedReads() const {
		return noUnal_;
	}

	bool passthrough() const {
		return print_xr_;
	}

protected:

	bool truncQname_;   // truncate QNAME to 255 chars?
	bool appendComment_;// Append FASTA/Q comment to SAM record
	bool omitsec_;      // omit secondary
	bool noUnal_;       // omit unaligned reads

	std::string pg_id_; // @PG ID: Program record identifier
	std::string pg_pn_; // @PG PN: Program name
	std::string pg_vn_; // @PG VN: Program version
	std::string pg_cl_; // @PG CL: Program command-line
	std::string rgs_;   // Read-group string to add to all records
	const StrList& refnames_; // reference sequence names
	const LenList& reflens_;  // reference sequence lengths

	// Which alignment flags to print?

	// Following are printed by BWA-SW
	bool print_as_; // AS:i: Alignment score generated by aligner
	bool print_xs_; // XS:i: Suboptimal alignment score
	bool print_xss_;// Xs:i: Best invalid alignment score found
	bool print_yn_; // YN:i:, Yn:i: minimum valid score and perfect score
	bool print_xn_; // XN:i: Number of ambiguous bases in the referenece

	// Following are printed by BWA
	bool print_x0_; // X0:i: Number of best hits
	bool print_x1_; // X1:i: Number of sub-optimal best hits
	bool print_xm_; // XM:i: Number of mismatches in the alignment
	bool print_xo_; // XO:i: Number of gap opens
	bool print_xg_; // XG:i: Number of gap extensions (incl. opens)
	bool print_nm_; // NM:i: Edit dist. to the ref, Ns count, clipping doesn't
	bool print_md_; // MD:Z: String for mms. [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*2

	// Following are Bowtie2-specific
	bool print_yf_; // YF:i: Read was filtered out?
	bool print_yi_; // YI:Z: Summary of inputs to MAPQ calculation
	bool print_ym_; // YM:i: Read was repetitive when aligned unpaired?
	bool print_yp_; // YP:i: Read was repetitive when aligned paired?
	bool print_yt_; // YT:Z: String representing alignment type
	bool print_ys_; // YS:i: Score of other mate
	bool print_zs_; // ZS:i: Pseudo-random seed

	bool print_xr_; // XR:Z: Original read string
	bool print_xt_; // XT:i: Time taken to align
	bool print_xd_; // XD:i: DP problems
	bool print_xu_; // XU:i: ungapped alignment
	bool print_ye_; // YE:i: streak of failed DPs at end
	bool print_yl_; // YL:i: longest streak of failed DPs
	bool print_yu_; // YU:i: index of last succeeded DP
	bool print_xp_; // XP:BI: seed hit information
	bool print_yr_; // YR:i: # redundant seed hits
	bool print_zb_; // ZB:i: # Ftab lookups
	bool print_zr_; // ZR:i: # redundant path checks
	bool print_zf_; // ZF:i: # FM Index ops
	bool print_zm_; // ZM:i: FM ops string for best-first search
	bool print_zi_; // ZI:i: # extend loop iters
	bool print_zp_; // ZP:i: Score of best/second-best paired-end alignment
	bool print_zu_; // ZU:i: Score of best/second-best unpaired alignment
	bool print_zt_; // ZT:Z: Extra features for MAPQ estimation
};

#endif /* SAM_H_ */
