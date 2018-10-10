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

#ifndef FORMATS_H_
#define FORMATS_H_

#include <iostream>

/**
 * File-format constants and names
 */

enum file_format {
	FASTA = 1,
	FASTA_CONT,
	FASTQ,
	BAM,
	TAB_MATE5,
	TAB_MATE6,
	RAW,
	CMDLINE,
	QSEQ,
	SRA_FASTA,
	SRA_FASTQ
};

static const std::string file_format_names[] = {
	"Invalid!",
	"FASTA",
	"FASTA sampling",
	"FASTQ",
	"BAM",
	"Tabbed 5-field",
	"Tabbed 6-field",
	"Raw",
	"Command line",
	"Qseq",
    "SRA Fasta",
    "SRA Fastq"
};

#endif /*FORMATS_H_*/
