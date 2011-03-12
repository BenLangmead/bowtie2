#ifndef FORMATS_H_
#define FORMATS_H_

#include <iostream>
#include <seqan/sequence.h>

/**
 * File-format constants and names
 */

enum file_format {
	FASTA = 1,
	FASTA_CONT,
	FASTQ,
	TAB_MATE,
	RAW,
	CMDLINE,
	INPUT_CHAIN,
	RANDOM,
	QSEQ
};

static const std::string file_format_names[] = {
	"Invalid!",
	"FASTA",
	"FASTA sampling",
	"FASTQ",
	"Tabbed mated",
	"Raw",
	"Command line",
	"Chain file",
	"Random",
	"Qseq"
};

/**
 * Print the given read information as a FASTA record.
 */
static inline void printFastaRecord(
		std::ostream& o,
		const seqan::String<char>& name,
		const seqan::String<seqan::Dna5>& seq)
{
	o << ">" << name << endl << seq << endl;
}

/**
 * Print the given read information as a FASTQ record.
 */
static inline void printFastqRecord(
		std::ostream& o,
		const seqan::String<char>& name,
		const seqan::String<seqan::Dna5>& seq,
		const seqan::String<char>& qual)
{
	o << "@" << name << endl << seq << endl << "+" << endl << qual << endl;
}

#endif /*FORMATS_H_*/
