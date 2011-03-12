#!/bin/sh

#
# Downloads sequence for H. sapiens (human) from NCBI.  This script was
# used to build the Bowtie index for H. sapiens 37.
#
# From README_CURRENT_BUILD:
# Organism: Homo sapiens (human)
# NCBI Build Number: 37    
# Version: 1
# Release date: 04 August 2009
#

BASE_CHRS="\
chr1 \
chr2 \
chr3 \
chr4 \
chr5 \
chr6 \
chr7 \
chr8 \
chr9 \
chr10 \
chr11 \
chr12 \
chr13 \
chr14 \
chr15 \
chr16 \
chr17 \
chr18 \
chr19 \
chr20 \
chr21 \
chr22 \
chrX \
chrY"

CHRS_TO_INDEX=$BASE_CHRS

FTP_BASE=ftp://ftp.ncbi.nih.gov/genomes/H_sapiens
FTP_ASM_BASE=$FTP_BASE/Assembled_chromosomes
FTP_MT_BASE=$FTP_BASE/CHR_MT

get() {
	file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o `basename $1` $1
		return $?
	else
		wget $1
		return $?
	fi
}

BOWTIE_BUILD_EXE=./bowtie-build
if [ ! -x "$BOWTIE_BUILD_EXE" ] ; then
	if ! which bowtie-build ; then
		echo "Could not find bowtie-build in current directory or in PATH"
		exit 1
	else
		BOWTIE_BUILD_EXE=`which bowtie-build`
	fi
fi

INPUTS=
BASE_NAME=hs_ref_GRCh37_
for c in $CHRS_TO_INDEX ; do
	c=${BASE_NAME}${c}
	if [ ! -f $c.fa ] ; then
		F=$c.fa.gz
		get ${FTP_ASM_BASE}/$F || (echo "Error getting $F" && exit 1)
		gunzip $F || (echo "Error unzipping $F" && exit 1)
	fi
	[ -n "$INPUTS" ] && INPUTS=$INPUTS,$c.fa
	[ -z "$INPUTS" ] && INPUTS=$c.fa
done

# Special case: get mitochondrial DNA from its home
if [ ! -f hs_ref_chrMT.fa ] ; then
	F=hs_ref_chrMT.fa.gz
	get ${FTP_MT_BASE}/$F || (echo "Error getting $F" && exit 1)
	gunzip $F || (echo "Error unzipping $F" && exit 1)
fi

INPUTS=$INPUTS,hs_ref_chrMT.fa

echo Running ${BOWTIE_BUILD_EXE} $* ${INPUTS} h_sapiens_37_asm
${BOWTIE_BUILD_EXE} $* ${INPUTS} h_sapiens_37_asm

if [ "$?" = "0" ] ; then
	echo "h_sapiens_37_asm index built:"
	echo "   h_sapiens_37_asm.1.ebwt h_sapiens_37_asm.2.ebwt"
	echo "   h_sapiens_37_asm.3.ebwt h_sapiens_37_asm.4.ebwt"
	echo "   h_sapiens_37_asm.rev.1.ebwt h_sapiens_37_asm.rev.2.ebwt"
	echo "You may remove hs_ref_chr*.fa"
else
	echo "Index building failed; see error message"
fi
