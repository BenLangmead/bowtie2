#!/bin/sh

#
# Downloads sequence for A. thaliana from TAIR.  This script was used
# to build the Bowtie index for A. thaliana.  The downloaded version
# was TAIR9.
#

GENOMES_MIRROR=ftp://ftp.arabidopsis.org/home/tair

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

for c in 1 2 3 4 5 C M ; do
	if [ ! -f chr$c.fas ] ; then
		F=${GENOMES_MIRROR}/Sequences/whole_chromosomes/chr$c.fas
		get $F || (echo "Error getting $F" && exit 1)
	fi
	
	if [ ! -f chr$c.fas ] ; then
		echo "Could not find chr$c.fas file!"
		exit 2
	fi
done

CMD="${BOWTIE_BUILD_EXE} $* chr1.fas,chr2.fas,chr3.fas,chr4.fas,chr5.fas,chrC.fas,chrM.fas  a_thaliana"
echo $CMD
if $CMD ; then
	echo "a_thaliana index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi
