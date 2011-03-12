#!/bin/sh

SEED="$1"
if [ -z "$SEED" ] ; then
	SEED=77
else
	shift
fi

if [ "$1" == "-c" ] ; then
	rm -f bowtie bowtie-debug bowtie-build-debug bowtie-inspect-debug
	shift
fi

make $* bowtie bowtie-debug bowtie-build-debug bowtie-inspect-debug

# Args: seed, outer, inner, tbase, trand, pbase, prand

echo "Short test emphasizing searching..."
perl scripts/random_bowtie_tests.pl $SEED 1000 200 300 200 8 30
if [ "$?" != "0" ] ; then echo "Error!" ; exit 1 ; fi
echo "Short test emphasizing searching with short patterns..."
perl scripts/random_bowtie_tests.pl $SEED 1000 200 300 200 6 6
if [ "$?" != "0" ] ; then echo "Error!" ; exit 1 ; fi
echo "Short test emphasizing searching with long patterns..."
perl scripts/random_bowtie_tests.pl $SEED 1000 200 300 200 30 20
if [ "$?" != "0" ] ; then echo "Error!" ; exit 1 ; fi
echo "Short test emphasizing building..."
perl scripts/random_bowtie_tests.pl $SEED 1200 10 300 200 8 30
if [ "$?" != "0" ] ; then echo "Error!" ; exit 1 ; fi

echo "Long test emphasizing searching..."
perl scripts/random_bowtie_tests.pl $SEED 5000 200 300 200 8 30
if [ "$?" != "0" ] ; then echo "Error!" ; exit 1 ; fi
echo "Long test emphasizing building..."
perl scripts/random_bowtie_tests.pl $SEED 6000 10 300 200 8 30
if [ "$?" != "0" ] ; then echo "Error!" ; exit 1 ; fi
