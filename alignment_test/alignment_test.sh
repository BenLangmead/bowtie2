#!/bin/bash

#PATH=$PATH:/data1/igm3/sw/bin; export PATH

ARGS=( "$@" )


#needs to be much more comprehensie
if [ $# -lt 1 ]
then
    echo 
    echo alignment_test.sh -h for help
    echo
    exit
elif [ $1 == '-h' ] || [ $1 == '-H' ] || [ $1 == '--help'  ]
then
    echo 
    echo To run this test:
    echo alignment_test.sh --runMason \<fasta_filename\> \<Index Base Name\> --runA \<path to run A\> \[options\] --runB \<path to run B\> \[options\] -x \<\# alignments to check\>
    echo 
    echo Note: If you already have a Mason-generated \(verbose\) fastq file, leave out --runMason and replace \<fasta_filename\> with \<fastq_filename\>
    echo
    exit
fi

count=1
if [ $1 == '--runMason' ]
then
    dotFASTA=$2
    dotFastaLength=${#dotFASTA}
    FASTA=${dotFASTA:0:dotFastaLength - 3}
    mason illumina -N 10000 -sq -i -o ./${FASTA}.fq $dotFASTA
    (( count++ ))
else
    dotFASTA=$1
    dotFastaLength=${#dotFASTA}
    FASTA=${dotFASTA:0:dotFastaLength - 3}
fi

idxName=${ARGS[$count]}


(( count += 2 ))
runApath=${ARGS[$count]}
runACommands=""
(( count++ ))
while [ $count -lt $# ]
do
    arg=${ARGS[$count]}
    if [ "$arg" == "--runB" ]
    then
	(( count++ ))
	break
    else
	runACommands+=" ${arg}"
	(( count++ ))
    fi
done

runBpath=${ARGS[$count]}
(( count++ ))
runBCommands=""
while [ $count -lt $# ]
do
    arg=${ARGS[$count]}
    if [ "$arg" == "-x" ]
    then
	(( count++ ))
	break
    else
	runBCommands+=" ${arg}"
	(( count++ ))
    fi
done

if [ ${ARGS[$count - 1]} == "-x" ]
then
    num_alignments_to_check=${ARGS[$count]}
else
    num_alignments_to_check=10000
fi


/usr/bin/time -v  -o runA.time $runApath $runACommands -x ./$idxName -U $FASTA.fq -S runA_alignments.sam

/usr/bin/time -v -o runB.time $runBpath $runBCommands -x ./$idxName -U $FASTA.fq -S runB_alignments.sam

awk 'NR==5' runA.time > runA_totals.time
awk 'NR==10' runA.time > runA_totals.mem
awk 'NR==5' runB.time > runB_totals.time
awk 'NR==10' runB.time > runB_totals.mem

/usr/bin/time -v  -o runA.time $runApath $runACommands -x ./$idxName -U $FASTA.fq -S runA_alignments.sam

/usr/bin/time -v -o runB.time $runBpath $runBCommands -x ./$idxName -U $FASTA.fq -S runB_alignments.sam

awk 'NR==5' runA.time >> runA_totals.time
awk 'NR==10' runA.time >> runA_totals.mem
awk 'NR==5' runB.time >> runB_totals.time
awk 'NR==10' runB.time >> runB_totals.mem

/usr/bin/time -v  -o runA.time $runApath $runACommands -x ./$idxName -U $FASTA.fq -S runA_alignments.sam

/usr/bin/time -v -o runB.time $runBpath $runBCommands -x ./$idxName -U $FASTA.fq -S runB_alignments.sam

awk 'NR==5' runA.time >> runA_totals.time
awk 'NR==10' runA.time >> runA_totals.mem
awk 'NR==5' runB.time >> runB_totals.time
awk 'NR==10' runB.time >> runB_totals.mem

python average_timeANDmem.py runA_totals.time runA_totals.mem runB_totals.time runB_totals.mem

rm runA.time
rm runB.time
rm runA_totals.time
rm runB_totals.time
rm runA_totals.mem
rm runB_totals.mem

python check_alignment_accuracy.py $FASTA.fq runA_alignments.sam runB_alignments.sam $num_alignments_to_check
