#!/bin/bash

# A prototype for scalability test.
# Run it from bowtie2 directory. 


SAMPLE_POINTS=2
MAX_THREADS=5
#CMD_SKEL="/usr/bin/time -f \"%U %S %E %PCPU \" ./bowtie2 -x ../data/genomes/hg19 -U ../data/reads/ERR000589.fastq -p "
CMD_SKEL="/usr/bin/time -f %U,%S,%E ./bowtie2 -x ecoli_s -U example/reads/longreads.fq -p "
DATA_DIR=$(pwd)
DATA_DIR+="/scripts/test/scalability"
DATA_FILE_PREFIX=$1

get_sample_points() {
    local cmd="$@"
    for ((i=0; i<SAMPLE_POINTS; i++))
    do
        echo $cmd >&2
        $($cmd)
    done
}

collect_all_data() {
    local no_threads=$1
    local cmd=$CMD_SKEL
    local data_file=$DATA_FILE_PREFIX
    local full_cmd
    data_file+="sample.dat"
    
    for ((t=1; t<MAX_THREADS; t++))
    do
        full_cmd="$cmd $t >/dev/null" 
        get_sample_points $full_cmd 2>>"${DATA_DIR}/data/${t}_${data_file}"
    done
}


collect_all_data

# format data definitely maybe
