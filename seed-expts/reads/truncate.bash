#!/bin/bash

for fn in *_local.fastq *_semiglobal.fastq ; do
	newfn=$(echo ${fn} | sed 's/\.fastq$/_300k.fastq/')
	head -n 1200000 ${fn} > ${newfn}
done
