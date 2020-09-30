#!/bin/bash

for fn in *_local.fastq *_semiglobal.fastq ; do
	newfn=$(echo ${fn} | sed 's/\.fastq$/_100k.fastq/')
	head -n 400000 ${fn} > ${newfn}
done
