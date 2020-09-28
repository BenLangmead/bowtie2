#!/usr/bin/env bash

set -ex

for mode in local semiglobal ; do
  while IFS=" " read -r genome index ; do
    READS_SUMM="${genome}_${mode}.reads.csv"
    CMDS_SUMM="${genome}_${mode}.cmds.csv"
    if [[ ! -f ${READS_SUMM} || ! -f ${CMDS_SUMM} ]] ; then
      python vassess_aligner.py \
          cmds_${mode}.txt \
          "reads/${genome}_${mode}.fastq" \
          "${index}" \
          "${READS_SUMM}" \
          "${CMDS_SUMM}"
    else
      echo "*** Found ${READS_SUMM} ${CMDS_SUMM} so skipping... ***"
    fi
  done < genome_index.txt
done
