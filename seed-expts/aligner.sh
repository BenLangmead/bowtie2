#!/usr/bin/env bash

set -ex

for ver in versions/* ; do
  ver="${ver:9}"
  echo "=== Version ${ver} ==="
  for mode in local semiglobal ; do
    echo "  --- Mode ${mode} ---"
    while IFS=" " read -r genome index ; do
      READS_SUMM="${genome}_${mode}_${ver}.reads.csv"
      CMDS_SUMM="${genome}_${mode}_${ver}.cmds.csv"
      echo "    ooo Genome ${genome} ooo"
      if [[ ! -f ${READS_SUMM} || ! -f ${CMDS_SUMM} ]] ; then
        python vassess_aligner.py \
            cmds_${mode}.txt \
            "${ver}" \
            "reads/${genome}_${mode}.fastq" \
            "${index}" \
            "${READS_SUMM}" \
            "${CMDS_SUMM}"
      else
        echo "*** Found ${READS_SUMM} ${CMDS_SUMM} so skipping... ***"
      fi
    done < genome_index.txt
  done
done
