#!/usr/bin/env python3

"""
vassess_reads.py

Given a Vargas SAM file, construct a corresponding FASTQ file that includes
all the reads with Vargas AS:i and XS:i in the read name.  Make many
replicates of each reads so we can assess how pseudo-randomness used by the
aligner affects the result.
"""

import os
import sys


def go(vargas_fn, trials):
    if not os.path.exists(vargas_fn):
        raise RuntimeException('Vargas SAM file "%s" does not exist' % vargas_fn)
    with open(vargas_fn, 'rt') as fh:
        for ln in fh:
            if ln[0] == '@':
                continue
            tokens = ln.rstrip().split('\t')
            assert len(tokens) >= 11
            name, seq, quality = tokens[0], tokens[9], tokens[10]
            asi, ssi = None, None
            for tok in tokens[11:]:
                if tok[:5] == 'AS:i:':
                    asi = int(tok[5:])
                if tok[:5] == 'ss:i:':
                    ssi = int(tok[5:])
                if asi is not None and ssi is not None:
                    break
            assert asi is not None and ssi is not None
            read_body = '\n'.join([seq, '+', quality])
            for i in range(trials):
                print('@%s:%d:%d:%d' % (name, asi, ssi, i))
                print(read_body)


if __name__ == '__main__':
    # TODO: make this a parameter
    trials = 20
    if len(sys.argv) < 2:
        raise ValueError('Expected argument: (1) Vargas SAM')
    go(sys.argv[1], trials)
