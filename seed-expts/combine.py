#!/usr/bin/env python3

"""
Takes a Vargas SAM file and an aligner-generated SAM file containing the same
reads, and adds all the Vargas annotations to the aligner-generated SAM.
"""

import sys
from collections import defaultdict

fn1 = sys.argv[1]
fn2 = sys.argv[2]

scrape = ['mc:i:', 'mp:i:', 'AS:i:', 'sc:i:', 'su:Z:', 'ss:i:', 'sp:i:']

# ERR239486.11882    0    chr7    58784192    255    *    *    0    0    CTCTGTTTGTAAAGTCTGCAAGTGGATATATGGAACTCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATTTCATGCTAGACAGAAGAATTCTCAGTA    DCEEFHEGGFFGGIHFGHGGHGFHHHDGGHGHCGHHGHGHGHFGHHHGGGFGGHFDEGHGGHFFHGDGGHDIHGGGHFIGCIGJGHHIHGC@CGFGHGBH    st:Z:fwd    sc:i:1327    su:Z:chr7    ss:i:-15    sp:i:58785313    RG:Z:VAUGRP    AS:i:-10    mp:i:58784291    mc:i:4

varg_results = defaultdict(dict)

with open(fn1, 'rt') as fh:
    for ln in fh:
        if ln[0] == '@':
            continue
        toks = ln.split('\t')
        assert len(toks) >= 11
        for tok in toks[12:]:
            for scr in scrape:
                if tok.startswith(scr):
                    varg_results[toks[0]][scr] = tok[len(scr):].strip()

with open(fn2, 'rt') as fh:
    for ln in fh:
        if ln[0] == '@':
            print(ln, end='')
            continue
        toks = ln.split('\t')
        assert toks[0] in varg_results
        print(ln.rstrip(), end='')
        for k, v in varg_results[toks[0]].items():
            print('\t' + k + v, end='')
        print()
