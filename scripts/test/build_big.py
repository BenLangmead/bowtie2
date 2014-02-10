#!/usr/bin/env python

import sys
import os
import random
import numpy

def gen_genome(fh, ref_id, tot_len, per_line=1000):
    fh.write('>%s\n' % ref_id)
    while tot_len > 0:
        ln = min(tot_len, per_line)
        assert ln > 0
        fh.write(''.join(numpy.random.choice(list('ACGT'), ln)) + '\n')
        tot_len -= ln

if __name__ == '__main__':
    tot_len = int(sys.argv[1])
    tmp_fasta = sys.argv[2]
    tmp_idx = sys.argv[3]
    bt2_build_exe = sys.argv[4]
    num_refs = 1
    if len(sys.argv) > 5:
        num_refs = int(sys.argv[5])
    with open(tmp_fasta, 'w') as fh:
        for i in xrange(num_refs):
            gen_genome(fh, 'ref%d' % i, tot_len // num_refs)
    ret = os.system('%s %s %s' % (bt2_build_exe, tmp_fasta, tmp_idx))
    print >> sys.stderr, "Exitlevel: %d" % ret
    exit(ret)
