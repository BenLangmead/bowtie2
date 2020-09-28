#!/usr/bin/env python3

"""
TODO: specify a particular setup as being the "reference" and compare against it
Then we can automate the process of checking if any of our experiments exceed old S/VS modes on 2 out of the 3 major criteria:

- lower n_unal
- lower as_i diff
- lower time
"""

import os
import sys
from collections import defaultdict

if len(sys.argv) < 6:
    raise ValueError('Requires arguments (1) aligner, (2) dataset, '
                     '(3) target series, (4) target command, (5) query series')

aligner = sys.argv[1]
dataset = sys.argv[2]
target_series = sys.argv[3]
target_cmd = sys.argv[4]
query_series = sys.argv[5]

target_fn = os.path.join(target_series, '%s_%s.cmds.csv' % (aligner, dataset))
if not os.path.exists(target_fn):
    raise RuntimeError('No such target file "%s"' % target_fn)
target = {}
with open(target_fn, 'rt') as fh:
    for ln in fh:
        tokens = ln.split(',')
        if tokens[0] == target_cmd:
            target['time'] = float(tokens[1])
            target[tokens[2]] = float(tokens[3])

assert 'time' in target
assert 'as_diff_mean' in target
assert 'n_unal' in target
assert 'xs_diff_mean' in target

query_fn = os.path.join(query_series, '%s_%s.cmds.csv' % (aligner, dataset))
if not os.path.exists(query_fn):
    raise RuntimeError('No such query file "%s"' % query_fn)
query = defaultdict(dict)
with open(query_fn, 'rt') as fh:
    # command,elapsed,statistic,mean,med,max,min
    # def0_vf,5.9200220108,as_diff_mean,18.2343866598,16.05,52.0,-1
    # def0_vf,5.9200220108,xs_diff_mean,17.1496019176,18.65,47.0,-12.0
    for ln in fh:
        tokens = ln.split(',')
        if tokens[0] == 'command':
            continue
        cmd, time, stat, mean, med, mx, mn = tokens
        query[cmd]['time'] = float(time)
        query[cmd][stat] = float(mean)

time_as = []
time_nunal = []
as_nunal = []
all3 = []

for cmd in query.keys():
    time_diff = query[cmd]['time'] - target['time']
    as_diff_mean_diff = query[cmd]['as_diff_mean'] - target['as_diff_mean']
    n_unal_diff = query[cmd]['n_unal'] - target['n_unal']
    if time_diff <= 0 and as_diff_mean_diff <= 0:
        time_as.append((cmd, time_diff, as_diff_mean_diff, n_unal_diff))
    if time_diff <= 0 and n_unal_diff <= 0:
        time_nunal.append((cmd, time_diff, as_diff_mean_diff, n_unal_diff))
    if as_diff_mean_diff <= 0 and n_unal_diff <= 0:
        as_nunal.append((cmd, time_diff, as_diff_mean_diff, n_unal_diff))
    if time_diff <= 0 and as_diff_mean_diff <= 0 and n_unal_diff <= 0:
        all3.append((cmd, time_diff, as_diff_mean_diff, n_unal_diff))


print('Time and AS:i better:')
for cmd, time_diff, as_diff_mean_diff, n_unal_diff in time_as:
    print(','.join(map(str, [cmd, time_diff, as_diff_mean_diff, n_unal_diff])))

print('Time and n_unal better:')
for cmd, time_diff, as_diff_mean_diff, n_unal_diff in time_nunal:
    print(','.join(map(str, [cmd, time_diff, as_diff_mean_diff, n_unal_diff])))

print('AS:i and n_unal better:')
for cmd, time_diff, as_diff_mean_diff, n_unal_diff in as_nunal:
    print(','.join(map(str, [cmd, time_diff, as_diff_mean_diff, n_unal_diff])))

print('All 3 better:')
for cmd, time_diff, as_diff_mean_diff, n_unal_diff in all3:
    print(','.join(map(str, [cmd, time_diff, as_diff_mean_diff, n_unal_diff])))
