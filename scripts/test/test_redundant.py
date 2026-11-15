"""Tests for the check that rejects redundant alignments.

When bowtie2 looks for more than one alignment per read (-k, -a), it
rejects a new alignment that shares a dynamic-programming cell with one it
has already found for the read: the same read position aligned to the same
reference position on the same strand. The check matters when two
separate dynamic-programming windows, framed around seed hits on
different diagonals, both reach the same cells: for example across a
long gap, or when a mate found by rescue was already found as an anchor.

The first test aligns reads to a reference holding several copies of one
segment, some with a gap or mismatches relative to the read and one
reverse complemented, and checks that each copy is reported exactly once.
The second checks, on the lambda phage example, that no two alignments
(or, for pairs, no two pairs) reported for a read share a cell.

Run from the repository root after `make`:

	python3 -m pytest scripts/test/test_redundant.py
"""

import random

import pytest

import bt2_helpers as h

SEG_LEN = 300     # length of the repeated segment
GAP_AT = 150      # offset in the segment of the gap in two of the copies
READ_LEN = 150
READ_OFFSETS = (40, 75, 110)  # read start offsets in the segment
SHORT_GAP = 3     # bases missing from the 'del' copy
LONG_GAP = 20     # extra bases in the 'ins' copy


def _copies(seg, extra):
	"""The copies of seg placed in the reference, as (name, sequence, kind).

	kind is how the copy relates to seg: 'same', 'del' (SHORT_GAP bases
	missing at GAP_AT, so reads spanning it have an insertion), 'ins'
	(LONG_GAP extra bases at GAP_AT, so reads have a deletion), 'mm' (two
	mismatches) or 'rc' (reverse complement). The deletion is long enough
	that seed hits on its two sides start separate DP windows."""
	return [
		('same', seg),
		('del', seg[:GAP_AT] + seg[GAP_AT + SHORT_GAP:]),
		('ins', seg[:GAP_AT] + extra + seg[GAP_AT:]),
		('mm', h.mutate(seg, [GAP_AT - 20, GAP_AT + 20])),
		('rc', h.revcomp(seg)),
	]


@pytest.fixture(scope='module')
def repeat_ref(bt2, tmp_path_factory):
	"""Index of one chromosome holding five copies of a random segment,
	separated by random spacers, plus the reads and the start of each copy."""
	rng = random.Random(20261024)
	seg = h.random_dna(SEG_LEN, rng)
	extra = h.random_dna(LONG_GAP, rng)
	chrom, starts = h.random_dna(300, rng), {}
	for kind, copy in _copies(seg, extra):
		starts[kind] = len(chrom)
		chrom += copy + h.random_dna(300, rng)
	d = tmp_path_factory.mktemp('repeat_ref')
	fasta = h.write_fasta(str(d / 'ref.fa'), [('rep', chrom)])
	prefix = str(d / 'ref')
	bt2.build_index(fasta, prefix, extra=['--seed', '0'])
	reads = []
	for off in READ_OFFSETS:
		r = seg[off:off + READ_LEN]
		reads.append(('fw%d' % off, r, 'I' * READ_LEN))
		reads.append(('rc%d' % off, h.revcomp(r), 'I' * READ_LEN))
	fq = h.write_fastq(str(d / 'reads.fq'), reads)
	return prefix, fq, starts


def expected_hits(name, starts):
	"""The (1-based leftmost position, reverse) of the read's alignment to
	each copy. Every read spans GAP_AT, so the gaps fall inside it."""
	off = int(name[2:])
	read_rc = name.startswith('rc')
	hits = set()
	for kind in ('same', 'del', 'ins', 'mm'):
		hits.add((starts[kind] + off + 1, read_rc))
	hits.add((starts['rc'] + SEG_LEN - off - READ_LEN + 1, not read_rc))
	return hits


def cells(rec):
	"""The DP cells of an alignment, as (row, column) pairs. Rows count read
	positions from the left end of SEQ; columns are reference positions.

	This follows the cells RedundantAlns records: a read position covers one
	column, an inserted read base (I) covers the column of the read base
	after it, and deleted reference bases (D) extend the read position
	before them to the right."""
	out = set()
	row, col = 0, rec.pos - 1
	for n, op in h.cigar_ops(rec.cigar):
		if op == 'S':
			row += n
		elif op in 'M=X':
			for _ in range(n):
				out.add((row, col))
				row += 1
				col += 1
		elif op == 'I':
			for _ in range(n):
				out.add((row, col))
				row += 1
		elif op == 'D':
			for _ in range(n):
				out.add((row - 1, col))
				col += 1
		else:
			raise AssertionError('unexpected CIGAR op %s in %s' % (op, rec.cigar))
	return out


def assert_no_shared_cells(recs):
	for i in range(len(recs)):
		for j in range(i + 1, len(recs)):
			a, b = recs[i], recs[j]
			if a.rname != b.rname or a.is_reverse != b.is_reverse:
				continue
			shared = cells(a) & cells(b)
			assert not shared, 'alignments share %d cells:\n%s\n%s' % (
				len(shared), a.line, b.line)


def gaps(cigar):
	"""The gaps of a CIGAR string, e.g. '108M3I39M' -> ((3, 'I'),). Where a
	gap sits can shift with the sequence around it, so only its length and
	type are compared."""
	return tuple((n, op) for n, op in h.cigar_ops(cigar) if op in 'ID')


# The gaps of each read's five alignments, sorted
EXPECTED_GAPS = sorted([(), (), (), ((SHORT_GAP, 'I'),), ((LONG_GAP, 'D'),)])


@pytest.mark.parametrize('align_mode', ['--end-to-end', '--local'])
@pytest.mark.parametrize('k_mode', ['-k 50', '-a'])
def test_one_alignment_per_copy(bt2, repeat_ref, align_mode, k_mode):
	"""Each read is reported exactly once per copy: once with no edits, once
	with a short insertion, once with a long deletion, once with two
	mismatches and once on the other strand. Finding the same alignment again
	from another seed, or a variant of it that shares a cell, is rejected."""
	prefix, fq, starts = repeat_ref
	_, recs, _ = bt2.align(['--quiet', align_mode, '-x', prefix, '-U', fq] +
		k_mode.split())
	groups = h.by_name(recs)
	assert len(groups) == 2 * len(READ_OFFSETS)
	for name, rs in groups.items():
		assert not any(r.is_unmapped for r in rs), name
		assert sorted((r.pos, r.is_reverse) for r in rs) == \
			sorted(expected_hits(name, starts)), (name, [r.line for r in rs])
		assert sorted(gaps(r.cigar) for r in rs) == EXPECTED_GAPS, name
		assert all(h.cigar_query_length(r.cigar) == READ_LEN for r in rs), name
		assert_no_shared_cells(rs)


def assert_no_shared_pair(recs):
	"""Concordant pairs are reported as consecutive mate 1, mate 2 records.
	A mate may appear in more than one pair, but no two pairs may share a
	cell in both mates. Among the alignments of one mate outside concordant
	pairs, no two may share a cell."""
	conc = [r for r in recs if r.flag & h.PROPER_PAIR]
	assert len(conc) % 2 == 0
	for mate in (h.FIRST_IN_PAIR, h.SECOND_IN_PAIR):
		assert_no_shared_cells([r for r in recs if not r.flag & h.PROPER_PAIR
			and r.flag & mate and not r.is_unmapped])
	pairs = [(conc[i], conc[i + 1]) for i in range(0, len(conc), 2)]
	for i in range(len(pairs)):
		for j in range(i + 1, len(pairs)):
			both = True
			for a, b in zip(pairs[i], pairs[j]):
				if a.is_unmapped or b.is_unmapped or a.rname != b.rname or \
						a.is_reverse != b.is_reverse or not cells(a) & cells(b):
					both = False
					break
			assert not both, 'pairs share cells:\n%s\n%s\n%s\n%s' % (
				pairs[i][0].line, pairs[i][1].line, pairs[j][0].line,
				pairs[j][1].line)


@pytest.mark.parametrize('align_mode', ['--end-to-end', '--local'])
@pytest.mark.parametrize('paired', [False, True], ids=['unpaired', 'paired'])
def test_lambda_no_shared_cells(bt2, lambda_index, example_reads, align_mode,
		paired):
	"""With -a on the lambda example, no two alignments of a read share a
	cell. Several reads there align across deletions of 13 to 17 bases, and
	in paired mode mate rescue finds mates that were also found from their
	own seeds."""
	reads = ['-1', example_reads['reads_1'], '-2', example_reads['reads_2']] \
		if paired else ['-U', example_reads['reads_1']]
	_, recs, _ = bt2.align(['--quiet', align_mode, '-a', '-x', lambda_index] +
		reads)
	for name, rs in h.by_name(recs).items():
		if paired:
			assert_no_shared_pair(rs)
		else:
			assert_no_shared_cells(rs)
