"""Alignment against a prebuilt large (.bt2l) index of more than 2^32 bases.

The rest of the suite checks the large-index code path on small genomes
forced into the .bt2l format with --large-index. Only an index whose joined
reference really exceeds 2^32 bases exercises offsets that do not fit in 32
bits, and building one takes hours and tens of GB of memory, so these tests
are skipped unless BOWTIE2_LARGE_INDEX names the prefix of such an index
(for example a human plus mouse index on a cluster).

Reads are cut from the index's own reference, obtained by streaming
`bowtie2-inspect` one sequence at a time, at positions past the 2^32 offset
in the joined reference (which leaves out runs of N) and straddling it. Each
read must align perfectly, and every read whose alignment is unique must
land where it was cut.

	BOWTIE2_LARGE_INDEX=/path/to/hs_mm python3 -m pytest scripts/test/test_large_genome.py

Optional settings:

	BOWTIE2_LARGE_INDEX_THREADS   -p for the alignment run (default 4)
	BOWTIE2_LARGE_INDEX_OFFSET    joined offset the reads must reach (default
	                              2^32); lowering it lets the test logic be
	                              checked against a small .bt2l index
"""

import os
import subprocess

import pytest

import bt2_helpers as h


INDEX = os.environ.get('BOWTIE2_LARGE_INDEX')
THRESHOLD = int(os.environ.get('BOWTIE2_LARGE_INDEX_OFFSET', 2 ** 32))
THREADS = os.environ.get('BOWTIE2_LARGE_INDEX_THREADS', '4')

READ_LEN = 100
FRAG_LEN = 300
FAR_READS = 150
CONTROL_READS = 5

pytestmark = pytest.mark.skipif(not INDEX,
	reason='set BOWTIE2_LARGE_INDEX to a .bt2l index of more than 2^32 bases')


def inspect_sequences(bt2, prefix):
	"""Yield (name, sequence as bytes) for each reference sequence, streaming
	bowtie2-inspect so that only one sequence is held in memory."""
	proc = subprocess.Popen([bt2.inspect, prefix], stdout=subprocess.PIPE)
	name, seq = None, bytearray()
	try:
		for line in proc.stdout:
			line = line.rstrip(b'\r\n')
			if line.startswith(b'>'):
				if name is not None:
					yield name, bytes(seq)
				name, seq = line[1:].decode().split()[0], bytearray()
			else:
				seq.extend(line)
		if name is not None:
			yield name, bytes(seq)
	finally:
		proc.stdout.close()
		assert proc.wait() == 0, 'bowtie2-inspect failed on %s' % prefix


def joined_offset(seq, start, p):
	"""Offset in the joined reference of position p of seq, where start is the
	joined offset of seq's first non-N base."""
	return start + p - seq.count(b'N', 0, p)


def position_at(seq, start, target):
	"""Smallest position p of seq whose joined offset is at least target."""
	lo, hi = 0, len(seq)
	while lo < hi:
		mid = (lo + hi) // 2
		if joined_offset(seq, start, mid) < target:
			lo = mid + 1
		else:
			hi = mid
	return lo


def is_clean(seq, q, length):
	return 0 <= q <= len(seq) - length and seq[q:q + length].translate(None, b'ACGT') == b''


def clean_window(seq, p, length, limit=10000):
	"""First position at or after p where seq[p:p+length] is all A/C/G/T, or None."""
	for q in range(p, min(p + limit, len(seq) - length + 1)):
		if is_clean(seq, q, length):
			return q
	return None


def straddling_window(seq, p, length):
	"""A position q with q < p < q+length and seq[q:q+length] all A/C/G/T, the
	one with p nearest the middle, or None."""
	for d in range(length // 2):
		for q in (p - length // 2 - d, p - length // 2 + d):
			if q < p < q + length and is_clean(seq, q, length):
				return q
	return None


@pytest.fixture(scope='module')
def large_index(bt2):
	assert os.path.exists(INDEX + '.1.bt2l'), \
		'BOWTIE2_LARGE_INDEX=%s: no %s.1.bt2l' % (INDEX, INDEX)
	return INDEX


@pytest.fixture(scope='module')
def cut_reads(bt2, large_index):
	"""Reads and pairs cut from the reference, with where they came from.

	Returns (unpaired, pairs, total joined length). unpaired holds
	(name, seq, rname, pos, joined offset) and pairs holds
	(name, mate1, mate2, rname, pos1, pos2, joined offset); positions are
	1-based as in SAM."""
	# Reads beyond the threshold are cut at evenly spaced joined offsets, so
	# most come from long sequences rather than from the many short, often
	# repetitive unplaced contigs an assembly may end with. The summary's
	# lengths include runs of N, so the spacing is an overestimate at worst.
	summary = bt2.run('bowtie2-inspect', ['-s', large_index]).stdout
	total_len = sum(int(l.split('\t')[-1]) for l in summary.splitlines()
		if l.startswith('Sequence-'))
	spacing = max(FRAG_LEN, (total_len - THRESHOLD) // FAR_READS)
	unpaired, pairs = [], []
	start, next_pick = 0, THRESHOLD
	for rname, seq in inspect_sequences(bt2, large_index):
		end = joined_offset(seq, start, len(seq))
		picks = []
		if not unpaired:
			picks += [len(seq) * (i + 1) // (CONTROL_READS + 1) for i in range(CONTROL_READS)]
		if start <= THRESHOLD < end:
			# A read straddling the threshold.
			q = straddling_window(seq, position_at(seq, start, THRESHOLD), READ_LEN)
			if q is not None:
				picks.append(q)
		while next_pick < end:
			picks.append(position_at(seq, start, next_pick))
			next_pick += spacing
		for p in picks:
			q = clean_window(seq, p, READ_LEN)
			if q is None:
				continue
			off = joined_offset(seq, start, q)
			name = 'r%d_%s_%d' % (len(unpaired), rname, q + 1)
			unpaired.append((name, seq[q:q + READ_LEN].decode(), rname, q + 1, off))
			f = clean_window(seq, p, FRAG_LEN)
			if f is not None and off >= THRESHOLD:
				frag = seq[f:f + FRAG_LEN].decode()
				pairs.append(('p%d_%s_%d' % (len(pairs), rname, f + 1), frag[:READ_LEN],
					h.revcomp(frag[-READ_LEN:]), rname, f + 1, f + FRAG_LEN - READ_LEN + 1,
					joined_offset(seq, start, f)))
		start = end
	return unpaired, pairs, start


def test_index_is_larger_than_threshold(cut_reads):
	"""The joined reference is longer than the offset the tests target, so
	reads beyond it exist."""
	_, _, total = cut_reads
	assert total > THRESHOLD, ('the index holds %d bases of non-N reference, not more '
		'than %d' % (total, THRESHOLD))


def test_large_offsets_align_where_cut(bt2, large_index, cut_reads, tmp_path):
	"""Unpaired reads cut beyond the threshold, and across it, align perfectly,
	and pairs cut beyond it align concordantly. Uniquely aligned reads and
	pairs land where they were cut, pairs with the right TLEN."""
	unpaired, pairs, _ = cut_reads
	far = [r for r in unpaired if r[4] >= THRESHOLD]
	assert far, 'no read could be cut beyond offset %d' % THRESHOLD
	assert any(r[4] < THRESHOLD < r[4] + READ_LEN for r in unpaired), \
		'no read straddles offset %d' % THRESHOLD
	assert pairs, 'no pair could be cut beyond offset %d' % THRESHOLD

	u = h.write_fastq(str(tmp_path / 'u.fq'), [(r[0], r[1]) for r in unpaired])
	m1 = h.write_fastq(str(tmp_path / 'm1.fq'), [(p[0], p[1]) for p in pairs])
	m2 = h.write_fastq(str(tmp_path / 'm2.fq'), [(p[0], p[2]) for p in pairs])
	# Unpaired and paired reads go in separate runs: mixing them in one run
	# with -p > 1 loses the pairs or hangs.
	recs = {}
	for name, reads in (('u', ['-U', u]), ('p', ['-1', m1, '-2', m2])):
		sam = str(tmp_path / ('%s.sam' % name))
		bt2.run('bowtie2', ['-x', large_index, '-p', THREADS, '-S', sam] + reads,
			timeout=7200)
		with open(sam) as f:
			recs.update(h.by_name(h.primary(h.parse_sam(f.read())[1])))

	def unique(rec):
		return 'XS' not in rec.tags or rec.tags['XS'] < rec.tags['AS']

	# Repeats are common in large genomes, so only a quarter of the reads are
	# required to align uniquely; every one that does must be in place.
	located = 0
	for name, _, rname, pos, off in unpaired:
		(rec,) = recs[name]
		assert not rec.is_unmapped and rec.tags['AS'] == 0, rec
		if unique(rec):
			assert (rec.rname, rec.pos) == (rname, pos), rec
			located += off >= THRESHOLD
	assert located * 4 >= len(far), \
		'only %d of %d reads beyond the threshold aligned uniquely' % (located, len(far))

	# For a pair in a repeat, bowtie2 may report a concordant placement that is
	# not the best one (one mate with mismatches and XS:i above AS:i), so a
	# perfect score is required only of pairs whose mates both align uniquely.
	located = 0
	for name, _, _, rname, pos1, pos2, _ in pairs:
		r1, r2 = sorted(recs[name], key=lambda r: not r.flag & h.FIRST_IN_PAIR)
		for r in (r1, r2):
			assert not r.is_unmapped and r.tags['YT'] == 'CP', r
		if unique(r1) and unique(r2):
			assert r1.tags['AS'] == 0 and r2.tags['AS'] == 0, (r1, r2)
			assert (r1.rname, r1.pos, r1.tlen) == (rname, pos1, FRAG_LEN), r1
			assert (r2.rname, r2.pos, r2.tlen) == (rname, pos2, -FRAG_LEN), r2
			located += 1
	assert located * 4 >= len(pairs), \
		'only %d of %d pairs aligned uniquely' % (located, len(pairs))
