"""Paired-end alignment tests for bowtie2.

Every pair is cut from the fixed random `genome` (conftest.py) at a known
position, orientation and fragment length, so the expected FLAG, POS, RNEXT,
PNEXT, TLEN and YT:Z of both mates are derived here before bowtie2 runs. The
same applies to the classification of overlapping, containing and dovetailing
mates (the MANUAL's own examples), the fragment-length bounds -I/-X, mixed
and discordant reporting, the reads written by --un-conc/--al-conc, and the
counts in the paired alignment summary.

Run from the repository root after `make`:

	python3 -m pytest scripts/test/test_paired.py
"""

import bz2
import gzip
import os
import random
import shutil
import struct
import subprocess
import tempfile
import zlib

import pytest

import bt2_helpers as h
from bt2_helpers import (PAIRED, PROPER_PAIR, UNMAPPED, MATE_UNMAPPED, REVERSE,
	MATE_REVERSE, FIRST_IN_PAIR, SECOND_IN_PAIR)


# ---------------------------------------------------------------------------
# Helpers (candidates for bt2_helpers.py)
# ---------------------------------------------------------------------------

class Fragment(object):
	"""Mates cut from both ends of ref[start:start + length].

	`orient` is the library orientation ('fr', 'rf' or 'ff', as in the
	MANUAL's --fr/--rf/--ff) and `strand` the strand the fragment came from
	('+' Watson, '-' Crick). For a Watson fragment mate 1 is upstream: with
	'fr' it reads forward and mate 2 is the reverse complement of the
	downstream end; 'rf' swaps the strands of the two mates; 'ff' reads both
	forward. A Crick fragment is the reverse complement of the Watson one, cut
	the same way. `pos`, `rev` and `tlen` hold the expected SAM POS (1-based),
	reverse-strand bit and TLEN of mate 1 and mate 2 when the pair aligns.
	"""

	def __init__(self, ref, start, length, len1, len2=None, orient='fr',
			strand='+', chrom='chr1'):
		len2 = len1 if len2 is None else len2
		frag = ref[start:start + length]
		flip = strand == '-'
		if flip:
			frag = h.revcomp(frag)
		left, right = frag[:len1], frag[length - len2:]
		if orient == 'fr':
			self.mate1, self.mate2 = left, h.revcomp(right)
		elif orient == 'rf':
			self.mate1, self.mate2 = h.revcomp(left), right
		elif orient == 'ff':
			self.mate1, self.mate2 = left, right
		else:
			raise ValueError(orient)
		self.chrom = chrom
		self.length = length
		self.rev = ((orient == 'rf') != flip, (orient == 'fr') != flip)
		if flip:
			self.pos = (start + length - len1 + 1, start + 1)
			self.tlen = (-length, length)
		else:
			self.pos = (start + 1, start + length - len2 + 1)
			self.tlen = (length, -length)

	def pair(self, name):
		return (name, self.mate1, self.mate2)


def pair_flags(rev1, rev2, proper):
	"""Expected FLAGs of mate 1 and mate 2 when both mates align."""
	base = PAIRED | (PROPER_PAIR if proper else 0)
	f1 = base | FIRST_IN_PAIR | (REVERSE if rev1 else 0) | (MATE_REVERSE if rev2 else 0)
	f2 = base | SECOND_IN_PAIR | (REVERSE if rev2 else 0) | (MATE_REVERSE if rev1 else 0)
	return f1, f2


def junk(length, seed):
	"""Random sequence that aligns nowhere in the 10 kbp test genome."""
	return h.random_dna(length, random.Random(seed))


def mismatching(seq):
	"""A sequence differing from `seq` at every position (soft-clip filler)."""
	return h.mutate(seq, range(len(seq)))


def write_pair_files(directory, pairs):
	"""Write (name, mate1, mate2[, qual1, qual2]) pairs as two FASTQ files."""
	r1, r2 = [], []
	for p in pairs:
		q1 = p[3] if len(p) > 3 else None
		q2 = p[4] if len(p) > 4 else None
		r1.append((p[0], p[1], q1))
		r2.append((p[0], p[2], q2))
	m1 = h.write_fastq(os.path.join(directory, 'reads_1.fq'), r1)
	m2 = h.write_fastq(os.path.join(directory, 'reads_2.fq'), r2)
	return m1, m2


def new_dir(tmp_path):
	return tempfile.mkdtemp(dir=str(tmp_path))


def align_pairs(bt2, index, tmp_path, pairs, args=(), check=True):
	"""Align pairs given with -1/-2 at -p 1; returns (records, proc)."""
	m1, m2 = write_pair_files(new_dir(tmp_path), pairs)
	_, records, proc = bt2.align(['-x', index, '-p', '1', '-1', m1, '-2', m2]
		+ list(args), check=check)
	return records, proc


def mates_by_name(records):
	"""qname -> [mate 1 record, mate 2 record], primary records only."""
	out = {}
	for r in h.primary(records):
		assert r.flag & PAIRED, r
		slot = 0 if r.flag & FIRST_IN_PAIR else 1
		assert out.setdefault(r.qname, [None, None])[slot] is None, r
		out[r.qname][slot] = r
	for name, (r1, r2) in out.items():
		assert r1 is not None and r2 is not None, name
	return out


def check_aligned_pair(r1, r2, frag, proper, yt):
	"""Both mates aligned where `frag` says, with mate fields filled in."""
	f1, f2 = pair_flags(frag.rev[0], frag.rev[1], proper)
	assert (r1.flag, r2.flag) == (f1, f2)
	assert (r1.rname, r2.rname) == (frag.chrom, frag.chrom)
	assert (r1.pos, r2.pos) == frag.pos
	assert (r1.rnext, r2.rnext) == ('=', '=')
	assert (r1.pnext, r2.pnext) == (frag.pos[1], frag.pos[0])
	assert (r1.tags['YT'], r2.tags['YT']) == (yt, yt)


def summary_block(stderr):
	"""The alignment summary lines of stderr, from 'N reads' to the rate."""
	lines = stderr.splitlines()
	start = [i for i, l in enumerate(lines) if l.endswith(' reads; of these:')]
	end = [i for i, l in enumerate(lines) if l.endswith(' overall alignment rate')]
	assert start and end, stderr
	return lines[start[-1]:end[-1] + 1]


def paired_summary(npairs, conc, disc, mates, discordant=True, mixed=True):
	"""Expected paired alignment summary, in the MANUAL's format.

	conc = pairs aligning concordantly (0 times, exactly 1 time, >1 times);
	disc = pairs aligning discordantly; mates = the mates of the remaining
	pairs aligning (0 times, exactly 1 time, >1 times).
	"""
	def pct(n, d):
		return '%.2f%%' % (100.0 * n / d if d else 0.0)
	lines = [
		'%d reads; of these:' % npairs,
		'  %d (%s) were paired; of these:' % (npairs, pct(npairs, npairs)),
		'    %d (%s) aligned concordantly 0 times' % (conc[0], pct(conc[0], npairs)),
		'    %d (%s) aligned concordantly exactly 1 time' % (conc[1], pct(conc[1], npairs)),
		'    %d (%s) aligned concordantly >1 times' % (conc[2], pct(conc[2], npairs))]
	if discordant:
		lines += [
			'    ----',
			'    %d pairs aligned concordantly 0 times; of these:' % conc[0],
			'      %d (%s) aligned discordantly 1 time' % (disc, pct(disc, conc[0]))]
	if mixed:
		rest = conc[0] - disc
		lines += [
			'    ----',
			'    %d pairs aligned 0 times concordantly or discordantly; of these:' % rest,
			'      %d mates make up the pairs; of these:' % (2 * rest),
			'        %d (%s) aligned 0 times' % (mates[0], pct(mates[0], 2 * rest)),
			'        %d (%s) aligned exactly 1 time' % (mates[1], pct(mates[1], 2 * rest)),
			'        %d (%s) aligned >1 times' % (mates[2], pct(mates[2], 2 * rest))]
	unpaired = mates[1] + mates[2] if mixed else 0
	aligned = 2 * (conc[1] + conc[2]) + 2 * disc + unpaired
	lines.append('%s overall alignment rate' % pct(aligned, 2 * npairs))
	return lines


def _bgzf_block(data):
	comp = zlib.compressobj(6, zlib.DEFLATED, -15)
	cdata = comp.compress(data) + comp.flush()
	header = struct.pack('<4BI2BH', 31, 139, 8, 4, 0, 0, 255, 6) + b'BC' + \
		struct.pack('<HH', 2, len(cdata) + 25)
	return header + cdata + struct.pack('<II', zlib.crc32(data) & 0xffffffff, len(data))


BGZF_EOF = bytes([0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42,
	0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0, 0, 0, 0, 0, 0, 0, 0])


def write_unaligned_bam(path, records):
	"""Write an unaligned BAM (no @SQ) from (name, flag, seq, qual) records,
	qual in Phred+33 text. Pure Python: BGZF blocks plus the EOF marker."""
	codes = dict((b, i) for i, b in enumerate('=ACMGRSVTWYHKDBN'))
	out = [b'BAM\x01', struct.pack('<i', 0), struct.pack('<i', 0)]
	for name, flag, seq, qual in records:
		qname = name.encode() + b'\x00'
		packed = bytearray((len(seq) + 1) // 2)
		for i, base in enumerate(seq):
			packed[i // 2] |= codes[base] << (4 * (1 - i % 2))
		body = struct.pack('<iiBBHHHiiii', -1, -1, len(qname), 0, 4680, 0, flag,
			len(seq), -1, -1, 0) + qname + bytes(packed) + \
			bytes(ord(c) - 33 for c in qual)
		out.append(struct.pack('<i', len(body)) + body)
	data = b''.join(out)
	with open(path, 'wb') as f:
		for i in range(0, len(data), 60000):
			f.write(_bgzf_block(data[i:i + 60000]))
		f.write(BGZF_EOF)
	return path


def fastq_records(text):
	"""Split FASTQ text into 4-line record strings."""
	lines = text.splitlines(True)
	assert len(lines) % 4 == 0, text
	return [''.join(lines[i:i + 4]) for i in range(0, len(lines), 4)]


def read_compressed(path, kind):
	"""Text of a file written by the wrapper with --*-gz/-bz2/-zst/-lz4."""
	if kind == 'gz':
		with gzip.open(path, 'rt') as f:
			return f.read()
	if kind == 'bz2':
		with bz2.open(path, 'rt') as f:
			return f.read()
	tool = {'zst': 'zstd', 'lz4': 'lz4'}[kind]
	proc = subprocess.run([tool, '-dc', path], capture_output=True, check=True)
	return proc.stdout.decode()


MAGIC = {'gz': b'\x1f\x8b', 'bz2': b'BZh', 'zst': b'\x28\xb5\x2f\xfd',
	'lz4': b'\x04\x22\x4d\x18'}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

REPEAT = ('chr2', 1000, 1400)   # copied into chr4 of `repeat_index`
REPEAT_COPY_START = 300         # offset of the copy in chr4


@pytest.fixture(scope='module')
def repeat_index(bt2, genome, tmp_path_factory):
	"""`genome` plus chr4, which holds a copy of chr2[1000:1400] between
	random flanks, so pairs from that region align concordantly twice."""
	g = dict(genome)
	rng = random.Random(4242)
	copy = g[REPEAT[0]][REPEAT[1]:REPEAT[2]]
	chr4 = h.random_dna(REPEAT_COPY_START, rng) + copy + h.random_dna(300, rng)
	d = tmp_path_factory.mktemp('repeat_index')
	fasta = h.write_fasta(str(d / 'repeat.fa'), list(genome) + [('chr4', chr4)])
	prefix = str(d / 'repeat')
	bt2.build_index(fasta, prefix, extra=['--seed', '0'])
	return prefix


# ---------------------------------------------------------------------------
# -1/-2, orientation (--fr/--rf/--ff), FLAG, RNEXT/PNEXT/TLEN, YT, YS
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('orient', ['fr', 'rf', 'ff'])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_concordant_pair_fields(bt2, genome_dict, genome_index, tmp_path, orient, strand):
	"""A pair built in --fr/--rf/--ff orientation aligns concordantly with that option: FLAG, POS, RNEXT/PNEXT, TLEN, YT:Z:CP, YS:i."""
	frag = Fragment(genome_dict['chr1'], 1000, 300, 50, 40, orient=orient, strand=strand)
	records, _ = align_pairs(bt2, genome_index, tmp_path, [frag.pair('p')],
		['--' + orient])
	r1, r2 = mates_by_name(records)['p']
	check_aligned_pair(r1, r2, frag, proper=True, yt='CP')
	assert (r1.tlen, r2.tlen) == frag.tlen
	assert (r1.cigar, r2.cigar) == ('50M', '40M')
	assert (r1.tags['YS'], r2.tags['YS']) == (0, 0)


def test_fr_is_default(bt2, genome_dict, genome_index, tmp_path):
	"""--fr is the default orientation: identical SAM records with and without it."""
	c = genome_dict['chr2']
	pairs = [Fragment(c, 100, 250, 50, orient=o, strand=s).pair('%s%s' % (o, s))
		for o in ('fr', 'rf', 'ff') for s in '+-']
	default, _ = align_pairs(bt2, genome_index, tmp_path, pairs)
	fr, _ = align_pairs(bt2, genome_index, tmp_path, pairs, ['--fr'])
	assert [r.line for r in default] == [r.line for r in fr]
	assert [r.tags['YT'] for r in default] == ['CP'] * 4 + ['DP'] * 8


@pytest.mark.parametrize('built,option', [('fr', 'rf'), ('fr', 'ff'), ('rf', 'fr'),
	('rf', 'ff'), ('ff', 'fr'), ('ff', 'rf')])
def test_wrong_orientation_is_discordant(bt2, genome_dict, genome_index, tmp_path, built, option):
	"""Mates in an orientation other than the one --fr/--rf/--ff asks for align discordantly: no 0x2, YT:Z:DP, TLEN still set on the same chromosome."""
	frag = Fragment(genome_dict['chr1'], 2000, 300, 50, orient=built)
	records, _ = align_pairs(bt2, genome_index, tmp_path, [frag.pair('p')],
		['--' + option])
	r1, r2 = mates_by_name(records)['p']
	check_aligned_pair(r1, r2, frag, proper=False, yt='DP')
	assert (r1.tlen, r2.tlen) == frag.tlen
	assert (r1.tags['YS'], r2.tags['YS']) == (0, 0)


def test_mate_names_lose_slash_suffix(bt2, genome_dict, genome_index, tmp_path):
	"""Mate names ending /1 and /2 are reported under one QNAME, with 0x40 on mate 1 and 0x80 on mate 2."""
	frag = Fragment(genome_dict['chr3'], 200, 220, 45, chrom='chr3')
	d = new_dir(tmp_path)
	m1 = h.write_fastq(os.path.join(d, 'r1.fq'), [('frag7/1 comment', frag.mate1)])
	m2 = h.write_fastq(os.path.join(d, 'r2.fq'), [('frag7/2 comment', frag.mate2)])
	_, records, _ = bt2.align(['-x', genome_index, '-p', '1', '-1', m1, '-2', m2])
	assert [r.qname for r in records] == ['frag7', 'frag7']
	assert records[0].flag & FIRST_IN_PAIR and not records[0].flag & SECOND_IN_PAIR
	assert records[1].flag & SECOND_IN_PAIR and not records[1].flag & FIRST_IN_PAIR
	check_aligned_pair(records[0], records[1], frag, proper=True, yt='CP')


def test_ys_is_opposite_mate_score(bt2, genome_dict, genome_index, tmp_path):
	"""YS:i of each mate equals AS:i of the other (1 and 3 mismatches at Phred 40: AS -6 and -18)."""
	frag = Fragment(genome_dict['chr1'], 1500, 300, 50)
	pair = ('p', h.mutate(frag.mate1, [45]), h.mutate(frag.mate2, [2, 25, 47]))
	records, _ = align_pairs(bt2, genome_index, tmp_path, [pair])
	r1, r2 = mates_by_name(records)['p']
	check_aligned_pair(r1, r2, frag, proper=True, yt='CP')
	assert (r1.tags['AS'], r2.tags['AS']) == (-6, -18)
	assert (r1.tags['YS'], r2.tags['YS']) == (-18, -6)


@pytest.mark.parametrize('orient', ['fr', 'rf', 'ff'])
@pytest.mark.parametrize('option', ['--nofw', '--norc'])
def test_nofw_norc_pertain_to_fragments(bt2, genome_dict, genome_index, tmp_path, orient, option):
	"""In paired mode --nofw keeps only fragments from the Crick strand and --norc only those from the Watson strand."""
	c = genome_dict['chr1']
	watson = Fragment(c, 1000, 300, 50, orient=orient, strand='+')
	crick = Fragment(c, 2000, 300, 50, orient=orient, strand='-')
	records, _ = align_pairs(bt2, genome_index, tmp_path,
		[watson.pair('watson'), crick.pair('crick')], ['--' + orient, option])
	mates = mates_by_name(records)
	kept, dropped = ('crick', 'watson') if option == '--nofw' else ('watson', 'crick')
	frag = crick if kept == 'crick' else watson
	check_aligned_pair(mates[kept][0], mates[kept][1], frag, proper=True, yt='CP')
	assert [r.flag for r in mates[dropped]] == [
		PAIRED | UNMAPPED | MATE_UNMAPPED | FIRST_IN_PAIR,
		PAIRED | UNMAPPED | MATE_UNMAPPED | SECOND_IN_PAIR]


# ---------------------------------------------------------------------------
# -I/--minins and -X/--maxins
# ---------------------------------------------------------------------------

def gap_pair(ref, start, gap, mate_len=20):
	"""Two mate_len-bp FR mates separated by `gap` bp (the MANUAL's -I/-X
	examples use 20 bp mates)."""
	return Fragment(ref, start, 2 * mate_len + gap, mate_len, chrom='chr2')


@pytest.mark.parametrize('option', ['-I', '--minins'])
def test_minins_boundary(bt2, genome_dict, genome_index, tmp_path, option):
	"""MANUAL example: with -I 60 two 20 bp mates 20 bp apart are concordant, 19 bp apart are not."""
	c = genome_dict['chr2']
	ok, short = gap_pair(c, 100, 20), gap_pair(c, 300, 19)
	records, _ = align_pairs(bt2, genome_index, tmp_path,
		[ok.pair('ok'), short.pair('short')], [option, '60'])
	mates = mates_by_name(records)
	check_aligned_pair(mates['ok'][0], mates['ok'][1], ok, proper=True, yt='CP')
	check_aligned_pair(mates['short'][0], mates['short'][1], short, proper=False, yt='DP')
	assert [r.tlen for r in mates['ok']] == [60, -60]
	assert [r.tlen for r in mates['short']] == [59, -59]


@pytest.mark.parametrize('option', ['-X', '--maxins'])
def test_maxins_boundary(bt2, genome_dict, genome_index, tmp_path, option):
	"""MANUAL example: with -X 100 two 20 bp mates 60 bp apart are concordant, 61 bp apart are not."""
	c = genome_dict['chr2']
	ok, wide = gap_pair(c, 500, 60), gap_pair(c, 700, 61)
	records, _ = align_pairs(bt2, genome_index, tmp_path,
		[ok.pair('ok'), wide.pair('wide')], [option, '100'])
	mates = mates_by_name(records)
	check_aligned_pair(mates['ok'][0], mates['ok'][1], ok, proper=True, yt='CP')
	check_aligned_pair(mates['wide'][0], mates['wide'][1], wide, proper=False, yt='DP')
	assert [r.tlen for r in mates['ok']] == [100, -100]
	assert [r.tlen for r in mates['wide']] == [101, -101]


def test_minins_maxins_together(bt2, genome_dict, genome_index, tmp_path):
	"""-I 120 -X 120 accepts only a 120 bp fragment; 119 and 121 bp are discordant."""
	c = genome_dict['chr1']
	frags = dict((n, Fragment(c, s, n, 40)) for n, s in ((119, 300), (120, 600), (121, 900)))
	records, _ = align_pairs(bt2, genome_index, tmp_path,
		[f.pair(str(n)) for n, f in frags.items()], ['-I', '120', '-X', '120'])
	mates = mates_by_name(records)
	for n, frag in frags.items():
		r1, r2 = mates[str(n)]
		check_aligned_pair(r1, r2, frag, proper=(n == 120), yt='CP' if n == 120 else 'DP')
		assert (r1.tlen, r2.tlen) == frag.tlen


def test_maxins_default_is_500(bt2, genome_dict, genome_index, tmp_path):
	"""Without -X a 500 bp fragment is concordant and a 501 bp fragment is not (default -X 500)."""
	c = genome_dict['chr1']
	ok, wide = Fragment(c, 1000, 500, 50), Fragment(c, 2500, 501, 50)
	records, _ = align_pairs(bt2, genome_index, tmp_path, [ok.pair('ok'), wide.pair('wide')])
	mates = mates_by_name(records)
	check_aligned_pair(mates['ok'][0], mates['ok'][1], ok, proper=True, yt='CP')
	check_aligned_pair(mates['wide'][0], mates['wide'][1], wide, proper=False, yt='DP')


@pytest.mark.parametrize('args,message', [
	(['-X', '0'], '-X arg must be at least 1'),
	(['-I', '-1'], '-I arg must be positive')])
def test_insert_size_bad_values(bt2, genome_dict, genome_index, tmp_path, args, message):
	"""-X below 1 and negative -I are rejected with bowtie2's own message."""
	frag = Fragment(genome_dict['chr1'], 1000, 300, 50)
	_, proc = align_pairs(bt2, genome_index, tmp_path, [frag.pair('p')], args, check=False)
	assert proc.returncode != 0
	assert message in proc.stderr


# ---------------------------------------------------------------------------
# --no-overlap, --no-contain, --dovetail (the MANUAL's example geometries)
# ---------------------------------------------------------------------------

# Offsets of mate 1 and mate 2 inside the MANUAL's 58 bp example reference,
# mate 1 expected on the left (FR): overlap, the first containment example
# (mates share their right end), and the dovetail example. 'abut' places the
# mates end to end without overlap.
GEOMETRY = {
	'overlap': ((0, 31), (27, 58)),
	'contain': ((0, 49), (27, 49)),
	'dovetail': ((13, 49), (8, 44)),
	'abut': ((0, 30), (30, 60)),
}

GEOMETRY_CONCORDANT = {
	None: {'overlap': True, 'contain': True, 'dovetail': False, 'abut': True},
	'--no-overlap': {'overlap': False, 'contain': False, 'dovetail': False, 'abut': True},
	'--no-contain': {'overlap': True, 'contain': False, 'dovetail': False, 'abut': True},
	'--dovetail': {'overlap': True, 'contain': True, 'dovetail': True, 'abut': True},
}


@pytest.mark.parametrize('option', [None, '--no-overlap', '--no-contain', '--dovetail'])
def test_overlap_contain_dovetail(bt2, genome_dict, genome_index, tmp_path, option):
	"""Overlapping and containing mates are concordant by default, dovetailing ones are not; --no-overlap, --no-contain and --dovetail change exactly that."""
	c = genome_dict['chr1']
	base = dict(zip(GEOMETRY, (3000, 3200, 3400, 3600)))
	pairs, expect = [], {}
	for name, ((s1, e1), (s2, e2)) in GEOMETRY.items():
		b = base[name]
		pairs.append((name, c[b + s1:b + e1], h.revcomp(c[b + s2:b + e2])))
		up, down = min(s1, s2), max(e1, e2)
		sign = 1 if s1 <= s2 else -1
		expect[name] = (b + s1 + 1, b + s2 + 1, sign * (down - up))
	records, _ = align_pairs(bt2, genome_index, tmp_path, pairs,
		[option] if option else [])
	mates = mates_by_name(records)
	for name in GEOMETRY:
		r1, r2 = mates[name]
		concordant = GEOMETRY_CONCORDANT[option][name]
		pos1, pos2, tlen = expect[name]
		f1, f2 = pair_flags(False, True, concordant)
		assert (r1.flag, r2.flag) == (f1, f2), name
		assert (r1.pos, r2.pos) == (pos1, pos2), name
		assert (r1.tlen, r2.tlen) == (tlen, -tlen), name
		yt = 'CP' if concordant else 'DP'
		assert (r1.tags['YT'], r2.tags['YT']) == (yt, yt), name


def test_contained_mate_sharing_left_end(bt2, genome_dict, genome_index, tmp_path):
	"""Mate 1 contained in a longer mate 2 starting at the same base is concordant by default and discordant with --no-contain."""
	c = genome_dict['chr2']
	pair = ('p', c[2000:2030], h.revcomp(c[2000:2060]))
	for args, proper in (([], True), (['--no-contain'], False)):
		records, _ = align_pairs(bt2, genome_index, tmp_path, [pair], args)
		r1, r2 = mates_by_name(records)['p']
		assert (r1.flag, r2.flag) == pair_flags(False, True, proper)
		assert (r1.pos, r2.pos) == (2001, 2001)
		assert (r1.cigar, r2.cigar) == ('30M', '60M')
		# Both mates start at 2001, so which one gets the positive TLEN is
		# arbitrary in the SAM specification; mate 1 pins the current output.
		assert (r1.tlen, r2.tlen) == (60, -60)


# ---------------------------------------------------------------------------
# Mixed mode, --no-mixed, discordant pairs, --no-discordant
# ---------------------------------------------------------------------------

def one_mate_pairs(genome_dict):
	"""A concordant pair, a pair whose mate 2 aligns nowhere, and a pair whose
	mate 1 aligns nowhere and whose mate 2 aligns to the reverse strand."""
	c = genome_dict['chr1']
	conc = Fragment(c, 1000, 300, 50)
	return conc, [conc.pair('conc'),
		('m1only', c[1500:1550], junk(50, 1)),
		('m2only', junk(50, 2), h.revcomp(c[2250:2300]))]


def test_mixed_mode_reports_aligned_mate(bt2, genome_dict, genome_index, tmp_path):
	"""When one mate fails, the other is reported alone (YT:Z:UP, 0x8) and the failed mate takes its RNAME/POS; TLEN is 0."""
	conc, pairs = one_mate_pairs(genome_dict)
	records, _ = align_pairs(bt2, genome_index, tmp_path, pairs)
	mates = mates_by_name(records)
	check_aligned_pair(mates['conc'][0], mates['conc'][1], conc, proper=True, yt='CP')

	al, un = mates['m1only']
	assert al.flag == PAIRED | MATE_UNMAPPED | FIRST_IN_PAIR
	assert un.flag == PAIRED | UNMAPPED | SECOND_IN_PAIR
	assert (al.rname, al.pos, al.cigar) == ('chr1', 1501, '50M')
	assert (al.rnext, al.pnext, al.tlen) == ('=', 1501, 0)
	assert (un.rname, un.pos, un.cigar, un.mapq) == ('chr1', 1501, '*', 0)
	assert (un.rnext, un.pnext, un.tlen) == ('=', 1501, 0)

	un, al = mates['m2only']
	assert al.flag == PAIRED | MATE_UNMAPPED | SECOND_IN_PAIR | REVERSE
	assert un.flag == PAIRED | UNMAPPED | FIRST_IN_PAIR | MATE_REVERSE
	assert (al.rname, al.pos, al.cigar) == ('chr1', 2251, '50M')
	assert (al.rnext, al.pnext, al.tlen) == ('=', 2251, 0)
	assert (un.rname, un.pos, un.rnext, un.pnext, un.tlen) == ('chr1', 2251, '=', 2251, 0)

	for name in ('m1only', 'm2only'):
		for r in mates[name]:
			assert r.tags['YT'] == 'UP'
			assert 'YS' not in r.tags


def test_no_mixed_suppresses_lone_mates(bt2, genome_dict, genome_index, tmp_path):
	"""--no-mixed reports a pair with one aligning mate as wholly unaligned; concordant pairs are unchanged."""
	conc, pairs = one_mate_pairs(genome_dict)
	default, _ = align_pairs(bt2, genome_index, tmp_path, pairs)
	records, _ = align_pairs(bt2, genome_index, tmp_path, pairs, ['--no-mixed'])
	mates = mates_by_name(records)
	assert [r.line for r in mates['conc']] == [r.line for r in mates_by_name(default)['conc']]
	for name in ('m1only', 'm2only'):
		r1, r2 = mates[name]
		assert r1.flag == PAIRED | UNMAPPED | MATE_UNMAPPED | FIRST_IN_PAIR
		assert r2.flag == PAIRED | UNMAPPED | MATE_UNMAPPED | SECOND_IN_PAIR
		for r in (r1, r2):
			assert (r.rname, r.pos, r.cigar, r.rnext, r.pnext, r.tlen) == ('*', 0, '*', '*', 0, 0)
			assert r.tags['YT'] == 'UP'


def discordant_pairs(genome_dict):
	"""Unique mates on two chromosomes, and unique FF mates on one chromosome
	(discordant under the default --fr)."""
	c1, c3 = genome_dict['chr1'], genome_dict['chr3']
	same = Fragment(c1, 3000, 250, 50, orient='ff')
	return same, [('diffchr', c1[2000:2050], h.revcomp(c3[1500:1540])), same.pair('samechr')]


def test_discordant_pair_fields(bt2, genome_dict, genome_index, tmp_path):
	"""Uniquely aligned mates that break the pair constraints are reported discordant: YT:Z:DP, no 0x2, mate fields set, TLEN only on one chromosome."""
	same, pairs = discordant_pairs(genome_dict)
	records, _ = align_pairs(bt2, genome_index, tmp_path, pairs)
	mates = mates_by_name(records)
	r1, r2 = mates['diffchr']
	assert (r1.flag, r2.flag) == pair_flags(False, True, False)
	assert (r1.rname, r1.pos, r1.rnext, r1.pnext) == ('chr1', 2001, 'chr3', 1501)
	assert (r2.rname, r2.pos, r2.rnext, r2.pnext) == ('chr3', 1501, 'chr1', 2001)
	assert (r1.tlen, r2.tlen) == (0, 0)
	assert (r1.tags['YT'], r2.tags['YT']) == ('DP', 'DP')
	assert (r1.tags['YS'], r2.tags['YS']) == (0, 0)
	r1, r2 = mates['samechr']
	check_aligned_pair(r1, r2, same, proper=False, yt='DP')
	assert (r1.tlen, r2.tlen) == same.tlen != (0, 0)


def test_no_discordant(bt2, genome_dict, genome_index, tmp_path):
	"""--no-discordant reports would-be discordant mates as unpaired alignments: YT:Z:UP, TLEN 0, no YS:i, mate fields kept."""
	same, pairs = discordant_pairs(genome_dict)
	records, _ = align_pairs(bt2, genome_index, tmp_path, pairs, ['--no-discordant'])
	mates = mates_by_name(records)
	r1, r2 = mates['diffchr']
	assert (r1.flag, r2.flag) == pair_flags(False, True, False)
	assert (r1.rnext, r1.pnext, r2.rnext, r2.pnext) == ('chr3', 1501, 'chr1', 2001)
	r1, r2 = mates['samechr']
	check_aligned_pair(r1, r2, same, proper=False, yt='UP')
	for r in mates['diffchr'] + mates['samechr']:
		assert r.tags['YT'] == 'UP'
		assert r.tlen == 0
		assert 'YS' not in r.tags


def test_no_discordant_no_mixed(bt2, genome_dict, genome_index, tmp_path):
	"""--no-discordant with --no-mixed leaves pairs that are not concordant wholly unaligned."""
	same, pairs = discordant_pairs(genome_dict)
	records, _ = align_pairs(bt2, genome_index, tmp_path, pairs,
		['--no-discordant', '--no-mixed'])
	mates = mates_by_name(records)
	assert sorted(mates) == ['diffchr', 'samechr']
	for r1, r2 in mates.values():
		assert r1.flag == PAIRED | UNMAPPED | MATE_UNMAPPED | FIRST_IN_PAIR
		assert r2.flag == PAIRED | UNMAPPED | MATE_UNMAPPED | SECOND_IN_PAIR


def test_discordant_needs_unique_mates(bt2, genome_dict, repeat_index, tmp_path):
	"""A pair is not called discordant when one mate aligns to two places; both mates are reported unpaired (YT:Z:UP)."""
	c2, c3 = genome_dict['chr2'], genome_dict['chr3']
	pair = ('p', c2[1100:1150], h.revcomp(c3[100:150]))
	records, _ = align_pairs(bt2, repeat_index, tmp_path, [pair])
	r1, r2 = mates_by_name(records)['p']
	assert (r1.flag, r2.flag) == pair_flags(False, True, False)
	assert (r1.tags['YT'], r2.tags['YT']) == ('UP', 'UP')
	assert (r1.cigar, r1.tags['AS']) == ('50M', 0)
	assert (r1.rname, r1.pos) in (('chr2', 1101), ('chr4', REPEAT_COPY_START + 101))
	assert (r2.rname, r2.pos) == ('chr3', 101)
	assert (r1.tlen, r2.tlen) == (0, 0)


def test_concordant_repeat_reports_one_pair(bt2, genome_dict, repeat_index, tmp_path):
	"""A pair inside a duplicated region is concordant at either copy, with XS:i equal to AS:i for both mates."""
	c2 = genome_dict['chr2']
	frag = Fragment(c2, 1020, 330, 50, chrom='chr2')
	records, _ = align_pairs(bt2, repeat_index, tmp_path, [frag.pair('p')])
	r1, r2 = mates_by_name(records)['p']
	assert (r1.flag, r2.flag) == pair_flags(False, True, True)
	assert r1.rname == r2.rname
	offset = {'chr2': 0, 'chr4': REPEAT_COPY_START - REPEAT[1]}[r1.rname]
	assert (r1.pos, r2.pos) == (frag.pos[0] + offset, frag.pos[1] + offset)
	assert (r1.tlen, r2.tlen) == frag.tlen
	for r in (r1, r2):
		assert r.tags['AS'] == r.tags['XS'] == 0


# ---------------------------------------------------------------------------
# MC:Z and MQ:i (off by default, enabled with --sam-opt-config)
# ---------------------------------------------------------------------------

def mate_tag_pairs(genome_dict):
	"""Concordant and discordant pairs whose mates differ in length, edits
	and alignment score, so their CIGARs and MAPQs differ. The ten mates
	have ten different CIGARs (with or without --xeq or --local), so an MC:Z
	equals the mate's CIGAR only if it was taken from that mate; one left
	over from an earlier record cannot match by chance."""
	c1, c2, c3 = genome_dict['chr1'], genome_dict['chr2'], genome_dict['chr3']
	return [
		('lengths', c1[1000:1052], h.revcomp(c1[1258:1300])),
		('mismatch', c2[1200:1246], h.revcomp(h.mutate(c2[1460:1500], [30]))),
		('deletion', c2[200:248], h.revcomp(c2[400:425] + c2[427:452])),
		('mapq', h.mutate(c3[100:150], [1, 3, 45, 47]), h.revcomp(c3[300:330])),
		('discordant', h.mutate(c1[2600:2644], [38, 41]), h.revcomp(c2[500:535])),
	]


def test_mc_mq_off_by_default(bt2, genome_dict, genome_index, tmp_path):
	"""MC:Z and MQ:i are not printed unless enabled with --sam-opt-config."""
	records, _ = align_pairs(bt2, genome_index, tmp_path, mate_tag_pairs(genome_dict))
	assert len(records) == 10
	for r in records:
		assert 'MC' not in r.tags and 'MQ' not in r.tags


@pytest.mark.parametrize('extra', [[], ['--xeq'], ['--local']])
def test_mc_mq_describe_the_mate(bt2, genome_dict, genome_index, tmp_path, extra):
	"""With --sam-opt-config MC,MQ each mate of a concordant or discordant pair carries the other mate's CIGAR (MC:Z) and MAPQ (MQ:i).

	In the current version the 'mapq' pair's two mates get different MAPQs,
	so there MQ:i is told apart from the record's own MAPQ. That difference
	is not required, since it depends on how MAPQ is computed.
	"""
	records, _ = align_pairs(bt2, genome_index, tmp_path, mate_tag_pairs(genome_dict),
		['--sam-opt-config', 'MC,MQ'] + extra)
	mates = mates_by_name(records)
	assert set(mates) == set(['lengths', 'mismatch', 'deletion', 'mapq', 'discordant'])
	for name, (r1, r2) in mates.items():
		expect_yt = 'DP' if name == 'discordant' else 'CP'
		assert (r1.tags['YT'], r2.tags['YT']) == (expect_yt, expect_yt), name
		assert r1.tags['MC'] == r2.cigar and r2.tags['MC'] == r1.cigar, name
		assert r1.tags['MQ'] == r2.mapq and r2.tags['MQ'] == r1.mapq, name
	assert 'D' in mates['deletion'][1].cigar
	cigars = [r.cigar for pair in mates.values() for r in pair]
	assert len(set(cigars)) == len(cigars), cigars
	if not extra:
		assert mates['lengths'][0].tags['MC'] == '42M'
		assert mates['lengths'][1].tags['MC'] == '52M'
	if extra == ['--xeq']:
		assert mates['mismatch'][0].tags['MC'] == '30=1X9='


def test_mc_mq_toggle_off(bt2, genome_dict, genome_index, tmp_path):
	"""--sam-opt-config MC,MQ,-MC,-MQ leaves the SAM output as without the option."""
	pairs = mate_tag_pairs(genome_dict)
	default, _ = align_pairs(bt2, genome_index, tmp_path, pairs)
	toggled, _ = align_pairs(bt2, genome_index, tmp_path, pairs,
		['--sam-opt-config', 'MC,MQ,-MC,-MQ'])
	assert len(default) == 2 * len(pairs)
	assert [r.line for r in toggled] == [r.line for r in default]


def test_mq_absent_when_mate_unaligned(bt2, genome_dict, genome_index, tmp_path):
	"""MQ:i is not printed on a mate whose partner did not align, nor on the unaligned mate.

	The SAM specification would also allow MQ:i on the unaligned record (its
	mate has a MAPQ), so that half pins the current output. MC:Z is checked
	only on the first pair of the file: a lone mate after an aligned pair
	carries a stale MC:Z in the current version, which is not pinned here.
	"""
	_, pairs = one_mate_pairs(genome_dict)
	records, _ = align_pairs(bt2, genome_index, tmp_path, pairs[1:],
		['--sam-opt-config', 'MC,MQ'])
	mates = mates_by_name(records)
	for name in ('m1only', 'm2only'):
		for r in mates[name]:
			assert 'MQ' not in r.tags, r
	# The first pair in the file has no earlier pair whose CIGAR could leak.
	assert 'MC' not in mates['m1only'][0].tags


# ---------------------------------------------------------------------------
# --soft-clipped-unmapped-tlen
# ---------------------------------------------------------------------------

def clipped_pair(genome_dict):
	"""FR pair on chr1: mate 1 = 6 non-matching bases + chr1[1000:1040];
	mate 2 = reverse complement of chr1[1200:1240] + 8 non-matching bases. In
	--local mode mate 1 is 6S40M at 1001 and mate 2 40M8S at 1201."""
	c = genome_dict['chr1']
	m1 = mismatching(c[994:1000]) + c[1000:1040]
	m2 = h.revcomp(c[1200:1240] + mismatching(c[1240:1248]))
	return ('clip', m1, m2)


@pytest.mark.parametrize('flag,tlen', [(False, 254), (True, 240)])
def test_soft_clipped_unmapped_tlen(bt2, genome_dict, genome_index, tmp_path, flag, tlen):
	"""In --local mode TLEN spans soft-clipped bases (994..1248 = 254) unless --soft-clipped-unmapped-tlen (1000..1240 = 240)."""
	args = ['--local'] + (['--soft-clipped-unmapped-tlen'] if flag else [])
	records, _ = align_pairs(bt2, genome_index, tmp_path, [clipped_pair(genome_dict)], args)
	r1, r2 = mates_by_name(records)['clip']
	assert (r1.flag, r2.flag) == pair_flags(False, True, True)
	assert (r1.cigar, r2.cigar) == ('6S40M', '40M8S')
	assert (r1.pos, r2.pos) == (1001, 1201)
	assert (r1.pnext, r2.pnext) == (1201, 1001)
	assert (r1.tlen, r2.tlen) == (tlen, -tlen)


def test_soft_clipped_unmapped_tlen_needs_local(bt2, genome_dict, genome_index, tmp_path):
	"""--soft-clipped-unmapped-tlen without --local is an error."""
	_, proc = align_pairs(bt2, genome_index, tmp_path, [clipped_pair(genome_dict)],
		['--soft-clipped-unmapped-tlen'], check=False)
	assert proc.returncode != 0
	assert '--soft-clipped-unmapped-tlen can only be set for local alignments' in proc.stderr


# ---------------------------------------------------------------------------
# --un-conc / --al-conc (+ -gz/-bz2/-zst/-lz4, % pattern), --un-mates
# ---------------------------------------------------------------------------

def conc_split_input(genome_dict, directory):
	"""Write four pairs as FASTQ with comments, '+name' lines and distinct
	qualities; returns (m1 path, m2 path, records by mate, names by class)."""
	c1, c2 = genome_dict['chr1'], genome_dict['chr2']
	conc = Fragment(c1, 1000, 300, 50)
	seqs = [
		('conc', conc.mate1, conc.mate2),
		('m1only', c1[1500:1550], junk(50, 3)),
		('none', junk(50, 4), junk(50, 5)),
		('disc', c1[2600:2650], h.revcomp(c2[500:550])),
	]
	text = {1: [], 2: []}
	for i, (name, s1, s2) in enumerate(seqs):
		for mate, seq in ((1, s1), (2, s2)):
			qual = ''.join(chr(43 + (7 * j + 3 * i + mate) % 31) for j in range(len(seq)))
			plus = '+%s/%d' % (name, mate) if i % 2 else '+'
			text[mate].append('@%s/%d extra words %d\n%s\n%s\n%s\n' % (
				name, mate, i, seq, plus, qual))
	paths = []
	for mate in (1, 2):
		p = os.path.join(directory, 'in_%d.fq' % mate)
		with open(p, 'w') as f:
			f.write(''.join(text[mate]))
		paths.append(p)
	return paths[0], paths[1], text, {'al': [0], 'un': [1, 2, 3]}


def run_with_read_files(bt2, index, directory, m1, m2, args):
	sam = os.path.join(directory, 'out.sam')
	proc = bt2.run('bowtie2', ['-x', index, '-p', '1', '-1', m1, '-2', m2, '-S', sam] + list(args))
	with open(sam) as f:
		return h.parse_sam(f.read())[1], proc


def test_un_conc_al_conc_split(bt2, genome_dict, genome_index, tmp_path):
	"""--al-conc gets the concordant pair and --un-conc the rest (mixed, unaligned, discordant), each record exactly as in the input."""
	d = new_dir(tmp_path)
	m1, m2, text, classes = conc_split_input(genome_dict, d)
	out = os.path.join(d, 'out')
	os.mkdir(out)
	run_with_read_files(bt2, genome_index, d, m1, m2,
		['--al-conc', os.path.join(out, 'al.fq'), '--un-conc', os.path.join(out, 'un.fq')])
	assert sorted(os.listdir(out)) == ['al.1.fq', 'al.2.fq', 'un.1.fq', 'un.2.fq']
	for kind in ('al', 'un'):
		for mate in (1, 2):
			with open(os.path.join(out, '%s.%d.fq' % (kind, mate))) as f:
				got = fastq_records(f.read())
			assert sorted(got) == sorted(text[mate][i] for i in classes[kind])


def test_conc_outputs_leave_sam_unchanged(bt2, genome_dict, genome_index, tmp_path):
	"""Adding --un-conc and --al-conc does not change the SAM records."""
	d = new_dir(tmp_path)
	m1, m2, _, _ = conc_split_input(genome_dict, d)
	plain, _ = run_with_read_files(bt2, genome_index, d, m1, m2, [])
	split, _ = run_with_read_files(bt2, genome_index, d, m1, m2,
		['--al-conc', os.path.join(d, 'al.fq'), '--un-conc', os.path.join(d, 'un.fq')])
	assert [r.line for r in split] == [r.line for r in plain]
	assert [r.tags['YT'] for r in plain] == ['CP', 'CP', 'UP', 'UP', 'UP', 'UP', 'DP', 'DP']


@pytest.mark.parametrize('name,files', [
	('reads.fq', ['reads.1.fq', 'reads.2.fq']),
	('reads.x.fastq', ['reads.x.1.fastq', 'reads.x.2.fastq']),
	('pair_%_mates.fq', ['pair_1_mates.fq', 'pair_2_mates.fq']),
	('noext', ['noext.1', 'noext.2']),
])
@pytest.mark.parametrize('kind', ['al', 'un'])
def test_conc_filenames(bt2, genome_dict, genome_index, tmp_path, name, files, kind):
	"""--al-conc/--un-conc put .1/.2 before the final dot, or replace % with 1/2; mate 1s go to the first file (no dot: .1/.2 appended)."""
	d = new_dir(tmp_path)
	m1, m2, text, classes = conc_split_input(genome_dict, d)
	out = os.path.join(d, 'out')
	os.mkdir(out)
	run_with_read_files(bt2, genome_index, d, m1, m2, ['--%s-conc' % kind, os.path.join(out, name)])
	assert sorted(os.listdir(out)) == sorted(files)
	for mate, fn in zip((1, 2), files):
		with open(os.path.join(out, fn)) as f:
			assert sorted(fastq_records(f.read())) == sorted(text[mate][i] for i in classes[kind])


@pytest.mark.parametrize('comp', ['gz', 'bz2', 'zst', 'lz4'])
@pytest.mark.parametrize('kind', ['al', 'un'])
def test_conc_compressed(bt2, genome_dict, genome_index, tmp_path, comp, kind):
	"""--al-conc-<c>/--un-conc-<c> write compressed per-mate files (out.fq.<c> -> out.fq.1.<c>, out.fq.2.<c>) holding the right reads."""
	tool = {'gz': 'gzip', 'bz2': 'bzip2', 'zst': 'zstd', 'lz4': 'lz4'}[comp]
	if not shutil.which(tool):
		pytest.skip('%s not installed (the wrapper compresses through it)' % tool)
	d = new_dir(tmp_path)
	m1, m2, text, classes = conc_split_input(genome_dict, d)
	out = os.path.join(d, 'out')
	os.mkdir(out)
	run_with_read_files(bt2, genome_index, d, m1, m2,
		['--%s-conc-%s' % (kind, comp), os.path.join(out, 'out.fq.' + comp)])
	files = ['out.fq.1.' + comp, 'out.fq.2.' + comp]
	assert sorted(os.listdir(out)) == files
	for mate, fn in zip((1, 2), files):
		path = os.path.join(out, fn)
		with open(path, 'rb') as f:
			assert f.read(len(MAGIC[comp])) == MAGIC[comp]
		got = fastq_records(read_compressed(path, comp))
		assert sorted(got) == sorted(text[mate][i] for i in classes[kind])


def test_conc_compressed_percent_pattern(bt2, genome_dict, genome_index, tmp_path):
	"""The % pattern also applies to compressed output: un_%.fq.gz -> un_1.fq.gz, un_2.fq.gz."""
	if not shutil.which('gzip'):
		pytest.skip('gzip not installed (the wrapper compresses through it)')
	d = new_dir(tmp_path)
	m1, m2, text, classes = conc_split_input(genome_dict, d)
	out = os.path.join(d, 'out')
	os.mkdir(out)
	run_with_read_files(bt2, genome_index, d, m1, m2,
		['--un-conc-gz', os.path.join(out, 'un_%.fq.gz')])
	assert sorted(os.listdir(out)) == ['un_1.fq.gz', 'un_2.fq.gz']
	for mate in (1, 2):
		got = fastq_records(read_compressed(os.path.join(out, 'un_%d.fq.gz' % mate), 'gz'))
		assert sorted(got) == sorted(text[mate][i] for i in classes['un'])


def test_al_conc_writes_repeat_pair_once(bt2, genome_dict, repeat_index, tmp_path):
	"""With -k 3 a pair with two concordant placements has secondary records in SAM but appears once in --al-conc."""
	c2 = genome_dict['chr2']
	frag = Fragment(c2, 1020, 330, 50, chrom='chr2')
	d = new_dir(tmp_path)
	m1, m2 = write_pair_files(d, [frag.pair('rep')])
	records, _ = run_with_read_files(bt2, repeat_index, d, m1, m2,
		['-k', '3', '--al-conc', os.path.join(d, 'al.fq')])
	assert len(records) == 4
	assert sum(1 for r in records if r.is_secondary) == 2
	for mate, seq in ((1, frag.mate1), (2, frag.mate2)):
		got = h.read_fastq(os.path.join(d, 'al.%d.fq' % mate))
		assert got == [('rep', seq, 'I' * len(seq))]


def test_un_mates(bt2, genome_dict, genome_index, tmp_path):
	"""--un-mates (wrapper option, not in the MANUAL) writes each unaligned mate of a pair to its mate's file; pins the current wrapper behaviour."""
	d = new_dir(tmp_path)
	m1, m2, text, _ = conc_split_input(genome_dict, d)
	out = os.path.join(d, 'out')
	os.mkdir(out)
	run_with_read_files(bt2, genome_index, d, m1, m2, ['--un-mates', os.path.join(out, 'um.fq')])
	assert sorted(os.listdir(out)) == ['um.1.fq', 'um.2.fq']
	# Unaligned: mate 2 of m1only (index 1) and both mates of none (index 2).
	with open(os.path.join(out, 'um.1.fq')) as f:
		assert fastq_records(f.read()) == [text[1][2]]
	with open(os.path.join(out, 'um.2.fq')) as f:
		assert sorted(fastq_records(f.read())) == sorted([text[2][1], text[2][2]])


# ---------------------------------------------------------------------------
# Paired alignment summary
# ---------------------------------------------------------------------------

def summary_pairs(genome_dict, with_discordant=True):
	"""Pairs of known class against `repeat_index`: 3 concordant once, 2
	concordant twice (inside the duplicated region), 2 discordant (other
	chromosome; FF orientation), and 3 that are neither: unique mate + junk,
	junk + junk, repeat mate + unique mate on another chromosome."""
	c1, c2, c3 = genome_dict['chr1'], genome_dict['chr2'], genome_dict['chr3']
	pairs = [
		Fragment(c1, 100, 300, 50).pair('conc1a'),
		Fragment(c1, 1100, 250, 50).pair('conc1b'),
		Fragment(c3, 500, 250, 50).pair('conc1c'),
		Fragment(c2, 1020, 330, 50).pair('conc2a'),
		Fragment(c2, 1050, 250, 50, strand='-').pair('conc2b'),
	]
	if with_discordant:
		pairs += [
			('disc_chr', c1[2000:2050], h.revcomp(c3[1500:1550])),
			Fragment(c1, 3000, 250, 50, orient='ff').pair('disc_ff'),
		]
	pairs += [
		('one_mate', c1[4000:4050], junk(50, 6)),
		('no_mate', junk(50, 7), junk(50, 8)),
		('rep_mate', c2[1100:1150], h.revcomp(c3[100:150])),
	]
	return pairs


@pytest.mark.parametrize('args,with_disc,expect', [
	([], True, dict(conc=(5, 3, 2), disc=2, mates=(3, 2, 1))),
	(['--no-discordant'], True, dict(conc=(5, 3, 2), disc=0, mates=(3, 6, 1), discordant=False)),
	(['--no-mixed'], False, dict(conc=(3, 3, 2), disc=0, mates=None, mixed=False)),
	(['--no-mixed', '--no-discordant'], True,
		dict(conc=(5, 3, 2), disc=0, mates=None, mixed=False, discordant=False)),
])
def test_paired_summary(bt2, genome_dict, repeat_index, tmp_path, args, with_disc, expect):
	"""The paired alignment summary counts concordant 0/1/>1, discordant, and mates 0/1/>1 times as the input was built; --no-discordant/--no-mixed drop their sections.

	The MANUAL shows only the full layout, so the dropped sections pin the
	current output. --no-mixed is run without the discordant pairs, which the
	current version wrongly suppresses.
	"""
	pairs = summary_pairs(genome_dict, with_disc)
	records, proc = align_pairs(bt2, repeat_index, tmp_path, pairs, args)
	expected = paired_summary(len(pairs), **expect)
	assert summary_block(proc.stderr) == expected
	yt = dict((name, r1.tags['YT']) for name, (r1, r2) in mates_by_name(records).items())
	assert sorted(n for n in yt if yt[n] == 'CP') == ['conc1a', 'conc1b', 'conc1c', 'conc2a', 'conc2b']
	assert sorted(n for n in yt if yt[n] == 'DP') == (
		['disc_chr', 'disc_ff'] if with_disc and expect['disc'] else [])


def test_paired_summary_manual_example_layout(bt2, genome_dict, genome_index, tmp_path):
	"""The summary for an all-concordant input follows the MANUAL's layout line by line."""
	c = genome_dict['chr1']
	pairs = [Fragment(c, s, 300, 50).pair('p%d' % s) for s in (200, 1200, 2200, 3200)]
	_, proc = align_pairs(bt2, genome_index, tmp_path, pairs)
	assert summary_block(proc.stderr) == [
		'4 reads; of these:',
		'  4 (100.00%) were paired; of these:',
		'    0 (0.00%) aligned concordantly 0 times',
		'    4 (100.00%) aligned concordantly exactly 1 time',
		'    0 (0.00%) aligned concordantly >1 times',
		'    ----',
		'    0 pairs aligned concordantly 0 times; of these:',
		'      0 (0.00%) aligned discordantly 1 time',
		'    ----',
		'    0 pairs aligned 0 times concordantly or discordantly; of these:',
		'      0 mates make up the pairs; of these:',
		'        0 (0.00%) aligned 0 times',
		'        0 (0.00%) aligned exactly 1 time',
		'        0 (0.00%) aligned >1 times',
		'100.00% overall alignment rate']


# ---------------------------------------------------------------------------
# Unequal mate counts
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('short_mate', [1, 2])
def test_unequal_mate_counts(bt2, genome_dict, genome_index, tmp_path, short_mate):
	"""A -1/-2 file pair with different read counts is an error naming the shorter file."""
	c = genome_dict['chr1']
	frags = [Fragment(c, 300 * i, 250, 40) for i in range(1, 6)]
	d = new_dir(tmp_path)
	r1 = [('p%d' % i, f.mate1) for i, f in enumerate(frags)]
	r2 = [('p%d' % i, f.mate2) for i, f in enumerate(frags)]
	if short_mate == 1:
		r1 = r1[:3]
	else:
		r2 = r2[:3]
	m1 = h.write_fastq(os.path.join(d, 'r1.fq'), r1)
	m2 = h.write_fastq(os.path.join(d, 'r2.fq'), r2)
	proc = bt2.run('bowtie2', ['-x', genome_index, '-p', '1', '-1', m1, '-2', m2,
		'-S', os.path.join(d, 'out.sam')], check=False, cwd=d)
	assert proc.returncode != 0
	other = 3 - short_mate
	assert ('Error, fewer reads in file specified with -%d than in file specified with -%d'
		% (short_mate, other)) in proc.stderr


# ---------------------------------------------------------------------------
# Same pairing from --interleaved, -b --align-paired-reads, and file lists
# ---------------------------------------------------------------------------

def input_pairs(genome_dict):
	c1, c2 = genome_dict['chr1'], genome_dict['chr2']
	return [
		Fragment(c1, 1000, 300, 50).pair('conc'),
		Fragment(c1, 2000, 300, 50, 40, strand='-').pair('crick'),
		('m1only', c1[1500:1550], junk(50, 9)),
		('disc', c1[2600:2650], h.revcomp(c2[500:550])),
		('none', junk(45, 10), junk(45, 11)),
	]


def test_interleaved_and_bam_pair_like_separate_files(bt2, genome_dict, genome_index, tmp_path):
	"""--interleaved FASTQ and -b --align-paired-reads give the same SAM records as -1/-2 for the same pairs."""
	pairs = input_pairs(genome_dict)
	d = new_dir(tmp_path)
	m1, m2 = write_pair_files(d, pairs)
	il = h.write_fastq(os.path.join(d, 'il.fq'),
		[r for name, s1, s2 in pairs for r in ((name, s1), (name, s2))])
	bam_records = []
	for name, s1, s2 in pairs:
		bam_records.append((name, PAIRED | UNMAPPED | MATE_UNMAPPED | FIRST_IN_PAIR, s1, 'I' * len(s1)))
		bam_records.append((name, PAIRED | UNMAPPED | MATE_UNMAPPED | SECOND_IN_PAIR, s2, 'I' * len(s2)))
	bam = write_unaligned_bam(os.path.join(d, 'reads.bam'), bam_records)
	base = ['-x', genome_index, '-p', '1']
	_, sep, _ = bt2.align(base + ['-1', m1, '-2', m2])
	_, inter, _ = bt2.align(base + ['--interleaved', il])
	_, frombam, _ = bt2.align(base + ['-b', bam, '--align-paired-reads'])
	assert [r.tags['YT'] for r in sep] == ['CP'] * 4 + ['UP'] * 2 + ['DP'] * 2 + ['UP'] * 2
	assert [r.line for r in inter] == [r.line for r in sep]
	assert [r.line for r in frombam] == [r.line for r in sep]


def test_mate_file_lists_pair_in_order(bt2, genome_dict, genome_index, tmp_path):
	"""-1 a,b -2 c,d pairs a with c and b with d: same SAM records as one file per mate."""
	pairs = input_pairs(genome_dict)
	d1, d2, d3 = new_dir(tmp_path), new_dir(tmp_path), new_dir(tmp_path)
	m1, m2 = write_pair_files(d1, pairs)
	a1, a2 = write_pair_files(d2, pairs[:2])
	b1, b2 = write_pair_files(d3, pairs[2:])
	base = ['-x', genome_index, '-p', '1']
	_, single, _ = bt2.align(base + ['-1', m1, '-2', m2])
	_, listed, _ = bt2.align(base + ['-1', a1 + ',' + b1, '-2', a2 + ',' + b2])
	assert [r.qname for r in single] == [p[0] for p in pairs for _ in (1, 2)]
	assert [r.line for r in listed] == [r.line for r in single]
