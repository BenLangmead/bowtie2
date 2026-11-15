"""Alignment, scoring, reporting and effort options of bowtie2.

Reads are cut from fixed random references at known positions, with
substitutions, Ns, insertions and deletions placed at known offsets, so the
expected POS, strand, CIGAR, MD, NM, XM, XO, XG and AS of each alignment are
computed here from the MANUAL's scoring rules before bowtie2 runs:

- end-to-end: mismatch penalty MN + floor((MX - MN) * min(Q, 40) / 40) with
  MX = 6, MN = 2; N penalty 1; read gap (deletion, CIGAR D) and reference gap
  (insertion, CIGAR I) of length n cost 5 + 3n; the best score is 0.
- local: as above, plus a bonus of 2 for each matching base.

Repeats come from a second, small reference (`rep`) with a segment present
three times exactly and once with two substitutions, and from a third
(`family`) holding 800 slightly diverged copies of one element, used where an
option only changes behaviour for reads with repetitive seeds (-D, -R).
"""

import random
import re

import pytest

import bt2_helpers as h


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_reads(tmp_path, reads, name='reads.fq'):
	return h.write_fastq(str(tmp_path / name), reads)


def align(bt2, index, fq, args=(), check=True):
	"""Align unpaired FASTQ `fq` with -p 1; returns the SAM records."""
	_, records, _ = bt2.align(['-p', '1', '-x', index, '-U', fq] + [str(a) for a in args],
		check=check)
	return records


def body(records):
	"""SAM lines without the header (whose @PG line holds the command line)."""
	return [r.line for r in records]


def same_records(a, b, nreads):
	"""Records a and b are identical and cover all `nreads` reads, so the
	comparison cannot pass on two empty or truncated outputs."""
	assert len(set(r.qname for r in a)) == nreads
	assert body(a) == body(b)


def one_per_read(records):
	"""Map qname -> its single record; fails if a read has several."""
	out = {}
	for r in records:
		assert r.qname not in out, 'more than one record for %s' % r.qname
		out[r.qname] = r
	return out


def other_base(*avoid):
	return [b for b in 'ACGT' if b not in avoid][0]


def gaps_unambiguous(ref, start, ops):
	"""True if every gap in `ops`, applied at `start`, has one placement only
	(it cannot slide left or right through equal bases)."""
	pos = start
	for op in ops:
		kind = op[0]
		if kind == 'D':
			n = op[1]
			if ref[pos - 1] == ref[pos + n - 1] or ref[pos] == ref[pos + n]:
				return False
			pos += n
		elif kind == 'I':
			pass  # Edited chooses inserted bases that cannot slide
		elif kind == '=':
			pos += op[1]
		else:
			pos += 1
	return True


def place(ref, start, ops):
	"""First offset >= start where the gaps in `ops` are unambiguous."""
	while not gaps_unambiguous(ref, start, ops):
		start += 1
	return start


class Edited(object):
	"""A read built from a reference and an edit script, with the alignment
	fields bowtie2 should report for it on the forward strand."""

	def __init__(self, ref, start, ops):
		seq, cigar, md = [], [], []
		self.sub_offsets = []    # read offsets of substitutions
		self.n_offsets = []      # read offsets of read Ns
		self.dels, self.inss = [], []
		self.matches = 0
		run, pos = 0, start

		def add(op, n):
			if cigar and cigar[-1][1] == op:
				cigar[-1] = (cigar[-1][0] + n, op)
			else:
				cigar.append((n, op))

		for op in ops:
			kind = op[0]
			if kind == '=':
				n = op[1]
				seq.append(ref[pos:pos + n])
				add('M', n)
				run += n
				pos += n
				self.matches += n
			elif kind in ('X', 'N'):
				offset = sum(len(s) for s in seq)
				if kind == 'X':
					seq.append(h.mutate(ref[pos], [0]))
					self.sub_offsets.append(offset)
				else:
					seq.append('N')
					self.n_offsets.append(offset)
				add('M', 1)
				md.append('%d%s' % (run, ref[pos]))
				run = 0
				pos += 1
			elif kind == 'D':
				n = op[1]
				add('D', n)
				md.append('%d^%s' % (run, ref[pos:pos + n]))
				run = 0
				pos += n
				self.dels.append(n)
			elif kind == 'I':
				# The first inserted base differs from the next reference base
				# and the last from the previous one, so the gap cannot slide.
				n = op[1]
				if n == 1:
					seq.append(other_base(ref[pos - 1], ref[pos]))
				else:
					seq.append(other_base(ref[pos]) + 'A' * (n - 2) + other_base(ref[pos - 1]))
				add('I', n)
				self.inss.append(n)
			else:
				raise ValueError(op)
		md.append(str(run))
		self.seq = ''.join(seq)
		self.pos = start + 1
		self.cigar = ''.join('%d%s' % (n, op) for n, op in cigar)
		self.md = ''.join(md)
		self.xm = len(self.sub_offsets) + len(self.n_offsets)
		self.xo = len(self.dels) + len(self.inss)
		self.xg = sum(self.dels) + sum(self.inss)
		self.nm = self.xm + self.xg

	def penalty(self, qual, mx=6, mn=2, np=1, rdg=(5, 3), rfg=(5, 3), ignore_quals=False):
		"""Sum of the MANUAL's penalties for this read's edits."""
		pen = 0
		for off in self.sub_offsets:
			pen += mx if ignore_quals else h.mismatch_penalty(qual[off], mx=mx, mn=mn)
		pen += np * len(self.n_offsets)
		pen += sum(rdg[0] + rdg[1] * n for n in self.dels)
		pen += sum(rfg[0] + rfg[1] * n for n in self.inss)
		return pen

	def e2e_score(self, qual, **kw):
		return -self.penalty(qual, **kw)

	def local_score(self, qual, ma=2, **kw):
		return ma * self.matches - self.penalty(qual, **kw)


def fastq_record(name, edited, qual, reverse):
	"""(name, seq, qual) as sequenced from the forward or reverse strand."""
	if reverse:
		return (name, h.revcomp(edited.seq), qual[::-1])
	return (name, edited.seq, qual)


def check_fields(rec, rname, edited, qual, reverse, score):
	"""Assert every alignment field of `rec` against the read's construction."""
	assert not rec.is_unmapped, rec
	assert rec.flag == (h.REVERSE if reverse else 0), rec
	assert rec.rname == rname
	assert rec.pos == edited.pos
	assert rec.cigar == edited.cigar
	assert rec.seq == edited.seq
	assert rec.qual == qual
	assert rec.tags['AS'] == score
	assert rec.tags['MD'] == edited.md
	assert rec.tags['NM'] == edited.nm
	assert rec.tags['XM'] == edited.xm
	assert rec.tags['XO'] == edited.xo
	assert rec.tags['XG'] == edited.xg
	assert rec.tags['XN'] == 0
	assert rec.tags['YT'] == 'UU'
	assert 'XS' not in rec.tags
	assert 'YF' not in rec.tags


def qual_with(length, changes):
	q = ['I'] * length
	for off, c in changes.items():
		q[off] = c
	return ''.join(q)


# ---------------------------------------------------------------------------
# Reads with known edits against the shared `genome`
# ---------------------------------------------------------------------------

# name: (chromosome, approximate start, edit script, {read offset: quality})
EDIT_READS = {
	'perfect': ('chr1', 100, [('=', 100)], {}),
	# Q40 ('I'), Q20 ('5') and Q10 ('+') substitutions: penalties 6, 4, 3.
	'subs': ('chr1', 1000, [('=', 20), ('X',), ('=', 39), ('X',), ('=', 24), ('X',),
		('=', 14)], {60: '5', 85: '+'}),
	'del3': ('chr1', 1500, [('=', 50), ('D', 3), ('=', 50)], {}),
	'ins2': ('chr2', 1000, [('=', 40), ('I', 2), ('=', 58)], {}),
	'combo': ('chr2', 2000, [('=', 8), ('X',), ('=', 16), ('X',), ('=', 20), ('D', 2),
		('=', 46), ('I', 1), ('=', 7)], {25: '5'}),
	'ns': ('chr3', 200, [('=', 30), ('N',), ('=', 19), ('N',), ('=', 49)], {}),
}


@pytest.fixture(scope='module')
def edit_reads(genome_dict):
	"""name -> (rname, Edited, qual) for EDIT_READS, gaps placed unambiguously."""
	out = {}
	for name, (rname, start, ops, quals) in EDIT_READS.items():
		ref = genome_dict[rname]
		edited = Edited(ref, place(ref, start, ops), ops)
		out[name] = (rname, edited, qual_with(len(edited.seq), quals))
	return out


@pytest.fixture(scope='module')
def edit_fastq(edit_reads, tmp_path_factory):
	"""FASTQ of every EDIT_READS read, from both strands (suffix _fw / _rc)."""
	reads = []
	for name, (rname, edited, qual) in sorted(edit_reads.items()):
		for strand, rev in (('fw', False), ('rc', True)):
			reads.append(fastq_record('%s_%s' % (name, strand), edited, qual, rev))
	path = tmp_path_factory.mktemp('edit_reads') / 'reads.fq'
	return h.write_fastq(str(path), reads)


@pytest.fixture(scope='module')
def edit_e2e(bt2, genome_index, edit_fastq):
	return one_per_read(align(bt2, genome_index, edit_fastq))


@pytest.mark.parametrize('name', sorted(EDIT_READS))
@pytest.mark.parametrize('strand', ['fw', 'rc'])
def test_fields_end_to_end(edit_reads, edit_e2e, name, strand):
	"""POS, FLAG, CIGAR, SEQ, QUAL, AS, MD, NM, XM, XO, XG, XN, YT match the read's construction."""
	rname, edited, qual = edit_reads[name]
	rec = edit_e2e['%s_%s' % (name, strand)]
	check_fields(rec, rname, edited, qual, strand == 'rc', edited.e2e_score(qual))


def test_fields_expected_values_spelled_out(edit_reads, edit_e2e):
	"""The combined read's fields, written out literally as a check on the helper."""
	rname, edited, qual = edit_reads['combo']
	rec = edit_e2e['combo_fw']
	# Q40 substitution (6) + Q20 substitution (4) + 2 bp deletion (5 + 2*3)
	# + 1 bp insertion (5 + 3).
	assert rec.tags['AS'] == -(6 + 4 + 11 + 8)
	assert rec.cigar == '46M2D46M1I7M'
	assert (rec.tags['XM'], rec.tags['XO'], rec.tags['XG'], rec.tags['NM']) == (2, 2, 3, 5)
	assert re.match(r'^8[ACGT]16[ACGT]20\^[ACGT]{2}53$', rec.tags['MD']), rec.tags['MD']


def test_end_to_end_is_default(bt2, genome_index, edit_fastq):
	"""--end-to-end gives the same SAM records as the default mode."""
	same_records(align(bt2, genome_index, edit_fastq, ['--end-to-end']),
		align(bt2, genome_index, edit_fastq), 2 * len(EDIT_READS))


def test_end_to_end_never_clips(edit_e2e):
	"""In --end-to-end mode every aligned read is aligned over its full length."""
	assert len(edit_e2e) == 2 * len(EDIT_READS)
	for rec in edit_e2e.values():
		assert 'S' not in rec.cigar
		assert h.cigar_query_length(rec.cigar) == len(rec.seq)
		assert rec.tags['AS'] <= 0


def test_ns_in_read_count_as_mismatches(edit_e2e):
	"""Read Ns cost --np each, appear in MD as the reference base and count in XM (pins current XM)."""
	rec = edit_e2e['ns_fw']
	assert re.match(r'^30[ACGT]19[ACGT]49$', rec.tags['MD']), rec.tags['MD']
	assert rec.tags['AS'] == -2
	assert rec.tags['XM'] == 2
	assert rec.tags['NM'] == 2


# ---------------------------------------------------------------------------
# Local mode
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('ma', [None, 3])
def test_local_scores_whole_reads(bt2, genome_index, edit_reads, edit_fastq, ma):
	"""--local (and --ma) adds the match bonus per matching base and keeps reads that align end to end unclipped."""
	args = ['--local'] + (['--ma', ma] if ma else [])
	recs = one_per_read(align(bt2, genome_index, edit_fastq, args))
	for name, (rname, edited, qual) in edit_reads.items():
		for strand in ('fw', 'rc'):
			rec = recs['%s_%s' % (name, strand)]
			check_fields(rec, rname, edited, qual, strand == 'rc',
				edited.local_score(qual, ma=ma or 2))


@pytest.fixture(scope='module')
def clip_reads(genome_dict):
	"""Reads whose 10-base ends differ from the reference at every base,
	around an 80-base core (exact, or with one Q40 substitution)."""
	ref = genome_dict['chr3']
	s = 1000
	left = h.mutate(ref[s - 10:s], range(10))
	right = h.mutate(ref[s + 80:s + 90], range(10))
	exact = Edited(ref, s, [('=', 80)])
	sub = Edited(ref, s, [('=', 40), ('X',), ('=', 39)])
	return {
		'clip': (left + exact.seq + right, exact),
		'clip_sub': (left + sub.seq + right, sub),
	}


@pytest.mark.parametrize('strand', ['fw', 'rc'])
def test_local_soft_clips_mismatched_ends(bt2, genome_index, clip_reads, tmp_path, strand):
	"""--local soft clips ends that do not match: CIGAR 10S..10S, POS at the core, AS = 2 * matches - penalties."""
	reads = []
	for name, (seq, _) in sorted(clip_reads.items()):
		reads.append((name, h.revcomp(seq) if strand == 'rc' else seq))
	fq = write_reads(tmp_path, reads)
	recs = one_per_read(align(bt2, genome_index, fq, ['--local']))
	for name, (seq, core) in clip_reads.items():
		rec = recs[name]
		assert rec.flag == (h.REVERSE if strand == 'rc' else 0)
		assert (rec.rname, rec.pos) == ('chr3', core.pos)
		assert rec.cigar == '10S%s10S' % core.cigar
		assert rec.seq == seq
		assert rec.tags['MD'] == core.md
		assert rec.tags['NM'] == core.nm
		assert rec.tags['AS'] == core.local_score('I' * 100)
	assert recs['clip'].tags['AS'] == 160
	assert recs['clip_sub'].tags['AS'] == 2 * 79 - 6
	# The same reads cannot align end to end: 20 mismatched ends cost far
	# more than the default minimum score allows.
	e2e = align(bt2, genome_index, fq)
	assert len(e2e) == 2
	for rec in e2e:
		assert rec.is_unmapped


@pytest.mark.parametrize('order', ['before', 'after'])
@pytest.mark.parametrize('preset', ['very-fast', 'fast', 'sensitive', 'very-sensitive'])
def test_local_with_preset_is_local_preset(bt2, lambda_index, example_reads, preset, order):
	"""--local combined with an end-to-end preset name equals the -local preset (MANUAL, --local)."""
	fq = example_reads['reads_1']
	args = ['--local', '--' + preset] if order == 'before' else ['--' + preset, '--local']
	same_records(align(bt2, lambda_index, fq, args + ['-u', 1000]),
		align(bt2, lambda_index, fq, ['--%s-local' % preset, '-u', 1000]), 1000)


# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------

PRESETS = {
	'very-fast': '-D 5 -R 1 -N 0 -L 22 -i S,0,2.50',
	'fast': '-D 10 -R 2 -N 0 -L 22 -i S,0,2.50',
	'sensitive': '-D 15 -R 2 -N 0 -L 22 -i S,1,1.15',
	'very-sensitive': '-D 20 -R 3 -N 0 -L 20 -i S,1,0.50',
	'very-fast-local': '-D 5 -R 1 -N 0 -L 25 -i S,1,2.00',
	'fast-local': '-D 10 -R 2 -N 0 -L 22 -i S,1,1.75',
	'sensitive-local': '-D 15 -R 2 -N 0 -L 20 -i S,1,0.75',
	'very-sensitive-local': '-D 20 -R 3 -N 0 -L 20 -i S,1,0.50',
}


def family(seed, elen, ncopies, ndiv, spacer):
	"""One random element, `ncopies` copies each with `ndiv` substitutions,
	separated by random spacers. Returns (sequence, copy start offsets)."""
	rng = random.Random(seed)
	element = h.random_dna(elen, rng)
	parts, origins, pos = [], [], 0
	for _ in range(ncopies):
		parts.append(h.random_dna(spacer, rng))
		pos += spacer
		parts.append(h.mutate(element, sorted(rng.sample(range(elen), ndiv))))
		origins.append(pos)
		pos += elen
	parts.append(h.random_dna(spacer, rng))
	return ''.join(parts), origins


@pytest.fixture(scope='module')
def family_index(bt2, tmp_path_factory):
	"""Index of 800 copies of a 200 bp element, 2 substitutions per copy."""
	seq, _ = family(11, 200, 800, 2, 20)
	d = tmp_path_factory.mktemp('family')
	fa = h.write_fasta(str(d / 'family.fa'), [('fam', seq)])
	prefix = str(d / 'family')
	bt2.build_index(fa, prefix, extra=['--seed', '0'])
	return prefix


@pytest.fixture(scope='module')
def family_fastq(tmp_path_factory):
	"""300 reads of 100 bp from random copies, 3 substitutions each, either strand.
	Their 22-mer seeds hit hundreds of copies, so they are repetitive reads."""
	seq, origins = family(11, 200, 800, 2, 20)
	rng = random.Random(3)
	reads = []
	for i in range(300):
		s = rng.choice(origins) + rng.randrange(0, 100)
		r = h.mutate(seq[s:s + 100], sorted(rng.sample(range(100), 3)))
		if rng.random() < 0.5:
			r = h.revcomp(r)
		reads.append(('f%d' % i, r))
	path = tmp_path_factory.mktemp('family_reads') / 'reads.fq'
	return h.write_fastq(str(path), reads)


@pytest.fixture(scope='module')
def preset_inputs(lambda_index, example_reads, family_index, family_fastq):
	"""(index, reads, extra args, number of reads): the shipped lambda reads
	(non-repetitive, sensitive to -N/-L/-i) and the family reads (sensitive
	to -D/-R)."""
	return {
		'lambda': (lambda_index, example_reads['reads_1'], ['-u', 3000], 3000),
		'family': (family_index, family_fastq, [], 300),
	}


@pytest.mark.parametrize('dataset', ['lambda', 'family'])
@pytest.mark.parametrize('preset', sorted(PRESETS))
def test_preset_equals_expansion(bt2, preset_inputs, preset, dataset):
	"""Each preset gives the same SAM records as its documented -D/-R/-N/-L/-i expansion."""
	index, fq, extra, nreads = preset_inputs[dataset]
	mode = ['--local'] if preset.endswith('-local') else []
	same_records(align(bt2, index, fq, ['--' + preset] + extra),
		align(bt2, index, fq, mode + PRESETS[preset].split() + extra), nreads)


@pytest.mark.parametrize('mode', ['end-to-end', 'local'])
def test_presets_differ_from_each_other(bt2, preset_inputs, mode):
	"""The four presets of a mode give four different outputs, so the equivalence tests are not vacuous."""
	index, fq, extra, _ = preset_inputs['family']
	names = [p for p in PRESETS if p.endswith('-local') == (mode == 'local')]
	outputs = set()
	for p in names:
		outputs.add(tuple(body(align(bt2, index, fq, ['--' + p] + extra))))
	assert len(outputs) == 4


@pytest.mark.parametrize('mode,preset', [('end-to-end', 'sensitive'), ('local', 'sensitive-local')])
def test_default_preset_is_sensitive(bt2, preset_inputs, mode, preset):
	"""With no preset, --end-to-end uses --sensitive and --local uses --sensitive-local."""
	for index, fq, extra, nreads in preset_inputs.values():
		same_records(align(bt2, index, fq, ['--' + mode] + extra),
			align(bt2, index, fq, ['--' + preset] + extra), nreads)


@pytest.mark.parametrize('preset', ['very-fast', 'very-sensitive', 'fast-local'])
def test_preset_option_equals_preset_flag(bt2, preset_inputs, preset):
	"""The undocumented --preset <name> equals --<name> (with --local for a -local name)."""
	index, fq, extra, nreads = preset_inputs['family']
	mode = ['--local'] if preset.endswith('-local') else []
	same_records(align(bt2, index, fq, mode + ['--preset', preset] + extra),
		align(bt2, index, fq, ['--' + preset] + extra), nreads)


# ---------------------------------------------------------------------------
# Seed options: -N, -L, -i, --no-1mm-upfront
# ---------------------------------------------------------------------------

def mapped(bt2, index, tmp_path, seq, args, qual=None):
	"""Align one read; returns its record (mapped or not)."""
	fq = write_reads(tmp_path, [('r', seq, qual)])
	recs = align(bt2, index, fq, args)
	assert len(recs) == 1
	return recs[0]


def test_seed_mismatches_N(bt2, genome_index, genome_dict, tmp_path):
	"""-N 1 finds a read whose every seed has one mismatch; -N 0 cannot."""
	# 50 bp, seeds of 20 every 20 bases (offsets 0 and 20), one substitution
	# in each seed; two substitutions in total, so no 1-mismatch end-to-end hit.
	edited = Edited(genome_dict['chr1'], 2000, [('=', 10), ('X',), ('=', 19), ('X',), ('=', 19)])
	seeds = ['-L', 20, '-i', 'C,20,0']
	assert mapped(bt2, genome_index, tmp_path, edited.seq, seeds + ['-N', 0]).is_unmapped
	rec = mapped(bt2, genome_index, tmp_path, edited.seq, seeds + ['-N', 1])
	assert (rec.rname, rec.pos, rec.cigar) == ('chr1', edited.pos, '50M')
	assert rec.tags['AS'] == -12


def test_seed_mismatches_N_range(bt2, genome_index, genome_dict, tmp_path):
	"""-N accepts only 0 or 1."""
	fq = write_reads(tmp_path, [('r', genome_dict['chr1'][:50])])
	proc = bt2.run('bowtie2', ['-x', genome_index, '-U', fq, '-N', 2], check=False)
	assert proc.returncode != 0
	assert '-N argument must be within the interval [0,1]' in proc.stderr


@pytest.mark.parametrize('length,aligns', [(10, True), (14, True), (15, False), (22, False)])
def test_seed_length_L(bt2, genome_index, genome_dict, tmp_path, length, aligns):
	"""With seeds at every offset, a read aligns iff -L <= its longest exact stretch (14)."""
	# 60 bp with substitutions every 15 bases: exact runs of 7, 14, 14, 14, 7.
	edited = Edited(genome_dict['chr1'], 2100, [('=', 7), ('X',), ('=', 14), ('X',), ('=', 14),
		('X',), ('=', 14), ('X',), ('=', 7)])
	rec = mapped(bt2, genome_index, tmp_path, edited.seq, ['-L', length, '-i', 'C,1,0'])
	assert rec.is_unmapped != aligns
	if aligns:
		assert (rec.pos, rec.tags['AS'], rec.tags['MD']) == (edited.pos, -24, edited.md)


@pytest.mark.parametrize('ival,aligns', [
	('C,1,0', True), ('C,5,0', True), ('C,30,0', False),
	('C,0,0', True),        # below 1, rounded up to 1
	('L,0,0.1', True), ('L,0,0.5', False),      # 6, 30
	('S,0,0.6', True), ('S,0,4', False),        # 4.6, 31.0
	('G,0,1', True), ('G,0,8', False),          # 4.1, 32.8
])
def test_seed_interval_i(bt2, genome_index, genome_dict, tmp_path, ival, aligns):
	"""-i sets the seed spacing: only the exact 20-mers at offsets 3-8 can seed this 60 bp read."""
	# Substitutions at 2, 28, 40 and 55: the only exact 20-base windows start
	# at offsets 3..8. Spacing <= 8 puts a seed there; spacing 30 or more puts
	# seeds only at 0 and 30, which both hold a substitution.
	edited = Edited(genome_dict['chr1'], 2200, [('=', 2), ('X',), ('=', 25), ('X',), ('=', 11),
		('X',), ('=', 14), ('X',), ('=', 4)])
	rec = mapped(bt2, genome_index, tmp_path, edited.seq, ['-L', 20, '-i', ival])
	assert rec.is_unmapped != aligns
	if aligns:
		assert (rec.pos, rec.tags['AS']) == (edited.pos, -24)


def test_no_1mm_upfront(bt2, genome_index, genome_dict, tmp_path):
	"""With -N 0 and -L equal to the read length, a 1-mismatch read aligns only via the 1-mismatch upfront search."""
	edited = Edited(genome_dict['chr3'], 1500, [('=', 16), ('X',), ('=', 15)])
	rec = mapped(bt2, genome_index, tmp_path, edited.seq, ['-N', 0, '-L', 32])
	assert (rec.pos, rec.cigar, rec.tags['AS'], rec.tags['MD']) == (edited.pos, '32M', -6, edited.md)
	assert mapped(bt2, genome_index, tmp_path, edited.seq,
		['-N', 0, '-L', 32, '--no-1mm-upfront']).is_unmapped
	rec = mapped(bt2, genome_index, tmp_path, edited.seq, ['-N', 1, '-L', 32, '--no-1mm-upfront'])
	assert (rec.pos, rec.tags['AS']) == (edited.pos, -6)


# ---------------------------------------------------------------------------
# N ceiling and N penalty
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def three_n_read(genome_dict):
	"""50 bp with Ns at offsets 4, 9 and 14 (seeds after offset 15 are clean)."""
	return Edited(genome_dict['chr1'], 3000, [('=', 4), ('N',), ('=', 4), ('N',), ('=', 4), ('N',),
		('=', 35)])


@pytest.mark.parametrize('func,filtered', [
	(None, False),                       # default L,0,0.15 -> 7.5
	('L,0,0.07', False), ('L,0,0.05', True),     # 3.5, 2.5
	('C,3,0', False), ('C,2,0', True),
	('S,0,0.45', False), ('S,0,0.4', True),     # 3.18, 2.83
	('G,0,0.8', False), ('G,0,0.7', True),      # 3.13, 2.74
])
def test_n_ceil(bt2, genome_index, three_n_read, tmp_path, func, filtered):
	"""A read with 3 Ns is filtered (YF:Z:NS) iff the --n-ceil function is below 3."""
	args = ['--n-ceil', func] if func else []
	rec = mapped(bt2, genome_index, tmp_path, three_n_read.seq, args)
	if filtered:
		assert rec.is_unmapped
		assert rec.tags['YF'] == 'NS'
		assert rec.tag_order == ['YT', 'YF']   # field order pins current output
	else:
		assert (rec.rname, rec.pos, rec.cigar) == ('chr1', three_n_read.pos, '50M')
		assert rec.tags['MD'] == three_n_read.md
		assert rec.tags['AS'] == -3
		assert 'YF' not in rec.tags


@pytest.mark.parametrize('np', [0, 1, 3])
def test_np(bt2, genome_index, three_n_read, tmp_path, np):
	"""--np sets the penalty per read N."""
	args = ['--np', np] if np != 1 else []
	rec = mapped(bt2, genome_index, tmp_path, three_n_read.seq, args)
	assert rec.tags['AS'] == -3 * np


@pytest.fixture(scope='module')
def ref_n_index(bt2, tmp_path_factory):
	"""Index of one 1000 bp sequence with NN at 500-501 and 10 Ns at 700-709."""
	rng = random.Random(5)
	seq = list(h.random_dna(1000, rng))
	seq[500:502] = ['N'] * 2
	seq[700:710] = ['N'] * 10
	seq = ''.join(seq)
	d = tmp_path_factory.mktemp('ref_n')
	prefix = str(d / 'refn')
	bt2.build_index(h.write_fasta(str(d / 'refn.fa'), [('chrN', seq)]), prefix, extra=['--seed', '0'])
	return prefix, seq


@pytest.mark.parametrize('local', [False, True])
def test_reference_ns_xn(bt2, ref_n_index, tmp_path, local):
	"""XN counts reference Ns under the alignment; each costs --np and shows as N in MD."""
	prefix, seq = ref_n_index
	reads = [
		('n2', seq[450:550].replace('N', 'A')),
		('n10', seq[650:760].replace('N', 'C')),
		('n0', seq[300:400]),
		('n0_adjacent', h.mutate(seq[600:695], [20, 50, 80])),
	]
	recs = one_per_read(align(bt2, prefix, write_reads(tmp_path, reads),
		['--local'] if local else []))
	for name, read in reads:
		rec = recs[name]
		start = {'n2': 450, 'n10': 650, 'n0': 300, 'n0_adjacent': 600}[name]
		nref = seq[start:start + len(read)].count('N')
		nsub = 3 if name == 'n0_adjacent' else 0
		assert (rec.pos, rec.cigar) == (start + 1, '%dM' % len(read))
		assert rec.tags['XN'] == nref
		assert rec.tags['MD'].count('N') == nref
		pen = nref + 6 * nsub
		assert rec.tags['AS'] == (2 * (len(read) - nref - nsub) - pen if local else -pen)
	assert recs['n10'].tags['MD'] == '50' + 'N0' * 9 + 'N50'


# ---------------------------------------------------------------------------
# Gaps: --rdg, --rfg, --dpad, --gbar
# ---------------------------------------------------------------------------

def test_rdg_rfg(bt2, genome_index, edit_reads, edit_fastq):
	"""--rdg prices deletions (read gaps, CIGAR D) and --rfg insertions (reference gaps, CIGAR I) as open + n * extend."""
	recs = one_per_read(align(bt2, genome_index, edit_fastq, ['--rdg', '4,2', '--rfg', '7,1']))
	for name in ('del3', 'ins2', 'combo'):
		rname, edited, qual = edit_reads[name]
		for strand in ('fw', 'rc'):
			check_fields(recs['%s_%s' % (name, strand)], rname, edited, qual, strand == 'rc',
				edited.e2e_score(qual, rdg=(4, 2), rfg=(7, 1)))
	assert recs['del3_fw'].tags['AS'] == -(4 + 2 * 3)
	assert recs['ins2_fw'].tags['AS'] == -(7 + 1 * 2)


@pytest.fixture(scope='module')
def long_deletion(genome_dict):
	"""100 bp read with a 20 bp deletion after 50 bases."""
	ops = [('=', 50), ('D', 20), ('=', 50)]
	ref = genome_dict['chr2']
	return Edited(ref, place(ref, 1000, ops), ops)


def test_long_deletion_needs_cheap_gaps(bt2, genome_index, long_deletion, tmp_path):
	"""A 20 bp deletion costs 5 + 60 by default, below the minimum score of a 100 bp read; --rdg 1,1 makes it 21."""
	assert mapped(bt2, genome_index, tmp_path, long_deletion.seq, []).is_unmapped
	rec = mapped(bt2, genome_index, tmp_path, long_deletion.seq, ['--rdg', '1,1'])
	assert (rec.pos, rec.cigar, rec.tags['MD']) == \
		(long_deletion.pos, '50M20D50M', long_deletion.md)
	assert rec.tags['AS'] == -21
	assert (rec.tags['XO'], rec.tags['XG'], rec.tags['NM']) == (1, 20, 20)


@pytest.mark.parametrize('dpad,aligns', [(0, False), (5, False), (20, True)])
def test_dpad(bt2, genome_index, long_deletion, tmp_path, dpad, aligns):
	"""--dpad bounds the gaps dynamic programming can open: none with 0 columns, a 20 bp deletion with 20."""
	rec = mapped(bt2, genome_index, tmp_path, long_deletion.seq, ['--rdg', '1,1', '--dpad', dpad])
	assert rec.is_unmapped != aligns
	if aligns:
		assert (rec.pos, rec.cigar, rec.tags['AS']) == (long_deletion.pos, '50M20D50M', -21)


def e2e_score_from_fields(rec):
	"""End-to-end AS implied by a record's XM, CIGAR gaps (default penalties,
	all bases Q40)."""
	gaps = [n for n, op in h.cigar_ops(rec.cigar) if op in 'ID']
	return -(6 * rec.tags['XM'] + sum(5 + 3 * n for n in gaps))


@pytest.mark.parametrize('end', ['start', 'end'])
def test_gbar(bt2, genome_index, genome_dict, tmp_path, end):
	"""A 2 bp deletion 8 bases from a read end is used by default (--gbar 4) but not with --gbar 10, which forces a costlier alignment with no gap near the ends."""
	if end == 'start':
		ops = [('=', 8), ('D', 2), ('=', 92)]
	else:
		ops = [('=', 92), ('D', 2), ('=', 8)]
	ref = genome_dict['chr3']
	edited = Edited(ref, place(ref, 600, ops), ops)
	rec = mapped(bt2, genome_index, tmp_path, edited.seq, [])
	assert (rec.pos, rec.cigar, rec.tags['AS'], rec.tags['MD']) == \
		(edited.pos, edited.cigar, -11, edited.md)
	rec = mapped(bt2, genome_index, tmp_path, edited.seq, ['--gbar', 10])
	assert not rec.is_unmapped
	assert rec.rname == 'chr3' and abs(rec.pos - edited.pos) <= 2
	assert h.cigar_query_length(rec.cigar) == 100
	assert rec.tags['AS'] < -11
	assert rec.tags['AS'] == e2e_score_from_fields(rec)
	if end == 'end':
		# No reported gap within 10 bases of either read end. Not checked for
		# 'start': there the current version reports 9M2D91M, a gap the
		# barrier forbids.
		off = 0
		for n, op in h.cigar_ops(rec.cigar):
			if op in 'ID':
				assert 10 <= off <= 90, rec.cigar
			if op in 'MIS=X':
				off += n


# ---------------------------------------------------------------------------
# Mismatch scoring: --mp, --ignore-quals
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def qual_read(genome_dict):
	"""60 bp, substitutions at a Q10 base (15), a Q20 base (45) and a Q40 base (52)."""
	edited = Edited(genome_dict['chr3'], 500, [('=', 15), ('X',), ('=', 29), ('X',), ('=', 6),
		('X',), ('=', 7)])
	return edited, qual_with(60, {15: '+', 45: '5'})


@pytest.mark.parametrize('args,mx,mn,ignore', [
	([], 6, 2, False),
	(['--ignore-quals'], 6, 2, True),
	(['--mp', '4'], 4, 2, False),            # MX only (not in the MANUAL); MN stays 2
	(['--mp', '8,4'], 8, 4, False),
	(['--mp', '8,4', '--ignore-quals'], 8, 4, True),
	(['--mp', '5,5'], 5, 5, False),
])
def test_mismatch_penalty(bt2, genome_index, qual_read, tmp_path, args, mx, mn, ignore):
	"""Mismatch penalty is MN + floor((MX - MN) * min(Q, 40) / 40), or MX with --ignore-quals."""
	edited, qual = qual_read
	rec = mapped(bt2, genome_index, tmp_path, edited.seq, args, qual=qual)
	assert (rec.pos, rec.tags['MD'], rec.tags['XM']) == (edited.pos, edited.md, 3)
	assert rec.tags['AS'] == edited.e2e_score(qual, mx=mx, mn=mn, ignore_quals=ignore)


def test_mismatch_penalty_spelled_out(bt2, genome_index, qual_read, tmp_path):
	"""Default penalties at Q10, Q20 and Q40 are 3, 4 and 6; with --mp 8,4 they are 5, 6 and 8."""
	edited, qual = qual_read
	assert mapped(bt2, genome_index, tmp_path, edited.seq, [], qual=qual).tags['AS'] == -(3 + 4 + 6)
	assert mapped(bt2, genome_index, tmp_path, edited.seq, ['--mp', '8,4'],
		qual=qual).tags['AS'] == -(5 + 6 + 8)


def test_mp_max_below_min_is_an_error(bt2, genome_index, qual_read, tmp_path):
	"""--mp with MX < MN is rejected."""
	edited, qual = qual_read
	fq = write_reads(tmp_path, [('r', edited.seq, qual)])
	proc = bt2.run('bowtie2', ['-x', genome_index, '-U', fq, '--mp', '3,5'], check=False)
	assert proc.returncode != 0
	assert 'Maximum mismatch penalty (3) is less than minimum penalty (5' in proc.stderr


# ---------------------------------------------------------------------------
# --score-min, --ma and the YF filters
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def score18_read(genome_dict):
	"""64 bp with three Q40 substitutions (AS -18); seeds at 10 and 20 are clean."""
	return Edited(genome_dict['chr2'], 2500, [('=', 5), ('X',), ('=', 39), ('X',), ('=', 12),
		('X',), ('=', 5)])


@pytest.mark.parametrize('func,aligns', [
	('C,-18,0', True), ('C,-17,0', False),
	('L,0,-0.28125', True), ('L,1,-0.28125', False),   # -18, -17
	('S,0,-2.25', True), ('S,0,-2.125', False),         # -18, -17
	('G,0,-5', True), ('G,0,-4', False),                # -20.8, -16.6
])
def test_score_min_end_to_end(bt2, genome_index, score18_read, tmp_path, func, aligns):
	"""A read scoring -18 is reported iff the --score-min function of its length is <= -18."""
	rec = mapped(bt2, genome_index, tmp_path, score18_read.seq, ['--score-min', func])
	assert rec.is_unmapped != aligns
	if aligns:
		assert (rec.pos, rec.tags['AS']) == (score18_read.pos, -18)
	else:
		assert 'YF' not in rec.tags


@pytest.mark.parametrize('func,aligns', [
	(None, True),                          # default G,20,8 -> 56.8
	('G,1,34', True), ('G,1,35', False),   # 157.6, 162.2
	('L,0,1.5', True), ('L,0,1.7', False),
	('C,160,0', True), ('C,161,0', False),
])
def test_score_min_local(bt2, genome_index, clip_reads, tmp_path, func, aligns):
	"""In --local mode the soft-clipped read scoring 160 is reported iff --score-min is <= 160."""
	seq, core = clip_reads['clip']
	args = ['--local'] + (['--score-min', func] if func else [])
	rec = mapped(bt2, genome_index, tmp_path, seq, args)
	assert rec.is_unmapped != aligns
	if aligns:
		assert (rec.pos, rec.tags['AS']) == (core.pos, 160)


@pytest.mark.parametrize('args,filtered', [
	(['--local', '--ma', '0'], True),                  # best possible 0 < 56.8
	(['--local', '--score-min', 'C,201,0'], True),     # best possible 200
	(['--local', '--score-min', 'C,200,0'], False),
])
def test_yf_sc_filter(bt2, genome_index, genome_dict, tmp_path, args, filtered):
	"""A read that cannot reach --score-min even with a perfect match is filtered with YF:Z:SC."""
	rec = mapped(bt2, genome_index, tmp_path, genome_dict['chr1'][100:200], args)
	if filtered:
		assert rec.is_unmapped
		assert rec.tags['YF'] == 'SC'
	else:
		assert (rec.pos, rec.cigar, rec.tags['AS']) == (101, '100M', 200)


def test_yf_ln_filter(bt2, genome_index, genome_dict, tmp_path):
	"""A read no longer than -N is filtered with YF:Z:LN; a longer one is not."""
	fq = write_reads(tmp_path, [('one', 'A'), ('two', genome_dict['chr1'][300:302])])
	recs = one_per_read(align(bt2, genome_index, fq, ['-N', 1]))
	assert recs['one'].is_unmapped
	assert recs['one'].tags['YF'] == 'LN'
	assert 'YF' not in recs['two'].tags


# ---------------------------------------------------------------------------
# Strands
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('args,fw_aligns,rc_aligns', [
	([], True, True),
	(['--nofw'], False, True),
	(['--norc'], True, False),
	(['--nofw', '--norc'], False, False),
])
def test_nofw_norc_unpaired(bt2, genome_index, genome_dict, tmp_path, args, fw_aligns, rc_aligns):
	"""--nofw skips the forward strand and --norc the reverse strand for unpaired reads."""
	seg = genome_dict['chr3'][1000:1060]
	fq = write_reads(tmp_path, [('fw', seg), ('rc', h.revcomp(seg))])
	recs = one_per_read(align(bt2, genome_index, fq, args))
	assert recs['fw'].is_unmapped != fw_aligns
	assert recs['rc'].is_unmapped != rc_aligns
	if fw_aligns:
		assert (recs['fw'].flag, recs['fw'].pos) == (0, 1001)
	if rc_aligns:
		assert (recs['rc'].flag, recs['rc'].pos) == (h.REVERSE, 1001)


@pytest.mark.parametrize('opt,watson,crick', [
	(None, True, True), ('--nofw', False, True), ('--norc', True, False),
])
def test_nofw_norc_paired_fragments(bt2, genome_index, genome_dict, tmp_path, opt, watson, crick):
	"""For pairs, --nofw/--norc apply to the fragment's strand: --nofw keeps only Crick fragments."""
	ref = genome_dict['chr1']
	frags = {'watson': ref[2000:2300], 'crick': h.revcomp(ref[3000:3300])}
	m1 = [(n, f[:50]) for n, f in sorted(frags.items())]
	m2 = [(n, h.revcomp(f[-50:])) for n, f in sorted(frags.items())]
	f1 = write_reads(tmp_path, m1, '1.fq')
	f2 = write_reads(tmp_path, m2, '2.fq')
	args = ['-p', '1', '-x', genome_index, '-1', f1, '-2', f2] + ([opt] if opt else [])
	_, recs, _ = bt2.align(args)
	by = h.by_name(recs)
	for name, expect in (('watson', watson), ('crick', crick)):
		assert len(by[name]) == 2
		for rec in by[name]:
			assert rec.tags['YT'] == ('CP' if expect else 'UP')
			assert rec.is_unmapped != expect
	if watson:
		assert sorted(r.pos for r in by['watson']) == [2001, 2251]
	if crick:
		assert sorted(r.pos for r in by['crick']) == [3001, 3251]


# ---------------------------------------------------------------------------
# Repeats: -k, -a, XS, MAPQ, -l, -d, --seed, --non-deterministic
# ---------------------------------------------------------------------------

REP_S_EXACT = {('rep1', 301, False), ('rep1', 1301, False), ('rep2', 501, True)}
REP_S_NEAR = ('rep2', 1201, False)
REP_T = {('rep1', 2001, False), ('rep1', 2601, False), ('rep2', 1701, False)}


@pytest.fixture(scope='module')
def rep():
	"""Two sequences with repeats of 100 bp segments S and T:
	S exactly at rep1:300 and rep1:1300 (forward) and rep2:500 (reverse
	complement), and with substitutions at 30 and 70 at rep2:1200; T exactly
	at rep1:2000, rep1:2600 and rep2:1700, all forward. U = rep1:700-800 is
	unique."""
	rng = random.Random(7)
	S = h.random_dna(100, rng)
	T = h.random_dna(100, rng)
	r1 = list(h.random_dna(3000, rng))
	r2 = list(h.random_dna(2000, rng))

	def put(chrom, off, seg):
		chrom[off:off + len(seg)] = list(seg)

	put(r1, 300, S)
	put(r1, 1300, S)
	put(r1, 2000, T)
	put(r1, 2600, T)
	put(r2, 500, h.revcomp(S))
	put(r2, 1200, h.mutate(S, [30, 70]))
	put(r2, 1700, T)
	seqs = [('rep1', ''.join(r1)), ('rep2', ''.join(r2))]
	return {'seqs': seqs, 'S': S, 'T': T, 'U': seqs[0][1][700:800],
		'Snear': h.mutate(S, [30, 70])}


@pytest.fixture(scope='module')
def rep_index(bt2, rep, tmp_path_factory):
	d = tmp_path_factory.mktemp('rep_index')
	fa = h.write_fasta(str(d / 'rep.fa'), rep['seqs'])
	prefix = str(d / 'rep')
	bt2.build_index(fa, prefix, extra=['--seed', '0'])
	return prefix


@pytest.fixture(scope='module')
def rep_fastq(rep, tmp_path_factory):
	reads = [('S', rep['S']), ('T', rep['T']), ('U', rep['U']), ('Snear', rep['Snear'])]
	path = tmp_path_factory.mktemp('rep_reads') / 'reads.fq'
	return h.write_fastq(str(path), reads)


def loc(rec):
	return (rec.rname, rec.pos, rec.is_reverse)


def test_repeat_default_mode(bt2, rep_index, rep_fastq):
	"""Default mode reports one best alignment; XS is the second-best score and absent for a unique read."""
	recs = one_per_read(align(bt2, rep_index, rep_fastq))
	s = recs['S']
	assert loc(s) in REP_S_EXACT
	assert (s.tags['AS'], s.tags['XS']) == (0, 0)
	assert loc(recs['T']) in REP_T
	assert (recs['T'].tags['AS'], recs['T'].tags['XS']) == (0, 0)
	near = recs['Snear']
	assert loc(near) == REP_S_NEAR
	assert (near.tags['AS'], near.tags['XS']) == (0, -12)
	assert near.tags['MD'] == '100'
	u = recs['U']
	assert (u.rname, u.pos, u.tags['AS']) == ('rep1', 701, 0)
	assert 'XS' not in u.tags
	for rec in recs.values():
		assert not rec.flag & h.SECONDARY


def test_repeat_mapq(bt2, rep_index, rep_fastq):
	"""MAPQ is 42 for a unique perfect hit, 1 for a tie and 30 when the second best is 12 lower (pins the current values).

	The MANUAL gives no MAPQ formula; the values follow BowtieMapq2 in
	unique.h for a 100 bp read, whose minimum score -60.6 truncates to -60
	(diff = 60): no second best and AS - min >= 0.8 * diff gives 42; AS ==
	XS with AS - min >= 0.67 * diff gives 1; AS 0 and XS -12 fall just below
	the 0.2f * diff band (12 < 12.0000002) into the 0.1 band, where AS ==
	perfect gives 30."""
	recs = one_per_read(align(bt2, rep_index, rep_fastq))
	assert recs['U'].mapq == 42
	assert recs['S'].mapq == 1
	assert recs['T'].mapq == 1
	assert (recs['Snear'].tags['AS'], recs['Snear'].tags['XS']) == (0, -12)
	assert recs['Snear'].mapq == 30


def check_k_records(recs, expected_locs, n):
	"""n distinct alignments from expected_locs, AS descending, secondary flag
	and MAPQ 255 on all but the first."""
	assert len(recs) == n
	locs = [loc(r) for r in recs]
	assert len(set(locs)) == n
	assert set(locs) <= set(expected_locs)
	scores = [r.tags['AS'] for r in recs]
	assert scores == sorted(scores, reverse=True)
	assert not recs[0].flag & h.SECONDARY
	for r in recs[1:]:
		assert r.flag & h.SECONDARY
		assert r.mapq == 255


@pytest.mark.parametrize('k', [1, 2, 3, 4, 10])
def test_k(bt2, rep_index, rep_fastq, rep, k):
	"""-k N reports min(N, number of valid alignments) distinct alignments in descending AS order."""
	by = h.by_name(align(bt2, rep_index, rep_fastq, ['-k', k]))
	s_all = REP_S_EXACT | {REP_S_NEAR}
	check_k_records(by['S'], s_all, min(k, 4))
	if k >= 4:
		near = by['S'][-1]
		assert loc(near) == REP_S_NEAR
		# MD gives the reference (near copy) bases at the two substitutions.
		assert (near.tags['AS'], near.tags['XM'], near.tags['MD']) == (-12, 2, '30%s39%s29' % (
			rep['Snear'][30], rep['Snear'][70]))
	check_k_records(by['T'], REP_T, min(k, 3))
	check_k_records(by['U'], {('rep1', 701, False)}, 1)
	for r in by['S']:
		# SEQ is reverse complemented for reverse-strand alignments, secondary or not.
		assert r.seq == (h.revcomp(rep['S']) if r.is_reverse else rep['S'])


@pytest.mark.parametrize('opt', ['-a', '--all'])
def test_all(bt2, rep_index, rep_fastq, opt):
	"""-a/--all reports every valid alignment: S at 4 places, T at 3, U at 1."""
	by = h.by_name(align(bt2, rep_index, rep_fastq, [opt]))
	check_k_records(by['S'], REP_S_EXACT | {REP_S_NEAR}, 4)
	assert by['S'][-1].tags['AS'] == -12
	check_k_records(by['T'], REP_T, 3)
	check_k_records(by['U'], {('rep1', 701, False)}, 1)
	for r in by['S'] + by['T']:
		assert 'XS' in r.tags


def test_all_alias(bt2, rep_index, rep_fastq):
	"""-a and --all give identical records."""
	same_records(align(bt2, rep_index, rep_fastq, ['-a']),
		align(bt2, rep_index, rep_fastq, ['--all']), 4)


NO_UPFRONT = ['--no-exact-upfront', '--no-1mm-upfront']


@pytest.mark.parametrize('lowseeds,t_aligns', [
	('0', True),            # 0 disables the cut
	('2', False), ('3', True),
	# Relative forms: the cut is ceil(R * value / divisor), where R is the
	# number of reference sequences in the index (2 here) and the divisor is
	# 100 for %, 1000 for m ("mili" in the usage text).
	('100%', False), ('101%', True), ('150%', True),
	('1000m', False), ('1001m', True),
	# n is "nano" in the usage text but the current version divides by
	# 10^6; these two cases come out the same with
	# 10^6 and 10^9: cuts of 2 and 1, then 2001 and 3.
	('1000000n', False), ('1000000001n', True),
])
def test_lowseeds(bt2, rep_index, rep_fastq, lowseeds, t_aligns):
	"""-l drops seeds whose hit range exceeds the cut: T's seeds hit 3 places, U's 1.

	-l is not in the MANUAL. The usage text says the relative forms are
	"relative to idx size"; that the base is the number of reference
	sequences (not the index length) pins the current output (bt2_search.cpp,
	lowseeds_ncut)."""
	recs = one_per_read(align(bt2, rep_index, rep_fastq, ['-l', lowseeds] + NO_UPFRONT))
	assert recs['T'].is_unmapped != t_aligns
	if t_aligns:
		assert loc(recs['T']) in REP_T
	assert (recs['U'].rname, recs['U'].pos) == ('rep1', 701)


def test_lowseeds_does_not_affect_upfront_search(bt2, rep_index, rep_fastq):
	"""With the exact upfront search on, -l 1 still finds T, whose end-to-end exact hit needs no seed."""
	recs = one_per_read(align(bt2, rep_index, rep_fastq, ['-l', '1']))
	assert loc(recs['T']) in REP_T


@pytest.mark.parametrize('opt', ['-d', '--deterministic-seeds'])
def test_deterministic_seeds(bt2, rep_index, rep_fastq, opt):
	"""-d/--deterministic-seeds with -a and both upfront searches off finds every alignment."""
	by = h.by_name(align(bt2, rep_index, rep_fastq, [opt, '-a'] + NO_UPFRONT))
	check_k_records(by['S'], REP_S_EXACT | {REP_S_NEAR}, 4)
	check_k_records(by['T'], REP_T, 3)
	check_k_records(by['U'], {('rep1', 701, False)}, 1)


def test_deterministic_seeds_alias(bt2, rep_index, rep_fastq):
	"""-d and --deterministic-seeds give identical records."""
	args = ['-a'] + NO_UPFRONT
	same_records(align(bt2, rep_index, rep_fastq, ['-d'] + args),
		align(bt2, rep_index, rep_fastq, ['--deterministic-seeds'] + args), 4)


@pytest.mark.parametrize('args', [['-a'], ['-a', '--no-1mm-upfront'], ['-a', '--no-exact-upfront']])
def test_deterministic_seeds_requires_no_upfront(bt2, rep_index, rep_fastq, args):
	"""-d without both --no-exact-upfront and --no-1mm-upfront is rejected, as the MANUAL requires."""
	proc = bt2.run('bowtie2', ['-x', rep_index, '-U', rep_fastq, '-d'] + args, check=False)
	assert proc.returncode != 0
	assert '-d must be used with --no-exact-upfront and --no-1mm-upfront' in proc.stderr


@pytest.fixture(scope='module')
def tied_reads(rep, tmp_path_factory):
	"""30 differently named copies of S, which has three equally good hits."""
	path = tmp_path_factory.mktemp('tied') / 'reads.fq'
	return h.write_fastq(str(path), [('s%02d' % i, rep['S']) for i in range(30)])


def test_seed(bt2, rep_index, tied_reads):
	"""--seed changes how ties are broken, the default is --seed 0, and a fixed seed reproduces the output."""
	runs = {}
	for seed in (None, 0, 1, 2):
		args = ['--seed', seed] if seed is not None else []
		runs[seed] = align(bt2, rep_index, tied_reads, args)
	assert body(runs[None]) == body(runs[0])
	assert body(align(bt2, rep_index, tied_reads, ['--seed', 1])) == body(runs[1])
	picks = {}
	for seed, recs in runs.items():
		assert len(recs) == 30
		for r in recs:
			assert loc(r) in REP_S_EXACT
			assert r.tags['AS'] == 0
		picks[seed] = [loc(r) for r in recs]
		assert len(set(picks[seed])) > 1      # ties go different ways across reads
	assert picks[0] != picks[1]
	assert picks[1] != picks[2]


def test_identical_reads_identical_alignments(bt2, rep_index, rep, tmp_path):
	"""Reads with the same name, sequence and quality get the same alignment among ties."""
	fq = write_reads(tmp_path, [('dup', rep['S'])] * 30)
	recs = align(bt2, rep_index, fq)
	assert len(recs) == 30
	assert len(set(loc(r) for r in recs)) == 1
	assert loc(recs[0]) in REP_S_EXACT


def test_non_deterministic(bt2, rep_index, rep, tmp_path):
	"""--non-deterministic breaks ties differently for identical reads."""
	fq = write_reads(tmp_path, [('dup', rep['S'])] * 30)
	recs = align(bt2, rep_index, fq, ['--non-deterministic'])
	assert len(recs) == 30
	locs = [loc(r) for r in recs]
	assert set(locs) <= REP_S_EXACT
	assert len(set(locs)) > 1
	assert all(r.tags['AS'] == 0 for r in recs)


# ---------------------------------------------------------------------------
# Effort: -D, -R
# ---------------------------------------------------------------------------

def total_score(records):
	assert all(not r.is_unmapped for r in records)
	return sum(r.tags['AS'] for r in records)


def test_D_more_extensions_find_better_alignments(bt2, family_index, family_fastq):
	"""On repetitive reads, raising -D finds better alignments (higher total AS); the size of the gain pins current output."""
	scores = [total_score(align(bt2, family_index, family_fastq, ['-D', d, '-R', 1]))
		for d in (1, 5, 15)]
	assert scores[0] < scores[1] < scores[2]


def test_D_irrelevant_for_unique_reads(bt2, lambda_index, example_reads):
	"""On the non-repetitive lambda reads, -D 1 and -D 40 give the same records."""
	fq = example_reads['reads_1']
	same_records(align(bt2, lambda_index, fq, ['-D', 1, '-u', 2000]),
		align(bt2, lambda_index, fq, ['-D', 40, '-u', 2000]), 2000)


def test_R_reseeding_helps_repetitive_reads(bt2, family_index, family_fastq):
	"""Reads whose seeds hit hundreds of places are re-seeded: -R 6 finds better alignments than -R 1 (gain pins current output)."""
	r1 = align(bt2, family_index, family_fastq, ['-R', 1])
	r6 = align(bt2, family_index, family_fastq, ['-R', 6])
	assert total_score(r1) < total_score(r6)


def test_R_irrelevant_for_unique_reads(bt2, lambda_index, example_reads):
	"""Reads without repetitive seeds are never re-seeded: -R 1 and -R 6 give the same records."""
	fq = example_reads['reads_1']
	same_records(align(bt2, lambda_index, fq, ['-R', 1, '-u', 2000]),
		align(bt2, lambda_index, fq, ['-R', 6, '-u', 2000]), 2000)


# ---------------------------------------------------------------------------
# Optional-field layout
# ---------------------------------------------------------------------------

def test_optional_field_order(bt2, rep_index, genome_index, rep_fastq, tmp_path):
	"""Default optional fields and their order for unpaired reads (pins current output)."""
	recs = one_per_read(align(bt2, rep_index, rep_fastq))
	assert recs['U'].tag_order == ['AS', 'XN', 'XM', 'XO', 'XG', 'NM', 'MD', 'YT']
	assert recs['S'].tag_order == ['AS', 'XS', 'XN', 'XM', 'XO', 'XG', 'NM', 'MD', 'YT']
	fq = write_reads(tmp_path, [('random', 'ACGT' * 25), ('short', 'A')])
	recs = one_per_read(align(bt2, genome_index, fq, ['-N', 1]))
	assert recs['random'].is_unmapped
	assert recs['random'].tag_order == ['YT']
	assert recs['short'].tag_order == ['YT', 'YF']
