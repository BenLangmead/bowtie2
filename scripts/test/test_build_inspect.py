"""Tests for bowtie2-build and bowtie2-inspect.

Most expected values are derived without bowtie2: the references are the
fixed random `genome` (with N runs, IUPAC codes, a gap and lowercase pasted
over known stretches), reads are cut from them at known places, and the
MANUAL's descriptions give what bowtie2-inspect prints and what the
--offrate and --ftabchars options change.

The central invariant: the options that choose HOW the suffix array is
built (--bmax, --bmaxdivn, --dcv, --packed, --noauto, --threads, --seed)
must not change the index at all, and the options that change the index's
shape (--offrate, --ftabchars, --large-index) must not change alignments.
"""

import gzip
import os
import re

import pytest

import bt2_helpers as h


SUFFIXES = ('1', '2', '3', '4', 'rev.1', 'rev.2')
READ_LEN = 60

# Every IUPAC ambiguity code plus X, in both cases. The MANUAL says
# bowtie2-inspect prints the reference "with all non-A/C/G/T characters
# converted to Ns".
IUPAC = 'RYKMSWBDHVNX'


def index_names(ext='bt2', suffixes=SUFFIXES):
	"""File-name endings of an index, e.g. '.rev.1.bt2'."""
	return sorted('.%s.%s' % (s, ext) for s in suffixes)


def files_written(prefix):
	"""Endings of every file in prefix's directory that starts with
	'<basename>.', e.g. ['.1.bt2', '.2.bt2', ...]."""
	d, base = os.path.split(prefix)
	return sorted(f[len(base):] for f in os.listdir(d) if f.startswith(base + '.'))


def read_bytes(path):
	with open(path, 'rb') as f:
		return f.read()


def assert_same_files(prefix_a, prefix_b, ext='bt2', suffixes=SUFFIXES):
	"""Each index file of prefix_a is byte-identical to prefix_b's."""
	for s in suffixes:
		a = '%s.%s.%s' % (prefix_a, s, ext)
		b = '%s.%s.%s' % (prefix_b, s, ext)
		assert read_bytes(a) == read_bytes(b), '%s differs from %s' % (a, b)


def inspected(seq):
	"""What bowtie2-inspect prints for a reference sequence, per the
	MANUAL: upper case, and every non-A/C/G/T character becomes N."""
	return ''.join(c if c in 'ACGT' else 'N' for c in seq.upper())


def fasta_text(records, width):
	"""FASTA text with `width` bases per line and no blank lines."""
	out = []
	for name, seq in records:
		out.append('>%s\n' % name)
		for i in range(0, len(seq), width):
			out.append(seq[i:i + width] + '\n')
	return ''.join(out)


def parse_fasta(text):
	"""Parse FASTA text into (name line, sequence) pairs."""
	records = []
	for line in text.splitlines():
		if line.startswith('>'):
			records.append([line[1:], []])
		elif line.strip():
			records[-1][1].append(line.strip())
	return [(name, ''.join(parts)) for name, parts in records]


def assert_same_text(got, expected):
	"""got == expected. A mismatch is reported as the first differing line:
	pytest's own diff of two multi-line strings that differ on every line
	takes minutes for the ~10,000-line outputs compared here."""
	if got == expected:
		return
	g, e = got.splitlines(), expected.splitlines()
	i = next((i for i, (a, b) in enumerate(zip(g, e)) if a != b), min(len(g), len(e)))
	pytest.fail('texts differ first at line %d: got %r, expected %r (%d vs %d lines)' % (
		i + 1, g[i] if i < len(g) else None, e[i] if i < len(e) else None, len(g), len(e)))


def summary_fields(text):
	"""bowtie2-inspect -s output as a list of tab-separated field lists."""
	return [line.split('\t') for line in text.splitlines() if line]


def summary_value(text, key):
	values = [f[1:] for f in summary_fields(text) if f[0] == key]
	assert len(values) == 1, (key, text)
	return values[0]


def sequence_lines(text):
	"""The 'Sequence-<i>' rows of bowtie2-inspect -s, as (i, name, len)."""
	rows = []
	for f in summary_fields(text):
		if f[0].startswith('Sequence-'):
			assert len(f) == 3, f
			rows.append((int(f[0][len('Sequence-'):]), f[1], int(f[2])))
	return rows


def logged_binary(stderr):
	"""Basename of the binary a wrapper run with --verbose says it ran."""
	m = re.search(r'Command: (\S+)', stderr)
	assert m, stderr
	return os.path.basename(m.group(1))


def md_tag(ref, mismatches):
	"""SAM MD:Z for an ungapped alignment against `ref` with substitutions
	at the given 0-based positions."""
	out, run = [], 0
	for i, base in enumerate(ref):
		if i in mismatches:
			out.append('%d%s' % (run, base))
			run = 0
		else:
			run += 1
	out.append(str(run))
	return ''.join(out)


def version():
	"""The version the binaries were built with: BOWTIE2_VERSION next to
	them, else in the repository."""
	for d in (h.BIN_DIR, h.ROOT):
		path = os.path.join(d, 'BOWTIE2_VERSION')
		if os.path.exists(path):
			with open(path) as f:
				return f.read().strip()
	pytest.skip('BOWTIE2_VERSION not found')


def align_lines(bt2, index, args):
	"""Header (minus @PG, which holds the command line) and alignment
	lines of a single-threaded bowtie2 run."""
	header, records, _ = bt2.align(['-p', '1', '-x', index] + list(args))
	return ([l for l in header if not l.startswith('@PG')],
		[r.line for r in records])


def arg_id(value):
	"""Test id for a parametrized list of command-line arguments."""
	if isinstance(value, list):
		return ' '.join(value) or 'default'
	return None


def splice(seq, start, text):
	return seq[:start] + text + seq[start + len(text):]


# ---------------------------------------------------------------------------
# A reference with N runs, IUPAC codes, a gap and lowercase pasted over known
# stretches of `genome`. Pasting keeps every coordinate, so reads cut from the
# untouched stretches have the same positions as in `genome`.

@pytest.fixture(scope='module')
def mixed_ref(genome_dict):
	chr1 = splice(genome_dict['chr1'], 4000, 'RYKMSWBDHV')   # IUPAC run
	chr1 = splice(chr1, 4990, 'N' * 10)                       # trailing N run
	chr2 = splice(genome_dict['chr2'], 1000, 'N' * 40)        # internal N run
	chr2 = splice(chr2, 2500, '--')                           # gap characters
	chr3 = genome_dict['chr3']
	chr3 = 'N' * 10 + chr3[10:500] + chr3[500:800].lower() + chr3[800:]
	return [('chr1', chr1), ('chr2', chr2), ('chr3', chr3)]


@pytest.fixture(scope='module')
def mixed_n(mixed_ref):
	"""Number of unambiguous (A/C/G/T) bases in mixed_ref."""
	return sum(1 for _, s in mixed_ref for c in s.upper() if c in 'ACGT')


@pytest.fixture(scope='module')
def mixed_fasta(mixed_ref, tmp_path_factory):
	path = tmp_path_factory.mktemp('mixed') / 'mixed.fa'
	return h.write_fasta(str(path), mixed_ref)


@pytest.fixture(scope='module')
def mixed_index(bt2, mixed_fasta, tmp_path_factory):
	"""mixed_ref built with default settings."""
	prefix = str(tmp_path_factory.mktemp('mixed_index') / 'mixed')
	bt2.build_index(mixed_fasta, prefix)
	return prefix


@pytest.fixture(scope='module')
def mixed_index_large(bt2, mixed_fasta, tmp_path_factory):
	"""mixed_ref built with --large-index."""
	prefix = str(tmp_path_factory.mktemp('mixed_index_large') / 'mixed')
	bt2.build_index(mixed_fasta, prefix, extra=['--large-index'])
	return prefix


# (name, chrom, 0-based start, strand, mismatch offsets along the reference)
READ_SPECS = [
	('fw_chr1_100', 'chr1', 100, '+', ()),
	('rc_chr1_2500', 'chr1', 2500, '-', (30,)),
	('fw_chr1_4100_after_iupac', 'chr1', 4100, '+', (5, 50)),
	('fw_chr2_200', 'chr2', 200, '+', ()),
	('rc_chr2_1100_after_nrun', 'chr2', 1100, '-', ()),
	('fw_chr2_2900', 'chr2', 2900, '+', (20,)),
	('fw_chr3_20_after_leading_n', 'chr3', 20, '+', ()),
	('fw_chr3_550_lowercase', 'chr3', 550, '+', ()),
	('rc_chr3_1900', 'chr3', 1900, '-', (10,)),
]


@pytest.fixture(scope='module')
def mixed_reads(genome_dict, tmp_path_factory):
	"""FASTQ of READ_SPECS and the expected alignment of each read."""
	reads, expect = [], {}
	for name, chrom, start, strand, mm in READ_SPECS:
		ref = genome_dict[chrom][start:start + READ_LEN]
		frag = h.mutate(ref, mm)
		reads.append((name, frag if strand == '+' else h.revcomp(frag)))
		expect[name] = dict(rname=chrom, pos=start + 1,
			flag=0 if strand == '+' else h.REVERSE, nm=len(mm),
			md=md_tag(ref, set(mm)))
	path = tmp_path_factory.mktemp('mixed_reads') / 'reads.fq'
	return h.write_fastq(str(path), reads), expect


@pytest.fixture(scope='module')
def mixed_pair(genome_dict, tmp_path_factory):
	"""One read pair from chr1: mate 1 forward at 1000, mate 2 the reverse
	complement of 1240..1300, so the fragment is 300 bp."""
	chr1 = genome_dict['chr1']
	d = tmp_path_factory.mktemp('mixed_pair')
	m1 = h.write_fastq(str(d / 'm1.fq'), [('pair', chr1[1000:1050])])
	m2 = h.write_fastq(str(d / 'm2.fq'), [('pair', h.revcomp(chr1[1240:1300]))])
	return m1, m2


@pytest.fixture(scope='module')
def mixed_baseline(bt2, mixed_index, mixed_reads, mixed_pair):
	"""Inspect output and alignments of the default-built mixed index."""
	fq, _ = mixed_reads
	m1, m2 = mixed_pair
	return {
		'fasta': bt2.run('bowtie2-inspect', [mixed_index]).stdout,
		'names': bt2.run('bowtie2-inspect', ['-n', mixed_index]).stdout,
		'summary': bt2.run('bowtie2-inspect', ['-s', mixed_index]).stdout,
		'unpaired': align_lines(bt2, mixed_index, ['-U', fq]),
		'paired': align_lines(bt2, mixed_index, ['-1', m1, '-2', m2]),
	}


def assert_matches_baseline(bt2, prefix, baseline, reads, pair):
	"""prefix inspects and aligns exactly like the default mixed index."""
	assert_same_text(bt2.run('bowtie2-inspect', [prefix]).stdout, baseline['fasta'])
	assert bt2.run('bowtie2-inspect', ['-n', prefix]).stdout == baseline['names']
	assert bt2.run('bowtie2-inspect', ['-s', prefix]).stdout == baseline['summary']
	assert align_lines(bt2, prefix, ['-U', reads[0]]) == baseline['unpaired']
	assert align_lines(bt2, prefix, ['-1', pair[0], '-2', pair[1]]) == baseline['paired']


# ---------------------------------------------------------------------------
# Files written

def test_build_writes_six_small_index_files(mixed_index):
	"""A default build writes exactly .1-.4.bt2, .rev.1.bt2 and .rev.2.bt2."""
	assert files_written(mixed_index) == index_names('bt2')


def test_build_large_index_writes_bt2l_files(mixed_index_large):
	"""--large-index writes the same six files with the .bt2l extension."""
	assert files_written(mixed_index_large) == index_names('bt2l')


def test_large_index_inspects_and_aligns_like_small(bt2, mixed_index_large,
		mixed_baseline, mixed_reads, mixed_pair):
	"""A .bt2l index of a reference inspects and aligns exactly like its .bt2."""
	assert_matches_baseline(bt2, mixed_index_large, mixed_baseline, mixed_reads,
		mixed_pair)


@pytest.mark.parametrize('flag', ['-r', '--noref'])
def test_build_noref_skips_3_and_4(bt2, mixed_fasta, mixed_index, mixed_baseline,
		tmp_path, flag):
	"""-r/--noref writes only .1/.2/.rev.1/.rev.2, identical to a full build's, still readable by -n and -s."""
	prefix = str(tmp_path / 'noref')
	bt2.build_index(mixed_fasta, prefix, extra=[flag])
	kept = ('1', '2', 'rev.1', 'rev.2')
	assert files_written(prefix) == index_names('bt2', kept)
	assert_same_files(prefix, mixed_index, suffixes=kept)
	assert bt2.run('bowtie2-inspect', ['-n', prefix]).stdout == mixed_baseline['names']
	assert bt2.run('bowtie2-inspect', ['-s', prefix]).stdout == mixed_baseline['summary']


@pytest.mark.parametrize('flag', ['-3', '--justref'])
def test_build_justref_writes_only_3_and_4(bt2, mixed_fasta, mixed_index,
		tmp_path, flag):
	"""-3/--justref writes only .3 and .4, identical to a full build's."""
	prefix = str(tmp_path / 'justref')
	bt2.build_index(mixed_fasta, prefix, extra=[flag])
	assert files_written(prefix) == index_names('bt2', ('3', '4'))
	assert_same_files(prefix, mixed_index, suffixes=('3', '4'))


def test_noref_plus_justref_is_the_full_index(bt2, mixed_fasta, mixed_index,
		mixed_baseline, mixed_reads, mixed_pair, tmp_path):
	"""--noref then --justref into one basename reassembles the default index."""
	prefix = str(tmp_path / 'split')
	bt2.build_index(mixed_fasta, prefix, extra=['--noref'])
	bt2.build_index(mixed_fasta, prefix, extra=['--justref'])
	assert files_written(prefix) == index_names('bt2')
	assert_same_files(prefix, mixed_index)
	assert_matches_baseline(bt2, prefix, mixed_baseline, mixed_reads, mixed_pair)


# ---------------------------------------------------------------------------
# Alignment against the mixed reference: expected values from READ_SPECS

def test_mixed_index_aligns_reads_to_known_positions(bt2, mixed_index, mixed_reads):
	"""Reads beside N runs, IUPAC codes and lowercase align at their source with the planted mismatches."""
	# AS is -6 per mismatch: MX = 6 at Phred 40 ('I').
	fq, expect = mixed_reads
	_, records, _ = bt2.align(['-p', '1', '-x', mixed_index, '-U', fq])
	assert sorted(r.qname for r in records) == sorted(expect)
	for r in records:
		e = expect[r.qname]
		assert (r.rname, r.pos, r.flag, r.cigar) == (e['rname'], e['pos'], e['flag'],
			'%dM' % READ_LEN), r
		assert r.tags['NM'] == e['nm'] and r.tags['XM'] == e['nm'], r
		assert r.tags['MD'] == e['md'], r
		assert r.tags['AS'] == -6 * e['nm'], r


def test_mixed_index_aligns_pair(bt2, mixed_index, mixed_pair):
	"""The planted pair aligns concordantly at 1001/1241 with TLEN +-300."""
	_, records, _ = bt2.align(['-p', '1', '-x', mixed_index, '-1', mixed_pair[0],
		'-2', mixed_pair[1]])
	assert len(records) == 2
	m1, m2 = sorted(records, key=lambda r: r.flag & h.SECOND_IN_PAIR)
	assert (m1.rname, m1.pos, m1.tlen) == ('chr1', 1001, 300)
	assert (m2.rname, m2.pos, m2.tlen) == ('chr1', 1241, -300)
	assert m1.flag == h.PAIRED | h.PROPER_PAIR | h.MATE_REVERSE | h.FIRST_IN_PAIR
	assert m2.flag == h.PAIRED | h.PROPER_PAIR | h.REVERSE | h.SECOND_IN_PAIR


# ---------------------------------------------------------------------------
# Index content invariance: SA construction settings change nothing

STRATEGY_VARIANTS = [
	['-f'],
	['-a'],
	['--noauto'],
	['-p'],
	['--packed'],
	['-a', '-p'],
	['-a', '--bmax', '64'],
	['--bmax', '500'],
	['-a', '--bmaxdivn', '16'],
	['--bmaxdivn', '1'],
	['-a', '--dcv', '16'],
	['-a', '--dcv', '4096'],
	['--dcv', '64'],
	['--threads', '4'],
	['--threads', '3', '-a', '--bmaxdivn', '9', '--dcv', '32'],
	['--seed', '7'],
	['--seed', '123456'],
	['-a', '--packed', '--bmax', '100', '--dcv', '128', '--threads', '2',
		'--seed', '99'],
]


@pytest.mark.parametrize('extra', STRATEGY_VARIANTS, ids=arg_id)
def test_build_strategy_does_not_change_index(bt2, mixed_fasta, mixed_index,
		mixed_baseline, mixed_reads, mixed_pair, tmp_path, extra):
	"""SA-construction settings, threads and seed give a byte-identical index that inspects and aligns the same."""
	prefix = str(tmp_path / 'variant')
	bt2.build_index(mixed_fasta, prefix, extra=extra)
	assert files_written(prefix) == index_names('bt2')
	assert_same_files(prefix, mixed_index)
	assert_matches_baseline(bt2, prefix, mixed_baseline, mixed_reads, mixed_pair)


@pytest.mark.parametrize('extra', [
	['--large-index', '-a', '-p'],
	['--large-index', '--bmax', '200', '--dcv', '32'],
	['--large-index', '--threads', '3', '--seed', '5'],
], ids=arg_id)
def test_large_build_strategy_does_not_change_index(bt2, mixed_fasta,
		mixed_index_large, tmp_path, extra):
	"""The same invariance for large indexes: byte-identical .bt2l files."""
	prefix = str(tmp_path / 'variant')
	bt2.build_index(mixed_fasta, prefix, extra=extra)
	assert files_written(prefix) == index_names('bt2l')
	assert_same_files(prefix, mixed_index_large, ext='bt2l')


@pytest.mark.parametrize('extra,message', [
	([], 'bmax according to bmaxDivN setting: {n4}'),
	(['--threads', '4'], 'bmax according to bmaxDivN setting: {n16}'),
	(['--threads', '4', '--bmaxdivn', '5'], 'bmax according to bmaxDivN setting: {n5}'),
	(['--bmaxdivn', '8', '--bmax', '300'], 'bmax according to bmax setting: 300'),
	(['--bmax', '300', '--bmaxdivn', '8'], 'bmax according to bmaxDivN setting: {n8}'),
	(['-a', '--dcv', '64'], '--dcv 64'),
	(['-a', '-p'], 'Strings: packed'),
	(['-a'], 'Strings: unpacked'),
], ids=arg_id)
def test_build_reports_block_settings(bt2, mixed_fasta, mixed_n, tmp_path, extra,
		message):
	"""The build log shows the block size from the MANUAL's --bmax/--bmaxdivn rules (pins the log wording)."""
	# Default n / (4 * threads); an explicit --bmaxdivn is not multiplied by
	# threads; of --bmax and --bmaxdivn the later one wins. n counts A/C/G/T only.
	n = mixed_n
	message = message.format(n4=n // 4, n16=n // 16, n5=n // 5, n8=n // 8)
	prefix = str(tmp_path / 'idx')
	proc = bt2.run('bowtie2-build', extra + [mixed_fasta, prefix])
	assert message in proc.stdout


# ---------------------------------------------------------------------------
# Index shape: --offrate and --ftabchars

@pytest.mark.parametrize('extra,offrate', [
	(['-o', '0'], 0), (['--offrate', '3'], 3), (['-o', '5'], 5), (['-o', '10'], 10),
], ids=arg_id)
@pytest.mark.parametrize('large', [False, True], ids=['small', 'large'])
def test_offrate_sets_sa_sample_density(bt2, mixed_fasta, mixed_n, tmp_path,
		extra, offrate, large):
	"""-o/--offrate k stores one SA offset per 2^k BWT rows and inspect -s reports 'SA-Sample 1 in 2^k'."""
	# n A/C/G/T bases give n+1 BWT rows; each .2 file is a 4-byte header plus
	# ceil((n+1)/2^k) offsets of 4 bytes (8 in a large index).
	prefix = str(tmp_path / 'idx')
	ext, width = ('bt2l', 8) if large else ('bt2', 4)
	bt2.build_index(mixed_fasta, prefix, extra=extra + (['--large-index'] if large else []))
	rows = mixed_n + 1
	samples = (rows + (1 << offrate) - 1) >> offrate
	for s in ('2', 'rev.2'):
		assert os.path.getsize('%s.%s.%s' % (prefix, s, ext)) == 4 + width * samples
	summary = bt2.run('bowtie2-inspect', ['-s', prefix]).stdout
	assert summary_value(summary, 'SA-Sample') == ['1 in %d' % (1 << offrate)]


@pytest.mark.parametrize('extra,chars', [
	(['-t', '1'], 1), (['--ftabchars', '5'], 5), (['-t', '8'], 8),
], ids=arg_id)
def test_ftabchars_sets_lookup_table_size(bt2, mixed_fasta, mixed_index, tmp_path,
		extra, chars):
	"""-t/--ftabchars k: inspect -s reports k and .1 shrinks by the MANUAL's ftab size, 4^11 - 4^(k+1) bytes."""
	# The default is k = 10; allow a few bytes of other k-dependent bookkeeping.
	prefix = str(tmp_path / 'idx')
	bt2.build_index(mixed_fasta, prefix, extra=extra)
	summary = bt2.run('bowtie2-inspect', ['-s', prefix]).stdout
	assert summary_value(summary, 'FTab-Chars') == [str(chars)]
	for s in ('1', 'rev.1'):
		shrink = (os.path.getsize('%s.%s.bt2' % (mixed_index, s))
			- os.path.getsize('%s.%s.bt2' % (prefix, s)))
		ftab = 4 ** 11 - 4 ** (chars + 1)
		assert ftab <= shrink <= ftab + 128, (s, shrink, ftab)
	# The reference itself (.3/.4) does not depend on the ftab.
	assert_same_files(prefix, mixed_index, suffixes=('3', '4'))


def test_default_ftabchars_is_10(bt2, mixed_index):
	"""Without -t, inspect -s reports the MANUAL's default of 10 ftab chars."""
	summary = bt2.run('bowtie2-inspect', ['-s', mixed_index]).stdout
	assert summary_value(summary, 'FTab-Chars') == ['10']


@pytest.fixture(scope='module')
def lambda_reads(example_reads, tmp_path_factory):
	"""The first 200 pairs of the example paired reads and the first 100
	example long reads."""
	d = tmp_path_factory.mktemp('lambda_reads')
	out = {}
	for key, n in (('reads_1', 800), ('reads_2', 800), ('longreads', 400)):
		with open(example_reads[key]) as f:
			lines = [next(f) for _ in range(n)]
		path = str(d / (key + '.fq'))
		with open(path, 'w') as f:
			f.writelines(lines)
		out[key] = path
	return out


@pytest.fixture(scope='module')
def lambda_shipped_alignments(bt2, lambda_index, lambda_reads):
	return lambda_alignments(bt2, lambda_index, lambda_reads)


def lambda_alignments(bt2, index, reads):
	return (align_lines(bt2, index, ['-1', reads['reads_1'], '-2', reads['reads_2']]),
		align_lines(bt2, index, ['-U', reads['longreads']]))


@pytest.mark.parametrize('extra', [
	[],
	['-o', '0'],
	['--offrate', '9'],
	['-t', '4'],
	['--ftabchars', '7'],
	['-o', '2', '-t', '6'],
	['--large-index'],
	['--large-index', '-o', '1', '-t', '5'],
], ids=arg_id)
def test_index_shape_does_not_change_alignments(bt2, lambda_fasta, lambda_index,
		lambda_reads, lambda_shipped_alignments, tmp_path, extra):
	"""Any --offrate/--ftabchars, small or large, aligns the example reads exactly like the shipped -o 7 index."""
	prefix = str(tmp_path / 'lambda')
	bt2.build_index(lambda_fasta, prefix, extra=extra)
	assert lambda_alignments(bt2, prefix, lambda_reads) == lambda_shipped_alignments
	assert_same_text(bt2.run('bowtie2-inspect', [prefix]).stdout,
		bt2.run('bowtie2-inspect', [lambda_index]).stdout)


def test_rebuilding_lambda_reproduces_shipped_index(bt2, lambda_fasta, lambda_index,
		tmp_path):
	"""bowtie2-build -o 7 on the lambda FASTA reproduces the six shipped example/index files byte for byte."""
	prefix = str(tmp_path / 'lambda_virus')
	bt2.build_index(lambda_fasta, prefix, extra=['-o', '7'])
	assert_same_files(prefix, lambda_index)


# ---------------------------------------------------------------------------
# Reference input: formats, file lists, -c

def write_text(path, text, compress=False):
	data = text.encode()
	if compress:
		with gzip.open(path, 'wb') as f:
			f.write(data)
	else:
		with open(path, 'wb') as f:
			f.write(data)
	return path


@pytest.mark.parametrize('layout', ['width13', 'one_line', 'crlf', 'gzip'])
def test_fasta_layout_does_not_change_index(bt2, mixed_ref, mixed_index, tmp_path,
		layout):
	"""Line width, CRLF line ends and gzip compression of the reference FASTA give a byte-identical index."""
	if layout == 'width13':
		path = write_text(str(tmp_path / 'ref.fa'), fasta_text(mixed_ref, 13))
	elif layout == 'one_line':
		path = write_text(str(tmp_path / 'ref.fa'), fasta_text(mixed_ref, 10 ** 6))
	elif layout == 'crlf':
		path = write_text(str(tmp_path / 'ref.fa'),
			fasta_text(mixed_ref, 60).replace('\n', '\r\n'))
	else:
		path = write_text(str(tmp_path / 'ref.fa.gz'), fasta_text(mixed_ref, 60),
			compress=True)
	prefix = str(tmp_path / 'idx')
	bt2.build_index(path, prefix)
	assert_same_files(prefix, mixed_index)


def split_reference(mixed_ref, tmp_path, compress_last=False):
	"""Write each mixed_ref record to its own file; returns the paths."""
	paths = []
	for i, rec in enumerate(mixed_ref):
		gz = compress_last and i == len(mixed_ref) - 1
		path = str(tmp_path / ('part%d.fa%s' % (i, '.gz' if gz else '')))
		paths.append(write_text(path, fasta_text([rec], 60), compress=gz))
	return paths


@pytest.mark.parametrize('compress_last', [False, True], ids=['plain', 'with_gz'])
def test_comma_separated_fasta_list(bt2, mixed_ref, mixed_index, tmp_path,
		compress_last):
	"""A comma-separated FASTA list (plain or gzipped) indexes like one multi-FASTA file in that order."""
	paths = split_reference(mixed_ref, tmp_path, compress_last)
	prefix = str(tmp_path / 'idx')
	bt2.build_index(','.join(paths), prefix)
	assert_same_files(prefix, mixed_index)


def test_space_separated_fasta_list(bt2, mixed_ref, mixed_index, tmp_path):
	"""FASTA files as separate arguments (undocumented wrapper alias) index like the comma list."""
	paths = split_reference(mixed_ref, tmp_path)
	prefix = str(tmp_path / 'idx')
	bt2.run('bowtie2-build', ['--quiet'] + paths + [prefix])
	assert_same_files(prefix, mixed_index)


def test_cmdline_sequences(bt2, genome_dict, tmp_path):
	"""-c indexes command-line sequences like a FASTA file of them named 0, 1, 2 (naming pins current output)."""
	seqs = [genome_dict['chr1'][100:400], genome_dict['chr2'][50:250].lower(),
		genome_dict['chr3'][1000:1180]]
	prefix = str(tmp_path / 'cmd')
	bt2.run('bowtie2-build', ['-q', '-c', ','.join(seqs), prefix])
	fasta = h.write_fasta(str(tmp_path / 'cmd.fa'),
		[(str(i), s) for i, s in enumerate(seqs)])
	bt2.build_index(fasta, str(tmp_path / 'fa'))
	assert_same_files(prefix, str(tmp_path / 'fa'))
	out = bt2.run('bowtie2-inspect', [prefix]).stdout
	assert parse_fasta(out) == [(str(i), s.upper()) for i, s in enumerate(seqs)]
	# A read from the second sequence aligns to record '1' at its offset.
	read = genome_dict['chr2'][100:160]
	fq = h.write_fastq(str(tmp_path / 'r.fq'), [('r', read)])
	_, records, _ = bt2.align(['-p', '1', '-x', prefix, '-U', fq])
	assert [(r.rname, r.pos, r.flag) for r in records] == [('1', 51, 0)]


# ---------------------------------------------------------------------------
# Odd reference records

@pytest.fixture(scope='module')
def odd_records(genome_dict):
	"""Records covering duplicate names, a name with a description, an
	all-N sequence, every IUPAC code in both cases, a gap, and an empty
	record (which bowtie2-build skips)."""
	g = genome_dict['chr1']
	return [
		('first', g[0:80]),
		('empty', ''),
		('dup', g[200:250]),
		('allN', 'N' * 12),
		('dup', g[300:340]),
		('with a description', g[400:470]),
		('codes', g[500:520] + IUPAC + IUPAC.lower() + '-' + g[520:550]),
		('last', g[600:633]),
	]


@pytest.fixture(scope='module')
def odd_index(bt2, odd_records, tmp_path_factory):
	d = tmp_path_factory.mktemp('odd')
	fasta = h.write_fasta(str(d / 'odd.fa'), odd_records)
	proc = bt2.build_index(fasta, str(d / 'odd'))
	return str(d / 'odd'), proc


def test_empty_record_is_skipped_with_warning(odd_index):
	"""An empty FASTA record is skipped with a warning (pins current behaviour; the MANUAL is silent)."""
	_, proc = odd_index
	assert 'Warning: Encountered empty reference sequence' in proc.stderr


def test_inspect_round_trips_odd_records(bt2, odd_records, odd_index):
	"""Default inspect prints each non-empty record under its full name, upper-cased, IUPAC/gap/N as N."""
	prefix, _ = odd_index
	expected = [(name, inspected(seq)) for name, seq in odd_records if seq]
	out = bt2.run('bowtie2-inspect', [prefix]).stdout
	assert_same_text(out, fasta_text(expected, 60))


def test_inspect_names_and_summary_of_odd_records(bt2, odd_records, odd_index):
	"""-n keeps input order and duplicate names; -s lengths count N, IUPAC and gap positions."""
	prefix, _ = odd_index
	kept = [(name, seq) for name, seq in odd_records if seq]
	names = bt2.run('bowtie2-inspect', ['-n', prefix]).stdout
	assert names.splitlines() == [name for name, _ in kept]
	summary = bt2.run('bowtie2-inspect', ['-s', prefix]).stdout
	assert sequence_lines(summary) == [(i + 1, name, len(seq))
		for i, (name, seq) in enumerate(kept)]


# ---------------------------------------------------------------------------
# bowtie2-inspect

def test_inspect_default_fasta(bt2, mixed_ref, mixed_index):
	"""Default inspect prints 60-column FASTA, upper-cased, with N runs, IUPAC codes and gaps as N."""
	out = bt2.run('bowtie2-inspect', [mixed_index]).stdout
	expected = [(name, inspected(seq)) for name, seq in mixed_ref]
	assert_same_text(out, fasta_text(expected, 60))


@pytest.mark.parametrize('flag', ['-a', '--across'])
@pytest.mark.parametrize('width', [1, 7, 60, 61, 100000])
def test_inspect_across(bt2, mixed_ref, mixed_index, flag, width):
	"""-a/--across N wraps the FASTA output every N bases."""
	out = bt2.run('bowtie2-inspect', [flag, width, mixed_index]).stdout
	expected = [(name, inspected(seq)) for name, seq in mixed_ref]
	assert_same_text(out, fasta_text(expected, width))


@pytest.mark.parametrize('flag', ['-n', '--names'])
def test_inspect_names(bt2, genome, genome_index, flag):
	"""-n/--names prints the reference names, one per line, and nothing else."""
	out = bt2.run('bowtie2-inspect', [flag, genome_index]).stdout
	assert out == ''.join(name + '\n' for name, _ in genome)


@pytest.mark.parametrize('flag', ['-s', '--summary'])
def test_inspect_summary(bt2, mixed_ref, mixed_fasta, tmp_path, flag):
	"""-s/--summary gives SA-Sample 1 in 2^offrate, FTab-Chars, and Sequence-<i> <name> <length> rows."""
	prefix = str(tmp_path / 'idx')
	bt2.build_index(mixed_fasta, prefix, extra=['-o', '3', '-t', '7'])
	out = bt2.run('bowtie2-inspect', [flag, prefix]).stdout
	assert summary_value(out, 'SA-Sample') == ['1 in 8']
	assert summary_value(out, 'FTab-Chars') == ['7']
	assert sequence_lines(out) == [(i + 1, name, len(seq))
		for i, (name, seq) in enumerate(mixed_ref)]


@pytest.mark.parametrize('flag', ['-o', '--output'])
def test_inspect_output_file(bt2, mixed_index, mixed_baseline, tmp_path, flag):
	"""-o/--output writes the FASTA to the file and nothing to stdout."""
	path = str(tmp_path / 'out.fa')
	proc = bt2.run('bowtie2-inspect', [flag, path, mixed_index])
	assert proc.stdout == ''
	with open(path) as f:
		assert_same_text(f.read(), mixed_baseline['fasta'])


def test_inspect_output_file_with_across(bt2, mixed_ref, mixed_index, tmp_path):
	"""-o combined with -a writes the re-wrapped FASTA to the file."""
	path = str(tmp_path / 'out.fa')
	bt2.run('bowtie2-inspect', ['-a', '25', '-o', path, mixed_index])
	with open(path) as f:
		assert_same_text(f.read(), fasta_text([(n, inspected(s)) for n, s in mixed_ref], 25))


def test_inspect_verbose(bt2, genome, genome_index):
	"""-v prints settings before the unchanged -n output (pins the first line's wording)."""
	quiet = bt2.run('bowtie2-inspect', ['-n', genome_index])
	proc = bt2.run('bowtie2-inspect', ['-v', '-n', genome_index])
	assert proc.stdout.startswith('Input ebwt file: "%s"' % genome_index)
	assert proc.stdout.endswith(quiet.stdout)
	assert len(proc.stdout) > len(quiet.stdout)


def test_inspect_wrapper_verbose_logs_binary(bt2, genome, genome_index):
	"""Wrapper --verbose logs the inspect binary it runs and leaves the -n output in place (pins the log wording)."""
	# Only the log line and the -n output are checked: the MANUAL's
	# "-v/--verbose" says --verbose should also print the -v settings, which
	# the current version does not do.
	proc = bt2.run('bowtie2-inspect', ['--verbose', '-n', genome_index])
	assert logged_binary(proc.stderr) == 'bowtie2-inspect-s'
	assert proc.stdout.endswith(''.join(name + '\n' for name, _ in genome))


def test_inspect_large_index_detected(bt2, mixed_index_large, mixed_baseline):
	"""A .bt2l-only index is found and inspects (FASTA, -n, -s) exactly like the .bt2 index."""
	assert_same_text(bt2.run('bowtie2-inspect', [mixed_index_large]).stdout, mixed_baseline['fasta'])
	assert bt2.run('bowtie2-inspect', ['-n', mixed_index_large]).stdout == mixed_baseline['names']
	assert bt2.run('bowtie2-inspect', ['-s', mixed_index_large]).stdout == mixed_baseline['summary']


def test_inspect_large_index_option(bt2, mixed_fasta, tmp_path):
	"""With .bt2 (-o 3) and .bt2l (-o 6) under one basename, --large-index selects the .bt2l."""
	prefix = str(tmp_path / 'both')
	bt2.build_index(mixed_fasta, prefix, extra=['-o', '3'])
	bt2.build_index(mixed_fasta, prefix, extra=['-o', '6', '--large-index'])
	small = bt2.run('bowtie2-inspect', ['--verbose', '-s', prefix])
	large = bt2.run('bowtie2-inspect', ['--verbose', '--large-index', '-s', prefix])
	assert summary_value(small.stdout, 'SA-Sample') == ['1 in 8']
	assert summary_value(large.stdout, 'SA-Sample') == ['1 in 64']
	assert logged_binary(small.stderr) == 'bowtie2-inspect-s'
	assert logged_binary(large.stderr) == 'bowtie2-inspect-l'


@pytest.mark.parametrize('kind', ['debug', 'sanitized'])
@pytest.mark.parametrize('large', [False, True], ids=['small', 'large'])
def test_inspect_debug_sanitized(bt2, mixed_index, mixed_index_large,
		mixed_baseline, kind, large):
	"""--debug/--sanitized run the matching inspect binary, whose output equals the release binary's."""
	binary = 'bowtie2-inspect-%s-%s' % ('l' if large else 's', kind)
	if not bt2.has_binary(binary):
		pytest.skip('%s not built' % binary)
	prefix = mixed_index_large if large else mixed_index
	extra = ['--' + kind] + (['--large-index'] if large else [])
	proc = bt2.run('bowtie2-inspect', extra + ['--verbose', prefix])
	assert logged_binary(proc.stderr) == binary
	assert_same_text(proc.stdout, mixed_baseline['fasta'])
	for flag, key in (('-n', 'names'), ('-s', 'summary')):
		assert bt2.run('bowtie2-inspect', extra + [flag, prefix]).stdout == mixed_baseline[key]


@pytest.mark.parametrize('flag', ['-n', '-s'])
def test_inspect_finds_index_through_bowtie2_indexes(bt2, mixed_index,
		mixed_baseline, tmp_path, flag):
	"""A basename not in the current directory is found through $BOWTIE2_INDEXES."""
	env = dict(os.environ)
	env['BOWTIE2_INDEXES'] = os.path.dirname(mixed_index)
	proc = bt2.run('bowtie2-inspect', [flag, os.path.basename(mixed_index)],
		cwd=str(tmp_path), env=env)
	assert proc.stdout == mixed_baseline['names' if flag == '-n' else 'summary']


def test_inspect_missing_index(bt2, tmp_path):
	"""Inspecting a basename with no index fails with bowtie2's message."""
	missing = str(tmp_path / 'nothing')
	proc = bt2.run('bowtie2-inspect', [missing], check=False)
	assert proc.returncode != 0
	assert 'Could not locate a Bowtie index corresponding to basename "%s"' % missing \
		in proc.stderr


@pytest.mark.parametrize('flag', ['-h', '--help'])
def test_inspect_help(bt2, flag):
	"""-h/--help prints usage listing every documented option and exits 0."""
	proc = bt2.run('bowtie2-inspect', [flag])
	assert proc.stdout.startswith('Bowtie 2 version %s' % version())
	assert 'Usage: bowtie2-inspect [options]* <bt2_base>' in proc.stdout
	for opt in ('--large-index', '--debug', '--sanitized', '--verbose',
			'-a/--across <int>', '-n/--names', '-s/--summary', '-o/--output',
			'-v/--verbose', '-h/--help'):
		assert opt in proc.stdout, opt


def test_inspect_version(bt2):
	"""--version prints the bowtie2-inspect version and exits 0."""
	proc = bt2.run('bowtie2-inspect', ['--version'])
	lines = proc.stdout.splitlines()
	assert re.search(r'bowtie2-inspect-s version %s$' % re.escape(version()), lines[0])
	assert '64-bit' in lines


# ---------------------------------------------------------------------------
# The shipped lambda index

def test_inspect_shipped_lambda_fasta(bt2, lambda_fasta, lambda_index):
	"""Inspecting the shipped lambda index gives back example/reference/lambda_virus.fa."""
	with open(lambda_fasta) as f:
		expected = parse_fasta(f.read())
	out = bt2.run('bowtie2-inspect', [lambda_index]).stdout
	assert parse_fasta(out) == [(n, inspected(s)) for n, s in expected]
	assert all(len(line) <= 60 for line in out.splitlines()[1:])


def test_inspect_shipped_lambda_names_and_summary(bt2, lambda_fasta, lambda_index):
	"""Shipped lambda index: -n gives the name line; -s gives length 48502, SA-Sample 1 in 128, 10 ftab chars."""
	with open(lambda_fasta) as f:
		[(name, seq)] = parse_fasta(f.read())
	assert len(seq) == 48502
	assert bt2.run('bowtie2-inspect', ['-n', lambda_index]).stdout == name + '\n'
	summary = bt2.run('bowtie2-inspect', ['-s', lambda_index]).stdout
	assert summary_value(summary, 'SA-Sample') == ['1 in 128']
	assert summary_value(summary, 'FTab-Chars') == ['10']
	assert sequence_lines(summary) == [(1, name, len(seq))]


# ---------------------------------------------------------------------------
# bowtie2-build wrapper options, help, version, quiet, errors

@pytest.mark.parametrize('extra,binary', [
	([], 'bowtie2-build-s'),
	(['--large-index'], 'bowtie2-build-l'),
	(['--debug'], 'bowtie2-build-s-debug'),
	(['--sanitized'], 'bowtie2-build-s-sanitized'),
	(['--debug', '--large-index'], 'bowtie2-build-l-debug'),
	(['--sanitized', '--large-index'], 'bowtie2-build-l-sanitized'),
], ids=arg_id)
def test_build_wrapper_binary_selection(bt2, mixed_fasta, mixed_index,
		mixed_index_large, tmp_path, extra, binary):
	"""--large-index/--debug/--sanitized pick the logged binary, and each writes the release index bytes."""
	if not bt2.has_binary(binary):
		pytest.skip('%s not built' % binary)
	prefix = str(tmp_path / 'idx')
	proc = bt2.build_index(mixed_fasta, prefix, extra=['--verbose'] + extra)
	assert logged_binary(proc.stderr) == binary
	if '--large-index' in extra:
		assert files_written(prefix) == index_names('bt2l')
		assert_same_files(prefix, mixed_index_large, ext='bt2l')
	else:
		assert files_written(prefix) == index_names('bt2')
		assert_same_files(prefix, mixed_index)


@pytest.mark.parametrize('flag', ['-q', '--quiet'])
def test_build_quiet(bt2, mixed_fasta, mixed_index, tmp_path, flag):
	"""-q/--quiet silences a clean build; without it the log names the input; the index is the same."""
	quiet = str(tmp_path / 'quiet')
	proc = bt2.run('bowtie2-build', [flag, mixed_fasta, quiet])
	assert proc.stdout == '' and proc.stderr == ''
	loud = str(tmp_path / 'loud')
	proc = bt2.run('bowtie2-build', [mixed_fasta, loud])
	assert '  %s\n' % mixed_fasta in proc.stdout
	assert_same_files(quiet, mixed_index)
	assert_same_files(loud, mixed_index)


BUILD_OPTIONS_IN_HELP = ('-f', '-c', '--large-index', '--debug', '--sanitized',
	'--verbose', '-a/--noauto', '-p/--packed', '--bmax <int>', '--bmaxdivn <int>',
	'--dcv <int>', '--nodc', '-r/--noref', '-3/--justref', '-o/--offrate <int>',
	'-t/--ftabchars <int>', '--threads <int>', '--seed <int>', '-q/--quiet',
	'--help', '--version')


@pytest.mark.parametrize('flag', ['-h', '--help'])
def test_build_help(bt2, tmp_path, flag):
	"""-h/--help prints usage listing every documented option, exits 0 and writes nothing."""
	proc = bt2.run('bowtie2-build', [flag], cwd=str(tmp_path))
	assert proc.stdout.startswith('Bowtie 2 version %s' % version())
	assert 'Usage: bowtie2-build [options]* <reference_in> <bt2_index_base>' in proc.stdout
	for opt in BUILD_OPTIONS_IN_HELP:
		assert opt in proc.stdout, opt
	assert os.listdir(str(tmp_path)) == []


def test_build_version(bt2, tmp_path):
	"""--version prints the bowtie2-build version, exits 0, writes nothing."""
	proc = bt2.run('bowtie2-build', ['--version'], cwd=str(tmp_path))
	lines = proc.stdout.splitlines()
	assert re.search(r'bowtie2-build-s version %s$' % re.escape(version()), lines[0])
	assert '64-bit' in lines
	assert os.listdir(str(tmp_path)) == []


@pytest.mark.parametrize('args,message', [
	(['-q'], 'No input sequence or sequence file specified!'),
	(['-q', '{fa}'], 'No output file specified!'),
	(['-q', '{missing}', '{prefix}'], 'Error: could not open'),
	(['-q', '{empty}', '{prefix}'], 'Warning: All fasta inputs were empty'),
	(['-q', '--bmax', '0', '{fa}', '{prefix}'], '--bmax arg must be at least 1'),
	(['-q', '--bmaxdivn', '0', '{fa}', '{prefix}'], '--bmaxdivn arg must be at least 1'),
	(['-q', '--dcv', '2', '{fa}', '{prefix}'], '--dcv arg must be at least 3'),
	(['-q', '-o', '-1', '{fa}', '{prefix}'], '-o/--offRate arg must be at least 0'),
	(['-q', '-t', '0', '{fa}', '{prefix}'], '-t/--ftabChars arg must be at least 1'),
	(['-q', '-t', '17', '{fa}', '{prefix}'], '-t/--ftabChars arg must not exceed 16'),
	(['-q', '--seed', '-1', '{fa}', '{prefix}'], '--seed arg must be at least 0'),
	(['-q', '--no-such-option', '{fa}', '{prefix}'], 'Usage: bowtie2-build'),
], ids=arg_id)
def test_build_errors(bt2, mixed_fasta, tmp_path, args, message):
	"""Bad invocations exit non-zero with bowtie2-build's message and leave no index files."""
	work = tmp_path / 'work'
	work.mkdir()
	empty = str(tmp_path / 'empty.fa')
	open(empty, 'w').close()
	subs = dict(fa=mixed_fasta, missing=str(tmp_path / 'missing.fa'), empty=empty,
		prefix=str(work / 'idx'))
	proc = bt2.run('bowtie2-build', [a.format(**subs) for a in args], check=False)
	assert proc.returncode != 0
	assert message in proc.stderr
	assert os.listdir(str(work)) == []
