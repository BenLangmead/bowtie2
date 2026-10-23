"""bowtie2 output options, performance options, the Perl wrapper and errors.

Covers -S, -t, --un/--al (plain and compressed), --quiet, the metrics options,
the SAM options (--no-unal, --no-hd, --no-sq, --rg-id, --rg, --omit-sec-seq,
--xeq, --sam-opt-config), header contents, the unpaired alignment summary,
-p/--reorder/--mm, OMP_NUM_THREADS, --version/-h, error exits, the wrapper-only
options (--debug, --sanitized, --verbose, --keep, --temp-directory,
--no-named-pipes, --log-file, --large-index, --ref-string, --bam) and index
auto-detection. Reads are cut from the fixed random `genome` (and from two
small references built here) at known places, so positions, strands, CIGARs,
MD strings and scores are known before bowtie2 runs.
"""

import bz2
import gzip
import os
import random
import re
import shutil
import stat
import subprocess
import sys
import threading
import time

import pytest

import bt2_helpers as h


# ---------------------------------------------------------------- helpers --

def version_string(bt2):
	"""The version bowtie2 was built with, read from BOWTIE2_VERSION."""
	for d in (bt2.bin_dir, h.ROOT):
		path = os.path.join(d, 'BOWTIE2_VERSION')
		if os.path.exists(path):
			with open(path) as f:
				return f.read().strip()
	pytest.skip('BOWTIE2_VERSION not found')


def body(text):
	"""SAM text without header lines."""
	return [line for line in text.splitlines() if line and not line.startswith('@')]


def headers(text):
	return [line for line in text.splitlines() if line.startswith('@')]


def opt_fields(line):
	"""Optional fields of a SAM line as (tag, type, value) triples. Values stay
	strings because some bowtie2 debug fields print NA in an i-typed field."""
	return [tuple(f.split(':', 2)) for f in line.split('\t')[11:] if f]


def opt_dict(line):
	return dict((t, v) for t, _, v in opt_fields(line))


def fastq_records(text):
	"""Split FASTQ text into raw 4-line records, keeping every line verbatim."""
	lines = text.splitlines()
	assert len(lines) % 4 == 0, text
	return [tuple(lines[i:i + 4]) for i in range(0, len(lines), 4)]


def fasta_records(text):
	lines = text.splitlines()
	assert len(lines) % 2 == 0, text
	return [tuple(lines[i:i + 2]) for i in range(0, len(lines), 2)]


def qname(read_name):
	"""SAM QNAME of a read: its name up to the first whitespace."""
	return read_name.split()[0]


def read_output(path, kind):
	"""Decompress a file written by --un/--al and friends."""
	if kind == 'gz':
		with gzip.open(path, 'rt') as f:
			return f.read()
	if kind == 'bz2':
		with bz2.open(path, 'rt') as f:
			return f.read()
	if kind in ('zst', 'lz4'):
		tool = {'zst': 'zstd', 'lz4': 'lz4'}[kind]
		return subprocess.run([tool, '-dc', path], capture_output=True, text=True,
			check=True).stdout
	with open(path) as f:
		return f.read()


def align_bin(bt2, name='bowtie2-align-s'):
	return os.path.join(bt2.bin_dir, name)


def env_without_omp(**extra):
	env = dict(os.environ)
	env.pop('OMP_NUM_THREADS', None)
	env.update(extra)
	return env


def metrics_rows(text):
	"""(header columns, data rows) of bowtie2 metrics output in `text`."""
	cols, rows = None, []
	for line in text.splitlines():
		if line.startswith('Time\tRead\tBase'):
			cols = line.rstrip('\t').split('\t')
		elif cols is not None and re.match(r'^\d+\t\d+\t', line):
			rows.append(line.rstrip('\t').split('\t'))
	return cols, rows


def pct(num, den):
	return '%.2f%%' % (100.0 * num / den if den else 0.0)


def unescape_passthrough(line):
	return re.sub(r'%([0-9A-Fa-f]{2})', lambda m: chr(int(m.group(1), 16)), line)


# ------------------------------------------------------------- read sets --

def varied_qual(n):
	"""A quality string cycling through Phred 2..40."""
	return ''.join(chr(33 + 2 + (i * 7) % 39) for i in range(n))


@pytest.fixture(scope='module')
def unp(genome_dict, tmp_path_factory):
	"""Unpaired reads with known placements, written to one FASTQ.

	Returns a dict with 'reads' (name, seq, qual) as written, 'fq' (path),
	'aligned' {qname: expected placement} and 'unaligned' (qnames)."""
	c1, c2, c3 = genome_dict['chr1'], genome_dict['chr2'], genome_dict['chr3']
	rng = random.Random(31337)
	reads = [
		('fw_exact first comment', c1[100:150], 'I' * 50),
		('rc_exact', h.revcomp(c2[500:560]), varied_qual(60)),
		('fw_mm', h.mutate(c1[2000:2060], [7, 50]), 'I' * 60),
		('rc_mm', h.revcomp(h.mutate(c3[1000:1060], [7, 50])), 'I' * 60),
		('unal_1 second comment', h.random_dna(50, rng), varied_qual(50)),
		('unal_2', h.random_dna(55, rng), 'I' * 55),
		('filt_ns', 'ACGT' * 5 + 'N' * 30, 'I' * 50),
	]
	aligned = {
		'fw_exact': dict(rname='chr1', pos=101, rev=False, length=50, md='50', nmm=0),
		'rc_exact': dict(rname='chr2', pos=501, rev=True, length=60, md='60', nmm=0),
		'fw_mm': dict(rname='chr1', pos=2001, rev=False, length=60,
			md='7%s42%s9' % (c1[2007], c1[2050]), nmm=2),
		'rc_mm': dict(rname='chr3', pos=1001, rev=True, length=60,
			md='7%s42%s9' % (c3[1007], c3[1050]), nmm=2),
	}
	path = str(tmp_path_factory.mktemp('unp') / 'unp.fq')
	h.write_fastq(path, reads)
	return dict(reads=reads, fq=path, aligned=aligned,
		unaligned=['unal_1', 'unal_2', 'filt_ns'])


@pytest.fixture(scope='module')
def pair(genome_dict, tmp_path_factory):
	"""One concordant --fr pair: mate 1 forward with one mismatch (AS -6),
	mate 2 reverse with two mismatches (AS -12), fragment 1001-1300 on chr1."""
	c1 = genome_dict['chr1']
	d = tmp_path_factory.mktemp('pair')
	m1 = h.write_fastq(str(d / 'p_1.fq'), [('pr', h.mutate(c1[1000:1050], [25]))])
	m2 = h.write_fastq(str(d / 'p_2.fq'),
		[('pr', h.revcomp(h.mutate(c1[1250:1300], [10, 40])))])
	return dict(m1=m1, m2=m2)


@pytest.fixture(scope='module')
def rep(bt2, tmp_path_factory):
	"""A one-chromosome reference with an 80-bp unit repeated three times and a
	50-bp segment followed by a copy with six mismatches in its middle (seed
	0-21 stays exact), plus reads.

	The FASTA name carries a description, which bowtie2 drops from @SQ SN."""
	rng = random.Random(42)
	unit = h.random_dna(80, rng)
	seg = h.random_dna(50, rng)
	bg = [h.random_dna(300, rng) for _ in range(6)]
	parts = [bg[0], unit, bg[1], unit, bg[2], unit, bg[3], seg, bg[4],
		h.mutate(seg, [24, 27, 30, 33, 36, 39]), bg[5]]
	seq = ''.join(parts)
	offsets = [sum(len(p) for p in parts[:i]) for i in range(len(parts))]
	unit_starts = [offsets[1], offsets[3], offsets[5]]
	d = tmp_path_factory.mktemp('rep')
	fasta = h.write_fasta(str(d / 'rep.fa'), [('rep extra description', seq)])
	prefix = str(d / 'rep')
	bt2.build_index(fasta, prefix, extra=['--seed', '0'])
	reads = [
		('rep_fw', unit[10:60]),
		('rep_rc', h.revcomp(unit[20:70])),
		('uni_fw', bg[0][100:150]),
		('uni_rc', h.revcomp(bg[3][50:100])),
		('near', seg),
		('unal', h.random_dna(50, random.Random(99))),
	]
	fq = h.write_fastq(str(d / 'rep.fq'), reads)
	return dict(index=prefix, fq=fq, seq=seq, reads=reads,
		rep_fw_pos=sorted(s + 10 + 1 for s in unit_starts),
		rep_rc_pos=sorted(s + 20 + 1 for s in unit_starts),
		near_pos=offsets[7] + 1)


@pytest.fixture(scope='module')
def many(genome, tmp_path_factory):
	"""1200 reads (1100 cut from the genome, some with a mismatch, half reverse
	complemented; 100 random) in shuffled order, for the threading tests. At
	this size -p 4 output without --reorder came out of input order in each of
	30 runs pinned to 1, 2 and 4 cores, so test_reorder notices a --reorder
	that does nothing (very likely, though not certainly)."""
	rng = random.Random(4242)
	reads = []
	for i in range(1100):
		name, seq = genome[rng.randrange(len(genome))]
		length = rng.choice([40, 50, 75, 100])
		start = rng.randrange(len(seq) - length)
		s = seq[start:start + length]
		if rng.random() < 0.5:
			s = h.mutate(s, [rng.randrange(2, 8)])
		if rng.random() < 0.5:
			s = h.revcomp(s)
		reads.append(('m%04d' % i, s))
	for i in range(100):
		reads.append(('n%03d' % i, h.random_dna(60, rng)))
	rng.shuffle(reads)
	path = str(tmp_path_factory.mktemp('many') / 'many.fq')
	h.write_fastq(path, reads)
	return dict(fq=path, names=[r[0] for r in reads])


def unp_args(genome_index, unp):
	return ['-x', genome_index, '-U', unp['fq'], '-p', '1']


# ------------------------------------------------------------ -S / -t --

@pytest.mark.parametrize('opt', ['-S', '--output'])
def test_S_writes_sam_to_file(bt2, genome_index, unp, tmp_path, opt):
	"""-S/--output sends the SAM to the file and nothing to stdout."""
	ref = bt2.run('bowtie2', unp_args(genome_index, unp)).stdout
	out = str(tmp_path / 'out.sam')
	proc = bt2.run('bowtie2', unp_args(genome_index, unp) + [opt, out])
	assert proc.stdout == ''
	with open(out) as f:
		text = f.read()
	assert body(text) == body(ref)
	assert [l for l in headers(text) if not l.startswith('@PG')] == \
		[l for l in headers(ref) if not l.startswith('@PG')]
	assert len(body(text)) == len(unp['reads'])


def test_time_prints_phase_times(bt2, genome_index, unp):
	"""-t/--time adds hh:mm:ss timings for index loading and the whole run;
	the labels pin the current output."""
	plain = bt2.run('bowtie2', unp_args(genome_index, unp))
	for opt in ('-t', '--time'):
		timed = bt2.run('bowtie2', unp_args(genome_index, unp) + [opt])
		lines = timed.stderr.splitlines()
		assert any(re.match(r'^Time loading .+: \d\d:\d\d:\d\d$', l) for l in lines), timed.stderr
		assert any(re.match(r'^Overall time: \d\d:\d\d:\d\d$', l) for l in lines), timed.stderr
		assert body(timed.stdout) == body(plain.stdout)
	assert not re.search(r'\d\d:\d\d:\d\d', plain.stderr)


# ---------------------------------------------------------- --un / --al --

COMPRESSION = [
	('', None, None),
	('-gz', 'gz', 'gzip'),
	('-bz2', 'bz2', 'bzip2'),
	('-zst', 'zst', 'zstd'),
	('-lz4', 'lz4', 'lz4'),
]


@pytest.mark.parametrize('suffix,kind,tool', COMPRESSION,
	ids=['plain', 'gz', 'bz2', 'zst', 'lz4'])
def test_un_al_split_reads_by_alignment(bt2, genome_index, unp, tmp_path, suffix, kind, tool):
	"""--un/--al (and their compressed forms) receive the unaligned/aligned
	reads exactly as they were in the input, and the SAM is unchanged."""
	if tool is not None and shutil.which(tool) is None:
		pytest.skip('%s not installed' % tool)
	un, al, sam = str(tmp_path / 'un.out'), str(tmp_path / 'al.out'), str(tmp_path / 'o.sam')
	bt2.run('bowtie2', unp_args(genome_index, unp) + ['--un' + suffix, un,
		'--al' + suffix, al, '-S', sam])
	records = fastq_records(h.fastq_text(unp['reads']))
	by_name = dict((qname(r[0][1:]), r) for r in records)
	got_un = fastq_records(read_output(un, kind))
	got_al = fastq_records(read_output(al, kind))
	assert sorted(got_un) == sorted(by_name[n] for n in unp['unaligned'])
	assert sorted(got_al) == sorted(by_name[n] for n in unp['aligned'])
	with open(sam) as f:
		assert body(f.read()) == body(bt2.run('bowtie2', unp_args(genome_index, unp)).stdout)


@pytest.mark.parametrize('variant', ['trim', 'phred64'])
def test_un_al_keep_original_record(bt2, genome_index, unp, tmp_path, variant):
	"""--un/--al write the untrimmed read with its original quality encoding."""
	reads, extra = unp['reads'], []
	if variant == 'trim':
		extra = ['-5', '3', '-3', '10']
	else:
		reads = [(n, s, ''.join(chr(ord(c) + 31) for c in q)) for n, s, q in reads]
		extra = ['--phred64']
	fq = h.write_fastq(str(tmp_path / 'in.fq'), reads)
	un, al = str(tmp_path / 'un.fq'), str(tmp_path / 'al.fq')
	_, recs, _ = bt2.align(['-x', genome_index, '-U', fq, '-p', '1', '--un', un, '--al', al] + extra)
	assert sorted(r.qname for r in recs if r.is_unmapped) == sorted(unp['unaligned'])
	records = dict((qname(r[0][1:]), r) for r in fastq_records(h.fastq_text(reads)))
	with open(un) as f:
		assert sorted(fastq_records(f.read())) == sorted(records[n] for n in unp['unaligned'])
	with open(al) as f:
		assert sorted(fastq_records(f.read())) == sorted(records[n] for n in unp['aligned'])


def test_un_al_fasta_input(bt2, genome_index, unp, tmp_path):
	"""With -f input, --un/--al write the FASTA records as they were input."""
	fa = str(tmp_path / 'in.fa')
	with open(fa, 'w') as f:
		for name, seq, _ in unp['reads']:
			f.write('>%s\n%s\n' % (name, seq))
	un, al = str(tmp_path / 'un.fa'), str(tmp_path / 'al.fa')
	bt2.run('bowtie2', ['-x', genome_index, '-f', '-U', fa, '-p', '1', '--un', un,
		'--al', al, '-S', str(tmp_path / 'o.sam')])
	with open(fa) as f:
		records = dict((qname(r[0][1:]), r) for r in fasta_records(f.read()))
	with open(un) as f:
		assert sorted(fasta_records(f.read())) == sorted(records[n] for n in unp['unaligned'])
	with open(al) as f:
		assert sorted(fasta_records(f.read())) == sorted(records[n] for n in unp['aligned'])


def test_un_al_directory_path(bt2, genome_index, unp, tmp_path):
	"""--un/--al given a directory write un-seqs/al-seqs inside it (wrapper
	behaviour; the file names pin the current output)."""
	d = tmp_path / 'outdir'
	d.mkdir()
	bt2.run('bowtie2', unp_args(genome_index, unp) + ['--un', str(d), '--al', str(d),
		'-S', str(tmp_path / 'o.sam')])
	assert sorted(os.listdir(str(d))) == ['al-seqs', 'un-seqs']
	with open(str(d / 'un-seqs')) as f:
		assert sorted(qname(r[0][1:]) for r in fastq_records(f.read())) == sorted(unp['unaligned'])
	with open(str(d / 'al-seqs')) as f:
		assert sorted(qname(r[0][1:]) for r in fastq_records(f.read())) == sorted(unp['aligned'])


# -------------------------------------------------------------- --quiet --

def test_quiet_suppresses_summary(bt2, genome_index, unp):
	"""--quiet leaves stderr empty on success and does not change the SAM."""
	plain = bt2.run('bowtie2', unp_args(genome_index, unp))
	quiet = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--quiet'])
	assert 'overall alignment rate' in plain.stderr
	assert quiet.stderr == ''
	assert body(quiet.stdout) == body(plain.stdout)


# ------------------------------------------------------------- metrics --

def check_metrics(text, unp):
	cols, rows = metrics_rows(text)
	assert cols is not None, text
	assert rows, text
	final = dict(zip(cols, rows[-1]))
	assert int(final['Read']) == len(unp['reads'])
	assert int(final['Base']) == sum(len(r[1]) for r in unp['reads'])
	assert int(final['UnfilteredRead']) == len(unp['reads']) - 1   # filt_ns


def test_met_file(bt2, genome_index, unp, tmp_path):
	"""--met-file writes a metrics table whose final row counts the reads and
	bases; nothing goes to stderr and the SAM is unchanged."""
	met = str(tmp_path / 'met.txt')
	proc = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--met-file', met])
	with open(met) as f:
		check_metrics(f.read(), unp)
	assert 'Time\tRead' not in proc.stderr
	assert body(proc.stdout) == body(bt2.run('bowtie2', unp_args(genome_index, unp)).stdout)


def test_met_stderr_and_file_together(bt2, genome_index, unp, tmp_path):
	"""--met-stderr prints the metrics table to stderr, also when --met-file
	is given (the two are not mutually exclusive)."""
	proc = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--met-stderr'])
	check_metrics(proc.stderr, unp)
	met = str(tmp_path / 'met.txt')
	proc = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--met-stderr', '--met-file', met])
	check_metrics(proc.stderr, unp)
	with open(met) as f:
		check_metrics(f.read(), unp)


def feed_copies(proc, data, stop, sent):
	"""Write `data` to proc's stdin again and again until `stop` is set, then
	close it; sent[0] counts the copies written in full."""
	try:
		while not stop.is_set():
			proc.stdin.write(data)
			proc.stdin.flush()
			sent[0] += 1
		proc.stdin.close()
	except (OSError, ValueError):
		pass   # the process died or was killed; the test reports it


def test_met_interval(bt2, genome_index, many, tmp_path):
	"""--met <int> sets the reporting period (default 1 s): fed reads for long
	enough that --met 1 and the default each print two interval rows before
	the final one, --met 100, fed the same way, prints only the final row,
	which counts every read."""
	# bowtie2 checks the clock only between reads, so an interval row needs
	# reads to be aligned at least <int> s after the search threads start. A
	# fixed sleep before sending the reads would depend on how fast bowtie2
	# starts on a loaded machine; instead each process is fed copies of the
	# reads until the --met 1 and default processes have both printed two
	# interval rows (a second or two of aligning), with a generous deadline.
	# Two rows are printed at different clock seconds, so the --met 100
	# process, fed alongside, has also been aligning across a second boundary.
	with open(many['fq'], 'rb') as f:
		data = f.read()
	ivals = ('1', 'default', '100')
	errs = dict((i, str(tmp_path / ('met_%s.txt' % i))) for i in ivals)
	procs, sent, feeders = {}, {}, []
	stop = threading.Event()

	def rows(ival):
		with open(errs[ival]) as f:
			return metrics_rows(f.read())[1]

	try:
		for ival in ivals:
			met = [] if ival == 'default' else ['--met', ival]
			with open(errs[ival], 'wb') as err:
				procs[ival] = subprocess.Popen([align_bin(bt2), '--wrapper', 'basic-0', '-x',
					genome_index, '-U', '-', '-p', '1', '--met-stderr', '--quiet', '-S',
					os.devnull] + met,
					stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=err,
					env=env_without_omp())
			sent[ival] = [0]
			t = threading.Thread(target=feed_copies, args=(procs[ival], data, stop, sent[ival]))
			t.daemon = True
			t.start()
			feeders.append(t)
		deadline = time.time() + 60
		while not all(len(rows(i)) >= 2 for i in ('1', 'default')):
			for ival, p in procs.items():
				assert p.poll() is None, 'bowtie2-align-s --met %s exited early' % ival
			assert time.time() < deadline, dict((i, len(rows(i))) for i in ivals)
			time.sleep(0.05)
		stop.set()
		for t in feeders:
			t.join(120)
		for ival, p in procs.items():
			assert p.wait(timeout=120) == 0, rows(ival)
	finally:
		stop.set()
		for p in procs.values():
			if p.poll() is None:
				p.kill()
				p.wait()
	final = dict((i, rows(i)) for i in ivals)
	assert len(final['100']) == 1
	assert len(final['1']) >= 3
	assert len(final['default']) >= 3
	for ival in ivals:
		assert sent[ival][0] >= 1
		assert int(final[ival][-1][1]) == sent[ival][0] * len(many['names'])


# ---------------------------------------------------------- SAM options --

def test_default_records_match_design(bt2, genome_index, unp):
	"""Baseline for the SAM option tests: each designed read lands where it
	was cut, with the documented default optional fields."""
	_, recs, _ = bt2.align(unp_args(genome_index, unp))
	assert [r.qname for r in recs] == [qname(r[0]) for r in unp['reads']]
	for r in recs:
		exp = unp['aligned'].get(r.qname)
		if exp is None:
			assert r.flag == h.UNMAPPED and r.rname == '*' and r.pos == 0
			assert sorted(r.tag_order) == (['YF', 'YT'] if r.qname == 'filt_ns' else ['YT'])
			continue
		assert (r.rname, r.pos, r.is_reverse) == (exp['rname'], exp['pos'], exp['rev'])
		assert r.cigar == '%dM' % exp['length']
		assert r.tags['MD'] == exp['md']
		assert r.tags['AS'] == -6 * exp['nmm']
		assert r.tags['XM'] == r.tags['NM'] == exp['nmm']
		assert sorted(r.tag_order) == ['AS', 'MD', 'NM', 'XG', 'XM', 'XN', 'XO', 'YT']
	assert recs[-1].tags['YF'] == 'NS'


@pytest.mark.parametrize('mode', ['wrapper', 'binary'])
def test_no_unal(bt2, genome_index, unp, tmp_path, mode):
	"""--no-unal drops every unaligned record and leaves the aligned ones as
	they were (the wrapper filters itself; the binary has its own path)."""
	ref = body(bt2.run('bowtie2', unp_args(genome_index, unp)).stdout)
	if mode == 'wrapper':
		out = str(tmp_path / 'o.sam')
		bt2.run('bowtie2', unp_args(genome_index, unp) + ['--no-unal', '-S', out])
		with open(out) as f:
			got = body(f.read())
	else:
		got = body(bt2.run(align_bin(bt2), unp_args(genome_index, unp) + ['--no-unal']).stdout)
	assert got == [l for l in ref if not int(l.split('\t')[1]) & h.UNMAPPED]
	assert sorted(l.split('\t')[0] for l in got) == sorted(unp['aligned'])


def test_no_unal_still_writes_un(bt2, genome_index, unp, tmp_path):
	"""--no-unal hides unaligned reads from the SAM but --un still gets them."""
	un = str(tmp_path / 'un.fq')
	proc = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--no-unal', '--un', un])
	assert sorted(l.split('\t')[0] for l in body(proc.stdout)) == sorted(unp['aligned'])
	with open(un) as f:
		assert sorted(qname(r[0][1:]) for r in fastq_records(f.read())) == sorted(unp['unaligned'])


NO_HD = ['--no-hd', '--no-head', '--no-HD', '--sam-no-hd', '--sam-no-head', '--sam-nohead',
	'--sam-noHD']
NO_SQ = ['--no-sq', '--no-SQ', '--sam-no-sq', '--sam-nosq', '--sam-noSQ']


@pytest.mark.parametrize('opt', NO_HD)
def test_no_hd(bt2, genome_index, unp, opt):
	"""--no-hd (and its aliases) removes every header line and nothing else."""
	ref = bt2.run('bowtie2', unp_args(genome_index, unp)).stdout
	out = bt2.run('bowtie2', unp_args(genome_index, unp) + [opt]).stdout
	assert headers(out) == []
	assert body(out) == body(ref)


@pytest.mark.parametrize('opt', NO_SQ)
def test_no_sq(bt2, genome_index, unp, opt):
	"""--no-sq (and its aliases) removes the @SQ lines but keeps @HD and @PG."""
	ref = bt2.run('bowtie2', unp_args(genome_index, unp)).stdout
	out = bt2.run('bowtie2', unp_args(genome_index, unp) + [opt]).stdout
	assert [l.split('\t')[0] for l in headers(out)] == ['@HD', '@PG']
	assert headers(out)[0] == headers(ref)[0]
	assert body(out) == body(ref)


def test_rg_id_and_rg(bt2, genome_index, unp):
	"""--rg-id makes an @RG line (ID first, then each --rg in order) and tags
	every record with RG:Z:<id>; @RG sitting between @SQ and @PG and RG:Z
	coming last pin the current output."""
	header, recs, _ = bt2.align(unp_args(genome_index, unp) + ['--rg-id', 'grp1',
		'--rg', 'SM:sample', '--rg', 'PL:ILLUMINA', '--rg', 'DS:two words'])
	assert [l.split('\t')[0] for l in header] == ['@HD', '@SQ', '@SQ', '@SQ', '@RG', '@PG']
	assert header[4] == '@RG\tID:grp1\tSM:sample\tPL:ILLUMINA\tDS:two words'
	assert len(recs) == len(unp['reads'])
	for r in recs:
		assert r.tags['RG'] == 'grp1', r
		assert r.tag_order[-1] == 'RG'


def test_rg_without_rg_id(bt2, genome_index, unp):
	"""--rg without --rg-id prints no @RG line and no RG:Z field."""
	ref = bt2.run('bowtie2', unp_args(genome_index, unp)).stdout
	out = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--rg', 'SM:x']).stdout
	assert not any(l.startswith('@RG') for l in headers(out))
	assert body(out) == body(ref)


@pytest.mark.parametrize('args', [
	['--rg', 'ID:grp1', '--rg', 'SM:s'],
	['--sam-rg-id', 'grp1', '--sam-rg', 'SM:s'],
	['--rg-id', 'grp1', '--sam-RG', 'SM:s'],
	['--rg-id', 'grp1', '--RG', 'SM:s'],
], ids=['rg-ID-form', 'sam-rg', 'sam-RG', 'RG'])
def test_rg_equivalent_spellings(bt2, genome_index, unp, args):
	"""--rg ID:<id> and the --sam-rg*/--RG aliases equal --rg-id/--rg."""
	ref = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--rg-id', 'grp1', '--rg', 'SM:s']).stdout
	out = bt2.run('bowtie2', unp_args(genome_index, unp) + args).stdout
	assert [l for l in headers(out) if not l.startswith('@PG')] == \
		[l for l in headers(ref) if not l.startswith('@PG')]
	assert body(out) == body(ref)


def test_no_hd_keeps_rg_field(bt2, genome_index, unp):
	"""--no-hd with --rg-id drops the @RG line but keeps RG:Z on records."""
	out = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--no-hd', '--rg-id', 'g']).stdout
	assert headers(out) == []
	assert all(opt_dict(l)['RG'] == 'g' for l in body(out))


@pytest.mark.parametrize('opt', ['--omit-sec-seq', '--sam-omit-sec-seq'])
def test_omit_sec_seq(bt2, rep, opt):
	"""--omit-sec-seq puts * in SEQ and QUAL of secondary (-k) records only."""
	args = ['-x', rep['index'], '-U', rep['fq'], '-p', '1', '-k', '3']
	_, full, _ = bt2.align(args)
	_, omit, _ = bt2.align(args + [opt])
	reads = dict(rep['reads'])
	for name, positions, rev in (('rep_fw', rep['rep_fw_pos'], False),
			('rep_rc', rep['rep_rc_pos'], True)):
		recs = h.by_name(full)[name]
		assert sorted(r.pos for r in recs) == positions
		assert [r.is_secondary for r in recs] == [False, True, True]
		assert all(r.is_reverse == rev for r in recs)
		seq = h.revcomp(reads[name]) if rev else reads[name]
		assert all(r.seq == seq and r.qual == 'I' * 50 for r in recs)
		orecs = h.by_name(omit)[name]
		assert orecs[0].seq == seq and orecs[0].qual == 'I' * 50
		assert [(r.seq, r.qual) for r in orecs[1:]] == [('*', '*'), ('*', '*')]
		assert [(r.flag, r.pos, r.cigar) for r in orecs] == [(r.flag, r.pos, r.cigar) for r in recs]
	# Reads without secondary alignments are untouched.
	assert [r.line for r in full if r.qname not in ('rep_fw', 'rep_rc')] == \
		[r.line for r in omit if r.qname not in ('rep_fw', 'rep_rc')]


def xeq_reads(genome_dict):
	"""Reads with mismatches and gaps placed so the left-aligned CIGAR is
	unique; returns (name, seq, rname, pos, xeq_cigar, m_cigar)."""
	c1, c3 = genome_dict['chr1'], genome_dict['chr3']
	out = [
		('x_fw', h.mutate(c1[2000:2060], [7, 50]), 'chr1', 2001, '7=1X42=1X9=', '60M'),
		('x_rc', h.revcomp(h.mutate(c3[1000:1060], [7, 50])), 'chr3', 1001, '7=1X42=1X9=', '60M'),
	]
	# 3-bp deletion at read offset a: unique when ref[a-1] != ref[a+2].
	st = 4000
	a = next(a for a in range(28, 40) if c1[st + a - 1] != c1[st + a + 2])
	out.append(('x_del', c1[st:st + a] + c1[st + a + 3:st + 63], 'chr1', st + 1,
		'%d=3D%d=' % (a, 60 - a), '%dM3D%dM' % (a, 60 - a)))
	# 2-bp insertion before ref offset b: unique when ins[-1] != ref[b-1] and
	# ins[0] != ref[b].
	st = 1500
	b = 30
	ins = next(x + y for x in 'ACGT' for y in 'ACGT'
		if y != c1[st + b - 1] and x != c1[st + b] and x != y)
	out.append(('x_ins', c1[st:st + b] + ins + c1[st + b:st + 58], 'chr1', st + 1,
		'%d=2I%d=' % (b, 28), '%dM2I%dM' % (b, 28)))
	return out


def test_xeq(bt2, genome_index, genome_dict, tmp_path):
	"""--xeq writes =/X instead of M at the placed mismatches; everything else
	in the record is unchanged."""
	design = xeq_reads(genome_dict)
	fq = h.write_fastq(str(tmp_path / 'x.fq'), [(d[0], d[1]) for d in design])
	args = ['-x', genome_index, '-U', fq, '-p', '1']
	_, plain, _ = bt2.align(args)
	_, xeq, _ = bt2.align(args + ['--xeq'])
	assert len(plain) == len(xeq) == len(design)
	for d, p, x in zip(design, plain, xeq):
		name, _, rname, pos, xcig, mcig = d
		assert (x.qname, x.rname, x.pos) == (name, rname, pos)
		assert x.cigar == xcig
		assert p.cigar == mcig
		assert x.line.split('\t')[:5] + x.line.split('\t')[6:] == \
			p.line.split('\t')[:5] + p.line.split('\t')[6:]


# ------------------------------------------------------- --sam-opt-config --

@pytest.fixture(scope='module')
def datasets(genome_index, unp, pair, rep):
	return {
		'unp': ['-x', genome_index, '-U', unp['fq'], '-p', '1'],
		'pair': ['-x', genome_index, '-1', pair['m1'], '-2', pair['m2'], '-p', '1'],
		'rep': ['-x', rep['index'], '-U', rep['fq'], '-p', '1'],
	}


_default_cache = {}


def default_lines(bt2, args):
	key = tuple(args)
	if key not in _default_cache:
		_default_cache[key] = body(bt2.run('bowtie2', args).stdout)
	return _default_cache[key]


def same_core(a, b):
	"""SAM lines a and b agree in the 11 mandatory fields."""
	assert a.split('\t')[:11] == b.split('\t')[:11]


@pytest.mark.parametrize('config,removed,data', [
	('-xn', ['XN'], 'unp'), ('-xm', ['XM'], 'unp'), ('-xo', ['XO'], 'unp'),
	('-xg', ['XG'], 'unp'), ('-nm', ['NM'], 'unp'), ('-md', ['MD'], 'unp'),
	('-yt', ['YT'], 'unp'), ('-yf', ['YF'], 'unp'), ('-ys', ['YS'], 'pair'),
	('-xs', ['XS'], 'rep'), ('-MD,-XM', ['MD', 'XM'], 'unp'),
])
def test_sam_opt_config_turns_off(bt2, datasets, config, removed, data):
	"""--sam-opt-config -<tag> removes that default field and nothing else."""
	ref = default_lines(bt2, datasets[data])
	out = body(bt2.run('bowtie2', datasets[data] + ['--sam-opt-config', config]).stdout)
	assert len(out) == len(ref)
	assert any(t in removed for l in ref for t, _, _ in opt_fields(l)), 'tag never printed'
	for a, b in zip(ref, out):
		same_core(a, b)
		assert opt_fields(b) == [f for f in opt_fields(a) if f[0] not in removed]


def check_yn(lines):
	for l in lines:
		f = opt_dict(l)
		length = len(l.split('\t')[9])
		assert int(f['YN']) == int(-0.6 + -0.6 * length), l
		assert f['Yn'] == '0', l


def check_present(*tags):
	def check(lines):
		for l in lines:
			f = opt_fields(l)
			for t in tags:
				assert [x for x in f if x[0] == t], (t, l)
	return check


def check_ints(*tags):
	def check(lines):
		for l in lines:
			f = opt_dict(l)
			for t in tags:
				assert re.match(r'^-?\d+$', f[t]), (t, l)
	return check


def aligned(lines):
	return [l for l in lines if not int(l.split('\t')[1]) & h.UNMAPPED]


def check_ym(lines):
	for l in aligned(lines):
		assert opt_dict(l)['YM'] == '0', l


def check_zu(lines):
	for l in aligned(lines):
		f = opt_dict(l)
		assert f['ZU'] == f['AS'] and 'Zu' in f, l
	for l in lines:
		if int(l.split('\t')[1]) & h.UNMAPPED:
			assert 'ZU' not in opt_dict(l), l


def check_zt(lines):
	for l in aligned(lines):
		f = opt_dict(l)
		assert f['ZT'].split(',')[0] == f['AS'], l


def check_xss(lines):
	near = [l for l in lines if l.split('\t')[0] == 'near']
	assert len(near) == 1
	f = opt_dict(near[0])
	assert f['AS'] == '0'
	assert f['Xs'] == '-36'   # the copy with six Q40 mismatches, below --score-min


def mates(lines):
	m1 = [l for l in lines if int(l.split('\t')[1]) & h.FIRST_IN_PAIR]
	m2 = [l for l in lines if int(l.split('\t')[1]) & h.SECOND_IN_PAIR]
	assert len(m1) == len(m2) == 1
	return m1[0], m2[0]


def check_mc(lines):
	m1, m2 = mates(lines)
	assert opt_dict(m1)['MC'] == m2.split('\t')[5] == '50M'
	assert opt_dict(m2)['MC'] == m1.split('\t')[5] == '50M'


def check_mq(lines):
	m1, m2 = mates(lines)
	assert opt_dict(m1)['MQ'] == m2.split('\t')[4]
	assert opt_dict(m2)['MQ'] == m1.split('\t')[4]


def check_yp(lines):
	for l in mates(lines):
		assert opt_dict(l)['YP'] == '0'


def check_zp(lines):
	m1, m2 = mates(lines)
	assert (opt_dict(m1)['AS'], opt_dict(m2)['AS']) == ('-6', '-12')
	for l in (m1, m2):
		assert opt_dict(l)['ZP'] == '-18'   # best concordant score: sum of mates' AS


def check_on_all(*tags):
	"""Each tag is on every record (aligned or not) with an integer value."""
	def check(lines):
		check_present(*tags)(lines)
		check_ints(*tags)(lines)
	return check


ADD_CASES = [
	('yn', ['YN', 'Yn'], 'unp', check_yn),
	('zs', ['ZS'], 'unp', check_on_all('ZS')),
	('ym', ['YM'], 'unp', check_ym),
	('xt', ['XT'], 'unp', check_on_all('XT')),
	('xd', ['XD', 'Xd'], 'unp', check_on_all('XD', 'Xd')),
	('xu', ['XU', 'Xu'], 'unp', check_on_all('XU', 'Xu')),
	('ye', ['YE', 'Ye'], 'unp', check_on_all('YE', 'Ye')),
	('yl', ['YL', 'Yl'], 'unp', check_on_all('YL', 'Yl')),
	('yu', ['YU', 'Yu'], 'unp', check_on_all('YU', 'Yu')),
	('xp', ['XP'], 'unp', check_present('XP')),
	('yr', ['YR'], 'unp', check_on_all('YR')),
	('zb', ['ZB'], 'unp', check_on_all('ZB')),
	('zr', ['ZR'], 'unp', check_present('ZR')),
	('zf', ['ZF', 'Zf'], 'unp', check_on_all('ZF', 'Zf')),
	('zm', ['ZM'], 'unp', check_present('ZM')),
	('zi', ['ZI'], 'unp', check_on_all('ZI')),
	('zu', ['ZU', 'Zu'], 'unp', check_zu),
	('zt', ['ZT'], 'unp', check_zt),
	('xss', ['Xs', 'Ys'], 'rep', check_xss),
	('mc', ['MC'], 'pair', check_mc),
	('mq', ['MQ'], 'pair', check_mq),
	('yp', ['YP'], 'pair', check_yp),
	('zp', ['ZP', 'Zp'], 'pair', check_zp),
	('YN,Zs', ['YN', 'Yn', 'ZS'], 'unp', check_present('YN', 'Yn', 'ZS')),
]


@pytest.mark.parametrize('config,added,data,check', ADD_CASES, ids=[c[0] for c in ADD_CASES])
def test_sam_opt_config_turns_on(bt2, datasets, config, added, data, check):
	"""--sam-opt-config <tag> adds that off-by-default field (values checked
	where they follow from the design) and changes nothing else."""
	ref = default_lines(bt2, datasets[data])
	out = body(bt2.run('bowtie2', datasets[data] + ['--sam-opt-config', config]).stdout)
	assert len(out) == len(ref)
	assert not any(t in added for l in ref for t, _, _ in opt_fields(l))
	for a, b in zip(ref, out):
		same_core(a, b)
		assert [f for f in opt_fields(b) if f[0] not in added] == opt_fields(a)
	check(out)


@pytest.mark.parametrize('config', ['x0', 'x1', 'yi', '-x0,-x1,-yi'])
def test_sam_opt_config_accepted_noop_tags(bt2, datasets, config):
	"""X0, X1 and YI are accepted but print nothing with the default MAPQ
	model (pins the current output)."""
	proc = bt2.run('bowtie2', datasets['unp'] + ['--sam-opt-config', config])
	assert 'not a valid SAM Optional flag' not in proc.stderr
	assert body(proc.stdout) == default_lines(bt2, datasets['unp'])


def test_sam_opt_config_case_insensitive(bt2, datasets):
	"""Tag names in --sam-opt-config are case-insensitive."""
	upper = bt2.run('bowtie2', datasets['unp'] + ['--sam-opt-config', '-MD,ZS']).stdout
	lower = bt2.run('bowtie2', datasets['unp'] + ['--sam-opt-config', '-md,zs']).stdout
	assert body(upper) == body(lower)
	assert all('MD' not in opt_dict(l) and 'ZS' in opt_dict(l) for l in body(upper))


def test_sam_opt_config_invalid_tag(bt2, datasets):
	"""An unknown tag in --sam-opt-config is reported on stderr."""
	proc = bt2.run('bowtie2', datasets['unp'] + ['--sam-opt-config', 'qq'], check=False)
	assert 'Error: qq is not a valid SAM Optional flag.' in proc.stderr


def test_sam_opt_config_xr_passthrough(bt2, datasets, unp):
	"""--sam-opt-config xr follows each record with the input record, %-escaped
	(the line the wrapper reads for --un/--al)."""
	out = bt2.run('bowtie2', datasets['unp'] + ['--sam-opt-config', 'xr', '--no-hd']).stdout
	lines = out.splitlines()
	assert len(lines) == 2 * len(unp['reads'])
	for i, (name, seq, qual) in enumerate(unp['reads']):
		assert lines[2 * i].split('\t')[0] == qname(name)
		assert unescape_passthrough(lines[2 * i + 1]) == '@%s\n%s\n+\n%s\n' % (name, seq, qual)


# --------------------------------------------------------------- header --

def test_header_lines(bt2, genome, genome_index, unp):
	"""@HD, one @SQ per reference (name, length, index order) and @PG with
	the program, version and command line; the @HD line (VN:1.5, SO, GO)
	pins the current output."""
	header, _, _ = bt2.align(unp_args(genome_index, unp))
	assert header[0] == '@HD\tVN:1.5\tSO:unsorted\tGO:query'
	assert header[1:4] == ['@SQ\tSN:%s\tLN:%d' % (n, len(s)) for n, s in genome]
	assert len(header) == 5 and header[4].startswith('@PG\t')
	pg = dict(f.split(':', 1) for f in header[4].split('\t')[1:])
	assert sorted(pg) == ['CL', 'ID', 'PN', 'VN']
	assert pg['ID'] == 'bowtie2' and pg['PN'] == 'bowtie2'
	assert pg['VN'] == version_string(bt2)
	cl = pg['CL']
	assert cl.startswith('"') and cl.endswith('"')
	assert 'bowtie2-align-s' in cl
	assert '-x %s' % genome_index in cl and unp['fq'] in cl and '-p 1' in cl


def test_sq_name_truncated_at_whitespace(bt2, rep):
	"""@SQ SN and RNAME use the FASTA name up to its first whitespace."""
	header, recs, _ = bt2.align(['-x', rep['index'], '-U', rep['fq'], '-p', '1'])
	assert [l for l in header if l.startswith('@SQ')] == ['@SQ\tSN:rep\tLN:%d' % len(rep['seq'])]
	assert set(r.rname for r in recs if not r.is_unmapped) == {'rep'}
	near = h.by_name(recs)['near'][0]
	assert near.pos == rep['near_pos']


# -------------------------------------------------------------- summary --

def test_unpaired_summary(bt2, rep):
	"""The unpaired summary counts reads aligned 0, 1 and >1 times as the
	MANUAL lays it out."""
	proc = bt2.run('bowtie2', ['-x', rep['index'], '-U', rep['fq'], '-p', '1'])
	n, n0, n1, nm = 6, 1, 3, 2   # unal; uni_fw, uni_rc, near; rep_fw, rep_rc
	expected = (
		'%d reads; of these:\n'
		'  %d (%s) were unpaired; of these:\n'
		'    %d (%s) aligned 0 times\n'
		'    %d (%s) aligned exactly 1 time\n'
		'    %d (%s) aligned >1 times\n'
		'%s overall alignment rate\n') % (n, n, pct(n, n), n0, pct(n0, n), n1, pct(n1, n),
		nm, pct(nm, n), pct(n1 + nm, n))
	assert proc.stderr == expected


def test_summary_zero_reads(bt2, genome_index, tmp_path):
	"""An empty read file gives a header-only SAM and a zero-read summary."""
	empty = tmp_path / 'empty.fq'
	empty.write_text('')
	proc = bt2.run('bowtie2', ['-x', genome_index, '-U', str(empty), '-p', '1'])
	assert body(proc.stdout) == []
	assert proc.stderr.splitlines()[-2:] == ['0 reads', '0.00% overall alignment rate']


# ---------------------------------------------------------- performance --

def test_threads_same_records(bt2, genome_index, many):
	"""-p 2, -p 4 and --threads 3 report the same records as -p 1 (as a
	multiset); -p 1 keeps input order."""
	args = ['-x', genome_index, '-U', many['fq']]
	ref = body(bt2.run('bowtie2', args + ['-p', '1']).stdout)
	assert [l.split('\t')[0] for l in ref] == many['names']
	for p in (['-p', '2'], ['-p', '4'], ['--threads', '3']):
		out = body(bt2.run('bowtie2', args + p).stdout)
		assert sorted(out) == sorted(ref), p


def test_reorder(bt2, genome_index, many):
	"""--reorder with -p 4 gives exactly the -p 1 output, in input order; with
	-p 1 it changes nothing."""
	args = ['-x', genome_index, '-U', many['fq']]
	ref = body(bt2.run('bowtie2', args + ['-p', '1']).stdout)
	assert body(bt2.run('bowtie2', args + ['-p', '4', '--reorder']).stdout) == ref
	assert body(bt2.run('bowtie2', args + ['-p', '1', '--reorder']).stdout) == ref


def test_mm(bt2, genome_index, unp):
	"""--mm (memory-mapped index) gives the same alignments."""
	ref = bt2.run('bowtie2', unp_args(genome_index, unp)).stdout
	assert body(bt2.run('bowtie2', unp_args(genome_index, unp) + ['--mm']).stdout) == body(ref)


def thread_counts(bt2, genome_index, runs, settled, timeout=30):
	"""Number of threads of bowtie2-align-s, waiting on stdin, for each
	(extra args, OMP_NUM_THREADS) in runs, read from /proc.

	A process starts its threads one after another, so a count read too
	early is too low, and a process held up on a loaded machine can sit at
	such a count for a while. The counts are therefore taken once every
	count has held for 0.5 s AND settled(counts) is true, and only if that
	never happens are the last counts returned after `timeout` s. A slow
	machine then delays the reading instead of failing the test, and a real
	failure still fails, after the timeout."""
	procs = []
	try:
		for args, omp in runs:
			env = env_without_omp() if omp is None else env_without_omp(OMP_NUM_THREADS=omp)
			procs.append(subprocess.Popen([align_bin(bt2), '--wrapper', 'basic-0', '-x',
				genome_index, '-U', '-', '-S', os.devnull] + args, stdin=subprocess.PIPE,
				stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env))
		last = [None] * len(procs)
		since = [time.time()] * len(procs)
		deadline = time.time() + timeout
		while True:
			now = time.time()
			for i, p in enumerate(procs):
				assert p.poll() is None, 'bowtie2-align-s %s exited early' % (runs[i],)
				n = len(os.listdir('/proc/%d/task' % p.pid))
				if n != last[i]:
					last[i], since[i] = n, now
			if all(now - s >= 0.5 for s in since) and settled(last):
				return last
			if now > deadline:
				return last
			time.sleep(0.01)
	finally:
		for p in procs:
			p.stdin.close()
		for p in procs:
			p.wait(timeout=60)


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='counts threads via /proc')
def test_omp_num_threads_sets_default_thread_count(bt2, genome_index):
	"""OMP_NUM_THREADS sets the number of search threads when -p is absent,
	and -p overrides it."""
	def settled(counts):
		p1, p3, omp3, omp3_p1 = counts
		return p3 - p1 == 2 and omp3 == p3 and omp3_p1 == p1

	p1, p3, omp3, omp3_p1 = thread_counts(bt2, genome_index, [
		(['-p', '1'], None), (['-p', '3'], None), ([], '3'), (['-p', '1'], '3')], settled)
	assert p3 - p1 == 2
	assert omp3 == p3
	assert omp3_p1 == p1


def test_omp_num_threads_same_records(bt2, genome_index, many):
	"""With OMP_NUM_THREADS=3 and no -p the records equal those of -p 1."""
	args = ['-x', genome_index, '-U', many['fq']]
	ref = body(bt2.run('bowtie2', args + ['-p', '1']).stdout)
	out = body(bt2.run('bowtie2', args, env=env_without_omp(OMP_NUM_THREADS='3')).stdout)
	assert sorted(out) == sorted(ref)


def test_omp_num_threads_not_a_number(bt2, genome_index, unp):
	"""A non-numeric OMP_NUM_THREADS is an error (no -p given)."""
	proc = bt2.run('bowtie2', ['-x', genome_index, '-U', unp['fq']], check=False,
		env=env_without_omp(OMP_NUM_THREADS='abc'))
	assert proc.returncode != 0
	assert 'Error' in proc.stderr
	assert body(proc.stdout) == []


# ----------------------------------------------------- --version / --help --

def test_version(bt2):
	"""--version names the aligner binary and the version, and exits 0.

	The binary is bowtie2-align-s, or bowtie2-align-s-v256 where the CPU
	dispatch selects the AVX2 build."""
	proc = bt2.run('bowtie2', ['--version'])
	first = proc.stdout.splitlines()[0]
	assert re.search(r'bowtie2-align-s(-v256)? version %s$' % re.escape(version_string(bt2)), first), first


@pytest.mark.parametrize('opt', ['-h', '--help'])
def test_help(bt2, opt):
	"""-h/--help prints the usage, listing the output and performance options,
	to stdout and exits 0."""
	proc = bt2.run('bowtie2', [opt])
	out = proc.stdout
	assert 'Usage:' in out
	assert 'bowtie2 [options]* -x <bt2-idx>' in out
	for o in ('-t/--time', '--un <path>', '--al <path>', '--quiet', '--met-file <path>',
			'--met-stderr', '--met <int>', '--no-unal', '--no-sq', '--rg-id <text>',
			'--rg <text>', '--omit-sec-seq', '--xeq', '--sam-opt-config', '-p/--threads',
			'--reorder', '--mm', '--version', '-h/--help'):
		assert o in out, o
	assert out == bt2.run('bowtie2', ['-h']).stdout


# --------------------------------------------------------------- errors --

def test_missing_index(bt2, unp, tmp_path):
	"""A -x that names no index fails with the wrapper's message."""
	prefix = str(tmp_path / 'nothing_here')
	proc = bt2.run('bowtie2', ['-x', prefix, '-U', unp['fq']], check=False)
	assert proc.returncode != 0
	assert '"%s" does not exist or is not a Bowtie 2 index' % prefix in proc.stderr
	assert proc.stdout == ''


def test_missing_index_binary(bt2, unp, tmp_path):
	"""bowtie2-align-s itself also fails on a missing index."""
	prefix = str(tmp_path / 'nothing_here')
	proc = bt2.run(align_bin(bt2), ['-x', prefix, '-U', unp['fq']], check=False)
	assert proc.returncode == 1
	assert 'Could not locate a Bowtie index corresponding to basename "%s"' % prefix in proc.stderr
	assert body(proc.stdout) == []


def test_missing_read_file(bt2, genome_index, tmp_path):
	"""A read file that does not exist is an error when it is the only one."""
	missing = str(tmp_path / 'missing.fq')
	proc = bt2.run('bowtie2', ['-x', genome_index, '-U', missing], check=False)
	assert proc.returncode == 1
	assert 'Could not open read file "%s"' % missing in proc.stderr
	assert 'No input read files were valid' in proc.stderr
	assert body(proc.stdout) == []


def test_missing_read_file_in_list_is_skipped(bt2, genome_index, unp, tmp_path):
	"""A missing file in a list is skipped with a warning and the others are
	aligned (pins the current behaviour)."""
	missing = str(tmp_path / 'missing.fq')
	proc = bt2.run('bowtie2', ['-x', genome_index, '-U', missing + ',' + unp['fq'], '-p', '1'])
	assert 'Could not open read file "%s"' % missing in proc.stderr
	assert body(proc.stdout) == body(bt2.run('bowtie2', unp_args(genome_index, unp)).stdout)


def test_unknown_option(bt2, genome_index, unp):
	"""An unknown option prints the usage and exits 1."""
	proc = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--no-such-option-xyz'], check=False)
	assert proc.returncode == 1
	assert 'Usage:' in proc.stderr
	assert 'bowtie2-align exited with value 1' in proc.stderr
	assert proc.stdout == ''


BAD_VALUES = [
	(['-k', '0'], '-k arg must be at least 1'),
	(['-k', '-1'], '-k arg must be at least 1'),
	(['-p', '0'], '-p/--threads arg must be at least 1'),
	(['--threads', '0'], '-p/--threads arg must be at least 1'),
	(['--met', '0'], '--metrics arg must be at least 1'),
	(['-u', '0'], '-u/--qupto arg must be at least 1'),
	(['-s', '-1'], '-s arg must be positive'),
	(['-3', '-1'], '-3/--trim3 arg must be at least 0'),
	(['-5', '-1'], '-5/--trim5 arg must be at least 0'),
	(['--trim-to', '4:10'], '--trim-to: trim position must be either 3 or 5'),
	(['-L', '0'], '-L argument must be within the interval [1,32]; was 0'),
	(['-L', '33'], '-L argument must be within the interval [1,32]; was 33'),
	(['-N', '2'], '-N argument must be within the interval [0,1]; was 2'),
	(['-i', 'X,1,2'], "Bad function type 'X'"),
	(['--score-min', 'X,1,2'], "Bad function type 'X'"),
	(['--score-min', 'foo'], "Bad function type 'foo'"),
	(['--local', '--score-min', 'L,-1,0'], 'the match penalty is greater than 0 (2) but the '
		'--score-min function can be less than or equal to zero'),
	(['--mp', '2,6'], 'Maximum mismatch penalty (2) is less than minimum penalty (6'),
	(['--gbar', '0'], '--gbar must be no less than 1'),
	(['--dpad', '-1'], '--dpad must be no less than 0'),
	(['--seed', '-1'], '--seed arg must be at least 0'),
	(['-I', '-1'], '-I arg must be positive'),
	(['-X', '0'], '-X arg must be at least 1'),
	(['--soft-clipped-unmapped-tlen'], '--soft-clipped-unmapped-tlen can only be set for local alignments'),
	(['--preserve-tags'], 'can only be used when aligning BAM reads'),
	(['--align-paired-reads'], '--align-paired-reads can only be used when aligning BAM reads'),
]


@pytest.mark.parametrize('args,message', BAD_VALUES, ids=[' '.join(a) for a, _ in BAD_VALUES])
def test_bad_option_value(bt2, genome_index, unp, args, message):
	"""Invalid option values exit 1 with bowtie2's message and no records."""
	proc = bt2.run('bowtie2', unp_args(genome_index, unp) + args, check=False)
	assert proc.returncode == 1
	assert message in proc.stderr
	assert body(proc.stdout) == []


def test_no_index_option(bt2, unp):
	"""Leaving out -x is an error."""
	proc = bt2.run('bowtie2', ['-U', unp['fq']], check=False)
	assert proc.returncode == 1
	assert 'No index, query, or output file specified!' in proc.stderr


def test_no_reads_option(bt2, genome_index):
	"""Leaving out every read input is an error."""
	proc = bt2.run('bowtie2', ['-x', genome_index], check=False)
	assert proc.returncode == 1
	assert 'Must specify at least one read input with -U/-1/-2' in proc.stderr


# -------------------------------------------------------------- wrapper --

def test_debug_binary(bt2, genome_index, unp):
	"""--debug runs bowtie2-align-s-debug (which announces debug mode) and
	gives the same alignments."""
	if not bt2.has_binary('bowtie2-align-s-debug'):
		pytest.skip('bowtie2-align-s-debug not built')
	ref = bt2.run('bowtie2', unp_args(genome_index, unp))
	proc = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--debug'])
	assert 'Running in debug mode' in proc.stderr
	assert 'Running in debug mode' not in ref.stderr
	assert body(proc.stdout) == body(ref.stdout)


def test_sanitized_binary(bt2, genome_index, unp, tmp_path):
	"""--sanitized runs bowtie2-align-s-sanitized and gives the same alignments."""
	if not bt2.has_binary('bowtie2-align-s-sanitized'):
		pytest.skip('bowtie2-align-s-sanitized not built')
	ref = bt2.run('bowtie2', unp_args(genome_index, unp))
	out = str(tmp_path / 'o.sam')
	proc = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--sanitized', '--verbose', '-S', out])
	assert 'bowtie2-align-s-sanitized"' in proc.stderr
	with open(out) as f:
		assert body(f.read()) == body(ref.stdout)


def test_verbose(bt2, genome_index, unp, tmp_path):
	"""--verbose makes the wrapper print (INFO): lines, including the binary
	command, without changing the alignments written with -S."""
	ref = bt2.run('bowtie2', unp_args(genome_index, unp))
	out = str(tmp_path / 'o.sam')
	proc = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--verbose', '-S', out])
	info = [l for l in proc.stderr.splitlines() if l.startswith('(INFO):')]
	assert any('bowtie2-align-s" --wrapper basic-0' in l for l in info), proc.stderr
	assert not any(l.startswith('(INFO):') for l in ref.stderr.splitlines())
	with open(out) as f:
		assert body(f.read()) == body(ref.stdout)


def bz2_reads(unp, tmp_path):
	"""unp's reads as .bz2, which the wrapper decompresses with the bzip2
	tool into a named pipe or a temporary file; skips without bzip2."""
	if shutil.which('bzip2') is None:
		pytest.skip('bzip2 not installed')
	path = str(tmp_path / 'reads.fq.bz2')
	h.write_fastq(path, unp['reads'], compress='bz2')
	return path


@pytest.mark.skipif(not hasattr(os, 'mkfifo'), reason='needs named pipes')
def test_keep_named_pipe_in_temp_directory(bt2, genome_index, unp, tmp_path):
	"""For .bz2 input the wrapper decompresses through a named pipe in
	--temp-directory; --keep leaves it there."""
	tmp = tmp_path / 'tmp'
	tmp.mkdir()
	fq = bz2_reads(unp, tmp_path)
	proc = bt2.run('bowtie2', ['-x', genome_index, '-U', fq, '-p', '1', '--temp-directory',
		str(tmp), '--keep'])
	assert body(proc.stdout) == body(bt2.run('bowtie2', unp_args(genome_index, unp)).stdout)
	left = os.listdir(str(tmp))
	assert len(left) == 1 and left[0].endswith('.unps'), left
	assert stat.S_ISFIFO(os.stat(str(tmp / left[0])).st_mode)


def test_temp_directory_cleaned_without_keep(bt2, genome_index, unp, tmp_path):
	"""Without --keep the wrapper removes its temporary files."""
	tmp = tmp_path / 'tmp'
	tmp.mkdir()
	fq = bz2_reads(unp, tmp_path)
	for extra in ([], ['--no-named-pipes']):
		proc = bt2.run('bowtie2', ['-x', genome_index, '-U', fq, '-p', '1', '--temp-directory',
			str(tmp)] + extra)
		assert body(proc.stdout) == body(bt2.run('bowtie2', unp_args(genome_index, unp)).stdout)
		assert os.listdir(str(tmp)) == []


def test_no_named_pipes_keep(bt2, genome_index, unp, tmp_path):
	"""--no-named-pipes decompresses .bz2 input to a regular file, which
	--keep leaves holding the decompressed reads."""
	tmp = tmp_path / 'tmp'
	tmp.mkdir()
	fq = bz2_reads(unp, tmp_path)
	proc = bt2.run('bowtie2', ['-x', genome_index, '-U', fq, '-p', '1', '--temp-directory',
		str(tmp), '--keep', '--no-named-pipes'])
	assert body(proc.stdout) == body(bt2.run('bowtie2', unp_args(genome_index, unp)).stdout)
	left = os.listdir(str(tmp))
	assert len(left) == 1 and left[0].endswith('.unps'), left
	path = str(tmp / left[0])
	assert stat.S_ISREG(os.stat(path).st_mode)
	with open(path) as f:
		assert f.read() == h.fastq_text(unp['reads'])


def test_log_file(bt2, genome_index, unp, tmp_path):
	"""--log-file sends what would go to stderr (the summary) to the file."""
	log = str(tmp_path / 'log.txt')
	ref = bt2.run('bowtie2', unp_args(genome_index, unp))
	proc = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--log-file', log])
	assert proc.stderr == ''
	with open(log) as f:
		assert f.read() == ref.stderr
	assert 'overall alignment rate' in ref.stderr
	assert body(proc.stdout) == body(ref.stdout)


@pytest.fixture(scope='module')
def dual_index(bt2, genome_index_large, tmp_path_factory):
	"""A prefix holding a small index of a 1500-bp 'other' reference and the
	large index of `genome`, so the @SQ lines show which one was used."""
	d = tmp_path_factory.mktemp('dual')
	src_dir, src_base = os.path.split(genome_index_large)
	for name in os.listdir(src_dir):
		if name.startswith(src_base + '.') and name.endswith('.bt2l'):
			shutil.copy(os.path.join(src_dir, name), str(d / ('dual' + name[len(src_base):])))
	other = h.write_fasta(str(d / 'other.fa'), [('other', h.random_dna(1500, random.Random(3)))])
	prefix = str(d / 'dual')
	bt2.build_index(other, prefix, extra=['--seed', '0'])
	assert os.path.exists(prefix + '.1.bt2') and os.path.exists(prefix + '.1.bt2l')
	return prefix


def sq_names(text):
	return [l.split('\t')[1][3:] for l in headers(text) if l.startswith('@SQ')]


def test_large_index_forces_bt2l(bt2, genome, dual_index, unp):
	"""With .bt2 and .bt2l at one prefix the wrapper uses the .bt2 unless
	--large-index is given, which selects the .bt2l and bowtie2-align-l."""
	args = ['-x', dual_index, '-U', unp['fq'], '-p', '1']
	small = bt2.run('bowtie2', args)
	assert sq_names(small.stdout) == ['other']
	large = bt2.run('bowtie2', args + ['--large-index'])
	assert sq_names(large.stdout) == [n for n, _ in genome]
	_, recs, _ = bt2.align(args + ['--large-index'])
	for r in recs:
		exp = unp['aligned'].get(r.qname)
		if exp:
			assert (r.rname, r.pos, r.is_reverse) == (exp['rname'], exp['pos'], exp['rev'])


def test_index_autodetect_bt2l(bt2, genome_index, genome_index_large, unp, tmp_path):
	"""A prefix with only a .bt2l index is used through bowtie2-align-l
	without --large-index, with the same alignments as the .bt2 index."""
	out = str(tmp_path / 'o.sam')
	proc = bt2.run('bowtie2', ['-x', genome_index_large, '-U', unp['fq'], '-p', '1', '--verbose',
		'-S', out])
	assert 'bowtie2-align-l" --wrapper basic-0' in proc.stderr
	with open(out) as f:
		assert body(f.read()) == body(bt2.run('bowtie2', unp_args(genome_index, unp)).stdout)
	small = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--verbose', '-S', out])
	assert 'bowtie2-align-s" --wrapper basic-0' in small.stderr


def test_large_index_missing_bt2l(bt2, genome_index, unp):
	"""--large-index with only a .bt2 index fails."""
	proc = bt2.run('bowtie2', unp_args(genome_index, unp) + ['--large-index'], check=False)
	assert proc.returncode != 0
	assert 'Cannot find the large index %s.1.bt2l' % genome_index in proc.stderr


def test_bowtie2_indexes_env(bt2, genome_index, unp, tmp_path):
	"""-x <name> is looked up in $BOWTIE2_INDEXES when it is not in the
	current directory."""
	d, name = os.path.split(genome_index)
	env = env_without_omp(BOWTIE2_INDEXES=d)
	proc = bt2.run('bowtie2', ['-x', name, '-U', unp['fq'], '-p', '1'], cwd=str(tmp_path), env=env)
	assert body(proc.stdout) == body(bt2.run('bowtie2', unp_args(genome_index, unp)).stdout)


def test_ref_string(bt2, genome_dict, tmp_path):
	"""--ref-string aligns against a reference given on the command line
	(named '1'); the wrapper needs bowtie2-build on PATH for it."""
	c1 = genome_dict['chr1']
	ref = c1[1000:1600]
	fq = h.write_fastq(str(tmp_path / 'r.fq'), [('a', ref[100:150]), ('b', h.revcomp(ref[300:350]))])
	tmp = tmp_path / 'tmp'
	tmp.mkdir()
	env = env_without_omp(PATH=bt2.bin_dir + os.pathsep + os.environ.get('PATH', ''))
	out = str(tmp_path / 'o.sam')
	bt2.run('bowtie2', ['--ref-string', ref, '-U', fq, '-p', '1', '-S', out,
		'--temp-directory', str(tmp)], env=env)
	with open(out) as f:
		header, recs = h.parse_sam(f.read())
	assert [l for l in header if l.startswith('@SQ')] == ['@SQ\tSN:1\tLN:600']
	assert [(r.qname, r.rname, r.pos, r.is_reverse) for r in recs] == \
		[('a', '1', 101, False), ('b', '1', 301, True)]
	assert os.listdir(str(tmp)) == []


def test_bam_output(bt2, genome_index, unp, tmp_path):
	"""--bam writes BAM (through samtools) holding the same records."""
	if shutil.which('samtools') is None:
		pytest.skip('samtools not installed')
	out = str(tmp_path / 'o.bam')
	bt2.run('bowtie2', unp_args(genome_index, unp) + ['--bam', '-S', out])
	with open(out, 'rb') as f:
		assert f.read(2) == b'\x1f\x8b'
	view = subprocess.run(['samtools', 'view', out], capture_output=True, text=True, check=True)
	assert body(view.stdout) == body(bt2.run('bowtie2', unp_args(genome_index, unp)).stdout)
