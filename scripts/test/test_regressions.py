"""Regression tests for specific bowtie2 issues.

Each test names the GitHub issue (or behaviour) it guards and uses the
lambda phage example shipped in example/. New regression tests belong here
when they do not fit naturally in one of the topic modules.

Run from the repository root after `make`:

	python3 -m pytest scripts/test/test_regressions.py
"""

import os
import random
import shutil
import subprocess

import pytest

import bt2_helpers as h


def sam_records(path):
	with open(path) as f:
		return h.parse_sam(f.read())[1]


def read_names(path):
	return [r[0] for r in h.read_fastq(path)]


def file_bytes(path):
	with open(path, 'rb') as f:
		return f.read()


def test_issue_288_un_with_no_unal(bt2, lambda_index, example_reads, tmp_path):
	"""--un writes the same reads with or without --no-unal (issue 288).

	With -a, every read that fails to align is one unaligned SAM record. Those
	reads are exactly the ones --un writes, and --no-unal removes exactly
	those records from the SAM and nothing else."""
	reads = example_reads['longreads']
	full_sam, full_un = str(tmp_path / 'full.sam'), str(tmp_path / 'full_un.fq')
	nounal_sam, nounal_un = str(tmp_path / 'nounal.sam'), str(tmp_path / 'nounal_un.fq')
	base = ['--quiet', '-x', lambda_index, '-U', reads, '-a']
	bt2.run('bowtie2', base + ['--un', full_un, '-S', full_sam])
	bt2.run('bowtie2', base + ['--un', nounal_un, '--no-unal', '-S', nounal_sam])
	full = sam_records(full_sam)
	unaligned = [r.qname for r in full if r.is_unmapped]
	assert unaligned, 'expected some unaligned long reads'
	assert [n.split()[0] for n in read_names(full_un)] == unaligned
	assert file_bytes(nounal_un) == file_bytes(full_un)
	assert [r.line for r in sam_records(nounal_sam)] == \
		[r.line for r in full if not r.is_unmapped]


def test_issue_279_un_conc_with_no_unal(bt2, lambda_index, example_reads, tmp_path):
	"""--un-conc and --al-conc write the same pairs with or without --no-unal (issue 279).

	Both mate files of each output hold the same number of reads, and
	--al-conc holds exactly the pairs the SAM reports as concordant."""
	outs = {}
	for name, extra in (('plain', []), ('nounal', ['--no-unal'])):
		d = tmp_path / name
		d.mkdir()
		sam = str(d / 'out.sam')
		bt2.run('bowtie2', ['--quiet', '-x', lambda_index, '-1', example_reads['reads_1'],
			'-2', example_reads['reads_2'], '--dovetail', '--un-conc', str(d / 'nconc.fq'),
			'--al-conc', str(d / 'conc.fq'), '-S', sam] + extra)
		outs[name] = d
		for kind in ('conc', 'nconc'):
			m1 = read_names(str(d / ('%s.1.fq' % kind)))
			m2 = read_names(str(d / ('%s.2.fq' % kind)))
			assert len(m1) == len(m2), (name, kind)
		if name == 'plain':
			concordant = sorted(set(r.qname for r in sam_records(sam) if r.tags.get('YT') == 'CP'))
			conc = sorted(set(n.split()[0].rsplit('/', 1)[0] for n in read_names(str(d / 'conc.1.fq'))))
			assert conc == concordant
			assert concordant, 'expected some concordant pairs'
	for f in ('conc.1.fq', 'conc.2.fq', 'nconc.1.fq', 'nconc.2.fq'):
		assert file_bytes(str(outs['nounal'] / f)) == file_bytes(str(outs['plain'] / f)), f


def test_issue_313_dots_in_output_paths(bt2, lambda_index, example_reads, tmp_path):
	"""--un-conc, --al-conc and --un work when the output path has a dot in a
	directory name, and when given a bare directory (issue 313).

	A bare directory gets the default names un-conc-mate.1/.2,
	al-conc-mate.1/.2 and un-seqs, and the outputs do not overwrite each
	other."""
	pairs = ['--quiet', '-x', lambda_index, '-1', example_reads['reads_1'],
		'-2', example_reads['reads_2'], '-S', os.devnull]
	plain_dir, dot_dir = tmp_path / 'nodotdir', tmp_path / 'v.12_.test'
	plain_dir.mkdir()
	dot_dir.mkdir()

	bt2.run('bowtie2', pairs + ['--un-conc', str(plain_dir / 'bname')])
	expected = file_bytes(str(plain_dir / 'bname.1'))
	assert read_names(str(plain_dir / 'bname.1')), 'expected some non-concordant pairs'

	bt2.run('bowtie2', pairs + ['--un-conc', str(dot_dir / 'bname')])
	assert file_bytes(str(dot_dir / 'bname.1')) == expected

	bt2.run('bowtie2', pairs + ['--un-conc', str(dot_dir)])
	assert file_bytes(str(dot_dir / 'un-conc-mate.1')) == expected

	bt2.run('bowtie2', pairs + ['--un-conc', str(dot_dir), '--al-conc', str(dot_dir),
		'--un', str(dot_dir)])
	for name in ('un-conc-mate.1', 'un-conc-mate.2', 'al-conc-mate.1', 'al-conc-mate.2',
			'un-seqs'):
		assert (dot_dir / name).is_file(), name
	assert file_bytes(str(dot_dir / 'un-conc-mate.1')) == expected


def test_index_requires_x(bt2, lambda_index, example_reads, tmp_path):
	"""An index given as a positional argument instead of with -x is an error,
	and the same command with -x succeeds."""
	sam = str(tmp_path / 'out.sam')
	reads = ['-U', example_reads['longreads'], '-S', sam]
	proc = bt2.run('bowtie2', [lambda_index] + reads, check=False)
	assert proc.returncode != 0
	bt2.run('bowtie2', ['-x', lambda_index] + reads)
	assert sam_records(sam)


def test_one_n_in_repeated_read(bt2, lambda_fasta, tmp_path):
	"""A read with one N aligns, end to end, to every copy of a repeat.

	Before seed search, bowtie2 looks for end-to-end hits with at most one
	mismatch, extending from one half of the read into the other. When the
	N lies in the half being extended and the range so far covers more than
	one place in the genome, that search must not index its per-base range
	arrays with the N (index 4 of a 4-element array). A reference made of two
	identical copies of the lambda genome puts every read in that situation.
	The sanitized aligner, when it was built, must run without reports."""
	with open(lambda_fasta) as f:
		lam = ''.join(l.strip() for l in f if not l.startswith('>'))
	fasta = h.write_fasta(str(tmp_path / 'lambda2.fa'), [('copyA', lam), ('copyB', lam)])
	prefix = str(tmp_path / 'lambda2')
	bt2.build_index(fasta, prefix, extra=['--seed', '0'])
	reads, expected = [], {}
	for off, npos in ((1000, 10), (1000, 49), (20000, 30), (30000, 5)):
		r = lam[off:off + 60]
		r = r[:npos] + 'N' + r[npos + 1:]
		for name, seq in (('fw%d_%d' % (off, npos), r), ('rc%d_%d' % (off, npos), h.revcomp(r))):
			reads.append((name, seq, 'I' * 60))
			expected[name] = (off + 1, '%d%s%d' % (npos, lam[off + npos], 59 - npos))
	fq = h.write_fastq(str(tmp_path / 'reads.fq'), reads)
	args = ['-x', prefix, '-U', fq]
	_, recs, _ = bt2.align(['--quiet', '-k', '2'] + args)
	groups = h.by_name(recs)
	assert sorted(groups) == sorted(expected)
	for name, rs in groups.items():
		pos, md = expected[name]
		assert sorted(r.rname for r in rs) == ['copyA', 'copyB'], name
		for r in rs:
			assert (r.pos, r.cigar, r.tags['NM'], r.tags['MD']) == (pos, '60M', 1, md), r.line
			assert r.is_reverse == name.startswith('rc'), r.line
	san = os.path.join(h.BIN_DIR, 'bowtie2-align-s-sanitized')
	if os.access(san, os.X_OK):
		env = dict(os.environ, ASAN_OPTIONS='detect_leaks=0',
			UBSAN_OPTIONS='halt_on_error=1')
		proc = bt2.run(san, args + ['-k', '2'], env=env, check=False)
		assert proc.returncode == 0, proc.stderr
		assert 'runtime error' not in proc.stderr and 'Sanitizer' not in proc.stderr, \
			proc.stderr


def shm_segments():
	"""(shmid, creator pid) of each SysV shared-memory segment, read from
	`ipcs -m -p`, whose columns differ between macOS and Linux."""
	cols, segs = None, []
	out = subprocess.run(['ipcs', '-m', '-p'], capture_output=True, text=True).stdout
	for line in out.splitlines():
		toks = line.split()
		low = [t.lower() for t in toks]
		if 'cpid' in low:
			cols = low
		elif cols and len(toks) == len(cols) and toks[cols.index('cpid')].isdigit():
			shmid = toks[cols.index('shmid' if 'shmid' in cols else 'id')]
			segs.append((shmid, int(toks[cols.index('cpid')])))
	return segs


def shm_modes():
	"""shmid -> permission bits of each SysV shared-memory segment, read from
	`ipcs -m` (macOS prints MODE as --rw-r--r--, Linux prints perms as 644)."""
	cols, modes = None, {}
	out = subprocess.run(['ipcs', '-m'], capture_output=True, text=True).stdout
	for line in out.splitlines():
		toks = line.split()
		low = [t.lower() for t in toks]
		if 'mode' in low or 'perms' in low:
			cols = low
		elif cols and len(toks) >= len(cols) - 1:  # Linux may leave status empty
			mode = toks[cols.index('mode' if 'mode' in cols else 'perms')]
			shmid = toks[cols.index('shmid' if 'shmid' in cols else 'id')]
			modes[shmid] = int(mode, 8) if mode.isdigit() else \
				sum(1 << (8 - i) for i, c in enumerate(mode[-9:]) if c != '-')
	return modes


@pytest.fixture
def shmem_align():
	"""Runs bowtie2-align-s --shmem directly (so the pid that creates the
	segments is known) and removes every segment those runs created.
	Skips unless the binary was built with BOWTIE_SHARED_MEM."""
	exe = os.path.join(h.BIN_DIR, 'bowtie2-align-s')
	version = subprocess.run([exe, '--version'], capture_output=True, text=True).stdout
	if '-DBOWTIE_SHARED_MEM' not in version:
		pytest.skip('bowtie2-align-s was built without BOWTIE_SHARED_MEM')
	pids = []

	def run(args, cwd, timeout=60):
		proc = subprocess.Popen([exe, '--shmem'] + [str(a) for a in args], cwd=cwd,
			stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		pids.append(proc.pid)
		try:
			out, err = proc.communicate(timeout=timeout)
		except subprocess.TimeoutExpired:
			proc.kill()
			proc.communicate()
			pytest.fail('bowtie2-align-s --shmem hung for %d s' % timeout)
		return subprocess.CompletedProcess(proc.args, proc.returncode, out, err)

	run.pids = pids
	yield run
	for shmid, cpid in shm_segments():
		if cpid in pids:
			subprocess.run(['ipcrm', '-m', shmid])


def test_shmem_does_not_share_another_index(bt2, shmem_align, tmp_path):
	"""--shmem must not attach to the segments of another index that was
	opened under the same relative name (-x idx) from another directory.
	Both genomes have the same length, so the segments have the same sizes."""
	rng = random.Random(26)
	refs = [h.random_dna(3000, rng), h.random_dna(3000, rng)]
	dirs = []
	for i, seq in enumerate(refs):
		d = tmp_path / ('d%d' % i)
		d.mkdir()
		bt2.build_index(h.write_fasta(str(d / 'ref.fa'), [('ref', seq)]),
			str(d / 'idx'), extra=['--seed', '0'])
		dirs.append(str(d))
	starts = range(0, 2900, 100)
	reads = h.write_fastq(str(tmp_path / 'reads.fq'),
		[('r%d' % p, refs[1][p:p + 50]) for p in starts])
	first = shmem_align(['-x', 'idx', '-U', reads], cwd=dirs[0])
	assert first.returncode == 0, first.stderr
	second = shmem_align(['-x', 'idx', '-U', reads], cwd=dirs[1])
	assert second.returncode == 0, second.stderr
	records = h.parse_sam(second.stdout)[1]
	assert [(r.qname, r.pos) for r in records] == [('r%d' % p, p + 1) for p in starts]


def test_shmem_failed_load_does_not_hang(bt2, shmem_align, lambda_index,
		example_reads, tmp_path):
	"""A --shmem run that dies while loading the index must not leave a
	segment that makes every later run wait for it forever."""
	for ext in ('1', '2', '3', '4', 'rev.1', 'rev.2'):
		shutil.copy('%s.%s.bt2' % (lambda_index, ext), str(tmp_path))
	ref = str(tmp_path / 'lambda_virus.4.bt2')
	with open(ref, 'r+b') as f:
		f.truncate(os.path.getsize(ref) // 2)
	for _ in range(2):
		proc = shmem_align(['-x', 'lambda_virus', '-U', example_reads['reads_1']],
			cwd=str(tmp_path), timeout=20)
		assert proc.returncode != 0
		assert 'Only read' in proc.stderr


def test_shmem_rejected_without_support(bt2, lambda_index, example_reads):
	"""A build without BOWTIE_SHARED_MEM rejects --shmem instead of crashing."""
	version = bt2.run(os.path.join(h.BIN_DIR, 'bowtie2-align-s'), ['--version']).stdout
	if '-DBOWTIE_SHARED_MEM' in version:
		pytest.skip('bowtie2-align-s was built with BOWTIE_SHARED_MEM')
	proc = bt2.run('bowtie2', ['--shmem', '-x', lambda_index, '-U',
		example_reads['reads_1']], check=False)
	assert proc.returncode == 1
	assert 'BOWTIE_SHARED_MEM' in proc.stderr


def test_shmem_segments_not_world_writable(bt2, shmem_align, lambda_index,
		example_reads, tmp_path):
	"""--shmem segments are writable only by their owner (they were 0666),
	and a second run, which attaches to them read-only, aligns the same."""
	args = ['-x', lambda_index, '-U', example_reads['reads_1']]
	first = shmem_align(args, cwd=str(tmp_path))
	assert first.returncode == 0, first.stderr
	second = shmem_align(args, cwd=str(tmp_path))
	assert second.returncode == 0, second.stderr
	assert second.stdout == first.stdout
	segs = shm_segments()
	assert not [shmid for shmid, cpid in segs if cpid == shmem_align.pids[1]]
	created = [shmid for shmid, cpid in segs if cpid == shmem_align.pids[0]]
	assert created
	modes = shm_modes()
	assert [oct(modes[shmid]) for shmid in created] == [oct(0o644)] * len(created)
