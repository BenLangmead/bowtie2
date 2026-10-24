"""Regression tests for specific bowtie2 issues.

Each test names the GitHub issue (or behaviour) it guards and uses the
lambda phage example shipped in example/. New regression tests belong here
when they do not fit naturally in one of the topic modules.

Run from the repository root after `make`:

	python3 -m pytest scripts/test/test_regressions.py
"""

import os

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
