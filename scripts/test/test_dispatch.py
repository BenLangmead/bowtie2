"""The x86-64 dispatch between the 128-bit and the -v256 (AVX2) aligners.

On x86-64, `make` builds bowtie2-align-s with SSE2 and bowtie2-align-s-v256
with -march=x86-64-v3. At startup the 128-bit binary checks CPUID and, if
the CPU has every x86-64-v3 feature, executes the -v256 binary of the same
name in its own directory, passing on its own argv[0]. These tests check
that it does so exactly when the CPU qualifies, that --version and @PG name
bowtie2-align-s either way, and that both binaries write the same SAM
records.

They run only on x86-64 Linux with the -v256 binary built, and skip
elsewhere.
"""

import os
import platform
import random
import shutil
import subprocess

import pytest

import bt2_helpers as h


# x86-64-v3 features as Linux names them in /proc/cpuinfo. Linux reports
# SSE3 as pni and LZCNT as abm, and it leaves out avx (and the features that
# need it) when the OS does not save the YMM registers.
X86_64_V3_FLAGS = {'pni', 'ssse3', 'fma', 'cx16', 'sse4_1', 'sse4_2', 'movbe',
	'popcnt', 'avx', 'f16c', 'avx2', 'bmi1', 'bmi2', 'abm'}


def cpu_flags():
	"""The CPU feature flags from /proc/cpuinfo, or None if unavailable."""
	try:
		with open('/proc/cpuinfo') as f:
			for line in f:
				if line.startswith('flags'):
					return set(line.split(':', 1)[1].split())
	except OSError:
		pass
	return None


@pytest.fixture(scope='module')
def v256(bt2):
	"""Path of bowtie2-align-s-v256 and whether this CPU can run it. Skips
	unless the platform is x86-64 Linux and the -v256 binary was built."""
	if platform.machine().lower() not in ('x86_64', 'amd64'):
		pytest.skip('the -v256 build exists only on x86-64')
	flags = cpu_flags()
	if flags is None:
		pytest.skip('/proc/cpuinfo not available')
	path = os.path.join(bt2.bin_dir, 'bowtie2-align-s-v256')
	if not os.access(path, os.X_OK):
		pytest.skip('bowtie2-align-s-v256 was not built')
	return path, X86_64_V3_FLAGS <= flags


@pytest.fixture(scope='module')
def align_128(bt2, v256, tmp_path_factory):
	"""A copy of bowtie2-align-s in a directory without the -v256 binary,
	so that it cannot dispatch and always runs the 128-bit code."""
	d = tmp_path_factory.mktemp('align128')
	path = str(d / 'bowtie2-align-s')
	shutil.copy(os.path.join(bt2.bin_dir, 'bowtie2-align-s'), path)
	return path


def version_lines(exe):
	"""The first line of `exe --version` and its Options: line."""
	proc = subprocess.run([exe, '--version'], capture_output=True, text=True)
	assert proc.returncode == 0, proc.stderr
	lines = proc.stdout.splitlines()
	options = [l for l in lines if l.startswith('Options:')]
	assert len(options) == 1, proc.stdout
	return lines[0], options[0]


def test_dispatch_follows_cpu_features(bt2, v256):
	"""bowtie2-align-s runs the -v256 build if and only if the CPU has every
	x86-64-v3 feature. The Options: line of --version shows which build ran,
	and the first line names bowtie2-align-s either way."""
	_, capable = v256
	first, options = version_lines(os.path.join(bt2.bin_dir, 'bowtie2-align-s'))
	assert ('-march=x86-64-v3' in options) == capable, options
	assert 'bowtie2-align-s version' in first and '-v256' not in first, first


def test_pg_names_program_run(bt2, v256, lambda_index, example_reads):
	"""The @PG line names bowtie2-align-s, not the -v256 binary, when the
	dispatch runs the -v256 build."""
	proc = subprocess.run([os.path.join(bt2.bin_dir, 'bowtie2-align-s'),
		'-x', lambda_index, '-U', example_reads['reads_1'], '-u', '10'],
		capture_output=True, text=True)
	assert proc.returncode == 0, proc.stderr
	pg = [l for l in proc.stdout.splitlines() if l.startswith('@PG')]
	assert len(pg) == 1, pg
	assert 'bowtie2-align-s ' in pg[0] and '-v256' not in pg[0], pg[0]


def test_no_dispatch_without_v256_binary(align_128):
	"""Without a -v256 binary beside it, bowtie2-align-s warns if the CPU
	qualifies and carries on with the 128-bit code."""
	first, options = version_lines(align_128)
	assert 'bowtie2-align-s version' in first, first
	assert '-march=x86-64-v3' not in options, options


def boundary_reads(fasta, n, length, seed):
	"""n reads of `length` cut from random places in the reference, each with
	up to four mismatches at random positions and random base qualities.
	With bowtie2's quality-aware mismatch penalty (2 to 6), their alignment
	scores land on and around the --score-min values tested below."""
	ref = ''.join(line.strip() for line in open(fasta) if not line.startswith('>'))
	rng = random.Random(seed)
	reads = []
	for i in range(n):
		pos = rng.randrange(len(ref) - length)
		seq = ref[pos:pos + length]
		seq = h.mutate(seq, rng.sample(range(length), rng.randint(0, 4)))
		if rng.random() < 0.5:
			seq = h.revcomp(seq)
		qual = ''.join(chr(33 + rng.randint(2, 40)) for _ in range(length))
		reads.append(('r%d' % i, seq, qual))
	return reads


def sam_without_pg(text):
	return [line for line in text.splitlines() if not line.startswith('@PG')]


def assert_same_sam(exes, args):
	"""Run each aligner binary with args and check that their SAM output,
	apart from @PG, is the same. On failure, show the first record that
	differs."""
	out = []
	for exe in exes:
		proc = subprocess.run([exe] + args, capture_output=True, text=True)
		assert proc.returncode == 0, proc.stderr
		out.append(sam_without_pg(proc.stdout))
	for a, b in zip(out[0], out[1]):
		assert a == b, '%s:\n%s\n%s:\n%s' % (exes[0], a, exes[1], b)
	assert len(out[0]) == len(out[1])


# Read lengths and --score-min functions for which C + L*x lands within one
# rounding of an integer, so that a fused multiply-add truncates to a
# different minimum score than a separate multiply and add. The built-in
# defaults are single-precision constants and hit no such length, but values
# given on the command line are parsed as doubles: L,-0.6,-0.6 does at x = 24,
# L,-0.1,-0.7 at 17, L,-0.7,-0.022 at 150 and, in local mode, L,0.6,0.6 at 24.
# With -ffp-contract=fast each of these cases changes some reads' output.
FP_CASES = [
	('L,-0.6,-0.6-24', 24, ['--score-min', 'L,-0.6,-0.6', '-L', '10']),
	('L,-0.1,-0.7-17', 17, ['--score-min', 'L,-0.1,-0.7', '-L', '10']),
	('L,-0.7,-0.022-150', 150, ['--score-min', 'L,-0.7,-0.022']),
	('local-L,0.6,0.6-24', 24, ['--local', '--score-min', 'L,0.6,0.6', '--ma', '1', '-L', '10']),
]


@pytest.mark.parametrize('name,length,opts', FP_CASES, ids=[c[0] for c in FP_CASES])
def test_v256_matches_128bit_at_score_boundaries(v256, align_128, lambda_fasta,
		lambda_index, tmp_path, name, length, opts):
	"""The -v256 and 128-bit builds write identical SAM records for reads
	whose minimum score depends on how C + L*x is rounded. The -v256 build
	is compiled with -ffp-contract=off so that it rounds the same way."""
	path, capable = v256
	if not capable:
		pytest.skip('this CPU cannot run the -v256 build')
	reads = str(h.write_fastq(str(tmp_path / 'reads.fq'),
		boundary_reads(lambda_fasta, 2000, length, seed=length)))
	assert_same_sam((path, align_128),
		['-x', lambda_index, '-U', reads, '-i', 'C,1'] + opts)


EXAMPLE_CASES = [
	('unpaired', ['-U', 'reads_1']),
	('paired', ['-1', 'reads_1', '-2', 'reads_2']),
	('local-paired', ['--local', '-1', 'reads_1', '-2', 'reads_2']),
	('longreads-local', ['--local', '-U', 'longreads']),
	('k5-paired', ['-k', '5', '-1', 'reads_1', '-2', 'reads_2']),
]


@pytest.mark.parametrize('name,args', EXAMPLE_CASES, ids=[c[0] for c in EXAMPLE_CASES])
def test_v256_matches_128bit_on_example(v256, align_128, lambda_index,
		example_reads, name, args):
	"""The -v256 and 128-bit builds write identical SAM records for the
	lambda example reads."""
	path, capable = v256
	if not capable:
		pytest.skip('this CPU cannot run the -v256 build')
	assert_same_sam((path, align_128),
		['-x', lambda_index] + [example_reads.get(a, a) for a in args])
