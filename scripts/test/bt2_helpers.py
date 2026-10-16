"""Helpers shared by the pytest suite in scripts/test (test_*.py).

Everything here is plain Python so that expected values in the tests can be
derived independently of bowtie2: references are generated from a fixed
seed, reads are cut from them at known positions, and SAM output is parsed
into records whose fields the tests compare against those positions.
"""

import gzip
import bz2
import os
import random
import shutil
import subprocess


ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN_DIR = os.path.abspath(os.environ.get('BOWTIE2_BIN_DIR', ROOT))
EXAMPLE_DIR = os.path.join(ROOT, 'example')

COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

# SAM FLAG bits
PAIRED = 0x1
PROPER_PAIR = 0x2
UNMAPPED = 0x4
MATE_UNMAPPED = 0x8
REVERSE = 0x10
MATE_REVERSE = 0x20
FIRST_IN_PAIR = 0x40
SECOND_IN_PAIR = 0x80
SECONDARY = 0x100
QC_FAIL = 0x200


def revcomp(seq):
	return seq.translate(COMPLEMENT)[::-1]


def random_dna(length, rng):
	return ''.join(rng.choice('ACGT') for _ in range(length))


def mutate(seq, positions):
	"""Substitute the base at each 0-based position with a different base.

	The replacement is the next base in ACGT order, so the result is
	deterministic and always differs from the original.
	"""
	seq = list(seq)
	for p in positions:
		seq[p] = 'ACGT'[('ACGT'.index(seq[p]) + 1) % 4]
	return ''.join(seq)


def write_fasta(path, records, width=60):
	"""Write (name, sequence) pairs to a FASTA file; returns path."""
	with open(path, 'w') as f:
		for name, seq in records:
			f.write('>%s\n' % name)
			for i in range(0, len(seq), width):
				f.write(seq[i:i + width] + '\n')
	return path


def fastq_text(reads):
	"""FASTQ text for (name, seq, qual) triples; qual defaults to all 'I'."""
	out = []
	for read in reads:
		name, seq = read[0], read[1]
		qual = read[2] if len(read) > 2 and read[2] is not None else 'I' * len(seq)
		out.append('@%s\n%s\n+\n%s\n' % (name, seq, qual))
	return ''.join(out)


def write_fastq(path, reads, compress=None):
	"""Write reads to FASTQ, optionally gzip or bzip2 compressed; returns path."""
	text = fastq_text(reads).encode()
	if compress == 'gz':
		with gzip.open(path, 'wb') as f:
			f.write(text)
	elif compress == 'bz2':
		with bz2.open(path, 'wb') as f:
			f.write(text)
	else:
		with open(path, 'wb') as f:
			f.write(text)
	return path


def write_fasta_reads(path, reads):
	"""Write (name, seq) read pairs as multi-FASTA, one line per sequence."""
	with open(path, 'w') as f:
		for read in reads:
			f.write('>%s\n%s\n' % (read[0], read[1]))
	return path


def read_fastq(path):
	"""Parse a FASTQ file (plain, .gz or .bz2) into (name, seq, qual) triples."""
	opener = open
	if path.endswith('.gz'):
		opener = gzip.open
	elif path.endswith('.bz2'):
		opener = bz2.open
	with opener(path, 'rt') as f:
		lines = [line.rstrip('\n') for line in f]
	reads = []
	for i in range(0, len(lines) - 3, 4):
		assert lines[i].startswith('@'), lines[i]
		reads.append((lines[i][1:], lines[i + 1], lines[i + 3]))
	return reads


class SamRecord(object):
	"""One SAM alignment line. Optional fields are in `tags`, keyed by name,
	with integer values for type i and strings otherwise."""

	__slots__ = ('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext',
		'pnext', 'tlen', 'seq', 'qual', 'tags', 'tag_order', 'line')

	def __init__(self, line):
		fields = line.rstrip('\n').split('\t')
		assert len(fields) >= 11, line
		self.line = line.rstrip('\n')
		self.qname = fields[0]
		self.flag = int(fields[1])
		self.rname = fields[2]
		self.pos = int(fields[3])
		self.mapq = int(fields[4])
		self.cigar = fields[5]
		self.rnext = fields[6]
		self.pnext = int(fields[7])
		self.tlen = int(fields[8])
		self.seq = fields[9]
		self.qual = fields[10]
		self.tags = {}
		self.tag_order = []
		for field in fields[11:]:
			name, typ, value = field.split(':', 2)
			self.tags[name] = int(value) if typ == 'i' else value
			self.tag_order.append(name)

	@property
	def is_unmapped(self):
		return bool(self.flag & UNMAPPED)

	@property
	def is_reverse(self):
		return bool(self.flag & REVERSE)

	@property
	def is_secondary(self):
		return bool(self.flag & SECONDARY)

	def __repr__(self):
		return 'SamRecord(%r)' % self.line


def parse_sam(text):
	"""Split SAM text into (header lines, list of SamRecord)."""
	header, records = [], []
	for line in text.splitlines():
		if not line:
			continue
		if line.startswith('@'):
			header.append(line)
		else:
			records.append(SamRecord(line))
	return header, records


def cigar_ops(cigar):
	"""'10M2I5M' -> [(10, 'M'), (2, 'I'), (5, 'M')]."""
	ops, num = [], ''
	for c in cigar:
		if c.isdigit():
			num += c
		else:
			ops.append((int(num), c))
			num = ''
	return ops


def cigar_ref_length(cigar):
	return sum(n for n, op in cigar_ops(cigar) if op in 'MDN=X')


def cigar_query_length(cigar):
	return sum(n for n, op in cigar_ops(cigar) if op in 'MIS=X')


def mismatch_penalty(qual_char, mx=6, mn=2, phred_offset=33):
	"""bowtie2's quality-aware mismatch penalty, from the MANUAL's --mp:
	MN + floor((MX - MN) * (min(Q, 40) / 40))."""
	q = min(ord(qual_char) - phred_offset, 40)
	return mn + (mx - mn) * q // 40


class Bowtie2(object):
	"""Runs the bowtie2, bowtie2-build and bowtie2-inspect wrappers."""

	def __init__(self, bin_dir=BIN_DIR):
		self.bin_dir = bin_dir
		self.bowtie2 = os.path.join(bin_dir, 'bowtie2')
		self.build = os.path.join(bin_dir, 'bowtie2-build')
		self.inspect = os.path.join(bin_dir, 'bowtie2-inspect')

	def has_binary(self, name):
		return os.access(os.path.join(self.bin_dir, name), os.X_OK)

	def run(self, prog, args, input=None, check=True, cwd=None, env=None,
			timeout=300):
		"""Run one of the three wrappers with a list of arguments.

		Returns the CompletedProcess with text stdout and stderr. With
		check=True a non-zero exit raises AssertionError showing stderr.
		"""
		exe = {'bowtie2': self.bowtie2, 'bowtie2-build': self.build,
			'bowtie2-inspect': self.inspect}.get(prog, prog)
		cmd = [exe] + [str(a) for a in args]
		proc = subprocess.run(cmd, input=input, capture_output=True,
			text=True, cwd=cwd, env=env, timeout=timeout)
		if check:
			assert proc.returncode == 0, 'exit %d from %s\nstderr:\n%s' % (
				proc.returncode, ' '.join(cmd), proc.stderr)
		return proc

	def align(self, args, input=None, check=True, cwd=None, env=None):
		"""Run bowtie2 with SAM on stdout; returns (header, records, proc)."""
		proc = self.run('bowtie2', args, input=input, check=check, cwd=cwd,
			env=env)
		header, records = parse_sam(proc.stdout) if proc.returncode == 0 else ([], [])
		return header, records, proc

	def build_index(self, fasta, prefix, extra=(), check=True):
		"""bowtie2-build a FASTA (or comma-separated list) into prefix."""
		return self.run('bowtie2-build', list(extra) + ['--quiet', fasta, prefix],
			check=check)


def by_name(records):
	"""Map qname -> list of records, in output order."""
	out = {}
	for r in records:
		out.setdefault(r.qname, []).append(r)
	return out


def primary(records):
	"""Records that are neither secondary nor supplementary."""
	return [r for r in records if not r.flag & (SECONDARY | 0x800)]


def which(name):
	return shutil.which(name)
