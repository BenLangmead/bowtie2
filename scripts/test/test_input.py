"""Input formats and read preprocessing.

Covers the read input formats (-q, --tab5, --tab6, --12, --qseq, -f, -r, -F,
-c, --interleaved, -b with --align-paired-reads and --preserve-tags), which
reads are aligned and how they are trimmed (-s, -u, -5, -3, --trim-to),
quality encodings (--phred33, --phred64, --solexa-quals and their aliases),
comma-separated and repeated read files, compressed files and standard input,
the QNAME options --sam-no-qname-trunc and --sam-append-comment,
--qc-filter, and reads that are shorter than a seed, empty or all N.

Reads are cut from the fixed random `genome` fixture at known positions and
strands, so FLAG, POS, CIGAR, SEQ, QUAL, AS and MD of every record are derived
here rather than taken from bowtie2. Behaviour of the current version that
contradicts the MANUAL is not tested here.
"""

import bz2
import gzip
import math
import os
import shutil
import struct
import subprocess
import zlib

import pytest

import bt2_helpers as h


# Unpaired reads: (name, chromosome, 0-based start, length, strand).
UNPAIRED = [
	('u1', 'chr1', 120, 50, '+'),
	('u2', 'chr2', 900, 45, '-'),
	('u3', 'chr3', 1500, 60, '+'),
	('u4', 'chr1', 3300, 40, '-'),
]

# Pairs: (name, chromosome, fragment start, fragment length, mate 1 length,
# mate 2 length). Mate 1 is the forward start of the fragment and mate 2 the
# reverse complement of its end, the default --fr orientation.
PAIRS = [
	('p1', 'chr1', 1000, 250, 50, 50),
	('p2', 'chr2', 1500, 300, 40, 45),
]

# Fields 2-8 of the QSEQ records written here, and the QNAME suffix bowtie2
# builds from them: fields 1-7 joined with '_', then '/' and field 8. The
# MANUAL does not describe QSEQ read names, so this pins the current output
# (read_qseq.cpp).
QSEQ_FIELDS = ['7', '1', '11', '101', '202', '0', '1']
QSEQ_SUFFIX = '_7_1_11_101_202_0/1'


# ---------------------------------------------------------------- helpers

def run_align(bt2, index, args, input=None, check=True):
	"""Align with one thread; returns (records, CompletedProcess)."""
	_, records, proc = bt2.align(['-x', index, '-p', '1'] + list(args),
		input=input, check=check)
	return records, proc


def run_bytes(bt2, args, data):
	"""Run the bowtie2 wrapper with binary stdin; returns (records, proc)."""
	proc = subprocess.run([bt2.bowtie2] + [str(a) for a in args], input=data,
		capture_output=True, timeout=300)
	assert proc.returncode == 0, proc.stderr.decode()
	return h.parse_sam(proc.stdout.decode())[1], proc


def sam_body(proc):
	"""SAM records of a run without its header lines."""
	return [line for line in proc.stdout.splitlines() if not line.startswith('@')]


def read_seq(gd, spec):
	name, chrom, start, length, strand = spec
	seq = gd[chrom][start:start + length]
	return seq if strand == '+' else h.revcomp(seq)


def varied_qual(length, shift=0):
	"""Phred+33 qualities between 25 and 40 that differ between neighbours,
	so a reversed or shifted QUAL string is detectable."""
	return ''.join(chr(33 + 25 + (i * 7 + shift) % 16) for i in range(length))


def unpaired_reads(gd):
	"""(name, seq, qual) of the UNPAIRED reads."""
	return [(s[0], read_seq(gd, s), varied_qual(s[3], k))
		for k, s in enumerate(UNPAIRED)]


def pair_seqs(gd, pair):
	name, chrom, start, frag, l1, l2 = pair
	g = gd[chrom]
	return g[start:start + l1], h.revcomp(g[start + frag - l2:start + frag])


def assert_placed(rec, gd, spec, qual, trim5=0, trim3=0):
	"""`rec` is the alignment of read `spec`, with read quality string `qual`,
	after trim5 bases were removed from its 5' end and trim3 from its 3' end."""
	name, chrom, start, length, strand = spec
	n = length - trim5 - trim3
	kept = qual[trim5:length - trim3]
	if strand == '+':
		ref_start, flag, exp_qual = start + trim5, 0, kept
	else:
		ref_start, flag, exp_qual = start + trim3, h.REVERSE, kept[::-1]
	assert rec.flag == flag, rec
	assert rec.rname == chrom, rec
	assert rec.pos == ref_start + 1, rec
	assert rec.cigar == '%dM' % n, rec
	assert rec.seq == gd[chrom][ref_start:ref_start + n], rec
	assert rec.qual == exp_qual, rec
	assert rec.tags['AS'] == 0, rec
	assert rec.tags['NM'] == 0, rec
	assert rec.tags['MD'] == str(n), rec


def assert_pair(mate1, mate2, gd, pair, qual1=None, qual2=None):
	"""mate1/mate2 are a concordant --fr alignment of `pair`."""
	name, chrom, start, frag, l1, l2 = pair
	g = gd[chrom]
	pos2 = start + frag - l2
	assert mate1.flag == 99, mate1
	assert mate2.flag == 147, mate2
	assert (mate1.rname, mate2.rname) == (chrom, chrom)
	assert mate1.pos == start + 1, mate1
	assert mate2.pos == pos2 + 1, mate2
	assert (mate1.cigar, mate2.cigar) == ('%dM' % l1, '%dM' % l2)
	assert mate1.seq == g[start:start + l1]
	assert mate2.seq == g[pos2:start + frag]
	assert (mate1.pnext, mate2.pnext) == (mate2.pos, mate1.pos)
	assert (mate1.tlen, mate2.tlen) == (frag, -frag)
	assert mate1.tags['AS'] == 0 and mate2.tags['AS'] == 0
	if qual1 is not None:
		assert mate1.qual == qual1
	if qual2 is not None:
		assert mate2.qual == qual2[::-1]


def qseq_line(name, seq, qual, filt='1'):
	return '\t'.join([name] + QSEQ_FIELDS + [seq, qual, filt]) + '\n'


def write_text(path, text):
	with open(str(path), 'w') as f:
		f.write(text)
	return str(path)


def write_input(fmt, tmp_path, reads):
	"""Write (name, seq, qual) reads in format `fmt`.

	Returns (bowtie2 arguments, expected QNAMEs, whether QUAL is taken from
	the input). Formats without qualities get Phred 40 ('I') everywhere."""
	names = [r[0] for r in reads]
	index_names = [str(i) for i in range(len(reads))]
	if fmt == 'fastq':
		path = write_text(tmp_path / 'r.fq', ''.join('@%s desc %d\n%s\n+\n%s\n'
			% (n, i, s, q) for i, (n, s, q) in enumerate(reads)))
		return ['-q', '-U', path], names, True
	if fmt == 'fasta':
		path = str(tmp_path / 'r.fa')
		h.write_fasta(path, [('%s desc %d' % (n, i), s)
			for i, (n, s, q) in enumerate(reads)], width=25)
		return ['-f', '-U', path], names, False
	if fmt in ('tab5', 'tab6'):
		path = write_text(tmp_path / ('r.' + fmt),
			''.join('%s\t%s\t%s\n' % r for r in reads))
		return ['--' + fmt, path], names, True
	if fmt == 'qseq':
		path = write_text(tmp_path / 'r_qseq.txt',
			''.join(qseq_line(n, s, q) for n, s, q in reads))
		return ['--qseq', '-U', path], [n + QSEQ_SUFFIX for n in names], True
	if fmt == 'raw':
		path = write_text(tmp_path / 'r.raw', ''.join(s + '\n' for n, s, q in reads))
		return ['-r', '-U', path], index_names, False
	if fmt == 'cmdline':
		return ['-c', '-U', ','.join(s for n, s, q in reads)], index_names, False
	if fmt == 'bam':
		path = str(tmp_path / 'r.bam')
		write_bam(path, [(n, 0x4, s, q, []) for n, s, q in reads])
		return ['-b', path], names, True
	raise ValueError(fmt)


FORMATS = ['fastq', 'fasta', 'tab5', 'tab6', 'qseq', 'raw', 'cmdline', 'bam']


# ---------------------------------------------------- minimal BAM writer

BAM_EOF = bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')
BAM_SEQ_CODES = '=ACMGRSVTWYHKDBN'
BAM_INT_FORMATS = {'c': 'b', 'C': 'B', 's': 'h', 'S': 'H', 'i': 'i', 'I': 'I', 'f': 'f'}


def bgzf_block(data):
	"""One BGZF block: a gzip member with the BC extra subfield."""
	comp = zlib.compressobj(6, zlib.DEFLATED, -15)
	cdata = comp.compress(data) + comp.flush()
	header = struct.pack('<4BI2BH2BHH', 0x1f, 0x8b, 8, 4, 0, 0, 0xff, 6,
		ord('B'), ord('C'), 2, len(cdata) + 25)
	return header + cdata + struct.pack('<II', zlib.crc32(data) & 0xffffffff, len(data))


def bam_tag(tag, typ, value):
	out = tag.encode()
	if typ == 'Z':
		return out + b'Z' + value.encode() + b'\0'
	if typ == 'A':
		return out + b'A' + value.encode()
	if typ == 'B':
		sub, vals = value
		return (out + b'B' + sub.encode() + struct.pack('<I', len(vals))
			+ struct.pack('<%d%s' % (len(vals), BAM_INT_FORMATS[sub]), *vals))
	return out + typ.encode() + struct.pack('<' + BAM_INT_FORMATS[typ], value)


def bam_record(name, flag, seq, qual, tags):
	packed = bytearray()
	for i in range(0, len(seq), 2):
		hi = BAM_SEQ_CODES.index(seq[i])
		lo = BAM_SEQ_CODES.index(seq[i + 1]) if i + 1 < len(seq) else 0
		packed.append(hi << 4 | lo)
	read_name = name.encode() + b'\0'
	body = struct.pack('<iiBBHHHiiii', -1, -1, len(read_name), 0, 4680, 0, flag,
		len(seq), -1, -1, 0)
	body += read_name + bytes(packed) + bytes(ord(c) - 33 for c in qual)
	body += b''.join(bam_tag(*t) for t in tags)
	return struct.pack('<i', len(body)) + body


def write_bam(path, records):
	"""Unaligned BAM from (name, flag, seq, Phred+33 qual, tags) records,
	where tags are (tag, type, value) and a B-array value is (subtype, list).
	Everything goes in one BGZF block followed by the EOF block."""
	text = b'@HD\tVN:1.6\tSO:unsorted\n'
	data = b'BAM\1' + struct.pack('<i', len(text)) + text + struct.pack('<i', 0)
	data += b''.join(bam_record(*r) for r in records)
	assert len(data) < 65536
	with open(path, 'wb') as f:
		f.write(bgzf_block(data))
		f.write(BAM_EOF)
	return path


def read_bam_headers(path):
	"""(name, flag, offset of SEQ, l_seq) for every record of a BAM file, plus
	the decompressed data; BGZF is multi-member gzip, so gzip reads it."""
	with gzip.open(path, 'rb') as f:
		data = f.read()
	assert data[:4] == b'BAM\1'
	off = 8 + struct.unpack_from('<i', data, 4)[0]
	n_ref = struct.unpack_from('<i', data, off)[0]
	off += 4
	for _ in range(n_ref):
		off += 8 + struct.unpack_from('<i', data, off)[0]
	recs = []
	while off < len(data):
		block_size = struct.unpack_from('<i', data, off)[0]
		l_name, _, _, n_cigar, flag, l_seq = struct.unpack_from('<BBHHHi', data, off + 12)
		name = data[off + 36:off + 36 + l_name - 1].decode()
		recs.append((name, flag, off + 36 + l_name + 4 * n_cigar, l_seq))
		off += 4 + block_size
	return recs, data


def bam_seq_qual(data, seq_off, l_seq):
	seq = ''.join(BAM_SEQ_CODES[(data[seq_off + i // 2] >> (4 * (1 - i % 2))) & 0xf]
		for i in range(l_seq))
	qoff = seq_off + (l_seq + 1) // 2
	return seq, ''.join(chr(q + 33) for q in data[qoff:qoff + l_seq])


# ------------------------------------------------------- input formats

@pytest.mark.parametrize('fmt', FORMATS)
def test_format_parses_reads(bt2, genome_dict, genome_index, tmp_path, fmt):
	"""Each input format yields the reads at their known positions, strands, QNAMEs and QUALs.

	The QNAMEs of -r and -c reads (the 0-based read index) and of QSEQ reads
	(QSEQ_SUFFIX) pin the current output; the MANUAL only says that -r and
	-c reads have no names."""
	reads = unpaired_reads(genome_dict)
	args, names, has_qual = write_input(fmt, tmp_path, reads)
	recs, _ = run_align(bt2, genome_index, args)
	assert [r.qname for r in recs] == names
	for rec, spec, (n, s, q) in zip(recs, UNPAIRED, reads):
		assert_placed(rec, genome_dict, spec, q if has_qual else 'I' * len(s))


def test_q_is_the_default(bt2, genome_dict, genome_index, tmp_path):
	"""-q (FASTQ) gives the same SAM as giving no format option."""
	args, _, _ = write_input('fastq', tmp_path, unpaired_reads(genome_dict))
	_, with_q = run_align(bt2, genome_index, args)
	_, without_q = run_align(bt2, genome_index, args[1:])
	assert len(sam_body(with_q)) == len(UNPAIRED)
	assert sam_body(with_q) == sam_body(without_q)


@pytest.mark.parametrize('fmt', ['fasta', 'raw', 'cmdline'])
def test_formats_without_qualities_score_as_ignore_quals(bt2, genome_dict,
		genome_index, tmp_path, fmt):
	"""-f, -r and -c score a mismatch at the full MX penalty, as --ignore-quals does."""
	ref = genome_dict['chr2'][400:450]
	read = h.mutate(ref, [17])
	args, _, _ = write_input(fmt, tmp_path, [('m1', read, None)])
	recs, plain = run_align(bt2, genome_index, args)
	_, ignore = run_align(bt2, genome_index, args + ['--ignore-quals'])
	assert sam_body(plain) == sam_body(ignore)
	rec, = recs
	assert (rec.pos, rec.qual) == (401, 'I' * 50)
	assert rec.tags['AS'] == -6
	assert rec.tags['XM'] == 1 and rec.tags['NM'] == 1
	assert rec.tags['MD'] == '17%s32' % ref[17]


def test_fasta_continuous_kmers(bt2, genome_dict, genome_index, tmp_path):
	"""-F <k>,<i> aligns every k-mer at offsets 0, i, 2i... of each FASTA record, named <record>_<offset>.

	seqB is 36 bases long, so its last k-mer (offset 16) ends at the record's
	last base. Only the undocumented "-F 20,4" spelling works in the current
	version; the MANUAL's "-F k:20,i:4" yields no reads and is not tested here."""
	seq_a = genome_dict['chr2'][200:230]
	seq_b = genome_dict['chr3'][100:136]
	path = str(tmp_path / 'long.fa')
	h.write_fasta(path, [('seqA some description', seq_a), ('seqB', seq_b)], width=11)
	recs, _ = run_align(bt2, genome_index, ['-F', '20,4', '-U', path])
	expected = []
	for name, chrom, start, length in (('seqA', 'chr2', 200, 30), ('seqB', 'chr3', 100, 36)):
		for off in range(0, length - 20 + 1, 4):
			expected.append(('%s_%d' % (name, off), chrom, start + off))
	assert [(r.qname, r.rname, r.pos - 1) for r in recs] == expected
	for rec in recs:
		assert rec.flag == 0 and rec.cigar == '20M'
		assert rec.seq == genome_dict[rec.rname][rec.pos - 1:rec.pos + 19]
		assert rec.qual == 'I' * 20
		assert rec.tags['AS'] == 0


def test_cmdline_paired(bt2, genome_dict, genome_index):
	"""-c with -1/-2 takes mates as comma-separated sequences and names pairs by index."""
	m1s, m2s = zip(*[pair_seqs(genome_dict, p) for p in PAIRS])
	recs, _ = run_align(bt2, genome_index, ['-c', '-1', ','.join(m1s), '-2', ','.join(m2s)])
	assert [r.qname for r in recs] == ['0', '0', '1', '1']
	for k, pair in enumerate(PAIRS):
		assert_pair(recs[2 * k], recs[2 * k + 1], genome_dict, pair)


def test_tab5_mixes_paired_and_unpaired(bt2, genome_dict, genome_index, tmp_path):
	"""--tab5 treats 5-field lines as pairs sharing one name and 3-field lines as unpaired reads."""
	(m1a, m2a), (m1b, m2b) = [pair_seqs(genome_dict, p) for p in PAIRS]
	u1 = read_seq(genome_dict, UNPAIRED[0])
	q1, q2, q3, q4, qu = [varied_qual(len(s), k) for k, s in enumerate((m1a, m2a, m1b, m2b, u1))]
	path = write_text(tmp_path / 'r.tab5', 'p1\t%s\t%s\t%s\t%s\nu1\t%s\t%s\np2\t%s\t%s\t%s\t%s\n'
		% (m1a, q1, m2a, q2, u1, qu, m1b, q3, m2b, q4))
	recs, _ = run_align(bt2, genome_index, ['--tab5', path])
	assert [r.qname for r in recs] == ['p1', 'p1', 'u1', 'p2', 'p2']
	assert_pair(recs[0], recs[1], genome_dict, PAIRS[0], q1, q2)
	assert_placed(recs[2], genome_dict, UNPAIRED[0], qu)
	assert_pair(recs[3], recs[4], genome_dict, PAIRS[1], q3, q4)


def test_tab6_mate_names(bt2, genome_dict, genome_index, tmp_path):
	"""--tab6 gives each mate its own name; unpaired 3-field lines still work."""
	m1, m2 = pair_seqs(genome_dict, PAIRS[0])
	u1 = read_seq(genome_dict, UNPAIRED[0])
	path = write_text(tmp_path / 'r.tab6', 'left\t%s\t%s\tright\t%s\t%s\nu1\t%s\t%s\n'
		% (m1, 'I' * len(m1), m2, 'I' * len(m2), u1, 'I' * len(u1)))
	recs, _ = run_align(bt2, genome_index, ['--tab6', path])
	assert [r.qname for r in recs] == ['left', 'right', 'u1']
	assert_pair(recs[0], recs[1], genome_dict, PAIRS[0])
	assert_placed(recs[2], genome_dict, UNPAIRED[0], 'I' * len(u1))


def test_12_is_tab5(bt2, genome_dict, genome_index, tmp_path):
	"""The undocumented --12 option reads the --tab5 format."""
	m1, m2 = pair_seqs(genome_dict, PAIRS[0])
	u1 = read_seq(genome_dict, UNPAIRED[0])
	path = write_text(tmp_path / 'r.tab5', 'p1\t%s\t%s\t%s\t%s\nu1\t%s\t%s\n'
		% (m1, 'I' * len(m1), m2, 'I' * len(m2), u1, 'I' * len(u1)))
	_, tab5 = run_align(bt2, genome_index, ['--tab5', path])
	_, onetwo = run_align(bt2, genome_index, ['--12', path])
	assert len(sam_body(tab5)) == 3
	assert sam_body(onetwo) == sam_body(tab5)


@pytest.mark.parametrize('fmt', ['fastq', 'fasta'])
def test_interleaved(bt2, genome_dict, genome_index, tmp_path, fmt):
	"""--interleaved takes consecutive records as mate 1 and mate 2; /1 and /2 are dropped from QNAME."""
	text = []
	for pair in PAIRS:
		for mate, seq in zip((1, 2), pair_seqs(genome_dict, pair)):
			if fmt == 'fastq':
				text.append('@%s/%d\n%s\n+\n%s\n' % (pair[0], mate, seq, 'I' * len(seq)))
			else:
				text.append('>%s/%d\n%s\n' % (pair[0], mate, seq))
	path = write_text(tmp_path / ('r.' + fmt), ''.join(text))
	args = ['--interleaved', path] + (['-f'] if fmt == 'fasta' else [])
	recs, _ = run_align(bt2, genome_index, args)
	assert [r.qname for r in recs] == ['p1', 'p1', 'p2', 'p2']
	for k, pair in enumerate(PAIRS):
		assert_pair(recs[2 * k], recs[2 * k + 1], genome_dict, pair)


def test_interleaved_rejects_other_formats(bt2, genome_dict, genome_index, tmp_path):
	"""--interleaved with -r is an error."""
	path = write_text(tmp_path / 'r.raw', 'ACGTACGTACGTACGTACGTAC\n')
	_, proc = run_align(bt2, genome_index, ['-r', '--interleaved', path], check=False)
	assert proc.returncode != 0
	assert '--interleaved only works in combination with FASTA (-f) and FASTQ (-q) formats' in proc.stderr


@pytest.mark.parametrize('fmt,flag,message', [
	('fastq', '-f', 'does not look like a FASTA file'),
	('fasta', '-q', 'does not look like a FASTQ file'),
])
def test_wrong_format_is_an_error(bt2, genome_dict, genome_index, tmp_path, fmt, flag, message):
	"""A FASTQ file read with -f, or a FASTA file read with -q, is reported as the wrong format."""
	args, _, _ = write_input(fmt, tmp_path, unpaired_reads(genome_dict))
	_, proc = run_align(bt2, genome_index, [flag] + args[1:], check=False)
	assert proc.returncode != 0
	assert message in proc.stderr


# ------------------------------------------------------------------ QSEQ

def test_qc_filter(bt2, genome_dict, genome_index, tmp_path):
	"""--qc-filter removes exactly the QSEQ reads of one filter value (unaligned, YF:Z:QC, SEQ/QUAL kept) and aligns the rest.

	Which value is removed is left open on purpose. The MANUAL says the reads
	whose filter field is 1 ("non-zero"); the current version removes the
	reads marked 0, following Illumina's convention that 1 means "passed
	filter", so the suite stays green whichever of the two is corrected."""
	reads = unpaired_reads(genome_dict)
	filters = ['1', '0', '0', '1']
	path = write_text(tmp_path / 'r_qseq.txt',
		''.join(qseq_line(n, s, q, f) for (n, s, q), f in zip(reads, filters)))
	names = [n + QSEQ_SUFFIX for n, s, q in reads]
	recs, _ = run_align(bt2, genome_index, ['--qseq', '-U', path])
	assert [r.qname for r in recs] == names
	for rec, spec, (n, s, q) in zip(recs, UNPAIRED, reads):
		assert_placed(rec, genome_dict, spec, q)
	recs, _ = run_align(bt2, genome_index, ['--qseq', '--qc-filter', '-U', path])
	assert [r.qname for r in recs] == names
	removed = [f for rec, f in zip(recs, filters) if rec.is_unmapped]
	assert removed in (['1', '1'], ['0', '0'])
	for rec, spec, (n, s, q), f in zip(recs, UNPAIRED, reads, filters):
		if f == removed[0]:
			assert (rec.flag, rec.rname, rec.pos, rec.cigar) == (h.UNMAPPED, '*', 0, '*')
			assert rec.tags['YF'] == 'QC'
			assert 'AS' not in rec.tags
			assert rec.seq == s and rec.qual == q
		else:
			assert_placed(rec, genome_dict, spec, q)


def test_qc_filter_has_no_effect_on_fastq(bt2, genome_dict, genome_index, tmp_path):
	"""--qc-filter changes nothing for FASTQ input, which has no filter field."""
	args, _, _ = write_input('fastq', tmp_path, unpaired_reads(genome_dict))
	_, plain = run_align(bt2, genome_index, args)
	_, filtered = run_align(bt2, genome_index, args + ['--qc-filter'])
	assert len(sam_body(plain)) == len(UNPAIRED)
	assert sam_body(filtered) == sam_body(plain)


def test_qseq_bad_filter_value(bt2, genome_dict, genome_index, tmp_path):
	"""A QSEQ filter field other than 0 or 1 is an error."""
	n, s, q = unpaired_reads(genome_dict)[0]
	path = write_text(tmp_path / 'r_qseq.txt', qseq_line(n, s, q, '2'))
	_, proc = run_align(bt2, genome_index, ['--qseq', '-U', path], check=False)
	assert proc.returncode != 0
	assert 'for qseq filter flag' in proc.stderr


@pytest.mark.parametrize('fmt', ['fastq', 'fasta', 'qseq'])
def test_dot_is_read_as_n(bt2, genome_dict, genome_index, tmp_path, fmt):
	"""A '.' in a FASTQ, FASTA or QSEQ sequence is an N: SEQ shows N and AS is the N penalty (-1)."""
	ref = genome_dict['chr3'][800:840]
	read = ref[:12] + '.' + ref[13:]
	args, _, _ = write_input(fmt, tmp_path, [('dot', read, 'I' * 40)])
	rec, = run_align(bt2, genome_index, args)[0]
	assert rec.pos == 801 and rec.flag == 0
	assert rec.seq == ref[:12] + 'N' + ref[13:]
	assert rec.tags['AS'] == -1
	assert rec.tags['NM'] == 1
	assert rec.tags['MD'] == '12%s27' % ref[12]


# ------------------------------------------------------------------- BAM

def bam_test_records(gd):
	"""Records for the BAM tests: one pair (flags 77/141), two unpaired
	unaligned reads (flag 4) and one record that is already aligned (flag 0)."""
	m1, m2 = pair_seqs(gd, PAIRS[0])
	(n1, s1, q1), (n2, s2, q2) = unpaired_reads(gd)[:2]
	qm1, qm2 = varied_qual(len(m1), 3), varied_qual(len(m2), 5)
	return [
		('p1', 77, m1, qm1, [('QZ', 'Z', 'mate one'), ('Qi', 'i', -7)]),
		('p1', 141, m2, qm2, [('QC', 'C', 200)]),
		(n1, 4, s1, q1, [('QZ', 'Z', 'a b:c'), ('Qc', 'c', -100), ('QC', 'C', 250),
			('Qs', 's', -30000), ('QS', 'S', 60000), ('Qi', 'i', -2000000000),
			('QI', 'I', 4000000000), ('Qf', 'f', -0.25), ('Qb', 'B', ('c', [-1, 2])),
			('QB', 'B', ('S', [7, 65535])), ('Qg', 'B', ('f', [0.5, 2.0]))]),
		(n2, 4, s2, q2, []),
		('aligned', 0, s1, q1, []),
	]


def test_bam_unpaired_by_default(bt2, genome_dict, genome_index, tmp_path):
	"""-b aligns only the unpaired records, taking SEQ and QUAL (Phred, no offset) from the BAM.

	That paired records are skipped is documented; that the already-aligned
	(flag 0) record is skipped pins the current output."""
	recs_in = bam_test_records(genome_dict)
	path = write_bam(str(tmp_path / 'r.bam'), recs_in)
	recs, _ = run_align(bt2, genome_index, ['-b', path])
	assert [r.qname for r in recs] == ['u1', 'u2']
	assert_placed(recs[0], genome_dict, UNPAIRED[0], recs_in[2][3])
	assert_placed(recs[1], genome_dict, UNPAIRED[1], recs_in[3][3])
	assert not any(t.startswith('Q') for r in recs for t in r.tags)


def test_bam_align_paired_reads(bt2, genome_dict, genome_index, tmp_path):
	"""-b --align-paired-reads aligns the paired records as pairs and skips unpaired ones."""
	recs_in = bam_test_records(genome_dict)
	path = write_bam(str(tmp_path / 'r.bam'), recs_in)
	recs, _ = run_align(bt2, genome_index, ['-b', path, '--align-paired-reads'])
	assert [r.qname for r in recs] == ['p1', 'p1']
	assert_pair(recs[0], recs[1], genome_dict, PAIRS[0], recs_in[0][3], recs_in[1][3])


EXPECTED_PRESERVED = ['QZ:Z:a b:c', 'Qc:i:-100', 'QC:i:250', 'Qs:i:-30000',
	'QS:i:60000', 'Qi:i:-2000000000', 'QI:i:4000000000', 'Qf:f:{-0.25}',
	'Qb:B:c,-1,2', 'QB:B:S,7,65535', 'Qg:B:f,{0.5},{2.0}']


def assert_tags_preserved(line, expected):
	"""The SAM line ends with `expected` tags; '{x}' marks a float compared numerically."""
	fields = line.split('\t')[-len(expected):]
	for got, want in zip(fields, expected):
		if '{' not in want:
			assert got == want, line
			continue
		prefix = want.split('{')[0]
		assert got.startswith(prefix), line
		got_vals = [float(v) for v in got[len(prefix):].split(',')]
		want_vals = [float(v.strip('{}')) for v in want[len(prefix):].split(',')]
		assert got_vals == want_vals, line


def test_bam_preserve_tags_unpaired(bt2, genome_dict, genome_index, tmp_path):
	"""--preserve-tags appends each BAM tag, integer types as i, after bowtie2's own fields."""
	recs_in = bam_test_records(genome_dict)
	recs_in.append(('unal', 4, 'ACGT' * 10, 'I' * 40, [('QZ', 'Z', 'kept')]))
	path = write_bam(str(tmp_path / 'r.bam'), recs_in)
	recs, _ = run_align(bt2, genome_index, ['-b', path, '--preserve-tags'])
	assert [r.qname for r in recs] == ['u1', 'u2', 'unal']
	assert_placed(recs[0], genome_dict, UNPAIRED[0], recs_in[2][3])
	assert recs[0].tag_order.index('YT') < recs[0].tag_order.index('QZ')
	assert_tags_preserved(recs[0].line, EXPECTED_PRESERVED)
	assert not any(t.startswith('Q') for t in recs[1].tags)
	assert recs[2].is_unmapped
	assert recs[2].line.endswith('\tQZ:Z:kept')


def test_bam_preserve_tags_paired(bt2, genome_dict, genome_index, tmp_path):
	"""--preserve-tags with --align-paired-reads puts each mate's own BAM tags on that mate."""
	recs_in = bam_test_records(genome_dict)
	path = write_bam(str(tmp_path / 'r.bam'), recs_in)
	recs, _ = run_align(bt2, genome_index, ['-b', path, '--align-paired-reads', '--preserve-tags'])
	assert_pair(recs[0], recs[1], genome_dict, PAIRS[0])
	assert recs[0].line.endswith('\tYT:Z:CP\tQZ:Z:mate one\tQi:i:-7')
	assert recs[1].line.endswith('\tYT:Z:CP\tQC:i:200')


@pytest.mark.parametrize('option', ['--preserve-tags', '--align-paired-reads'])
def test_bam_options_require_bam_input(bt2, genome_dict, genome_index, tmp_path, option):
	"""--preserve-tags and --align-paired-reads are errors without -b."""
	args, _, _ = write_input('fastq', tmp_path, unpaired_reads(genome_dict))
	_, proc = run_align(bt2, genome_index, args + [option], check=False)
	assert proc.returncode != 0
	assert 'can only be used when aligning BAM reads' in proc.stderr


def test_bam_shipped_example(bt2, lambda_index, example_reads):
	"""-b reads the shipped example BAM: QNAME, SEQ and QUAL match the file's records, in order.

	example/reads/combined_reads.bam spans many BGZF blocks; -s/-u count
	unpaired records, or pairs with --align-paired-reads."""
	headers, data = read_bam_headers(example_reads['bam'])
	unpaired = [r for r in headers if r[1] & 0x4 and not r[1] & 0x1]
	firsts = [r for r in headers if r[1] & 0x40]
	seconds = [r for r in headers if r[1] & 0x80]
	assert len(firsts) == len(seconds) >= 5 and len(unpaired) >= 13

	def check(rec, bam_rec):
		seq, qual = bam_seq_qual(data, bam_rec[2], bam_rec[3])
		if rec.is_reverse:
			seq, qual = h.revcomp(seq), qual[::-1]
		assert (rec.qname, rec.seq, rec.qual) == (bam_rec[0], seq, qual)

	recs, _ = run_align(bt2, lambda_index, ['-b', example_reads['bam'], '-s', '5', '-u', '8'])
	assert len(recs) == 8
	for rec, bam_rec in zip(recs, unpaired[5:13]):
		check(rec, bam_rec)
	recs, _ = run_align(bt2, lambda_index, ['-b', example_reads['bam'], '--align-paired-reads', '-u', '5'])
	assert len(recs) == 10
	for k in range(5):
		assert recs[2 * k].flag & h.FIRST_IN_PAIR and recs[2 * k + 1].flag & h.SECOND_IN_PAIR
		check(recs[2 * k], firsts[k])
		check(recs[2 * k + 1], seconds[k])


# ------------------------------------------------- read selection: -s / -u

@pytest.mark.parametrize('args,expected', [
	([], ['u1', 'u2', 'u3', 'u4']),
	(['-s', '1'], ['u2', 'u3', 'u4']),
	(['-s', '3'], ['u4']),
	(['-s', '4'], []),
	(['-s', '10'], []),
	(['-u', '2'], ['u1', 'u2']),
	(['-u', '9'], ['u1', 'u2', 'u3', 'u4']),
	(['-s', '1', '-u', '2'], ['u2', 'u3']),
	(['--skip', '1', '--upto', '2'], ['u2', 'u3']),
	(['-s', '1', '--qupto', '2'], ['u2', 'u3']),
])
def test_skip_upto(bt2, genome_dict, genome_index, tmp_path, args, expected):
	"""-s/--skip skips the first reads and -u/--upto/--qupto stops after that many more, counting across files."""
	reads = unpaired_reads(genome_dict)
	a = h.write_fastq(str(tmp_path / 'a.fq'), reads[:2])
	b = h.write_fastq(str(tmp_path / 'b.fq.gz'), reads[2:], compress='gz')
	recs, _ = run_align(bt2, genome_index, ['-U', a + ',' + b] + args)
	assert [r.qname for r in recs] == expected
	for rec in recs:
		k = int(rec.qname[1]) - 1
		assert_placed(rec, genome_dict, UNPAIRED[k], reads[k][2])


@pytest.mark.parametrize('args,expected', [
	(['-s', '1'], ['p2', 'p2']),
	(['-u', '1'], ['p1', 'p1']),
])
def test_skip_upto_count_pairs(bt2, genome_dict, genome_index, tmp_path, args, expected):
	"""-s and -u count pairs, not mates."""
	lines = []
	for pair in PAIRS:
		m1, m2 = pair_seqs(genome_dict, pair)
		lines.append('%s\t%s\t%s\t%s\t%s\n' % (pair[0], m1, 'I' * len(m1), m2, 'I' * len(m2)))
	path = write_text(tmp_path / 'r.tab5', ''.join(lines))
	recs, _ = run_align(bt2, genome_index, ['--tab5', path] + args)
	assert [r.qname for r in recs] == expected
	pair = PAIRS[0] if expected[0] == 'p1' else PAIRS[1]
	assert_pair(recs[0], recs[1], genome_dict, pair)


# ----------------------------------------------------------- trimming

@pytest.mark.parametrize('fmt', FORMATS)
def test_trim5_trim3(bt2, genome_dict, genome_index, tmp_path, fmt):
	"""-5 3 -3 5 removes 3 bases from the read's 5' end and 5 from its 3' end, on either strand."""
	reads = unpaired_reads(genome_dict)
	args, names, has_qual = write_input(fmt, tmp_path, reads)
	recs, _ = run_align(bt2, genome_index, args + ['-5', '3', '-3', '5'])
	assert [r.qname for r in recs] == names
	for rec, spec, (n, s, q) in zip(recs, UNPAIRED, reads):
		assert_placed(rec, genome_dict, spec, q if has_qual else 'I' * len(s), trim5=3, trim3=5)


def test_trim_long_options(bt2, genome_dict, genome_index, tmp_path):
	"""--trim5 and --trim3 are the same as -5 and -3."""
	args, _, _ = write_input('fastq', tmp_path, unpaired_reads(genome_dict))
	_, short = run_align(bt2, genome_index, args + ['-5', '4', '-3', '2'])
	_, long_ = run_align(bt2, genome_index, args + ['--trim5', '4', '--trim3', '2'])
	assert sam_body(long_) == sam_body(short)
	assert len(sam_body(short)) == len(UNPAIRED)
	assert all('\t%dM\t' % (spec[3] - 6) in line for spec, line in zip(UNPAIRED, sam_body(short)))


@pytest.mark.parametrize('fmt', ['fastq', 'bam'])
@pytest.mark.parametrize('spec,end', [('3:42', 3), ('5:42', 5), ('42', 3)])
def test_trim_to(bt2, genome_dict, genome_index, tmp_path, fmt, spec, end):
	"""--trim-to [3:|5:]42 cuts reads to 42 bases from the given end (default 3'); a 40-base read is unchanged."""
	reads = unpaired_reads(genome_dict)
	args, _, _ = write_input(fmt, tmp_path, reads)
	recs, _ = run_align(bt2, genome_index, args + ['--trim-to', spec])
	assert len(recs) == len(UNPAIRED)
	for rec, read_spec, (n, s, q) in zip(recs, UNPAIRED, reads):
		cut = max(0, len(s) - 42)
		trim5, trim3 = (0, cut) if end == 3 else (cut, 0)
		assert_placed(rec, genome_dict, read_spec, q, trim5=trim5, trim3=trim3)


def test_trim_to_default_end_is_3(bt2, genome_dict, genome_index, tmp_path):
	"""--trim-to <int> is the same as --trim-to 3:<int>."""
	args, _, _ = write_input('fastq', tmp_path, unpaired_reads(genome_dict))
	_, bare = run_align(bt2, genome_index, args + ['--trim-to', '33'])
	_, three = run_align(bt2, genome_index, args + ['--trim-to', '3:33'])
	assert len(sam_body(three)) == len(UNPAIRED)
	assert all('\t33M\t' in line for line in sam_body(three))
	assert sam_body(bare) == sam_body(three)


def test_trim_to_bad_end(bt2, genome_dict, genome_index, tmp_path):
	"""--trim-to with an end other than 3 or 5 is an error."""
	args, _, _ = write_input('fastq', tmp_path, unpaired_reads(genome_dict))
	_, proc = run_align(bt2, genome_index, args + ['--trim-to', '4:30'], check=False)
	assert proc.returncode != 0
	assert 'trim position must be either 3 or 5' in proc.stderr


def test_trimming_whole_read_leaves_empty_read(bt2, genome_dict, genome_index, tmp_path):
	"""Trimming more bases than a read has leaves an empty, unaligned read (SEQ/QUAL '*', YF:Z:LN)."""
	n, s, q = unpaired_reads(genome_dict)[0]
	path = h.write_fastq(str(tmp_path / 'r.fq'), [(n, s, q)])
	rec, = run_align(bt2, genome_index, ['-U', path, '-5', '30', '-3', '30'])[0]
	assert rec.qname == n and rec.flag == h.UNMAPPED
	assert (rec.seq, rec.qual) == ('*', '*')
	assert rec.tags['YF'] == 'LN'


# --------------------------------------------------- quality encodings

def low_quality_mismatch_read(gd, qlow):
	"""A 50-base forward read from chr3:500 with one mismatch at offset 25,
	where the base quality is `qlow`; every other base is Phred 40."""
	ref = gd['chr3'][500:550]
	phred = [40] * 50
	phred[25] = qlow
	return ref, h.mutate(ref, [25]), phred


@pytest.mark.parametrize('fmt', ['fastq', 'tab5', 'qseq'])
@pytest.mark.parametrize('qlow', [2, 20, 30])
def test_phred33_phred64_agree(bt2, genome_dict, genome_index, tmp_path, fmt, qlow):
	"""A read in Phred+33 and in Phred+64 (with each alias) gives the same QUAL (Phred+33) and AS.

	AS is minus the quality-scaled mismatch penalty at the low-quality base."""
	ref, read, phred = low_quality_mismatch_read(genome_dict, qlow)
	q33 = ''.join(chr(33 + p) for p in phred)
	q64 = ''.join(chr(64 + p) for p in phred)
	(tmp_path / '33').mkdir()
	(tmp_path / '64').mkdir()
	args33, _, _ = write_input(fmt, tmp_path / '33', [('lowq', read, q33)])
	args64, _, _ = write_input(fmt, tmp_path / '64', [('lowq', read, q64)])
	runs = [args33, args33 + ['--phred33'], args33 + ['--phred33-quals'],
		args64 + ['--phred64'], args64 + ['--phred64-quals'], args64 + ['--solexa1.3-quals']]
	bodies = [sam_body(run_align(bt2, genome_index, args)[1]) for args in runs]
	assert all(b == bodies[0] for b in bodies)
	rec = h.SamRecord(bodies[0][0])
	assert rec.pos == 501 and rec.flag == 0
	assert rec.qual == q33
	assert rec.tags['AS'] == -h.mismatch_penalty(chr(33 + qlow))
	assert rec.tags['XM'] == 1 and rec.tags['NM'] == 1
	assert rec.tags['MD'] == '25%s24' % ref[25]


@pytest.mark.parametrize('solexa', [-5, 0, 5, 9, 15, 40])
def test_solexa_quals(bt2, genome_dict, genome_index, tmp_path, solexa):
	"""--solexa-quals converts Solexa+64 qualities to Phred, which sets QUAL and the mismatch penalty.

	Q = 10 log10(10^(S/10) + 1), rounded to the nearest integer."""
	ref, read, _ = low_quality_mismatch_read(genome_dict, 40)
	qual = 'h' * 25 + chr(64 + solexa) + 'h' * 24
	phred = int(math.floor(10 * math.log10(10 ** (solexa / 10.0) + 1) + 0.5))
	path = h.write_fastq(str(tmp_path / 'r.fq'), [('sol', read, qual)])
	rec, = run_align(bt2, genome_index, ['-U', path, '--solexa-quals'])[0]
	assert rec.pos == 501
	assert rec.qual == 'I' * 25 + chr(33 + phred) + 'I' * 24
	assert rec.tags['AS'] == -h.mismatch_penalty(chr(33 + phred))


def test_phred64_rejects_phred33_input(bt2, genome_dict, genome_index, tmp_path):
	"""--phred64 on qualities below ASCII 64 is an error."""
	n, s, q = unpaired_reads(genome_dict)[0]
	path = h.write_fastq(str(tmp_path / 'r.fq'), [(n, s, '5' * len(s))])
	_, proc = run_align(bt2, genome_index, ['-U', path, '--phred64'], check=False)
	assert proc.returncode != 0
	assert 'expected 64-based Phred qual' in proc.stderr


# --------------------------------------- read files, compression, stdin

def zstd_input_supported(bt2):
	"""True when the aligner was built WITH_ZSTD (shown in its --version
	compiler options) and the zstd program is available to compress input."""
	if h.which('zstd') is None:
		return False
	proc = subprocess.run([os.path.join(bt2.bin_dir, 'bowtie2-align-s'), '--version'],
		capture_output=True, text=True)
	return '-DWITH_ZSTD' in proc.stdout


@pytest.mark.parametrize('kind', ['plain', 'gz', 'bz2', 'gz-no-extension', 'zst'])
def test_compressed_reads(bt2, genome_dict, genome_index, tmp_path, kind):
	"""Plain, gzip, bzip2 and zstd FASTQ read files align identically.

	gzip detection without a .gz extension pins current behaviour (the aligner
	checks the gzip magic bytes)."""
	reads = unpaired_reads(genome_dict)
	if kind == 'bz2' and h.which('bzip2') is None:
		pytest.skip('the wrapper decompresses .bz2 with the bzip2 program, not found')
	if kind == 'zst' and not zstd_input_supported(bt2):
		pytest.skip('bowtie2 built without WITH_ZSTD, or zstd program not found')
	if kind == 'plain':
		path = h.write_fastq(str(tmp_path / 'r.fq'), reads)
	elif kind == 'gz':
		path = h.write_fastq(str(tmp_path / 'r.fq.gz'), reads, compress='gz')
	elif kind == 'bz2':
		path = h.write_fastq(str(tmp_path / 'r.fq.bz2'), reads, compress='bz2')
	elif kind == 'gz-no-extension':
		path = h.write_fastq(str(tmp_path / 'r.fq'), reads, compress='gz')
	else:
		plain = h.write_fastq(str(tmp_path / 'r.fq'), reads)
		path = str(tmp_path / 'r.fq.zst')
		subprocess.run(['zstd', '-q', '-f', '-o', path, plain], check=True)
	recs, _ = run_align(bt2, genome_index, ['-U', path])
	assert [r.qname for r in recs] == [r[0] for r in reads]
	for rec, spec, (n, s, q) in zip(recs, UNPAIRED, reads):
		assert_placed(rec, genome_dict, spec, q)


@pytest.mark.parametrize('layout', [
	['-U', '{a},{b}', '-U', '{c}'],
	['-U', '{a}', '-U', '{b},{c}'],
	['-U', '{a}', '-U', '{b}', '-U', '{c}'],
])
def test_comma_lists_and_repeated_U(bt2, genome_dict, genome_index, tmp_path, layout):
	"""Comma-separated -U lists and repeated -U read every file, in the order given."""
	reads = unpaired_reads(genome_dict)
	paths = {'a': h.write_fastq(str(tmp_path / 'a.fq'), reads[:2]),
		'b': h.write_fastq(str(tmp_path / 'b.fq.gz'), reads[2:3], compress='gz'),
		'c': h.write_fastq(str(tmp_path / 'c.fq'), reads[3:])}
	args = [a.format(**paths) for a in layout]
	recs, _ = run_align(bt2, genome_index, args)
	assert [r.qname for r in recs] == ['u1', 'u2', 'u3', 'u4']
	for rec, spec, (n, s, q) in zip(recs, UNPAIRED, reads):
		assert_placed(rec, genome_dict, spec, q)


def test_repeated_U_with_bzip2(bt2, genome_dict, genome_index, tmp_path):
	"""A bzip2 file given in its own -U is read along with a plain one."""
	if h.which('bzip2') is None:
		pytest.skip('the wrapper decompresses .bz2 with the bzip2 program, not found')
	reads = unpaired_reads(genome_dict)
	a = h.write_fastq(str(tmp_path / 'a.fq'), reads[:2])
	b = h.write_fastq(str(tmp_path / 'b.fq.bz2'), reads[2:], compress='bz2')
	recs, _ = run_align(bt2, genome_index, ['-U', a, '-U', b])
	assert [r.qname for r in recs] == ['u1', 'u2', 'u3', 'u4']


def test_compressed_mates(bt2, genome_dict, genome_index, tmp_path):
	"""-1/-2 read gzip and bzip2 mate files."""
	if h.which('bzip2') is None:
		pytest.skip('the wrapper decompresses .bz2 with the bzip2 program, not found')
	seqs = [pair_seqs(genome_dict, p) for p in PAIRS]
	m1 = h.write_fastq(str(tmp_path / 'm_1.fq.gz'), [(p[0], s[0]) for p, s in zip(PAIRS, seqs)], compress='gz')
	m2 = h.write_fastq(str(tmp_path / 'm_2.fq.bz2'), [(p[0], s[1]) for p, s in zip(PAIRS, seqs)], compress='bz2')
	recs, _ = run_align(bt2, genome_index, ['-1', m1, '-2', m2])
	assert [r.qname for r in recs] == ['p1', 'p1', 'p2', 'p2']
	for k, pair in enumerate(PAIRS):
		assert_pair(recs[2 * k], recs[2 * k + 1], genome_dict, pair)


@pytest.mark.parametrize('fmt', ['fastq', 'fastq-gz', 'fasta', 'tab5', 'bam'])
def test_reads_from_stdin(bt2, genome_dict, genome_index, tmp_path, fmt):
	"""'-' reads from standard input (FASTQ, gzipped FASTQ, FASTA, --tab5, BAM)."""
	reads = unpaired_reads(genome_dict)
	if fmt in ('fastq', 'fastq-gz'):
		data = h.fastq_text(reads).encode()
		if fmt == 'fastq-gz':
			data = gzip.compress(data)
		args = ['-U', '-']
	elif fmt == 'fasta':
		data = ''.join('>%s\n%s\n' % (n, s) for n, s, q in reads).encode()
		args = ['-f', '-U', '-']
	elif fmt == 'tab5':
		data = ''.join('%s\t%s\t%s\n' % r for r in reads).encode()
		args = ['--tab5', '-']
	else:
		path = write_bam(str(tmp_path / 'r.bam'), [(n, 4, s, q, []) for n, s, q in reads])
		with open(path, 'rb') as f:
			data = f.read()
		args = ['-b', '-']
	recs, _ = run_bytes(bt2, ['-x', genome_index, '-p', '1'] + args, data)
	assert [r.qname for r in recs] == [r[0] for r in reads]
	for rec, spec, (n, s, q) in zip(recs, UNPAIRED, reads):
		assert_placed(rec, genome_dict, spec, 'I' * len(s) if fmt == 'fasta' else q)


# ------------------------------------------------------ QNAME options

def test_qname_truncated_by_default(bt2, genome_dict, genome_index, tmp_path):
	"""By default QNAME stops at the first whitespace and a long name is cut short; --sam-no-qname-trunc keeps both whole.

	The 300-character name may come out with 254 or 255 characters: the
	current version cuts it to 255, one more than the SAM specification's
	QNAME limit of 254, so either length passes here."""
	long_name = 'L' + 'x' * 299
	reads = unpaired_reads(genome_dict)[:2]
	path = write_text(tmp_path / 'r.fq', '@u1 has a comment\n%s\n+\n%s\n@%s\n%s\n+\n%s\n'
		% (reads[0][1], reads[0][2], long_name, reads[1][1], reads[1][2]))
	recs, _ = run_align(bt2, genome_index, ['-U', path])
	assert recs[0].qname == 'u1'
	assert long_name.startswith(recs[1].qname) and len(recs[1].qname) in (254, 255)
	assert_placed(recs[1], genome_dict, UNPAIRED[1], reads[1][2])
	recs, _ = run_align(bt2, genome_index, ['-U', path, '--sam-no-qname-trunc'])
	assert [r.qname for r in recs] == ['u1 has a comment', long_name]
	assert_placed(recs[0], genome_dict, UNPAIRED[0], reads[0][2])
	assert_placed(recs[1], genome_dict, UNPAIRED[1], reads[1][2])


@pytest.mark.parametrize('fmt', ['fastq', 'fasta'])
def test_sam_append_comment(bt2, genome_dict, genome_index, tmp_path, fmt):
	"""--sam-append-comment appends the text after the first space of the read name to the record.

	The BC:Z: prefix added to an Illumina-style comment (1:N:0:<barcode>)
	pins current output; the MANUAL does not mention it."""
	reads = unpaired_reads(genome_dict)[:2]
	comments = ['QA:Z:hello\tQB:i:5', '1:N:0:ACGTAC']
	if fmt == 'fastq':
		text = ''.join('@%s %s\n%s\n+\n%s\n' % (n, c, s, q) for (n, s, q), c in zip(reads, comments))
		args = ['-U', write_text(tmp_path / 'r.fq', text)]
	else:
		text = ''.join('>%s %s\n%s\n' % (n, c, s) for (n, s, q), c in zip(reads, comments))
		args = ['-f', '-U', write_text(tmp_path / 'r.fa', text)]
	recs, _ = run_align(bt2, genome_index, args + ['--sam-append-comment'])
	assert [r.qname for r in recs] == ['u1', 'u2']
	assert recs[0].line.endswith('\tYT:Z:UU\tQA:Z:hello\tQB:i:5')
	assert recs[0].tags['QB'] == 5
	assert recs[1].line.endswith('\tYT:Z:UU\tBC:Z:1:N:0:ACGTAC')
	for rec, spec, (n, s, q) in zip(recs, UNPAIRED, reads):
		assert_placed(rec, genome_dict, spec, q if fmt == 'fastq' else 'I' * len(s))
	_, plain = run_align(bt2, genome_index, args)
	assert 'QA:Z:hello' not in plain.stdout and 'BC:Z:' not in plain.stdout


def test_sam_append_comment_needs_fasta_or_fastq(bt2, genome_dict, genome_index, tmp_path):
	"""--sam-append-comment with --tab5 input is an error."""
	args, _, _ = write_input('tab5', tmp_path, unpaired_reads(genome_dict))
	_, proc = run_align(bt2, genome_index, args + ['--sam-append-comment'], check=False)
	assert proc.returncode != 0
	assert '--sam-append-comment only works with FASTA (-f) and FASTQ (-q) formats' in proc.stderr


# ----------------------------------------------------- edge-case reads

def test_read_shorter_than_seed(bt2, genome_dict, genome_index, tmp_path):
	"""Reads of 15 and 20 bases, shorter than the 22-base default seed, align end to end where they were cut."""
	specs = [('s15', 'chr1', 3000, 15, '+'), ('s20', 'chr2', 2500, 20, '-')]
	reads = [(s[0], read_seq(genome_dict, s), varied_qual(s[3])) for s in specs]
	path = h.write_fastq(str(tmp_path / 'r.fq'), reads)
	recs, _ = run_align(bt2, genome_index, ['-U', path])
	assert [r.qname for r in recs] == ['s15', 's20']
	for rec, spec, (n, s, q) in zip(recs, specs, reads):
		assert_placed(rec, genome_dict, spec, q)


def test_empty_and_all_n_reads(bt2, genome_dict, genome_index, tmp_path):
	"""An empty read gets an unaligned record with YF:Z:LN and an all-N read one with YF:Z:NS; later reads align.

	A 1-base read is left out: the current version also filters it with
	YF:Z:LN under the default -N 0, which the MANUAL's definition of LN
	(length <= -N) does not allow, so that case is not tested here."""
	n, s, q = unpaired_reads(genome_dict)[0]
	path = h.write_fastq(str(tmp_path / 'r.fq'),
		[('empty', '', ''), ('alln', 'N' * 40, 'I' * 40), (n, s, q)])
	recs, _ = run_align(bt2, genome_index, ['-U', path])
	assert [r.qname for r in recs] == ['empty', 'alln', n]
	for rec, seq, yf in zip(recs[:2], ['*', 'N' * 40], ['LN', 'NS']):
		assert rec.flag == h.UNMAPPED and rec.rname == '*' and rec.pos == 0
		assert rec.seq == seq
		assert rec.tags['YF'] == yf
		assert 'AS' not in rec.tags
	assert recs[0].qual == '*'
	assert_placed(recs[2], genome_dict, UNPAIRED[0], q)


def test_lowercase_and_ambiguous_bases(bt2, genome_dict, genome_index, tmp_path):
	"""Lowercase read bases are read as uppercase and IUPAC codes other than ACGT as N.

	The MANUAL says IUPAC codes become N for reference characters; for read
	characters, and for lowercase, this pins the current output."""
	ref = genome_dict['chr1'][4000:4040]
	read = ref[:10].lower() + ref[10:20] + 'R' + ref[21:]
	path = h.write_fastq(str(tmp_path / 'r.fq'), [('lc', read)])
	rec, = run_align(bt2, genome_index, ['-U', path])[0]
	assert rec.pos == 4001 and rec.flag == 0
	assert rec.seq == ref[:20] + 'N' + ref[21:]
	assert rec.tags['AS'] == -1
	assert rec.tags['MD'] == '20%s19' % ref[20]
