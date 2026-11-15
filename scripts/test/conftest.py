"""Fixtures for the bowtie2 command-line test suite (scripts/test/test_*.py).

Build first, then run from the repository root:

	make -j4 allall
	python3 -m pytest scripts/test

Set BOWTIE2_BIN_DIR to test binaries that were built somewhere else.
"""

import os
import random

import pytest

import bt2_helpers as h


def pytest_report_header(config):
	return 'bowtie2 binaries: %s' % h.BIN_DIR


@pytest.fixture(scope='session')
def bt2():
	runner = h.Bowtie2()
	for prog in ('bowtie2', 'bowtie2-build', 'bowtie2-inspect', 'bowtie2-align-s',
			'bowtie2-build-s', 'bowtie2-inspect-s'):
		if not runner.has_binary(prog):
			pytest.fail('%s not found in %s: build bowtie2 first (make), or set '
				'BOWTIE2_BIN_DIR' % (prog, runner.bin_dir), pytrace=False)
	return runner


@pytest.fixture(scope='session')
def genome():
	"""Three random chromosomes, fixed by seed, as an ordered list of
	(name, sequence). At these lengths every 20-mer is unique with high
	probability, so a read cut from one place aligns only there."""
	rng = random.Random(20261005)
	return [('chr1', h.random_dna(5000, rng)), ('chr2', h.random_dna(3000, rng)),
		('chr3', h.random_dna(2000, rng))]


@pytest.fixture(scope='session')
def genome_dict(genome):
	return dict(genome)


@pytest.fixture(scope='session')
def genome_fasta(genome, tmp_path_factory):
	path = tmp_path_factory.mktemp('genome') / 'genome.fa'
	return str(h.write_fasta(str(path), genome))


@pytest.fixture(scope='session')
def genome_index(bt2, genome_fasta, tmp_path_factory):
	"""Small (.bt2) index of `genome`."""
	prefix = str(tmp_path_factory.mktemp('genome_index') / 'genome')
	bt2.build_index(genome_fasta, prefix, extra=['--seed', '0'])
	return prefix


@pytest.fixture(scope='session')
def genome_index_large(bt2, genome_fasta, tmp_path_factory):
	"""Large (.bt2l) index of `genome`, built with --large-index."""
	prefix = str(tmp_path_factory.mktemp('genome_index_large') / 'genome')
	bt2.build_index(genome_fasta, prefix, extra=['--large-index', '--seed', '0'])
	return prefix


@pytest.fixture(scope='session')
def lambda_fasta():
	return os.path.join(h.EXAMPLE_DIR, 'reference', 'lambda_virus.fa')


@pytest.fixture(scope='session')
def lambda_index():
	"""The index shipped in example/index."""
	return os.path.join(h.EXAMPLE_DIR, 'index', 'lambda_virus')


@pytest.fixture(scope='session')
def example_reads():
	d = os.path.join(h.EXAMPLE_DIR, 'reads')
	return {'reads_1': os.path.join(d, 'reads_1.fq'), 'reads_2': os.path.join(d, 'reads_2.fq'),
		'longreads': os.path.join(d, 'longreads.fq'),
		'bam': os.path.join(d, 'combined_reads.bam')}
