#!/usr/bin/perl -w

#
# Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
#
# This file is part of Bowtie 2.
#
# Bowtie 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Bowtie 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
#

package Sim;
use strict;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use DNA;
use Test;
use RandDNA;
use SampleRead;
use Mutate;
use AlignmentCheck;
use Math::Random;
use List::Util qw(max min);
use POSIX;

##
# Replacement for "die" that additionally writes error message to file so that
# run.pl can read it later.
#
sub mydie($) {
	my $fn = ".run.pl.child.$$";
	open(EO, ">$fn") || die "Could not open $fn for writing";
	print EO "$_[0]\n";
	close(EO);
	confess $_[0];
}

# Generates random printable strings of a given length
sub randStr($) {
	my $len = shift;
	my @chars = ('a'..'z', 'A'..'Z', '0'..'9', '_');
	my $str = "";
	foreach (1..$len) {
		$str .= $chars[int(rand(scalar(@chars)))];
	}
	return $str;
}

##
# Default random generator for number of reference per test case.
#
sub defaultRefNumGen() { return int(Math::Random::random_exponential(1, 8))+1; }

##
# Default random generator for reference length.
#
sub defaultRefLenGen() {
	return int(Math::Random::random_exponential(1, 50000))+1;
}

##
# Default random generator for number of reference per test case.
#
sub defaultReadNumGen() {
	return int(Math::Random::random_exponential(1, 10000))+1;
}

##
# Default random generator for read length.
#
sub defaultFragLenGen() {
	return int(Math::Random::random_normal(1, 200, 40))+1;
}

##
# Default random generator for reference length.
#
sub defaultReadLenGen() {
	my $r = int(rand(3));
	if($r == 0) {
		return int(Math::Random::random_exponential(1, 60))+1;
	} elsif($r == 1) {
		return int(Math::Random::random_exponential(1, 20))+1;
	} else {
		return int(Math::Random::random_exponential(1, 150))+1;
	}
}

##
# Default random generator for fraction of reference characters = N.
#
sub defaultNGen() {
	return Math::Random::random_uniform(1, 0, 0.05);
}

##
# Default random generator for fraction of reference characters = an
# ambiguous IUPAC code.
#
sub defaultIupacGen() {
	return Math::Random::random_uniform(1, 0, 0.01);
}

##
# Default random generator for AT/ACGT fraction.
#
sub defaultAtGen() {
	return min(max(Math::Random::random_normal(1, 0.5, 0.18), 0), 1);
}

##
# Default random generator for A/AT fraction.
#
sub defaultAGen() {
	return min(max(Math::Random::random_normal(1, 0.5, 0.18), 0), 1);
}

##
# Default random generator for C/CG fraction.
#
sub defaultCGen() {
	return min(max(Math::Random::random_normal(1, 0.5, 0.18), 0), 1);
}

##
# Default SNP rate generator.  Doesn't generate the SNP per se, just
# the rate.
#
sub defaultSNPGen() {
	return Math::Random::random_uniform(1, 0, 0.05);
}

##
# Default read gap rate generator.  Doesn't generate the gaps or
# lengths, just the rate.
#
sub defaultRdGapGen() {
	return Math::Random::random_uniform(1, 0, 0.005);
}

##
# Default reference gap rate generator.  Doesn't generate the gaps or
# lengths, just the rate.
#
sub defaultRfGapGen() {
	return Math::Random::random_uniform(1, 0, 0.005);
}

##
# Default rearrangement rate generator.
#
sub defaultRearrGen() {
	return Math::Random::random_uniform(1, 0, 0.005);
}

##
# Default function for generating gap lengths when introducing a gap.
#
sub defaultGapLenGen($) {
	return int(Math::Random::random_exponential(1, 3))+1;
}

##
# Default function for generating random sequence to insert into a gap.
#
sub defaultSeqGen($) {
	my $len = shift;
	($len == int($len) && $len > 0) ||
		mydie("Bad length for sequence generator: $len");
	my $ret = "";
	for (1..$len) {
		$ret .= substr("ACGT", int(rand(4)), 1);
	}
	return $ret;
}

##
# Default sequencing miscall rate generator.
#
sub defaultSeqMmGen() {
	return Math::Random::random_uniform(1, 0, 0.1);
}

##
# Create a new test case simulator
#
sub new {
	my (
		$class,
		$name,       # name of simulator
		$rfnumgen,   # number of reference sequences
		$rflengen,   # reference length
		$rdnumgen,   # number of read sequences per run
		$rdlengen,   # read length generator
		$fraglengen, # fragment length generator
		$ngen,       # N fraction
		$iupacgen,   # Non-A/C/G/T/N IUPAC fraction (after N fraction removed)
		$atgen,      # AT fraction (after N/IUPAC fractions removed)
		$agen,       # A fraction (of AT)
		$cgen,       # C fraction (of CG)
		$snpgen,     # SNP rate gen
		$rdgapgen,   # read gap generator
		$rfgapgen,   # ref gap generator
		$rearrgen,   # rearrangement generator
		$gaplengen,  # gap length generator
		$seqgen,     # gap filler sequence generator
		$seqmm,      # sequencing error generator
	) = @_;
	$rfnumgen   = \&defaultRefNumGen  unless defined($rfnumgen);
	$rflengen   = \&defaultRefLenGen  unless defined($rflengen);
	$rdnumgen   = \&defaultReadNumGen unless defined($rdnumgen);
	$rdlengen   = \&defaultReadLenGen unless defined($rdlengen);
	$fraglengen = \&defaultFragLenGen unless defined($fraglengen);
	$ngen       = \&defaultNGen       unless defined($ngen);
	$iupacgen   = \&defaultIupacGen   unless defined($iupacgen);
	$atgen      = \&defaultAtGen      unless defined($atgen);
	$agen       = \&defaultAGen       unless defined($agen);
	$cgen       = \&defaultCGen       unless defined($cgen);
	$snpgen     = \&defaultSNPGen     unless defined($snpgen);
	$rdgapgen   = \&defaultRdGapGen   unless defined($rdgapgen);
	$rfgapgen   = \&defaultRfGapGen   unless defined($rfgapgen);
	$rearrgen   = \&defaultRearrGen   unless defined($rearrgen);
	$gaplengen  = \&defaultGapLenGen  unless defined($gaplengen);
	$seqgen     = \&defaultSeqGen     unless defined($seqgen);
	$seqmm      = \&defaultSeqMmGen   unless defined($seqmm);
	$name = "noname" unless defined($name);
	return bless {
		_name       => $name,
		_rfnumgen   => $rfnumgen,
		_rflengen   => $rflengen,
		_rdnumgen   => $rdnumgen,
		_rdlengen   => $rdlengen,
		_fraglengen => $fraglengen,
		_ngen       => $ngen,
		_iupacgen   => $iupacgen,
		_atgen      => $atgen,
		_agen       => $agen,
		_cgen       => $cgen,
		_snpgen     => $snpgen,
		_rdgapgen   => $rdgapgen,
		_rfgapgen   => $rfgapgen,
		_rearrgen   => $rearrgen,
		_gaplengen  => $gaplengen,
		_seqgen     => $seqgen,
		_seqmm      => $seqmm,
	}, $class;
}
sub rfnumgen   { return $_[0]->{_rfnumgen}   }
sub rflengen   { return $_[0]->{_rflengen}   }
sub rdnumgen   { return $_[0]->{_rdnumgen}   }
sub rdlengen   { return $_[0]->{_rdlengen}   }
sub fraglengen { return $_[0]->{_fraglengen} }
sub ngen       { return $_[0]->{_ngen}       }
sub iupacgen   { return $_[0]->{_iupacgen}   }
sub atgen      { return $_[0]->{_atgen}      }
sub agen       { return $_[0]->{_agen}       }
sub cgen       { return $_[0]->{_cgen}       }
sub snpgen     { return $_[0]->{_snpgen}     }
sub rdgapgen   { return $_[0]->{_rdgapgen}   }
sub rfgapgen   { return $_[0]->{_rfgapgen}   }
sub rearrgen   { return $_[0]->{_rearrgen}   }
sub gaplengen  { return $_[0]->{_gaplengen}  }
sub seqgen     { return $_[0]->{_seqgen}     }
sub seqmm      { return $_[0]->{_seqmm}      }

##
# Generate DNA generator.
#
sub genDNAgen {
	my $self = shift;
	my $nfrac     = $self->ngen->();
	my $iupacfrac = $self->iupacgen->();
	my $atfrac    = $self->atgen->();
	my $afrac     = $self->agen->();
	my $cfrac     = $self->cgen->();
	my $refdnagen = RandDNA->new(
		"Sim.pm gen",
		$nfrac,
		$iupacfrac,
		$atfrac,
		$afrac,
		$cfrac);
	printf STDERR "Created DNA generator\n";
	printf STDERR "  N frac: %0.3f\n", $nfrac;
	printf STDERR "  IUPAC frac: %0.3f\n", $iupacfrac;
	printf STDERR "  AT/ACGT frac: %0.3f\n", $atfrac;
	printf STDERR "  A/AT frac: %0.3f\n", $afrac;
	printf STDERR "  C/CG frac: %0.3f\n", $cfrac;
	return $refdnagen;
}

##
# Generate and print reference sequences to file of given name.  Also,
# install reference sequences into hash ref $ref.  To allow for
# "overhang" (alignment that hang off the end of the reference), we
# actually write out a little bit less than the full reference sequence
# for each sequence.
#
sub genRef {
	my ($self, $ref, $refdnagen, $conf, $tmpfn) = @_;
	# Get a generator for reference length
	my $reflen = $self->rflengen;
	# Generate the number of references
	my $refnum = $self->rfnumgen->();
	$refnum = sqrt($refnum) if $conf->{small};
	$refnum = 1 if $refnum <= 0;
	$refnum = sqrt($refnum) if $conf->{small};
	$refnum = 1 if $refnum <= 0;
	$refnum = ceil($refnum);
	$refnum = $conf->{numrefs} if defined($conf->{numrefs});
	# Open output file
	open (FA, ">$tmpfn") ||
		mydie("Could not open temporary fasta file '$tmpfn' for writing");
	my %ccnt = ();
	print STDERR "Generating $refnum references\n";
	for (1..$refnum) {
		# Randomly generate length
		my $len = $reflen->();
		$len = sqrt($len) if $conf->{small};
		$len = 1 if $len <= 0;
		$len = ceil($len);
		my $seq = $refdnagen->nextSeq($len);
		length($seq) >= $len || die;
		my $name = "Sim.pm.$_";
		$ref->{$name} = $seq;
		# Select amount to trim from LHS
		my $trimleft = int(Math::Random::random_exponential(1, 200));
		# Select amount to trim from RHS
		my $trimright = int(Math::Random::random_exponential(1, 200));
		# Make sure we're leaving some sequence after trimming
		while($trimleft + $trimright > $len) {
			if(int(rand(2))) {
				$trimleft = int($trimleft*0.5);
			} else {
				$trimright = int($trimright*0.5);
			}
		}
		# Trim the sequence
		substr($seq, 0, $trimleft) = "";
		$seq = substr($seq, 0, length($seq)-$trimright);
		my $trimlen = length($seq);
		$trimlen == $len - $trimleft - $trimright ||
			mydie("Unexpected trim combo: $len, $trimleft, $trimright, $trimlen");
		print STDERR "  Generated reference '$name' of untrimmed length $len, trimmed length $trimlen\n";
		print FA ">$name\n";
		my $buf = "";
		length($seq) >= $trimlen || die;
		for my $i (1..$trimlen) {
			my $c = substr($seq, $i-1, 1);
			defined($c) || die;
			$ccnt{$c}++;
			$buf .= $c;
			$ref->{$name} .= $c;
			if($i % 60 == 0) {
				print FA "$buf\n";
				$buf = "";
			}
		}
		print FA "$buf\n" if $buf ne "";
	}
	close(FA);
	print STDERR "Wrote references to $tmpfn\n";
	for my $k (sort keys %ccnt) {
		print STDERR "  $k: $ccnt{$k}\n";
	}
}

##
# Generate a hash of key/value arguments to pass to bowtie2.
#
sub genBuildArgs {
	my ($self, $large_index) = @_;
	my %args = ();
	my $r1 = int(rand(3));
	if($r1 == 0) {
		$args{"--bmaxdivn"} = int(Math::Random::random_exponential(1, 4))+1;
	} elsif($r1 == 1) {
		$args{"--bmax"} = int(Math::Random::random_exponential(1, 10000))+100;
	}
	my $r2 = int(rand(2));
	if($r2 == 0) {
		$args{"--dcv"} = 2 ** (int(rand(10))+4);
	}
	my $r3 = int(rand(5));
	if($r3 == 0) {
		$args{"--packed"} = "";
	}
	my $r4 = int(rand(3));
	if($r4 == 0) {
		$args{"--offrate"} = int(rand(8))+1;
	}
	$args{"--large-index"} = "" if $large_index;
	return \%args;
}

##
# Given a fasta filename, an index basename, and a path to the
# bowtie2-build executable, build nucleotide-space and colorpace
# indexes for the sequences in the fasta file.
#
sub build {
	my ($self, $fa, $idx, $conf, $args) = @_;
	my $argstr = "";
	for (keys %$args) {
		$argstr .= " $_";
		if($args->{$_} ne "") {
			$argstr .= " ".$args->{$_};
		}
	}
	$argstr .= " --sanity";
	# Build nucleotide index
	my $cmd = "$conf->{bowtie2_build_debug} $argstr $fa $idx";
	print STDERR "$cmd\n";
	system($cmd);
	$? == 0 || mydie("Error running '$cmd'; exitlevel=$?");
	print STDERR "Built nucleotide index '$idx'\n";
	# Build colorspace index
	unless($conf->{no_color}) {
		$cmd = "$conf->{bowtie2_build_debug} $argstr -C $fa ${idx}.c";
		print STDERR "$cmd\n";
		system($cmd);
		$? == 0 || mydie("Error running '$cmd'; exitlevel=$?");
		print STDERR "Built colorspace index '$idx'\n";
	}
}

##
# Given a hash of sequences, flatten all IUPAC codes into unambiguous
# nucleotides.
#
sub flattenIUPAC() {
	my ($self, $h) = @_;
	for my $c (keys %$h) {
		my $len = length($h->{$c});
		for my $i (0..$len-1) {
			my $ch = uc substr($h->{$c}, $i, 1);
			my $nc = $ch;
			if(DNA::isIUPAC($ch) || $ch eq "N") {
				if(rand() < $self->snpgen->()) {
					$nc = DNA::pickIncompat($ch);
					defined($nc) || mydie("Couldn't find incompatible base for $ch");
				} else {
					$nc = DNA::pickCompat($ch);
					defined($nc) || mydie("Couldn't find compatible base for $ch");
				}
			}
			if($ch ne $nc) {
				substr($h->{$c}, $i, 1) = $nc;
			}
		}
	}
}

##
# Mutate reference genome into a subject genome.
#
sub mutate() {
	my ($self, $refs) = @_;
	my %subj = %$refs;
	$self->flattenIUPAC(\%subj);
	print STDERR "Flattened IUPAC characters\n";
	my $mutator = Mutate->new(
		"Sim.pm mutator",
		$self->snpgen,
		$self->rdgapgen,
		$self->rfgapgen,
		$self->rearrgen,
		$self->gaplengen,
		$self->seqgen);
	my ($nsnp, $nrfgap, $nrdgap, $nrearr) = (0, 0, 0, 0);
	for(keys %subj) {
		print STDERR "  Mutating sequence $_\n";
		my ($nsnp_, $nrfgap_, $nrdgap_, $nrearr_) = $mutator->mutateSeq($_, \%subj);
		$nsnp   += $nsnp_;
		$nrfgap += $nrfgap_;
		$nrdgap += $nrdgap_;
		$nrearr += $nrearr_;
	}
	print STDERR "Mutated reference genome to subject genome\n";
	print STDERR "  SNPs introduced: $nsnp\n";
	print STDERR "  Reference gaps introduced: $nrfgap\n";
	print STDERR "  Read gaps introduced: $nrdgap\n";
	print STDERR "  Rearrangements introduced: $nrearr\n";
	return \%subj;
}

sub dumpFastq {
	my ($self, $input, $fh1, $fh2) = @_;
	for (1..scalar(@{$input->{seq1s}})) {
		my $seq1 = $input->{seq1s}->[$_-1];
		my $qual1 = $input->{qual1s}->[$_-1];
		print {$fh1} "\@$_\n";
		print {$fh1} "$seq1\n";
		print {$fh1} "+$_\n";
		print {$fh1} "$qual1\n";
		if($input->{paired}) {
			my $seq2 = $input->{seq2s}->[$_-1];
			my $qual2 = $input->{qual2s}->[$_-1];
			print {$fh2} "\@$_\n";
			print {$fh2} "$seq2\n";
			print {$fh2} "+$_\n";
			print {$fh2} "$qual2\n";
		}
	}
}

sub dumpQseq {
	my ($self, $input, $fh1, $fh2) = @_;
	for (1..scalar(@{$input->{seq1s}})) {
		my $seq1 = $input->{seq1s}->[$_-1];
		my $qual1 = $input->{qual1s}->[$_-1];
		print {$fh1} "R\t1\t1\t1\t$_\t$_\t1\t1\t$seq1\t$qual1\t1\n";
		if($input->{paired}) {
			my $seq2 = $input->{seq2s}->[$_-1];
			my $qual2 = $input->{qual2s}->[$_-1];
			print {$fh2} "R\t1\t1\t1\t$_\t$_\t1\t1\t$seq2\t$qual2\t1\n";
		}
	}
}

sub dumpFasta {
	my ($self, $input, $fh1, $fh2) = @_;
	for (1..scalar(@{$input->{seq1s}})) {
		my $seq1 = $input->{seq1s}->[$_-1];
		print {$fh1} ">$_\n";
		print {$fh1} "$seq1\n";
		if($input->{paired}) {
			my $seq2 = $input->{seq2s}->[$_-1];
			print {$fh2} ">$_\n";
			print {$fh2} "$seq2\n";
		}
	}
}

sub dumpRaw {
	my ($self, $input, $fh1, $fh2) = @_;
	for (1..scalar(@{$input->{seq1s}})) {
		my $seq1 = $input->{seq1s}->[$_-1];
		print {$fh1} "$seq1\n";
		if($input->{paired}) {
			my $seq2 = $input->{seq2s}->[$_-1];
			print {$fh2} "$seq2\n";
		}
	}
}

##
# Generate the input (reads plus paired/fragment information)
#
sub genInput {
	my ($self, $refs, $conf) = @_;
	# Select whether we're doing colorspace
	my $color = int(rand(2));
	$color = 0 if $conf->{no_color};
	# Select whether we're doing unpaired or paired-end.
	my $paired = int(rand(2));
	$paired = 0 if $conf->{no_paired};
	# Select format for read file
	my @formats    = ("fastq",   "qseq", "fasta", "raw");
	my @format_arg = (   "-q", "--qseq",    "-f",  "-r");
	my $formati = int(rand(scalar(@formats)));
	my $format     = $formats[$formati];
	my $format_arg = $format_arg[$formati];
	my $tmprdfn1 = "$conf->{tempdir}/Sim.pm.$conf->{randstr}_1.$format";
	my $tmprdfn2 = "$conf->{tempdir}/Sim.pm.$conf->{randstr}_2.$format";
	# Generate reads from the subject genome; no sequencing error yet
	my %input = (
		seq1s      => [],
		seq2s      => [],
		qual1s     => [],
		qual2s     => [],
		mate1fw    => 1,
		mate2fw    => 0,
		paired     => $paired,
		color      => $color,
		format     => $format,
		format_arg => $format_arg,
		file1      => $tmprdfn1,
		file2      => $tmprdfn2 );
	my $read_sampler = SampleRead->new(
		"Sim.pm read sampler",
		$self->fraglengen,
		$self->rdlengen,
		$self->rdlengen);
	print STDERR "Created read sampler\n";
	my $numreads = $self->rdnumgen->();
	$numreads = ceil(sqrt($numreads)) if $conf->{small};
	$numreads == int($numreads) || mydie("numreads $numreads not a number");
	my $tmp = int(rand(3));
	if($tmp == 0) {
		$input{mate2fw} = 1;
	} elsif($tmp == 1) {
		$input{mate1fw} = 0;
		$input{mate2fw} = 1;
	}
	print STDERR "Sampling $numreads reads\n";
	ref($refs) eq "HASH" || mydie("Reference input must be hash ref");
	if($paired) { 
		$read_sampler->genReadPairs(
			$numreads,       # number of reads/fragments to generate
			$input{color},   # colorize?
			$refs,           # hash ref holding reference sequences
			$input{mate1fw}, # orientation of mate 1 when fragment comes from Watson strand
			$input{mate2fw}, # orientation of mate 2 when fragment comes from Watson strand
			$input{seq1s},   # put generated mate1 sequences here
			$input{seq2s},   # put generated mate2 sequences here
			$input{qual1s},  # put generated mate1 quality sequences here
			$input{qual2s}); # put generated mate2 quality sequences here
	} else {
		$read_sampler->genReads(
			$numreads,       # number of reads/fragments to generate
			$input{color},   # colorize?
			$refs,           # hash ref holding reference sequences
			$input{seq1s},   # put generated sequences here
			$input{qual1s}); # put generated quality sequences here
	}
	# TODO: with some probability, sort the reads
	print STDERR "Dumping reads to temporary files $tmprdfn1 & $tmprdfn2\n";
	# Dump reads to output file
	my ($fh1, $fh2);
	open($fh1, ">$tmprdfn1") || mydie("Could not open '$tmprdfn1' for writing");
	open($fh2, ">$tmprdfn2") || mydie("Could not open '$tmprdfn2' for writing");
	if($format eq "fastq") {
		$self->dumpFastq(\%input, $fh1, $fh2);
	} elsif($format eq "qseq") {
		$self->dumpQseq(\%input, $fh1, $fh2);
	} elsif($format eq "fasta") {
		$self->dumpFasta(\%input, $fh1, $fh2);
	} elsif($format eq "raw") {
		$self->dumpRaw(\%input, $fh1, $fh2);
	}
	close($fh1);
	close($fh2);
	return \%input;
}

##
# Mutate reads according to sequencing error model.
#
sub mutateSeq {
	my ($self, $input) = @_;
	return $input;
}

##
# Generate a setting for MA (match bonus).
#
sub genPolicyMA($) {
	my $local = shift;
	return undef if $local;
	return Math::Random::random_uniform(1, 1, 40)
}

##
# Generate a setting for MMP (mismatch penalty).
#
sub genPolicyMMP() {
	my $op = substr("CQR", int(rand(3)), 1);
	if($op eq "C") {
		$op .= Math::Random::random_uniform(1, 1, 40);
	}
	return $op;
}

##
# Generate a setting for NP (penalty for a mismatch involving an N).
#
sub genPolicyNP() {
	my $op = substr("CQR", int(rand(3)), 1);
	if($op eq "C") {
		$op .= int(Math::Random::random_exponential(1, 3))+1;
	}
	return $op;
}

##
# Generate a setting for RDG (read gap open and extend penalties).
#
sub genPolicyRDG() {
	my $op = Math::Random::random_uniform(1, 1, 50);
	if(int(rand(2)) == 0) {
		$op .= ",";
		$op .= Math::Random::random_uniform(1, 1, 20);
	}
	return "$op";
}

##
# Generate a setting for RFG (ref gap open and extend penalties).
#
sub genPolicyRFG() {
	my $op = Math::Random::random_uniform(1, 1, 50);
	if(int(rand(2)) == 0) {
		$op .= ",";
		$op .= Math::Random::random_uniform(1, 1, 20);
	}
	return "$op";
}

##
# Generate a setting for MIN (function determining minimum acceptable score).
#
sub genPolicyMIN($) {
	my $local = shift;
	return undef if ($local || int(rand(2)) == 0);
	my $xx = Math::Random::random_uniform(1, 1, 10);
	my $yy = Math::Random::random_uniform(1, 1, 10);
	if(!$local) {
		$xx = -$xx if int(rand(2)) == 0;
		$yy = -$yy;
	}
	return "L,$xx,$yy";
}

##
# Generate a setting for NCEIL (function determining maximum number of Ns
# allowed).
#
sub genPolicyNCEIL() {
	return undef if int(rand(2)) == 0;
	my $xx = Math::Random::random_uniform(1, 0, 1.5);
	my $yy = Math::Random::random_uniform(1, 0, 1.5);
	return "$xx,$yy";
}

##
# Generate a setting for SEED (# mismatches, length, interval).
#
sub genPolicySEED() {
	return undef if int(rand(2)) == 0;
	# Pick a number of mismatches
	my $sd = substr("012", int(rand(2)), 1);
	if(rand() < 0.9) {
		# Length
		$sd .= ",".int(Math::Random::random_uniform(1, 12, 32));
	}
	return $sd;
}

##
# Generate a setting for -D (# DP fails in a row).
#
sub genPolicyFailStreak() {
	return undef if int(rand(2)) == 0;
	return int(Math::Random::random_uniform(1, 2, 50));
}

##
# Generate a setting for -R (# seeding rounds).
#
sub genPolicySeedRounds() {
	return undef if int(rand(2)) == 0;
	return int(Math::Random::random_uniform(1, 1, 5));
}

##
# Generate a setting for IVAL.  Interval between seeds is a function of the
# read length OR sqaure root of read length OR cube root of read length.
#
sub genPolicyIVAL() {
	return undef if int(rand(2)) == 0;
	# Pick a number of mismatches
	my $iv = substr("LSC", int(rand(3)), 1);
	if($iv eq "L") {
		if(rand() < 0.9) {
			# Multiplier
			$iv .= ",".Math::Random::random_uniform(1, 0.0, 0.5);
		}
		if(rand() < 0.3) {
			# Offset
			$iv .= ",".Math::Random::random_uniform(1, 0.0, 4.0);
		}
	} elsif($iv eq "S") {
		if(rand() < 0.9) {
			# Multiplier
			$iv .= ",".Math::Random::random_uniform(1, 0.0, 3.0);
		}
		if(rand() < 0.3) {
			# Offset
			$iv .= ",".Math::Random::random_uniform(1, 0.0, 7.0);
		}
	} elsif($iv eq "C") {
		if(rand() < 0.9) {
			# Multiplier
			$iv .= ",".Math::Random::random_uniform(1, 0.0, 5.0);
		}
		if(rand() < 0.3) {
			# Offset
			$iv .= ",".Math::Random::random_uniform(1, 0.0, 14.0);
		}
	}
	return $iv;
}

##
# Generate a hash of key/value arguments to pass to bowtie2.
#
sub genAlignArgs {
	my ($self, $input, $color, $large_index, $conf) = @_;
	my %args = ();
	my $local = int(rand(2)) == 0;
	$args{"-u"}         = $conf->{maxreads} if defined($conf->{maxreads});
	$args{"--mm"}       = "" if int(rand(2)) == 0;
	$args{"--trim3"}    = int(rand(10)) if int(rand(2)) == 0;
	$args{"--trim5"}    = int(rand(10)) if int(rand(2)) == 0;
	$args{"--nofw"}     = "" if int(rand(4)) == 0;
	$args{"--norc"}     = "" if int(rand(4)) == 0;
	$args{"--col-keepends"} = "" if ($color && int(rand(3)) == 0);
	$args{"--gbar"}     = int(Math::Random::random_exponential(1, 3))+1 if int(rand(4)) == 0;
	$args{"--local"}    = "" if $local;
	my $rep = int(rand(5));
	if($rep == 0) {
		$args{"-a"} = "";
	} elsif($rep == 1) {
		$args{"-k"} = int(Math::Random::random_exponential(1, 3))+2;
	} elsif($rep == 2) {
		$args{"-M"} = int(Math::Random::random_exponential(1, 3))+2;
	}
	$args{"--rdg"} = genPolicyRDG() if rand() < 0.5;
	$args{"--rfg"} = genPolicyRFG() if rand() < 0.5;
	$args{"--score-min"} = genPolicyMIN($local) if rand() < 0.5;
	$args{"--n-ceil"} = genPolicyNCEIL() if rand() < 0.5;
	$args{"-N"} = genPolicySEED() if rand() < 0.5;
	$args{"-D"} = genPolicyFailStreak() if rand() < 0.5;
	$args{"-R"} = genPolicySeedRounds() if rand() < 0.5;
	$args{"--ma"} = genPolicyMA($local) if rand() < 0.5;
	$args{"--mp"} = genPolicyMMP() if rand() < 0.5;
	$args{"--np"} = genPolicyNP() if rand() < 0.5;
	$args{"-i"} = genPolicyIVAL() if rand() < 0.5;
	$args{"--large-index"} = "" if $large_index;
	$args{"--cp-min"} = int(Math::Random::random_exponential(1, 3)) + 2;
	$args{"--cp-ival"} = int(Math::Random::random_exponential(1, 1)) + 1;
	return \%args;
}

##
# Align the given input set against the given index using the given
# bowtie2 binary and arguments.  Sanity-check the SAM output.
#
sub align {
	my ($self, $fa, $idx, $input, $conf, $args) = @_;
	my $argstr = "";
	for (keys %$args) {
		if(defined($args->{$_})) {
			$argstr .= " $_";
			if($args->{$_} ne "") {
				$argstr .= " ".$args->{$_};
			}
		}
	}
	$argstr .= " -C" if $input->{color};
	$argstr .= " ".$input->{format_arg};
	$idx .= ".c" if $input->{color};
	my $inputfn;
	if($input->{paired}) {
		$inputfn = "-1 $input->{file1} -2 $input->{file2}";
	} else {
		$inputfn = $input->{file1};
	}
	# Create object that will help us sanity-check alignments
	my $ac = AlignmentCheck->new(
		"Sim.pm alignment checker", # name
		[ $fa ],                    # fasta
		"sam",                      # SAM-formatted alignments
		0,                          # no bis-C
		0                           # no bis-CpG
	);
	$ac->nrefs() > 0 || mydie("No references");
	# Run normal (non-debug) Bowtie
	defined($conf->{tempdir}) || mydie("No tmp dir");
	my $als       = "$conf->{tempdir}/Sim.pm.$conf->{randstr}.als";
	my $als_debug = "$conf->{tempdir}/Sim.pm.$conf->{randstr}.debug.als";
	my $als_px    = "$conf->{tempdir}/Sim.pm.$conf->{randstr}.px.als";
	my $als_px_reord = "$conf->{tempdir}/Sim.pm.$conf->{randstr}.px.reord.als";
	my $cmd = "$conf->{bowtie2_debug} $argstr -x $idx $inputfn";
	print "$cmd\n";
	open(ALSDEB, ">$als_debug") || mydie("Could not open '$als_debug' for writing");
	open(ALSDEBCMD, "$cmd |") || mydie("Could not open pipe '$cmd |'");
	my $ival = 50;
	my $nals = 0;
	my @lines = ();
	while(<ALSDEBCMD>) {
		# Remove @PG line because CL: tag can legitimately differ
		print ALSDEB $_ unless /^\@PG/;
		push @lines, $_;
		$nals++;
		print STDERR "  Read $nals alignments...\n" if ($nals % $ival) == 0;
	}
	close(ALSDEBCMD);
	$ac->checkAlignments(\@lines, 0);
	$? == 0 || mydie("bowtie2-align-debug exited with exitlevel $?:\n$cmd\n");
	close(ALSDEB);
	$ac->printSummary();
	# With some probability, also run debug Bowtie and check that
	# results are identical
	if(int(rand(3)) == 0) {
		print STDERR "ALSO checking that bowtie2 and bowtie2-align-debug match up\n";
		# Remove @PG line because CL: tag can legitimately differ
		$cmd = "$conf->{bowtie2} $argstr -x $idx $inputfn | grep -v '^\@PG' > $als";
		print "$cmd\n";
		system($cmd);
		$? == 0 ||
			mydie("Command '$cmd' failed with exitlevel $?");
		$cmd = "diff -uw $als $als_debug";
		print "$cmd\n";
		system($cmd);
		$? == 0 ||
			mydie("diff found a difference between bowtie2 and bowtie2-align-debug ".
			      "output for same input (above)\n");
	}
	# With some probability, also run debug Bowtie in -p X mode with  X > 1 and
	# without the --reorder argument and check that results are identical
	if(int(rand(3)) == 0) {
		print STDERR "ALSO checking that bowtie2 and bowtie2 -p X w/ X > 1 match up\n";
		my $p = int(rand(3))+2;
		$cmd = "$conf->{bowtie2} $argstr -p $p -x $idx $inputfn | grep -v '^\@PG' > $als_px";
		print "$cmd\n";
		system($cmd);
		$? == 0 ||
			mydie("Command '$cmd' failed with exitlevel $?");
		# Sort the $als_px and $als_debug files to guarantee that reads and
		# alignments for a given read appear in the same order in both
		$cmd = "sort -k 1,1 -n -k 2,2 -k 3,3 -k 4,4 < $als_px | grep -v '^\@PG' > $als_px.sorted";
		print "$cmd\n";
		system($cmd);
		$? == 0 ||
			mydie("Failed to sort alignment file $als_px\n");
		# Sort the $als_px and $als_debug files to guarantee that reads and
		# alignments for a given read appear in the same order in both
		$cmd = "sort -k 1,1 -n -k 2,2 -k 3,3 -k 4,4 < $als_debug | grep -v '^\@PG' > $als_debug.sorted";
		print "$cmd\n";
		system($cmd);
		$? == 0 ||
			mydie("Failed to sort alignment file $als_debug\n");
		$cmd = "diff -uw $als_debug.sorted $als_px.sorted";
		print "$cmd\n";
		system($cmd);
		$? == 0 ||
			mydie("diff found a difference between bowtie2-align-debug and bowtie2 ".
			      "-p output for same input (above)\n");
	}

	# With some probability, also run debug Bowtie in -p X mode with X > 1 and
	# with the --reorder argument and check that results are identical
	if(int(rand(3)) == 0) {
		print STDERR "ALSO checking that bowtie2 and bowtie2 -p X --reorder w/ X > 1 match up\n";
		my $p = int(rand(3))+2;
		$cmd = "$conf->{bowtie2} $argstr -p $p -x $idx --reorder $inputfn | grep -v '^\@PG' > $als_px_reord";
		print "$cmd\n";
		system($cmd);
		$? == 0 || mydie("Command '$cmd' failed with exitlevel $?");
		$cmd = "diff -uw $als_debug $als_px_reord";
		print "$cmd\n";
		system($cmd);
		$? == 0 ||
			mydie("diff found a difference between bowtie2-align-debug and bowtie2 ".
			      "-p --reorder output for same input (above)\n");
	}
}

##
# Generate a new test case
#
# Possible key/value pairs in $conf hash:
#
# 1. bowtie2_build:       path to bowtie2-build binary
# 2. bowtie2:             path to bowtie2 binary
# 3. bowtie2_build_debug: path to bowtie2-build-debug binary
# 4. bowtie2_debug:       path to bowtie2-debug binary
# 5. tempdir:             temporary directory for reference/reads/index
# 6. no_paired:           defined & non-0 -> don't generate paired-end datasets
# 7. no_color:            defined & non-0 -> don't generate colorspace datasets
# 8. single_thread:       defined & non-0 -> don't use -p X where X > 1
#
sub nextCase {
	my ($self, $conf) = @_;

	$conf->{bowtie2_build}       = "bowtie2-build"         unless defined($conf->{bowtie2_build});
	$conf->{bowtie2}             = "bowtie2-align"         unless defined($conf->{bowtie2});
	$conf->{bowtie2_build_debug} = "bowtie2-build --debug" unless defined($conf->{bowtie2_build_debug});
	$conf->{bowtie2_debug}       = "bowtie2-align --debug" unless defined($conf->{bowtie2_debug});
	$conf->{tempdir}             = "/tmp"                  unless defined($conf->{tempdir});
	srand(time ^ $$);
	$conf->{randstr} = randStr(8);

	print "*** TEST CASE ***\n";
	
	# Build a large index?
	my $large_index = int(rand(2)) == 0;
	
	# Generate the references
	my $refdnagen = $self->genDNAgen();
	# Generate references and write them to a temporary fasta file
	my $tmpfn = "$conf->{tempdir}/Sim.pm.$conf->{randstr}.fa";
	my %refs = ();
	$self->genRef(\%refs, $refdnagen, $conf, $tmpfn);
	# Run bowtie2-build
	my $tmpidxfn = "$conf->{tempdir}/Sim.pm.$conf->{randstr}";
	my $buildArgs = $self->genBuildArgs($large_index);
	$self->build($tmpfn, $tmpidxfn, $conf, $buildArgs);
	my $numruns = 10;
	$numruns *= 10 if $conf->{small}; # Lots of short runs
	# For each batch of reads / bowtie options
	for(1..$numruns) {
		print "*** Run $_ of $numruns\n";
		# Generate mutated version of the reference as our subject genome
		my $subj = $self->mutate(\%refs);
		# Generate all the input, including reads, pairedness,
		# fragment information, whether it's colorspace, etc
		my $input = $self->genInput($subj, $conf);
		# Mutate the input
		my $mutinput = $self->mutateSeq($input);
		# Select Bowtie arguments
		my $args = $self->genAlignArgs($mutinput, $input->{color}, $large_index, $conf);
		$self->align($tmpfn, $tmpidxfn, $mutinput, $conf, $args);
		# Sanity check output.  Possible sanity checks are:
		# 1. Check alignments & edits against reference
		# 2. Compare bowtie2 and bowtie2-debug
		# 3. Compare -p X>1 and -p 1
	}
}

if($0 =~ /Sim\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
}

1;
