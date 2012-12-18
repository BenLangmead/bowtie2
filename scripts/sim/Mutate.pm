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

package Mutate;
use strict;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use DNA;
use Test;
use List::Util qw(max min);
use Math::Random;

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
sub defaultGapLenGen() {
	return int(Math::Random::random_exponential(1, 3))+1;
}

##
# Default function for generating random sequence to insert into a gap.
#
sub defaultSeqGen($) {
	my $len = shift;
	($len == int($len) && $len > 0) ||
		die "Bad length for sequence generator: $len";
	my $ret = "";
	for (1..$len) {
		$ret .= substr("ACGT", int(rand(4)), 1);
	}
	return $ret;
}

##
# Create a new DNA mutator
#
sub new {
	my (
		$class,
		$name,   # name
		$snp,    # SNP rate
		$rdgap,  # read gap rate
		$rfgap,  # ref gap rate
		$rearr,  # rearrangement rate
		$gaplen, # gap length
		$seqgen, # DNA generator
	) = @_;
	$name = "noname" unless defined($name);
	$snp    = \&defaultSNPGen    unless defined($snp);
	$rdgap  = \&defaultRdGapGen  unless defined($rdgap);
	$rfgap  = \&defaultRfGapGen  unless defined($rfgap);
	$rearr  = \&defaultRearrGen  unless defined($rearr);
	$gaplen = \&defaultGapLenGen unless defined($gaplen);
	$seqgen = \&defaultSeqGen    unless defined($seqgen);
	return bless {
		_name   => $name,
		_snp    => $snp,
		_rdgap  => $rdgap,
		_rfgap  => $rfgap,
		_rearr  => $rearr,
		_gaplen => $gaplen,
		_seqgen => $seqgen,
	}, $class;
}
sub snp    { return $_[0]->{_snp}    }
sub rdgap  { return $_[0]->{_rdgap}  }
sub rfgap  { return $_[0]->{_rfgap}  }
sub rearr  { return $_[0]->{_rearr}  }
sub gaplen { return $_[0]->{_gaplen} }
sub seqgen { return $_[0]->{_seqgen} }

##
# Given a sequence (i.e. a key $srcchr into the reference hash),
# mutate that string.  Note that rearrangement mutations can affect
# more than one sequence at a time.
#
# Returns a list containing counts for:
#
# 1: number of SNPs added
# 2: number of read gaps added
# 3: number of ref gaps added
# 4: number of rearrangements added
#
sub mutateSeq {
	my ($self, $srcchr, $ref) = @_;
	my ($nsnp, $nrfgap, $nrdgap, $nrearr) = (0, 0, 0, 0);
	my $mutseq = $ref->{$srcchr};
	# Calculate # SNPs to add
	my $len = length($mutseq);
	my $snpRate   = $self->snp->();
	my $rfgapRate = $self->rfgap->();
	my $rdgapRate = $self->rdgap->();
	my $rearrRate = $self->rearr->();
	$nsnp   = Math::Random::random_binomial(1, $len, $snpRate);
	$nrfgap = Math::Random::random_binomial(1, $len, $rfgapRate);
	$nrdgap = Math::Random::random_binomial(1, $len, $rdgapRate);
	$nrearr = Math::Random::random_binomial(1, $len, $rearrRate);
	print STDERR "    Introducing $nsnp SNPs, $nrfgap/$nrdgap ref/read gaps, and $nrearr rearrangements\n";
	$nsnp = min($nsnp, $len);
	# Add the SNPs
	for (1..$nsnp) {
		my $off = int(rand($len)); # where to mutate
		my $add = int(rand(3))+1;  # how to mutate
		my $c = substr($mutseq, $off, 1);
		$c eq "A" || $c eq "C" || $c eq "G" || $c eq "T" || $c eq "N" || die "Bad char '$c' in:\n$ref->{$srcchr}";
		substr($mutseq, $off, 1) = DNA::plus(substr($mutseq, $off, 1), $add);
	}
	print STDERR "    Finished SNPs\n";
	# Calculate # ref gaps to add
	for (1..$nrfgap) {
		my $off = int(rand($len));      # where to mutate
		my $gaplen = $self->gaplen->(); # how many gap positions in ref
		# Insert characters into the subject genome
		my $insseq = $self->seqgen->($gaplen);
		substr($mutseq, $off, 0) = $insseq;
		$len = length($mutseq);
	}
	print STDERR "    Finished ref gaps\n";
	# Calculate # read gaps to add
	for (1..$nrdgap) {
		my $off = int(rand($len));      # where to mutate
		my $gaplen = $self->gaplen->(); # how many gap positions in ref
		# Delete characters from subject genome
		substr($mutseq, $off, $gaplen) = "";
		$len = length($mutseq);
	}
	print STDERR "    Finished read gaps\n";
	$ref->{$srcchr} = $mutseq;
	return ($nsnp, $nrfgap, $nrdgap, $nrearr);
	
	my $totlen = 0;
	for (keys %$ref) { $totlen += length($ref->{$_}); }
	# Calculate # rearrangements to add
	for (1..$nrearr) {
		# Pick two loci, at least one on this reference sequence and
		# then cross them over somehow
		my $off     = int(rand($len));
		my @refkeys = keys %$ref;
		my $ochr    = $refkeys[int(rand(scalar(@refkeys)))];
		my $oseq    = $ref->{$ochr};
		my $ooff    = int(rand(length($oseq)));
		my $srcleft = int(rand(2));
		my $dstleft = int(rand(2));
		my $srcrc   = int(rand(2));
		my $dstrc   = int(rand(2));
		# Check that the source and dest don't overlap
		next if $srcchr eq $ochr;
		# Get the sequence to move
		my $presrclen = length($mutseq);
		my $predstlen = length($oseq);
		my $srcseq;
		if($srcleft) {
			$srcseq = substr($mutseq, 0, $off);
		} else {
			$srcseq = substr($mutseq, $off);
		}
		my $dstseq;
		if($dstleft) {
			$dstseq = substr($oseq, 0, $ooff);
		} else {
			$dstseq = substr($oseq, $ooff);
		}
		# Delete the sequence from the source
		length($srcseq) <= length($mutseq) || die;
		length($dstseq) <= length($oseq) || die;
		if($srcleft) {
			substr($mutseq, 0, length($srcseq)) = "";
		} else {
			substr($mutseq, -length($srcseq)) = "";
		}
		if($dstleft) {
			substr($oseq, 0, length($dstseq)) = "";
		} else {
			substr($oseq, -length($dstseq)) = "";
		}
		# Possibly reverse the pieces we broke off
		my $len1 = length($srcseq);
		my $len2 = length($dstseq);
		$srcseq = DNA::revcomp($srcseq) if $srcrc;
		$dstseq = DNA::revcomp($dstseq) if $dstrc;
		length($srcseq) == $len1 || die "$srcseq";
		length($dstseq) == $len2 || die;
		# Mutate the current chromosome
		if($srcleft) {
			$mutseq = $dstseq . $mutseq;
		} else {
			$mutseq = $mutseq . $dstseq;
		}
		# Mutate the other chromosome
		if($dstleft) {
			$oseq = $srcseq . $oseq;
		} else {
			$oseq = $oseq . $srcseq;
		}
		my $postsrclen = length($mutseq);
		my $postdstlen = length($oseq);
		($presrclen + $presrclen) == ($postsrclen + $postsrclen) ||
			die "from $srcchr to $ochr: $presrclen + $presrclen != $postsrclen + $postsrclen";
		$ref->{$srcchr} = $mutseq;
		$ref->{$ochr} = $oseq;
		my $ntotlen = 0;
		for (keys %$ref) { $ntotlen += length($ref->{$_}); }
		$totlen == $ntotlen || die "Total length changed after rearrangements from $srcchr to $ochr ($totlen -> $ntotlen)";
	}
	print STDERR "    Finished rearrangements\n";
	$ref->{$srcchr} = $mutseq;
	return ($nsnp, $nrfgap, $nrdgap, $nrearr);
}

sub test1 {
	my $mut = Mutate->new("UnitTest mutator");
	my %refs = (
		"r1" => "TATGACGGTCGAAACCAGGCGA",
		"r2" => "TATATTTAGTCTCGTCTGGCTGTCTCGGCTGCGCGCGAGTAAAGACCGGCCTGATC");
	$mut->mutateSeq("r1", \%refs);
	$mut->mutateSeq("r2", \%refs);
	return 1;
}

sub test2 {
	my $mut = Mutate->new(
		"UnitTest mutator",
		\&defaultSNPGen,
		\&defaultRdGapGen,
		\&defaultRfGapGen,
		sub { return 0.1 },
		\&defaultGapLenGen,
		\&defaultSeqGen);
	my %refs = (
		"r1" => "TATGACGGTCGAAACCAGGCGA",
		"r2" => "TATATTTAGTCTCGTCTGGCTGTCTCGGCTGCGCGCGAGTAAAGACCGGCCTGATC",
		"r3" => "TATATTTAGTCTCGTCTGGCTGTCTCGGCTGCGCGCGAGTAAAGACCGGCCTGATC".
				"ATTGGTGTCGCGGCGCGCGTATATATATATATATATAGCCTGCTACGTCAGCTAGC",
		"r4" => "TATATTTAGTCTCGTCTGGCTGTCTCGGCTGCGCGCGAGTAAAGACCGGCCTGATC".
				"ATTGGTGTCGCGGCGCGCGTATATATATATATATATAGCCTGCTACGTCAGCTAGC".
				"ATATAACAAAAAAACCCCACACGACGCGGACTCTAGCACTATCGGACTATCATCGG");
	$mut->mutateSeq("r1", \%refs);
	$mut->mutateSeq("r2", \%refs);
	$mut->mutateSeq("r3", \%refs);
	$mut->mutateSeq("r4", \%refs);
	return 1;
}

if($0 =~ /[^0-9a-zA-Z_]?Mutate\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
	Test::shouldSucceed("test1", \&test1);
	Test::shouldSucceed("test2", \&test2);
}

1;
