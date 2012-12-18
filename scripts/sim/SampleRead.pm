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

package SampleRead;
use strict;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use DNA;
use Test;
use List::Util qw(max min);
use Math::Random;

##
# Default sequencing miscall rate generator.
#
sub defaultSeqMmGen() {
	return Math::Random::random_uniform(1, 0, 0.1);
}

##
# Default random generator for read length.
#
sub defaultFragLenGen() {
	return int(Math::Random::random_normal(1, 200, 40))+1;
}

##
# Default random generator for read length.
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
# Create a new read sampler
#
sub new {
	my (
		$class,
		$name,       # name of simulator
		$fraglengen, # paired-end fragment length generator
		$m1lengen,   # random mate1 length generator
		$m2lengen,   # random mate2 length generator
	) = @_;
	$fraglengen = \&defaultFragLenGen unless defined($fraglengen);
	$m1lengen   = \&defaultReadLenGen unless defined($m1lengen);
	$m2lengen   = \&defaultReadLenGen unless defined($m2lengen);
	$name = "noname" unless defined($name);
	return bless {
		_name       => $name,
		_fraglengen => $fraglengen,
		_m1lengen   => $m1lengen,
		_m2lengen   => $m2lengen,
	}, $class;
}
sub name       { return $_[0]->{_name}       }
sub fraglengen { return $_[0]->{_fraglengen} }
sub m1lengen   { return $_[0]->{_m1lengen}   }
sub m2lengen   { return $_[0]->{_m2lengen}   }

##
# Generate a set of reads from a subject genome encoded in a hash ref.
#
sub genReads {
	my (
		$self,
		$num,          # number of reads/fragments to generate
		$color,        # colorize?
		$refs,         # hash ref holding reference sequences
		$seqs,         # put generated read sequences here
		$quals,        # put generated quality sequences here
		$lengen) = @_; # length generator

	ref($refs)  eq "HASH"  || die "Reference input must be hash ref, is ".ref($refs);
	ref($seqs)  eq "ARRAY" || die "seqs input must be array ref, is ".ref($seqs);
	ref($quals) eq "ARRAY" || die "quals input must be array ref, is".ref($quals);
	$lengen = $self->m1lengen() unless defined($lengen);
	my $totreflen = 0;
	my @keys = keys %$refs;
	for (@keys) { $totreflen += length($refs->{$_}); }
	for(1..$num) {
		if(rand() < 0.05 && scalar(@$seqs) > 0) {
			# Clone a previous read
			my $ci = int(rand(scalar(@$seqs)));
			push @$seqs, $seqs->[$ci];
			push @$quals, $quals->[$ci];
		} else {
			while(1) {
				my $ro = int(rand($totreflen));
				my $len = $lengen->();
				$len = 1 if $len < 1;
				my $key = undef;
				my $rflen = 0;
				for (@keys) {
					$rflen = length($refs->{$_});
					if($ro < $rflen) {
						$key = $_;
						last;
					}
					$ro -= $rflen;
				}
				defined($key) || die;
				$rflen > 0 || die;
				# If we are overhanging the end, discard and try again
				next if $ro + $len > $rflen;
				my $rfseq = substr($refs->{$key}, $ro, $len);
				length($rfseq) == $len || die;
				my $rc = int(rand(2));
				# Possibly reverse-complement it
				$rfseq = DNA::revcomp($rfseq) if $rc == 1;
				# Possible colorize
				if($color) {
					my $cseq = "";
					for(0..$len-2) {
						my ($c1, $c2) = (substr($rfseq, $_, 1), substr($rfseq, $_+1, 1));
						my $col = DNA::dinucToColor($c1, $c2);
						$cseq .= $col;
					}
					$rfseq = $cseq;
					$len = length($rfseq);
				}
				push @$seqs, $rfseq;
				# TODO: generate interesting qualities
				push @$quals, "I" x $len;
				last;
			}
		}
		# Simulate next read
	}
}

##
# Generate a set of read pairs from a subject genome encoded in a hash
# ref.  First we extract unpaired fragments, then take sequences from
# either end to make the mates.
#
sub genReadPairs {
	my (
		$self,
		$num,          # number of reads/fragments to generate
		$color,        # colorize?
		$refs,         # hash ref holding reference sequences
		$m1fw,         # orientation of mate 1 when fragment comes from Watson strand
		$m2fw,         # orientation of mate 2 when fragment comes from Watson strand
		$seq1s,        # put generated mate1 sequences here
		$seq2s,        # put generated mate2 sequences here
		$qual1s,       # put generated mate1 quality sequences here
		$qual2s) = @_; # put generated mate2 quality sequences here
	
	# First simulate fragments
	ref($refs)   eq "HASH"  || die "Reference input must be hash ref";
	ref($seq1s)  eq "ARRAY" || die "seq1s input must be array ref";
	ref($seq2s)  eq "ARRAY" || die "seq2s input must be array ref";
	ref($qual1s) eq "ARRAY" || die "qual1s input must be array ref";
	ref($qual2s) eq "ARRAY" || die "qual2s input must be array ref";
	my @fragseqs = ();
	my @fragquals = ();
	$self->genReads(
		$num,
		$color,
		$refs,
		\@fragseqs,
		\@fragquals,
		$self->fraglengen);
	scalar(@fragseqs) == scalar(@fragquals) || die;
	# For each fragment
	for (1..scalar(@fragseqs)) {
		# Take mates from either end
		my $m1len = $self->m1lengen->();
		my $m2len = $self->m2lengen->();
		$m1len = min($m1len, length($fragseqs[$_-1]));
		$m2len = min($m2len, length($fragseqs[$_-1]));
		my $m1seq  = substr($fragseqs [$_-1], 0, $m1len);
		my $m2seq  = substr($fragseqs [$_-1], -$m2len);
		my $m1qual = substr($fragquals[$_-1], 0, $m1len);
		my $m2qual = substr($fragquals[$_-1], -$m2len);
		if(!$m1fw) {
			$m1seq  = DNA::revcomp($m1seq);
			$m1qual = reverse $m1qual;
		}
		if(!$m2fw) {
			$m2seq  = DNA::revcomp($m2seq);
			$m2qual = reverse $m2qual;
		}
		# Commit new pair to the list
		push @$seq1s,  $m1seq;
		push @$seq2s,  $m2seq;
		push @$qual1s, $m1qual;
		push @$qual2s, $m2qual;
		# Simulate next pair
	}
}

sub test1 {
	my $samp = SampleRead->new("UnitTest read sampler");
	my %refs = (
		"r1" => "TATGACGGTCGAAACCAGGCGA",
		"r2" => "TATATTTAGTCTCGTCTGGCTGTCTCGGCTGCGCGCGAGTAAAGACCGGCCTGATC");
	my @seqs = ();
	my @quals = ();
	$samp->genReads(10, 0, \%refs, \@seqs, \@quals, \&defaultReadLenGen);
	scalar(@seqs) == 10 || die;
	scalar(@quals) == 10 || die;
	return 1;
}

sub test2 {
	return 1;
}

if($0 =~ /[^0-9a-zA-Z_]?SampleRead\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
	Test::shouldSucceed("test1", \&test1);
	Test::shouldSucceed("test2", \&test2);
}

1;
