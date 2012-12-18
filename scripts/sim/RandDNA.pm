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

package RandDNA;
use strict;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use DNA;
use Test;
use Math::Random;

##
# Create a new random DNA generator
#
sub new {
	my (
		$class,
		$name,  # name of generator
		$n,     # N fraction
		$iupac, # Non-A/C/G/T/N IUPAC fraction (after N fraction removed)
		$at,    # AT fraction (after N/IUPAC fractions removed)
		$a,     # A fraction (of AT)
		$c,     # C fraction (of CG)
	) = @_;
	$name = "noname" unless defined($name);
	return bless {
		_name   => $name,
		_n      => $n     || croak("No N frac"),
		_iupac  => $iupac || croak("No IUPAC frac"),
		_at     => $at    || rand(),
		_a      => $a     || rand(),
		_c      => $c     || rand()
	}, $class;
}
sub name      { return $_[0]->{_name}  }
sub nFrac     { return $_[0]->{_n}     }
sub iupacFrac { return $_[0]->{_iupac} }
sub atFrac    { return $_[0]->{_at}    }
sub aFrac     { return $_[0]->{_a}     }
sub cFrac     { return $_[0]->{_c}     }

##
# Return a random IUPAC character.
#
sub randIUPAC() {
	my @iu  = (
		"R",
		"Y",
		"M",
		"K",
		"S",
		"W",
		"B",
		"V",
		"H",
		"D"
	);
	my $iui = int(rand(scalar(@iu)));
	defined($iu[$iui]) || die "Bad index $iui";
	return $iu[$iui];
}

##
# Given parameters controlling character frequencies, build a palette
# of short sequences that can be composed to make longer sequences.
#
sub genBuildingBlocks {
	my ($self, $arr, $num, $bbrnd) = @_;
	defined($arr) || die;
	$num = 30 unless defined($num);
	# Random generator for length
	$bbrnd = sub { return int(rand(100))+1 } unless defined($bbrnd);
	for my $i (1..$num) {
		my $seq = "";
		# Generate length
		my $len = $bbrnd->();
		$len > 0 || die "Bad length: $len";
		# Generate characters
		for my $j (1..$len) {
			my $c = "";
			if(rand() < $self->nFrac()) {
				$c = "N";
			} elsif(rand() < $self->iupacFrac()) {
				$c = randIUPAC();
				defined($c) || die;
			} else {
				$c     = (rand() < $self->atFrac()) ? "AT" : "CG";
				if($c eq "AT") {
					$c = (rand() < $self->aFrac())  ? "A"  : "T";
				} else {
					$c = (rand() < $self->cFrac())  ? "C"  : "G";
				}
			}
			$seq .= $c;
		}
		# Add to return list
		push @$arr, $seq;
	}
}

##
# Use this generator to generate another random sequence.
#
sub nextSeq {
	my ($self, $len, $bbsr, $runrnd) = @_;
	my $seq = "";
	# Generate building blocks
	my @bbs = ();
	if(defined($bbsr)) {
		@bbs = @$bbsr;
	} else {
		$self->genBuildingBlocks(\@bbs);
	}
	scalar(@bbs) > 0 || die;
	# Random generator for run length
	my $defaultRnd = sub {
		# Mean of exp is 1/lambda
		return int(Math::Random::random_exponential(1, 2))+1;
	};
	$runrnd = $defaultRnd unless defined($runrnd);
	# Build the sequence by repeatedly inserting runs of building blocks
	while(length($seq) < $len * 5) {
		# Choose building block
		my $bbi = int(rand(scalar(@bbs)));
		# Choose how many times to add it
		my $runlen = $runrnd->();
		# Choose insert point
		my $insat = int(rand(length($seq)));
		# Repeatedly insert building lock
		for my $i (1..$runlen) {
			substr($seq, $insat, 0) = $bbs[$bbi];
		}
	}
	# Return chopped out piece
	return substr($seq, $len * 2, $len);
}

sub test1 {
	my $rd = RandDNA->new(
		"randtest1", # name
		0.02,        # n frac
		0.01,        # IUPAC frac
		0.4,         # AT frac
		0.45,        # A/AT frac
		0.35);       # C/CG frac
	my $seq = $rd->nextSeq(200, undef);
	length($seq) == 200 || die;
	return 1;
}

sub test2 {
	my @bb = ("AAA", "CCCC");
	my $rd = RandDNA->new(
		"randtest2", # name
		0.02,        # n frac
		0.01,        # IUPAC frac
		0.4,         # AT frac
		0.45,        # A/AT frac
		0.35);       # C/CG frac
	my $seq = $rd->nextSeq(300, \@bb);
	length($seq) == 300 || die;
	return 1;
}

if($0 =~ /[^0-9a-zA-Z_]?RandDNA\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
	Test::shouldSucceed("test1", \&test1);
	Test::shouldSucceed("test2", \&test2);
}

1;
