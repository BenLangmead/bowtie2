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

package DNA;
use strict;
use warnings;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use Test;

##
# Set up uppercase IUPAC characters minus N
#
my %iupac_u_nn = (
	"R" => 1,
	"Y" => 1,
	"M" => 1,
	"K" => 1,
	"S" => 1,
	"W" => 1,
	"B" => 1,
	"V" => 1,
	"H" => 1,
	"D" => 1
);

##
# Return true iff arg is an IUPAC code and not ACGT or N.
#
sub isIUPAC($) {
	return defined($iupac_u_nn{$_[0]});
}

##
# Replace IUPAC characters with Ns.
#
sub iupacSubN($) {
	$_[0] =~ tr/RYMKSWBVHDrymkswbvhd/NNNNNNNNNNnnnnnnnnnn/;
	return $_[0];
}

my %compat = (
	"A" => "A",
	"T" => "T",
	"C" => "C",
	"G" => "G",
	"R" => "AG",
	"Y" => "CT",
	"M" => "AC",
	"K" => "GT",
	"S" => "CG",
	"W" => "AT",
	"B" => "CGT",
	"V" => "ACG",
	"H" => "ACT",
	"D" => "AGT",
	"N" => "ACGT"
);

##
# Pick a random character that's compatible with input character.
#
sub pickCompat($) {
	my $c = uc $_[0];
	defined($compat{$c}) || die "Bad input $c";
	my $cc = $compat{$c};
	if(length($cc) == 1) {
		return $cc;
	} else {
		length($cc) > 1 || die;
		return substr($cc, int(rand(length($cc))), 1);
	}
}

my %incompat = (
	"A" => "CGT",
	"T" => "ACG",
	"C" => "AGT",
	"G" => "ACT",
	"R" => "CT",
	"Y" => "AG",
	"M" => "GT",
	"K" => "AC",
	"S" => "AT",
	"W" => "CG",
	"B" => "A",
	"V" => "T",
	"H" => "G",
	"D" => "C",
	"N" => "ACGT"
);

##
# Pick a random character that's compatible with input character.
#
sub pickIncompat($) {
	my $c = uc $_[0];
	defined($incompat{$c}) || die "Bad input $c";
	my $cc = $incompat{$c};
	if(length($cc) == 1) {
		return $cc;
	} else {
		length($cc) > 1 || die;
		return substr($cc, int(rand(length($cc))), 1);
	}
}

##
# Set up lowercase IUPAC characters minus N
#
my %iupac_l_nn = (
	"r" => 1,
	"y" => 1,
	"m" => 1,
	"k" => 1,
	"s" => 1,
	"w" => 1,
	"b" => 1,
	"v" => 1,
	"h" => 1,
	"d" => 1
);

my %revcompMap = (
	"A" => "T",
	"T" => "A",
	"C" => "G",
	"G" => "C",
	"R" => "Y",
	"Y" => "R",
	"M" => "K",
	"K" => "M",
	"S" => "S",
	"W" => "W",
	"B" => "V",
	"V" => "B",
	"H" => "D",
	"D" => "H",
	"N" => "N",
	"a" => "t",
	"t" => "a",
	"c" => "g",
	"g" => "c",
	"r" => "y",
	"y" => "r",
	"m" => "k",
	"k" => "m",
	"s" => "s",
	"w" => "w",
	"b" => "v",
	"v" => "b",
	"h" => "d",
	"d" => "h",
	"n" => "n"
);

my %unambigSet = (
	"A" => 1, "a" => 1,
	"C" => 1, "c" => 1,
	"G" => 1, "g" => 1,
	"T" => 1, "t" => 1
);

##
# Return the complement, incl. if it's IUPAC.
#
sub comp($) {
	my $s = uc shift;
	return $s =~ tr/ACGTRYMKSWBVHDN/TGCAYRKMSWVBDHN/;
}

##
# Return the reverse complement of a string.
#
sub revcomp($) {
	my $s = uc shift;
	$s = reverse $s;
	$s =~ tr/ACGTRYMKSWBVHDN/TGCAYRKMSWVBDHN/;
	return $s;
}

##
# Return true iff it's unambiguous.
#
sub unambig($) {
	return $unambigSet{$_[0]};
}

##
# Manipulate DNA in an integer-indexed fashion.
#
sub plus($$) {
	my ($c, $amt) = @_;
	my %ctoi = ("A" => 0, "C" => 1, "G" => 2, "T" => 3);
	my %itoc = (0 => "A", 1 => "C", 2 => "G", 3 => "T");
	$c = uc $c;
	defined($ctoi{$c}) || die "Not an unambiguous nucleotide: $c";
	return $itoc{($ctoi{$c}+$amt) % 4};
}

my %dinucToColorMap = (
	"AA" => "0",
	"AC" => "1",
	"AG" => "2",
	"AT" => "3",
	"CC" => "0",
	"CG" => "3",
	"CT" => "2",
	"GG" => "0",
	"GT" => "1",
	"TT" => "0",
	
	"AN" => ".",
	"CN" => ".",
	"GN" => ".",
	"NT" => ".",
	"NN" => ".",
);

sub dinucToColor($$) {
	my ($n1, $n2) = @_;
	if(ord($n2) < ord($n1)) {
		my $tmp = $n1;
		$n1 = $n2;
		$n2 = $tmp;
	}
	ord($n1) <= ord($n2) || die;
	defined($dinucToColorMap{"$n1$n2"}) ||
		die "Bad nucleotide dinuc: '$n1$n2'";
	return $dinucToColorMap{"$n1$n2"};
}

sub test1 {
	plus("A", 1) eq "C" || die;
	plus("C", 1) eq "G" || die;
	plus("G", 1) eq "T" || die;
	plus("T", 1) eq "A" || die;
	return 1;
}

sub test2 {
	plus("A", 2) eq "G" || die;
	plus("C", 2) eq "T" || die;
	plus("G", 2) eq "A" || die;
	plus("T", 2) eq "C" || die;
	return 1;
}

sub test3 {
	revcomp("ACGT") eq "ACGT" || die;
	return 1;
}

sub test4 {
	revcomp("ACGTYR") eq "YRACGT" || die;
	return 1;
}

if($0 =~ /[^0-9a-zA-Z_]?DNA\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
	Test::shouldSucceed("test1", \&test1);
	Test::shouldSucceed("test2", \&test2);
	Test::shouldSucceed("test3", \&test3);
	Test::shouldSucceed("test4", \&test4);
}

1;
