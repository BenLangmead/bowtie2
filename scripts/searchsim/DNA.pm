#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 6/21/2010
#
#

package DNA;
use strict;
use warnings;
use Carp;

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
	"N" => "N"
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
	my $ret = $revcompMap{$_[0]} || die;
	return $ret;
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
	defined($ctoi{$c}) || die;
	return $itoc{($ctoi{$c}+$amt) % 4};
}

1;
