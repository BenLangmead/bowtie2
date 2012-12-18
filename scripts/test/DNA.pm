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

my %revcompMap = (
	"A" => "T", "a" => "t",
	"T" => "A", "t" => "a",
	"C" => "G", "c" => "g",
	"G" => "C", "g" => "c",
	"R" => "Y", "r" => "y",
	"Y" => "R", "y" => "r",
	"M" => "K", "m" => "k",
	"K" => "M", "k" => "m",
	"S" => "S", "s" => "s",
	"W" => "W", "w" => "w",
	"B" => "V", "b" => "v",
	"V" => "B", "v" => "b",
	"H" => "D", "h" => "d",
	"D" => "H", "d" => "h",
	"N" => "N", "n" => "n"
);

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
	"N" => "N"
);

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
	my $ret = $revcompMap{$_[0]} || die "Can't reverse-complement '$_[0]'";
	return $ret;
}

##
# Return the complement, incl. if it's IUPAC.
#
sub revcomp {
	my ($ret, $color) = @_;
	$ret = reverse $ret;
	unless($color) {
		for(my $i = 0; $i < length($ret); $i++) {
			substr($ret, $i, 1) = comp(substr($ret, $i, 1));
		}
	}
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
