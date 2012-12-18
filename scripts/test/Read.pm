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

package Read;
use strict;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use DNA;
use Test;

sub new {
	my ($class, $name, $seq, $qual, $color, $fw, $orig) = @_;
	$name = "noname" unless defined($name);
	return bless {
		_name   => $name,
		_seq    => $seq   || croak("No sequence"),
		_qual   => $qual  || croak("No qualities"),
		_color  => $color || 0,
		_fw     => $fw    || croak("No orientation"),
		_orig   => $orig  || croak("No original read string")
	}, $class;
}
sub name  { return $_[0]->{_name}  }
sub seq   { return $_[0]->{_seq}   }
sub qual  { return $_[0]->{_qual}  }
sub color { return $_[0]->{_color} }
sub fw    { return $_[0]->{_fw}    }
sub orig  { return $_[0]->{_orig}  }
sub len   { return length($_[0]->seq()) }

##
# Obtain a character from the read.
#
sub at {
	my ($self, $off, $ori) = @_;
	my ($c, $q) = "";
	if($ori eq "RtL") {
		$c = uc substr($self->seq(), -$off-1, 1);
		$q = substr($self->qual(), -$off-1, 1);
	} else {
		$c = uc substr($self->seq(), $off, 1);
		$q = substr($self->qual(), $off, 1);
	}
	length($c) == 1 || die;
	return ($c, $q);
}

##
# Load a set of FASTQ reads into the given reads array.
#
sub fromFastq {
	my ($fh, $color, $reads) = @_;
	$reads = [] unless defined($reads);
	while(<$fh>) {
		my $l1 = $_;
		my $l2 = <$fh>; defined($l2) || croak("Name line followed by EOF");
		my $l3 = <$fh>; defined($l3) || croak("Sequence line followed by EOF");
		my $l4 = <$fh>; defined($l4) || croak("Name2 line followed by EOF");
		my $orig = "$l1$l2$l3$l4";
		chomp($l1); chomp($l2); chomp($l3); chomp($l4);
		push @{$reads}, Read->new(substr($l1, 1), $l2, $l4, $color, "FW", $orig);
	}
	return $reads;
}

##
# Load a set of FASTQ reads into the given reads array.
#
sub fromFastqs {
	my ($fqs, $color, $reads) = @_;
	$reads = [] unless defined($reads);
	for my $f (@$fqs) {
		my $fqfh;
		open($fqfh, $f =~ /\.gz$/ ? "gzip -dc $f |" : "$f") || croak("Could not open $f for reading");
		fromFastq($fqfh, $color, $reads);
		close($fqfh);
	}
	return $reads;
}

##
# Load a set of FASTA reads into the given reads array.
#
sub fromFasta {
	my ($fh, $color, $reads) = @_;
	$reads = [] unless defined($reads);
	while(<$fh>) {
		my $l1 = $_;
		my $l2 = <$fh>; defined($l2) || croak("Name line followed by EOF");
		my $orig = "$l1$l2";
		chomp($l1); chomp($l2);
		my $qual = "I" x length($l2);
		push @{$reads}, Read->new(substr($l1, 1), $l2, $qual, $color, "FW", $orig);
	}
	return $reads;
}

##
# Load a set of FASTA reads into the given reads array.
#
sub fromFastas {
	my ($fas, $color, $reads) = @_;
	$reads = [] unless defined($reads);
	for my $f (@$fas) {
		my $fafh;
		open($fafh, $f =~ /\.gz$/ ? "gzip -dc $f |" : "$f") || croak("Could not open $f for reading");
		fromFasta($fafh, $color, $reads);
		close($fafh);
	}
	return $reads;
}

##
# Load a set of FASTA reads into the given reads array.
#
sub fromStrings {
	my ($strs, $color, $reads) = @_;
	$reads = [] unless defined($reads);
	my $idx = 0;
	for my $str (@$strs) {
		my $qual = "I" x length($str);
		push @{$reads}, new Read($idx, $str, $qual, $color, "FW", $str);
		$idx++;
	}
	return $reads;
}

sub test1 {
	my $r = new Read("blah", "TTACGAACCACAACGTATCG", "I"x20, 0, "FW", "?");
	my ($c, $q) = $r->at(0, "LtR");
	($c eq "T" && $q eq "I") || croak("Expected (T, I), got ($c, $q)\n");
	($c, $q) = $r->at(0, "RtL");
	($c eq "G" && $q eq "I") || croak("Expected (G, I), got ($c, $q)\n");
	($c, $q) = $r->at(1, "LtR");
	($c eq "T" && $q eq "I") || croak("Expected (T, I), got ($c, $q)\n");
	($c, $q) = $r->at(1, "RtL");
	($c eq "C" && $q eq "I") || croak("Expected (C, I), got ($c, $q)\n");
	return 1;
}

sub test2 {
	my $rs = fromStrings(["ACGATGCTACG", "TGACGATGCTAG"], 0);
	$rs->[0]->seq()  eq "ACGATGCTACG"  || croak($rs->[0]->seq());
	$rs->[0]->qual() eq "IIIIIIIIIII"  || croak($rs->[0]->qual());
	$rs->[0]->name() eq "0"            || croak($rs->[0]->name());
	$rs->[1]->seq()  eq "TGACGATGCTAG" || croak($rs->[1]->seq());
	$rs->[1]->qual() eq "IIIIIIIIIIII" || croak($rs->[1]->qual());
	$rs->[1]->name() eq "1"            || croak($rs->[1]->name());
	return 1;
}

if($0 =~ /Read\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
	print "Test \"test1\"..."; test1(); print "PASSED\n";
	print "Test \"test2\"..."; test2(); print "PASSED\n";
}

1;
