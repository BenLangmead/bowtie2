#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 6/24/2010
#
# Routines for loading and indexing references.
#

package Reference;
use strict;
use warnings;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use Test;

sub new {
	my ($class, $name, $seq) = @_;
	return bless {
		_name => $name || "noname",
		_seq  => $seq  || croak("No sequence")
	}, $class;
}
sub name { return $_[0]->{_name} }
sub seq  { return $_[0]->{_seq}  }

##
# Read all of the sequences out of a fasta file, storing each in the
# given array ref.
#
sub newFromFasta($$) {
	my ($fasta, $refs) = @_;
	my $name = "";
	open(FA, $fasta) || croak("Couldn't open $fasta\n");
	while(<FA>) {
		chomp;
		if(/^>/) {
			$name = substr($_, 1);
			push @{$refs}, Reference->new($name, "");
		} else {
			$name ne "" || croak("No name for first record in fasta file $fasta\n");
			$refs->[-1]->seq() .= $_;
		}
	}
	close(FA);
}

##
# Take each key/value pair to be a name/seq pair and populate the
# 'refs' array ref with corresponding Reference objects.
#
sub newFromHash($$) {
	my ($hash, $refs) = @_;
	for my $k (keys %{$hash}) {
		push @{$refs}, Reference->new($k, $hash->{$k});
	}
}

##
# Take each key/value pair to be a name/seq pair and populate the
# 'refs' array ref with corresponding Reference objects.
#
sub newFromString($$) {
	my ($str, $refs) = @_;
	my $name = scalar(@{$refs});
	push @{$refs}, Reference->new($name, $str);
}

##
# Take a comma-separated list that may be a mix of filenames, raw
# strings, and name/sequence pairs and try to interpret each correctly.
#
sub newFromList {
	my ($list, $refs) = @_;
	$refs = [] unless defined($refs);
	for my $i (@$list) {
		my $fastaLike = ($i =~ /\.(fna|fasta|fa|fas|mfa)$/);
		newFromFasta ($i, $refs) if $fastaLike;
		newFromString($i, $refs) unless $fastaLike;
	}
	return $refs;
}

#
# The following routines return functions (or closures) that answer the
# question "does the reference contain substring x?".  Some do this
# with the aid of an index, others don't.  Some are guarnteed to return
# the right answer every time, while others are prone to returning
# false positives (they'll sometimes indicate that a substring is
# present even when it's not.)
#

##
# Given a ref to a list of references, return a closure that uses the
# Perl index function to return 1 iff the query is truly a reference
# substring. 
#
sub perlIndexStrset($) {
	my ($refs) = @_;
	return sub {
		return 1 if $_[0] eq "";
		for my $r (@{$refs}) { return 1 if index($r->seq(), $_[0]) != -1 }
		return 0;
	}
}

##
# Given a ref to a list of references, return a closure that uses a
# pre-configured Perl hash to return 1 iff the query is truly a
# reference substring.
#
sub perlHashStrset($) {
	my ($refs) = @_;
	my %hash = ();
	for my $r (@{$refs}) {
		for(my $i = 1; $i <= length($r->seq()); $i++) {
			for(my $j = 0; $j <= length($r->seq())-$i; $j++) {
				my $str = substr($r->seq(), $j, $i);
				$hash{$str} = 1;
			}
		}
	}
	$hash{""} = 1;
	return sub { return defined($hash{$_[0]}); }
}

sub test1 {
	my @refs = ();
	newFromHash({ "seq1" => "ACGT"x5 }, \@refs);
	my $strset1 = Reference::perlIndexStrset(\@refs);
	my $strset2 = Reference::perlIndexStrset(\@refs);
	for my $in ("ACGT", "CGTA", "GTAC", "TACG", "ACGT"x5, "CGTA"x4, "A", "C", "G", "T", "") {
		$strset1->($in) || croak("$in should have been in index\n");
		$strset2->($in) || croak("$in should have been in index\n");
	}
	for my $notin ("AAA", "ACGTT", "B") {
		$strset1->($notin) && croak("$notin should not have been in index\n");
		$strset2->($notin) && croak("$notin should not have been in index\n");
	}
	return 1;
}

if($0 =~ /Reference\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
	Test::shouldSucceed("test1", \&test1);
}

1;
