#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 9/5/2009
#
# A SearchScheme is a collection of SearchRoots.
#

package SearchScheme;
use strict;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use Edit;
use Constraint;
use SearchRoot;

sub new {
	my ($class, $name, $roots, $instantiated) = @_;
	return bless {
		_name  => $name  || "noname",
		_roots => $roots || [ ],
		_instantiated => $instantiated || 0
	}, $class;
}
sub name         { return $_[0]->{_name} }
sub roots        { return $_[0]->{_roots} }
sub instantiated { return $_[0]->{_instantiated} }

##
# Parse a policy string of this format:
#
# <seed len>,<seed spacing>,<seed edits>,<edits or mismatches>,<>
#
sub parsePolicy {
	my ($policy, $delim) = @_;
	$delim = "," unless defined($delim);
	$policy =~ s/\s//g;
	my @ps = split(/$delim/, $policy);
	scalar(@ps) == 7 ||
		croak("Bad policy tuple; expected 7 elements, got: \"$policy\"\n");
	my ($len, $rep, $edits, $mmOnly, $inclFw, $inclRc, $allow) = @ps;
	$edits <= 1 || croak("Max edits allowed in a seed is 1; got $edits\n");
	my $c0 = Constraint::zero();
	my $c1 = $mmOnly ? Constraint::allow_mms(1) : Constraint::allow_edits(1);
	my $cons = { "0" => $c0, "1" => $c1 };
	my ($z1, $z2) = ("01", "10");
	my @roots = ();
	if($inclFw) {
		push @roots, new SearchRoot(
			"${policy}_LtR_FW_Off", undef, "FW", 0, $len, $rep, $z1, $cons,
			Constraint::allowance($allow), 0, 0, 0);
		push @roots, new SearchRoot(
			"${policy}_RtL_FW_Off", undef, "FW", 0, $len, $rep, $z2, $cons,
			Constraint::allowance($allow), 0, 0, 0);
	}
	if($inclRc) {
		push @roots, new SearchRoot(
			"${policy}_LtR_RC_Off", undef, "RC", 0, $len, $rep, $z1, $cons,
			Constraint::allowance($allow), 0, 0, 0);
		push @roots, new SearchRoot(
			"${policy}_RtL_RC_Off", undef, "RC", 0, $len, $rep, $z2, $cons,
			Constraint::allowance($allow), 0, 0, 0);
	}
}

##
# Return a SearchScheme consisting of 1 SearchRoot that implements a
# silly non-half-and-half strategy for 1-mm or 1-edit search.  For
# 1-mm, pass 1 as the only argument; otherwise, pass 0.
#
sub scheme_one_simple {
	my ($mm, $inclFw, $inclRc, $off, $len, $expand, $retract, $allowance) = @_;
	$off = 0 unless defined($off);
	$len = 2 unless defined($len);
	$expand = 1 unless defined($expand);
	$retract = 1 unless defined($retract);
	my $c1 = $mm ? Constraint::allow_mms(1) : Constraint::allow_edits(1);
	my $cons = { "1" => $c1 };
	my $z1 = "1";
	my $name = $mm ? "OneMismatchSilly" : "OneEditSilly";
	my @roots = ();
	push @roots, new SearchRoot(
		"${name}_1_FW", undef, "FW", $off, $len, 0, $z1, $cons,
		Constraint::allowance($allowance), $expand, $retract, 0) if $inclFw;
	push @roots, new SearchRoot(
		"${name}_1_RC", undef, "RC", $off, $len, 0, $z1, $cons,
		Constraint::allowance($allowance), $expand, $retract, 0) if $inclRc;
	my $scheme = new SearchScheme($name, \@roots);
	for my $r (@roots) { $r->{_scheme} = $scheme; }
 	return $scheme;
}

##
# Return a SearchScheme consisting of 2 SearchRoots that implement a
# simple half-and-half strategy for 1-mm or 1-edit search.  For 1-mm,
# pass 1 as the only argument; otherwise, pass 0.
#
sub scheme_one {
	my ($mm, $inclFw, $inclRc, $off, $len, $expand, $retract, $allowance) = @_;
	$off = 0 unless defined($off);
	$len = 2 unless defined($len);
	$expand = 1 unless defined($expand);
	$retract = 1 unless defined($retract);
	my $c0 = Constraint::zero();
	my $c1 = $mm ? Constraint::allow_mms(1) : Constraint::allow_edits(1);
	my $cons = { "0" => $c0, "1" => $c1 };
	my $z1 = "01";
	my $z2 = "10";
	my $name = $mm ? "OneMismatchHalfAndHalf" : "OneEditHalfAndHalf";
	my @roots = ();
	if($inclFw) {
		push @roots, new SearchRoot(
			"${name}_LtR_FW", undef, "FW", $off, $len, 0, $z1, $cons,
			Constraint::allowance($allowance), $expand, $retract, 0);
		push @roots, new SearchRoot(
			"${name}_RtL_FW", undef, "FW", $off, $len, 0, $z2, $cons,
			Constraint::allowance($allowance), $expand, $retract, 0);
	}
	if($inclRc) {
		push @roots, new SearchRoot(
			"${name}_LtR_RC", undef, "RC", $off, $len, 0, $z1, $cons,
			Constraint::allowance($allowance), $expand, $retract, 0);
		push @roots, new SearchRoot(
			"${name}_RtL_RC", undef, "RC", $off, $len, 0, $z2, $cons,
			Constraint::allowance($allowance), $expand, $retract, 0);
	}
	my $scheme = new SearchScheme($name, \@roots);
	for my $r (@roots) { $r->{_scheme} = $scheme; }
	return $scheme;
}

##
# Instatiate a SearchScheme with respect to a particular read.
#
sub instantiate($$) {
	my ($self, $read) = @_;
	$self->instantiated() && croak("Tried to instantiate an instantiated SearchScheme");
	# Instantiated roots will be placed in here:
	my @instRoots = ();
	# For each root:
	for my $r (@{$self->roots()}) {
		defined($r->zone_str()) || die;
		push @instRoots, $r->instantiate($read);
	}
	return undef if scalar(@instRoots) == 0;
	return SearchScheme->new($self->name(), \@instRoots, 1);
}

##
# Return an abstract Bowtie-like scheme with the given values for N, L and E.
#
sub maqLike($$$$$$) {
	my ($n, $l, $e, $mm, $inclFw, $inclRc) = @_;
	$n = 2 unless defined($n);
	$l = 28 unless defined($l);
	$e = 70 unless defined($e);
	my $c0 = Constraint::zero();
	my $c1 = $mm ? Constraint::allow_mms($n) : Constraint::allow_edits($n);
	my $cons = { "0" => $c0, "1" => $c1 };
	my $z1 = "01";
	my $z2 = "10";
	my $name = $mm ? "Maq,n=${n},l=${l},$e=${e},ungapped" : "Maq,n=${n},l=${l},$e=${e},gapped";
	my @roots = ();
	if($inclFw) {
		push @roots, new SearchRoot(
			"${name}_LtR_FW", undef, "FW", 0, 2, $z1, $cons, Constraint::none(), 1, 1, 0);
		push @roots, new SearchRoot(
			"${name}_RtL_FW", undef, "FW", 0, 2, $z2, $cons, Constraint::none(), 1, 1, 0);
	}
	if($inclRc) {
		push @roots, new SearchRoot(
			"${name}_LtR_RC", undef, "RC", 0, 2, $z1, $cons, Constraint::none(), 1, 1, 0);
		push @roots, new SearchRoot(
			"${name}_RtL_RC", undef, "RC", 0, 2, $z2, $cons, Constraint::none(), 1, 1, 0);
	}
	my $scheme = new SearchScheme($name, \@roots);
	for my $r (@roots) { $r->{_scheme} = $scheme; }
	return $scheme;
}

##
# Return an abstract BWA-like scheme with the given values for N, L and E.
# Usage:   bwa aln [options] <prefix> <in.fq>
#
# Options: -n NUM    max #diff (int) or missing prob under 0.02 err rate (float) [0.04]
#          -o INT    maximum number or fraction of gap opens [1]
#          -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]
#          -i INT    do not put an indel within INT bp towards the ends [5]
#          -d INT    maximum occurrences for extending a long deletion [10]
#          -l INT    seed length [32]
#          -k INT    maximum differences in the seed [2]
#          -m INT    maximum entries in the queue [2000000]
#          -t INT    number of threads [1]
#          -M INT    mismatch penalty [3]
#          -O INT    gap open penalty [11]
#          -E INT    gap extension penalty [4]
#          -R INT    stop searching when there are >INT equally best hits [30]
#          -q INT    quality threshold for read trimming down to 35bp [0]
#          -c        input sequences are in the color space
#          -L        log-scaled gap penalty for long deletions
#          -N        non-iterative mode: search for all n-difference hits (slooow)
#
#
sub scheme_bwa_like() {
	my ($i, $O, $E, $M) = @_;
}

##
# Check that the SearchScheme is internally consistent.
#
sub repOkTest($) {
	my $self = shift;
	my $inst = $self->instantiated();
	for my $r (@{$self->roots()}) {
		$r->instantiated() == $inst || return 0;
	}
	return 1;
}

sub repOk($) { repOkTest($_[0]) || die }

sub test1() {
	my @s = (
		new SearchRoot("blah", undef, "FW", 0, 28, 0, "01", undef, undef, 1, 1, 0),
		new SearchRoot("blah", undef, "FW", 0, 32, 0, "01", undef, undef, 1, 1, 0),
		new SearchRoot("blah", undef, "FW", 0, 32, 0, "01", undef, undef, 0, 0, 0),
		new SearchRoot("blah", undef, "FW", 0, 28, 0, "01", undef, undef, 0, 0, 0));
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my $sch = SearchScheme->new(undef, \@s);
	my $sch2 = $sch->instantiate($r);
	scalar(@{$sch2->roots()}) == 3 || croak("Expected 3 roots, got ".scalar(@{$sch2->roots()}));
	$sch2->roots()->[0]->off() == 0  || die "off was ".$sch2->roots()->[0]->off()."\n";
	$sch2->roots()->[0]->len() == 30 || die "len was ".$sch2->roots()->[0]->len()."\n";
	$sch2->roots()->[1]->off() == 0  || die "off was ".$sch2->roots()->[1]->off()."\n";
	$sch2->roots()->[1]->len() == 30 || die "len was ".$sch2->roots()->[1]->len()."\n";
	$sch2->roots()->[2]->off() == 0  || die "off was ".$sch2->roots()->[2]->off()."\n";
	$sch2->roots()->[2]->len() == 28 || die "len was ".$sch2->roots()->[2]->len()."\n";
	my $ex_zone_str15 = ("0"x15).("1"x15);
	my $ex_zone_str14 = ("0"x14).("1"x14);
	$sch2->roots()->[0]->zone_str() eq $ex_zone_str15 || die;
	$sch2->roots()->[1]->zone_str() eq $ex_zone_str15 || die;
	$sch2->roots()->[2]->zone_str() eq $ex_zone_str14 || die;
	return 1;
}

sub test2() {
	my $s = scheme_one(1, 1, 1);
	my $r = Read->new("blah", "ACGTA"x4, "IIIII"x4, 0, "FW", "?");
	my $sch = $s->instantiate($r);
	$sch->roots()->[0]->off() == 0 || die;
	$sch->roots()->[1]->off() == 0 || die;
	$sch->roots()->[0]->len() == 20 || die;
	$sch->roots()->[1]->len() == 20 || die;
	$sch->roots()->[0]->zone_str() eq "00000000001111111111" || die;
	$sch->roots()->[1]->zone_str() eq "11111111110000000000" || die;
	return 1;
}

if($0 =~ /SearchScheme\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
	print "Test \"test1\"..."; test1(); print "PASSED\n";
	print "Test \"test2\"..."; test2(); print "PASSED\n";
}

1;
