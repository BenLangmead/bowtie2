#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 9/5/2009
#
# A constraint for a search zone.  Constraints have components
# describing:
#  1. Number of edits (of any kind) allowed
#  2. Number of mismatches allowed
#  3. Number of insertions allowed
#  4. Number of deletions allowed
#

package Constraint;
use strict;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use Test;

sub new {
	my ($class, $edits, $mms, $ins, $dels, $cost, $allow) = @_;
	scalar(@_) == 7 || croak("Expected 7 args, got ".scalar(@_)."\n");
	$edits = 2**30 unless defined($edits);
	$mms   = 2**30 unless defined($mms);
	$ins   = 2**30 unless defined($ins);
	$dels  = 2**30 unless defined($dels);
	$cost  = 2**30 unless defined($cost);
	return bless {
		_edits => $edits,
		_mms   => $mms,
		_ins   => $ins,
		_dels  => $dels,
		_cost  => $cost,
		_allow => $allow
	}, $class;
}
sub edits { return $_[0]->{_edits} }
sub mms   { return $_[0]->{_mms}   }
sub ins   { return $_[0]->{_ins}   }
sub dels  { return $_[0]->{_dels}  }
sub cost  { return $_[0]->{_cost}  }
sub allow { return $_[0]->{_allow} }
sub costDiff { return abs($_[1]->{_cost} - $_[0]->{_cost}) }

##
# Instantiate a Constraint with respect to a read.
#
sub instantiate($$) {
	my ($self, $read) = @_;
	return $self unless defined($self->allow);
	return new Constraint(
		$self->edits,
		$self->mms,
		$self->ins,
		$self->dels,
		int($self->allow * $read->len + 0.5),
		undef
	);
}

##
# Make a copy of this Constraint and return it.
#
sub clone($) {
	my $self = shift;
	my $copy = { %$self };
	bless $copy, ref $self;
}

##
# Return true iff the given qdepth is 
#
sub mustMatch($) {
	my ($self) = @_;
	return $self->edits() == 0 ||
	       ($self->mms() == 0 && $self->dels() == 0 && $self->ins() == 0) ||
	       $self->cost() == 0;
}

##
# Return true if the constraint doesn't rule out the possibility of
# another mismatch.
#
sub canMismatch($) {
	return !$_[0]->mustMatch() && $_[0]->mms() > 0;
}

##
# Return true if the constraint doesn't rule out the possibility of
# another insertion.
#
sub canInsert($) {
	return !$_[0]->mustMatch() && $_[0]->ins() > 0;
}

##
# Return true if the constraint doesn't rule out the possibility of
# another deletion.
#
sub canDelete($) {
	return !$_[0]->mustMatch() && $_[0]->dels() > 0;
}

##
# Return true iff the constraint doesn't rule out the possibility of
# another gap (insertion or deletion).
#
sub canGap($) {
	return $_[0]->canInsert() || $_[0]->canDelete();
}

##
# Modify the Constraint to reflect an additional mismatch.
#
sub chargeMismatch($$) {
	my ($self, $cost) = @_;
	$self->canMismatch() || croak("Bad call to chargeMismatch");
	$self->{_edits}--;
	$self->{_mms}--;
	$self->{_cost} -= $cost;
}

##
# Modify the Constraint to reflect an additional insertion.
#
sub chargeInsert($$) {
	my ($self, $cost) = @_;
	$self->canInsert() || croak("Bad call to chargeInsert");
	$self->{_edits}--;
	$self->{_ins}--;
	$self->{_cost} -= $cost;
}

##
# Modify the Constraint to reflect an additional mismatch.
#
sub chargeDelete($$) {
	my ($self, $cost) = @_;
	$self->canDelete() || croak("Bad call to chargeDelete");
	$self->{_edits}--;
	$self->{_dels}--;
	$self->{_cost} -= $cost;
}

##
# Generate a constraint that allows anything.
#
sub none() {
	return new Constraint(undef, undef, undef, undef, undef, undef);
}

##
# Generate a constraint that allows anything up to some per-base
# allowance, which must be instantiated on a per-read basis.
#
sub allowance($) {
	return new Constraint(undef, undef, undef, undef, undef, $_[0]);
}

##
# Generate a constraint that forbids any edits whatsoever.
#
sub zero() {
	return new Constraint(0, undef, undef, undef, undef, undef);
}

##
# Generate a constraint that constrains to a given number of
# mismatches, but constrains nothing else.
#
sub allow_mms($) {
	my $mms = shift;
	return new Constraint(undef, $mms, 0, 0, undef, undef);
}

##
# Generate a constraint that constrains to a given number of
# edits, but constrains nothing else.
#
sub allow_edits($) {
	my $edits = shift;
	return new Constraint($edits, undef, undef, undef, undef, undef);
}

##
# Generate a constraint that constraints score to the given ceiling,
# and constrains nothing else.
#
sub score($) {
	return Constraint->new(1, undef, undef, undef, $_[0], undef);
}

sub test1 {
	my $c = Constraint::zero();
	$c->mustMatch() || die;
	$c->canMismatch() && die;
	$c->canDelete() && die;
	$c->canInsert() && die;
	$c->canGap() && die;
	return 1;
}

sub test2 {
	my $c = Constraint::allow_mms(1);
	$c->mustMatch() && die;
	$c->canMismatch() || die;
	$c->canDelete() && die;
	$c->canInsert() && die;
	$c->canGap() && die;
	$c->chargeMismatch(10);
	$c->mustMatch() || die;
	$c->canMismatch() && die;
	$c->canDelete() && die;
	$c->canInsert() && die;
	$c->canGap() && die;
	return 1;
}

sub test3 {
	my $c = Constraint::allow_edits(1);
	$c->mustMatch() && die;
	$c->canMismatch() || die;
	$c->canDelete() || die;
	$c->canInsert() || die;
	$c->canGap() || die;
	$c->chargeInsert(10);
	$c->mustMatch() || die;
	$c->canMismatch() && die;
	$c->canDelete() && die;
	$c->canInsert() && die;
	$c->canGap() && die;
	return 1;
}

sub test4 {
	use Read;
	my $c  = Constraint::allowance(1);
	my $r  = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my $cc = $c->instantiate($r);
	$cc->cost == 30 || die;
	defined($cc->allow) && die;
	return 1;
}

if($0 =~ /Constraint\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
	print "Test \"test1\"..."; test1(); print "PASSED\n";
	print "Test \"test2\"..."; test2(); print "PASSED\n";
	print "Test \"test3\"..."; test3(); print "PASSED\n";
	print "Test \"test4\"..."; test4(); print "PASSED\n";
}

1;
