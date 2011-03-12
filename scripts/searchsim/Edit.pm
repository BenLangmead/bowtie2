#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 6/21/2010
#

package Edit;
use strict;
use Carp;

sub new($$$$$) {
	my ($class, $name, $off, $readch, $refch) = @_;
	defined($off) || croak("No offset");
	return bless {
		_name   => $name   || "noname",
		_off    => $off,
		_readch => $readch || croak("No read char"),
		_refch  => $refch  || croak("No reference char")
	}, $class;
}
sub name   { return $_[0]->{_name}   }
sub off    { return $_[0]->{_off}    }
sub readch { return $_[0]->{_readch} }
sub refch  { return $_[0]->{_refch}  }

##
# Return a string version of this Edit.
#
sub toString($) {
	my $self = shift;
	return $self->off().":".$self->refch().">".$self->readch();
}

##
# Return a string version of an array ref of Edits.
#
sub listToString($) {
	my $list = shift;
	my $ret = "";
	if(scalar(@$list) > 0) {
		for my $e (@$list) {
			$ret .= "," if $ret ne "";
			$ret .= $e->toString();
		}
	} else {
		$ret = "-";
	}
	return $ret;
}

sub is_insert($) {
	return $_[0]->{_readch} eq "-";
}

sub is_delete($) {
	return $_[0]->{"_refch"} eq "-";
}

sub is_mismatch($) {
	return !$_[0]->is_insert() && !$_[0]->is_delete();
}

##
# Utility function that takes an array ref of edits and reverses it
# w/r/t an interval of length $len.
#
sub flip($$) {
	my ($elist, $len) = @_;
	@$elist = reverse @$elist;
	for my $e (@$elist) {
		if($e->is_insert()) {
			$e->{_off} = $len - $e->{_off};
		} else {
			$e->{_off} = $len - $e->{_off} - 1;
		}
	}
}

sub test1() {
	my $e = Edit->new("blah", 4, "-", "a");
	$e->is_insert() || die;
	$e->is_delete() && die;
	$e->is_mismatch() && die;
}

sub test2() {
	my $e = Edit->new("blah", 4, "a", "-");
	$e->is_insert() && die;
	$e->is_delete() || die;
	$e->is_mismatch() && die;
}

sub test3() {
	my $e = Edit->new("blah", 4, "a", "c");
	$e->is_insert() && die;
	$e->is_delete() && die;
	$e->is_mismatch() || die;
}

sub test4() {
	my $e = [ Edit->new("blah", 4, "a", "c"),
	          Edit->new("blah", 6, "-", "c"),
	          Edit->new("blah", 8, "a", "-") ];
	Edit::flip($e, 20);
	$e->[0]->off() == 11   || croak();
	$e->[0]->is_delete()   || croak();
	$e->[1]->off() == 14   || croak();
	$e->[1]->is_insert()   || croak();
	$e->[2]->off() == 15   || croak();
	$e->[2]->is_mismatch() || croak();
}

if($0 =~ /Edit\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
	print "Test \"test1\"..."; test1(); print "PASSED\n";
	print "Test \"test2\"..."; test2(); print "PASSED\n";
	print "Test \"test3\"..."; test3(); print "PASSED\n";
	print "Test \"test4\"..."; test4(); print "PASSED\n";
}

1;
