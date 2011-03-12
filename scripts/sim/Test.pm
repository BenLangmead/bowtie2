#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 6/22/2010
#
# Module for test routines.
#

package Test;
use strict;
use Carp;

sub test($$) {
	my ($name, $f) = @_;
	print "Test \"$name\"...";
	$f->();
	print "PASSED\n";
}

sub shouldSucceed($$) {
	my ($name, $f) = @_;
	print "Test \"$name\"...";
	$f->() || die "";
	print "PASSED\n";
}

sub shouldFail($$) {
	my ($name, $f) = @_;
	print "Test \"$name\"...";
	$f->() && die "";
	print "PASSED\n";
}

1;
