#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 6/29/2010
#
# Encapsulates the search results for a given read.  For each valid
# alignment found, we record a map from cost to (a) reference string,
# and (b) edits.
#

package ResultAln;
use strict;
use warnings;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use Edit;

sub new {
	my ($class, $cost, $rstr, $fw, $read, $edits) = @_;
	defined($cost) || croak("No cost");
	return bless {
		_cost  => $cost,
		_rstr  => $rstr  || croak("No rstr"),
		_fw    => $fw    || croak("No fw"),
		_read  => $read  || croak("No read"),
		_edits => $edits || [ ]
	}, $class;
}
sub cost  { return $_[0]->{_cost}  }
sub rstr  { return $_[0]->{_rstr}  }
sub fw    { return $_[0]->{_fw}    }
sub read  { return $_[0]->{_read}  }
sub edits { return $_[0]->{_edits} }

1;
