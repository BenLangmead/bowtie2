#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 6/30/2010
#
# Encapsulates output formats.
#

package Formatter;
use strict;
use warnings;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use Edit;

##
# Bowtie-like output formatter.
#
sub bowtie {
	my ($fh, $rdname, $fw, $chr, $off, $seq, $qual, $mms, $edits) = @_;
	defined($rdname) || croak("No read name");
	defined($fw)     || croak("No fw");
	defined($chr)    || croak("No chr");
	defined($off)    || croak("No off");
	defined($seq)    || croak("No seq");
	defined($qual)   || croak("No qual");
	defined($mms)    || croak("No mms");
	defined($edits)  || croak("No edits");
	$fw = $fw eq "FW" ? "+" : "-";
	$edits = Edit::listToString($edits);
	print {$fh} "$rdname\t$fw\t$chr\t$off\t$seq\t$qual\t$mms\t$edits\n";
}

1;
