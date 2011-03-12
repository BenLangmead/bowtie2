#!/usr/bin/perl -w

#
#
#
#

use strict;
use warnings;
use Heap::Priority;
use FindBin qw($Bin);
use lib $Bin;
use Edit;
use PartialAln;

my $h = new Heap::Priority;

##
# Encapsulates everything needed to resume a search from a particular
# point (i.e. partial alignment).
#
sub new {
	my ($class, $name, $seqs, $snps) = @_;
	my $self = {
		# _seqs: name -> sequence
		_name => defined($name) ? $name : "noname",
		# _seqs: name -> sequence
		_seqs => defined($seqs) ? $seqs : { },
		# _snps: name,pos -> alleles
		_snps => defined($snps) ? $snps : { }
	};
	bless $self, $class;
	return $self;
}
