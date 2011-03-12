#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 6/29/2010
#
# Encapsulates the search results for a given read.  For each valid
# alignment found, we record a map from cost to (a) reference string,
# and (b) edits.
#

package ResultAlns;
use strict;
use warnings;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use Test;
use Edit;
use Reference;
use Read;
use ResultAln;

sub new {
	my ($class, $refs, $results) = @_;
	my $self = bless {
		_refs  => $refs || croak("No refs"),
		_results => {},
		_minCost => 2 ** 30,
		_maxCost => 0,
		_num => 0,
		_disallowRepeats => 0
	}, $class;
	# Deep copy $rstre into $self->
	if(defined($results)) {
		for my $k (keys %{$results}) {
			$self->add($results->{$k});
		}
	}
	return $self;
}

##
# Clear state as though no alignments have been accumulated.
#
sub clearAlignments {
	my $self = shift;
	$self->{_results} = {};
	$self->{_minCost} = 2 ** 30;
	$self->{_maxCost} = 0;
	$self->{_num} = 0;
}

sub results { return $_[0]->{_results} }
sub refs    { return $_[0]->{_refs}    }
sub num     { return $_[0]->{_num}     }
sub minCost { return $_[0]->{_minCost} }
sub maxCost { return $_[0]->{_maxCost} }
sub disallowRepeats { return $_[0]->{_disallowRepeats} }
sub empty { return scalar(keys %{$_[0]->rstre()}) == 0 }

##
# Add a new result
#
sub add($$) {
	my ($self, $raln) = @_;
	$self->disallowRepeats() &&
		defined($self->results()->{$raln->rstr()}) &&
			croak("ResultAlns got same rstr (\"".$raln->rstr()."\") more than once\n");
	if(!$self->disallowRepeats() && defined($self->results()->{$raln->rstr()})) {
		$raln->cost() >= $self->results()->{$raln->rstr()}->cost() || croak("Better rstr came later\n");
		return;
	}
	$self->results()->{$raln->rstr()} = $raln;
	$self->{_minCost} = $raln->cost() if $raln->cost() < $self->{_minCost};
	$self->{_maxCost} = $raln->cost() if $raln->cost() > $self->{_maxCost};
	$self->minCost() <= $raln->cost() || croak();
	$self->maxCost() >= $raln->cost() || croak();
	$self->{_num}++;
}

sub test1 {
	my @refs = ();
	Reference::newFromHash({ "seq1" => "ACGT"x5 }, \@refs);
	my $res = new ResultAlns(\@refs);
	my $r = new Read("blah", "TTACGAACCCAACGTATCG", "I"x20, 0, "FW", "?");
	$res->add(new ResultAln(0, "ACGTACGT", "FW", $r, []));
	$res->minCost() == 0 || croak();
	$res->maxCost() == 0 || croak();
	$res->add(new ResultAln(3, "CCGTACGT", "FW", $r, [ new Edit("blah", 0, "A", "C") ]));
	$res->minCost() == 0 || croak();
	$res->maxCost() == 3 || croak();
	$res->num() == 2 || croak();
	return 1;
}

if($0 =~ /ResultAlns\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
	print "Test \"test1\"..."; test1(); print "PASSED\n";
}

1;
