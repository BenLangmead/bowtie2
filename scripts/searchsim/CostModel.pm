#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 6/21/2010
#
#

package CostModel;
use strict;
use warnings;
use Carp;

sub new {
	my ($class, $name, $mmfunc, $readGapo, $readGape, $refGapo, $refGape) = @_;
	return bless {
		_name     => $name   || "noname",
		_mmfunc   => $mmfunc,
		_readGapo => $readGapo,
		_readGape => $readGape,
		_refGapo  => $refGapo,
		_refGape  => $refGape
	}, $class;
}
sub name     { return $_[0]->{_name}     }
sub mmfunc   { return $_[0]->{_mmfunc}   }
sub readGapo { return $_[0]->{_readGapo} }
sub readGape { return $_[0]->{_readGape} }
sub refGapo  { return $_[0]->{_refGapo}  }
sub refGape  { return $_[0]->{_refGape}  }

sub readGap($$) { return $_[1] ? $_[0]->readGape() : $_[0]->readGapo() }
sub refGap($$)  { return $_[1] ? $_[0]->refGape()  : $_[0]->refGapo()  }

sub qualIdentityFunc($) {
	my $q = ord($_[0]) - 33;
	$q >= 0 || croak("Bad quality char: $_[0]\n");
	return $q;
}

my %maqRound = (
	 0 =>  0,  1 =>  0,  2 =>  0,  3 =>  0,  4 =>  0,
	 5 => 10,  6 => 10,  7 => 10,  8 => 10,  9 => 10,
	10 => 10, 11 => 10, 12 => 10, 13 => 10, 14 => 10,
	15 => 20, 16 => 20, 17 => 20, 18 => 20, 19 => 20,
	20 => 20, 21 => 20, 22 => 20, 23 => 20, 24 => 20,
	25 => 30, 26 => 30, 27 => 30, 28 => 30, 29 => 30,
	30 => 30, 31 => 30, 32 => 30, 33 => 30, 34 => 30,
	35 => 30, 36 => 30, 37 => 30, 38 => 30, 39 => 30,
	40 => 30, 41 => 30, 42 => 30, 43 => 30, 44 => 30
);

sub qualMaqRoundFunc($) {
	my $q = ord($_[0]) - 33;
	defined($maqRound{$q}) || die;
	return $maqRound{$q};
}

##
# Return a Maq-like cost model.
#
sub maqLike() {
	my ($round) = @_;
	my $qf = $round ? sub { qualMaqRoundFunc($_[0]) } :
	                  sub { qualIdentityFunc($_[0]) };
	return CostModel->new("Maq-like", $qf, 1, 1, 1, 1);
}

sub bwaDefault() {
	# 3, 11, 4
	return CostModel->new("BWA-like", sub { return 3 }, 11, 4, 11, 4);
}

sub allOnes() {
	# 3, 11, 4
	return CostModel->new("All-1s", sub { return 1 }, 1, 1, 1, 1);
}

1;
