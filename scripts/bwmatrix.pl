#!/usr/bin/perl -w

##
# bwmatrix.pl
# Author: Ben Langmead
#
# Simple, naive script that constructs a Burrows-Wheeler matrix from a
# string input.  fmap and bmap are a simple way of ensuring that $ is
# treated as > other characters, as in Bowtie.
#

use strict;
use warnings;

# Replace $s with ~s
my %fmap = ( "\$" => "~" );

# Re-replace ~s with $s before printing result
my %bmap = ( "~" => "\$" );

my @mat = ();
for my $s (@ARGV) {
	@mat = ();
	# Forward map
	for(my $i = 0; $i < length($s); $i++) {
		my $c = substr($s, $i, 1);
		substr($s, $i, 1) = $fmap{$c} if defined($fmap{$c});
	}
	# Produce all cyclic rotations
	for(my $i = 0; $i < length($s); $i++) {
		push @mat, $s;
		# Rotate s
		my $c = substr($s, 0, 1);
		$s = substr($s, 1) . $c;
	}
	@mat = sort @mat;
	# Backward map
	my $i = 0;
	for my $ss (@mat) {
		for(my $i = 0; $i < length($ss); $i++) {
			my $c = substr($ss, $i, 1);
			substr($ss, $i, 1) = $bmap{$c} if defined($bmap{$c});
		}
		printf "%03d $ss\n", $i;
		$i++;
	}
}
