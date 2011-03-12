#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 6/24/2010
#
# Routines for loading and indexing references.
#

package BloomIndex;
use strict;
use warnings;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use Test;
use Bloom::Faster;

##
# Given a ref to a list of references, return a closure that uses the
# Perl index function to return 1 iff the query is truly a reference
# substring. 
#
sub perlIndexBloom {
	my ($refs, $m, $k, $merLens, $msg) = @_;
	my $bloom = new Bloom::Faster({m => $m, k => $k});
	my $ents = 0;
	my $dups = 0;
	for my $r (@{$refs}) {
		for my $i (@{$merLens}) {
			for(my $j = 0; $j <= length($r->seq())-$i; $j++) {
				my $str = substr($r->seq(), $j, $i);
				$dups += $bloom->add($str);
				$ents++;
			}
		}
	}
	defined($msg) &&
		$msg->("Bloom filled with $ents entries, $dups (%.2f%%) duplicates",
		       $dups*100.0/($ents || 1));
	@{$merLens} = sort @{$merLens};
	return sub {
		my @strs = ();
		for(my $i = 0; $i < scalar(@{$merLens}); $i++) {
			if(length($_[0]) <= $merLens->[$i]) {
				# We could grab more than just the leftmost substring
				# here
				push @strs, substr($_[0], 0, $merLens->[$i]);
			}
		}
		for(my $i = 0; $i < scalar(@strs); $i++) {
			# If the substring of the query is not in the bloom, the
			# query definitely isn't
			return 0 unless($bloom->check($strs[$i]));
		}
		return 1; # may be a false positive
	}
}

1;
