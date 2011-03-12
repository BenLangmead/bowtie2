#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 6/29/2010
#

package RefSearch;
use strict;
use warnings;
use Algorithm::AhoCorasick qw(find_all);;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use Test;
use Read;
use Reference;
use ResultAlns;
use ResultAln;

##
# Takes a Reference, a list of alignments results (a set of reference
# strings per read), and ceilings on the number of alignments that (a)
# can be reported (-k), and (b) can be detected before the read is
# considered repetitive (-m).
#
sub findOffsets($$$$) {
	my ($refs, $ress, $k, $m) = @_;
	# Build an array of all the rstrs
	my @rstrs = ();
	for my $res (@$ress) { push @rstrs, keys %{$res->results()}; }
	return {} if scalar(@rstrs) == 0;
	my $found = {};
	# Make a batch of keywords from all of the different reference
	# strings found
	for my $r (@$refs) {
		my $acfound = find_all($r->seq(), @rstrs);
		next unless defined($acfound);
		# For each pos
		for my $pos (keys %$acfound) {
			# For each result at this pos
			for my $pres (@{$acfound->{$pos}}) {
				# Map the rstr to chr:pos
				push @{$found->{$pres}}, ($r->name(), $pos);
				#print "Found \"$pres\" at ".$r->name().":$pos\n";
			}
		}
	}
	return $found;
}

sub test1() {
	my @refs = ();
	Reference::newFromHash(
		{ "seq1" => "AGATACTGCGGAGAGGAAAACCCTTCTGAATCGAGCTACGATTTATGCG",
	#                                  AACCCTTCTGAATC
	#                0123456789012345678
		  "seq2" => "GAGGACCCATATCTACTTACGAACCACAACGTATCGATCGAGTTATGCGTATTCGCGG" },
	#                                 TACGAACCACAACGTATC
	#                012345678901234567
		\@refs
	);
	my $r = new Read("blah", "TTACGAACCCAACGTATCG", "I"x20, 0, "FW", "?");
	my $res = new ResultAlns(\@refs);
	$res->add(new ResultAln(0, "TACGAACCACAACGTATC", "FW", $r, []));
	my $res2 = new ResultAlns(\@refs);
	$res2->add(new ResultAln(3, "AACCCTTCTGAATC", "FW", $r, [ new Edit("blah", 0, "C", "A") ]));
	my $found = findOffsets(\@refs, [$res, $res2], 1, 1);
	scalar(keys %$found) == 2 || croak("Expected 2 results, got ".scalar(keys %$found));
	defined($found->{TACGAACCACAACGTATC}) || croak("Expected a result for rstr \"TACGAACCACAACGTATC\" but didn't\n");
	defined($found->{AACCCTTCTGAATC})     || croak("Expected a result for rstr \"AACCCTTCTGAATC\" but didn't\n");
	scalar(@{$found->{TACGAACCACAACGTATC}}) == 2 || croak("Expected length 2, got ".scalar(@{$found->{TACGAACCACAACGTATC}})."\n");
	scalar(@{$found->{AACCCTTCTGAATC}}) == 2 || croak("Expected length 2, got ".scalar(@{$found->{AACCCTTCTGAATC}})."\n");
	return 1;
}

if($0 =~ /RefSearch\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
	print "Test \"test1\"..."; test1(); print "PASSED\n";
}

1;
