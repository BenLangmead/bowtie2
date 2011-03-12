#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 6/29/2010
#
# Module for reporting alignments given a heap of PartialAlns
# representing the (prioritized) alignments and a Reference.
#
#
# Accepts alignments
# Tells the caller when it can stop searching for more alignments
# Buffers alignments for previous reads until the buffer reaches a certain point, then purges them
#
#

package Reporter;
use strict;
use warnings;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use Test;
use ResultAlns;
use RefSearch;
use Formatter;

sub new {
	my ($class, $refs, $k, $m, $formatter, $flushIval, $reported) = @_;
	$k = 1 unless defined($k);
	$m = 2**30 unless defined($m);
	return bless {
		_refs => $refs || croak("No ref"),
		_k   => $k,
		_m   => $m,
		_formatter => $formatter || croak("No formatter"),
		_reported => $reported || 0,
		_pending => [ ],
		_flushIval => $flushIval || 1,
		_flushCur => 0
	}, $class;
}
sub refs      { return $_[0]->{_refs}      }
sub k         { return $_[0]->{_k}         }
sub m         { return $_[0]->{_m}         }
sub formatter { return $_[0]->{_formatter} }
sub flushIval { return $_[0]->{_flushIval} }
sub flushCur  { return $_[0]->{_flushCur}  }
sub pending   { return $_[0]->{_pending}   }
sub reported  { return $_[0]->{_reported}  }

##
# Reset state as though no alignments have been reported so far.
#
sub clear($) {
	my $self = shift;
	$self->{_pending} = [ ];
	$self->{_reported} = 0;
	$self->{_flushCur} = 0;
}

##
# Takes a ResultAlns object and either flushes it immediately or stores
# it in the "pending" list for subsequent flushing.
#
sub report($$) {
	my ($self, $pa) = @_;
	push @{$self->pending()}, $pa;
	if(++$self->{_flushCur} == $self->flushIval()) {
		$self->{_flushCur} = 0;
		$self->flush();
	}
	$self->{_reported}++;
}

##
# Flush all pending reports.  TODO: maybe give RefSearch some
# polymorphism and add a sanity-check method.
#
sub flush($) {
	my ($self) = @_;
	my $found = RefSearch::findOffsets(
		$self->refs(), $self->pending(), $self->k(), $self->m());
	# Now, for each ResultAlns, we iterate through the rstrs in
	# best-to-worst (in terms of cost) order and start accumulating
	# results until we bump up against a k or m limit, or until we're
	# done
	ALNS: for my $ra (@{$self->pending()}) {
		# Examine this ResultAlns's rstrs in best to worst order
		my @results = ();
		my $numResults = 0;
		my @sortedRs = sort { $ra->results()->{$a}->cost() <=>
		                      $ra->results()->{$b}->cost() } keys %{$ra->results()};
		RS: for my $rs (@sortedRs) {
			# Were there any matches for this ref str?
			if(defined($found->{$rs})) {
				# Premptively check if the number of candidate
				# placements for $rs pushes us above the -m ceiling
				my $placements = scalar(@{$found->{$rs}})/2;
				if($numResults + $placements > $self->m()) {
					# Exceeded -m ceiling
					next ALNS;
				}
				# There were matches, keep adding them to @results
				# until we run up against the k and/or m limit.
				my $ral = $ra->results()->{$rs}; # Get ResultAln
				defined($ral->edits()) || croak();
				for(my $i = 0; $i < scalar(@{$found->{$rs}}); $i += 2) {
					my ($chr, $off) = ($found->{$rs}->[$i+0], $found->{$rs}->[$i+1]);
					if(scalar(@results) < $self->k()) {
						push @results, [
							$ral->read()->name(), $ral->fw(), $chr, $off, $ral->read()->seq(), $ral->read()->qual(), 0, $ral->edits() ];
					}
					scalar(@results) <= $self->m() || croak();
				}
			}
		}
		# Finished compiling results for this Read
		for my $res (@results) { $self->formatter()->(@$res); }
	}
}

1;
