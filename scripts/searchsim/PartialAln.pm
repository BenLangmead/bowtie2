#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 6/21/2010
#
#

package PartialAln;
use strict;
use warnings;
use Heap::Priority;
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval nanosleep
                    clock stat );
use FindBin qw($Bin); 
use lib $Bin;
use Edit;
use Carp;
use Test;
use CostModel;
use SearchRoot;
use SearchScheme;
use Constraint;
use Read;
use Reference;
use ResultAlns;

sub new {
	scalar(@_) == 15 || croak("Expected 15 arguments, got ".scalar(@_)."\n");
	my ($class, $name, $orient, $qdepth, $qrootdepth, $rdepth, $cons,
	    $ocons, $read, $root, $edits, $rstr, $ralstr, $qalstr, $pingPongs) = @_;
	defined($qdepth) || croak("No read depth");
	defined($rdepth) || croak("No ref depth");
	defined($pingPongs) || croak("No # ping pongs");
	my $self = bless {
		_name   => $name   || "noname",
		_orient => $orient || croak("No orientation"),
		_qdepth => $qdepth,
		_qrootdepth => $qrootdepth,
		_rdepth => $rdepth,
		_cons   => $cons   || croak("No remaining per-zone constraints"),
		_ocons  => $ocons  || croak("No remaining overall constraints"),
		_read   => $read   || croak("No Read"),
		_root   => $root   || croak("No SearchRoot"),
		_edits  => $edits  || croak("No edits"),
		_rstr   => $rstr   || "",
		_ralstr => $ralstr || "",
		_qalstr => $qalstr || "",
		_pingPongs => $pingPongs
	}, $class;
	defined($self->pingPongs()) || croak();
	return $self;
}
sub name   { return $_[0]->{_name}   }
sub orient { return $_[0]->{_orient} }
sub qdepth { return $_[0]->{_qdepth} }
sub qrootdepth { return $_[0]->{_qrootdepth} }
sub rdepth { return $_[0]->{_rdepth} }
sub cons   { return $_[0]->{_cons}   }
sub ocons  { return $_[0]->{_ocons}  }
sub read   { return $_[0]->{_read}   }
sub root   { return $_[0]->{_root}   }
sub rstr   { return $_[0]->{_rstr}   }
sub ralstr { return $_[0]->{_ralstr} }
sub qalstr { return $_[0]->{_qalstr} }
sub edits  { return $_[0]->{_edits}  }
sub pingPongs { return $_[0]->{_pingPongs} }
sub color  { return $_[0]->read()->color(); }
sub fw     { return $_[0]->root()->fw();    }
sub cost   { return $_[0]->ocons()->cost(); }
sub costDiff { return $_[0]->ocons()->costDiff($_[0]->root()->ocons()); }

sub toString {
	my $self = shift;
	return      $self->orient()
	       .":".$self->fw()
	       .":".$self->qdepth()
	       .":".$self->qrootdepth()
	       .":".$self->ralstr()
	       .":".$self->qalstr()
	       #.":".$self->root()->name()
	        ;
}

##
# Called when a PartialAln runs into the edge of the read and needs to
# return to the root position to extend in the opposite direction.
#
sub pingPong($) {
	my $self = shift;
	# Flip the orientation; done!
	$self->{_orient} = ($self->{_orient} eq "RtL") ? "LtR" : "RtL";
	$self->{_qrootdepth} = 0;
	Edit::flip($self->{_edits}, $self->qdepth());
	$self->{_pingPongs}++;
	$self->{_pingPongs} <= 1 || croak("Bad number of ping-pongs: $self->{_pingPongs}\n");
}

##
# Print the stacked version of this partial alignment, with reference
# on top and query on bottom.
#
sub print($$) {
	my ($self, $fh) = @_;
	print {$fh} $self->ralstr()."\n";
	for(my $i = 0; $i < length($self->ralstr()); $i++) {
		my $r = substr($self->ralstr(), $i, 1);
		my $q = substr($self->qalstr(), $i, 1);
		print {$fh} ($r eq $q) ? "|" : " ";
	}
	print {$fh} $self->qalstr()."\n";
}

##
# For a given SearchRoot, return a partial alignment representing the
# root of the search space.
#
sub newFromSearchRoot {
	my ($class, $root, $read) = @_;
	defined($root) || croak("No SearchRoot\n");
	$root->instantiated() || croak("SearchRoot not instantiated\n");
	my $first = substr($root->zone_str(), 0, 1);
	my $last  = substr($root->zone_str(), -1);
	my $ori = "LtR";
	my $lhs = $root->off();
	my $rhs = $lhs + $root->len();
	($lhs >= 0 && $lhs <  $read->len()) || croak("Bad LHS for seed $lhs\n");
	($rhs >= 0 && $rhs <= $read->len()) || croak("Bad RHS for seed $rhs\n");
	my $rootDepth = $lhs;
	if(ord($first) > ord($last)) {
		# Start search in right-to-left direction
		$ori = "RtL";
		$rootDepth = $read->len() - $rhs;
	}
	return bless {
		_name   => $root->name(),
		_orient => $ori,
		_qdepth => 0,
		_qrootdepth => $rootDepth,
		_rdepth => 0,
		_cons   => $root->cons(),
		_ocons  => $root->ocons(),
		_read   => $read || croak("No Read"),
		_root   => $root,
		_edits  => [ ],
		_rstr   => "",
		_ralstr => "",
		_qalstr => "",
		_pingPongs => 0
	}, $class;
}

##
# Return a clone of the current PartialAln, with perhaps a few fields
# tweaked.
#
sub clone {
	my ($self, $qdepth, $rdepth, $nrstr, $nralstr, $nqalstr) = @_;
	defined($self) || croak("No self\n");
	# Copy zone constraints
	my $ncons = { };
	for my $con (keys %{$self->cons()}) {
		$ncons->{$con} = $self->cons()->{$con}->clone();
	}
	# Copy overall constraints
	my $nocons = $self->ocons()->clone();
	# Either copy or override these fields
	$qdepth  = $self->qdepth() unless defined($qdepth);
	$rdepth  = $self->rdepth() unless defined($rdepth);
	$nrstr   = $self->rstr()   unless defined($nrstr);
	$nralstr = $self->ralstr() unless defined($nralstr);
	$nqalstr = $self->qalstr() unless defined($nqalstr);
	defined($self->pingPongs()) || croak();
	# Construct and return the clone
	my $ret = PartialAln->new(
		$self->name(),
		$self->orient(),
		$qdepth,
		$self->qrootdepth(),
		$rdepth,
		$ncons,
		$nocons,
		$self->read(),
		$self->root(),
		[],
		$nrstr,
		$nralstr,
		$nqalstr,
		$self->pingPongs()
	);
	for my $e (@{$self->edits()}) {
		push @{$ret->edits()}, $e;
	}
	return $ret;
}

##
# Create a root PartialAln for the given root/read and add it to the
# given heap.
#
sub addSearchRootToHeap($$$) {
	my ($heap, $root, $read) = @_;
	my $aln = PartialAln->newFromSearchRoot($root, $read);
	$aln->repOk(1) || croak();
	$heap->add($aln, $aln->cost());
}

##
# Create a set of root PartialAlns for teh given scheme/read and add
# them to the given heap.
#
sub addSearchSchemeToHeap($$$) {
	my ($heap, $scheme, $read) = @_;
	for my $sr (@{$scheme->roots()}) {
		addSearchRootToHeap($heap, $sr, $read);
	}
}

##
# Return true iff introducing a reference gap would be a gap extension
# rather than a gap open.
#
sub extendRef($)  {
	my $self = shift;
	return 0 if scalar(@{$self->edits()}) == 0;
	$self->edits()->[-1]->off() <= $self->qdepth() || croak();
	return $self->edits()->[-1]->off() == $self->qdepth()-1 &&
	       $self->edits()->[-1]->is_delete();
}

##
# Return true iff introducing a read gap would be a gap extension
# rather than a gap open.
#
sub extendRead($) {
	my $self = shift;
	return 0 if scalar(@{$self->edits()}) == 0;
	$self->edits()->[-1]->off() <= $self->qdepth() || croak();
	return $self->edits()->[-1]->off() == $self->qdepth() &&
	       $self->edits()->[-1]->is_insert();
}

##
# Get the zone that applies to the next edit.
#
sub zone($) {
	my $self = shift;
	my $offset = $self->qdepth();
	$offset < length($self->root()->zone_str()) || return undef;
	$offset = -$offset-1 if $self->orient() eq "RtL";
	my $z = uc substr($self->root()->zone_str(), $offset, 1);
	return $z;
}

## Small utility function; returns true iff arg is a lowercase letter
sub islower($) { return $_[0] =~ /[a-z]/; }

##
# Get the zone that applies to the next insertion; insertions require
# a special rule that avoids redundant alignments with insertions in
# the cusp between two zones.
#
sub zoneInsert($) {
	my $self = shift;
	my $offset = $self->qdepth();
	$offset < length($self->root()->zone_str()) || return undef;
	$offset = -$offset-1 if $self->orient() eq "RtL";
	my $z = substr($self->root()->zone_str(), $offset, 1);
	if(islower($z)) {
		# Go with the previous zone
		$offset += ($self->orient() eq "RtL") ? -1 : 1;
		($offset >= 0 && $offset < length($self->root()->zone_str())) || die;
		$z = substr($self->root()->zone_str(), $offset, 1);
		islower($z) && die;
	}
	return $z;
}

##
# Return non-zero iff the search is still on the seed zone.
#
sub withinSeed($) {
	return $_[0]->qdepth() < $_[0]->root()->len();
}

##
# Check whether this partial alignment is at a point where it *has* to
# match.
#
sub mustMatch($) {
	my ($self) = @_;
	my $zone = $self->zone();
	defined($zone) || return 0;
	defined($self->cons()->{$zone}) || die "Bad zone $zone\n";
	return $self->cons()->{$zone}->mustMatch();
}

##
# Do the constraints categorically allow a mismatch as a way of
# extending this partial alignment.  A particular mismatch may still be
# disallowed because of a cost constraint.
#
sub canMismatch($) {
	defined($_[0]->zone()) || return $_[0]->ocons()->canMismatch();
	defined($_[0]->cons()->{$_[0]->zone()}) || croak("No such zone as ".$_[0]->zone()."\n");
	return $_[0]->cons()->{$_[0]->zone()}->canMismatch() &&
	       $_[0]->ocons()->canMismatch();
}

##
# Do the constraints categorically allow a gap as a way of extending
# this partial alignment.  A particular mismatch may still be
# disallowed because of a cost constraint.
#
sub canGap($) {
	defined($_[0]->zone()) || return $_[0]->ocons()->canGap();
	defined($_[0]->cons()->{$_[0]->zone()}) || croak("No such zone as ".$_[0]->zone()."\n");
	return $_[0]->cons()->{$_[0]->zone()}->canGap() &&
	       $_[0]->ocons()->canGap();
}

##
# Do the constraints categorically allow an insertion (gap in the read)
# as a way of extending this partial alignment.  A particular insertion
# may still be disallowed because of, say, a cost constraint.
#
sub canInsert($) {
	defined($_[0]->zone()) || return $_[0]->ocons()->canInsert();
	defined($_[0]->cons()->{$_[0]->zoneInsert()}) || croak("No such zone as ".$_[0]->zoneInsert()."\n");
	return $_[0]->cons()->{$_[0]->zoneInsert()}->canInsert() &&
	       $_[0]->ocons()->canInsert();
}

##
# Do the constraints categorically allow a deletion (gap in the ref)
# as a way of extending this partial alignment.  A particular deletion
# may still be disallowed because of, say, a cost constraint.
#
sub canDelete($) {
	defined($_[0]->zone()) || return $_[0]->ocons()->canDelete();
	defined($_[0]->cons()->{$_[0]->zone()}) || croak("No such zone as ".$_[0]->zone()."\n");
	return $_[0]->cons()->{$_[0]->zone()}->canDelete() &&
	       $_[0]->ocons()->canDelete();
}

##
# Do the constraints allow us to incur the given marginal cost?
#
sub canTakeCost($$) {
	my ($self, $cost) = @_;
	defined($self->zone()) || return $self->cost() >= $cost;
	defined($self->cons()->{$self->zone()}) || croak("No such zone as ".$self->zone()."\n");
	$cost > 0 || croak("Bad cost: $cost");
	return $self->cons()->{$self->zone()}->cost() >= $cost &&
	       $self->cost() >= $cost;
}

##
# Append a character on the "far" end of the rstr (far = farther from
# the root, rstr = reference string aligned to so far).
#
sub appendToRstr($$) {
	my ($self, $c) = @_;
	return ($self->orient() eq "RtL") ? $c.$self->rstr() : $self->rstr().$c;
}

##
# Append a character on the "far" end of the qstr (far = farther from
# the root, qstr = read string aligned to so far).
#
sub appendToRAlStr($$) {
	my ($self, $c) = @_;
	return ($self->orient() eq "RtL") ? $c.$self->ralstr() : $self->ralstr().$c;
}
sub appendToQAlStr($$) {
	my ($self, $c) = @_;
	return ($self->orient() eq "RtL") ? $c.$self->qalstr() : $self->qalstr().$c;
}

##
# Pops the top PartialAln off the top of the head and calls advance on
# it.
#
sub popAndAdvance {
	my ($costModel, $heap, $strset, $results, $log, $redunCheck, $progress, $cfg) = @_;
	$cfg = { } unless defined($cfg);
	$log = sub { } unless defined($log);
	return 0 if $heap->get_heap() == 0;
	my $pa = $heap->pop();
	$pa->advance($costModel, $heap, $strset, $results, $log, $redunCheck, $progress, $cfg);
	return scalar($heap->get_heap());
}

sub readCharsLeft($) {
	my $self = shift;
	my $left = $self->read()->len() - ($self->qdepth() + $self->qrootdepth());
	$left >= 0 || croak("Bad number of read chars left: $left\n");
	return $left;
}

##
# Takes a partial alignment (self), a priority heap of partial
# alignments, and a reference to a function that, given a hypothetical
# substring, will return 1 if that substring occurs in the reference or
# 0 if it does not.
#
sub advance {
	my ($self, $costModel, $heap, $strset, $results, $log, $redunCheck, $progress, $cfg) = @_;
	$cfg = { } unless defined($cfg);
	$log = sub { } unless defined($log);
	defined($redunCheck) || ($redunCheck = sub { return 0 });
	$redunCheck->($self) &&
		print STDERR "Failed non-redundant partial alignment check: \"".$self->toString()."\"\n";
	defined($progress) && $progress->($heap, $self);
	$log->("Entering advance() with ".scalar($heap->get_heap())." partials on the heapl rstr=".$self->rstr());
	if($self->readCharsLeft() == 0) {
		# Hit the end of the read; either stop here and report the
		# alignment or ping-pong back to the search root and extend in
		# the opposite direction.
		if($self->qrootdepth() > 0) {
			$log->("  extended to edge of read; ping-ponging back to the root");
			$self->pingPong();
			$self->readCharsLeft() > 0 || croak();
			# change orientation
			# keep qdepth
		} else {
			$log->("  no chars left; done!");
			$self->repOk(1) || croak();
			my $naln = new ResultAln(
				$self->costDiff(), $self->rstr(), $self->fw(),
				$self->read(), $self->edits());
			#print "XXX ResultAln generated\n";
			$results->add($naln);
			return;
		}
	}
	my ($c, $q) = $self->read()->at(
		$self->qdepth() + $self->qrootdepth(), $self->orient(), $self->color());
	my @strsetResults = (undef, undef, undef, undef);
	my $nrstr = $self->appendToRstr($c);
	$log->("  read at offset ".$self->qdepth()." = ($c, $q), orient=".$self->orient());
	if($strset->($nrstr)) {
		# We can proceed on a match
		my $npa = $self->clone(
			$self->qdepth()+1,
			$self->rdepth()+1,
			$nrstr,
			$self->appendToRAlStr($c),
			$self->appendToQAlStr($c));
		# Prioritize by cost budget remaining
		$npa->repOk(1) || croak();
		#$redunCheck->($npa) ||
		#	print STDERR "Failed non-redundant partial alignment check: \"".$npa->toString()."\"\n";
		$heap->add($npa, $npa->cost());
		$log->("  added match partial");
		$strsetResults[0] = 1;
	} else {
		$log->("  match partial \"$nrstr\" didn't occur in reference");
		$strsetResults[0] = 0;
	}
	if($self->canMismatch()) {
		# Try to proceed on each of the three different mismatches.
		# Check whether we still have a reference and whether we've
		# exceeded.
		for(my $i = 1; $i <= 3; $i++) {
			my $mmc = DNA::plus($c, $i);
			$log->("    trying mismatch partial for $mmc");
			$nrstr = $self->appendToRstr($mmc);
			# Is it acceptable according to the cost model?
			my $margMmCost = $costModel->mmfunc()->($q);
			if($self->canTakeCost($margMmCost)) {
				# Haven't exceeded cost budget; now check to see if the
				# hypothetical substring exists in the reference
				if($strset->($nrstr)) {
					# Bingo, a valid outgoing path
					# Construct a new constraints hash
					my $npa = $self->clone(
						$self->qdepth()+1,
						$self->rdepth()+1,
						$nrstr,
						$self->appendToRAlStr($mmc),
						$self->appendToQAlStr($c));
					my $zone = $self->zone();
					defined($zone) && $npa->cons()->{$zone}->chargeMismatch($margMmCost);
					$npa->ocons()->chargeMismatch($margMmCost);
					push @{$npa->edits()}, Edit->new($self->name(), $self->qdepth(), $c, $mmc);
					$npa->repOk(1) || croak();
					#$redunCheck->($npa) ||
					#	print STDERR "Failed non-redundant partial alignment check: \"".$npa->toString()."\"\n";
					$heap->add($npa, $npa->cost());
					$log->("      added partial alignment");
					$strsetResults[$i] = 1;
				} else {
					# The mismatch was legal according to the
					# constraints, but it didn't correspond to a
					# substring of the reference.
					$log->("      rejected because \"$nrstr\" not present in reference");
					$strsetResults[$i] = 0;
				}
			} else {
				# This mismatch would have been allowed categorically,
				# but this particular mismatch cost too much
				$log->("      rejected because no more mismatches allowed either overall or in zone ".$self->zone());
			}
		}
	}
	# Is there still at least 1 edit (or cost=1) allowed at this depth?
	if($self->canDelete()) {
		$log->("    trying delete");
		# Consider reference gap (deletion) first
		my $margDelCost = $costModel->refGap($self->extendRef());
		if($self->canTakeCost($margDelCost)) {
			$log->("      delete is cost-allowed");
			my $npa = $self->clone(
				$self->qdepth()+1,
				$self->rdepth(),
				$self->rstr(), # Don't add anything to this
				$self->appendToRAlStr("-"),
				$self->appendToQAlStr($c));
			my $zone = $self->zone();
			defined($zone) && $npa->cons()->{$zone}->chargeDelete($margDelCost);
			$npa->ocons()->chargeDelete($margDelCost);
			push @{$npa->edits()}, Edit->new($self->name(), $self->qdepth(), $c, "-");
			$npa->repOk(1) || croak();
			#$redunCheck->($npa) ||
			#	print STDERR "Failed non-redundant partial alignment check: \"".$npa->toString()."\"\n";
			$heap->add($npa, $npa->cost());
			$log->("        added partial alignment");
		} else {
			$log->("      delete is disallowed due to cost");
		}
	}
	if($self->canInsert()) {
		my $margInsCost = $costModel->readGap($self->extendRead());
		$log->("    trying insertions");
		if($self->canTakeCost($margInsCost)) {
			$log->("      insert is cost-allowed");
			# Try to proceed on each of the four different insertions.
			for(my $i = 0; $i <= 3; $i++) {
				my $insc = DNA::plus($c, $i);
				$log->("        trying insertion of $insc");
				$nrstr = $self->appendToRstr($insc);
				$strsetResults[$i] = $strset->($nrstr) unless defined($strsetResults[$i]);
				if($strsetResults[$i]) {
					# The mismatch was legal according to the
					# constraints and it corresponds to an actual
					# substring of the reference.
					$log->((" "x10)."accepted; \"$nrstr\" is a reference substring");
					my $npa = $self->clone(
						$self->qdepth(),
						$self->rdepth()+1,
						$nrstr, # Don't add anything to this
						$self->appendToRAlStr($insc),
						$self->appendToQAlStr("-"));
					my $zone = $self->zoneInsert();
					defined($zone) && $npa->cons()->{$zone}->chargeInsert($margInsCost);
					$npa->ocons()->chargeInsert($margInsCost);
					push @{$npa->edits()}, Edit->new($self->name(), $self->qdepth(), "-", $insc);
					$npa->repOk(1) || croak();
					#$redunCheck->($npa) ||
					#	print STDERR "Failed non-redundant partial alignment check: \"".$npa->toString()."\"\n";
					$heap->add($npa, $npa->cost());
					$log->("        added partial alignment");
				} else {
					# The mismatch was legal according to the
					# constraints, but it didn't correspond to a
					# substring of the reference.
					$log->("        rejected because \"$nrstr\" is not a reference substring");
				}
			}
		}
	}
	$log->("Exiting advance() with ".scalar($heap->get_heap())." partials on the heap");
	return;
}

##
# Check that the partial alignment is internally consistent.
#
sub repOk($$) {
	my ($self, $croak) = @_;
	if(!defined($self->read())) {
		$croak && croak("Read not defined\n");
		return 0;
	}
	if($self->qdepth() > $self->read()->len()) {
		$croak && croak("qdepth (".$self->qdepth().") > read length (".$self->read()->len().")\n");
		return 0;
	}
	if($self->orient() ne "RtL" && $self->orient() ne "LtR") {
		$croak && croak("Bad orientation: ".$self->orient()."\n");
		return 0;
	}
	# Make sure the difference between qdepth and rdepth is consistent
	# with the Edits
	my $qMinusR = 0;
	for my $e (@{$self->edits()}) {
		if($e->off > $self->qdepth()) {
			$croak && croak("Bad edit offset (".$e->off."), qdepth=".$self->qdepth()."\n");
			return 0;
		}
		if($e->off() == $self->qdepth()) {
			unless($e->is_insert()) {
				$croak && croak("Bad non-insert edit offset (".$e->off."), qdepth=".$self->qdepth().", read char=\"".$e->readch()."\", ref char=\"".$e->refch()."\"\n");
				return 0;
			}
		}
		$qMinusR-- if $e->is_insert();
		$qMinusR++ if $e->is_delete();
	}
	if($self->qdepth() - $self->rdepth() != $qMinusR) {
		$croak && croak("qdepth (".$self->qdepth().") - rdepth (".$self->rdepth().") != qMinusR ($qMinusR)\n");
		return 0; # rep not OK
	}
	if($self->rdepth() != length($self->rstr())) {
		$croak && croak("rdepth (".$self->rdepth().") != length(self->rstr) (".length($self->rstr())."); rstr=\"".$self->rstr()."\"\n");
		return 0; # rep not OK
	}
	my $qstr = $self->qalstr();
	$qstr =~ s/-//g;
	if($self->qdepth() != length($qstr)) {
		$croak && croak("qdepth (".$self->qdepth().") != length(qstr); qstr=\"$qstr\"\n");
		return 0; # rep not OK
	}
	my $ralstr = $self->ralstr();
	$ralstr =~ s/-//g;
	if($ralstr ne $self->rstr()) {
		$croak && croak("ralstr \"".$self->ralstr()."\" with gaps removed \"$ralstr\" doesn't match rstr=\"".$self->rstr()."\"\n");
		return 0; # rep not OK
	}
	return 1;
}

my %testRedun = ();
my $testRedunCheck = sub {
	my $ret = defined($testRedun{$_[0]->toString()});
	$testRedun{$_[0]->toString()}++;
	$testRedun{$_[0]->toString()} < 2 || croak();
	return $ret;
};

my $progressCnt = 0;
my $progressCntWithinSeed = 0;
my $progressT0 = [ gettimeofday ];
my %progressDepthHash = ();
sub progressFunc($) {
	my $ival = shift;
	return sub {
		my ($heap, $self) = @_;
		$progressCntWithinSeed++ if $self->withinSeed();
		$progressDepthHash{$self->qdepth()}++;
		if((++$progressCnt % $ival) == 0) {
			my $elapsed = tv_interval($progressT0, [ gettimeofday ]);
			$progressT0 = [ gettimeofday ];
			my ($heapSz, $curCost) = ("?", "?");
			$heapSz = scalar($heap->get_heap()) if defined($heap);
			$curCost = $self->costDiff() if defined($self);
			printf "  searched $progressCnt partials (%0.2f pps, within seed: $progressCntWithinSeed, heap sz: $heapSz, min cost: $curCost)\n", ($ival/$elapsed);
			#for my $k (sort {$a <=> $b} keys %progressDepthHash) {
			#	print "    $k: $progressDepthHash{$k}\n";
			#}
		}
	};
}

sub test1() {
	my $sra = new SearchRoot(
		"blah", undef, "FW", 0, 28, 0, "01",
		{ 0 => Constraint::none(), 1 => Constraint::none() },
		Constraint::none(),
		1, 1, 0);
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my $sr = $sra->instantiate($r);
	my $e = PartialAln->new(
		"blah",
		"RtL",
		10, 0, 11,
		{ 0 => Constraint::zero() },
		Constraint::score(100),
		$r, $sr,
		[ Edit->new("blah", 4, "-", "a") ],
		"aaaaaaaaaaa",
		"aaaaaaaaaaa",
		"aaaa-aaaaaa",
		0);
	$e->repOk(1);
	return 1;
}

sub test2() {
	my $sra = new SearchRoot(
		"blah", undef, "FW", 0, 28, 0, "01",
		{ 0 => Constraint::none(), 1 => Constraint::none() },
		Constraint::none(),
		1, 1, 0);
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my @srs = $sra->instantiate($r);
	scalar(@srs) == 1 || croak();
	my $sr = $srs[0];
	my $e = PartialAln->new(
		"blah",
		"RtL",
		10, 0, 11,
		{ 0 => Constraint::zero() },
		Constraint::score(100),
		$r, $sr,
		[ ],
		"aaaaaaaaaaa",
		"aaaaaaaaaaa",
		"aaaa-aaaaaa",
		0);
	$e->repOk(0) && die;
	return 1;
}

sub test3() {
	my $sra = new SearchRoot(
		"blah", undef, "FW", 0, 28, 0, "01",
		{ 0 => Constraint::none(), 1 => Constraint::none() },
		Constraint::none(),
		1, 1, 0);
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my @srs = $sra->instantiate($r);
	scalar(@srs) == 1 || croak();
	my $sr = $srs[0];
	my $e = PartialAln->new(
		"blah",
		"RtL",
		10, 0, 11,
		{ 0 => Constraint::zero() },
		Constraint::score(100),
		$r, $sr,
		[ Edit->new("blah", 4, "-", "a"),
		  Edit->new("blah", 6, "a", "-"),
		  Edit->new("blah", 8, "-", "a") ],
		"aaaaaaaaaaa",
		"aaaaaaaaaaa",
		"aaaaa-aaaaa",
		0);
	$e->repOk(1);
	return 1;
}

sub test4() {
	my $sra = new SearchRoot(
		"blah", undef, "FW", 0, 28, 0, "01",
		{ 0 => Constraint::none(), 1 => Constraint::none() },
		Constraint::none(),
		1, 1, 0);
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my @srs = $sra->instantiate($r);
	scalar(@srs) == 1 || croak();
	my $sr = $srs[0];
	my $e = PartialAln->new(
		"blah",
		"RtL",
		11, 0, 11,
		{ 0 => Constraint::zero() },
		Constraint::score(100),
		$r, $sr,
		[ Edit->new("blah", 4, "-", "a"),
		  Edit->new("blah", 6, "a", "-") ],
		"aaaaaaaaaaa",
		"aaaaaaaaaaa",
		"aaaaaaaaaaa",
		0);
	$e->repOk(1);
	return 1;
}

sub test5() {
	my $sra = new SearchRoot(
		"blah", undef, "FW", 0, 28, 0, "01",
		{ 0 => Constraint::none(), 1 => Constraint::none() },
		Constraint::none(),
		1, 1, 0);
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my @srs = $sra->instantiate($r);
	scalar(@srs) == 1 || croak();
	my $sr = $srs[0];
	my $e = PartialAln->new(
		"blah",
		"RtL",
		5, 0, 5,
		{ 0 => Constraint::zero() },
		Constraint::score(100),
		$r, $sr,
		[ Edit->new("blah", 4, "-", "a"),
		  Edit->new("blah", 6, "a", "-") ],
		"aaaaa",
		"aaaaa",
		"aaaaa",
		0);
	$e->repOk(0) && die;
	return 1;
}

sub test6() {
	my $s = new SearchRoot(
		"blah", undef, "FW", 0, 28, 0, "01",
		{ 0 => Constraint::none(), 1 => Constraint::none() },
		Constraint::none(),
		0, 0, 0);
	$s->repOkTest() || return 0;
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my @srs = $s->instantiate($r);
	scalar(@srs) == 1 || croak();
	my $s2 = $srs[0];
	my $pa = PartialAln->newFromSearchRoot($s2, $r);
	$pa->qdepth() == 0 || die;
	$pa->qdepth() == 0 || die;
	$pa->orient() eq "LtR" || die;
	$pa->read()->len() == 30 || die;
	return 1;
}

sub test7() {
	my @refs = ();
	Reference::newFromHash(
		{ "seq1" => "AGATACTGCGGAGAGGAAAACCCTTCTGAATCGAGCTACGATTTATGCG",
		  "seq2" => "GAGGACCCATATCTACTTACGAACCACAACGTATCGATCGAGTTATGCGTATTCGCGG" },
		\@refs
	);
	my $r = Read->new("blah", "TTACGAACCACAACGTATCG", "I"x20, 0, "FW", "?");
	my $h = new Heap::Priority;
	my $results = new ResultAlns(\@refs);
	my $cm = CostModel::bwaDefault();
	for my $strset (Reference::perlHashStrset(\@refs), Reference::perlIndexStrset(\@refs)) {
		for my $ssa (SearchScheme::scheme_one(1, 1, 0, undef, undef, undef, undef, 0.4),
		             SearchScheme::scheme_one_simple(1, 1, 0, undef, undef, undef, undef, 0.4))
		{
			$h->get_heap() == 0 || die;
			$results->num() == 0 || die;
			my $ss = $ssa->instantiate($r);
			my $results = new ResultAlns(\@refs);
			addSearchSchemeToHeap($h, $ss, $r);
			if($ssa->name() eq "OneMismatchSilly") {
				$h->get_heap() == 1 ||
					croak("Expected 1 items on PartialAln heap, got ".scalar($h->get_heap())."\n");
			} else {
				$h->get_heap() == 2 ||
					croak("Expected 2 items on PartialAln heap, got ".scalar($h->get_heap())."\n");
			}
			my $advances = 0;
			my $log = sub { }; # sub { print STDERR "$_[0]\n" };
			%testRedun = ();
			while(popAndAdvance($cm, $h, $strset, $results, $log, $testRedunCheck, progressFunc(100))) { $advances++ }
			$advances > 0 || croak("Expected > 0 advances, but got $advances\n");
			$results->num() == 1 ||
				croak("Expected 1 result, got ".$results->num()." after $advances advances\n");
			scalar(keys %{$results->results()}) == 1 ||
				croak("Expected 1 result rstr, got ".scalar(keys %{$results->results()})."\n");
			for my $k (keys %{$results->results()}) {
				$k eq "TTACGAACCACAACGTATCG" || croak();
			}
			$results->minCost() == 0 || croak("Expected minCost 0, got ".$results->minCost());
			$results->maxCost() == 0 || croak("Expected maxCost 0, got ".$results->maxCost());
			$results->clearAlignments();
		}
	}
	return 1;
}

sub test8() {
	my @refs = ();
	Reference::newFromHash(
		{ "seq1" => "AGATACTGCGGAGAGGAAAACCCTTCTGAATCGAGCTACGATTTATGCG",
		  "seq2" => "GAGGACCCATATCTACTTACGAACCACAACGTATCGATCGAGTTATGCGTATTCGCGG" },
		\@refs
	);
	my $r = Read->new("blah", "TTACGAACCtCAACGTATCG", "I"x20, 0, "FW", "?");
	my $h = new Heap::Priority;
	my $results = new ResultAlns(\@refs);
	my $cm = CostModel::bwaDefault();
	for my $strset (Reference::perlHashStrset(\@refs), Reference::perlIndexStrset(\@refs)) {
		for my $ssa (SearchScheme::scheme_one(1, 1, 0), SearchScheme::scheme_one_simple(1, 1, 0)) {
			$h->get_heap() == 0 || die;
			$results->num() == 0 || die;
			my $ss = $ssa->instantiate($r);
			addSearchSchemeToHeap($h, $ss, $r);
			if($ssa->name() eq "OneMismatchSilly") {
				$h->get_heap() == 1 ||
					croak("Expected 1 items on PartialAln heap, got ".scalar($h->get_heap())."\n");
			} else {
				$h->get_heap() == 2 ||
					croak("Expected 2 items on PartialAln heap, got ".scalar($h->get_heap())."\n");
			}
			my $advances = 0;
			my $log = sub { }; # sub { print STDERR "$_[0]\n" };
			%testRedun = ();
			while(popAndAdvance($cm, $h, $strset, $results, $log, $testRedunCheck, progressFunc(100))) { $advances++ }
			$advances > 0 || croak("Expected > 0 advances, but got $advances\n");
			scalar(keys %{$results->results()}) == 1 ||
				croak("Expected 1 result rstr, got ".scalar(keys %{$results->results()})."\n");
			for my $k (keys %{$results->results()}) {
				$k eq "TTACGAACCACAACGTATCG" || croak();
			}
			$results->minCost() == 3 || croak();
			$results->maxCost() == 3 || croak();
			$results->clearAlignments();
		}
	}
	return 1;
}

sub test9 {
	my $verbose = shift;
	my @refs = ();
	Reference::newFromHash(
		{ "seq1" => "AGATACTGCGGAGAGGAAAACCCTTCTGAATCGAGCTACGATTTATGCG",
		  "seq2" => "GAGGACCCATATCTACTTACGAACCACAACGTATCGATCGAGTTATGCGTATTCGCGG" },
		\@refs
	);
	my $r = Read->new("Read1", "TTACGAACCtCAACGTATCG", "I"x20, 0, "FW", "?");
	# TTACGAACCtCAACGTATCG
	#   ----seed----
	my $h = new Heap::Priority;
	my $results = new ResultAlns(\@refs);
	my $cm = CostModel::bwaDefault();
	for my $strset (Reference::perlHashStrset(\@refs), Reference::perlIndexStrset(\@refs)) {
		# Create a search scheme that puts the seed in the middle of
		# the read and doesn't expand or retract.
		for my $ssa (SearchScheme::scheme_one       (1, 1, 0, 2, 12, 0, 1, 0.4),
		             SearchScheme::scheme_one_simple(1, 1, 0, 2, 12, 0, 1, 0.4))
		{
			$h->get_heap() == 0 || die;
			$results->num() == 0 || die;
			my $ss = $ssa->instantiate($r);
			addSearchSchemeToHeap($h, $ss, $r);
			if($ssa->name() eq "OneMismatchSilly") {
				$h->get_heap() == 1 ||
					croak("Expected 1 items on PartialAln heap, got ".scalar($h->get_heap())."\n");
			} else {
				$h->get_heap() == 2 ||
					croak("Expected 2 items on PartialAln heap, got ".scalar($h->get_heap())."\n");
			}
			my $advances = 0;
			my $log = $verbose ? sub { print STDERR "$_[0]\n" } : sub { };
			%testRedun = ();
			while(popAndAdvance($cm, $h, $strset, $results, $log, $testRedunCheck, progressFunc(100))) { $advances++ }
			$advances > 0 || croak("Expected > 0 advances, but got $advances\n");
			my @rk = keys %{$results->results()};
			scalar(@rk) == 1 ||
				croak("Expected 1 result rstr, got ".scalar(@rk)."\n:@rk\n");
			for my $k (keys %{$results->results()}) {
				$k eq "TTACGAACCACAACGTATCG" || croak();
			}
			$results->minCost() == 3 || croak();
			$results->maxCost() == 3 || croak();
			$results->clearAlignments();
		}
	}
	return 1;
}

sub testHeapExpect($$$) {
	my ($heap, $num, $sname) = @_;
	if($heap->get_heap() != $num) {
		print STDERR "Expected $num results in heap, got ".scalar($heap->get_heap())."\n";
		print STDERR "Using search scheme $sname\n";
		print STDERR "Results:\n";
		while($heap->get_heap() > 0) {
			my $r = $heap->pop();
			my $se;
			open($se, ">/dev/stderr") || die;
			$r->print($se);
		}
		croak();
	}
}

sub test10() {
	my @refs = ();
	Reference::newFromHash(
		{ "seq1" => "AGATACTGCGGAGAGGAAAACCCTTCTGAATCGAGCTACGATTTATGCG",
		  "seq2" => "GAGGACCCATATCTACTTACGAACCACAACGTATCGATCGAGTTATGCGTATTCGCGG" },
		\@refs
	);
	#  Ref: TTACGAACCACAACGTATCG
	# Read: TTACGAACCC-AACGTATCG
	#       01234567890123456789
	my $r = Read->new("blah", "TTACGAACCCAACGTATCG", "I"x20, 0, "FW", "?");
	my $h = new Heap::Priority;
	my $results = new ResultAlns(\@refs);
	my $cm = CostModel::bwaDefault();
	for my $strset (Reference::perlHashStrset(\@refs), Reference::perlIndexStrset(\@refs)) {
		for my $ssa (SearchScheme::scheme_one(0, 1, 0), SearchScheme::scheme_one_simple(0, 1, 0)) {
			$h->get_heap() == 0 || die;
			$results->num() == 0 || die;
			my $ss = $ssa->instantiate($r);
			addSearchSchemeToHeap($h, $ss, $r);
			if($ssa->name() eq "OneEditSilly") {
				testHeapExpect($h, 1, $ssa->name());
			} else {
				testHeapExpect($h, 2, $ssa->name());
			}
			my $advances = 0;
			my $log = sub { }; # sub { print STDERR "$_[0]\n" };
			%testRedun = ();
			while(popAndAdvance($cm, $h, $strset, $results, $log, $testRedunCheck, progressFunc(100))) { $advances++ }
			$advances > 0 || croak("Expected > 0 advances, but got $advances\n");
			
			scalar(keys %{$results->results()}) == 1 ||
				croak("Expected 1 result rstr, got ".scalar(keys %{$results->results()})."\n");
			for my $k (keys %{$results->results()}) {
				$k eq "TTACGAACCACAACGTATCG" || croak();
			}
			$results->minCost() == 11 || croak();
			$results->maxCost() == 11 || croak();
			$results->clearAlignments();
		}
	}
	return 1;
}

if($0 =~ /PartialAln\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
	print "Test \"test1\"...";  test1();  print "PASSED\n";
	print "Test \"test2\"...";  test2();  print "PASSED\n";
	print "Test \"test3\"...";  test3();  print "PASSED\n";
	print "Test \"test4\"...";  test4();  print "PASSED\n";
	print "Test \"test5\"...";  test5();  print "PASSED\n";
	print "Test \"test6\"...";  test6();  print "PASSED\n";
	print "Test \"test7\"...";  test7();  print "PASSED\n";
	print "Test \"test8\"...";  test8();  print "PASSED\n";
	print "Test \"test9\"...";  test9();  print "PASSED\n";
	print "Test \"test10\"..."; test10(); print "PASSED\n";
}

1;
