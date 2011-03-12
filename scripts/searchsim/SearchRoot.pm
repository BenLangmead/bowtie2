#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 9/5/2009
#
# Refsim encapsulates a reference sequence along with several known-SNP
# sites.  The routine simulateSubject allows the caller to simuate a
# new reference with some known SNPs and some novel SNPs.
#

package SearchRoot;
use strict;
use Carp;
use FindBin qw($Bin); 
use lib $Bin;
use Edit;
use Read;
use Test;
use Constraint;

# A search root defines a particular set of paths through the search
# space for a read.  One or more search roots together might constitute
# an exhaustive search w/r/t a particular alignment policy; e.g., two
# roots can be combined to make an efficient and exhaustive 1-edit
# policy.

# A set of search roots may or may not define a set of paths s.t. some
# paths are redundant with other paths.  I.e. two paths corresponding
# to the exact same alignment but occurring in two different roots.

# What happens when the read is shorter than 1 seed?  And what happens
# when the final seed doesn't fit entirely on the read?

sub new {
	scalar(@_) == 13 || croak("Expected 13 arguments, got ".scalar(@_)."\n");
	my ($class, $name, $scheme, $fw, $off, $len, $rep, $zone_str,
	    $cons, $ocons, $expands, $retracts, $instantiated) = @_;
	defined($off) || croak("SearchRoot->new: must specify offset");
	(defined($len) && $len > 0) || croak("SearchRoot->new: must specify length >0");
	defined($zone_str) || croak("No zone string");
	defined($cons)  || ($cons = { });
	defined($ocons) || ($ocons = Constraint::none());
	my $self = bless {
		_name   => $name   || "noname",
		_scheme => $scheme,
		_fw     => $fw     || croak("No fw"),
		# _off: offset from 5' end of read of first character in the
		# seed
		_off => $off,
		# _len: length of the seed
		_len => $len,
		# _rep: repeat this root every 'rep' bases
		_rep => $rep,
		# _zone_str: defines the constraint zones
		_zone_str => $zone_str,
		# _cons: defines the zone constraints
		_cons => $cons || croak("No zone constraints"),
		# _ocons: defines the overall constraints
		_ocons => $ocons || croak("No overall constraints"),
		# If a seed expands, then a read larger than the seed causes
		# the _zone_str to be lengthened to be the same length as the
		# read.  This is useful for implementing, e.g., an end-to-end
		# 1-edit policy.
		_expands => $expands || 0,
		# If a seed retracts, then a read that fails to entirely
		# overlap the seed causes the _zone_str to be shortened to be
		# the same length as the read.  If a seed does not retract,
		# then a read that fails to entirely overlap the seed causes
		# the seed to be excluded from consideration.
		_retracts => $retracts || 0,
		# When a root is instantiated (i.e. tailored to a particular
		# input read), this becomes 1
		_instantiated => (defined($instantiated) ? $instantiated : 1)
	}, $class;
	$self->ocons eq "?" && croak();
	$self->ocons->isa("Constraint") || croak();
	for my $c (values %{$self->cons}) {
		$c->isa("Constraint") || croak();
	}
	return $self;
}
sub name         { return $_[0]->{_name}         }
sub scheme       { return $_[0]->{_scheme}       }
sub fw           { return $_[0]->{_fw}           }
sub off          { return $_[0]->{_off}          }
sub len          { return $_[0]->{_len}          }
sub rep          { return $_[0]->{_rep}          }
sub zone_str     { return $_[0]->{_zone_str}     }
sub cons         { return $_[0]->{_cons}         }
sub ocons        { return $_[0]->{_ocons}        }
sub expands      { return $_[0]->{_expands}      }
sub retracts     { return $_[0]->{_retracts}     }
sub instantiated { return $_[0]->{_instantiated} }

##
# Helper that takes a zone string and resizes it to match a new seed
# length.
#
sub resizeZoneStr($$) {
	my ($str, $len) = @_;
	while($len < length($str)) {
		# Remove characters
		if((length($str) % 2) == 0) {
			# Even, so remove from right
			$str = substr($str, 0, -1);
		} else {
			# Odd, so remove from left
			$str = substr($str, 1);
		}
	}
	while($len > length($str)) {
		# Add characters
		if((length($str) % 2) == 0) {
			# Even, so add to left
			$str = substr($str, 0, 1) . $str;
		} else {
			# Odd, so add to right
			$str = $str . substr($str, -1, 1);
		}
	}
	return $str;
}

##
# Test whether a SearchRoot is internally consistent.
#
sub repOkTest($) {
	my $self = shift;
	$self->len() > 0 || return 0;
	$self->off() >= 0 || return 0;
}

##
# Die if SearchRoot isn't internally consistent.
#
sub repOk($) { $_[0]->repOkTest() || die }

##
# Instatiate a SearchRoot with respect to a particular read.  This may
# result in multiple instantiated SearchRoots, if 'rep' is set.
#
sub instantiate($$) {
	my ($self, $read) = @_;
	$self->instantiated() && croak("Tried to instantiate an instantiated SearchRoot");
	$self->repOk();
	my @iroots = ();
	# Parameters for the Root we'll eventually instantiate
	my ($len, $zone_str, $cons, $ocons) =
		($self->len(), $self->zone_str(), $self->cons(),
		 $self->ocons()->instantiate($read));
	my $rep = $self->rep;
	$rep = $read->len if $rep == 0;
	for(my $off = $self->off; $off < $read->len; $off += $rep) {
		my $retracts = ($self->retracts && $off == $self->off);
		# Does the seed fall entirely within the bounds of the read?
		if($retracts && $read->len < $off + $len) {
			# Retract the instantiated seed
			$len = $read->len - $off;
		}
		if($self->expands && $read->len > $off + $len) {
			# Expand the instantiated seed
			$len = $read->len - $self->off;
		}
		next if $off + $len > $read->len;
		$zone_str = resizeZoneStr($zone_str, $len);
		push @iroots, new SearchRoot(
			$self->name, undef, "FW", $off, $len, 0, $zone_str, $cons,
			$ocons, $self->expands, $retracts, 1);
	}
	return @iroots;
}

sub test1() {
	my $s = SearchRoot->new(
		"blah", undef, "FW", 0, 28, 0, "01",
		{ 0 => Constraint::none(), 1 => Constraint::none() },
		Constraint::none(),
		1, 1, 0);
	$s->repOkTest() || return 0;
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my @ss2 = $s->instantiate($r);
	scalar(@ss2) == 1 || croak("Expected 1 instantiated root, got ".scalar(@ss2));
	my $s2 = $ss2[0];
	$s2->off() == 0 || die "off was ".$s2->off()."\n";
	$s2->len() == 30 || die "len was ".$s2->len()."\n";
	my $ex_zone_str = ("0"x15).("1"x15);
	$s2->zone_str() eq $ex_zone_str || die "expected $ex_zone_str, got ".$s2->zone_str()."\n";
	return 1;
}

sub test2() {
	my $s = SearchRoot->new(
		"blah", undef, "FW", 0, 32, 0, "01",
		{ 0 => Constraint::none(), 1 => Constraint::none() },
		Constraint::none(),
		1, 1, 0);
	$s->repOkTest() || return 0;
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my @ss2 = $s->instantiate($r);
	scalar(@ss2) == 1 || croak("Expected 1 instantiated root, got ".scalar(@ss2));
	my $s2 = $ss2[0];
	$s2->off() == 0 || die "off was ".$s2->off()."\n";
	$s2->len() == 30 || die "len was ".$s2->len()."\n";
	my $ex_zone_str = ("0"x15).("1"x15);
	$s2->zone_str() eq $ex_zone_str || die "expected $ex_zone_str, got ".$s2->zone_str()."\n";
	return 1;
}

sub test3() {
	my $s = SearchRoot->new(
		"blah", undef, "FW", 0, 32, 0, "01",
		{ 0 => Constraint::none(), 1 => Constraint::none() },
		Constraint::none(),
		0, 0, 0);
	$s->repOkTest() || return 0;
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my @ss2 = $s->instantiate($r);
	scalar(@ss2) == 0 || croak();
	return 1;
}

sub test4() {
	my $s = SearchRoot->new(
		"blah", undef, "FW", 0, 28, 0, "01",
		{ 0 => Constraint::none(), 1 => Constraint::none() },
		Constraint::none(),
		0, 0, 0);
	$s->repOkTest() || return 0;
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my @ss2 = $s->instantiate($r);
	scalar(@ss2) == 1 || croak();
	my $s2 = $ss2[0];
	defined($s2) || die;
	$s2->off() == 0 || die "off was ".$s2->off()."\n";
	$s2->len() == 28 || die "len was ".$s2->len()."\n";
	my $ex_zone_str = ("0"x14).("1"x14);
	$s2->zone_str() eq $ex_zone_str || die "expected $ex_zone_str, got ".$s2->zone_str()."\n";
	return 1;
}

sub test5() {
	my $s = new SearchRoot(
		"blah", undef, "FW", 0, 30, 1, "01",
		{ 0 => Constraint::none(), 1 => Constraint::none() },
		Constraint::none(),
		0, 1, 0);
	$s->repOkTest() || return 0;
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my @ss2 = $s->instantiate($r);
	scalar(@ss2) == 1 || croak("Expected 1 instantiated root, got ".scalar(@ss2));
	$ss2[0]->off() == 0 || die "off was ".$ss2[0]->off()."\n";
	$ss2[0]->len() == 30 || die "len was ".$ss2[0]->len()."\n";
	$ss2[0]->retracts || die;
	return 1;
}

sub test6() {
	my $s = SearchRoot->new(
		"blah", undef, "FW", 0, 10, 1, "01",
		{ 0 => Constraint::none(), 1 => Constraint::none() },
		Constraint::none(),
		0, 1, 0);
	$s->repOkTest() || return 0;
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my @ss2 = $s->instantiate($r);
	scalar(@ss2) == 21 || croak();
	for(my $i = 0; $i < scalar(@ss2); $i++) {
		$ss2[$i]->off() == $i || die "off was ".$ss2[$i]->off()."\n";
		$ss2[$i]->len() == 10 || die "len was ".$ss2[$i]->len()."\n";
	}
	$ss2[0]->retracts || die;
	!$ss2[1]->retracts || die;
	return 1;
}

sub test7() {
	my $s = SearchRoot->new(
		"blah", undef, "FW", 0, 10, 2, "01",
		{ 0 => Constraint::none(), 1 => Constraint::none() },
		Constraint::none(),
		0, 1, 0);
	$s->repOkTest() || return 0;
	my $r = Read->new("blah", "ACGTA"x6, "IIIII"x6, 0, "FW", "?");
	my @ss2 = $s->instantiate($r);
	scalar(@ss2) == 11 || croak();
	for(my $i = 0; $i < scalar(@ss2); $i++) {
		$ss2[$i]->off() == $i*2 || die "off was ".$ss2[$i]->off()."\n";
		$ss2[$i]->len() == 10 || die "len was ".$ss2[$i]->len()."\n";
	}
	$ss2[0]->retracts || die;
	!$ss2[1]->retracts || die;
	return 1;
}

if($0 =~ /SearchRoot\.pm$/) {
	print "Running unit tests\n";
	# Run unit tests
	print "Test \"test1\"..."; test1(); print "PASSED\n";
	print "Test \"test2\"..."; test2(); print "PASSED\n";
	print "Test \"test3\"..."; test3(); print "PASSED\n";
	print "Test \"test4\"..."; test4(); print "PASSED\n";
	print "Test \"test5\"..."; test5(); print "PASSED\n";
	print "Test \"test6\"..."; test6(); print "PASSED\n";
	print "Test \"test7\"..."; test7(); print "PASSED\n";
}

1;
