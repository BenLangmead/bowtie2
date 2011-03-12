#!/usr/bin/perl -w

##
# Test ideas for how to extend partial, gapped alignments such that no
# redundant alignments are considered.
#

use strict;
use warnings;

my $ref;  # Current reference
my $read; # Current read

my $nrefs = 3;
my $nreadsPerRef = 10;

my $refLenLo = 200;
my $refLenDelt = 100;

my $readLenLo = 40;
my $readLenDelt = 30;

my $mmsLo = 0;
my $mmsDelt = 7;

my $gapsLo = 0;
my $gapsDelt = 7;

my $mmpen = 30;
my $rdOpenPen = 40;
my $rfOpenPen = 40;
my $rdExtendPen = 15;
my $rfExtendPen = 15;
my $maxPenPerPos = 3;
my $maxPen;

my $userref  = shift @ARGV;
my $userread = shift @ARGV;

sub generateRef($) {
	my ($len) = @_;
	if(defined($userref)) {
		$ref = $userref;
		return;
	}
	$ref = "";
	for(my $i = 0; $i < $len; $i++) {
		my $c = substr("ACGT", int(rand(4)), 1);
		$ref .= $c;
	}
	length($ref) == $len || die;
}

sub generateRead($$$) {
	my ($len, $mms, $gaps) = @_;
	if(defined($userread)) {
		$read = $userread;
		return;
	}
	$read = substr($ref, int(rand(length($ref) - $len)), $len);
	length($read) == $len || die "Length: ".length($read).", len: $len\n";
	for(my $i = 0; $i < $mms; $i++) {
		substr($read, int(rand(length($read))), 1) = substr("ACGT", int(rand(4)), 1);
	}
	for(my $i = 0; $i < $gaps; $i++) {
		if(int(rand(2)) == 0) {
			# Insertion
			substr($read, int(rand(length($read))), 0) = substr("ACGT", int(rand(4)), 1);
		} else {
			# Deletion
			substr($read, int(rand(length($read))), 1) = "";
		}
	}
}

sub nextDna($) {
	my $c = shift;
	return 'C' if $c eq 'A';
	return 'G' if $c eq 'C';
	return 'T' if $c eq 'G';
	return 'A' if $c eq 'T';
	die;
}

# Map from partial alignment keys (combination of $rdoff, $rfoff, and
# $parrf with gaps removed) to $parrd, $parrf, $pen triplets.
my %parmap = ();
my %parredundant = ();
my %pathmap = ();

sub extendPartialAlignment($$$$$);

##
# Recursive function that takes a partial alignment and extends it in
# all possible ways.
#
sub extendPartialAlignment($$$$$) {
	my ($parrd, $parrf, $rdoff, $rfoff, $pen) = @_;
	defined($pathmap{"$parrd$parrf"}) && die "Visited more than once:\n$parrd\n$parrf\n";
	$pathmap{"$parrd$parrf"} = 1;
	if($rdoff == length($read) || $rfoff == length($ref)) {
		print "Finished one; penalty=$pen:\n$parrf\n$parrd\n";
		return;
	}
	# Bail if 
	my $parrfGapless = $parrf;
	$parrfGapless =~ tr/-//;
	return if index($ref, $parrfGapless) == -1;
	my $key = "$parrfGapless,$rdoff,$rfoff";
	my $val = "$parrd,$parrf,$pen";
	if(defined($parmap{$key})) {
		# Redundant
		$parmap{$key} .= ";$val";
		$parredundant{$key} = 1;
	} else {
		$parmap{$key} = $val;
	}
	my $lastrdc = substr($parrd, -1, 1);
	my $lastrfc = substr($parrf, -1, 1);
	my $rdc = substr($read, $rdoff, 1);
	my $rfc = substr($ref,  $rfoff, 1);
	length($rdc) == 1 || die;
	length($rfc) == 1 || die;
	# Try introducing a gap in the read
	my $rdIsOpen = ($lastrdc ne '-' && $lastrdc ne '=');
	my $rdgpc = $rdIsOpen ? '=' : '-';
	my $rdgp = $rdIsOpen ? $rdOpenPen : $rdExtendPen;
	my $rfsame = $rfc eq $lastrfc;
	$rfsame = 0 if $rdoff == 1;
	if(!$rfsame && $pen + $rdgp <= $maxPen) {
		extendPartialAlignment($parrd.$rdgpc, $parrf.$rfc, $rdoff, $rfoff+1, $pen + $rdgp);
	}
	# Try introducing a gap in the reference
	my $rfIsOpen = ($lastrfc ne '-' && $lastrfc ne '=');
	my $rfgpc = $rfIsOpen ? '=' : '-';
	my $rfgp = $rfIsOpen ? $rfOpenPen : $rfExtendPen;
	if($pen + $rfgp <= $maxPen) {
		extendPartialAlignment($parrd.$rdc, $parrf.$rfgpc, $rdoff+1, $rfoff, $pen + $rfgp);
	}
	# Try introducing a mismatch; try each separately
	if($pen + $mmpen <= $maxPen) {
		my $mmc = $rdc;
		for(my $i = 0; $i < 3; $i++) {
			$mmc = nextDna($mmc);
			extendPartialAlignment($parrd.$rdc, $parrf.$mmc, $rdoff+1, $rfoff+1, $pen + $mmpen);
		}
	}
	extendPartialAlignment($parrd.$rdc, $parrf.$rdc, $rdoff+1, $rfoff+1, $pen);
}

sub alignRead() {
	# Force first read char to align to first ref char
	extendPartialAlignment(substr($read, 0, 1), substr($ref, 0, 1), 1, 1, 0);
}

if(defined($userref) && defined($userread)) {
	$nrefs = 1;
	$nreadsPerRef = 1;
}

for(my $i = 0; $i < $nrefs; $i++) {
	generateRef(int(rand($refLenDelt)) + $refLenLo);
	for(my $j = 0; $j < $nreadsPerRef; $j++) {
		generateRead(int(rand($readLenDelt)) + $readLenLo,
		             int(rand($mmsDelt))     + $mmsLo,
		             int(rand($gapsDelt))    + $gapsLo);
		$maxPen = length($read) * $maxPenPerPos;
		print "$ref\n$read\n";
		alignRead();
		# Now pick out all the equivalence classes of alignments
		for my $k (keys %parredundant) {
			print "Redundant alignments for $k:\n";
			for my $al (split(/;/, $parmap{$k})) {
				my ($parrd, $parrf, $pen) = split(/,/, $al);
				print " Ref: $parrf\n";
				print "Read: $parrd\n";
				print " Pen: $pen\n\n";
			}
		}
		%parredundant = ();
		%parmap = ();
		%pathmap = ();
	}
}
