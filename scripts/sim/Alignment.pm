#!/usr/bin/perl -w

##
# Alignment.pm
#
# Read in fasta files containing reference sequences that might be
# aligned to, then read in alignment files, checking each alignment to
# be sure it's sane and consistent with the reference sequence it
# aligns to.
#

package Alignment;
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin); 

##
# Parse a fasta file into the %ref hash
#
sub parseFasta($$) {
	my ($fa, $ref) = shift;
	print STDERR "Parsing FASTA file $fa...\n";
	my $fapipe = "$fa";
	$fapipe = "gzip -dc $fa |" if $fa =~ /\.gz$/;
	open(FA, $fapipe) || die "Could not open '$fapipe' for reading";
	my $name = "";
	my $bases = 0;
	while(<FA>) {
		chomp;
		if(/^>/) {
			$name = substr($_, 1);
			$name =~ s/\s.*//;
			print STDERR "  parsing sequence $name...\n";
			$ref->{std}{$name} = "";
		} else {
			$name ne "" || die "sequence before name";
			$bases += length($_);
			$ref->{std}{$name} .= $_;
		}
	}
	print STDERR "Parsed $bases reference bases\n";
	close(FA);
}

##
# Take all the references in %ref and make both Watson and Crick
# versions where the sequence is in-silico bisulfite treated.
#
sub bisulfiteC($) {
	my $ref = shift;
	for(keys %{$ref->{std}}) {
		$ref->{bisc_fw}{$_} = $ref->{std}{$_};
		$ref->{bisc_fw}{$_} = s/C/Y/g;
		$ref->{bisc_rc}{$_} = $ref->{std}{$_};
		$ref->{bisc_rc}{$_} = s/G/R/g;
	}
}

##
# Take all the references in %ref and make both Watson and Crick
# versions where the sequence is in-silico bisulfite treated.
#
sub bisulfiteCpG($) {
	my $ref = shift;
	for(keys %{$ref->{std}}) {
		$ref->{biscpg_fw}{$_} = $ref->{std}{$_};
		$ref->{biscpg_fw}{$_} =~ s/CG/YG/g;
		$ref->{biscpg_fw}{$_} =~ s/C/T/g;
		$ref->{biscpg_rc}{$_} = $ref->{std}{$_};
		$ref->{biscpg_rc}{$_} =~ s/CG/CR/g;
		$ref->{biscpg_rc}{$_} =~ s/G/A/g;
	}
}

##
# Given a sequence that represents the read oriented s.t. the Watson-
# upstream end is on the left, and given a set of edits to apply to the
# read, also oriented assuming that Watson-upstream is on the left,
# return the string corresponding to the read mutated to match the
# reference.
#
my $nedits = 0;
sub applyEdits($$) {
	my ($seq, $edits) = @_;
	my $rfseq = $seq;
	my $lpos = length($seq)+1;
	$nedits += scalar(@$edits);
	foreach (reverse @$edits) {
		next unless defined($_);
		# Apply the edit
		$_->{pos} <= $lpos || die;
		if($_->{qchr} eq "-") {
			# Insert
			substr($rfseq, $_->{pos}, 0) = $_->{chr};
		} elsif($_->{chr} eq "-") {
			# Deletion
			substr($rfseq, $_->{pos}, 1) eq $_->{qchr} ||
				die "Edit: $_->{pos}:$_->{chr}>$_->{qchr}\n$rfseq";
			substr($rfseq, $_->{pos}, 1) = "";
		} else {
			# Mismatch
			substr($rfseq, $_->{pos}, 1) eq $_->{qchr} ||
				die "Edit: $_->{pos}:$_->{chr}>$_->{qchr}\n$rfseq";
			substr($rfseq, $_->{pos}, 1) = $_->{chr};
		}
	}
	return $rfseq;
}

##
# Given a list of Bowtie edits, invert them by reversing the list and
# changing the poss to be with respect to the other end of the read.
#
sub invertEdits($$) {
	my ($edits, $len) = @_;
	@$edits = reverse @$edits;
	for (@$edits) {
		next unless defined($_);
		defined($_->{qchr}) || die;
		if($_->{qchr} eq "-") {
			$_->{pos} = $len - $_->{pos};
		} else {
			$_->{pos} = $len - $_->{pos} - 1;
		}
	}
}

##
# Given an edit string, parses it into a list of hashes and returns
# the list.
#
sub parseEdits($) {
	my $editstr = shift;
	return undef if (!defined($editstr) || $editstr eq "-" || $editstr eq "");
	my @edits = ();
	# For each edit
	for (split(/,/, $editstr)) {
		# Parse pos
		my ($pos, $ed) = split(/:/);
		defined($ed) || die;
		# Parse chr, qchr
		my ($chr, $qchr) = split(/>/, $ed);
		push @edits, { pos => $pos, chr => $chr, qchr => $qchr };
	}
	return \@edits;
}

##
# Given an edit string, possibly with 4 semicolon-delimited fields,
# parse it into a set of 4 lists of hashes and return the set as an
# array ref.
#
sub parseAllEdits($) {
	my $editstr = shift;
	return [ undef, undef, undef, undef ] if ($editstr eq "-" || $editstr eq "");
	my @editls = split(/;/, $editstr, -1);
	if(scalar(@editls) > 1) {
		scalar(@editls) == 4 || die;
		return [
			parseEdits($editls[0]),
			parseEdits($editls[1]),
			parseEdits($editls[2]),
			parseEdits($editls[3]) ];
	} else {
		scalar(@editls) == 1 || die;
		return [
			parseEdits($editls[0]),
			undef,
			undef,
			undef ];
	}
}

##
# Given a Bowtie orientation string, return true iff the 5' end of the
# read is at the left of the printed sequence.
#
sub fivePrimeLeft($) {
	return ($_[0] eq "+" || $_[0] eq "W" || $_[0] eq "CR");
}

##
# Given the orientation of the read and the state of the global
# bisulfite variables, determine which version of the reference to
# compare against.
#
sub calcRefType($) {
	my $orient = shift;
	if($bisC || $bisCpG) {
		if($orient eq "W" || $orient eq "WR" || $orient eq "+") {
			return $bisC ? "bisc_fw" : "biscpg_fw";
		} else {
			$orient eq "C" || $orient eq "CR" || $orient eq "-" || die;
			return $bisC ? "bisc_rc" : "biscpg_rc";
		}
	} else {
		return "std";
	}
}

##
# Given array refs for two lists of edits, one corresponding to the
# nucleotide edit list and the other corresponding to the resolved
# ambiguous base list, eliminate any edits that appear in both lists.
# Really this shouldn't happen, but I observe that merman does report
# an edit in both categories if the reference base is being resolved to
# an incompatible nucleotide, e.g. R>C.
#
sub removeDups($$) {
	my ($ned, $aed) = @_;
	return unless (defined($ned) && defined($aed));
	for my $i (0..scalar(@$ned)-1) {
		next unless defined($ned->[$i]);
		for my $j (0..scalar(@$aed)-1) {
			next unless defined($aed->[$j]);
			if($ned->[$i]->{qchr} eq $aed->[$j]->{qchr} &&
			   $ned->[$i]->{chr}  eq $aed->[$j]->{chr} &&
			   $ned->[$i]->{pos}  eq $aed->[$j]->{pos})
			{
				#print STDERR "  Eliminated a duplicate edit\n";
				$aed->[$j] = undef;
			}
		}
	}
}

##
# Parse lines from a Bowtie alignment file and check that the
# alignments are sane and consistent with the reference.
#
sub parseBowtie($) {
	my $fh = shift;
	# Print a summary every this many alignments
	my $alIval = 20000;
	while(<$fh>) {
		chomp;
		my ($rdname, $orient, $refname, $off, $seq, $qual, $oms, $editstr) = split(/\t/, $_, -1);
		next if $refname eq "*";
		defined($editstr) || die "Failed to get 8 tokens from line:\n$_";
		$off == int($off) || die "Offset field (col 4) must be an integer:\n$_";
		$oms == int($oms) || die "OMS field (col 7) must be an integer:\n$_";
		my $reftype = calcRefType($orient);
		defined($ref->{$reftype}{$refname}) || die "No such refname as $refname";
		my $edits4 = parseAllEdits($editstr);
		my ($ned, $aed) = ($edits4->[0], $edits4->[1]);
		removeDups($ned, $aed);
		my $fpl = fivePrimeLeft($orient);
		invertEdits($ned, length($seq)) unless ($fpl || !defined($ned));
		invertEdits($aed, length($seq)) unless ($fpl || !defined($aed));
		my $rfseq = $seq;
		$rfseq = applyEdits($rfseq, $ned) if defined($ned);
		$rfseq = applyEdits($rfseq, $aed) if defined($aed);
		my $trueRfseq = substr($ref->{$reftype}{$refname}, $off, length($rfseq));
		$rfseq eq $trueRfseq || die "Did not match:\n$seq\n$rfseq\n$trueRfseq\n";
		$nals++;
		if(($nals % $alIval) == 0) {
			print STDERR "  Checked $nals alignments, $nedits edits\n";
		}
	}
}

##
# Parse lines from a SAM alignment file and check that alignments are
# sane and consistent with the reference.
#
sub parseSam($) {
	while(<$_[0]>) {
		chomp;
		my ($rdname, $orient, $refname, $off, $seq, $qual, $oms, $edits) = split(/\t/);
		defined($edits) || die;
	}
}

##
#
#
sub checkAlignments() {
	my ($fas, $als, $parser, $bisC, $bisCpG) = @_;
	my $nals = 0;
	my %ref = ();
	foreach (@$fas) { parseFasta($_, \%ref); }
	if($bisC) {
		print STDERR "Generating all-C bisulfite-treated references\n";
		bisulfiteC(\%ref);
	}
	if($bisCpG) {
		print STDERR "Generating all-CpG bisulfite-treated references\n";
		bisulfiteCpG(\%ref);
	}
	foreach (@$als) {
		my $alnpipe = $_;
		print STDERR "Processing alignment file '$_'\n";
		$alnpipe = "gzip -dc $_ |" if $_ =~ /\.gz$/;
		my $alnfh = undef;
		open($alnfh, $alnpipe) || die "Could not open '$alnpipe' for reading";
		$parser->($alnfh);
		close($alnfh);
	}

	print STDERR "--- Summary ---\n";
	print STDERR "Read ".scalar(keys %ref)." reference strings\n";
	print STDERR "Checked $nals alignments, $nedits edits\n";
	print STDERR "---------------\n";
	print STDERR "PASSED\n";
}

if($0 =~ /^Alignment\.pm$/) {

	my @fas = ();
	my @als = ();
	my $format = "bowtie";
	my $bisC = 0;
	my $bisCpG = 0;

	GetOptions (
		"fasta|ref=s"           => \@fas,
		"als|aln=s"             => \@als,
		"format=s"              => \$format,
		"bis-C|bisulfite-C"     => \$bisC,
		"bis-CpG|bisulfite-CpG" => \$bisCpG,
		"bowtie"                => sub {$format = "bowtie"},
		"sam"                   => sub {$format = "sam"}) || die;

	my %parsers = (
		"bowtie" => \&parseBowtie,
		"sam"    => \&parseSam);
	my $parser = $parsers{$format};
	defined($parser) || die "No parser for format '$format'";
	
	checkAlignments(\@fas, $
		my ($fas, $parser, $bisC, $bisCpG) = @_;
}

1;
