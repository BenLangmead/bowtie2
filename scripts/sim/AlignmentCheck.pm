#!/usr/bin/perl -w

#
# Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
#
# This file is part of Bowtie 2.
#
# Bowtie 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Bowtie 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
#

##
# AlignmentCheck.pm
#
# Read in fasta files containing reference sequences that might be
# aligned to, then read in alignment files, checking each alignment to
# be sure it's sane and consistent with the reference sequence it
# aligns to.
#

package AlignmentCheck;
use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin;
use DNA;
use Data::Dumper;

##
# Parse a fasta file into the %ref hash
#
sub parseFasta($$) {
	my ($fa, $ref) = @_;
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
# Create a new alignment checker
#
sub new {
	my (
		$class,
		$name,   # name of checker
		$fas,    # list of fasta files containing reference sequences
		         # or hash of the references themselves (former if it's
		         # an array ref latter if it's a hash ref)
		$format, # alignment format
		$bisC,   # whether alignment was w/r/t ref w/ all Cs converted to Ys
		$bisCpG, # whether alignment was w/r/t ref w/ all CpGs converted to YpGs
	) = @_;
	(defined($fas)) || die "Must specify non-empty list of fasta files";
	# Set defaults for format, bisC, bisCpG, name
	$format = "bowtie" unless defined($format);
	$bisC   = 0 unless defined($bisC);
	$bisCpG = 0 unless defined($bisCpG);
	$name = "noname" unless defined($name);
	# Parse all the fasta files into the ref hash
	my %ref = ();
	if(ref($fas) eq "HASH") {
		for (keys %$fas) {
			$ref{std}{$_} = $fas->{$_}
		}
	} else {
		ref($fas) eq "ARRAY" || die;
		foreach (@$fas) { parseFasta($_, \%ref); }
	}
	return bless {
		_name       => $name,
		_fas        => $fas,
		_format     => $format,
		_bisC       => $bisC,
		_bisCpG     => $bisCpG,
		_refs       => \%ref,
		_nals       => 0,
		_nedits     => 0
	}, $class;
}
sub name      { return $_[0]->{_name}      }
sub fas       { return $_[0]->{_fas}       }
sub format    { return $_[0]->{_format}    }
sub bisC      { return $_[0]->{_bisC}      }
sub bisCpG    { return $_[0]->{_bisCpG}    }
sub refs      { return $_[0]->{_refs}      }
sub nals      { return $_[0]->{_nals}      }
sub nedits    { return $_[0]->{_nedits}    }
sub nrefs     { return scalar(keys %{$_[0]->{_refs}}); }

##
# Given a sequence that represents the read oriented s.t. the Watson-
# upstream end is on the left, and given a set of edits to apply to the
# read, also oriented assuming that Watson-upstream is on the left,
# return the string corresponding to the read mutated to match the
# reference.
#
my $nedits = 0;
sub applyEdits($$$) {
	my ($seq, $edits, $line) = @_;
	my $rfseq = $seq;
	my $lpos = length($seq)+1;
	$nedits += scalar(@$edits);
	foreach (reverse @$edits) {
		next unless defined($_);
		#print STDERR "Applying edit at $_->{pos}\n";
		# Apply the edit
		$_->{pos} <= $lpos || die "Edit position $_->{pos} was greater than previous $lpos";
		if($_->{qchr} eq "-") {
			# Insert
			$_->{pos} <= length($rfseq) || die "Read gap pos $_->{pos} was not <= string len ".length($rfseq)."\n$line";
			substr($rfseq, $_->{pos}, 0) = $_->{chr};
		} elsif($_->{chr} eq "-") {
			# Deletion
			$_->{pos} < length($rfseq) || die "Ref gap pos $_->{pos} was not < string len ".length($rfseq)."\n$line";
			my $dc = substr($rfseq, $_->{pos}, 1);
			$dc eq $_->{qchr} ||
				die "Edit: $_->{pos}:$_->{chr}>$_->{qchr} but ref char was $dc".
				    "\n$rfseq\n$line";
			substr($rfseq, $_->{pos}, 1) = "";
		} else {
			# Mismatch
			$_->{pos} < length($rfseq) || die "Mismatch pos $_->{pos} was not < string len ".length($rfseq)."\n$line";
			substr($rfseq, $_->{pos}, 1) eq $_->{qchr} ||
				die "Edit: $_->{pos}:$_->{chr}>$_->{qchr}\n$rfseq\n$line";
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
			length($_->{chr}) >= 1 || die;
			$_->{chr} = reverse $_->{chr};
		} else {
			$_->{pos} = $len - $_->{pos} - 1;
			length($_->{chr}) == 1 || die;
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
sub calcRefType {
	my ($self, $orient) = @_;
	if($self->bisC || $self->bisCpG) {
		if($orient eq "W" || $orient eq "WR" || $orient eq "+") {
			return $self->bisC ? "bisc_fw" : "biscpg_fw";
		} else {
			$orient eq "C" || $orient eq "CR" || $orient eq "-" || die;
			return $self->bisC ? "bisc_rc" : "biscpg_rc";
		}
	} else {
		return "std";
	}
}

##
# Parse a CIGAR string into parallel arrays of CIGAR operations (M, I, D)
#
sub cigarParse($$$) {
	my ($cigar, $ops, $runs) = @_;
	my $i = 0;
	while($i < length($cigar)) {
		substr($cigar, $i) =~ /^([0-9]+)/;
		defined($1) || die "Could not parse number at pos $i: '$cigar'";
		$i += length($1);
		$i < length($cigar) || die;
		push @$runs, $1;
		my $op = substr($cigar, $i, 1);
		defined($op) || die "Could not parse operation at pos $i: '$cigar'";
		push @$ops, $op;
		$i++;
	}
}

##
# Trim a read sequence according to the soft clipping in the CIGAR string.
#
sub cigarTrim($$) {
	my ($seq, $cigar) = @_;
	my @ops = ();
	my @runs = ();
	cigarParse($cigar, \@ops, \@runs);
	my ($trimup, $trimdn) = (0, 0);
	if($ops[0] eq 'S') {
		$runs[0] < length($seq) || die "Soft clipped entire alignment!";
		$seq = substr($seq, $runs[0]);
		$trimup = $runs[0];
	}
	if(scalar(@ops) > 1 && $ops[-1] eq 'S') {
		$runs[-1] < length($seq) || die "Soft clipped entire alignment!";
		$seq = substr($seq, 0, -$runs[-1]);
		$trimdn = $runs[-1];
	}
	return ($seq, $trimup, $trimdn);
}

##
# Parse a CIGAR string into a string of operators.  Operators are expanded into
# runs where appropriate.  = and X are collapsed into M.
#
sub parse_cigar($) {
	my ($cigar) = @_;
	my $ret = "";
	my $i = 0;
	my ($rdlen, $rflen) = (0, 0);
	while($i < length($cigar)) {
		substr($cigar, $i) =~ /^([0-9]+)/;
		defined($1) || die "Could not parse number at pos $i: '$cigar'";
		my $runlen = $1;
		$i += length($1);
		$i < length($cigar) || die;
		my $op = substr($cigar, $i, 1);
		defined($op) || die "Could not parse operation at pos $i: '$cigar'";
		if($op eq "X" || $op eq "=") {
			$op = "M";
		}
		$rdlen += $runlen if $op ne "D";
		$rflen += $runlen if $op ne "I";
		$ret .= ($op x $runlen);
		$i++;
	}
	return ($ret, $rdlen, $rflen);
}

##
# Parse an MD:Z string into a string with length equal to query length.  Each
# position contains either a space, if the read matches the reference at that
# position, or a character, if the reference contains a character that doesn't
# match its opposite in the alignment.  In the latter case, the character in
# the string is the reference character.
#
sub parse_md($) {
	my ($md) = @_;
	my $i = 0;
	my $ret = "";
	while($i < length($md)) {
		# Starts with a number?
		my $ch = substr($md, $i, 1);
		if($ch =~ /[0-9]/) {
			# Parse the number off the beginning
			substr($md, $i) =~ /^([0-9]+)/;
			defined($1) || die "Could not parse number at pos $i: '$md'";
			my $runlen = $1;
			$ret .= (" " x $runlen) if $runlen > 0;
			$i += length($runlen);
		} elsif($ch eq "^") {
			# Read gap
			$i++;
			substr($md, $i) =~ /^([A-Za-z]+)/;
			defined($1) || die "Could not parse read gap at pos $i: '$md'";
			my $chrs = $1;
			$i += length($chrs);
			$ret .= $chrs;
		} else {
			# DNA character
			$ch =~ /[ACGTN.]/i || die "Bad char '$ch' at pos $i: '$md'";
			$ret .= $ch;
			$i++;
		}
	}
	return $ret;
}

##
# Given a read sequence (with characters in upstream-to-downstream order with
# respect to the reference - NOT necessarily 5'-to-3') and a CIGAR string and
# an MD:Z string, build the alignment strings.  The alignment strings will only
# contain the portion of the read that aligned.  Any portions that were either
# hard-trimmed or soft-trimmed are trimmed from this function's result.
#
# For now, I'm assuming that the MD:Z string only describes aligned characters,
# i.e. *after* trimming.
#
sub _read_md_cigar_to_al($$$) {
	my ($seq, $md, $cigar) = @_;
	my $alread = "";
	my $alref = "";
	$cigar ne "*" || die "CIGAR string was star!";
	$seq ne "" || die "Empty sequence given to _read_md_cigar_to_al";
	my $parsed_md  = parse_md($md);
	my ($parsed_cig, $cigar_rdlen, $cigar_rflen) = parse_cigar($cigar);
	my ($rdoff, $mdoff) = (0, 0);
	my ($htriml, $striml, $htrimr, $strimr) = (0, 0, 0, 0);
	my $nonsh = 0; # have I seen a non-S, non-H CIGAR op?
	my $nonh = 0;  # have I seen a non-H CIGAR op?
	for(my $i = 0; $i < length($parsed_cig); $i++) {
		my $cigop = substr($parsed_cig, $i, 1);
		$nonh++ if $cigop ne "H";
		$nonsh++ if ($cigop ne "H" && $cigop ne "S");
		if($cigop eq "S") {
			if($nonsh) {
				$strimr++;
			} else {
				$striml++;
			}
			$rdoff++;
			next;
		}
		if($cigop eq "H") {
			if($nonh) {
				$htrimr++;
			} else {
				$htriml++;
			}
			next;
		}
		$cigop = "M" if $cigop eq "=" || $cigop eq "X";
		if($cigop eq "P") {
			# Padding
			$alread .= "-";
			$alref .= "-";
		} elsif($cigop eq "M") {
			my $rdc = substr($seq, $rdoff, 1);
			$mdoff < length($parsed_md) ||
				die "Bad mdoff ($mdoff)\nlength(parsed_md)=".length($parsed_md)."\nseq:\n$seq\ncigar:\n$cigar\nmd:\n$md\nparsed md:\n$parsed_md";
			my $rfc = substr($parsed_md, $mdoff, 1);
			$rfc = $rdc if $rfc eq " ";
			$alread .= $rdc;
			$alref .= $rfc;
			$rdoff++;
			$mdoff++;
		} elsif($cigop eq "D") {
			# Read gap
			#  Read: AAA-AAA
			#   Ref: AAAAAAA
			my $rfc = substr($parsed_md, $mdoff, 1);
			$rfc ne " " ||
				die "Must have a ref character opposite a gap in the read:\n".
				    "cig: $parsed_cig ($i)\nmd:  $parsed_md ($mdoff)\n";
			$alread .= "-";
			$alref .= $rfc;
			$mdoff++;
		} else {
			# Reference gap
			#  Read: AAAAAAA
			#   Ref: AAA-AAA
			$cigop eq "I" || die "Unsupported cigop: $cigop in cigar $cigar";
			my $rdc = substr($seq, $rdoff, 1);
			$alread .= $rdc;
			$alref .= "-";
			$rdoff++;
			# $mdoff SHOULD NOT be incremented here
		}
		$rdoff <= length($seq) ||
			die "Bad rdoff:$rdoff for seq '$seq' cigop=$cigop\nseq: $seq\ncigar=$cigar\nmd=$md";
	}
	return ($alread, $alref, $htriml, $striml, $htrimr, $strimr);
}

##
# Parse a line from a Bowtie alignment file and check that the
# alignment is sane and consistent with the reference.
#
sub parseBowtieLines {
	my ($self, $lines) = @_;
	for my $line (@$lines) {
		chomp($line);
		my ($rdname, $orient, $refname, $off, $seq, $qual, $oms, $editstr,
		    $flags) = split(/\t/, $line, -1);
		next if $refname eq "*";
		$flags =~ /XC:([^,\s]+)/;
		my $cigar = $1;
		defined($cigar) ||
			die "Could not parse CIGAR string from flags: '$flags'";
		defined($editstr) || die "Failed to get 8 tokens from line:\n$_";
		$off == int($off) || die "Offset field (col 4) must be an integer:\n$_";
		$oms == int($oms) || die "OMS field (col 7) must be an integer:\n$_";
		my $reftype = $self->calcRefType($orient);
		defined($self->refs->{$reftype}{$refname}) ||
			die "No such refname as $refname for reftype $reftype:\n".
			Dumper($self->refs->{$reftype});
		my $edits4 = parseAllEdits($editstr);
		my ($ned, $aed) = ($edits4->[0], $edits4->[1]);
		removeDups($ned, $aed);
		my $fpl = fivePrimeLeft($orient);
		# Trim seq according to CIGAR string
		my $rfseq = $seq;
		my ($trimup, $trimdn);
		($rfseq, $trimup, $trimdn) = cigarTrim($rfseq, $cigar);
		invertEdits($ned, length($rfseq)) unless ($fpl || !defined($ned));
		invertEdits($aed, length($rfseq)) unless ($fpl || !defined($aed));
		$rfseq = applyEdits($rfseq, $ned, $line) if defined($ned);
		$rfseq = applyEdits($rfseq, $aed, $line) if defined($aed);
		# Check if our alignment falls off the end of the reference, in
		# which case we need to pad the reference string with Ns
		my $exoff = $off;
		my $padleft = "";
		my $exlen = length($rfseq);
		my $tlen = length($self->refs->{$reftype}{$refname});
		if($exoff < 0) {
			# Alignment hangs off LHS; pad it
			my $npad = -$exoff;
			$padleft = "N" x $npad;
			$exlen += $exoff;
			$exlen >= 0 ||
				die "Read was entirely off the LHS of the reference\n".
					"Referemce len=$tlen\n".
					"Alignment referemce len=$tlen\n".
					"$line\n";
			$exoff = 0;
		}
		my $padright = "";
		my $roverhang = $off + length($rfseq) - $tlen;
		if($roverhang > 0) {
			$padright = "N" x $roverhang;
			$exlen -= $roverhang;
			$exlen >= 0 ||
				die "Read was entirely off the RHS of the reference\n".
					"Referemce len=$tlen\n".
					"Alignment referemce len=$tlen\n".
					"$line\n";
		}
		my $refsub = substr($self->refs->{$reftype}{$refname}, $exoff, $exlen);
		length($refsub) == $exlen ||
			die "Tried to extract ref substring of length $exlen, got ".
			    "\"$refsub\" from \"".$self->refs->{$reftype}{$refname}."\"".
				"\n$line\n".
				"\noff=$off, rfseq=$rfseq\n";
		$refsub = DNA::iupacSubN($refsub);
		my $trueRfseq = $padleft . $refsub . $padright;
		length($trueRfseq) == length($rfseq) ||
			die "Different lengths for edited read and ref:\n".
			"       Read: $seq\n".
			"Edited read: $rfseq\n".
			"        Ref: $trueRfseq\n";
		$rfseq eq $trueRfseq ||
			die "Did not match:\n".
			"       Read: $seq\n".
			"Edited read: $rfseq\n".
			"        Ref: $trueRfseq\n";
		$self->{_nals}++;
	}
}

##
# Parse a line from a SAM alignment file and check that the
# alignment is sane and consistent with the reference.
#
sub parseSamLines {
	my ($self, $lines) = @_;
	my ($lastseq, $lastqual) = ("", "");
	my $idx = 0;
	for my $line (@$lines) {
		$idx++;
		print STDERR "Processing line...\n";
		chomp($line);
		next if $line =~ /^\@/;
		my @toks = split(/\t/, $line, -1);
		my (
			$qname, #1
			$flag,  #2
			$rname, #3
			$pos,   #4
			$mapq,  #5
			$cigar, #6
			$rnext, #7
			$pnext, #8
			$tlen,  #9
			$seq,   #10
			$qual) = @toks;
		defined($qual) || die "Not enough SAM tokens:\n$line\n";
		my @opt_flags_list = @toks[11..$#toks];
		my %opt_flags = ();
		next if $cigar eq "*"; # Failed to align
		# Get the read sequence & qualities from previous record if necessary
		if($seq eq "*") {
			$lastseq ne "" || die "Line #$idx:\n$line";
			$seq = $lastseq;
			$qual = $lastqual;
		} else {
			$lastseq = $seq;
			$lastqual = $qual;
		}
		$seq ne "*" || die;
		my ($parsed_cigar, $rdlen_cigar, $rflen_cigar) = parse_cigar($cigar);
		length($seq) == $rdlen_cigar ||
			die "Sequence length and parsed cigar string length ($rdlen_cigar) mismatch:\n".
			    "$seq\n$parsed_cigar\nLine:\n$line";
		# Stick optional flags into a hash
		for my $fl (@opt_flags_list) {
			my @fs = split(/:/, $fl, -1);
			scalar(@fs) > 2 || die "Bad optional flag: $fl\n$line\n";
			$opt_flags{$fs[0]}{type} = $fs[1];
			$opt_flags{$fs[0]}{value} = join(":", @fs[2..$#fs]);
		}
		defined($opt_flags{"MD"}) || die "No MD:Z flag:\n$line\n";
		$opt_flags{"MD"}{type} eq "Z" || die "Bad type for MD:Z flag\n$line\n";
		my $md = $opt_flags{"MD"}{value};
		$pos   == int($pos)   || die "POS field (col 4) must be an int:\n$line\n";
		$pnext == int($pnext) || die "PNEXT field (col 8) must be an int:\n$line\n";
		$tlen  == int($tlen)  || die "TLEN field (col 9) must be an int:\n$line\n";
		$mapq  == int($mapq)  || die "MAPQ field (col 5) must be an int:\n$line\n";
		# TODO: deal with bisulfite strands??
		my $fw = (($flag & 0x10) == 0);
		my $orient = $fw ? "+" : "-";
		my $reftype = $self->calcRefType($orient);
		defined($self->refs->{$reftype}{$rname}) ||
			die "No such refname as $rname for reftype $reftype:\n$line\n".
			Dumper($self->refs->{$reftype});
		my $exoff = $pos-1; # expected 0-based reference offset
		my ($alread, $alref, $htriml, $striml, $htrimr, $strimr) =
			_read_md_cigar_to_al($seq, $md, $cigar);
		print STDERR "$alread\n$alref\n";
		my $rfseq = $alref;
		$rfseq =~ s/[ -]//g; # remove spaces & gaps
		my $exlen = length($rfseq);
		my $refsub = substr($self->refs->{$reftype}{$rname}, $exoff, $exlen);
		length($refsub) == $exlen ||
			die "Tried to extract ref substring of length $exlen from:\n".
			    $self->refs->{$reftype}{$rname}.
				"\ngot string of length ".length($refsub).":\n".
				$refsub.
				"\nfrom:\n".
				$line.
				"\nexlen is the length of:\n$rfseq\npos=$pos, rfseq=$rfseq\n";
		$refsub = DNA::iupacSubN($refsub);
		my $trueRfseq = $refsub;
		length($trueRfseq) == length($rfseq) ||
			die "Different lengths for edited read and ref:\n".
			"       Read: $seq\n".
			"Edited read: $rfseq\n".
			"        Ref: $trueRfseq\n";
		$rfseq eq $trueRfseq ||
			die "Did not match:\n".
			"       Read: $seq\n".
			"Edited read: $rfseq\n".
			"        Ref: $trueRfseq\n".
			"$line";
		$self->{_nals}++;
	}
}

##
# Parse lines from a Bowtie alignment file and check that the
# alignments are sane and consistent with the reference.
#
sub parseBowtie {
	my ($self, $fh) = @_;
	while(<$fh>) {
		$self->parseBowtieLines([$_]);
	}
}

##
# Parse lines from a SAM alignment file and check that alignments are
# sane and consistent with the reference.
#
sub parseSam {
	my ($self, $fh) = @_;
	my @lines = ();
	while(<$fh>) { push @lines, $_; }
	$self->parseSamLines(\@lines);
}

##
# Parse lines from an alignment file of the type given by self->format
#
sub parseLines {
	my ($self, $lines) = @_;
	if($self->format eq "bowtie") {
		$self->parseBowtieLines($lines);
	} else {
		$self->format eq "sam" || die;
		$self->parseSamLines($lines);
	}
}

##
# Parse lines from an alignment file of the type given by self->format
#
sub parse {
	my ($self, $fh) = @_;
	if($self->format eq "bowtie") {
		$self->parseBowtie($fh);
	} else {
		$self->format eq "sam" || die;
		$self->parseSam($fh);
	}
}

##
# Print summary of how many alignments and edits were checked.
#
sub printSummary {
	my $self = shift;
	print STDERR "--- Summary ---\n";
	print STDERR "Read ".scalar(keys %{$self->refs})." reference strings\n";
	print STDERR "Checked $self->{_nals} alignments, $self->{_nedits} edits\n";
	print STDERR "---------------\n";
	print STDERR "PASSED\n";
}

##
# Check the given batch of alignments.  We check that they're
# internally consistent in some basic ways, and we check that the
# sequence and edits are consistent with the reference.
#
# The $als argument is either a list of (possibly compressed) filenames
# of files containing alignments, or a list of alignment strings.  If
# the former, $filenames is non-zero.
#
sub checkAlignments {
	my ($self, $als, $filenames) = @_;
	if($self->bisC) {
		print STDERR "Generating all-C bisulfite-treated references\n";
		bisulfiteC($self->refs);
	}
	if($self->bisCpG) {
		print STDERR "Generating all-CpG bisulfite-treated references\n";
		bisulfiteCpG($self->refs);
	}
	if($filenames) {
		foreach (@$als) {
			my $alnpipe = $_;
			print STDERR "Processing alignment file '$_'\n";
			$alnpipe = "gzip -dc $_ |" if ($_ =~ /\.gz$/);
			my $alnfh = undef;
			open($alnfh, $alnpipe) || die "Could not open '$alnpipe' for reading";
			$self->parse($alnfh);
			close($alnfh);
		}
	} else {
		$self->parseLines($als);
	}
}

##
# Check simple alignments
#
sub test1 {
	my $ac = AlignmentCheck->new(
		"AlignmentCheck.pm test1 checker",
		{ "r1" => "TTGTTCGT" },
		"bowtie",
		0,
		0
	);
	$ac->checkAlignments([
		"0\t+\tr1\t1\tTGTTCGT\tIIIIIII\t40\t-",
		"1\t+\tr1\t0\tTTGTTCG\tIIIIIII\t40\t-",
		"2\t+\tr1\t2\tGTTCGTA\tIIIIIII\t40\t6:N>A",
		"3\t+\tr1\t-1\tATTGTTC\tIIIIIII\t40\t0:N>A"], 0);
	return 1;
}

##
# Check simple alignments from files
#
sub test2 {
	open(TMP, ">/tmp/.AlignmentCheck.pm.fa") || die;
	print TMP ">r1\n";
	print TMP "TTGTTCGT\n";
	close(TMP);
	my $ac = AlignmentCheck->new(
		"AlignmentCheck.pm test1 checker",
		[ "/tmp/.AlignmentCheck.pm.fa" ],
		"bowtie",
		0,
		0
	);
	$ac->checkAlignments([
		"0\t+\tr1\t1\tTGTTCGT\tIIIIIII\t40\t-",
		"1\t+\tr1\t0\tTTGTTCG\tIIIIIII\t40\t-",
		"2\t+\tr1\t2\tGTTCGTA\tIIIIIII\t40\t6:N>A",
		"3\t+\tr1\t-1\tATTGTTC\tIIIIIII\t40\t0:N>A"], 0);
	return 1;
}

if($0 =~ /[^0-9a-zA-Z_]?AlignmentCheck\.pm$/) {
	my @fas = ();
	my @als = ();
	my $format = "sam";
	my $bisC = 0;
	my $bisCpG = 0;
	my $test = 0;
	
	use Getopt::Long;
	GetOptions (
		"test"                  => \$test,
		"fasta|ref=s"           => \@fas,
		"als|aln=s"             => \@als,
		"format=s"              => \$format,
		"bis-C|bisulfite-C"     => \$bisC,
		"bis-CpG|bisulfite-CpG" => \$bisCpG,
		"bowtie"                => sub {$format = "bowtie"},
		"sam"                   => sub {$format = "sam"}) || die;

	if($test) {
		use Test;
		print "Running unit tests\n";
		# Run unit tests
		Test::shouldSucceed("test1", \&test1);
		Test::shouldSucceed("test2", \&test2);
		exit 0;
	}

	my $ac = AlignmentCheck->new(
		"AlignmentCheck.pm checker",
		\@fas,
		$format,
		$bisC,
		$bisCpG);
	$ac->checkAlignments(\@als, 1);
}

1;
