#!/usr/bin/perl -w

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
# Parse a line from a SAM alignment file and check that the
# alignment is sane and consistent with the reference.
#
sub parseSamLines {
	my ($self, $lines) = @_;
	for my $line (@$lines) {
	}
}

##
# Parse lines from a SAM alignment file and check that alignments are
# sane and consistent with the reference.
#
sub parseSam {
	my ($self, $fh) = @_;
	while(<$fh>) { $self->parseSamLines([$_]); }
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
	my $format = "bowtie";
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
