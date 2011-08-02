#!/usr/bin/perl -w

# 
# reconcile_alignments.pl
#
#  Author: Ben Langmead
#    Date: 6/14/2009
#
# TODO: Support QSEQ reads
#
# Reconcile and sanity-check an input read file, a set of input
# reference files, and an output hits file (Bowtie verbose format), and
# an output unaligned-read file (generated with the --un option) to be
# sure they're consistent with each other.  If the --max option was
# used during alignment, then the --max file should also be specified.
# If the --al option was also used, specify the --al file as the last
# argument to include it in the sanity check.
#
# Usage: perl reconcile_alignments.pl \
#            [-k <int>] [-a] [-m <int>] \
#            [-u <int>] [-f|-q|-r] \
#            <input read file> \
#            <hits file> \
#            <--un file> \
#            <--max file> \
#            <--al file>
#

use strict;
use warnings;
use Getopt::Long;
Getopt::Long::Configure ("no_ignore_case");

my $k = undef;
my $m = undef;
my $M = undef;
my $u = undef;
my $f = undef;
my $r = undef;
my $a = undef;
my $q = undef;
my $C = undef;
my $t = undef;
my $qseq = undef;
my $colCseq = undef;
my $colCqual = undef;
my $colCedit = undef;
my $verbose = undef;

# For reference fasta files
my @reffn = ();
my %ref = ();
sub parseFasta {
	my $fa = shift;
	my $fapipe = ($fa =~ /\.gz$/ ? "gzip -dc $fa |" : $fa);
	open(FA, $fapipe) || die "Could not open '$fapipe'";
	my $curname = "";
	while(<FA>) {
		chomp;
		if(/^>/) {
			# Name line
			$curname = substr($_, 1);
			$curname =~ s/\s.*//;
			$ref{$curname} = "";
		} else {
			# Seq line
			$ref{$curname} .= $_ if $curname ne "";
		}
	}
	close(FA);
}

# Bowtie executables
my $bowtie = "";
my $bowtie_build = "";
my $index = "";

GetOptions (
	"bowtie=s"  => \$bowtie,
	"bowtie-build=s" => \$bowtie_build,
	"index=s"   => \$index,
	"k=i"       => \$k,
	"m=i"       => \$m,
	"M=i"       => \$M,
	"u=i"       => \$u,
	"q"         => \$q,
	"f"         => \$f,
	"qseq"      => \$qseq,
	"a"         => \$a,
	"C"         => \$C,
	"t"         => \$t,
	"fa|ref=s"  => \@reffn,
	"col-cseq"  => \$colCseq,
	"col-cedit" => \$colCedit,
	"col-cqual" => \$colCqual,
	"verbose"   => \$verbose,
	"r"         => \$r) || die "One or more errors parsing script arguments";

for my $f (@reffn) { parseFasta($f); }
my $checkFasta = scalar(@reffn) > 0;

##
# Check a nucleotide-space alignment to ensure it's consistent with the
# reference sequence aligned to.
#
sub checkNucAlignment {
	$checkFasta || die;
	my $al = shift;
	my ($rdname, $or, $refname, $refoff, $seq, $qual, $oms, $mms) = split(/\t/, $al, -1);
	my $origSeq = $seq;
	if($mms ne "-") {
		my @eds = split(/,/, $mms);
		if($or eq "+") {
			@eds = reverse @eds;
		}
		my $len = length($seq);
		# From RHS of read to LHS
		my $lastpos = 99999;
		for(my $ei = 0; $ei < scalar(@eds); $ei++) {
			my ($pos, $ed) = split(/:/, $eds[$ei]);
			$pos = $len-$pos-1 if $or eq "-";
			$pos <= $lastpos || die "$pos not <= $lastpos";
			$lastpos = $pos;
			my ($rfc, $rdc) = split(/>/, $ed);
			if($rfc eq "-") {
				# Reference gap; delete char out of read
				substr($seq, $pos, 1) = "";
			} elsif($rdc eq "-") {
				# Read gap; add characters to read
				substr($seq, $pos, 0) = $rfc;
			} else {
				# Mismatch; mutate read char
				substr($seq, $pos, 1) = $rfc;
			}
		}
	}
	my $rfstr = substr($ref{$refname}, $refoff, length($seq));
	$rfstr eq $seq ||
		die "Read sequence:\n".
		    "$origSeq\n".
		    "And edits:\n".
		    "$mms\n".
		    "Make this sequence:\n".
		    "$seq\n".
		    "But reference is:\n".
		    "$rfstr\n";
}

# Set khits
my $khits = 1;
$khits = int($k) if defined($k);
$khits = 999999 if $a;

# Set maxhits
my $maxhits = 999999;
$maxhits = int($m) if defined($m);
$maxhits = int($M) if defined($M);

# Set max # reads
my $num_reads = -1;
$num_reads = $u if defined($u);

# Set read format
my $format = "fastq";
$format = "fasta" if $f;
$format = "raw" if $r;

# Utility function that returns the reverse complement of its argument
sub reverseComp($) {
	my $r = shift;
	$r = reverse($r);
	$r =~ tr/aAcCgGtT/tTgGcCaA/;
	return $r;
}

##
# Encode colors as proxy nucleotides
#
sub nucencode($) {
	my $s = shift;
	my %nmap = ("0" => "A", "1" => "C", "2" => "G", "3" => "T", "." => "N",
	            "A" => "A", "C" => "C", "G" => "G", "T" => "T", "N" => "N");
	my $ret = "";
	for(my $i = 0; $i < length($s); $i++) {
		my $c = uc substr($s, $i, 1);
		defined($nmap{$c}) || die;
		$ret .= $nmap{$c};
	}
	return $ret;
}

defined($ARGV[0]) || die "Must specify input read file as first arg";
defined($ARGV[1]) || die "Must specify alignment file as second arg";

my $read_files  = $ARGV[0]; # reads
my $algn_file   = $ARGV[1]; # alignments
my $un_file     = "";
$un_file = $ARGV[2] if defined($ARGV[2]); # unaligned reads
my $max_file    = "";
$max_file = $ARGV[3] if defined($ARGV[3]); # non-uniquely aligned reads
my $al_file     = "";
$al_file = $ARGV[4] if defined($ARGV[4]);  # reads with >= 1 alignment

# hits_hash stores the read
my %hits_hash = ();
# al_hash stores the
my %al_hash = ();
my %un_hash = ();
my %max_hash = ();

##
# If the alignment file specified is -, this indicates that we should
# run Bowtie to get the alignments.
#
if($algn_file eq "-") {
	# Generate the alignment file by running Bowtie
	if(! -x $bowtie) {
		die "Alignment file was '-', but bowtie executable '$bowtie' failed the -x test";
	}
	# If --index wasn't specified, build index on the fly
	if($index eq "") {
		# Must be able to locate a bowtie-build
		if(! -x $bowtie_build) {
			die "Alignment file was '-' and --index was not specified, ".
			    "but bowtie-build executable '$bowtie_build' failed the ".
			    "-x test";
		}
		# Must have fasta files specified
		if(scalar(@reffn) == 0) {
			die "Alignment file was '-' and --index was not specified, ".
			    "but reference FASTA files were not specified with --ref";
		}
		my $refstr = join(",", @reffn);
		# Build the index
		my $cmd = "$bowtie_build $refstr .reconcile_alignments.index";
		print STDERR "$cmd\n";
		system($cmd);
		if($? != 0) {
			die "Bad exitlevel $? from bowtie-build command";
		}
		$index = ".reconcile_alignments.index";
	}
	# Index must exist
	if(! -f "$index.1.bt2") {
		die "Index basename $index.1.bt2 was not found";
	}
	my $btargs = "";
	$btargs .= " -k $k" if defined($k);
	$btargs .= " -m $m" if defined($m);
	$btargs .= " -M $M" if defined($M);
	$btargs .= " -u $u" if defined($u);
	$btargs .= " -q" if $q;
	$btargs .= " -f" if $f;
	$btargs .= " --qseq" if $qseq;
	$btargs .= " -a" if $a;
	$btargs .= " -r" if $r;
	$btargs .= " -C" if $C;
	$btargs .= " -t" if $t;
	$btargs .= " --col-cseq" if $colCseq;
	$btargs .= " --col-cedit" if $colCedit;
	$btargs .= " --col-cqual" if $colCqual;
	$btargs .= " --verbose" if $verbose;
	$btargs .= " $index"; # index
	$btargs .= " $read_files"; # reads
	$algn_file = ".reconcile_alignments.als";
	$un_file   = ".reconcile_alignments.un";
	$al_file   = ".reconcile_alignments.al";
	$max_file  = ".reconcile_alignments.max";
	my $cmd = "$bowtie $btargs $algn_file --un $un_file --max $max_file --al $al_file";
	print STDERR "$cmd\n";
	system($cmd);
	if($? != 0) {
		die "Bad exitlevel $? from bowtie command";
	}
}

# Go through Bowtie-produced alignments
my $hits = 0;
my $distinctHits = 0;
if(-f $algn_file) {
	open(ALGN, $algn_file) || die "Could not open alignment file $algn_file\n";
	my %num_hash = (); # for checking -k/-m
	while(<ALGN>) {
		chomp;
		checkNucAlignment($_) if $checkFasta;
		my @s = split /\t/;
		my $name = $s[0];
		my $add_read = 1;
		$hits++;
		if($khits > 1) {
			if(defined($num_hash{$name})) {
				$num_hash{$name}++;
				$add_read = 0; # already added this one
			} else {
				$num_hash{$name} = 1;
				$distinctHits++;
			}
			if($num_hash{$name} > $khits) {
				die "Number of alignments for read $name exceeds -k limit of $khits";
			} elsif($num_hash{$name} > $maxhits) {
				die "Number of alignments for read $name exceeds -m limit of $maxhits";
			}
		} else {
			defined($hits_hash{$name}{seq}) && 
				die "Read w/ name $name appears twice in alignment file but khits is 1";
			$distinctHits++;
		}
		if($add_read) {
			my $fw = ($s[1] eq "+");
			my $seq = $s[4];
			my $qual = $s[5];
			if(!$fw) {
				# Reverse / reverse-comp
				if($C && $colCseq) {
					$seq = reverse $seq;
				} else {
					$seq = reverseComp($seq);
				}
				$qual = reverse $qual;
			}
			$hits_hash{$name}{seq} = $seq;
			$hits_hash{$name}{qual} = (($format eq "fasta")? "" : $qual);
		}
	}
	close(ALGN);
}

##
# Parse a read from the given filehandle.  The format to parse is in
# $format.
#
sub get_read($) {
	my $fh = shift;
	my ($name, $seq, $qual) = ("", "", "");
	if($format eq "qseq") {
		die "qseq not supported";
		my $line = <$fh>;
		my @s = split(/\t/, $line, -1);
		my ($machname, $runnum, $lanenum,
		    $tilenum, $xcoord, $ycoord,
		    $index, $matenum, $seq1,
		    $qual1, $filter) = split(/\t/, $line, -1);
		$seq = $seq1;
		$qual = $qual1;
		$name = "${machname}_${runnum}_${lanenum}_${tilenum}_${xcoord}_${ycoord}_${index}/${matenum}";
	} elsif($format eq "fastq") {
		$name = <$fh>;
		return ("", "", "") unless defined($name);
		chomp($name);
		$name = substr($name, 1);
		$seq = <$fh>;
		chomp($seq);
		my $tmp = <$fh>;
		$qual = <$fh>;
		chomp($qual);
	} elsif($format eq "fasta") {
		$name = <$fh>;
		return ("", "", "") unless defined($name);
		chomp($name);
		$name = substr($name, 1);
		$seq = <$fh>;
		chomp($seq);
	} else {
		$format eq "raw" || die;
		die "Raw format not supported in reconcile_alignment.pl; read names required";
	}
	return ($name, $seq, $qual);
}

# Go through entries of the FASTQ file for the unaligned reads
my $uns = 0;
if($un_file ne "" && -f $un_file) {
	my $UN;
	open $UN, $un_file || die "Could not open unaligned-read file $un_file\n";
	while(1) {
		my ($name, $seq, $qual) = get_read($UN);
		last if $name eq "";
		$uns++;
		unless($M) {
			defined($hits_hash{$name}) &&
				die "Read $name appears both in hits file $algn_file and in --un file $un_file";
		} elsif(defined($hits_hash{$name})) {
			$distinctHits--;
			delete $hits_hash{$name};
		}
		defined($un_hash{$name}) &&
			die "Read $name appears more than once in --un file $un_file";
		$un_hash{$name}{seq} = $seq;
		$un_hash{$name}{qual} = $qual;
	}
	close($UN);
}

my $maxs = 0;
if($max_file ne "" && -f $max_file) {
	my $MAX;
	open $MAX, $max_file || die "Could not open maxed-read file $max_file\n";
	# Go through entries of the MAX file for the unaligned reads
	while(1) {
		my ($name, $seq, $qual) = get_read($MAX);
		last if $name eq "";
		$maxs++;
		if($M) {
			defined($hits_hash{$name}) ||
				die "Read $name appears in --max file $max_file but not in alignment file $algn_file";
			$distinctHits--;
			delete $hits_hash{$name};
		} else {
			defined($hits_hash{$name}) &&
				die "Read $name appears both in hits file $algn_file and in --max file $max_file";
		}
		defined($un_hash{$name})   &&
			die "Read $name appears both in --un file $un_file and in --max file $max_file";
		defined($max_hash{$name})  &&
			die "Read $name appears in --max file $max_file more than once";
		$max_hash{$name}{seq} = $seq;
		$max_hash{$name}{qual} = $qual;
	}
	close($MAX);
}

my $als = 0;
if($al_file ne "") {
	my $AL;
	open $AL, $al_file;
	# Go through entries of the MAX file for the unaligned reads
	while(1) {
		my ($name, $seq, $qual) = get_read($AL);
		last if $name eq "";
		$als++;
		defined($hits_hash{$name}) ||
			die "Read $name appears --al file $al_file but not in hits file $algn_file";
		defined($un_hash{$name})   &&
			die "Read $name appears both in --un file $un_file and in --al file $al_file";
		defined($max_hash{$name})  &&
			die "Read $name appears both in --max file $max_file and in --al file $al_file";
		defined($al_hash{$name})  &&
			die "Read $name appears in --al file $al_file more than once";
		$al_hash{$name}{seq} = $seq;
		$al_hash{$name}{qual} = $qual;
	}
	close($AL);
}

my @read_list = split(/,/, $read_files);
my $reads = 0;
for my $read_file (@read_list) {
	# Go through entries of the FASTQ file for the input reads and make
	# sure that each entry is mirrored by an entry either in the alignment
	# file or in the unaligned-read file.
	my $READ;
	open $READ, $read_file;
	my $patid = 0;
	while(1) {
		my ($name, $seq, $qual) = get_read($READ);
		last if $name eq "";
		$reads++;
		if(defined($hits_hash{$name})) {
			my $alseq = $seq;
			if($C && $colCseq) {
				$alseq = nucencode($seq);
			}
			if(!$C || $colCseq) {
				$hits_hash{$name}{seq} eq $alseq ||
					die "Read $name in hits file $algn_file has different sequence ".
					    "from input read.\nHit: \"$hits_hash{$name}{seq}\"\nInput: \"$alseq\"";
			}
			# Qualities can be legitimately different
			#$hits_hash{$name}{qual} eq $qual ||
			#	die "Read $name in hits file $algn_file has different sequence ".
			#	    "from input read.\nHit: $hits_hash{$name}{qual}\nInput: $qual";
		}
		elsif(defined($un_hash{$name})) {
			$un_hash{$name}{seq} eq $seq ||
				die "Read $name in --un file $un_file has different sequence ".
				    "from input read.\nHit: \"$un_hash{$name}{seq}\"\nInput: \"$seq\"";
			$un_hash{$name}{qual} eq $qual ||
				die "Read $name in --un file $un_file has different sequence ".
				    "from input read.\nHit: \"$un_hash{$name}{qual}\"\nInput: \"$qual\"";
		}
		elsif(defined($max_hash{$name})) {
			$max_hash{$name}{seq} eq $seq ||
				die "Read $name in --max file $max_file has different sequence ".
				    "from input read.\nHit: \"$max_hash{$name}{seq}\"\nInput: \"$seq\"";
			$max_hash{$name}{qual} eq $qual ||
				die "Read $name in --max file $max_file has different sequence ".
				    "from input read.\nHit: \"$max_hash{$name}{qual}\"\nInput: \"$qual\"";
		}
		elsif($un_file ne "") {
			die "Read with name $name appears in input, but not in any of the output files";
		}
		if(defined($al_hash{$name})) {
			$al_hash{$name}{seq} eq $seq ||
				die "Read $name in --al file $al_file has different sequence ".
				    "from input read.\nHit: \"$al_hash{$name}{seq}\"\nInput: \"$seq\"";
			$al_hash{$name}{qual} eq $qual ||
				die "Read $name in --al file $al_file has different sequence ".
				    "from input read.\nHit: \"$al_hash{$name}{qual}\"\nInput: \"$qual\"";
		}
		$patid++;
		last if $patid == $num_reads;
	}
	close($READ);
}

if($al_file ne "" && $un_file ne "") {
	$als == $distinctHits || die "number of --al $als does not match number of distinct hits $distinctHits";
}
if($un_file ne "") {
	$distinctHits + $uns + $maxs == $reads ||
		die "distinct hits $distinctHits, --un $uns, --max $maxs doesn't add to reads $reads";
}

# Print a summary of what was checked
print STDERR "Checked:\n";
print STDERR "  Reads: $reads\n";
print STDERR "  Alignments: $hits ($distinctHits from distinct reads)\n";
print STDERR "  Aligned reads: $als\n";
print STDERR "  Unaligned reads: $uns\n";
print STDERR "  Non-uniquely aligned reads: $maxs\n";
print STDERR "Read format: $format\n";
print STDERR "Checked that khits of $khits was not violated\n";
print STDERR "Checked that mhits of $maxhits was not violated\n";
if($un_file ne "") {
	print STDERR "Checked that all reads ended up in some output\n";
} else {
	print STDERR "DID NOT check that all reads ended up in some output; specify at least on an --un file for this\n";
}
if($checkFasta) {
	print STDERR "Checked all alignments against reference sequence\n";
} else {
	print STDERR "DID NOT check alignments against reference sequence; specify reference FASTAs with --ref\n";
}
print STDERR "PASSED; processed $hits hits, $reads reads, $uns --un, $maxs --max, $als --al\n";
