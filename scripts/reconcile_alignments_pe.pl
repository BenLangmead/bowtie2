#!/usr/bin/perl -w

# 
# reconcile_alignments_pe.pl
#
#  Author: Ben Langmead
#    Date: 6/14/2009
#
# Reconcile and sanity-check an pair of input files, an output
# alignment file (the default, verbose kind), and a pair of output
# unaligned-read files (generated with the --un option) to be sure
# they're consistent with each other.  If the --max option was used
# during alignment, then the --max files should also be specified.  If
# the --al option was also used, specify the --al files as the last
# arguments to include it in the sanity check.
#
# Usage: perl reconcile_alignments_pe.pl \
#            [-k <int>] [-a] [-m <int>] \
#            [-u <int>] [-v] [-q | -f | -r] \
#            <input mate1 file> \
#            <input mate2 file> \
#            <alignment-file> \
#            <--un file 1> \
#            <--un file 2> \
#            <--max file 1> \
#            <--max file 2> \
#            <--al file 1> \
#            <--al file 2>
#

use strict;
use warnings;
use Getopt::Long;
Getopt::Long::Configure ("no_ignore_case");

my $khits = 1;
my $allHits = 0;
my $m = 999999;
my $M = undef;
my $u = -1;
my $verbose = 0;
my $fastq = 0;
my $fasta = 0;
my $raw = 0;
my $C = undef;
my $colCseq = undef;
my $colCqual = undef;
my $colCedit = undef;

GetOptions("m=i" => \$m,
           "M=i" => \$M,
           "k=i" => \$khits,
           "u=i" => \$u,
           "a"   => \$allHits,
           "v"   => \$verbose,
           "q"   => \$fastq,
           "f"   => \$fasta,
           "C"   => \$C,
           "col-cseq" => \$colCseq,
           "col-cedit" => \$colCedit,
           "col-cqual" => \$colCqual,
           "e"   => \$raw) || die "Bad arguments";

$m = $M if $M;

$fastq = 1 if ($fastq + $fasta + $raw == 0);

$khits = 999999 if $allHits;

print "khits: $khits\nmaxhits: $m\nnum_reads: $u\n" if $verbose;

# Utility function that returns the reverse complement of its argument
sub reverseComp($) {
	my $r = shift;
	$r = reverse($r);
	$r =~ tr/aAcCgGtT/tTgGcCaA/;
	return $r;
}

defined($ARGV[0]) || die "Must specify input read file 1 as first arg";
defined($ARGV[1]) || die "Must specify input read file 2 as first arg";
defined($ARGV[2]) || die "Must specify alignment file as second arg";
defined($ARGV[3]) || die "Must specify unaligned-read file 1 as third arg";
defined($ARGV[4]) || die "Must specify unaligned-read file 2 as third arg";

my $read1_file   = $ARGV[0];
my $read2_file   = $ARGV[1];
my $hits_file    = $ARGV[2];
my $un1_file = $ARGV[3];
my $un2_file = $ARGV[4];
my $max1_file = "";
$max1_file = $ARGV[5] if defined($ARGV[5]);
my $max2_file = "";
$max2_file = $ARGV[6] if defined($ARGV[6]);
my $al1_file = "";
$al1_file  = $ARGV[7] if defined($ARGV[7]);
my $al2_file = "";
$al2_file = $ARGV[8] if defined($ARGV[8]);

if($verbose) {
	print "read1: $read1_file\nread2: $read2_file\nalignments: $hits_file\n";
	print "unaligned reads 1: $un1_file\nunaligned reads 2: $un2_file\n";
	print "max reads 1: $max1_file\nmax reads 2: $max2_file\n";
}

sub get_read($) {
	my $fh = shift;
	my ($name, $seq, $qual) = ("", "", "");
	if($fastq) {
		$name = <$fh>;
		return ("", "", "") unless defined($name);
		chomp($name);
		$name = substr($name, 1);
		$seq = <$fh>;
		chomp($seq);
		my $tmp = <$fh>;
		$qual = <$fh>;
		chomp($qual);
	} elsif($fasta) {
		$name = <$fh>;
		return ("", "", "") unless defined($name);
		chomp($name);
		$name = substr($name, 1);
		$seq = <$fh>;
		chomp($seq);
	} else {
		$raw || die;
		die "Raw format not supported in reconcile_alignment_pe.pl; read names required";
	}
	return ($name, $seq, $qual);
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

my %hits_hash = ();
my %un_hash = ();
my %max_hash = ();
my %al_hash = ();

# Go through Bowtie-produced paired-end alignments
my $hits = 0;
my $distinctHits = 0;
open(HITS, $hits_file);
{
	my %num_hash = ();
	while(<HITS>) {
		# Read two alignments at a time, cleave off the trailing /1 or
		# /2 and assert that they have the same base name.
		my $l1 = $_;
		chomp($l1);
		my @s1 = split(/[\t]/, $l1);
		my $l2 = <HITS>;
		chomp($l2);
		my @s2 = split(/[\t]/, $l2);
		my $name1 = $s1[0];
		my $name2 = $s2[0];
		$hits++;
		my @name1s = split(/\//, $name1);
		my @name2s = split(/\//, $name2);
		$#name1s == $#name2s || die "Read names formatted differently: $name1, $name2";
		$name1s[$#name1s] eq "1" || $name1s[$#name1s] eq "2" || die "Bad suffix on read name: $name1";
		$name2s[$#name2s] eq "1" || $name2s[$#name2s] eq "2" || die "Bad suffix on read name: $name2";
		$name1s[$#name1s] ne $name2s[$#name2s] || die "Read names have same suffix: $name1, $name2";
		$name1s[0] eq $name2s[0] || die "Names don't match: $name1,$name2";
		my $swap = 0;
		$swap = 1 if $name1s[$#name1s] eq "2";
		my $read_name = $name1s[0];
		my $add_read = 1;
		# Check that the number of times that the basename has appeared
		# in a paired alignment does not violate the -k or -m policies.
		if($khits > 1) {
			# Update num_hash
			if(defined($num_hash{$read_name})) {
				$num_hash{$read_name}++;
				$add_read = 0; # already added this one
			} else {
				$num_hash{$read_name} = 1;
				$distinctHits++;
			}
			# Check whether num_hash violates -k
			if($num_hash{$read_name} > $khits) {
				die "Number of alignments for read $read_name exceeds -k limit of $khits";
			}
			# Check whether num_hash violates -m
			if($num_hash{$read_name} > $m) {
				die "Number of alignments for read $read_name exceeds -m limit of $m";
			}
		} else {
			defined($hits_hash{$read_name}) &&
				die "Read w/ name $read_name appears twice in alignment file";
			$distinctHits++;
		}
		if($add_read) {
			my $fw1 = ($s1[1] eq "+");
			my $seq1 = $s1[4];
			my $qual1 = $s1[5];
			if(!$fw1) {
				# Reverse / reverse-comp
				if($C && $colCseq) {
					$seq1 = reverse $seq1;
				} else {
					$seq1 = reverseComp($seq1);
				}
				$qual1 = reverse $qual1;
			}
			my $fw2 = ($s2[1] eq "+");
			my $seq2 = $s2[4];
			my $qual2 = $s2[5];
			if(!$fw2) {
				# Reverse / reverse-comp
				if($C && $colCseq) {
					$seq2 = reverse $seq2;
				} else {
					$seq2 = reverseComp($seq2);
				}
				$qual2 = reverse $qual2;
			}
			if(!$swap) {
				$hits_hash{$read_name}{seq1} = $seq1;
				$hits_hash{$read_name}{qual1} = $fasta ? "" : $qual1;
				$hits_hash{$read_name}{seq2} = $seq2;
				$hits_hash{$read_name}{qual2} = $fasta ? "" : $qual2;
			} else {
				$hits_hash{$read_name}{seq1} = $seq2;
				$hits_hash{$read_name}{qual1} = $fasta ? "" : $qual2;
				$hits_hash{$read_name}{seq2} = $seq1;
				$hits_hash{$read_name}{qual2} = $fasta ? "" : $qual1;
			}
		}
	}
	# Drop %num_hash - don't need it anymore
}
close(HITS);

# Go through entries of the FASTQ file for the unaligned reads
my $uns = 0;
my ($UN1, $UN2);
open $UN1, $un1_file;
open $UN2, $un2_file;
while(1) {
	my ($name1, $seq1, $qual1) = get_read($UN1);
	my ($name2, $seq2, $qual2) = get_read($UN2);
	if($name1 eq "") {
		$name2 eq "" || die;
		last;
	}
	my @ls1 = split(/\//, $name1);
	my @ls2 = split(/\//, $name2);
	$uns++;
	$#ls1 == $#ls2 || die "Differently-formatted names: $name1, $name2";
	$ls1[0] eq $ls2[0] || die "Different names for paired alignments: $ls1[0], $ls2[0]";
	my $name = $ls1[0];
	# Where have we seen it before?  Nowhere, hopefully
	if($M) {
		if(defined($hits_hash{$name})) {
			delete $hits_hash{$name};
			$distinctHits--;
		}
	} else {
		defined($hits_hash{$name}) &&
			die "Read $name appears both in hits file $hits_file and in --un file $un1_file/$un2_file";
	}
	defined($un_hash{$name}) &&
		die "Read $name appears in --un file $un1_file/$un2_file more than once";
	# Insert summary of the pair that didn't align
	$un_hash{$name}{seq1} = $seq1;
	$un_hash{$name}{qual1} = $qual1;
	$un_hash{$name}{seq2} = $seq2;
	$un_hash{$name}{qual2} = $qual2;
}
close($UN1);
close($UN2);

#
# Go through entries of the FASTQ file for the reads with a number of
# reportable alignments that exceeded the -m limit.
#
my $maxs = 0;
if($max1_file ne "") {
	$max2_file ne "" || die;
	my ($MAX1, $MAX2);
	my $bad = 0;
	open ($MAX1, $max1_file) || ($bad = 1);
	open ($MAX2, $max2_file) || ($bad = 1);
	while(!$bad) {
		my ($name1, $seq1, $qual1) = get_read($MAX1);
		my ($name2, $seq2, $qual2) = get_read($MAX2);
		if($name1 eq "") {
			$name2 eq "" || die;
			last;
		}
		$maxs++;
		my @ls1 = split(/\//, $name1);
		my @ls2 = split(/\//, $name2);
		$#ls1 == $#ls2 || die "Differently-formatted names: $name1, $name2";
		$ls1[0] eq $ls2[0] || die "Different names for paired alignments: $ls1[0], $ls2[0]";
		my $name = $ls1[0];
		# Where have we seen it before?  Nowhere, hopefully
		if($M) {
			defined($hits_hash{$name}) ||
				die "Read $name appears in --max file $max1_file/$max2_file but not in hits file $hits_file in -M mode";
			delete $hits_hash{$name};
			$distinctHits--;
		} else {
			defined($hits_hash{$name}) &&
				die "Read $name appears in hits file $hits_file and in --max file $max1_file/$max2_file";
		}
		defined($un_hash{$name}) &&
			die "Read $name appears in --un file $un1_file/$un2_file and in --max file $max1_file/$max2_file";
		defined($max_hash{$name}) &&
			die "Read $name appears in --max file $max1_file/$max2_file more than once";
		# Insert summary of the pair that didn't align
		$max_hash{$name}{seq1} = $seq1;
		$max_hash{$name}{qual1} = $qual1;
		$max_hash{$name}{seq2} = $seq2;
		$max_hash{$name}{qual2} = $qual2;
	}
	close($MAX1);
	close($MAX2);
}

#
# Go through entries of the 
#
my $als = 0;
if($al1_file ne "") {
	$al2_file ne "" || die;
	my ($AL1, $AL2);
	my $bad = 0;
	open ($AL1, $al1_file) || ($bad = 1);
	open ($AL2, $al2_file) || ($bad = 1);
	while(!$bad) {
		my ($name1, $seq1, $qual1) = get_read($AL1);
		my ($name2, $seq2, $qual2) = get_read($AL2);
		if($name1 eq "") {
			$name2 eq "" || die;
			last;
		}
		$als++;
		my @ls1 = split(/\//, $name1);
		my @ls2 = split(/\//, $name2);
		$#ls1 == $#ls2 || die "Differently-formatted names: $name1, $name2";
		$ls1[0] eq $ls2[0] || die "Different names for paired alignments: $ls1[0], $ls2[0]";
		my $name = $ls1[0];
		defined($hits_hash{$name}) ||
			die "Read $name appears in --al file $al1_file/$al2_file but not in hits file $hits_file";
		defined($un_hash{$name}) &&
			die "Read $name appears in --un file $un1_file/$un2_file and in --al file $al1_file/$al2_file";
		defined($max_hash{$name}) &&
			die "Read $name appears in --max file $max1_file/$max2_file and in --al file $al1_file/$al2_file";
		defined($al_hash{$name}) &&
			die "Read $name appears in --al file $al1_file/$al2_file more than once";
		# Insert summary of the pair that didn't align
		$al_hash{$name}{seq1} = $seq1;
		$al_hash{$name}{qual1} = $qual1;
		$al_hash{$name}{seq2} = $seq2;
		$al_hash{$name}{qual2} = $qual2;
	}
	close($AL1);
	close($AL2);
}

# Go through entries of the FASTQ file for the input reads and make
# sure that each entry is mirrored by an entry either in the alignment
# file or in the unaligned-read file.
my ($READ1, $READ2);
open $READ1, $read1_file;
open $READ2, $read2_file;
my $patid = 0;
my $reads = 0;
while(1) {
	my ($name1, $seq1, $qual1) = get_read($READ1);
	my ($name2, $seq2, $qual2) = get_read($READ2);
	if($name1 eq "") {
		$name2 eq "" || die;
		last;
	}
	$reads++;
	my @ls1 = split(/\//, $name1);
	my @ls2 = split(/\//, $name2);
	$#ls1 == $#ls2 || die "Differently-formatted names: $name1, $name2";
	$ls1[0] eq $ls2[0] || die "Different names for paired alignments: $ls1[0], $ls2[0]";
	my $name = $ls1[0];
	if(defined($hits_hash{$name})) {
		my $alseq1 = $seq1;
		my $alseq2 = $seq2;
		if($C && $colCseq) {
			$alseq1 = nucencode($seq1);
			$alseq2 = nucencode($seq2);
		}
		if(!$C || $colCseq) {
			$hits_hash{$name}{seq1} eq $alseq1 ||
				die "Read $name in alignment file has different _1 sequence than in ".
			 	    "input read file.\nAlgn: \"$hits_hash{$name}{seq1}\"\nInput: \"$alseq1\"";
		}
		# Qualities can be legitimately different
		#$hits_hash{$name}{qual1} eq $qual1 ||
		#	die "Read $name in alignment file has different _1 quals than in ".
		#	    "input read file.\nAlgn: $hits_hash{$name}{qual1}\nInput: $qual1";
		if(!$C || $colCseq) {
			$hits_hash{$name}{seq2} eq $alseq2 ||
				die "Read $name in alignment file has different _2 sequence than in ".
				    "input read file.\nAlgn: \"$hits_hash{$name}{seq2}\"\nInput: \"$alseq2\"";
		}
		#$hits_hash{$name}{qual2} eq $qual2 ||
		#	die "Read $name in alignment file has different _2 quals than in ".
		#	    "input read file.\nAlgn: $hits_hash{$name}{qual2}\nInput: $seq2";
	}
	elsif(defined($un_hash{$name})) {
		$un_hash{$name}{seq1} eq $seq1 ||
			die "Read $name in --un file has different _1 sequence than in ".
			    "input read file.\nAlgn: $un_hash{$name}{seq1}\nInput: $seq1";
		$un_hash{$name}{qual1} eq $qual1 ||
			die "Read $name in --un file has different _1 quals than in ".
			    "input read file.\nAlgn: $un_hash{$name}{qual1}\nInput: $qual1";
		$un_hash{$name}{seq2} eq $seq2 ||
			die "Read $name in --un file has different _2 sequence than in ".
			    "input read file.\nAlgn: $un_hash{$name}{seq2}\nInput: $seq2";
		$un_hash{$name}{qual2} eq $qual2 ||
			die "Read $name in --un file has different _2 quals than in ".
			    "input read file.\nAlgn: $un_hash{$name}{qual2}\nInput: $seq2";
	}
	elsif(defined($max_hash{$name})) {
		$max_hash{$name}{seq1} eq $seq1 ||
			die "Read $name in --max file has different _1 sequence than in ".
			    "input read file.\nAlgn: $max_hash{$name}{seq1}\nInput: $seq1";
		$max_hash{$name}{qual1} eq $qual1 ||
			die "Read $name in --max file has different _1 quals than in ".
			    "input read file.\nAlgn: $max_hash{$name}{qual1}\nInput: $qual1";
		$max_hash{$name}{seq2} eq $seq2 ||
			die "Read $name in --max file has different _2 sequence than in ".
			    "input read file.\nAlgn: $max_hash{$name}{seq2}\nInput: $seq2";
		$max_hash{$name}{qual2} eq $qual2 ||
			die "Read $name in --max file has different _2 quals than in ".
			    "input read file.\nAlgn: $max_hash{$name}{qual2}\nInput: $seq2";
	}
	else {
		die "Read $name does not appear in hits, --un or ---max file";
	}
	
	if(defined($al_hash{$name})) {
		$al_hash{$name}{seq1} eq $seq1 ||
			die "Read $name in --al file has different _1 sequence than in ".
			    "input read file.\nAlgn: $al_hash{$name}{seq1}\nInput: $seq1";
		$al_hash{$name}{qual1} eq $qual1 ||
			die "Read $name in --al file has different _1 quals than in ".
			    "input read file.\nAlgn: $al_hash{$name}{qual1}\nInput: $qual1";
		$al_hash{$name}{seq2} eq $seq2 ||
			die "Read $name in --al file has different _2 sequence than in ".
			    "input read file.\nAlgn: $al_hash{$name}{seq2}\nInput: $seq2";
		$al_hash{$name}{qual2} eq $qual2 ||
			die "Read $name in --al file has different _2 quals than in ".
			    "input read file.\nAlgn: $al_hash{$name}{qual2}\nInput: $seq2";
	}
	
	$patid++;
	last if $patid == $u;
}
close($READ1);
close($READ2);

if($al1_file ne "") {
	$als == $distinctHits || die "number of --al $als does not match number of distinct hits $distinctHits";
}
$distinctHits + $uns + $maxs == $reads ||
	die "distinct hits $distinctHits, --un $uns, --max $maxs doesn't add to reads $reads";

print "PASSED; processed $hits hits, $reads reads, $uns --un, $maxs --max, $als --al\n";
