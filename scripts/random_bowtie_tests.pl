#!/usr/bin/perl -w

#
# Generate and run a series of random (but non-trivial) test cases for
# the Bowtie suite of tools, including bowtie, bowtie-build and
# bowtie-inspect.  Most of the problems turned up this way are in the
# form of assertions in the tools.  However, we also do some sanity-
# checking of the results; e.g. we use the pe_verify.pl script to
# verify that paired-end alignments are consistent with matched-up
# single-end alignments.
#
# Usage: perl random_tester.pl [rand seed] \
#                              [# outer iters] \
#                              [# inner iters] \
#                              [min # text bases] \
#                              [max text bases to add to min] \
#                              [min # read bases] \
#                              [max read bases to add to min]
#
# Options:
#   -n         don't attempt to compile binaries; use existing binaries
#   -p "<pol>" use only the specified alignment policy (e.g. "-n 3")
#

use List::Util qw[min max];
use Getopt::Std;

$| = 1; # flush immediately

my %options=();
getopts("mhnop:we",\%options);

if(defined $options{h}) {
	print "Usage: perl random_bowtie_tests.pl seed outer inner tbase trand pbase prand\n";
	exit 0;
}

# Let user specify a policy that will always override pickPolicy()
my $setPolicy;
$setPolicy = $options{p} if defined($options{p});

my $seed = 0;
$seed = int $ARGV[0] if defined($ARGV[0]);
srand $seed;

# make all the relevant binaries, unless we were asked not to
unless(defined $options{n}) {
	run("make bowtie bowtie-debug bowtie-build-debug ".
	       "bowtie-inspect-debug") == 0 || die "Error building";
}

# Alignment policies
my @policies = (
	"-n 3",
	"-n 2",
	"-n 1",
	"-n 0",
	"-v 3",
	"-v 2",
	"-v 1",
	"-v 0"
);

sub pickPolicy {
	my $pe = shift;
	my $r = int(rand($#policies + 1));
	my $pol = $policies[$r];
	$pol = $setPolicy if defined($setPolicy);
	if($pe && int(rand(2)) == 0) {
		my $min = int(rand(200))+10;
		$pol .= " -I $min";
		$pol .= " -X ".int(rand(30)+$min);
	}
	if($pol =~ /-n/ && int(rand(2)) == 0) {
		$pol .= " --nomaqround";
	}
	if($pol =~ /-n/ && int(rand(2)) == 0) {
		$pol .= " -l ".int(rand(30)+8);
	}
	if($pol =~ /-n/ && int(rand(2)) == 0) {
		$pol .= " -e ".int(rand(120)+40);
	}
	$pol .= " -g ";
	return $pol;
}

my $outer = 5000;
$outer = int $ARGV[1] if defined($ARGV[1]);
my $limit = $outer;

my $inner = 100;
$inner = int $ARGV[2] if defined($ARGV[2]);

my $tbase = 100;
$tbase = int $ARGV[3] if defined($ARGV[3]);
my $trand = 300;
$trand = int $ARGV[4] if defined($ARGV[4]);

my $pbase = 10;
$pbase = int $ARGV[5] if defined($ARGV[5]);
my $prand = 30;
$prand = int $ARGV[6] if defined($ARGV[6]);

my $ibase = 50;
$ibase = int $ARGV[7] if defined($ARGV[7]);
my $irand = 250;
$irand = int $ARGV[8] if defined($ARGV[8]);

my $verbose = 0;
my $exitOnFail = 1;
my @dnaMap = ('A', 'T', 'C', 'G',
              'N',
              'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B');

sub nonACGTtoN {
	my $t = shift;
	$t =~ tr/-NnMmRrWwSsYyKkVvHhDdBbXx/N/;
	$t =~ /[ACGTN]+/ || die "Bad N-ized DNA string: $t";
	return $t;
}

sub randGap() {
	my $or = int(rand(4));
	my $gap = "";
	if($or == 0) {
		my $ir = int(rand(100))+1;
		if(($ir & 1) == 1) {
			for(my $i = 0; $i < $ir; $i++) {
				$gap .= 'N';
			}
		} else {
			for(my $i = 0; $i < $ir; $i++) {
				$gap .= '-';
			}
		}
	}
	return $gap;
}

# Generates a random DNA string of the given length
sub randDna($) {
	my $num = shift;
	my $i;
	my $t = '';
	my $noAmbigs = int(rand(3)) == 0;
	for($i = 0; $i < $num; $i++) {
		my $or = int(rand(50));
		if($or == 0 && !$noAmbigs) {
			# Add a random, possibly ambiguous character
			$t .= $dnaMap[int(rand($#dnaMap+1))];
		} elsif($or == 1 && !$noAmbigs) {
			# Add a random-length streak of Ns (max: 20)
			my $streak = int(rand(20))+1;
			for(my $j = $i; $j < $num && $j < $streak; $j++) {
				$t .= 'N';
			}
		} else {
			# Add a random non-ambiguous character
			$t .= $dnaMap[int(rand(4))];
		}
	}
	return $t;
}

sub randColor($) {
	my %cmap = ("A" => "0", "C" => "1", "G" => "2", "T" => "3", "N" => ".");
	my $s = nonACGTtoN(randDna(shift));
	for(my $i = 0; $i < length($s); $i++) {
		defined($cmap{substr($s, $i, 1)}) || die "No entry for \"".substr($s, $i, 1)."\"\n";
		substr($s, $i, 1) = $cmap{substr($s, $i, 1)};
	}
	return $s;
}

# Utility function that returns the reverse complement of its argument
sub revcomp {
	my $r = shift;
	my $c = shift; # may be undef
	$r = reverse($r);
	$r =~ tr/aAcCgGtT/tTgGcCaA/ unless (defined($c) && $c);
	return $r;
}

# Utility function that returns the complement of its argument
sub comp {
	my $r = shift;
	my $c = shift;
	$r =~ tr/aAcCgGtT/tTgGcCaA/ unless (defined($c) && $c);
	return $r;
}

# Insert $src just before the $off'th character of $dst
sub ins($$$$$$) {
	my ($dst, $src, $off, $origoff, $origlen, $barrier) = @_;
	$origoff+1 > $barrier || die "Offset $off is within barrier ($barrier) from beginning";
	$origoff < $origlen   || die "Offset is not within read; length = $origlen";
	($origlen - $origoff) > ($barrier-1) || die "Offset $off is within barrier ($barrier-1) from end; length = $origlen";
	$origoff < length($dst) || die "ins offset $off is not less than dst length ".length($dst)."\n";
	substr($dst, $off, 0) = $src;
	return $dst;
}

# Delete the $off'th character of $dst
sub del($$$$$) {
	my ($dst, $off, $origoff, $origlen, $barrier) = @_;
	$origoff+1 > $barrier || die "Offset $off is within barrier ($barrier) from beginning";
	$origoff < $origlen   || die "Offset is not within read; length = $origlen";
	($origlen - $origoff) > $barrier || die "Offset $off is within barrier ($barrier) from end; length = $origlen";
	$origoff < length($dst) || die "del offset $off is not less than dst length ".length($dst)."\n";
	substr($dst, $off, 1) = "";
	return $dst;
}

# Substitute $src for the $off'th character of $dst
sub subst($$$) {
	my ($dst, $src, $off) = @_;
	substr($dst, $off, 1) = $src;
	return $dst;
}

##
# Sanity check whether a read and a list of edits corresponds correctly
# to a substring of the reference.
#
sub checkAlignmentRef {
	my ($ref, $read, $fw, $off, $edits, $barrier, $alnuc) = @_;
	my $orig = $read;
	my $origlen = length($read);
	my $origr = $fw ? $orig : revcomp($orig, $alnuc);
	$read = $origr;
	$off-- unless $alnuc;
	if($edits ne '-') {
		my $adjust = 0;
		my @es = split(/[,]/, $edits);
		for my $e (@es) {
			my $colonoff = index($e, ":");
			my $caratoff = index($e, ">");
			my $pos = substr($e, 0, $colonoff);
			my $ref = substr($e, $colonoff+1, $caratoff-$colonoff-1);
			my $qry = substr($e, $caratoff+1);
			length($qry) == 1 || die "Query char in edit $e, \"$qry\", isn't 1 char\n";
			$ref = comp($ref, $alnuc) unless $fw;
			if($qry eq "-") {
				# insertion
				my $inspos = $pos + $adjust;
				$read = ins($read, $ref, $inspos, $pos, $origlen, $barrier);
				$adjust += length($ref);
			} elsif($ref eq "-") {
				# deletion
				my $delpos = $pos + $adjust;
				$read = del($read, $delpos, $pos, $origlen, $barrier);
				$adjust--;
			} else {
				# mismatch
				my $mmpos = $pos + $adjust;
				$read = subst($read, $ref, $mmpos);
			}
		}
	}
	my $readr = $read;
	$read = revcomp($read, $alnuc) unless $fw;
	my $rstr = substr($ref, $off, length($read));
	$fw = $fw ? '+' : '-';
	$read eq $rstr || die "\n\nAlignment failed to pass sanity check:\n".
	                      "FW: $fw, Off: $off, Edits: $edits\n".
	                      "Orig:   $orig\n".
	                      "Orig5L: $origr\n".
	                      "Qry5L:  $readr\n".
	                      "Qry:    $read\n".
	                      "Ref:    $rstr\n";
}

##
# Given a string in nucleotide space, convert to colorspace.
#
sub colorize($$) {
	my ($s, $nucs) = @_;
	defined($s) || die;
	my %cmap = (
		"AA" => "0", "CC" => "0", "GG" => "0", "TT" => "0",
		"AC" => "1", "CA" => "1", "GT" => "1", "TG" => "1",
		"AG" => "2", "GA" => "2", "CT" => "2", "TC" => "2",
		"AT" => "3", "TA" => "3", "CG" => "3", "GC" => "3",
		"NA" => ".", "NC" => ".", "NG" => ".", "NT" => ".",
		"AN" => ".", "CN" => ".", "GN" => ".", "TN" => ".",
		"NN" => "."
	);
	my %nmap = ("0" => "A", "1" => "C", "2" => "G", "3" => "T", "." => "N");
	my $ret = "";
	for(my $i = 0; $i < length($s)-1; $i++) {
		my $di = uc substr($s, $i, 2);
		$di =~ tr/-NnMmRrWwSsYyKkVvHhDdBbXx/N/;
		defined($cmap{$di}) || die "Bad dinuc: $di\n";
		$ret .= ($nucs ? $nmap{$cmap{$di}} : $cmap{$di});
	}
	return $ret;
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

# Add some random quality values to encourage excercising the
# backtracking code
sub addQual($) {
	my $r = shift;
	my $len = length($r);
	$r .= ":";
	for(my $i = 0; $i < $len; $i++) {
		my $c = "-";
		while(not $c =~ /[0-9A-Z\/=@%]/) {
			$c = chr(33 + int(rand(41)));
		}
		$r .= $c;
	}
	return $r;
}

# Trim whitespace from a string argument
sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub run {
	my $cmd = shift;
	open(CMDTMP, ">.randtmp$seed.cmd") || die;
	print CMDTMP "$cmd\n";
	close(CMDTMP);
	return system($cmd);
}

sub runBacktick {
	my $cmd = shift;
	open(CMDTMP, ">.randtmp$seed.cmd") || die;
	print CMDTMP "$cmd\n";
	close(CMDTMP);
	return `$cmd`;
}

# Build an Ebwt based on given arguments
sub build {
	my($t, $color, $offRate, $ftabChars) = @_;
	my $ret = 0;
	
	my $file1 = "-c";
	my $file2 = "\"$t\"";
	if(substr($t, 0, 1) eq '-') {
		# Add backslash to escape first dash
		$file2 = "\"\\$t\"";
	}
	
	$color = ($color ? "-C" : "");
	
	# Write reference sequences to a FASTA file
	open(FA, ">.randtmp$seed.fa") || die "Could not open temporary fasta file";
	my @seqs = split(/,/, $t);
	for(my $i = 0; $i <= $#seqs; $i++) {
		print FA ">$i\n";
		print FA "$seqs[$i]\n";
	}
	close(FA);
	
	# Make a version of the FASTA file where all non-A/C/G/T characters
	# are Ns.  This is useful if we'd like to compare to the output of
	# bowtie-inspect.
	open(FAN, ">.randtmp$seed.ns.fa") || die "Could not open temporary fasta file";
	for(my $i = 0; $i <= $#seqs; $i++) {
		print FAN ">$i\n";
		my $t = nonACGTtoN($seqs[$i]);
		defined($t) || die;
		$t = colorize($t, 1) if $color;
		print FAN "$t\n";
	}
	close(FAN);
	
	my $fasta = int(rand(2)) == 0;
	if($fasta) {
		# Use the FASTA file as input
		$file1 = "-f";
		$file2 = ".randtmp$seed.fa";
	}
	
	my $bucketArg = "";
	my $bucketRand = int(rand(3));
	if($bucketRand == 0) {
		$bucketArg = "--bmaxdivn ";
		$bucketArg .= (int(rand(30))+1);
	} elsif($bucketRand == 1) {
		$bucketArg = "-a ";
	}

	$offRate   = "--offrate $offRate"   if $offRate ne "";
	$ftabChars = "--ftabchars $ftabChars" if $ftabChars ne "";
	my $noauto = "";
	$noauto = "-a" if $offRate eq "" && $ftabChars eq "";
	my $args = "-q --sanity $color $file1 $noauto $offRate $ftabChars $bucketArg $file2";
	
	# Do unpacked version
	my $cmd = "./bowtie-build-debug $args .tmp$seed";
	run("echo \"$cmd\" > .tmp$seed.build.cmd");
	print "$cmd\n";
	my $out = trim(runBacktick("$cmd 2>&1"));
	if($out eq "") {
		$ret++;
	} else {
		print "$out\n";
		if($exitOnFail) {
			exit 1;
		}
	}
	
	# Use bowtie-inspect to compare the output of bowtie-build to the
	# original reference sequences
	$cmd = "./bowtie-inspect-debug -a -1 .tmp$seed > .tmp$seed.inspect.ref";
	print "$cmd\n";
	run($cmd) == 0 || die "$cmd - failed";
	$cmd = "diff .randtmp$seed.ns.fa .tmp$seed.inspect.ref";
	print "$cmd\n";
	run($cmd) == 0 || die "$cmd - failed";

	# Do packed version and assert that it matches unpacked version
	# (sometimes, but not all the time because it takes a while)
	if(int(rand(4)) == 0) {
		$cmd = "./bowtie-build-debug -a -p $args .tmp$seed.packed";
		print "$cmd\n";
		$out = trim(runBacktick("$cmd 2>&1"));
		if($out eq "") {
			if(run("diff .tmp$seed.1.ebwt .tmp$seed.packed.1.bt2") != 0) {
				die if $exitOnFail;
			} elsif(run("diff .tmp$seed.2.ebwt .tmp$seed.packed.2.bt2") != 0) {
				die if $exitOnFail;
			} else {
				$ret++;
			}
		} else {
			print "$out\n";
			if($exitOnFail) {
				exit 1;
			}
		}
	}
	
	return $ret;
}

sub deleteReadParts {
	run("rm -f .tmp.un$seed". ".* .tmp.un$seed". "_1.* .tmp.un$seed". "_2.*");
	run("rm -f .tmp.max$seed".".* .tmp.max$seed"."_1.* .tmp.max$seed"."_2.*");
	run("rm -f .tmp.al$seed". ".* .tmp.al$seed". "_1.* .tmp.al$seed". "_2.*");
}

sub search($$$$$$$$$$) {
	my($tstr, $cstr, $pe, $color, $p1, $p2, $policy, $oneHit, $requireResult, $offRate) = @_;
	my $ret = doSearch($tstr, $cstr, $pe, $color, $p1, $p2, $policy, $oneHit, $requireResult, $offRate);
	deleteReadParts();
	return $ret;
}

##
# Make sure that a verbose alignment is consistent with the content of
# the reference
#
# $alnuc: Alignments are in nucleotides (not colors, as with -C --col-cseq)
#
sub checkRefVerbose($$$$$) {
	my ($l, $textsr, $ctextsr, $barrier, $alnuc) = @_;
	my @texts = @{$textsr};
	my @ctexts = @{$ctextsr};
	scalar(@texts) > 0 || die;
	my @ls = split(/\t/, $l);
	my ($fw, $ref, $off, $seq, $mms) = ($ls[1] eq '+', $ls[2], $ls[3], $ls[4], $ls[7]);
	defined($ref) || die "Malformed verbose output line: $l\n";
	$ref == int($ref) || die;
	$ref <= $#texts || die "Ref idx $ref exceeds number of texts ".($#texts+1)."\n";
	my $orig = $seq;
	print STDERR "Checking alignemnt: $l\n";
	checkAlignmentRef($alnuc ? $texts[$ref] : $ctexts[$ref],
	                  $orig, $fw, $off, $mms, $barrier, $alnuc);
	return 1;
}

# Search for a pattern in an existing Ebwt
sub doSearch($$$$$$$$$$) {
	my($nstr, $cstr, $pe, $color, $p1, $p2, $policy, $oneHit, $requireResult, $offRate) = @_;
	
	my @nts = split(/,/, $nstr); # nucleotide texts
	my @cts = split(/,/, $cstr); # color texts
	
	my $patarg = "-c";
	my $patstr = "\"$p1\"";
	$patstr = "-1 \"$p1\" -2 \"$p2\"" if $pe;
	my $outfile = ".tmp$seed.out";

	my $alnuc = 1;
	if($color) {
		$color = "-C";
		if(int(rand(3))) {
			$color .= " --col-cseq --col-cedit";
			$alnuc = 0;
		}
		$color .= " --col-cqual" if int(rand(3)) == 0;
	}
	
	my $ext = ".fq";
	my $format = int(rand(5));
	if($format == 0) {
		# FASTA
		open(FA, ">.randread$seed.fa") || die "Could not open temporary fasta file";
		$ext = ".fa";
		my @cs = split /[,]/, $p1;
		my $idx = 0;
		foreach my $c (@cs) {
			my @cms = split /[:]/, $c;
			print FA ">$idx\n$cms[0]\n";
			$idx++;
		}
		close FA;
		$patarg = "-f";
		$patstr = ".randread$seed.fa";
		if($pe) {
			run("mv .randread$seed.fa .randread$seed"."_1.fa");
			$patstr = "-1 .randread$seed"."_1.fa";
			open(FA, ">.randread$seed"."_2.fa") || die "Could not open temporary fasta file";
			@cs = split /[,]/, $p2;
			my $idx = 0;
			foreach my $c (@cs) {
				my @cms = split /[:]/, $c;
				print FA ">$idx\n$cms[0]\n";
				$idx++;
			}
			close FA;
			$patstr .= " -2 .randread$seed"."_2.fa";
		}
	} elsif($format == 1) {
		# FASTQ with ASCII qualities
		open(FQ, ">.randread$seed.fq") || die "Could not open temporary fastq file";
		$ext = ".fq";
		my @cs = split /[,]/, $p1;
		my $idx = 0;
		foreach my $c (@cs) {
			my @cms = split /[:]/, $c;
			$#cms == 1 || die "Should be pair with : separating seq from quals: $c\n$p1\n";
			print FQ "\@$idx\n$cms[0]\n+\n$cms[1]\n";
			$idx++;
		}
		close FQ;
		$patarg = "-q";
		$patstr = ".randread$seed.fq";
		if($pe) {
			run("mv .randread$seed.fq .randread$seed"."_1.fq");
			$patstr = "-1 .randread$seed"."_1.fq";
			open(FQ, ">.randread$seed"."_2.fq") || die "Could not open temporary fastq file";
			@cs = split /[,]/, $p2;
			my $idx = 0;
			foreach my $c (@cs) {
				my @cms = split /[:]/, $c;
				$#cms == 1 || die "Should be pair with : separating seq from quals: $c\n$p2\n";
				print FQ "\@$idx\n$cms[0]\n+\n$cms[1]\n";
				$idx++;
			}
			close FQ;
			$patstr .= " -2 .randread$seed"."_2.fq";
		}
	} elsif($format == 2) {
		# FASTQ with integer qualities
		open(FQ, ">.randread$seed.integer.fq") || die "Could not open temporary solexa fastq file";		my $ext = ".fa";
		$ext = ".fq";
		my @cs = split /[,]/, $p1;
		my $idx = 0;
		foreach my $c (@cs) {
			my @cms = split /[:]/, $c;
			print FQ "\@$idx\n$cms[0]\n+\n";
			for(my $i = 0; $i < length($cms[1]); $i++) {
				my $q = substr($cms[1], $i, 1);
				$q = ord($q) - 33;
				print FQ "$q ";
			}
			print FQ "\n";
			$idx++;
		}
		close FQ;
		$patarg = "-q --integer-quals";
		$patstr = ".randread$seed.integer.fq";
		if($pe) {
			run("mv .randread$seed.integer.fq .randread$seed.integer_1.fq");
			$patstr = "-1 .randread$seed.integer_1.fq";
			open(FQ, ">.randread$seed.integer_2.fq") || die "Could not open temporary fastq file";
			@cs = split /[,]/, $p2;
			my $idx = 0;
			foreach my $c (@cs) {
				my @cms = split /[:]/, $c;
				print FQ "\@$idx\n$cms[0]\n+\n";
				for(my $i = 0; $i < length($cms[1]); $i++) {
					my $q = substr($cms[1], $i, 1);
					$q = ord($q) - 33;
					print FQ "$q ";
				}
				$idx++;
				print FQ "\n";
			}
			close FQ;
			$patstr .= " -2 .randread$seed.integer_2.fq";
		}
	} elsif($format == 3) {
		# Raw
		open(RAW, ">.randread$seed.raw") || die "Could not open temporary raw file";
		$ext = ".raw";
		my @cs = split /[,]/, $p1;
		foreach my $c (@cs) {
			my @cms = split /[:]/, $c;
			print RAW "$cms[0]\n";
		}
		close RAW;
		$patarg = "-r";
		$patstr = ".randread$seed.raw";
		if($pe) {
			run("mv .randread$seed.raw .randread$seed"."_1.raw");
			$patstr = "-1 .randread$seed"."_1.raw";
			open(RAW, ">.randread$seed"."_2.raw") || die "Could not open temporary fastq file";
			@cs = split /[,]/, $p2;
			foreach my $c (@cs) {
				my @cms = split /[:]/, $c;
				print RAW "$cms[0]\n";
			}
			close RAW;
			$patstr .= " -2 .randread$seed"."_2.raw";
		}
	}
	
	# Perhaps dump unaligned reads using --un argument
	my $unalignArg = "";
	my $unalignReconArg = "";
	my $unalign = int(rand(5));
	if($unalign == 0 || $unalign == 2) {
		$unalignArg .= "--un .tmp.un$seed$ext ";
		if($pe) {
			$unalignReconArg .= " .tmp.un$seed"."_1$ext .tmp.un$seed"."_2$ext";
		} else {
			$unalignReconArg .= " .tmp.un$seed$ext";
		}
	}
	if($unalign == 1 || $unalign == 2) {
		$unalignArg .= "--max .tmp.max$seed$ext ";
		if($unalign == 2) {
			if($pe) {
				$unalignReconArg .= " .tmp.max$seed"."_1$ext .tmp.max$seed"."_2$ext";
			} else {
				$unalignReconArg .= " .tmp.max$seed$ext";
			}
		}
	}
	if($unalign == 2 || $unalign == 3) {
		$unalignArg .= "--al .tmp.al$seed$ext ";
		if($unalign == 2) {
			if($pe) {
				$unalignReconArg .= " .tmp.al$seed"."_1$ext .tmp.al$seed"."_2$ext";
			} else {
				$unalignReconArg .= " .tmp.al$seed$ext";
			}
		}
	}
	
	my $khits = "-k 1";
	my $mhits = 0;
	if(int(rand(2)) == 0) {
		$khits = "-a";
	} else {
		$khits = "-k " . (int(rand(20))+2);
	}
	if(int(rand(3)) == 0) {
		$khits .= " --strata";
	}
	if(int(rand(2)) == 0) {
		$requireResult = 0;
		$mhits = (int(rand(20))+2);
	}
	
	if($mhits > 0) {
		if(int(rand(2)) == 0) {
			$khits .= " -m $mhits";
		} else {
			$khits .= " -M $mhits";
		}
	}
	
	my $strand = "";
	if(int(rand(4)) == 0) {
		if(int(rand(2)) == 0) {
			# Reverse-complement reference strand only
			$strand = "--nofw";
		} else {
			# Forward reference strand only
			$strand = "--norc";
		}
	}
	
	if($oneHit || 1) {
		$oneHit = "";
	} else {
		$oneHit = "-a";
	}
	my $offRateStr = "";
	if($offRate ne "" && int(rand(3)) == 0) {
		$offRate == int($offRate) || die "Bad offrate: $offRate\n";
		$offRateStr = "--offrate " . ($offRate + 1 + int(rand(4)));
	}
	defined($policy) || die;
	defined($color) || die;
	defined($strand) || die;
	defined($unalignArg) || die;
	defined($khits) || die;
	defined($offRateStr) || die;
	defined($nstr) || die;
	my $cmd = "./bowtie-debug $policy $color $strand $unalignArg $khits $offRateStr --cost --orig \"$nstr\" $oneHit --sanity $patarg .tmp$seed $patstr";
	print "$cmd\n";
	my $out = trim(runBacktick("$cmd 2>.tmp$seed.stderr | tee .tmp$seed.stdout"));
	
	# Bad exitlevel?
	if($? != 0) {
		print "Exitlevel: $?\n";
		if($exitOnFail) {
			my $err = `cat .tmp$seed.stderr 2> /dev/null`;
			print "Stdout:\n$out\nStderr:\n$err\n";
			exit 1;
		}
		return 0;
	}
	my $err = `cat .tmp$seed.stderr 2> /dev/null`;

	# No output?
	if($out eq "" && $requireResult) {
		print "Expected results but got \"No Results\"\n";
		if($exitOnFail) {
			print "Stdout:\n$out\nStderr:\n$err\n";
			exit 1;
		}
		return 0;
	}
	
	# Parse output to see if any of it is bad
	print $out;
	my @outlines = split(/[\r\n]+/, $out);
	my %outhash = ();
	my %readcount = ();
	my %readStratum = (); # for unpaired reads
	my %readStratum1 = (); # for mate 1
	my %readStratum2 = (); # for mate 2
	my $lastread = "";
	my $barrier = 4;
	
	for(my $i = 0; $i <= $#outlines; $i++) {
		# Get the alignment for the first mate (or the unpaired read)
		my $l = $outlines[$i];
		chomp($l);
		checkRefVerbose($l, \@nts, \@cts, $barrier, $alnuc);
		my $key = "$l";
		my $l2 = "";
		if($pe) {
			# Get the alignment for the second mate
			$l2 = $outlines[++$i];
			defined($l2) || die "Odd number of output lines";
			chomp($l2);
			checkRefVerbose($l2, \@nts, \@cts, $barrier, $alnuc);
			$key .= ", $l2";
		}
		print "$key\n";
		# No two results should be the same
		if(!$pe) {
			!defined($outhash{$key}) || die "Result $key appears in output twice";
		}
		$outhash{$key} = 1;
		# Parse out the read id
		my @ls = split(/[\t]/, $l);
		my ($mate1, $fw1, $stratum1) = ($ls[0], $ls[1], $ls[8]);
		my ($mate2, $fw2, $stratum2) = (undef, undef, undef);
		if($pe) {
			my @l2s = split(/[\t]/, $l2);
			($mate2, $fw2, $stratum2) = ($l2s[0], $l2s[1], $l2s[8]);
		}
		my $peor = ($color ? "ff" : "fr");
		my $m1fw = ($peor eq 'fr' || $peor eq 'ff') ? '+' : '-';
		my $m2fw = ($peor eq 'ff') ? '+' : '-';
		if($pe && $peor eq 'ff') {
			$fw1 eq $fw2 || die "Saw different orientations for mates\n";
		} elsif($pe) {
			$fw1 ne $fw2 || die "Saw same orientation for mates\n";
		}
		if($strand =~ /nofw/) {
			if($pe) {
				my $m = substr($mate1, index($mate1, "/")+1);
				if($peor eq 'ff') {
					($fw1 eq '-' && $fw2 eq '-') ||
						die "Saw non-rc alignment with --nofw specified\n";
				} else {
					"$m$fw1" eq "2+" ||
						die "Saw a forward alignment on line ".($i+1)." when --nofw was specified";
				}
			} else {
				$fw1 eq "-" ||
					die "Saw a forward alignment on line ".($i+1)." when --nofw was specified";
			}
		} elsif($strand =~ /norc/) {
			if($pe) {
				if($peor eq 'ff') {
					($fw1 eq '+' && $fw2 eq '+') ||
						die "Saw non-fw alignment with --nofw specified\n";
				} else {
					my $m = substr($mate1, index($mate1, "/")+1);
					"$m$fw1" eq "1+" ||
						die "Saw a rev-comp alignment on line ".($i+1)." when --norc was specified";
				}
			} else {
				$fw1 eq "+" ||
					die "Saw a rev-comp alignment on line ".($i+1)." when --norc was specified";
			}
		}
		if($mate1 ne $lastread && !$pe) {
			die "Read $mate1 appears multiple times non-consecutively" if defined($readcount{$mate1});
		}
		if(!$pe && $policy =~ /--strata /) {
			!defined($readStratum{$mate1}) ||
				$readStratum{$mate1} == $stratum1 ||
				die "Incompatible strata: old: $readStratum{$mate1}, cur: $stratum1";
			$readStratum{$mate1} = $stratum1;
		} elsif($policy =~ /--strata /) {
			# Paired-end case; use minimum of 2 mates
			if($mate2 =~ /\/1$/) {
				# Swticheroo
				my $tmp = $mate1; $mate1 = $mate2; $mate2 = $tmp;
				$tmp = $stratum1; $stratum1 = $stratum2; $stratum2 = $tmp;
			}
			if(defined($readStratum1{$mate1})) {
				# One or ther other mate has to have the same stratum as we recorded previously
				defined($readStratum2{$mate2}) || die "$mate1 was seen before, but not $mate2";
				$readStratum1{$mate1} == $stratum1 ||
				$readStratum2{$mate2} == $stratum2 ||
				$readStratum1{$mate1} == $stratum2 ||
				$readStratum2{$mate2} == $stratum1 ||
					die "Incompatible strata: old: <$readStratum1{$mate1}, $readStratum2{$mate2}> cur: <$stratum1, $stratum2>";
			}
			$readStratum1{$mate1} = $stratum1;
			$readStratum2{$mate2} = $stratum2;
		}
		$lastread = $mate1;
		$readcount{$mate1}++;
		if($mhits > 0) {
			$readcount{$mate1} <= $mhits || die "Read $mate1 matched more than $mhits times";
		}
	}

	{
		$cmd = "./bowtie $policy $color $strand $unalignArg $khits $offRateStr --cost --orig \"$nstr\" $oneHit --sanity $patarg .tmp$seed $patstr";
		print "$cmd\n";
		my $out2 = trim(runBacktick("$cmd 2>.tmp$seed.stderr"));
		$out2 eq $out || die "Normal bowtie output did not match debug bowtie output";

		$cmd = "./bowtie $policy $color $strand $unalignArg $khits --orig \"$nstr\" $oneHit --sanity $patarg .tmp$seed $patstr $outfile";
		print "$cmd\n";
		my $out3 = trim(runBacktick("$cmd 2>.tmp$seed.stderr"));

		$cmd = "./bowtie --mm $policy $color $strand $unalignArg $khits --orig \"$nstr\" $oneHit --sanity $patarg .tmp$seed $patstr $outfile";
		print "$cmd\n";
		my $out4 = trim(runBacktick("$cmd 2>.tmp$seed.stderr"));
		$out3 eq $out4 || die "Normal bowtie output did not match memory-mapped bowtie output";
	}
	
	# Now do another run with verbose output so that we can check the
	# mismatch strings
	$cmd = "./bowtie-debug $policy $color $strand $unalignArg $khits $offRateStr --orig \"$nstr\" $oneHit --sanity $patarg .tmp$seed $patstr";
	print "$cmd\n";
	$out = trim(runBacktick("$cmd 2>.tmp$seed.stderr"));
	# Parse output to see if any of it is bad
	@outlines = split('\n', $out);
	for my $l (@outlines) {
		checkRefVerbose($l, \@nts, \@cts, $barrier, $alnuc);
	}
	
	if($pe) {
		# If search was for paired-end alignments, then verify the
		# outcome using pe_verify.pl
		my $pol = $policy;
		$pol =~ s/--.*//;
		my $patstr2 = $patstr;
		$patstr2 =~ s/-1//;
		$patstr2 =~ s/-2//;
		my $col = ($color ? "-C" : "");
		open(TMP, ">.tmp$seed.pe_verify.cmd") || die;
		$cmd = "perl scripts/pe_verify.pl --reference=.randtmp$seed.fa --args=\"--quiet\" $col -d $pol .tmp$seed $patstr2";
		print "$cmd\n";
		print TMP "$cmd\n";
		$out = trim(runBacktick("$cmd 2>.tmp$seed.pe_verify.stderr"));
		close(TMP);
		# Bad exitlevel?
		if($? != 0) {
			print "scripts/pe_verify.pl exitlevel: $?\n";
			if($exitOnFail) {
				my $err = `cat .tmp$seed.pe_verify.stderr 2> /dev/null`;
				print "Stdout:\n$out\nStderr:\n$err\n";
				exit 1;
			}
			return 0;
		}
	}
	
	if(!$pe) {
		my $col = ($color ? "-C" : "");
		$cmd = "perl scripts/best_verify.pl -d $policy $col .tmp$seed $patstr";
		print "$cmd\n";
		$out = trim(runBacktick("$cmd 2>.tmp$seed.best_verify.stderr"));
		
		if($? != 0) {
			print "scripts/best_verify.pl exitlevel: $?\n";
			if($exitOnFail) {
				my $err = `cat .tmp$seed.best_verify.stderr 2> /dev/null`;
				print "Stdout:\n$out\nStderr:\n$err\n";
				exit 1;
			}
			return 0;
		}
	}
		
	# Run the reconciler to ensure that --un, --max, and --al had
	# sane output
	# .tmp$seed.verbose.out
	if($format >= 0 && $format <= 2 && $unalignReconArg ne "") {
		deleteReadParts();
		my $cmd = "./bowtie $policy $color $strand $unalignArg $khits $offRateStr --orig \"$nstr\" $oneHit --sanity $patarg .tmp$seed $patstr .tmp$seed.verbose.out";
		print "$cmd\n";
		run($cmd) == 0 || die "Error performing reconciler run\n";
		$khits =~ s/--strata//;
		$patarg =~ s/--integer-quals//;
		$unalignReconArg =~ /\.tmp\.un/ || die;
		if($pe) {
			$patstr =~ s/-1//;
			$patstr =~ s/-2//;
			$cmd = "perl scripts/reconcile_alignments_pe.pl $color $patarg $khits $patstr .tmp$seed.verbose.out $unalignReconArg";
			print "$cmd\n";
			run($cmd) == 0 || die "Failed to reconcile";
		} else {
			$cmd = "perl scripts/reconcile_alignments.pl $color $patarg $khits $patstr .tmp$seed.verbose.out $unalignReconArg";
			print "$cmd\n";
			run($cmd) == 0 || die "Failed to reconcile";
		}
	}
	
	# Success
	return 1;
}

my $pass = 0;
my $tests = 0;
my $fail = 0;

for(; $outer > 0; $outer--) {

	# Generate random parameters
	my $offRate = "";
	if(int(rand(2)) == 0) {
		$offRate = int(rand(16));
	}
	my $ftabChars = "";
	if(int(rand(2)) == 0) {
		$ftabChars = 1 + int(rand(8));
	}

	# Generate random text(s)
	my $nt = int(rand(10)) + 1;
	my $tstr = '', $cstr = '';
	my @nts;
	my @cts;
	
	# For each reference sequence...
	for(my $i = 0; $i < $nt; $i++) {
		my $tlen = $tbase + int(rand($trand));
		my $tt = randDna($tlen);             # add text meat
		push(@nts, $tt);
		push(@cts, colorize($tt, 1));
		my $newt = randGap() . $tt . randGap();
		$tstr .= $newt;
		$cstr .= colorize($newt, 1);
		if($i < $nt-1) {
			$tstr .= ",";
			$cstr .= ",";
		}
	}
	
	my $color = (int(rand(2)) == 0);
	
	# Run the command to build the Ebwt from the random text
	$pass += build($tstr, $color, $offRate, $ftabChars);
	last if(++$tests > $limit);

	my $in = $inner;
	for(; $in >= 0; $in--) {
		# Paired-end?
		my $pe = (int(rand(2)) == 0);
		$pe = 1 if defined($options{e});
		# Generate random pattern(s) based on text
		my $pfinal1 = '';
		my $pfinal2 = '';
		# Decide number of patterns to generate
		my $np = int(rand(30)) + 1;
		for(my $i = 0; $i < $np; $i++) {
			# Pick a text
			my $tt = $nts[int(rand($#nts))];
			# Pick a length
			my $pl;
			my $plen;
			my $pr;
			my $p1;
			my $p2;
			if($pe) {
				# Pick a left-hand offset for the insert
				my $il = int(rand(length($tt))) - 10;
				# Pick a length for the insert
				my $ilen = int(rand($irand)) + $ibase;
				$il = min($il, length($tt)-$ilen);
				my $ir = min($il + $ilen, length($tt));
				my $is = substr $tt, $il, $ir - $il;
				# Pick a length for the #1 mate
				my $plen1 = int(rand($prand)) + $pbase;
				# Pick a length for the #2 mate
				my $plen2 = int(rand($prand)) + $pbase;
				$p1 = substr $is, 0, $plen1;
				$is = reverse $is;
				$p2 = substr $is, 0, $plen2;
				$p2 = reverse $p2;
				$p2 = revcomp($p2);
			} else {
				# Pick a length for the read
				$plen = int(rand($prand)) + $pbase;
				$pl = int(rand(length($tt))) - 10;
				$pl = max($pl, $color ? 5 : 4);
				$pl = min($pl, length($tt));
				$pr = min($pl + $plen, length($tt));
				$p1 = substr $tt, $pl, $pr - $pl;
			}
			# Check for empty pattern or pattern that spans a comma
			if(length($p1) < ($color ? 5 : 4) || index($p1, ",") != -1) {
				$i--; next;
			}
			if($pe && (length($p2) < ($color ? 5 : 4) || index($p2, ",") != -1)) {
				$i--; next;
			}
			# Optionally add nucleotide changes to pattern
			if($i > 0) {
				my $nummms = int(rand($color ? 3 : 5));
				for(my $j = 0; $j < $nummms; $j++) {
					substr($p1, int(rand(length($p1))), 1) = randDna(1);
				}
				if($pe) {
					$nummms = int(rand(4));
					for(my $j = 0; $j < $nummms; $j++) {
						substr($p2, int(rand(length($p2))), 1) = randDna(1);
					}
				}
			}
			# Possibly reverse complement it
			if((int(rand(2)) == 0)) {
				$p1 = revcomp($p1);
				if($pe) {
					$p2 = revcomp($p2);
					my $ptmp = $p1;
					$p1 = $p2;
					$p2 = $ptmp;
				}
			}
			$p1 =~ tr/MRWSYKVHDBX/N/;
			if($color) {
				defined($p1) || die;
				!$pe || defined($p2) || die;
				my ($nize1, $nize2) = (int(rand(2)) == 0, int(rand(2)) == 0);
				$p1 = colorize($p1, $nize1);
				$p2 = colorize($p2, $nize2) if $pe;
				# Optionally add color changes to pattern
				if($i > 0) {
					my $nummms = int(rand(3));
					for(my $j = 0; $j < $nummms; $j++) {
						my $r = randColor(1);
						$r = nucencode($r) if $nize1;
						substr($p1, int(rand(length($p1))), 1) = $r;
					}
					if($pe) {
						$nummms = int(rand(4));
						for(my $j = 0; $j < $nummms; $j++) {
							my $r = randColor(1);
							$r = nucencode($r) if $nize2;
							substr($p2, int(rand(length($p2))), 1) = $r;
						}
					}
				}
			}
			# Add valid random quality values
			$p1 = addQual($p1);
			$pfinal1 .= $p1;
			$pfinal1 .= "," if($i < $np-1);
			if($pe) {
				$p2 =~ tr/MRWSYKVHDBX/N/;
				$p2 = addQual($p2);
				$pfinal2 .= $p2;
				$pfinal2 .= "," if($i < $np-1);
			}
		}
		
		# Run the command to search for the pattern from the Ebwt
		my $oneHit = (int(rand(3)) == 0);
		my $policy = pickPolicy($pe);
		my $expectResult = 1;
		if(!$pe) {
			for(my $i = 0; $i < length($pfinal1); $i++) {
				my $c = substr($pfinal1, $i, 1);
				last if($c eq ',');
				if($c ne 'A' && $c ne 'C' && $c ne 'G' && $c ne 'T') {
					$expectResult = 0;
					last;
				}
			}
		} else {
			$expectResult = 0;
		}
		$pass += search($tstr, $cstr, $pe, $color, $pfinal1, $pfinal2, $policy, $oneHit, $expectResult, $offRate); # require 1 or more results
		last if(++$tests > $limit);
	}

	$in = $inner;
	for(; $in >= 0; $in--) {
		my $pe = (int(rand(2)) == 0);
		# Generate random pattern *not* based on text
		my $pfinal1 = '';
		my $pfinal2 = '';
		my $np = int(rand(10)) + 1;
		for(my $i = 0; $i < $np; $i++) {
			my $plen = int(rand($prand)) + $pbase;
			my $p1 = randDna($plen);
			$plen = int(rand($prand)) + $pbase;
			my $p2 = randDna($plen);
			$p1 =~ tr/MRWSYKVHDBX/N/;
			$p1 = colorize($p1, int(rand(2)) == 0) if $color;
			$p1 = addQual($p1);
			$p2 =~ tr/MRWSYKVHDBX/N/ if $pe;
			$p2 = colorize($p2, int(rand(2)) == 0) if $color && $pe;
			$p2 = addQual($p2) if $pe;
			$pfinal1 .= $p1;
			$pfinal2 .= $p2 if $pe;
			if($i < $np-1) {
				$pfinal1 .= ",";
				$pfinal2 .= "," if $pe;
			}
		}

		# Run the command to search for the pattern from the Ebwt
		my $oneHit = (int(rand(3)) == 0);
		my $policy = pickPolicy($pe);
		$pass += search($tstr, $cstr, $pe, $color, $pfinal1, $pfinal2, $policy, $oneHit, 0, $offRate); # do not require any results
		last if(++$tests > $limit);
	}
}

print "$pass tests passed, $fail failed\n";
exit 1 if $fail > 0;
