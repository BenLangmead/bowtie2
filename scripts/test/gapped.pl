#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my $bwa = "";
my $novo = "";
my $upto = 0;

GetOptions (
	"bwa:s" => \$bwa,
	"novo:s" => \$novo,
	"upto:i" => \$upto
) || die "Bad options\n";

my $btbuild = "./bowtie-build";
my $btalign = "./bowtie";
my $btalign_d = "./bowtie-debug";
system("make $btbuild $btalign $btalign_d") == 0 || die;

# Utility function that returns the reverse complement of its argument
sub revcomp($) {
	my $r = shift;
	$r = reverse($r);
	$r =~ tr/aAcCgGtT/tTgGcCaA/;
	return $r;
}

# Utility function that returns the complement of its argument
sub comp($) {
	my $r = shift;
	$r =~ tr/aAcCgGtT/tTgGcCaA/;
	return $r;
}

# Insert $src just before the $off'th character of $dst
sub ins {
	my ($dst, $src, $off) = @_;
	$off < length($dst) || die "ins offset $off is not less than dst length ".length($dst)."\n";
	my $s = substr($dst, 0, $off) . $src . substr($dst, $off);
	return $s;
}

# Delete the $off'th character of $dst
sub del {
	my ($dst, $off) = @_;
	my $s = substr($dst, 0, $off).substr($dst, $off+1);
	return $s;
}

# Substitute $src for the $off'th character of $dst
sub subst {
	my ($dst, $src, $off) = @_;
	my $s = substr($dst, 0, $off).$src.substr($dst, $off+1);
	return $s;
}

##
# Sanity check whether a read and a list of edits corresponds correctly
# to a substring of the reference.
#
sub checkAlignmentRef {
	my ($ref, $read, $decread, $fw, $off, $edits, $color) = @_;
	$read = $decread if $color;
	my $orig = $read;
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
			$ref = comp($ref) unless $fw;
			if($qry eq "-") {
				# insertion
				#my $inspos = ($fw ? $pos + $adjust : length($read) - $pos - $adjust);
				my $inspos = $pos + $adjust;
				$read = ins($read, $ref, $inspos);
				$adjust += length($ref);
			} elsif($ref eq "-") {
				# deletion
				#my $delpos = ($fw ? $pos + $adjust : length($read) - $pos - $adjust);
				my $delpos = $pos + $adjust;
				$read = del($read, $delpos);
				$adjust--;
			} else {
				# mismatch
				#my $mmpos = ($fw ? $pos + $adjust : length($read) - $pos - $adjust);
				my $mmpos = $pos + $adjust;
				$read = subst($read, $ref, $mmpos);
			}
		}
	}
	$read = revcomp($read) unless $fw;
	my $rstr = substr($ref, $off, length($read));
	$read eq $rstr || die "FW: $fw, Off: $off, Edits: $edits\nOrig read:  $orig\nXform read: $read\nRef:        $rstr\n";
}

my @unp_cases = (

	#
	# Small inserts
	#

	# Small insert too close to the beginning, but seeded strategy is
	# required to align it (regression test to make sure --gbar is
	# working)
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCTCGATCGTATCTCTACGCATGGCG:1:-n 1 -l 20 -a:0:0:0",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCAATCGATCGTATCTCTACGCATGGCG:1:-n 1 -l 20 -a:0:0:0",
	
	# Regressions 2/9/10:
	"GCGAGCTGATGGAGAATGACTAAAGGTTGTGCACGAGTGGTACCCGCACACT:GGTGCGTCAGGGATCGCGTAAGACACCCTCCTG:0:-C -n 1 -a -e 100:1:0:0",
	"GCGAGCTGATGGAGAATGACTAAAGGTTGTGCACGAGTGGTACCCGCACACT:GGTGCGTCAGGGATCGCGTAAGACACCCTCCTG:0:-C -n 2 -a -e 100:1:0:0",
	"GCGAGCTGATGGAGAATGACTAAAGGTTGTGCACGAGTGGTACCCGCACACT:GGTGCGTCAGGGATCGCGTAAGACACCCTCCTG:0:-C -n 3 -a -e 100:1:0:0",

	# Regression 6/3/10:
	"ACCCAAGGATAGACACACCGGGCACATGCCACCTTCGAGATCCTTATCTATCAACCACCTTAAGACGGTAAACGCACTTGATCTGGCAGG:CAACAGAGTTGGC:0:-v 3 -C -a -m 15:2:0:0",
	"ACCCAAGGATAGACACACCGGGCACATGCCACCTTCGAGATCCTTATCTATCAACCACCTTAAGACGGTAAACGCACTTGATCTGGCAGG:CAACAGAGTTGGC:0:-v 3 --gbar 5 -C -a -m 15:1:0:0",

	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -l 9 -e 210 -O 40 -E 30:1",

	# Trying to fool the redundant path detection by putting a
	# homopolymer spanning the halfway point
	"CCCCCCCCCCCCCCCAAAAAAAAAACCCCCCCCCCCCCCCCC:CCCCCCCCCCAAAAAAAAACCCCCCCCCC:1:-v 1 -a -O 40:2:1:30",
	"CCCCCCCCCCCCCCCAAAAAAAAAACCCCCCCCCCCCCCCCC:CCCCCCCCCCAAAAAAAAACCCCCCCCCC:1:-v 1 -a -O 40 --redundant 2:4:1:30",
	"CCCCCCCCCCCCCCCAAAAAAAAAACCCCCCCCCCCCCCCCC:CCCCCCCCCCAAAAAAAAACCCCCCCCCC:1:-v 1 -a -O 40 --redundant 3:12:1:30",

	# Read is missing the 10th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCTATCTCTACGCATGGCG:1:-v 0 -a:0:0:0",
	# Read is missing the 10th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCTATCTCTACGCATGGCG:1:-v 1 -a --gbar 11:0:0:0",

	# Read is missing the 5th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCACGATCGTATCTCTACGCATGGCG:1:-v 1 -a:1:1:100",
	# Read is missing the 5th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCACGATCGTATCTCTACGCATGGCG:1:-v 1 -a --gbar 6:0:0:0",
	# Read is missing the 5th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCACGATCGTATCTCTACGCATGGCG:1:-v 1 -a --gbar 4:1:1:100",

	# Read is missing the 4th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCTCGATCGTATCTCTACGCATGGCG:1:-v 1 -a:0:0:0",
	# Read is missing the 4th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCTCGATCGTATCTCTACGCATGGCG:1:-v 1 -a --gbar 6:0:0:0",
	# Read is missing the 4th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCTCGATCGTATCTCTACGCATGGCG:1:-v 1 -a --gbar 5:0:0:0",
	# Read is missing the 4th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCTCGATCGTATCTCTACGCATGGCG:1:-v 1 -a --gbar 4:0:0:0",
	# Read is missing the 4th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCTCGATCGTATCTCTACGCATGGCG:1:-v 1 -a --gbar 3:1:1:100",

	# Read is missing the 2nd char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:ACATCGATCGTATCTCTACGCATGGCG:1:-v 1 -a:1:1:30",
	# Read is missing the 10th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCTATCTCTACGCATGGCG:1:-v 1 -a:1:1:100",
	# Read is missing the 10th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCTATCTCTACGCATGGCG:1:-v 2:1:1:100",
	# Read is missing the 10th and 11th chars from the reference, so that's a size-2 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCATCTCTACGCATGGCG:1:-v 1 -a:0:0:0",
	# Read is missing the 10th and 11th chars from the reference, so that's a size-2 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCATCTCTACGCATGGCG:1:-v 2 -a:1:2:140",
	# Read is missing 10th, 11th 12th chars
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCTCTCTACGCATGGCG:1:-v 2 -a:0:0:0",
	# Read is missing 10th, 11th 12th chars
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCTCTCTACGCATGGCG:1:-v 3 -a --redundants 3:2:3:180",
	# Read is missing 10th, 11th 12th chars
	"AGCATCGATCAAATCTCTACGCATGGCG:AGCATCGATCTCTCTACGCATGGCG:1:-v 3 -a:1:3:180",
	# Read is missing the 10th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCTATCTCTACGCATGGCG:0:-n 0:0:0:0",
	# Read is missing the 10th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCTATCTCTACGCATGGCG:0:-n 1:0:0:0",
	# Read is missing the 10th char from the reference, so that's a size-1 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCTATCTCTACGCATGGCG:0:-n 1 -O 40 -E 40:1:1:40",
	# Read is missing the 10th and 11th chars from the reference, so that's a size-2 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCATCTCTACGCATGGCG:0:-n 1 -O 40 -E 40:0:0:0",
	# Read is missing the 10th and 11th chars from the reference, so that's a size-2 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCATCTCTACGCATGGCG:0:-n 2 -O 40 -E 40:0:0:0",
	# Read is missing the 10th and 11th chars from the reference, so that's a size-2 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCATCTCTACGCATGGCG:0:-n 2 -O 40 -E 30:1:2:70",
	# Read is missing the 10th and 11th chars from the reference, so that's a size-2 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCATCTCTACGCATGGCG:0:-n 2 -O 30 -E 40:1:2:60",
	# Read is missing the 10th and 11th chars from the reference, so that's a size-2 insert
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCATCTCTACGCATGGCG:0:-n 2 -O 20 -E 10:1:2:30",

	#
	# Small deletions
	#

	# Inserted "ACA" after 10 chars
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:1:-v 0:0:0:0",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:1:-v 1:0:0:0",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:1:-v 2:0:0:0",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:1:-v 3:1:3:180",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:1:-v 3 --gbar 10:1:3:180",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:1:-v 3 --gbar 11:0:0:0",

	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -e 90 -O 40 -E 30:0:0:0",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -e 99 -O 40 -E 30:0:0:0",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -e 100 -O 40 -E 30:1:3:100",
	#"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -e 110 -O 40 -E 30:1",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -l 9 -e 90 -O 40 -E 30:0",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -l 9 -e 100 -O 40 -E 30:1",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -l 10 -e 90 -O 40 -E 30:0",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -l 10 -e 100 -O 40 -E 30:1",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -l 11 -e 90 -O 40 -E 30:0",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -l 11 -e 100 -O 40 -E 30:1",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -l 12 -e 90 -O 40 -E 30:0",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -l 12 -e 100 -O 40 -E 30:1",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -l 13 -e 90 -O 40 -E 30:0",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 3 -l 13 -e 100 -O 40 -E 30:1",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 0 -l 9 -e 90 -O 40 -E 30:0",
	"AGCATCGATCGTATCTCTACGCATGGCG:AGCATCGATCACAGTATCTCTACGCATGGCG:0:-n 0 -l 9 -e 100 -O 40 -E 30:1",
	
"" );

my $lastcolor = 0;
my $lastref = "";
my $idx = 0;
for my $case (@unp_cases) {
	$idx++;
	last if $upto != 0 && $idx > $upto;
	next if $case eq "";
	my ($ref, $read, $revcomp, $args, $alns, $strat, $cost) = split(/:/, $case);
	my $color = ($args =~ /-C/ ? "-C" : "");
	if($ref ne $lastref || $color ne $lastcolor) {
		my $cmd = "$btbuild -c $color $ref .gapped.pl.tmp >/dev/null";
		print "$cmd\n";
		system($cmd) == 0 || die "";
	}
	for(my $i = 0; $i <= $revcomp; $i++) {
		# Try with --redundant 0 (or some other --redundant, if it's specified)
		$read = revcomp($read) if $i == 1;
		my $it = 0;
		my @exes = ($btalign_d, $btalign_d, $btalign);
		for (@exes) {
			my $exe = $_;
			$it++;
			my $myargs = $args;
			if($it == 2) {
				if($myargs =~ /--redundant/) {
					print "Skipping --redundant modification because it's already specified\n";
					next;
				} else {
					$myargs .= " --redundant 1";
				}
			}
			my $cmd = "$exe --sanity --orig $ref --cost -g -c $myargs .gapped.pl.tmp $read";
			print "$cmd\n";
			open BTIE, "$cmd 2>/dev/null |" || die "Could not open pipe for '$cmd'\n";
			my $als = 0;
			my @lines;
			while(<BTIE>) {
				push @lines, $_;
			}
			close(BTIE);
			$? == 0 || die "Command quit with exitlevel $?\n";
			for (@lines) {
				# Check that alignment is valid
				print $_;
				chomp;
				my @s = split(/[\t]/);
				if($als == 0) {
					# Parse stratum and cost
					my $mystrat = $s[-2];
					my $mycost = $s[-1];
					!defined($strat) || $mystrat == $strat || die "Expected stratum $strat, got $mystrat\n$case\n";
					!defined($cost)  || $mycost  == $cost  || die "Expected cost $cost, got $mycost\n$case\n";
				}
				checkAlignmentRef($ref, $read, $s[4], $s[1] eq '+', $s[3], $s[7], $color eq "-C");
				$als++;
			}
			$als == $alns || die "Expected $alns alignments, got $als\n$case\n";
			print "PASSED case: $case\n";
		}
		if($bwa ne "") {
			-x $bwa || die;
			my $bwacolor = $color eq "-C" ? "-c" : "";
			open(BWAFA, ">.bwa.fa") || die;
			print BWAFA ">0\n$ref\n";
			close(BWAFA);
			open(BWAFQ, ">.bwa.fq") || die;
			print BWAFQ "\@0\n$read\n+\n".("I" x length($read))."\n";
			close(BWAFQ);
			my $idx_cmd = "$bwa index $bwacolor .bwa.fa";
			print STDERR "$idx_cmd\n";
			system("$idx_cmd 2>/dev/null") == 0 || die;
			my $aln_cmd = "$bwa aln $bwacolor .bwa.fa .bwa.fq > .bwa.sai";
			print STDERR "$aln_cmd\n";
			system("$aln_cmd 2>/dev/null") == 0 || die;
			my $sam_cmd = "$bwa samse .bwa.fa .bwa.sai .bwa.fq";
			print STDERR "$sam_cmd\n";
			open(BWAOUT, "$sam_cmd 2>/dev/null |") || die;
			print "BWA output:\n";
			while(<BWAOUT>) {
				next if /^@/;
				my @s = split(/\t/);
				print $_ if $s[2] ne "*";
			}
		}
		if($novo ne "" && $color eq "") {
			-x "${novo}index" || die;
			-x "${novo}align" || die;
			open(NOVOFA, ">.novo.fa") || die;
			print NOVOFA ">0\n$ref\n";
			close(NOVOFA);
			open(NOVOFQ, ">.novo.fq") || die;
			print NOVOFQ "\@0\n$read\n+\n".("I" x length($read))."\n";
			close(NOVOFQ);
			my $idx_cmd = "${novo}index .novo.idx .novo.fa";
			print STDERR "$idx_cmd\n";
			system("$idx_cmd 2>/dev/null >/dev/null") == 0 || die;
			my $aln_cmd = "${novo}align -d .novo.idx -f .novo.fq -o SAM";
			print STDERR "$aln_cmd\n";
			open(NOVOOUT, "$aln_cmd 2>/dev/null |") || die;
			print "Novo output:\n";
			while(<NOVOOUT>) {
				next if /^@/;
				my @s = split(/\t/);
				print $_ if $s[2] ne "*";
			}
		}
	}
	$lastref = $ref;
	$lastcolor = $color;
}
print "\nPASSED\n\n";
