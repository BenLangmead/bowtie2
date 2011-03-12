#!/usr/bin/perl -w

#
# Verifies that Bowtie's paired-end mode gives alignments that are
# consistent with the alignments produced in single-end mode with -a
# and --nostrata options.
#
# Run this from the Bowtie directory.
#
# Usage: perl ../scripts/pe_verify.pl [mode] [index] [mates1] [nates2]
#
# Defaults are:
#  [mode] = -n 2
#  [index] = e_coli
#  [mates1] = reads/e_coli_1000_1.fq
#  [mates2] = reads/e_coli_1000_2.fq
#
# E.g.: perl ../scripts/pe_verify.pl -v 0
#       perl ../scripts/pe_verify.pl -v 0 e_coli reads/e_coli_1000_1.fq reads/e_coli_100000_e1_2.fq
#

use warnings;
use strict;
use List::Util qw[min max];
use Getopt::Long;
Getopt::Long::Configure ("no_ignore_case");

my $l = undef;
my $n = undef;
my $e = undef;
my $v = undef;
my $I = undef;
my $X = undef;
my $d = undef;
my $C = undef;
my $g = undef;
my $s = undef;
my $u = undef;
my $colCseq = undef;
my $colCqual = undef;
my $args = "";
my $verbose = undef;
my $m1fw = 1;
my $m2fw = 0;
my $reffa = undef;
my %rf = ();

GetOptions ("I=i" => \$I,
            "X=i" => \$X,
            "v=i" => \$v,
            "n=i" => \$n,
            "e=i" => \$e,
            "l=i" => \$l,
            "d"   => \$d,
            "g"   => \$g,
            "s=i" => \$s,
            "u=i" => \$u,
            "C"   => \$C,
            "reference:s" => \$reffa,
            "verbose" => \$verbose,
            "col-cseq" => \$colCseq,
            "col-cqual" => \$colCqual,
            "args:s" => \$args,
            ) || die "One or more errors parsing script arguments";

my $inner = 0;
my $outer = 250;
$inner = $I if ${I};
$outer = $X if ${X};

defined($reffa) || die "Must specify fasta file for reference";
open (FA, $reffa) || die "Could not open -reference file $reffa\n";
my $name = "";
while(<FA>) {
	chomp;
	if(/^>/) {
		$name = substr($_, 1);
		$rf{$name} = "";
	} else {
		$name ne "" || die;
		$rf{$name} .= $_;
	}
}
close(FA);

# Utility function that returns the complement of its argument
sub comp {
	my $r = shift;
	my $c = shift;
	$r =~ tr/aAcCgGtT/tTgGcCaA/ unless (defined($c) && $c);
	return $r;
}

# Insert $src just before the $off'th character of $dst
sub ins {
	my ($dst, $src, $off) = @_;
	$off < length($dst) || die "ins offset $off is not less than dst length ".length($dst)."\n";
	substr($dst, $off, 0) = $src;
	return $dst;
}

# Delete the $off'th character of $dst
sub del {
	my ($dst, $off) = @_;
	$off < length($dst) || die "del offset $off is not less than dst length ".length($dst)."\n";
	substr($dst, $off, 1) = "";
	return $dst;
}

# Substitute $src for the $off'th character of $dst
sub subst {
	my ($dst, $src, $off) = @_;
	substr($dst, $off, 1) = $src;
	return $dst;
}

# Utility function that returns the reverse complement of its argument
sub revcomp {
	my $r = shift;
	my $c = shift; # may be undef
	$r = reverse($r);
	$r =~ tr/aAcCgGtT/tTgGcCaA/ unless (defined($c) && $c);
	return $r;
}

##
# Sanity check whether a read and a list of edits corresponds correctly
# to a substring of the reference.
#
sub checkAlignmentRef {
	my ($ref, $read, $fw, $off, $edits) = @_;
	my $orig = $read;
	my $origr = $fw ? $orig : revcomp($orig);
	$read = $origr;
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
				my $inspos = $pos + $adjust;
				$read = ins($read, $ref, $inspos);
				$adjust += length($ref);
			} elsif($ref eq "-") {
				# deletion
				my $delpos = $pos + $adjust;
				$read = del($read, $delpos);
				$adjust--;
			} else {
				# mismatch
				my $mmpos = $pos + $adjust;
				$read = subst($read, $ref, $mmpos);
			}
		}
	}
	my $reado = $read;
	$read = revcomp($read) unless $fw;
	my $rstr = substr($ref, $off, length($read));
	$fw = $fw ? '+' : '-';
	$read eq $rstr || die "\n\nAlignment failed to pass sanity check:\n".
	"FW: $fw, Off: $off, Edits: $edits\n".
	"OrigFW: $orig\n".
	"Orig:   $origr\n".
	"Qry:    $reado\n".
	"QryFW:  $read\n".
	"Ref:    $rstr\n";
}


my $extra_args = "";
$extra_args = $e if defined($e);
my $match_mode = "-n 2";
$match_mode = "-v " . $v if defined($v);
$match_mode = "-n " . $n if defined($n);
$match_mode .= " -l " . $l if defined($l);
$match_mode .= " -e " . $e if defined($e);
$match_mode .= " -C" if defined($C);
$match_mode .= " -g" if defined($g);
$match_mode .= " -u $u " if defined($u);
$match_mode .= " -s $s " if defined($s);
$match_mode .= " --col-cseq" if defined($colCseq);
$m2fw = 1 if $C;

print "Using match mode: $match_mode\n";

my $bowtie_dir = ".";
my $bowtie_exe = "bowtie";
$bowtie_exe .= "-debug" if $d;

my $index  = "e_coli";
$index = $ARGV[0] if defined($ARGV[0]);
my $reads1 = "reads/e_coli_1000_1.fq";
$reads1 = $ARGV[1] if (defined($ARGV[1]));
my $reads2 = "reads/e_coli_1000_2.fq";
$reads2 = $ARGV[2] if (defined($ARGV[2]));

# Infer input type so we can provide Bowtie with appropriate option
if($reads1 =~ /\.fa/) {
	$reads2 =~ /\.fa/ || die "Reads files $reads1 and $reads2 have different extensions";
	$extra_args .= " -f ";
} elsif($reads1 =~ /\.raw/) {
	$reads2 =~ /\.raw/ || die "Reads files $reads1 and $reads2 have different extensions";
	$extra_args .= " -r ";
} elsif(!($reads1 =~ /\.fq/)) {
	(!($reads2 =~ /\.fq/)) || die "Reads files $reads1 and $reads2 have different extensions";
	$extra_args .= " -c ";
}

my $sesMustHavePes = 0; # force se pairs to have corresponding pe pairs 

system("make -C $bowtie_dir $bowtie_exe") == 0 || die;

# Run Bowtie not in paired-end mode for first mate file
my $bowtie_se_cmd1 = "$bowtie_dir/$bowtie_exe $match_mode $args -y $extra_args -a --refidx $index $reads1";

# Run Bowtie not in paired-end mode for second mate file
my $bowtie_se_cmd2 = "$bowtie_dir/$bowtie_exe $match_mode $args -y $extra_args -a --refidx $index $reads2";

# Run Bowtie in paired-end mode
my $bowtie_pe_cmd = "$bowtie_dir/$bowtie_exe $match_mode $args -I $inner -X $outer -y $extra_args -a --refidx $index -1 $reads1 -2 $reads2";
print "$bowtie_pe_cmd --chunkmbs 1024\n";
open BOWTIE_PE, "$bowtie_pe_cmd --chunkmbs 1024 |";

if($C) {
	$inner = max(0, $inner-1);
	$outer = max(0, $outer-1);
}

##
# Given leftmost reference char, mismatch string, and length of read,
# calculate the rightmost character covered by the alignment.
#
sub calcRightExtent($$$) {
	my ($left, $mm, $len) = @_;
	my $adj = 0;
	for my $m (split(/,/, $mm)) {
		if($m =~ /->/) {
			$adj--;
		} elsif($m =~ />-/) {
			$m =~ /:(.*)>/;
			#print "  subtracting ".length($1)." (length($1))\n";
			$adj += length($1);
		}
	}
	my $ret = $left + $len + $adj;
	#print "$left, $mm, $len -> $ret\n";
	return $ret;
}

##
# Given left and right extends for two alignments, calculate whether
# one is contained in the other.
#
sub contained($$$$) {
	my ($l1, $r1, $l2, $r2) = @_;
	#print "[$l1, $r1], [$l2, $r2]\n";
	return ($l1 >= $l2 && $r1 <= $r2) || ($l2 >= $l1 && $r2 <= $r1);
}

my $pes = 0;
my $pesFw = 0;
my $pesRc = 0;
my %peHash = ();
while(<BOWTIE_PE>) {
	chomp;
	my $l1 = $_;
	my $l2 = <BOWTIE_PE>;
	chomp($l2);
	print "$l1\n$l2\n";
	my @l1s = split(/[\t]/, $l1);
	my @l2s = split(/[\t]/, $l2);
	$#l1s >= 5 || die "Paired-end alignment not formatted correctly: $l1";
	$#l2s >= 5 || die "Paired-end alignment not formatted correctly: $l2";
	# Split the read name
	my @l1rs = split(/\//, $l1s[0]);
	my @l2rs = split(/\//, $l2s[0]);
	$#l1rs >= 1 || die "Read not formatted correctly: $l1s[0]";
	$#l2rs >= 1 || die "Read not formatted correctly: $l2s[0]";
	$l1rs[0] eq $l2rs[0] || die "Before-/ parts of read names don't match: $l1rs[0], $l2rs[0]";
	my $mate1 = ($l1rs[$#l1rs] eq "1");
	my $mate1str = $mate1 ? "1" : "0";
	my $basename = $l1rs[0];
	my $lleft = int($l1s[3]);
	my $rleft = int($l2s[3]);
	my $lright = calcRightExtent($lleft, $l1s[7], length($l1s[4]));
	my $rright = calcRightExtent($rleft, $l2s[7], length($l2s[4]));
	if(contained($lleft, $lright, $rleft, $rright)) {
		print "Skipping above due to containment\n";
		next;
	}
	my $insLen = $rleft - $lleft + length($l2s[4]);
	$insLen > length($l1s[4]) || die "Insert length $insLen did not exceed first mate length ".length($l1s[4]);
	$insLen > length($l1s[5]) || die "Insert length $insLen did not exceed first mate length ".length($l1s[5]);
	$insLen >= $inner || die "Insert length was $insLen < $inner\n";
	$insLen <= $outer || die "Insert length was $insLen > $outer\n";
	my $read1Short = $l1s[4];
	my $qual1Short = $l1s[5];
	my $read2Short = $l2s[4];
	my $qual2Short = $l2s[5];
	my $read1Mms = "";
	$read1Mms = $l1s[7] if defined($l1s[7]);
	$read1Mms = "-" if $read1Mms eq "";
	my $read2Mms = "";
	$read2Mms = $l2s[7] if defined($l2s[7]);
	$read2Mms = "-" if $read2Mms eq "";
	my $content = "$read1Short:$qual1Short $read2Short:$qual2Short";
	$content = "" if $C && !$colCseq;
	my $mcont = "$read1Mms $read2Mms";
	$mcont = "" if $C;
	my $val = "$basename $mate1str $l1s[2] $l1s[1] $l1s[3] $l2s[1] $l2s[3] $content $mcont";
	#defined ($peHash{$basename}) && die "Already saw paired-end alignment for basename $basename";
	$peHash{$basename} = $val;
	if($mate1) { $pesFw++; } else { $pesRc++; }
	$pes++;
	{
		my ($fw, $ref, $off, $seq, $mms) = ($l1s[1] eq '+', $l1s[2], $l1s[3], $l1s[4], $l1s[7]);
		checkAlignmentRef($rf{$ref}, $seq, $fw, $off, $mms, 1);
		($fw, $ref, $off, $seq, $mms) = ($l2s[1] eq '+', $l2s[2], $l2s[3], $l2s[4], $l2s[7]);
		checkAlignmentRef($rf{$ref}, $seq, $fw, $off, $mms, 1);
	}
}
close(BOWTIE_PE);

my $ses = 0;
my %seHash = ();
my %seMatedHash = ();
my %unmatchedSe = ();
my $unmatchedSes = 0;
my $unmatchedSeReads = 0;

print "$bowtie_se_cmd1 --redundants 2 --chunkmbs 1024\n";
open BOWTIE_SE1, "$bowtie_se_cmd1 --redundants 2 --chunkmbs 1024 |";
while(<BOWTIE_SE1>) {
	print "$_";
	chomp;
	$ses++;
	my @ls = split(/[\t]/, $_);
	my @lrs = split(/\//, $ls[0]);
	my $basename = $lrs[0];
	my $ref = $ls[2];
	my $off = int($ls[3]);
	my $len = length($ls[4]);
	my $readShort = $ls[4];
	my $qualShort = $ls[5];
	my $mms = "";
	$mms = $ls[7] if defined($ls[7]);
	$mms = "-" if $mms eq "";
	my $content = "$readShort $qualShort $mms";
	my $key = "$ref $ls[1] $off $len $content";
	push @{ $seHash{$basename}{$ref} }, $key;
	{
		my ($fw, $ref, $off, $seq, $mms) = ($ls[1] eq '+', $ls[2], $ls[3], $ls[4], $ls[7]);
		checkAlignmentRef($rf{$ref}, $seq, $fw, $off, $mms, 1);
	}
}
close(BOWTIE_SE1);

print "$bowtie_se_cmd2 --redundants 2 --chunkmbs 1024\n";
open BOWTIE_SE2, "$bowtie_se_cmd2 --redundants 2 --chunkmbs 1024 |";
open UNMATCHED_SE, ">.unmatched.se";
while(<BOWTIE_SE2>) {
	print "$_";
	chomp;
	$ses++;
	my @ls = split(/[\t]/, $_);
	my @lrs = split(/\//, $ls[0]);
	my $basename = $lrs[0];
	my $ref = $ls[2];
	my $off = int($ls[3]);
	my $len = length($ls[4]);
	my $fw = $ls[1] eq "+" ? 1 : 0;
	my $readShort = $ls[4];
	my $qualShort = $ls[5];
	my $mms = "";
	$mms = $ls[7] if defined($ls[7]);
	$mms = "-" if $mms eq "";
	my $content = "$readShort $qualShort";
	my $mcont = "$mms";
	my $key = "$ref $ls[1] $off $len $content $mcont";
	{
		checkAlignmentRef($rf{$ref}, $readShort, $fw, $off, $mms, 1);
	}
	# Is the other mate already aligned?
	if(defined($seHash{$basename}{$ref})) {
		# Get all of the alignments for the mate
		for my $om (@{ $seHash{$basename}{$ref} }) {
			my @oms = split(/ /, $om);
			$#oms == 6 || die "Wrong number of elements for oms: $#oms";
			my ($oref, $ofw, $ooff, $olen, $oreadShort, $oqualShort) = @oms;
			$ofw = $oms[1] eq "+" ? 1 : 0;
			$ooff = int($oms[2]);
			$olen = int($oms[3]);
			my $omms = "";
			print "Trying $ref:$off and $oref:$ooff\n" if $verbose;
			$omms = $oms[6] if defined($oms[6]);
			$oref eq $ref || die "Refs don't match: $oref, $ref";
			my $diff;
			my $peKey = "$basename ";
			if($ooff > $off) {
				# The #1 mate is on the right
				my $my_m1fw = $m1fw ? 0 : 1;
				my $my_m2fw = $m2fw ? 0 : 1;
				$diff = $ooff - $off + $olen;
				if ($diff <= $olen || $diff <= $len) {
					print "diff $diff is <= $olen and $len\n" if $verbose;
					next;
				}
				# upstream mate contains downstream one?
				#next if $off + $len >= $ooff + $olen;
				# mates are at the same position?
				if($ooff == $off) {
					print "overlapping offsets: $ooff, $off\n" if $verbose;
					next;
				}
				if ($diff < $inner || $diff > $outer) {
					print "diff $diff is outside of inner/outer: [$inner, $outer]\n" if $verbose;
					next;
				}
				if($ofw != $my_m1fw) {
					print "orientation of other $ofw doesn't match expected $my_m1fw\n" if $verbose;
					next;
				}
				if($fw != $my_m2fw) {
					print "orientation of anchor $fw doesn't match expected $my_m2fw\n" if $verbose;
					next;
				}
				$content = "$readShort:$qualShort $oreadShort:$oqualShort";
				$content = "" if $C && !$colCseq;
				$mcont = "$mms $omms";
				$mcont = "" if $C;
				my $fwc = $fw ? "+":"-";
				my $ofwc = $ofw ? "+":"-";
				my $right = calcRightExtent($off, $mms, $len);
				my $oright = calcRightExtent($ooff, $omms, $olen);
				$peKey .= "0 $ref $fwc $off $ofwc $ooff $content $mcont";
				if(contained($off, $right, $ooff, $oright)) {
					print "Skipping [$peKey] due to containment\n";
					next;
				}
			} else {
				# The #1 mate is on the left
				$diff = $off - $ooff + $len;
				if ($diff <= $olen || $diff <= $len) {
					print "diff $diff is <= $olen and $len\n" if $verbose;
					next;
				}
				# upstream mate contains downstream one?
				#next if $ooff + $olen >= $off + $len;
				# mates are at the same position?
				if($ooff == $off) {
					print "overlapping offsets: $ooff, $off\n" if $verbose;
					next;
				}
				if ($diff < $inner || $diff > $outer) {
					print "diff $diff is outside of inner/outer: [$inner, $outer]\n" if $verbose;
					next;
				}
				if($ofw != $m1fw) {
					print "orientation of other $ofw doesn't match expected $m1fw\n" if $verbose;
					next;
				}
				if($fw != $m2fw) {
					print "orientation of anchor $fw doesn't match expected $m2fw\n" if $verbose;
					next;
				}
				$content = "$oreadShort:$oqualShort $readShort:$qualShort";
				$content = "" if $C && !$colCseq;
				$mcont = "$omms $mms";
				$mcont = "" if $C;
				my $fwc = $fw ? "+":"-";
				my $ofwc = $ofw ? "+":"-";
				$peKey .= "1 $ref $ofwc $ooff $fwc $off $content $mcont";
				my $right = calcRightExtent($off, $mms, $len);
				my $oright = calcRightExtent($ooff, $omms, $olen);
				if(contained($ooff, $oright, $off, $right)) {
					print "Skipping [$peKey] due to containment\n";
					next;
				}
			}
			# Found a legitimate paired-end alignment using a pair of
			# single-end alignments
			if($seMatedHash{$basename}) {
				print "already found corresponding paired-end\n" if $verbose;
				next;
			}
			if($sesMustHavePes) {
				defined($peHash{$basename}) ||
					die "Found single-end alignment for $basename, but no paired-end";
			} else {
				if(!defined($peHash{$basename})) {
					if(!defined($unmatchedSe{$basename})) {
						$unmatchedSe{$basename} = 0;
						$unmatchedSeReads++;
					}
					$unmatchedSe{$basename}++;
					$unmatchedSes++;
					print UNMATCHED_SE "Read $basename:\n$om\n$key\n";
				}
			}
			if(defined($peHash{$basename}) && $peHash{$basename} eq $peKey) {
				delete $peHash{$basename};
				$seMatedHash{$basename} = 1;
			} else {
				print "No matchup:\n$peHash{$basename}\n$peKey\n" if $verbose;
			}
			#print "Found alignment for mate $otherMate of $ls[0]; diff: $diff\n";
		}
	}
}
close(BOWTIE_SE2);
close(UNMATCHED_SE);

my $die = 0;
for my $peKey (keys %peHash) {
	print "Paired-end $peKey has a paired-end alignment without single-end support\n";
	print "[ $peHash{$peKey} ]\n";
	$die++;
}

if($unmatchedSes > 0) {
	print "Total of $unmatchedSes unmatched single-end alignments found for $unmatchedSes distinct pairs\n";
	print "Ref orientation off len seq quals mms\n";
	system("cat .unmatched.se");
	$die++;
}

$die && die "Found $die paired-end reads with no corresponding single-end mates";
print "PASSED; analyzed $pes paired-end ($pesFw fw, $pesRc rc) and $ses single-end alignments\n";
