#!/usr/bin/perl -w

#
# Verifies that Bowtie's --better and --best alignment modes really
# produce alignments that are the best possible in terms of score.
#
# Run this from the Bowtie directory.
#
# E.g.: perl ../scripts/best_verify.pl -v 0
#       perl ../scripts/best_verify.pl -v 0 e_coli genomes/NC_008253.fna reads/e_coli_1000.fq
#

use warnings;
use strict;
use Getopt::Long;
Getopt::Long::Configure ("no_ignore_case");

my %options=();
my $varg;
my $narg;
my $g = undef;
my $l = undef;
my $e = undef;
my $C = undef;
my $O = undef;
my $E = undef;
my $debug = 0;
my $round = "";
my $nomaqround = 0;
my $result = GetOptions("v=i" => \$varg,
                        "n=i" => \$narg,
                        "l=i" => \$l,
                        "e=i" => \$e,
                        "d"   => \$debug,
                        "g"   => \$g,
                        "O"   => \$O,
                        "E"   => \$E,
                        "C"   => \$C,
                        "nomaqround" => \$nomaqround) || die "Bad options";

my $match_mode = "-n 2";
$match_mode = "-v " . $varg if defined($varg);
$match_mode = "-n " . $narg if defined($narg);
$match_mode .= " -l $l " if defined($l);
$match_mode .= " -e $e " if defined($e);
$match_mode .= " -g " if defined($g);
$match_mode .= " -C " if defined($C);
$match_mode .= " -O $O " if defined($O);
$match_mode .= " -E $E " if defined($E);
$round = "--nomaqround" if $nomaqround;

# These should match Bowtie defaults
$O = 100 unless defined($O);
$E = 40  unless defined($E);

print "Maq-like rounding is: ".($nomaqround ? "off" : "on") . "\n";
print "Using match mode: $match_mode\n";

if(defined($C)) {
	print "Aborting with exitlevel 0 because colorspace alignment was selected\n";
	exit 0;
}

my $bowtie_dir = ".";
my $bowtie_exe = "bowtie";
$bowtie_exe .= "-debug" if $debug;

my $index  = "e_coli";
$index = $ARGV[0] if defined($ARGV[0]);
my $reads = "reads/e_coli_1000.fq";
$reads = $ARGV[1] if defined($ARGV[1]);

if($reads =~ /\.fa$/) {
	$reads = " -f $reads ";
} elsif($reads =~ /\.raw$/) {
	$reads = " -r $reads ";
} elsif($reads =~ /\.integer.fq$/) {
	$reads = " -q --integer-quals $reads ";
} elsif($reads =~ /\.fq$/) {
	$reads = " -q $reads ";
} else {
	$reads = " -c $reads ";
}

my $seedLen = ($l ? $l : 28);
my $vmode = ($match_mode =~ /[-]v/);

system("make -C $bowtie_dir $bowtie_exe") == 0 || die;

# Run Bowtie to get best alignments
my $bowtie_best_cmd = "$bowtie_dir/$bowtie_exe -y --cost -l $seedLen $round $match_mode --refidx $index $reads";

# Run Bowtie to get all alignments
my $bowtie_all_cmd = "$bowtie_dir/$bowtie_exe -y --cost -l $seedLen $round $match_mode -a --refidx $index $reads";

sub costFromMms($$$) {
	my ($mmstr, $fw, $quals) = @_;
	my $cost = 0;
	my $len = length($quals);
	my @mms = split(/,/, $mmstr);
	my $lastRefGap = -2;
	for my $mm (@mms) {
		my @mmss = split(/:/, $mm);
		my $mmoff = int($mmss[0]);
		my $affectsStratum = 0;
		if($mmoff < $seedLen) {
			$affectsStratum = 1;
			$cost += (1 << 14);
		}
		if(!$fw) {
			$mmoff = $len - $mmoff - 1;
		}
		my @mmchrs = split(/>/, $mmss[1]);
		$#mmchrs == 1 || die;
		if($mmchrs[0] eq '-') {
			# Reference has a gap
			if(abs($lastRefGap - int($mmss[0])) == 1) {
				$cost += $E; # extension of reference gap
			} else {
				$cost += $O; # opening of reference gap
			}
			$lastRefGap = int($mmss[0]);
		} elsif($mmchrs[1] eq '-') {
			$lastRefGap = -2;
			$cost += $O;
			if(length($mmchrs[0]) > 1) {
				my $extensions = length($mmchrs[0]) - 1;
				$cost += ($E * $extensions);
				$cost += ((1 << 14) * $extensions) if $affectsStratum;
			}
		} else {
			# Grab the quality from the quality string
			$lastRefGap = -2;
			my $q = substr($quals, $mmoff, 1);
			$q = int(ord($q) - 33);
			if(!$nomaqround) {
				$q = int(int($q + 5) / 10);
				$q = 3 if $q > 3;
				$q *= 10;
			}
			my $qcost += $q;
			$cost += $qcost;
		}
	}
	return $cost;
}

print "$bowtie_best_cmd\n";
open BOWTIE_BEST, "$bowtie_best_cmd |";
my %nameToBestScore = ();
my %nameToBestAlignment = ();
my $bestAls = 0;
while(<BOWTIE_BEST>) {
	next if /^Reported/;
	chomp;
	my $line = $_;
	my @ls = split(/[\t]/, $line);
	$#ls >= 5 || die "Alignment not formatted correctly: $line";
	my $name = $ls[0];
	defined($nameToBestAlignment{$name}) && die "Read with name $name appeared more than once in best-hit output";
	defined($nameToBestScore{$name}) && die "Read with name $name appeared more than once in best-hit output";
	my $len = length($ls[4]);
	my $quals = $ls[5];
	my $mmstr = "";
	$mmstr = $ls[7] if defined($ls[7]);
	my $cost = 0;
	my $fw = $ls[1] eq "+";
	$seedLen = $len if $vmode;
	$cost = costFromMms($mmstr, $fw, $quals) if $mmstr ne "-";
	$ls[-2] == int($ls[-2]) || die "Bad stratum in second-to-last column: $ls[-2]\n";
	$ls[-1] == int($ls[-1]) || die "Bad cost in last column: $ls[-1]\n";
	my $btcost = $ls[-1] + ($ls[-2] << 14);
	print "$line: $cost\n";
	$cost == $btcost || die "best_verify.pl calculated cost of $cost, but bowtie calculated $btcost\n";
	$nameToBestAlignment{$name} = $line;
	$nameToBestScore{$name} = $cost;
	$bestAls++;
}
close(BOWTIE_BEST);
if($? != 0) {
	die "bowtie --best process failed with exitlevel $?\n";
}

print "$bowtie_all_cmd\n";
open BOWTIE_ALL, "$bowtie_all_cmd |";
my $allAls = 0;
while(<BOWTIE_ALL>) {
	next if /^Reported/;
	chomp;
	my $line = $_;
	my @ls = split(/[\t]/, $line);
	$#ls >= 5 || die "Alignment not formatted correctly: $line";
	my $name = $ls[0];
	my $len = length($ls[4]);
	my $quals = $ls[5];
	my $mmstr = "";
	$mmstr = $ls[7] if defined($ls[7]);
	my $cost = 0;
	my $fw = $ls[1] eq "+";
	$seedLen = $len if $vmode;
	$cost = costFromMms($mmstr, $fw, $quals) if $mmstr ne "-";
	$ls[-2] == int($ls[-2]) || die "Bad stratum in second-to-last column: $ls[-2]\n";
	$ls[-1] == int($ls[-1]) || die "Bad cost in last column: $ls[-1]\n";
	my $btcost = $ls[-1] + ($ls[-2] << 14);
	print "$line: $cost\n";
	$cost == $btcost || die "best_verify.pl calculated cost of $cost, but bowtie calculated $btcost\n";
	defined($nameToBestAlignment{$name}) ||
		die "Read with alignment:\n$line\nhas no corresponding alignment in best-hit mode\n";
	int($cost) >= int($nameToBestScore{$name}) ||
		die "Alignment:\n$line\n$cost\nis better than:\n$nameToBestAlignment{$name}\n$nameToBestScore{$name}\n";
	$allAls++;
}
close(BOWTIE_ALL);
if($? != 0) {
	die "bowtie -a process failed with exitlevel $?\n";
}

print "Checked $bestAls best-alignments against $allAls all-alignments\n";
print "PASSED\n";
