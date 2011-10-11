#!/usr/bin/env perl

##
# Copyright 2011, Ben Langmead <blangmea@jhsph.edu>
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

#
# The paired-end data is made by (a) changing to the reads subdirectory and (b)
# running 'perl simulate.pl --ref=../reference/lambda_virus.fa'.
#

#
# The long-read data is made by (a) changing to the reads subdirectory and (b)
# running 'perl simulate.pl --ref=../reference/lambda_virus.fa --long
# --unpaired --prefix=longreads'.
#

use strict;
use warnings;
use Carp;
use Math::Random qw(random_normal random_exponential);
use Getopt::Long;
use List::Util qw(max min);

my @fa_fn = ();       # files with reference FASTA
my $rf = "";          # reference sequence
my $long = 0;         # 1 -> generate long reads
my $paired = 1;       # 1 -> generate paired-end reads
my $prefix = "reads"; # output files start with this string

GetOptions (
	"fasta|reference=s" => \@fa_fn,
	"long"              => \$long,
	"unpaired"          => sub { $paired = 0; },
	"prefix=s"          => \$prefix
) || die "Bad option";

for my $fn (@fa_fn) {
	open(FN, $fn) || confess;
	my $name = "";
	while(<FN>) {
		chomp;
		$rf .= $_ if substr($_, 0, 1) ne ">";
	}
	close(FN);
}

my %revcompMap = (
	"A" => "T", "T" => "A", "a" => "t", "t" => "a",
	"C" => "G", "G" => "C", "c" => "g", "g" => "c",
	"R" => "Y", "Y" => "R", "r" => "y", "y" => "r",
	"M" => "K", "K" => "M", "m" => "k", "k" => "m",
	"S" => "S", "W" => "W", "s" => "s", "w" => "w",
	"B" => "V", "V" => "B", "b" => "v", "v" => "b",
	"H" => "D", "D" => "H", "h" => "d", "d" => "h",
	"N" => "N", "." => ".", "n" => "n" );

sub comp($) {
	my $ret = $revcompMap{$_[0]} || confess "Can't reverse-complement '$_[0]'";
	return $ret;
}

sub revcomp {
	my ($ret) = @_;
	$ret = reverse $ret;
	for(0..length($ret)-1) { substr($ret, $_, 1) = comp(substr($ret, $_, 1)); }
	return $ret;
}

my $nreads    = 10000;  # number of reads/end to generate
my $rdlen_av  = 75;     # average to use when drawing from exponential
my $rdlen_min = 40;     # minimum read length (added to exponential draw)
my $frag_av   = 250;    # mean fragment len
my $frag_sd   = 45;     # s.d. to use when drawing frag len from normal dist
my @fraglens  = ();     # fragment lengths (for paired)
my @readlens  = ();     # read/end lengths

if($long) {
	$nreads = 6000;
	$rdlen_av = 300;
	$rdlen_min = 40;
}

sub rand_dna($) {
	my $ret = "";
	for(1..$_[0]) { $ret .= substr("ACGT", int(rand(4)), 1); }
	return $ret;
}

#
# Mutate the reference
#

# Add SNPs
my $nsnp = 0;
for(0..length($rf)-1) {
	if(rand() < 0.0012) {
		my $oldc = substr($rf, $_, 1);
		substr($rf, $_, 1) = substr("ACGT", int(rand(4)), 1);
		$nsnp++ if substr($rf, $_, 1) ne $oldc;
	}
}

# Add microindels
my $newrf = "";
my $microgap = 0;
for(0..length($rf)-1) {
	$newrf .= substr($rf, $_, 1);
	if(rand() < 0.0010) {
		my $len = int(random_exponential(1, 3))+1;
		if(int(rand()) == 0) {
			# Remove a bit of the reference
			$len = min($len, length($newrf));
			$newrf = substr($newrf, 0, length($newrf)-$len);
		} else {
			# Add a bit to the reference
			$newrf .= rand_dna($len);
		}
		$microgap++;
	}
}
$rf = $newrf;

# Add some larger rearrangements
my $nrearr = int(random_exponential(1, 3)+1);
for(0..$nrearr) {
	my $break = int(rand(length($rf)));
	my $before = substr($rf, 0, $break);
	my $after = substr($rf, $break);
	$after = revcomp($after) if int(rand()) == 0;
	$rf = $after.$before;
}

print STDERR "Added $nsnp SNPs\n";
print STDERR "Added $microgap Microindels\n";
print STDERR "Added $nrearr Rearrangements\n";

#
# Simulate reads
#

if($paired) {
	# Pick random fragment and read lengths
	@fraglens = random_normal($nreads, $frag_av, $frag_sd);
	@readlens = random_exponential($nreads*2, $rdlen_av);
} else {
	# Pick random read lengths
	@readlens = random_exponential($nreads, $rdlen_av);
}
@fraglens = map int, @fraglens;
@readlens = map { int($_ + $rdlen_min) } @readlens;

sub rand_quals($) {
	my $ret = "";
	my $upper = (rand() < 0.2 ? 11 : 40);
	$upper = 4 if rand() < 0.02;
	for(1..$_[0]) {
		$ret .= chr(33+int(rand($upper)));
	}
	return $ret;
}

sub add_seq_errs($$) {
	my($rd, $qu) = @_;
	for(0..length($rd)-1) {
		my $c = substr($rd, $_, 1);
		my $q = substr($qu, $_, 1);
		$q = ord($q)-33;
		my $p = 10 ** (-0.1 * $q);
		if(rand() < $p) {
			$c = substr("ACGTNNNNNN", int(rand(10)), 1);
		}
		substr($rd, $_, 1) = $c;
		substr($qu, $_, 1) = $q;
	}
	return $rd;
}

# Now simulate 
my $rflen = length($rf);
if($paired) {
	open(RD1, ">${prefix}_1.fq") || die;
	open(RD2, ">${prefix}_2.fq") || die;
	for my $i (0..$#fraglens) {
		# Extract fragment
		my $flen = $fraglens[$i];
		my $off = int(rand($rflen - ($flen-1)));
		my $fstr = substr($rf, $off, $flen);
		# Possibly reverse complement
		$fstr = revcomp($fstr) if (int(rand(2)) == 0);
		# Get reads 1 and 2
		my $rdlen1 = min($readlens[2*$i], $flen);
		my $rdlen2 = min($readlens[2*$i+1], $flen);
		my $rd1 = substr($fstr, 0, $rdlen1);
		my $rd2 = substr($fstr, length($fstr)-$rdlen2);
		length($rd2) == $rdlen2 || die "Got ".length($rd2)." expected $rdlen2";
		# Reverse complement 2 to simulate --fr orientation
		$rd2 = revcomp($rd2);
		# Generate random quality values
		my $qu1 = rand_quals($rdlen1);
		$rd1 = add_seq_errs($rd1, $qu1);
		length($rd1) == length($qu1) || die;
		my $qu2 = rand_quals($rdlen2);
		$rd2 = add_seq_errs($rd2, $qu2);
		length($rd2) == length($qu2) || die;
		# Print
		print RD1 "\@r".($i+1)."\n$rd1\n+\n$qu1\n";
		print RD2 "\@r".($i+1)."\n$rd2\n+\n$qu2\n";
	}
	close(RD1);
	close(RD2);
	print STDERR "Made pairs: reads_1.fq/reads_2.fq\n";
} else {
	open(RD1, ">${prefix}.fq") || die;
	for my $i (0..$#readlens) {
		# Extract fragment
		my $rdlen = $readlens[$i];
		my $off = int(rand($rflen - ($rdlen-1)));
		my $rd = substr($rf, $off, $rdlen);
		# Possibly reverse complement
		$rd = revcomp($rd) if int(rand(2)) == 0;
		# Generate random quality values
		my $qu = rand_quals($rdlen);
		$rd = add_seq_errs($rd, $qu);
		length($rd) == length($qu) || die;
		# Print
		print RD1 "\@r".($i+1)."\n$rd\n+\n$qu\n";
	}
	close(RD1);
	print STDERR "Made unpaired reads: $prefix.fq\n";
}

print STDERR "DONE\n";
