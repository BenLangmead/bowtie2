#!/usr/bin/perl -w

##
# fastq_to_tabbed.pl
#
# Convert a pair of FASTQ files to a single file in a one-pair-per-line
# format, where each line has these five tab-delimited fields:
#
# 1. Pair name
# 2. Sequence of mate 1
# 3. Qualities of mate 1
# 4. Sequence of mate 2
# 5. Qualities of mate 2
#
# This script comes in handy if (a) you'd just like to store your
# paired-end data in a less awkward format than the usual pair of
# parallel FASTQ files, or (b) you'd like to use Bowtie with gzipped
# paired-end input without unzipping it first.  In that case, you can
# pipe the output of this script (which handles gzipped inputs by
# piping them through gzip -dc) to Bowtie and use '--12 -'.
#
# Note that this script can also handle unpaired input (with -u), which
# Bowtie handles approrpaitely in --12 mode even when it's intermingled
# with paired-end input.
#

use strict;
use warnings;
use Getopt::Long;

# Courtesy: http://www.perlmonks.org/?node_id=625977
sub shuffle (@) {
	my @a = \(@_);
	my $n;
	my $i = @_;
	map {
		$n = rand($i--);
		(${$a[$n]}, $a[$n] = $a[$i])[0];
	} @_;
}

my $unpaired = "";
my $mate1 = "";
my $mate2 = "";
my $shuffle = 0;

GetOptions ("u=s"     => \$unpaired,
            "1=s"     => \$mate1,
            "2=s"     => \$mate2,
            "shuffle" => \$shuffle) || die "Bad option";

my @output;
my @unpaireds = split(/,/, $unpaired);

for my $f (@unpaireds) {
	open UNP, "$f" || die;
	while(<UNP>) {
		chomp;
		my $name = $_;
		$name = substr($name, 1);
		my $seq = <UNP>;
		chomp($seq);
		my $name2 = <UNP>;
		my $qual = <UNP>;
		chomp($qual);
		push @output, "$name\t$seq\t$qual\n";
	}
	close(UNP);
}

my @mate1s = split(/,/, $mate1);
my @mate2s = split(/,/, $mate2);

for(my $i = 0; $i <= $#mate1s; $i++) {
	if($mate1s[$i] =~ /\.gz$/) {
		open M1, "gzip -dc $mate1s[$i] |" || die;
	} else {
		open M1, "$mate1s[$i]" || die;
	}
	if($mate2s[$i] =~ /\.gz$/) {
		open M2, "gzip -dc $mate2s[$i] |" || die;
	} else {
		open M2, "$mate2s[$i]" || die;
	}
	while(<M1>) {
		my $name1 = $_;
		chomp($name1);
		$name1 = substr($name1, 1, -2);
		my $name2 = <M2>;
		chomp($name2);
		$name2 = substr($name2, 1);
		my $seq1 = <M1>;
		chomp($seq1);
		my $seq2 = <M2>;
		chomp($seq2);
		my $tmp = <M1>;
		$tmp = <M2>;
		my $qual1 = <M1>;
		chomp($qual1);
		my $qual2 = <M2>;
		chomp($qual2);
		print "$name1\t$seq1\t$qual1\t$seq2\t$qual2\n" unless $shuffle;
		push @output, "$name1\t$seq1\t$qual1\t$seq2\t$qual2\n" if $shuffle;
	}
	close(M1);
	close(M2);
}

if($shuffle) {
	@output = shuffle(@output);
	print join("", @output);
}
