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
# args.pl
#
# Basic tests to ensure that bad combinations of arguments are rejected
# and good ones are accepted.
#

my $bowtie2 = "./bowtie2";
my $bowtie2_d = "./bowtie2-debug";
if(system("$bowtie2 --version") != 0) {
	print STDERR "Could not execute ./bowtie2; looking in PATH...\n";
	$bowtie2 = `which bowtie2`;
	chomp($bowtie2);
	if(system("$bowtie2 --version") != 0) {
		die "Could not find bowtie2 in current directory or in PATH\n";
	}
}
if(system("$bowtie2_d --version") != 0) {
	print STDERR "Could not execute ./bowtie2-debug; looking in PATH...\n";
	$bowtie2_d = `which bowtie2-debug`;
	chomp($bowtie2_d);
	if(system("$bowtie2_d --version") != 0) {
		die "Could not find bowtie2-debug in current directory or in PATH\n";
	}
}

if(! -f "e_coli_c.1.ebwt") {
	print STDERR "Making colorspace e_coli index\n";
	my $bowtie2_build = "./bowtie2-build";
	if(system("$bowtie2_build --version") != 0) {
		print STDERR "Could not execute ./bowtie2-build; looking in PATH...\n";
		$bowtie2_build = `which $bowtie2_build`;
		chomp($bowtie2_build);
		if(system("$bowtie2_build --version") != 0) {
			die "Could not find bowtie2-build in current directory or in PATH\n";
		}
	}
	system("$bowtie2_build -C genomes/NC_008253.fna e_coli_c") && die;
} else {
	print STDERR "Colorspace e_coli index already present...\n";
}

open TMP, ">.args.pl.1.fa" || die;
print TMP ">\nT0120012002012030303023\n";
close(TMP);
open TMP, ">.args.pl.1.qv" || die;
print TMP ">\n10 11 12 10 10 11 12 10 10 12 10 22 33 23 13 10 12 23 24 25 26 27\n";
close(TMP);
open TMP, ">.args.pl.2.qv" || die;
print TMP ">\n9 10 11 12 10 10 11 12 10 10 12 10 22 33 23 13 10 12 23 24 25 26 27\n";
close(TMP);

my @bad = (
	"-N 6",
	"-N 5",
	"-N 4",
	"-N 3"
);

my @badEx = (
	"e_coli -f .args.pl.1.fa -Q .args.pl.1.qv",
	"e_coli_c -f .args.pl.1.fa -Q .args.pl.1.qv"
);

my @good = (
	"-N 0",
	"-N 1",
	"-N 2"
);

my @goodEx = (
	"-C e_coli_c -f .args.pl.1.fa -Q .args.pl.1.qv",
	"-C e_coli_c -f .args.pl.1.fa -Q .args.pl.2.qv"
);

sub run($) {
	my $cmd = shift;
	print "$cmd\n";
	return system($cmd);
}

print "Bad:\n";
for my $a (@bad) {
	run("$bowtie2 $a e_coli reads/e_coli_1000.fq /dev/null") != 0 ||
		die "bowtie2 should have rejected: \"$a\"\n";
	run("$bowtie2_d $a e_coli reads/e_coli_1000.fq /dev/null") != 0 ||
		die "bowtie2-debug should have rejected: \"$a\"\n";
	print "PASSED: bad args \"$a\"\n";
}
print "\nBadEx:\n";
for my $a (@badEx) {
	run("$bowtie2 $a /dev/null") != 0 ||
		die "bowtie2 should have rejected: \"$a\"\n";
	run("$bowtie2_d $a /dev/null") != 0 ||
		die "bowtie2-debug should have rejected: \"$a\"\n";
	print "PASSED: bad args \"$a\"\n";
}
print "\nGood:\n";
for my $a (@good) {
	run("$bowtie2 $a e_coli reads/e_coli_1000.fq /dev/null") == 0 ||
		die "bowtie2 should have accepted: \"$a\"\n";
	run("$bowtie2_d $a e_coli reads/e_coli_1000.fq /dev/null") == 0 ||
		die "bowtie2-debug should have accepted: \"$a\"\n";
	print "PASSED: good args \"$a\"\n";
}
print "\nGoodEx:\n";
for my $a (@goodEx) {
	run("$bowtie2 $a /dev/null") == 0 ||
		die "bowtie2 should have accepted: \"$a\"\n";
	run("$bowtie2_d $a /dev/null") == 0 ||
		die "bowtie2-debug should have accepted: \"$a\"\n";
	print "PASSED: good args \"$a\"\n";
}
