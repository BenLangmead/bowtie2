#!/usr/bin/perl -w

##
# long_read.pl
#
# Basic tests to ensure that long reads are handled properly.
#

use strict;
use warnings;

my $bowtie   = "./bowtie";
my $bowtie_d = "./bowtie-debug";
if(system("$bowtie --version") != 0) {
	$bowtie = `which bowtie`;
	chomp($bowtie);
	if(system("$bowtie --version") != 0) {
		die "Could not find bowtie in current directory or in PATH\n";
	}
}
if(system("$bowtie_d --version") != 0) {
	$bowtie_d = `which bowtie-debug`;
	chomp($bowtie_d);
	if(system("$bowtie_d --version") != 0) {
		die "Could not find bowtie-debug in current directory or in PATH\n";
	}
}

if(! -f "e_coli_c.1.bt2") {
	print STDERR "Making colorspace e_coli index\n";
	system("make bowtie-build") && die;
	system("bowtie-build -C genomes/NC_008253.fna e_coli_c") && die;
} else {
	print STDERR "Colorspace e_coli index already present...\n";
}

my @reads = (
	# 70:
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGACCACAAA:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????",
	# 140:
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGACCACAAA".
	"CTGGAACGGCACGTTCGGTCAGCATTTGTACGGTTTCCACCAGCCACTCACCGCCTTCAATTTTGACCAT:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????".
	"??????????=====;96666/.......;=>>>?DDDDDDDDDDDDDDDDDDDDDBBBBDDDDDCCCDD",
	# 210:
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGACCACAAA".
	"CTGGAACGGCACGTTCGGTCAGCATTTGTACGGTTTCCACCAGCCACTCACCGCCTTCAATTTTGACCAT".
	"GTTGGCTCCGGCACGCATTACCGTTGCGGCGTTTTCAAAGGCTTGTTCCGGCGTGGCATACGCCATAAAT:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????".
	"??????????=====;96666/.......;=>>>?DDDDDDDDDDDDDDDDDDDDDBBBBDDDDDCCCDD".
	"DDDDDDDDC>>>CDDDDDD>>?DDDDDDDDDDDDDCC???CCCCCCDDDDDDDDDDDDDDDDDDDDDDDD",
	# 280:
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGACCACAAA".
	"CTGGAACGGCACGTTCGGTCAGCATTTGTACGGTTTCCACCAGCCACTCACCGCCTTCAATTTTGACCAT".
	"GTTGGCTCCGGCACGCATTACCGTTGCGGCGTTTTCAAAGGCTTGTTCCGGCGTGGCATACGCCATAAAT".
	"GGCAGGTCAGCCAGCAGAAGGCAGTTTGGCGCGCCGCGACGCACAGCGGCGGTGTGATAGGCGATATCGG:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????".
	"??????????=====;96666/.......;=>>>?DDDDDDDDDDDDDDDDDDDDDBBBBDDDDDCCCDD".
	"DDDDDDDDC>>>CDDDDDD>>?DDDDDDDDDDDDDCC???CCCCCCDDDDDDDDDDDDDDDDDDDDDDDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD",
	# 350:
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGACCACAAA".
	"CTGGAACGGCACGTTCGGTCAGCATTTGTACGGTTTCCACCAGCCACTCACCGCCTTCAATTTTGACCAT".
	"GTTGGCTCCGGCACGCATTACCGTTGCGGCGTTTTCAAAGGCTTGTTCCGGCGTGGCATACGCCATAAAT".
	"GGCAGGTCAGCCAGCAGAAGGCAGTTTGGCGCGCCGCGACGCACAGCGGCGGTGTGATAGGCGATATCGG".
	"CAACGGTAACCGGTAGCGTGGAGTCATGCCCCTGGACCGTCATGCCCAGTGAATCGCCTACCAGCATGAC:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????".
	"??????????=====;96666/.......;=>>>?DDDDDDDDDDDDDDDDDDDDDBBBBDDDDDCCCDD".
	"DDDDDDDDC>>>CDDDDDD>>?DDDDDDDDDDDDDCC???CCCCCCDDDDDDDDDDDDDDDDDDDDDDDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDDDDDDDDDDDDDDDAAADDDDDDD;;;;DDDEFFFDEGGGGDDDDDDDDDCCCDDDDDD",
	""
);

my @args1 = (
	"",
	"-m 1",
	"--strata -m 1"
);

my @args2 = (
	"-v 0",
	"-v 1",
	"-v 2",
	"-v 3",
	"-n 0",
	"-n 1",
	"-n 2",
	"-n 3"
);

sub btrun {
	my ($read, $args, $color) = @_;
	$args .= $color ? " -C" : "";
	my $cmd = "$bowtie $args -c e_coli \"$read\"";
	print "\n$cmd\n\n";
	system($cmd) && die;
	$cmd = "$bowtie_d $args -c e_coli \"$read\"";
	print "\n$cmd\n\n";
	system($cmd) && die;
	print "PASSED: \"$args\"\n";
}

for my $r (@reads) {
	next if $r eq "";
	for my $a1 (@args1) {
		for my $a2 (@args2) {
			btrun($r, "$a1 $a2", 0);
		}
	}
}
