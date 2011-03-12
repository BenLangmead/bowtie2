#!/usr/bin/perl -w

##
# cs_trim.pl
#
# Basic tests to ensure that colorspace primer trimming is working as
# expected.
#

use strict;
use warnings;

my $bowtie = "./bowtie";
if(system("$bowtie --version") != 0) {
	$bowtie = `which bowtie`;
	chomp($bowtie);
	if(system("$bowtie --version") != 0) {
		die "Could not find bowtie in current directory or in PATH\n";
	}
}

my $bowtie_d = "./bowtie-debug";
if(system("$bowtie_d --version") != 0) {
	$bowtie_d = `which bowtie-debug`;
	chomp($bowtie_d);
	if(system("$bowtie_d --version") != 0) {
		die "Could not find bowtie-debug in current directory or in PATH\n";
	}
}

if(! -f "e_coli_c.1.ebwt") {
	print STDERR "Making colorspace e_coli index\n";
	system("make bowtie-build") && die;
	system("bowtie-build -C genomes/NC_008253.fna e_coli_c") && die;
} else {
	print STDERR "Colorspace e_coli index already present...\n";
}

sub readToFastq {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	system("rm -f $fname");
	open TMPFQ, ">$fname" || die "Could not open $fname for writing\n";
	print TMPFQ "\@r\n$r[0]\n+\n$r[1]\n";
	close(TMPFQ);
}

sub readToBFASTFastq {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	system("rm -f $fname");
	open TMPFQ, ">$fname" || die "Could not open $fname for writing\n";
	my $q = $r[1];
	if(substr($r[0], 0, 1) =~ /[ACGT]/ &&
	   substr($r[0], 1, 1) =~ /[0123\.]/)
	{
		$q = substr($q, 1);
	}
	print TMPFQ "\@r\n$r[0]\n+\n$q\n";
	close(TMPFQ);
}

sub readToFasta {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	system("rm -f $fname");
	open TMPFA, ">$fname" || die "Could not open $fname for writing\n";
	print TMPFA ">r\n$r[0]\n";
	close(TMPFA);
}

sub readToRaw {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	system("rm -f $fname");
	open TMPR, ">$fname" || die "Could not open $fname for writing\n";
	print TMPR "$r[0]\n";
	close(TMPR);
}

sub readToTabbed {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	system("rm -f $fname");
	open TMPT, ">$fname" || die "Could not open $fname for writing\n";
	print TMPT "r\t$r[0]\t$r[1]\n";
	close(TMPT);
}

sub readToTabbed2 {
	my ($rstr1, $rstr2, $fname) = @_;
	my @r1 = split(/[:]/, $rstr1);
	my @r2 = split(/[:]/, $rstr2);
	system("rm -f $fname");
	open TMPT, ">$fname" || die "Could not open $fname for writing\n";
	print TMPT "r\t$r1[0]\t$r1[1]\t$r2[0]\t$r2[1]\n";
	close(TMPT);
}

my $evenodd = 0;
sub readToFastaQV {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	$evenodd++;
	system("rm -f $fname $fname.qv");
	open TMPT, ">$fname" || die "Could not open $fname for writing\n";
	print TMPT "# Some\n#\t annoting\n#comments\n>r\n$r[0]\n";
	print TMPT "\n\n>r\nNNNN\n";
	close(TMPT);
	open TMPQ, ">$fname.qv" || die "Could not open $fname.qv for writing\n";
	print TMPQ "#annoying\n;comment\n>r\n$r[1]".((($evenodd % 2) == 0) ? " " : "")."\n";
	print TMPQ ">r\n10 10 10 10\n".((($evenodd % 2) == 0) ? " " : "")."\n";
	close(TMPQ);
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
		$di =~ s/[URYMKWSBDHV]/N/gi;
		defined($cmap{$di}) || die "Bad dinuc: $di\n";
		$ret .= ($nucs ? $nmap{$cmap{$di}} : $cmap{$di});
	}
	return $ret;
}

# Utility function that returns the reverse complement of its argument
sub reverseComp($) {
	my $r = shift;
	$r = reverse($r);
	$r =~ tr/aAcCgGtT/tTgGcCaA/;
	return $r;
}

sub btrun {
	my ($name, $args, $num, $char1, $qual1, $char2, $qual2) = @_;
	my $case = "$name, $args, $num, $char1, $qual1, $char2, $qual2";
	$char1 =~ tr/0123./ACGTN/ if defined($char1);
	$char2 =~ tr/0123./ACGTN/ if defined($char2);
	for my $bt ($bowtie, $bowtie_d) {
		my $cmd = "$bt $args";
		print "$cmd\n";
		open BTIE, "$cmd |" || die;
		while(<BTIE>) {
			my @s = split;
			my ($seq, $quals) = ($s[4], $s[5]);
			my $sl = length($seq);
			my $el = $num;
			$sl == $el || die "Expected seq length $el, got $sl\n";
			my $ql = length($quals);
			$ql == $el || die "Expected qual length $el, got $ql\n";
		}
		if($args =~ /-C/) {
			$cmd .= " --col-cseq --col-cqual";
			print "$cmd\n";
			open BTIE, "$cmd |" || die;
			my $line = 0;
			while(<BTIE>) {
				my @s = split;
				my ($seq, $quals) = ($s[4], $s[5]);
				my $fc = substr($seq, 0, 1);
				my $fq = substr($quals, 0, 1);
				if($line == 0) {
					!defined($char1) || $fc eq $char1 || die "Expected first char on line 1 $char1, got $fc\nCase: $case\nBowtie: $_";
					!defined($qual1) || $fq eq $qual1 || die "Expected first qual on line 1 $qual1, got $fq\nCase: $case\nBowtie: $_";
				} elsif($line == 1) {
					!defined($char2) || $fc eq $char2 || die "Expected first char $char2, got $fc\nCase: $case\nBowtie: $_";
					!defined($qual2) || $fq eq $qual2 || die "Expected first qual $qual2, got $fq\nCase: $case\nBowtie: $_";
				}
				$line++;
			}
			close(BTIE);
		}
		$? == 0 || die "$bt returned non-zero status $?\n";
	}
	print "PASSED $name\n";
}

my $m1 = colorize("CTGACTGCAACGGGCAATATGTCTCTGTGTGGA", 0);
my $m2 = colorize("CATCACCATTACCACAGGTAACGGTGCGGGCTG", 0);

my $m1n = colorize("CTGACTGCAACGGGCAATATGTCTCTGTGTGGA", 1);
my $m2n = colorize("CATCACCATTACCACAGGTAACGGTGCGGGCTG", 1);

my $n1 = "CTGACTGCAACGGGCAATATGTCTCTGTGTGGA";
my $n2 = reverseComp("CATCACCATTACCACAGGTAACGGTGCGGGCTG");

my $q = "ABCDEFGHIJKLMNOPQRSTUVWXYZZZZZZZZZ";

sub intize {
	my $read = shift;
	my $start = shift;
	$start = 0 unless defined($start);
	my @rs = split(/:/, $read);
	my $q = $rs[1];
	my $qint = "";
	for(my $i = $start; $i < length($q); $i++) {
		my $qc = substr($q, $i, 1);
		$qint .= " " if $qint ne "";
		$qint .= (ord($qc)-33);
	}
	return "$rs[0]:$qint";
}

my @trimReads = (
	# Trim me
	"T0$m1:$q",
	# Trim me
	"A1$m2:$q",
	# Trim me
	"G.$m1:$q",
	# Trim me
	"C.$m2:$q"
);

my @nucReads = (
	# Don't trim me
	"$n1:".(substr($q, 2)),
	# Don't trim me
	"$n2:".(substr($q, 2))
);

my @untrimReads = (
	# Don't trim me
	"TT$m1:$q",
	# Don't trim me
	"CG$m2:$q",
	# Don't trim me
	"31$m1:$q",
	# Don't trim me
	"20$m2:$q",
	# Don't trim me
	"AC$m1n:$q",
	# Don't trim me
	"CC$m2n:$q"
);

sub ca {
	my ($str, $off) = @_;
	return substr($str, $off, 1);
}

for(my $i = 0; $i < scalar(@trimReads); $i += 2) {
	
	my $r1 = $trimReads[$i];
	my $r2 = $trimReads[$i+1];

	btrun("trim/-c/paired",   "-C -c e_coli_c -1 $r1 -2 $r2",
	      length($m1)-1,
	      ca($m1, 0), ca($q, 2),
	      ca($m2, 0), ca($q, 2));
	btrun("trim/-c/unpaired", "-C -c e_coli_c $r1,$r2",
	      length($m1)-1,
	      ca($m1, 0), ca($q, 2),
	      ca($m2, 0), ca($q, 2));

	for(my $j = 0; $j < 4; $j++) {
		my ($fq1, $fq2) = (".tmp1.fq", ".tmp2.fq");
		my $args = "";
		my $peargs = "";
		my $name = "";
		if($j == 0) {
			readToFastq($r1, $fq1);
			readToFastq($r2, $fq2);
			$peargs = "-q -1 $fq1 -2 $fq2";
			$args = "-q $fq1,$fq2";
			$name = "fastq";
		} elsif($j == 1) {
			$fq1 = ".tmp1.bfast.fq";
			$fq2 = ".tmp2.bfast.fq";
			readToBFASTFastq($r1, $fq1);
			readToBFASTFastq($r2, $fq2);
			$peargs = "-q -1 $fq1 -2 $fq2";
			$args = "-q $fq1,$fq2";
			$name = "fastq-bfast";
		} elsif($j == 2) {
			$fq1 = ".tmp1.fa";
			$fq2 = ".tmp2.fa";
			readToFastaQV(intize($r1), $fq1);
			readToFastaQV(intize($r2), $fq2);
			$peargs = "-f -1 $fq1 -2 $fq2 --Q1 $fq1.qv --Q2 $fq2.qv";
			$args = "-f $fq1,$fq2 -Q $fq1.qv,$fq2.qv";
			$name = "fastq-qv";
		} else {
			$fq1 = ".tmp1.fa";
			$fq2 = ".tmp2.fa";
			readToFastaQV(intize($r1, 1), $fq1);
			readToFastaQV(intize($r2, 1), $fq2);
			$peargs = "-f -1 $fq1 -2 $fq2 --Q1 $fq1.qv --Q2 $fq2.qv";
			$args = "-f $fq1,$fq2 -Q $fq1.qv,$fq2.qv";
			$name = "fastq-qv";
		}
	
		btrun("no-trim/$name/paired", "-C e_coli_c $peargs",
		      length($m1)-1,
		      ca($m1, 0), ca($q, 2),
		      ca($m2, 0), ca($q, 2));
		btrun("no-trim/$name/unpaired", "-C e_coli_c $args",
		      length($m1)-1,
		      ca($m1, 0), ca($q, 2),
		      ca($m2, 0), ca($q, 2));
		btrun("trim5/$name/unpaired", "-5 3 -C e_coli_c $args",
		      length($m1)-1-3,
		      ca($m1, 3), ca($q, 5),
		      ca($m2, 3), ca($q, 5));
		btrun("trim35/$name/unpaired", "-5 5 -3 3 -C e_coli_c $args",
		      length($m1)-1-8,
		      ca($m1, 5), ca($q, 7),
		      ca($m2, 5), ca($q, 7));
	}
	
	readToFasta($r1, ".tmp1.fa");
	readToFasta($r2, ".tmp2.fa");
	btrun("no-trim/fasta/paired",   "-C e_coli_c -f -1 .tmp1.fa -2 .tmp2.fa",
	      length($m1)-1,
	      ca($m1, 0), "I",
	      ca($m2, 0), "I");
	btrun("no-trim/fasta/unpaired", "-C e_coli_c -f .tmp1.fa,.tmp2.fa",
	      length($m1)-1,
	      ca($m1, 0), "I",
	      ca($m2, 0), "I");
	btrun("trim5/fasta/unpaired", "-5 3 -C e_coli_c -f .tmp1.fa,.tmp2.fa",
	      length($m1)-1-3,
	      ca($m1, 3), "I",
	      ca($m2, 3), "I");
	btrun("trim35/fasta/unpaired", "-5 5 -3 3 -C e_coli_c -f .tmp1.fa,.tmp2.fa",
	      length($m1)-1-8,
	      ca($m1, 5), "I",
	      ca($m2, 5), "I");
	
	readToFasta($r1, ".tmp1.fa");
	readToFasta($r2, ".tmp2.fa");
	
	readToRaw($r1, ".tmp1.raw");
	readToRaw($r2, ".tmp2.raw");
	btrun("no-trim/raw/paired",   "-C e_coli_c -r -1 .tmp1.raw -2 .tmp2.raw",
	      length($m1)-1,
	      ca($m1, 0), "I",
	      ca($m2, 0), "I");
	btrun("no-trim/raw/unpaired", "-C e_coli_c -r .tmp1.raw,.tmp2.raw",
	      length($m1)-1,
	      ca($m1, 0), "I",
	      ca($m2, 0), "I");
	btrun("trim5/raw/unpaired", "-5 3 -C e_coli_c -r .tmp1.raw,.tmp2.raw",
	      length($m1)-1-3,
	      ca($m1, 3), "I",
	      ca($m2, 3), "I");
	btrun("trim35/raw/unpaired", "-5 5 -3 3 -C e_coli_c -r .tmp1.raw,.tmp2.raw",
	      length($m1)-1-8,
	      ca($m1, 5), "I",
	      ca($m2, 5), "I");
}

for(my $i = 0; $i < scalar(@nucReads); $i += 2) {
	my $r1 = $nucReads[$i];
	my $r2 = $nucReads[$i+1];
	btrun("trim/-c/paired",   "-c e_coli -1 $r1 -2 $r2",
	      length($n1), undef, undef, undef, undef);
	btrun("trim/-c/unpaired", "-c e_coli $r1,$r2",
	      length($n1), undef, undef, undef, undef);
}

for(my $i = 0; $i < scalar(@untrimReads); $i += 2) {
	my $r1 = $untrimReads[$i];
	my $r2 = $untrimReads[$i+1];
	my $colri = index($r1, ":");
	my $colrp = index($r2, ":");
	btrun("no-trim/-c/paired",   "-C -c e_coli_c -1 $r1 -2 $r2",
	      length($m1)+1,
	      ca($r1, 0), ca($r1, $colri+1+0),
	      ca($r2, 0), ca($r2, $colrp+1+0));
	btrun("no-trim/-c/unpaired", "-C -c e_coli_c $r1,$r2",
	      length($m1)+1,
	      ca($r1, 0), ca($r1, $colri+1+0),
	      ca($r2, 0), ca($r2, $colri+1+0));
	btrun("trim5/-c/unpaired", "-5 3 -C -c e_coli_c $r1,$r2",
	      length($m1)+1-3,
	      ca($r1, 3), ca($r1, $colri+1+3),
	      ca($r2, 3), ca($r2, $colri+1+3));
	btrun("trim35/-c/unpaired", "-5 5 -3 3 -C -c e_coli_c $r1,$r2",
	      length($m1)+1-8,
	      ca($r1, 5), ca($r1, $colri+1+5),
	      ca($r2, 5), ca($r2, $colri+1+5));
	
	readToFastq($r1, ".tmp1.fq");
	readToFastq($r2, ".tmp2.fq");
	btrun("no-trim/fastq/paired", "-C e_coli_c -q -1 .tmp1.fq -2 .tmp2.fq",
	      length($m1)+1,
	      ca($r1, 0), ca($r1, $colri+1+0),
	      ca($r2, 0), ca($r2, $colri+1+0));
	btrun("no-trim/fastq/unpaired", "-C e_coli_c -q .tmp1.fq,.tmp2.fq",
	      length($m1)+1,
	      ca($r1, 0), ca($r1, $colri+1+0),
	      ca($r2, 0), ca($r2, $colri+1+0));
	btrun("trim5/fastq/unpaired", "-5 3 -C e_coli_c -q .tmp1.fq,.tmp2.fq",
	      length($m1)+1-3,
	      ca($r1, 3), ca($r1, $colri+1+3),
	      ca($r2, 3), ca($r2, $colri+1+3));
	btrun("trim35/fastq/unpaired", "-5 5 -3 3 -C e_coli_c -q .tmp1.fq,.tmp2.fq",
	      length($m1)+1-8,
	      ca($r1, 5), ca($r1, $colri+1+5),
	      ca($r2, 5), ca($r2, $colri+1+5));

	readToFasta($r1, ".tmp1.fa");
	readToFasta($r2, ".tmp2.fa");
	btrun("no-trim/fasta/paired",   "-C e_coli_c -f -1 .tmp1.fa -2 .tmp2.fa",
	      length($m1)+1,
	      ca($r1, 0), "I",
	      ca($r2, 0), "I");
	btrun("no-trim/fasta/unpaired", "-C e_coli_c -f .tmp1.fa,.tmp2.fa",
	      length($m1)+1,
	      ca($r1, 0), "I",
	      ca($r2, 0), "I");
	btrun("trim5/fasta/unpaired", "-5 3 -C e_coli_c -f .tmp1.fa,.tmp2.fa",
	      length($m1)+1-3,
	      ca($r1, 3), "I",
	      ca($r2, 3), "I");
	btrun("trim35/fasta/unpaired", "-5 5 -3 3 -C e_coli_c -f .tmp1.fa,.tmp2.fa",
	      length($m1)+1-8,
	      ca($r1, 5), "I",
	      ca($r2, 5), "I");

	readToRaw($r1, ".tmp1.raw");
	readToRaw($r2, ".tmp2.raw");
	btrun("no-trim/raw/paired",   "-C e_coli_c -r -1 .tmp1.raw -2 .tmp2.raw",
	      length($m1)+1,
	      ca($r1, 0), "I",
	      ca($r2, 0), "I");
	btrun("no-trim/raw/unpaired", "-C e_coli_c -r .tmp1.raw,.tmp2.raw",
	      length($m1)+1,
	      ca($r1, 0), "I",
	      ca($r2, 0), "I");
	btrun("trim5/raw/unpaired", "-5 3 -C e_coli_c -r .tmp1.raw,.tmp2.raw",
	      length($m1)+1-3,
	      ca($r1, 3), "I",
	      ca($r2, 3), "I");
	btrun("trim35/raw/unpaired", "-5 5 -3 3 -C e_coli_c -r .tmp1.raw,.tmp2.raw",
	      length($m1)+1-8,
	      ca($r1, 5), "I",
	      ca($r2, 5), "I");

	readToTabbed($r1, ".tmp1.tabbed");
	readToTabbed($r2, ".tmp2.tabbed");
	readToTabbed2($r1, $r2, ".tmp.tabbed");
	btrun("no-trim/tabbed/paired",   "-C e_coli_c --12 .tmp.tabbed",
	      length($m1)+1,
	      ca($r1, 0), ca($r1, $colri+1+0),
	      ca($r2, 0), ca($r2, $colri+1+0));
	btrun("no-trim/tabbed/unpaired", "-C e_coli_c --12 .tmp1.tabbed,.tmp2.tabbed",
	      length($m1)+1,
	      ca($r1, 0), ca($r1, $colri+1+0),
	      ca($r2, 0), ca($r2, $colri+1+0));
	btrun("trim5/tabbed/unpaired", "-5 3 -C e_coli_c --12 .tmp1.tabbed,.tmp2.tabbed",
	      length($m1)+1-3,
	      ca($r1, 3), ca($r1, $colri+1+3),
	      ca($r2, 3), ca($r2, $colri+1+3));
	btrun("trim35/tabbed/unpaired", "-5 5 -3 3 -C e_coli_c --12 .tmp1.tabbed,.tmp2.tabbed",
	      length($m1)+1-8,
	      ca($r1, 5), ca($r1, $colri+1+5),
	      ca($r2, 5), ca($r2, $colri+1+5));
}
