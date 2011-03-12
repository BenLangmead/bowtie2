#!/usr/bin/perl -w

##
# Give simple tests with known results to bowtie2.
#

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin); 
use lib $Bin;
use List::Util qw(max min);
use Data::Dumper;
use DNA;
use Clone qw(clone);

my $bowtie2 = "";
my $bowtie2_build = "";

GetOptions(
	"bowtie2=s"       => \$bowtie2,
	"bowtie2-build=s" => \$bowtie2_build) || die "Bad options";

if(! -x $bowtie2 || ! -x $bowtie2_build) {
	my $bowtie2_dir = `dirname $bowtie2`;
	my $bowtie2_exe = `basename $bowtie2`;
	my $bowtie2_build_exe = `basename $bowtie2_build`;
	chomp($bowtie2_dir);
	chomp($bowtie2_exe);
	chomp($bowtie2_build_exe);
	system("make -C $bowtie2_dir $bowtie2_exe $bowtie2_build_exe") && die;
}

(-x $bowtie2)       || die "Cannot run '$bowtie2'";
(-x $bowtie2_build) || die "Cannot run '$bowtie2_build'";

#
# Penalties for mismatches
#   MMP={Cxx|Q|RQ}
#     Cxx = constant, where constant is integer xx
#     Q   = equal to quality
#     R   = equal to maq-rounded quality value (rounded to nearest
#           10, can't be greater than 30)
#
# Penalties for mismatches where read char=N
#   NP={Cxx|Q|RQ}
#     Cxx = constant, where constant is integer xx
#     Q   = equal to quality
#     R   = equal to maq-rounded quality value (rounded to nearest
#           10, can't be greater than 30)
#
# Penalties for read gaps
#   RDG=xx,yy,zz
#     xx = read gap open penalty
#     yy = read gap extension penalty constant coefficient
#          (defaults to open penalty)
#     zz = read gap extension penalty linear coefficient
#          (defaults to 0)
#
# Penalties for reference gaps
#   RFG=xx,yy,zz
#     xx = ref gap open penalty
#     yy = ref gap extension penalty constant coefficient
#          (defaults to open penalty)
#     zz = ref gap extension penalty linear coefficient
#          (defaults to 0)
#
# Per-read penalty ceiling as a function of read length
#   CEIL=xx,yy
#     xx = cost ceiling constant coefficient
#     yy = cost ceiling linear coefficient (set to 0 if
#          unspecified)
#
# Per-read N ceiling as a function of read length
#   NCEIL=xx,yy
#     xx = N ceiling constant coefficient
#     yy = N ceiling linear coefficient (set to 0 if unspecified)
#
# Multiseed parameters
#   SEED=xx,yy,zz
#     xx = mismatches allowed in seed family; can be 0, 1, 2
#     yy = length of seed family (set to 28 if unspecified)
#     zz = period between seed family placements (seed to 14 if
#          unspecified)
#

my @cases = (

	#
	# Alignment involving ambiguous reference character
	#

	# First read has non-compatible unambiguous charcacter (G for Y),
	# second read has compatible one
	{ ref    => [ "TTGTTYGT" ],
	  reads  => [ "TTGTTGGT", "TTGTTCGT" ],
	  args   => "",
	  report => "-a -P \"SEED=0,5,1;NCEIL=2,0\"",
	  hits   => [ { 0 => 1 }, { 0 => 1 } ],
	  norc   => 1,
	  edits  => [ "5:Y>G", "5:Y>C" ] },

	#
	# Alignment with multi-character read gap
	#

	# Reads 1 and 2 don't have overhang, reads 3 and 4 overhang opposite ends
	{ ref    => [ "ATATGCCCCATGCCCCCCTCCG" ],
	  reads  => [ "ATATGCCCCCCCCCCTCCG" ],
	  #                     ^
	  #                     9:ATG>- 
	  args   => "",
	  report => "-a --overhang -P \"RDG=5,5;SEED=0,8,1\"",
	  hits   => [ { 0 => 1 } ],
	  edits  => [ "9:ATG>-" ],
	  norc   => 1 },

	# Reads 1 and 2 don't have overhang, reads 3 and 4 overhang opposite ends
	{ ref    => [ "ATATGCCCCATGCCCCCCTCCG" ],
	  reads  => [ "CGGAGGGGGGGGGGCATAT" ],
	  #            ATATGCCCCCCCCCCTCCG
	  #                     ^         
	  #                     10:GTA>- 
	  args   => "",
	  report => "-a --overhang -P \"RDG=5,5;SEED=0,8,1\"",
	  hits   => [ { 0 => 1 } ],
	  edits  => [ "10:GTA>-" ],
	  norc   => 1 },

	#
	# Alignment with overhang
	#

	# Reads 1 and 2 don't have overhang, reads 3 and 4 overhang opposite ends
	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TGTTCGT", "TTGTTCG", "GTTCGTA", "ATTGTTC" ],
	  args   => "",
	  report => "-a --overhang -P \"SEED=0,2,1;NCEIL=2,0\"",
	  hits   => [ { 1 => 1 }, { 0 => 1 }, { 2 => 1 }, { -1 => 1 } ] },

	# Same as previous case but --overhang not specified
	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TGTTCGT", "TTGTTCG", "GTTCGTA", "ATTGTTC" ],
	  args   => "",
	  report => "-a -P \"SEED=0,2,1;NCEIL=2,0\"",
	  hits   => [ { 1 => 1 }, { 0 => 1 } ] }, # only the internal hits

	#
	# Colorspace alignment
	#

	# 1 SNP and 1 N
	{ ref      => [ "TTGTTC" ],
	  #         Ref: TTGTTCn
	  #               TGTACA
	  #                  ^snp
	  reads    => [ "11311" ],
	  dec_seq  => [ "TGTACA" ],
	  dec_qual => [ "IqqqqI" ],
	  args     =>   "--col-keepends --norc --overhang -P \"NCEIL=1,0;CEIL=10,0;SEED=0,2,1;MMP=C5;SNP=9\"",
	  report   =>   "-a",
	  hits     => [ { 1 => 1 } ],
	  color    => 1 },
	
	{ ref      => [ "TTGTTC" ],
	  reads    => [ "01102"  ],
	  dec_seq  => [ "TTGTTC" ],
	  dec_qual => [ "IqqqqI" ],
	  args     =>   "--col-keepends",
	  report   =>   "-a",
	  hits     => [ { 0 => 1 } ],
	  color    => 1 },

	{ ref      => [ "TTGTTC" ],
	  reads    => [ "01102"  ],
	  dec_seq  => [ "TTGTTC" ],
	  dec_qual => [ "IqqqqI" ],
	  args     =>   "--col-keepends -P \"SEED=0,3,1\"",
	  report   =>   "-a",
	  hits     => [ { 0 => 1 } ],
	  color    => 1 },

	{ ref      => [ "TTGTTC" ],
	  reads    => [ "01102" ],
	  dec_seq  => [ "TGTT" ],
	  dec_qual => [ "qqqq" ],
	  args     =>   "-P \"SEED=0,3,1\"",
	  report   =>   "-a",
	  hits     => [ { 1 => 1 } ],
	  color    => 1 },

	{ ref      => [ "TTGTTC" ],
	  reads    => [ "01112" ],
	  #                 ^mm
	  args     =>   "--col-keepends -P \"SEED=0,4,1;MMP=C5\"",
	  report   =>   "-a",
	  hits     => [ { } ],
	  color    => 1 },

	{ ref      => [ "TTGTTC" ],
	  reads    => [ "01112" ],
	  #                 ^mm
	  dec_seq  => [ "TTGTTC" ],
	  dec_qual => [ "Iqq!!I" ],
	  args     =>   "--col-keepends -P \"SEED=0,3,1;MMP=C5\"",
	  report   =>   "-a",
	  hits     => [ { 0 => 1 } ],
	  color    => 1 },

	# Two color mismatches carry a smaller penalty than 1 SNP
	{ ref      => [ "TTGTTC" ],
	  #              TTGTAC
	  #                  ^snp
	  reads    => [ "01131" ],
	  dec_seq  => [ "TTGTTC" ],
	  dec_qual => [ "Iqq!!!" ],
	  args     =>   "--col-keepends -P \"SEED=0,3,1;MMP=C5;SNP=11\"",
	  report   =>   "-a",
	  hits     => [ { 0 => 1 } ],
	  color    => 1 },

	# 1 SNP carries a smaller penalty than 2 color mismatches
	{ ref      => [ "TTGTTC" ],
	  #              TTGTAC
	  #                  ^snp
	  reads    => [ "01131" ],
	  dec_seq  => [ "TTGTAC" ],
	  dec_qual => [ "IqqqqI" ],
	  args     =>   "--col-keepends -P \"SEED=0,3,1;MMP=C5;SNP=9\"",
	  report   =>   "-a",
	  hits     => [ { 0 => 1 } ],
	  color    => 1 },

	# Check that pseudo-random generation is always the same for
	# same-sequence, same-name reads

	{ ref    => [ "AAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCC" ],
	  reads  => [ "AA", "AA", "AA", "AA", "CC", "CC", "CC", "CC", "AA", "AA", "AA", "AA", "CC", "CC", "CC", "CC" ],
	  names  => [ "r1", "r1", "r1", "r1", "r2", "r2", "r2", "r2", "r3", "r3", "r3", "r3", "r4", "r4", "r4", "r4" ],
	  args   => "",
	  check_random => 1,
	  report => "-k 1" },

	# Check basic -m funtionality

	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TTGTTCGT", "TTGTTCGT", "TTGTTCGT", "TTGTTCGT", "TTGTTCGT" ],
	  args   => "",
	  report => "-m 1 -a",
	  hits   => [ { 0 => 1 }, { 0 => 1 }, { 0 => 1 }, { 0 => 1 }, { 0 => 1 } ],
	  edits  => [ "-",        "-",        "-",        "-",        "-",       ] },

	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TTGTTCGT", "TAAAACGT", "TTGTTCGT", "TAATTCGT",    "TTGTTCGT" ],
	  args   => "",
	  report => "-m 1 -a -P \"SEED=0,3,1;MMP=C9\"",
	  hits   => [ { 0 => 1 }, {        }, { 0 => 1 }, { 0 => 1 },    { 0 => 1 } ],
	  edits  => [ "-",        undef,      "-",        "1:T>A,2:G>A", "-",       ],
	  norc   => 1 },

	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => "",
	  report => "-m 3 -a",
	  hits   => [ { 0 => 1, 8 => 1 } ] },

	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => "",
	  report => "-m 2 -a",
	  hits   => [ { 0 => 1, 8 => 1 } ] },

	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => "",
	  report => "-m 1 -a",
	  hits   => [ { } ] },

	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => "",
	  report => "-M 1",
	  hits   => [ { 0 => 1, 8 => 1 } ],
	  hits_are_superset => [ 1 ] },

	# Mess with arguments
	
	# Default should be 1-mismatch, so this shouldn't align
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTATTAGT" ],
	  args   => "",
	  report => "-a",
	  hits   => [ {  } ] },

	# Shouldn't align with 0 mismatches either
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTATTAGT" ],
	  args   => "-P SEED=0",
	  report => "-a",
	  hits   => [ { } ] },

	# Should align with 2 mismatches, provided the mismatches are in the right spots
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TAGTTCAT" ],
	  args   => "-P \"SEED=2;MMP=C1\"",
	  report => "-a",
	  hits   => [ { 0 => 1, 8 => 1 } ] },

	# Should align with 2 mismatches, provided the mismatches are in the right spots
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTATTAGT" ],
	  args   =>   "-P \"SEED=2;MMP=C1\"",
	  report =>   "-a",
	  hits   => [ { } ] },

	# Should align with 0 mismatches if we can wedge a seed into the 2
	# matching characters between the two mismatches.  Here we fail to
	# wedge a length-3 seed in (there's no room)
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTATTAGT" ],
	  args   => "-P \"SEED=0,3,1;MMP=C1\"",
	  report => "-a",
	  hits   => [ { } ] },

	# Should align with 0 mismatches if we can wedge a seed into the 2
	# matching characters between the two mismatches.  Here we wedge a
	# length-2 seed in
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTATTAGT" ],
	  args   => "-P \"SEED=0,2,1;MMP=C1\"",
	  report => "-a",
	  hits   => [ { 0 => 1, 3 => 1, 4 => 1, 5 => 1, 7 => 1, 8 => 1} ] },

	# Following cases depend on this being the case:
	#
	# static const float DEFAULT_CEIL_CONST = 3.0f;
	# static const float DEFAULT_CEIL_LINEAR = 3.0f;

	# Just enough budget for hits, so it should align
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCAT" ], # budget = 3 + 8 * 3 = 27
	  args   => "-P \"SEED=2;MMP=C27\"", # penalty = 27
	  report => "-a",
	  hits   => [ { 0 => 1, 8 => 1 } ] },

	# Not quite enough budget for hits, so it should NOT align
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCAT" ], # budget = 3 + 8 * 3 = 27
	  args   =>   "-P \"SEED=2;MMP=C28\"", # penalty = 28
	  report =>   "-a",
	  hits   => [ { } ] },

	# Alignment with 1 read gap
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCTTTGTT" ], # budget = 3 + 12 * 3 = 39
	  args   =>   "-P \"SEED=0,3,1;RDG=39\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ] },

	# Alignment with 1 read gap, but not enough budget
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCTTTGTT" ], # budget = 3 + 12 * 3 = 39
	  args   =>   "-P \"SEED=0,3,1;RDG=40\"",
	  report =>   "-a",
	  hits   => [ { } ] },

	# Alignment with 1 reference gap
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGATTTGTT" ], # budget = 3 + 14 * 3 = 45
	  args   =>   "-P \"SEED=0,3,1;RFG=45\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ] },

	# Alignment with 1 reference gap, but not enough budget
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGATTTGTT" ], # budget = 3 + 14 * 3 = 45
	  args   =>   "-P \"SEED=0,3,1;RFG=46\"",
	  report =>   "-a",
	  hits   => [ { } ] },

	# Alignment with 1 reference gap and 1 read gap
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTGTTTGATTCGT" ], # budget = 3 + 16 * 3 = 51
	  args   =>   "-P \"SEED=0,3,1;RFG=25;RDG=26\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ] },

	# Alignment with 1 reference gap and 1 read gap, but not enough budget
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTGTTTGATTCGT" ], # budget = 3 + 16 * 3 = 51
	  args   =>   "-P \"SEED=0,3,1;RFG=26;RDG=26\"",
	  report =>   "-a",
	  hits   => [ { } ] },
);

##
# Take a list of reference sequences and write them to a temporary
# FASTA file of the given name.
#
sub writeFasta($$) {
	my ($l, $fa) = @_;
	open(FA, ">$fa") || die "Could not open $fa for writing";
	my $idx = 0;
	for(my $i = 0; $i < scalar(@$l); $i++) {
		print FA ">$idx\n".$l->[$i]."\n";
		$idx++;
	}
	close(FA);
}

##
# Take a lists of named reads/mates and write them to appropriate
# files.
#
sub writeReads($$$$$$) {
	my ($reads, $mate1s, $mate2s, $names, $fq1, $fq2) = @_;
	open(FQ1, ">$fq1") || die "Could not open '$fq1' for writing";
	open(FQ2, ">$fq2") || die "Could not open '$fq2' for writing";
	my $pe = (defined($mate1s) && $mate1s ne "");
	if($pe) {
		for(my $i = 0; $i < scalar(@$mate1s); $i++) {
			my $m1 = $mate1s->[$i];
			my $m2 = $mate2s->[$i];
			my $nm = $names->[$i];
			print FQ1 "\@$nm/1\n$m1\n+\n".("I" x length($m1))."\n";
			print FQ2 "\@$nm/2\n$m2\n+\n".("I" x length($m2))."\n";
		}
	} else {
		for(my $i = 0; $i < scalar(@$reads); $i++) {
			my $read = $reads->[$i];
			my $nm = $names->[$i];
			print FQ1 "\@$nm\n$read\n+\n".("I" x length($read))."\n";
		}
	}
	close(FQ1);
	close(FQ2);
}

##
# Run bowtie2 with given arguments
#
sub runbowtie2($$$$$$$$$$) {
	my (
		$args,
		$color,
		$fa,
		$reportargs,
		$reads,
		$mate1s,
		$mate2s,
		$names,
		$ls,
		$rawls) = @_;
	$args .= " --quiet";
	$reportargs = $reportargs || "-a";
	$args .= " -C" if $color;
	$args .= " $reportargs";
	# Write the reference to a fasta file
	my $build_args = ($color ? "-C" : "");
	print "References:\n";
	open(FA, $fa) || die;
	while(<FA>) { print $_; }
	close(FA);
	my $cmd = "$bowtie2_build --quiet $build_args $fa .simple_tests.tmp";
	print "$cmd\n";
	system($cmd);
	($? == 0) || die "Bad exitlevel from bowtie2-build: $?";
	my $pe = (defined($mate1s) && $mate1s ne "");
	my $mate1arg;
	my $mate2arg;
	my $readarg;
	my $formatarg = "-c";
	my $readstr = join(",", @$reads);
	if(defined($names)) {
		writeReads($reads, $mate1s, $mate2s, $names, ".simple_tests.1.fq", ".simple_tests.2.fq");
		$mate1arg = ".simple_tests.1.fq";
		$mate2arg = ".simple_tests.2.fq";
		$formatarg = "-q";
		$readarg = $mate1arg;
	} else {
		$mate1arg = $mate1s;
		$mate2arg = $mate2s;
		$readarg = $readstr;
	}
	if($pe) {
		# Paired-end case
		$cmd = "$bowtie2 $args .simple_tests.tmp $formatarg -1 $mate1arg -2 $mate2arg";
		print "$cmd\n";
		open(BT, "$cmd |") || die "Could not open pipe '$cmd |'";
		while(<BT>) {
			my $m1 = $_;
			my $m2 = <BT>;
			defined($m2) || die;
			print $m1;
			print $m2;
			chomp($m1);
			chomp($m2);
			push @$ls,    [ split(/\t/, $m1, -1) ];
			push @$ls,    [ split(/\t/, $m2, -1) ];
			push @$rawls, $m1;
			push @$rawls, $m2;
		}
		close(BT);
	} else {
		# Unpaired case
		$cmd = "$bowtie2 $args .simple_tests.tmp $formatarg $readarg";
		print "$cmd\n";
		open(BT, "$cmd |") || die "Could not open pipe '$cmd |'";
		while(<BT>) {
			print $_;
			chomp;
			push @$ls,    [ split(/\t/, $_, -1) ];
			push @$rawls, $_;
		}
		close(BT);
	}
	($? == 0) || die "bowtie2 exited with level $?\n";
}

my $tmpfafn = ".simple_tests.pl.fa";
for (my $ci = 0; $ci < scalar(@cases); $ci++) {
	my $c = $cases[$ci];
	writeFasta($c->{ref}, $tmpfafn);
	# For each set of arguments...
	my $a = $c->{args};
	for(my $fwi = 0; $fwi <= 1; $fwi++) {
		my $fw = ($fwi == 0);
		next if !$fw && $c->{norc};
		# Run bowtie2
		my @lines = ();
		my @rawlines = ();
		print $c->{name}."\n" if defined($c->{name});
		my $color = 0;
		$color = $c->{color} if defined($c->{color});
		my $reads = $c->{reads};
		if(!$fw) {
			# Reverse-complement the reads
			my @s = @$reads;
			for(my $i = 0; $i < scalar(@s); $i++) {
				if($color) {
					$s[$i] = reverse $s[$i];
				} else {
					$s[$i] = DNA::revcomp($s[$i]);
				}
			}
			$reads = \@s;
		}
		runbowtie2(
			"$a",
			$color,
			$tmpfafn,
			$c->{report},
			$reads,
			$c->{mate1s},
			$c->{mate2s},
			$c->{names},
			\@lines,
			\@rawlines);
		my $pe = defined($c->{mate1s}) && $c->{mate1s} ne "";
		my ($lastchr, $lastoff) = ("", -1);
		# Keep temporary copies of hits and pairhits so that we can
		# restore for the next orientation
		my $hitstmp = [];
		$hitstmp = clone($c->{hits}) if defined($c->{hits});
		my $pairhitstmp = [];
		$pairhitstmp = clone($c->{pairhits}) if defined($c->{pairhits});
		# Record map from already-seen read name, read sequence and
		# quality to the place on the reference where it's reported.
		# This allows us to check that the pseudo-random generator
		# isn't mistakenly yielding different alignments for identical
		# reads.
		my %seenNameSeqQual = ();
		for my $li (0 .. scalar(@lines)-1) {
			my $l = $lines[$li];
			scalar(@$l) == 8 || die "Bad number of fields; expected 8 got ".scalar(@$l).":\n$rawlines[$li]\n";
			my ($readname, $orient, $chr, $off, $seq, $qual) = @$l;
			if($c->{check_random}) {
				my $rsqKey = "$readname\t$orient\t$seq\t$qual";
				my $rsqVal = "$chr\t$off";
				if(defined($seenNameSeqQual{$rsqKey})) {
					$seenNameSeqQual{$rsqKey} eq $rsqVal ||
						die "Two hits for read/seq/qual:\n$rsqKey\n".
							"had different alignments:\n".
							"$seenNameSeqQual{$rsqKey}\n$rsqVal\n";
				}
				$seenNameSeqQual{$rsqKey} = $rsqVal;
			}
			my $rdi = $readname;
			if($rdi !~ /^[0-9]+$/) {
				# Read name has non-numeric characters.  Figure out
				# what number it is by scanning the names list.
				defined($c->{names}) || die "Non-numeric read name for case with no names specified";
				my $found = 0;
				for(my $i = 0; $i < scalar(@{$c->{names}}); $i++) {
					if($c->{names}->[$i] eq $readname) {
						$rdi = $i;
						$found = 1;
						last;
					}
				} 
				$found || die "No specified name matched reported name $readname";
			}
			# Check that the sequence printed in the alignment is sane
			if($color) {
				# It's a decoded nucleotide sequence
				my $dseq = $c->{dec_seq}->[$rdi];
				if(defined($dseq)) {
					$seq eq $dseq || die "Expected decoded sequence '$seq' from alignment to match '$dseq'";
				}
				my $dqual = $c->{dec_qual}->[$rdi];
				if(defined($dqual)) {
					$qual eq $dqual || die "Expected decoded qualities '$qual' from alignment to match '$dqual'";
				}
			} else {
				
			}
			# Make simply-named copies of some portions of the test case
			# 'hits'
			my %hits = ();
			%hits = %{$c->{hits}->[$rdi]} if defined($c->{hits}->[$rdi]);
			# 'pairhits'
			my %pairhits = ();
			%pairhits = %{$c->{pairhits}->[$rdi]} if defined($c->{pairhits}->[$rdi]);
			# 'hits_are_superset'
			my $hits_are_superset = 0;
			$hits_are_superset = $c->{hits_are_superset}->[$rdi] if defined($ci);
			# edits
			my $edits = undef;
			$edits = $c->{edits}->[$rdi] if defined($c->{edits}->[$rdi]);
			next if $orient eq '*'; # read not aligned
			if($pe && $lastchr ne "") {
				my $offkey = min($lastoff, $off).",".max($lastoff, $off);
				if(defined($c->{pairhits}->[$rdi])) {
					defined($pairhits{$offkey}) ||
						die "No such paired off as $offkey in pairhits list: ".Dumper(\%pairhits)."\n";
					$c->{pairhits}->[$rdi]->{$off}--;
					delete $c->{pairhits}->[$rdi]->{$off} if $c->{pairhits}->[$rdi]->{$off} == 0;
					%pairhits = %{$c->{pairhits}->[$rdi]};
				}
				($lastchr, $lastoff) = ("", -1);
			} elsif($pe) {
				($lastchr, $lastoff) = ($chr, $off);
			} else {
				if(defined($c->{hits}->[$rdi])) {
					defined($hits{$off}) ||
						die "No such off as $off in hits list: ".Dumper(\%hits)."\n";
					$c->{hits}->[$rdi]->{$off}--;
					delete $c->{hits}->[$rdi]->{$off} if $c->{hits}->[$rdi]->{$off} == 0;
					%hits = %{$c->{hits}->[$rdi]};
				}
			}
			if(defined($edits)) {
				my $eds = $l->[-1];
				$eds eq $edits || die "For edit string, expected \"$edits\" got \"$eds\"\n";
			}
		}
		# Go through all the per-read 
		my $klim = scalar(@{$c->{hits}});
		for (my $k = 0; $k < $klim; $k++) {
			my %hits     = %{$c->{hits}->[$k]}     if defined($c->{hits}->[$k]);
			my %pairhits = %{$c->{pairhits}->[$k]} if defined($c->{pairhits}->[$k]);
			my $hits_are_superset = $c->{hits_are_superset}->[$k];
			# Check if there are any hits left over
			my $hitsLeft = scalar(keys %hits);
			if($hitsLeft != 0 && !$hits_are_superset) {
				print Dumper(\%hits);
				die "Had $hitsLeft hit(s) left over";
			}
			if($pe && $lastchr eq "") {
				my $pairhitsLeft = scalar(keys %pairhits);
				$pairhitsLeft == 0 || die "Had $pairhitsLeft hit(s) left over";
			}
		}
		
		$c->{hits} = $hitstmp;
		$c->{pairhits} = $pairhitstmp;
	}
}
print "PASSED\n";
