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

my @cases = (


	# Alignment with 1 read gap
	{ name   => "Gap penalties 1",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCTTTGTT" ], # budget = 3 + 12 * 3 = 39
	  args   =>   "-P \"SEED=0,3,1;RDG=29,10\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ],
	  flags  => [ "XT:UU" ] },

	# Alignment with 1 read gap - colorspace
	{ name   => "Gap penalties 1 (colorspace)",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "01102200110" ], # budget = 3 + 12 * 3 = 39
	  #           "TTGTTCTTTGTT"
	  args   =>   "-P \"SEED=0,3,1;RDG=26,10\"",
	  report =>   "-a",
	  hits   => [ { 1 => 1 } ],
	  color  => 1,
	  flags  => [ "XT:UU" ] },

	# Alignment with 1 read gap, but not enough budget
	{ name   => "Gap penalties 2",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCTTTGTT" ], # budget = 3 + 12 * 3 = 39
	  args   =>   "-P \"SEED=0,3,1;RDG=30,10\"",
	  report =>   "-a",
	  hits   => [ { } ],
	  flags  => [ "XT:UU" ] },

	# Alignment with 1 read gap, but not enough budget - colorspace
	{ name   => "Gap penalties 2 (colorspace)",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "01102200110" ], # budget = 3 + 11 * 3 = 36
	  #           "TTGTTCTTTGTT"
	  args   =>   "-P \"SEED=0,3,1;RDG=27,10\"",
	  report =>   "-a",
	  hits   => [ { } ],
	  color  => 1,
	  flags  => [ "XT:UU" ] },

	# Alignment with 1 reference gap
	{ name   => "Gap penalties 3",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGATTTGTT" ], # budget = 3 + 14 * 3 = 45
	  args   =>   "-P \"SEED=0,3,1;RFG=30,15\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ],
	  flags  => [ "XT:UU" ] },

	# Alignment with 1 reference gap - colorspace
	{ name   => "Gap penalties 3 (colorspace)",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "0110232300110" ], # budget = 3 + 13 * 3 = 42
	  #           "TTGTTCGATTTGTT"
	  args   =>   "-P \"SEED=0,3,1;RFG=27,15\"",
	  report =>   "-a",
	  hits   => [ { 1 => 1 } ],
	  color  => 1,
	  flags  => [ "XT:UU" ] },

	# Alignment with 1 reference gap, but not enough budget
	{ name   => "Gap penalties 4",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGATTTGTT" ], # budget = 3 + 14 * 3 = 45
	  args   =>   "-P \"SEED=0,3,1;RFG=30,16\"",
	  report =>   "-a",
	  hits   => [ { } ],
	  flags  => [ "XT:UU" ] },

	# Alignment with 1 reference gap, but not enough budget - colorspace
	{ name   => "Gap penalties 4 (colorspace)",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "0110232300110" ], # budget = 3 + 13 * 3 = 42
	  #           "TTGTTCGATTTGTT"
	  args   =>   "-P \"SEED=0,3,1;RFG=27,16\"",
	  report =>   "-a",
	  hits   => [ { } ],
	  color  => 1,
	  flags  => [ "XT:UU" ] },

	# Alignment with 1 reference gap, but not enough budget
	{ name   => "Gap penalties 5",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGATTTGTT" ], # budget = 3 + 14 * 3 = 45
	  args   =>   "-P \"SEED=0,3,1;RFG=31,15\"",
	  report =>   "-a",
	  hits   => [ { } ],
	  flags  => [ "XT:UU" ] },

	# Alignment with 1 reference gap, but not enough budget - colorspace
	{ name   => "Gap penalties 5 (colorspace)",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "0110232300110" ], # budget = 3 + 13 * 3 = 42
	  #           "TTGTTCGATTTGTT"
	  args   =>   "-P \"SEED=0,3,1;RFG=28,15\"",
	  report =>   "-a",
	  hits   => [ { } ],
	  color  => 1,
	  flags  => [ "XT:UU" ] },

	# Alignment with 1 reference gap and 1 read gap
	{ name   => "Gap penalties 6",
	  ref    => [ "ATTGTTCGTTTGTTCGTA" ],
	  reads  => [ "ATTGTTGTTTGATTCGTA" ], # budget = 3 + 18 * 3 = 57
	  args   =>   "-P \"SEED=0,3,1;RFG=18,10;RDG=19,10\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ],
	  flags  => [ "XT:UU" ] },

	# Alignment with 1 reference gap and 1 read gap - colorspace
	{ name   => "Gap penalties 6 (colorspace)",
	  ref    => [ "ATTGTTCGTTTGTTCGTA" ],
	  reads  => [ "30110110012302313" ], # budget = 3 + 17 * 3 = 54
	  #           "ATTGTTGTTTGATTCGTA"
	  args   =>   "-P \"SEED=0,3,1;RFG=18,10;RDG=16,10\"",
	  report =>   "-a",
	  hits   => [ { 1 => 1 } ],
	  color  => 1,
	  flags  => [ "XT:UU" ] },

	# Alignment with 1 reference gap and 1 read gap, but not enough budget
	{ name   => "Gap penalties 7",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTGTTTGATTCGT" ], # budget = 3 + 16 * 3 = 51
	  args   =>   "-P \"SEED=0,3,1;RFG=16,10;RDG=16,10\"",
	  report =>   "-a",
	  hits   => [ { } ],
	  flags  => [ "XT:UU" ] },

	# Alignment with 1 reference gap and 1 read gap, but not enough budget -
	# colorspace
	{ name   => "Gap penalties 7 (colorspace)",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "011011001230231" ], # budget = 3 + 15 * 3 = 48
	  #           "TTGTTGTTTGATTCGT"
	  args   =>   "-P \"SEED=0,3,1;RFG=16,10;RDG=13,10\"",
	  report =>   "-a",
	  hits   => [ { } ],
	  color  => 1,
	  flags  => [ "XT:UU" ] },

	# Experiment with N filtering
	
	{ name => "N filtering 1",
	  ref      => [ "GAGACTTTATACGCATCGAACTATCGCTCTA" ],
	  reads    => [         "ATACGCATCGAAC" ],
	  #              0123456789012345678901234567890
	  #                        1         2         3
	  args     =>   "-P \"NCEIL=0,0\"",
	  report   =>   "-a",
	  hits     => [ { 8 => 1 } ],
	  flags => [ "XT:UU" ] },

	{ name => "N filtering 2",
	  ref      => [ "GAGACTTTATNCGCATCGAACTATCGCTCTA" ],
	  reads    => [         "ATACGCATCGAAC" ],
	  #              0123456789012345678901234567890
	  #                        1         2         3
	  args     =>   "-P \"NCEIL=0,0\"",
	  report   =>   "-a",
	  hits     => [ { } ] },

	{ name => "N filtering 3",
	  ref      => [ "GAGACTTTATACGCATCGAANTATCGCTCTA" ],
	  reads    => [         "ATACGCATCGAAC" ],
	  #              0123456789012345678901234567890
	  #                        1         2         3
	  args     =>   "-P \"NCEIL=0,0\"",
	  report   =>   "-a",
	  hits     => [ { } ] },

	{ name => "N filtering 4",
	  ref      => [ "GAGACTTTNTACGCATCGAACTATCGCTCTA" ],
	  reads    => [         "ATACGCATCGAAC" ],
	  #              0123456789012345678901234567890
	  #                        1         2         3
	  args     =>   "-P \"NCEIL=0,0\"",
	  report   =>   "-a",
	  hits     => [ { } ] },

	{ name => "N filtering 5",
	  ref      => [ "GAGACTTTATNCGCATCGAACTATCGCTCTA" ],
	  reads    => [         "ATACGCATCGAAC" ],
	  #              0123456789012345678901234567890
	  #                        1         2         3
	  args     =>   "-P \"NCEIL=0,0.1;SEED=0,10,1\"",
	  report   =>   "-a",
	  hits     => [ { 8 => 1 } ] },

	{ name => "N filtering 6",
	  ref      => [ "GAGACTTTNTACGCATCGAANTATCGCTCTA" ],
	  reads    => [         "ATACGCATCGAAC" ],
	  #              0123456789012345678901234567890
	  #                        1         2         3
	  args     =>   "-P \"NCEIL=0,0.1;SEED=0,10,1\"",
	  report   =>   "-a",
	  hits     => [ { } ] },

	# No discordant alignment because one mate is repetitive.

	{ name => "Simple paired-end 15",
	  ref    => [ "TTTATAAAAATATTTTTTATAAAAATATTTTCCCCCCCCCCCGCCGGCGCGCCCCCCGCCGGCGCGCCCCC" ],
	#                 ATAAAAATAT     ATAAAAATAT             CGCCGGCGCG     CGCCGGCGCG
	#              01234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7
	#                 -------------------------------------------
	#                 0123456789012345678901234567890123456789012
	  mate1s => [ "ATAAAAATAT" ],
	  mate2s => [ "CGCCGGCGCG" ],
	  args   =>   "--ff -I 0 -X 70 -m 2 --print-placeholders",
	  report =>   "-a",
	  # Not really any way to flag an alignment as discordant
	  pairhits => [ { 3 => 1, 18 => 1, 41 => 1, 56 => 1 } ],
	  flags  => [ "XP:1,XT:UP" # Pair aligns repetitively
	                           ] },

	# No discordant alignment because one mate is repetitive.

	{ name => "Simple paired-end 14",
	  ref    => [ "TTTATAAAAATATTTCCCCCCCGCCCGCCCGCCCCCCGCCCGCCCGCCCCC" ],
	#                 ATAAAAATAT        CGCCCGCCCG     CGCCCGCCCG
	#              012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5
	#                 -------------------------------------------
	#                 0123456789012345678901234567890123456789012
	  mate1s => [ "ATAAAAATAT" ],
	  mate2s => [ "CGCCCGCCCG" ],
	  args   =>   "--ff -I 0 -X 50 -m 1",
	  report =>   "-a",
	  # Not really any way to flag an alignment as discordant
	  pairhits => [ { 3 => 1 } ],
	  flags  => [ "XP:1,XT:UP" # Pair aligns repetitively
	                           ] },

	# 1 discordant alignment and one concordant alignment.  Discordant because
	# the fragment is too long.

	{ name => "Simple paired-end 13",
	  ref    => [ "TTTATAAAAATATTTCCCCCCCCCCCCCCTGTCGCTACCGCCCCCCCCCCC" ],
	#                 ATAAAAATAT                 GTCGCTACCG
	#                 ATAAAAATAT                TGTCGCTACC
	#                 ATAAAAATAT               CTGTCGCTAC
	#                 ATAAAAATAT              CCTGTCGCTA
	#                  TAAAAATATT                GTCGCTACCG
	#                  TAAAAATATT               TGTCGCTACC
	#                  TAAAAATATT              CTGTCGCTAC
	#                  TAAAAATATT             CCTGTCGCTA
	#              012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5
	#                 -----------------------------------
	#                 012345678901234567890123456789012345678901234567
	#                 0         1         2         3         4
	  mate1s   => [ "ATAAAAATAT", "ATAAAAATAT", "ATAAAAATAT", "ATAAAAATAT",
	                "TAAAAATATT", "TAAAAATATT", "TAAAAATATT", "TAAAAATATT", ],
	  mate2s   => [ "GTCGCTACCG", "TGTCGCTACC", "CTGTCGCTAC", "CCTGTCGCTA",
	                "GTCGCTACCG", "TGTCGCTACC", "CTGTCGCTAC", "CCTGTCGCTA" ],
	  args     =>   "--ff -I 0 -X 35",
	  report   =>   "-a",
	  # Not really any way to flag an alignment as discordant
	  pairhits => [ { "3,30" => 1 }, { "3,29" => 1 }, { "3,28" => 1 }, { "3,27" => 1 },
	                { "4,30" => 1 }, { "4,29" => 1 }, { "4,28" => 1 }, { "4,27" => 1 } ],
	  flags    => [ "XT:DP",      "XT:DP",      "XT:CP",      "XT:CP",
	                "XT:DP",      "XT:CP",      "XT:CP",      "XT:CP" ] },

	# 1 discordant alignment and one concordant alignment.  Discordant because
	# the fragment is too long.

	{ name => "Simple paired-end 12",
	  ref    => [ "TTTATAAAAATATTTCCCCCCCCCCCCCCGGGCCCGCCCGCCCCCCCCCCC" ],
	#                 ATAAAAATAT                 GGCCCGCCCG
	#                 ATAAAAATAT              CCGGGCCCGC
	#              012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5
	#                 -------------------------------------
	#                 012345678901234567890123456789012345678901234567
	  mate1s => [ "ATAAAAATAT", "ATAAAAATAT" ],
	  mate2s => [ "GGCCCGCCCG", "CCGGGCCCGC" ],
	  args   =>   "--ff -I 0 -X 36",
	  report =>   "-a",
	  # Not really any way to flag an alignment as discordant
	  pairhits => [ { "3,30" => 1 }, { "3,27" => 1 } ],
	  flags => [ "XT:DP", "XT:CP" ] },

	# 1 discordant alignment.  Discordant because the fragment is too long.

	{ name => "Simple paired-end 11",
	  ref    => [ "TTTATAAAAATATTTCCCCCCCCCCCCCCCCGATCGCCCGCCCCCCCCCCC" ],
	#                 ATAAAAATAT                 CGATCGCCCG
	#              012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5
	#                 -------------------------------------
	#                 012345678901234567890123456789012345678901234567
	  mate1s => [ "ATAAAAATAT" ],
	  mate2s => [ "CGATCGCCCG" ],
	  args   =>   "--ff -I 0 -X 36",
	  report =>   "-a",
	  # Not really any way to flag an alignment as discordant
	  pairhits => [ { "3,30" => 1 } ],
	  flags => [ "XT:DP" ] },

	# 1 discordant alignment.  Discordant because the fragment is too short.

	{ name => "Simple paired-end 10",
	  ref    => [ "TTTATAAAAATATTTCCCCCCGATCGCCCGCCCCCCCCCCC" ],
	#                 ATAAAAATAT       CGATCGCCCG
	#              01234567890123456789012345678901234567890
	#              0         1         2         3         4
	#                 ---------------------------
	#                 012345678901234567890123456
	  mate1s => [ "ATAAAAATAT" ],
	  mate2s => [ "CGATCGCCCG" ],
	  args   =>   "--ff -I 28 -X 80",
	  report =>   "-a",
	  # Not really any way to flag an alignment as discordant
	  pairhits => [ { "3,20" => 1 } ],
	  flags => [ "XT:DP" ] },

	# Like 6, but with -M limit

	{ name => "Simple paired-end 9",
	  ref    => [ "CCCATATATATATCCTCCCATATATATATCCCTCCCCATATATATATCCCTTTTCCTTTCGCGCGCGCGTTTCCCCCCCCC" ],
	#                 ATATATATAT      ATATATATAT        ATATATATAT            CGCGCGCGCG
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CGCGCGCGCG" ],
	  args   =>   "--fr -I 0 -X 80",
	  report =>   "-a -M 2",
	  lines  => 2,
	  pairhits => [ { "3,59" => 1, "19,59" => 1, "37,59" => 1 } ],
	  hits_are_superset => [ 1 ],
	  flags  => [ "XM:1,XP:1,XT:CP", "XM:1,XP:1,XT:CP" ] },

	# Like 6, but without -m limit

	{ name => "Simple paired-end 8",
	  ref    => [ "CCCATATATATATCCTCCCATATATATATCCCTTCCCATATATATATCCCTTTTTTTTTCGCGCGCGCGTTTCCCCCCCCC" ],
	#                 ATATATATAT      ATATATATAT        ATATATATAT            CGCGCGCGCG
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CGCGCGCGCG" ],
	  args   =>   "--fr -I 0 -X 80",
	  report =>   "-a",
	  pairhits => [ { "3,59" => 1, "19,59" => 1, "37,59" => 1 } ],
	  flags  => [ "XT:CP" ] },

	# Like 6, but with lower -m limit

	{ name => "Simple paired-end 7",
	  ref    => [ "TTTATATATATATTTTTTTATATATATATTTTTTTTTATATATATATTTTTTTTCCCCCCGCGCGCGCGCCCCCCCCCCCC" ],
	#                 ATATATATAT      ATATATATAT        ATATATATAT            CGCGCGCGCG
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CGCGCGCGCG" ],
	  args   =>   "--fr -I 0 -X 80 -m 1 --print-placeholders",
	  report =>   "-a",
	  pairhits => [ { "*,*" => 1 } ],
	  flags  => [ "XM:1,XP:1,XT:CP", "XM:1,XP:1,XT:CP" ] },

	# Same but with mates reversed, first mate aligns 3 times

	{ name => "Simple paired-end 6",
	  ref    => [ "CCCATATATATATCCCCCCATATATATATCCCTTCCCATATATATATCCCTTTTCCTTTCGCGCGCGCGTTTCCCCCCCCC" ],
	#                 ATATATATAT      ATATATATAT        ATATATATAT            CGCGCGCGCG
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CGCGCGCGCG" ],
	  args   =>   "--fr -I 0 -X 80 -m 2",
	  report =>   "-a",
	  pairhits => [ { 59 => 2 } ],
	  flags  => [ "XP:1,XT:UP" ] },

	# Same but with mates reversed, second mate doesn't align

	{ name => "Simple paired-end 5",
	  ref    => [ "CCCATATATATATCCCTCCATATATATATCCCTTTCCATATATATATCCCTTTTCCCTTCGCGCGCGCGTTTCCCCCCCCC" ],
	#                 ATATATATAT      ATATATATAT        ATATATATAT            CGCGCGCGCG
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CCCCCGGGGG" ],
	  args   =>   "--fr -I 0 -X 80 -m 2 --print-placeholders",
	  report =>   "-a",
	  pairhits => [ { } ],
	  flags  => [ "XM:1,XT:UP", "XT:UP" ] },

	# Paired-end read, but only the first mate aligns within the -m 2 ceiling.
	# Second mate aligns 3 places.

	{ name => "Simple paired-end 4",
	  ref    => [ "CCCATATATATATCCCTTTTTTTCCCCCCCCCTTTCGCGCGCGCGTTTCTTTCGCGCGCGCGTTTCCCTTTCGCGCGCGCG" ],
	#                 ATATATATAT                      CGCGCGCGCG       CGCGCGCGCG         CGCGCGCGCG
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CGCGCGCGCG" ],
	  args   =>   "--ff -I 0 -X 80 -m 2",
	  report =>   "-a",
	  pairhits => [ { 3 => 2 } ],
	  flags  => [ "XP:1,XT:UP" ] },

	# Paired-end read, but only the first mate aligns within the -m 2 ceiling.
	# Second mate doesn't align at all.

	{ name => "Simple paired-end 3",
	  ref    => [ "CCCATATATATATCCCTTTTTTTCCCCCCCCCTTTCGCGCGCGCGTTTCCTTCGCGCGCGCGTTTTCCTTTCGCGCGCGCG" ],
	#                 ATATATATAT                      CGCGCGCGCG       CGCGCGCGCG         CGCGCGCGCG
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CCCCCGGGGG" ],
	  args   =>   "--ff -I 0 -X 80 -m 2",
	  report =>   "-a",
	  pairhits => [ { 3 => 2 } ],
	  flags => [ "XT:UP" ] },

	# Paired-end read, but only one mate aligns

	{ name => "Simple paired-end 2; no --no-mixed",
	  ref    => [ "CCCATATATATATCCCTTTTTTTCCCCCCCCCCTTCGCGCGCGCGTTTCCCCC" ],
	#                 ATATATATAT                      CGCGCGCGCG
	#              01234567890123456789012345678901234567890123456789012
	#              0         1         2         3         4         5
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CCCCCGGGGG" ],
	  args   =>   "--ff -I 0 -X 50",
	  report =>   "-a",
	  pairhits => [ { 3 => 2 } ],
	  flags => [ "XT:UP" ] },

	{ name => "Simple paired-end 2; --no-mixed",
	  ref    => [ "CCCATATATATATCCCTTTTTTTCCCCCCCCTTTTCGCGCGCGCGTTTCCCCC" ],
	#                 ATATATATAT                      CGCGCGCGCG
	#              01234567890123456789012345678901234567890123456789012
	#              0         1         2         3         4         5
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CCCCCGGGGG" ],
	  args   =>   "--ff -I 0 -X 50 --no-mixed",
	  report =>   "-a",
	  pairhits => [ ] },

	# Simple paired-end alignment
	
	{ name => "Simple paired-end 1",
	  ref    => [ "CCCATATATATATCCCTTTTTTTCCCCCCCCTTTTCGCGCGCGCGTTTTCCCC" ],
	#                 ATATATATAT                      CGCGCGCGCG
	#              01234567890123456789012345678901234567890123456789012
	#              0         1         2         3         4         5
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CGCGCGCGCG" ],
	  args   =>   "--ff -I 0 -X 50",
	  report =>   "-a",
	  pairhits => [ { "3,35" => 1 } ],
	  flags => [ "XT:CP" ] },

	#
	# Alignment with overhang
	#

	# A simple case that should align with or without overhang, with or without
	# a special NCEIL setting.
	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TTGTTCG" ],
	  args   => "",
	  report => "-a --overhang -P \"SEED=0,2,1;NCEIL=2,0\"",
	  hits   => [ { 0 => 1 } ],
	  flags => [ "XT:UU" ] },

	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TTGTTCG" ],
	  args   => "",
	  report => "-a -P \"SEED=0,2,1;NCEIL=2,0\"",
	  hits   => [ { 0 => 1 } ],
	  flags => [ "XT:UU" ] },
	
	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TTGTTCG" ],
	  args   => "",
	  report => "-a --overhang",
	  hits   => [ { 0 => 1 } ],
	  flags => [ "XT:UU" ] },

	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TTGTTCG" ],
	  args   => "",
	  report => "-a",
	  hits   => [ { 0 => 1 } ],
	  flags => [ "XT:UU" ] },

	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TGTTCGT", "TTGTTCG" ],
	  args   => "",
	  report => "-a --overhang",
	  hits   => [ { 1 => 1 }, { 0 => 1 } ],
	  flags => [ "XT:UU", "XT:UU" ] },

	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TGTTCGT", "TTGTTCG" ],
	  args   => "",
	  report => "-a",
	  hits   => [ { 1 => 1 }, { 0 => 1 } ],
	  flags => [ "XT:UU", "XT:UU" ] },

	# Reads 1 and 2 don't have overhang, reads 3 and 4 overhang opposite ends
	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TGTTCGT", "TTGTTCG", "GTTCGTA", "ATTGTTC" ],
	  args   => "",
	  report => "-a --overhang -P \"SEED=0,2,1;NCEIL=2,0\"",
	  hits   => [ { 1 => 1 }, { 0 => 1 }, { 2 => 1 }, { -1 => 1 } ],
	  flags => [ "XT:UU", "XT:UU", "XT:UU", "XT:UU" ] },

	# Same as previous case but --overhang not specified
	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TGTTCGT", "TTGTTCG", "GTTCGTA", "ATTGTTC" ],
	  args   => "",
	  report => "-a -P \"SEED=0,2,1;NCEIL=2,0\"",
	  hits   => [ { 1 => 1 }, { 0 => 1 } ], # only the internal hits
	  flags => [ "XT:UU", "XT:UU" ] },

	# Alignment with 1 reference gap
	{ ref    => [ "TTTTGTTCGTTTG" ],
	  reads  => [ "TTTTGTTCGATTTG" ], # budget = 3 + 14 * 3 = 45
	  args   =>   "-P \"SEED=0,8,1;RFG=25,20\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ],
	  flags => [ "XT:UU" ] },

	# Alignment with 1 reference gap
	{ ref    => [ "TTGTTCGTTTGTT" ],
	  reads  => [ "TTGTTCGATTTGTT" ], # budget = 3 + 14 * 3 = 45
	  args   =>   "-P \"SEED=0,3,1;RFG=25,20\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ],
	  flags => [ "XT:UU" ] },

	{ ref    => [ "ACNCA" ],
	  reads  => [ "CA" ],
	  args   => "",
	  report => "-a -P \"SEED=0,2,1;NCEIL=0,0\"",
	  hits   => [ { 3 => 1 } ],
	  edits  => [ ],
	  flags => [ "XT:UU" ] },

	{ ref    => [ "ACNCA" ],
	  reads  => [ "AC" ],
	  args   => "",
	  report => "-a -P \"SEED=0,2,1;NCEIL=0,0\"",
	  hits   => [ { 0 => 1 } ],
	  edits  => [ ],
	  flags => [ ] },

	{ ref    => [ "ACNCANNNNNNNNCGNNNNNNNNCG" ],
	#              0123456789012345678901234
	#              0         1         2
	  reads  => [ "CG" ],
	  args   => "",
	  report => "-a -P \"SEED=0,2,1;NCEIL=0,0\"",
	  hits   => [ { 13 => 2, 23 => 2 } ],
	  edits  => [ ],
	  flags => [ "XT:UU", "XT:UU" ] },

	{ ref    => [ "ACNCANNNNNNAACGNNNNNNNACGAANNNNCGAAAN" ],
	#              0123456789012345678901234567890123456
	#              0         1         2         3
	  reads  => [ "CG" ],
	  args   => "",
	  report => "-a -P \"SEED=0,2,1;NCEIL=0,0\"",
	  hits   => [ { 13 => 2, 23 => 2, 31 => 2 } ],
	  edits  => [ ],
	  flags => [ "XT:UU", "XT:UU", "XT:UU" ] },

	{ ref    => [ "ACNCANNNNNNAACGNNNNNNNACGAANNNNCGAAAN" ],
	#              0123456789012345678901234567890123456
	#              0         1         2         3
	  reads  => [ "CG" ],
	  args   => "",
	  report => "-a -P \"SEED=0,1,1;NCEIL=0,0\"",
	  hits   => [ { 13 => 2, 23 => 2, 31 => 2 } ],
	  edits  => [ ],
	  flags => [ "XT:UU", "XT:UU", "XT:UU" ] },

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
	  edits  => [ "5:N>G", "5:N>C" ],
	  flags => [ "XT:UU", "XT:UU" ] },

	#
	# Alignment with multi-character read gap
	#

	# Reads 1 and 2 don't have overhang, reads 3 and 4 overhang opposite ends
	{ ref    => [ "ATATGCCCCATGCCCCCCTCCG" ],
	  reads  => [ "ATATGCCCCCCCCCCTCCG" ], # 3 * 19 + 3 = 60
	  #                     ^
	  #                     9:ATG>- 
	  args   => "",
	  report => "-a --overhang -P \"RDG=5,5;SEED=0,8,1\"",
	  hits   => [ { 0 => 1 } ],
	  edits  => [ "9:ATG>-" ],
	  norc   => 1,
	  flags => [ "XT:UU" ] },

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
	  norc   => 1,
	  flags => [ "XT:UU" ] },

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
	  args     =>   "--col-keepends --norc --overhang -P \"NCEIL=1,0;MIN=-10,0;SEED=0,2,1;MMP=C5;SNP=9\"",
	  report   =>   "-a",
	  hits     => [ { 1 => 1 } ],
	  color    => 1,
	  flags => [ "XT:UU" ] },

	{ ref      => [ "TTGTTC" ],
	  reads    => [ "01102"  ],
	  dec_seq  => [ "TTGTTC" ],
	  dec_qual => [ "IqqqqI" ],
	  args     =>   "--col-keepends",
	  report   =>   "-a",
	  hits     => [ { 0 => 1 } ],
	  color    => 1,
	  flags => [ "XT:UU" ] },

	{ ref      => [ "TTGTTC" ],
	  reads    => [ "01102"  ],
	  dec_seq  => [ "TTGTTC" ],
	  dec_qual => [ "IqqqqI" ],
	  args     =>   "--col-keepends -P \"SEED=0,3,1\"",
	  report   =>   "-a",
	  hits     => [ { 0 => 1 } ],
	  color    => 1,
	  flags => [ "XT:UU" ] },

	{ ref      => [ "TTGTTC" ],
	  reads    => [ "01102" ],
	  dec_seq  => [ "TGTT" ],
	  dec_qual => [ "qqqq" ],
	  args     =>   "-P \"SEED=0,3,1\"",
	  report   =>   "-a",
	  hits     => [ { 1 => 1 } ],
	  color    => 1,
	  flags => [ "XT:UU" ] },

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
	  color    => 1,
	  flags => [ "XT:UU" ] },

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
	  color    => 1,
	  flags => [ "XT:UU" ] },

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
	  color    => 1,
	  flags => [ "XT:UU" ] },

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
	  edits  => [ "-",        "-",        "-",        "-",        "-",       ],
	  flags  => [ "XT:UU",    "XT:UU",    "XT:UU",    "XT:UU",    "XT:UU" ]  },

	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TTGTTCGT", "TAAAACGT", "TTGTTCGT", "TAATTCGT",    "TTGTTCGT" ],
	  args   => "",
	  report => "-m 1 -a -P \"SEED=0,3,1;MMP=C9\" --print-placeholders",
	  hits   => [ { 0 => 1 }, { '*' => 1 }, { 0 => 1 }, { 0 => 1 },    { 0 => 1 } ],
	  edits  => [ "-",        undef,        "-",        "1:T>A,2:G>A", "-",        ],
	  flags  => [ "XT:UU",    "XT:UU",      "XT:UU",    "XT:UU",       "XT:UU"     ],
	  norc   => 1 },

	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => "",
	  report => "-m 3 -a",
	  hits   => [ { 0 => 1, 8 => 1 } ],
	  flags  => [ "XT:UU", "XT:UU" ] },

	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => "",
	  report => "-m 2 -a",
	  hits   => [ { 0 => 1, 8 => 1 } ],
	  flags  => [ "XT:UU", "XT:UU" ] },

	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => "",
	  report => "-m 1 -a",
	  hits   => [ { } ],
	  flags  => [ "XM:1,XT:UU" ] },

	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => "",
	  report => "-M 1",
	  hits   => [ { 0 => 1, 8 => 1 } ],
	  flags  => [ "XM:1,XT:UU" ],
	  hits_are_superset => [ 1 ] },

	# Mess with arguments
	
	# Default should be 1-mismatch, so this shouldn't align
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTATTAGT" ],
	  args   => "",
	  report => "-a",
	  hits   => [ {  } ],
	  flags  => [ "XT:UU" ] },

	# Shouldn't align with 0 mismatches either
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTATTAGT" ],
	  args   => "-P SEED=0",
	  report => "-a",
	  hits   => [ { } ],
	  flags  => [ "XT:UU" ] },

	# Should align with 2 mismatches, provided the mismatches are in the right spots
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TAGTTCAT" ],
	  args   => "-P \"SEED=2;MMP=C1\"",
	  report => "-a",
	  hits   => [ { 0 => 1, 8 => 1 } ],
	  flags  => [ "XT:UU" ] },

	# Should align with 2 mismatches, provided the mismatches are in the right spots
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTATTAGT" ],
	  args   =>   "-P \"SEED=2;MMP=C1\"",
	  report =>   "-a",
	  hits   => [ { } ],
	  flags  => [ "XT:UU" ] },

	# Should align with 0 mismatches if we can wedge a seed into the 2
	# matching characters between the two mismatches.  Here we fail to
	# wedge a length-3 seed in (there's no room)
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTATTAGT" ],
	  args   => "-P \"SEED=0,3,1;MMP=C1\"",
	  report => "-a",
	  hits   => [ { } ],
	  flags  => [ "XT:UU" ] },

	# Should align with 0 mismatches if we can wedge a seed into the 2
	# matching characters between the two mismatches.  Here we wedge a
	# length-2 seed in
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [      "TTATTAGT" ],
	  args   => "-P \"SEED=0,2,1;MMP=C1\"",
	  report => "-a",
	  hits   => [ { 0 => 1, 3 => 1, 4 => 1, 5 => 1, 7 => 1, 8 => 1} ],
	  flags  => [ "XT:UU" ] },

	# Following cases depend on this being the case:
	#
	# static const float DEFAULT_CEIL_CONST = 3.0f;
	# static const float DEFAULT_CEIL_LINEAR = 3.0f;

	# Just enough budget for hits, so it should align
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCAT" ], # budget = 3 + 8 * 3 = 27
	  args   => "-P \"SEED=2;MMP=C27\"", # penalty = 27
	  report => "-a",
	  hits   => [ { 0 => 1, 8 => 1 } ],
	  flags  => [ "XT:UU" ] },

	# Not quite enough budget for hits, so it should NOT align
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCAT" ], # budget = 3 + 8 * 3 = 27
	  args   =>   "-P \"SEED=2;MMP=C28\"", # penalty = 28
	  report =>   "-a",
	  hits   => [ { } ],
	  flags  => [ "XT:UU" ] },
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
		for (0..scalar(@$mate1s)-1) {
			my $m1 = $mate1s->[$_];
			my $m2 = $mate2s->[$_];
			my $nm = $names->[$_];
			print FQ1 "\@$nm/1\n$m1\n+\n".("I" x length($m1))."\n";
			print FQ2 "\@$nm/2\n$m2\n+\n".("I" x length($m2))."\n";
		}
	} else {
		for (0..scalar(@$reads)-1) {
			my $read = $reads->[$_];
			my $nm = $names->[$_];
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
	$args .= " --print-flags";
	$reportargs = $reportargs || "-a";
	$args .= " -C" if $color;
	$args .= " $reportargs";
	# Write the reference to a fasta file
	my $build_args = ($color ? "-C" : "");
	print "References:\n";
	open(FA, $fa) || die;
	while(<FA>) { print $_; }
	close(FA);
	my $cmd = "$bowtie2_build --quiet --sanity $build_args $fa .simple_tests.tmp";
	print "$cmd\n";
	system($cmd);
	($? == 0) || die "Bad exitlevel from bowtie2-build: $?";
	my $pe = (defined($mate1s) && $mate1s ne "");
	my $mate1arg;
	my $mate2arg;
	my $readarg;
	my $formatarg = "-c";
	my ($readstr, $m1str, $m2str) = (undef, undef, undef);
	$readstr = join(",", @$reads)  if defined($reads);
	$m1str   = join(",", @$mate1s) if defined($mate1s);
	$m2str   = join(",", @$mate2s) if defined($mate2s);
	if(defined($names)) {
		writeReads($reads, $mate1s, $mate2s, $names, ".simple_tests.1.fq", ".simple_tests.2.fq");
		$mate1arg = ".simple_tests.1.fq";
		$mate2arg = ".simple_tests.2.fq";
		$formatarg = "-q";
		$readarg = $mate1arg;
	} else {
		$mate1arg = $m1str;
		$mate2arg = $m2str;
		$readarg = $readstr;
	}
	if($pe) {
		# Paired-end case
		$cmd = "$bowtie2 $args .simple_tests.tmp $formatarg -1 $mate1arg -2 $mate2arg";
	} else {
		# Unpaired case
		$cmd = "$bowtie2 $args .simple_tests.tmp $formatarg $readarg";
	}
	print "$cmd\n";
	open(BT, "$cmd |") || die "Could not open pipe '$cmd |'";
	while(<BT>) {
		print $_;
		chomp;
		push @$ls,    [ split(/\t/, $_, -1) ];
		push @$rawls, $_;
	}
	close(BT);
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
		my $m1s   = $c->{mate1s};
		my $m2s   = $c->{mate2s};
		if(!$fw) {
			# Reverse-complement the reads
			my @s = (); @s = @$reads if defined($reads);
			# Reverse-complement mates and switch mate1 with mate2
			my @m1 = (); @m1 = @$m1s if defined($m1s);
			my @m2 = (); @m2 = @$m2s if defined($m2s);
			for(0..scalar(@s )-1) { $s [$_] = DNA::revcomp($s [$_], $color); }
			for(0..scalar(@m1)-1) { $m1[$_] = DNA::revcomp($m1[$_], $color); }
			for(0..scalar(@m2)-1) { $m2[$_] = DNA::revcomp($m2[$_], $color); }
			$reads = \@s if defined($reads);
			$m1s = \@m2 if defined($m1s);
			$m2s = \@m1 if defined($m2s);
		}
		runbowtie2(
			"$a",
			$color,
			$tmpfafn,
			$c->{report},
			$reads,
			$m1s,
			$m2s,
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
		if(defined($c->{lines})) {
			my $l = scalar(@lines);
			$l == $c->{lines} || die "Expected $c->{lines} lines, got $l";
		}
		for my $li (0 .. scalar(@lines)-1) {
			my $l = $lines[$li];
			scalar(@$l) == 9 || die "Bad number of fields; expected 9 got ".scalar(@$l).":\n$rawlines[$li]\n";
			my ($readname, $orient, $chr, $off, $seq, $qual, $oms, $editstr, $flagstr) = @$l;
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
			$rdi = substr($rdi, 1) if substr($rdi, 0, 1) eq "r";
			my $mate = 0;
			if($readname =~ /\//) {
				($rdi, $mate) = split(/\//, $readname);
				defined($rdi) || die;
			}
			if($rdi != int($rdi)) {
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
			# 'flags'
			my $flags = undef;
			$flags = $c->{flags}->[$rdi] if defined($c->{flags}->[$rdi]);
			# 'pairhits'
			my %pairhits = ();
			%pairhits = %{$c->{pairhits}->[$rdi]} if defined($c->{pairhits}->[$rdi]);
			# 'pairflags'
			my %pairflags = ();
			%pairflags = %{$c->{pairflags}->[$rdi]} if defined($c->{pairflags}->[$rdi]);
			# 'hits_are_superset'
			my $hits_are_superset = 0;
			$hits_are_superset = $c->{hits_are_superset}->[$rdi] if defined($ci);
			# edits
			my $edits = undef;
			$edits = $c->{edits}->[$rdi] if defined($c->{edits}->[$rdi]);
			# flags
			if(defined($flags)) {
				$flagstr eq $flags ||
					die "Expected flags=\"$flags\", got \"$flagstr\"";
			}
			if($pe && $lastchr ne "") {
				my $offkey = $lastoff.",".$off;
				if($lastoff ne "*") {
					$offkey = min($lastoff, $off).",".max($lastoff, $off);
				}
				if(defined($c->{pairhits}->[$rdi])) {
					defined($pairhits{$offkey}) ||
						die "No such paired off as $offkey in pairhits list: ".Dumper(\%pairhits)."\n";
					$c->{pairhits}->[$rdi]->{$offkey}--;
					delete $c->{pairhits}->[$rdi]->{$offkey} if $c->{pairhits}->[$rdi]->{$offkey} == 0;
					%pairhits = %{$c->{pairhits}->[$rdi]};
				}
				($lastchr, $lastoff) = ("", -1);
			} elsif($pe) {
				# Found an unpaired alignment from aligning a pair
				my $foundSe =
					defined($c->{pairhits}->[$rdi]) &&
					$c->{pairhits}->[$rdi]->{$off};
				if($foundSe) {
					$c->{pairhits}->[$rdi]->{$off}--;
					delete $c->{pairhits}->[$rdi]->{$off}
						if $c->{pairhits}->[$rdi]->{$off} == 0;
					%pairhits = %{$c->{pairhits}->[$rdi]};
				} else {
					($lastchr, $lastoff) = ($chr, $off);
				}
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
				my $eds = $l->[7];
				$eds eq $edits || die "For edit string, expected \"$edits\" got \"$eds\"\n";
			}
		}
		# Go through all the per-read 
		my $klim = 0;
		$klim = scalar(@{$c->{hits}}) if defined($c->{hits});
		$klim = max($klim, scalar(@{$c->{pairhits}})) if defined($c->{pairhits});
		for (my $k = 0; $k < $klim; $k++) {
			# For each read
			my %hits     = %{$c->{hits}->[$k]}     if defined($c->{hits}->[$k]);
			my %pairhits = %{$c->{pairhits}->[$k]} if defined($c->{pairhits}->[$k]);
			my $hits_are_superset = $c->{hits_are_superset}->[$k];
			# Check if there are any hits left over
			my $hitsLeft = scalar(keys %hits);
			if($hitsLeft != 0 && !$hits_are_superset) {
				print Dumper(\%hits);
				die "Had $hitsLeft hit(s) left over";
			}
			my $pairhitsLeft = scalar(keys %pairhits);
			if($pairhitsLeft != 0 && !$hits_are_superset) {
				print Dumper(\%pairhits);
				die "Had $pairhitsLeft hit(s) left over";
			}
		}
		
		$c->{hits} = $hitstmp;
		$c->{pairhits} = $pairhitstmp;
	}
}
print "PASSED\n";
