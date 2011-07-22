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
use Test::Deep;

my $bowtie2 = "";
my $bowtie2_build = "";
my $skipColor = 1;

GetOptions(
	"bowtie2=s"       => \$bowtie2,
	"bowtie2-build=s" => \$bowtie2_build,
	"skip-color"      => \$skipColor) || die "Bad options";

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

	# Some tricky SAM FLAGS field tests
	
	{ name     => "SAM paired-end where both mates align 1",
	  ref      => [ "GCACTATCTACGCTTCGGCGTCGGCGAAAAAACGCACGACCGGGTGTGTGACAATCATATATAGCGCGC" ],
	  #              012345678901234567890123456789012345678901234567890123456789012345678
	  #              0         1         2         3         4         5         6
	  mate1s   => [    "CTATCTACGCTTCGGCGTCGGTGA" ],
	  mate2s   =>                                    [ "GATTGTCACACACCCGGTCGT" ],
	  #                 -----------------------------------------------------
	  #                 01234567890123456789012345678901234567890123456789012
	  #                 0         1         2         3         4         5
	  #  0x1    template having multiple fragments in sequencing
	  #  0x2    each fragment properly aligned according to the aligner
	  #  0x4    fragment unmapped
	  #  0x8    next fragment in the template unmapped
	  # 0x10    SEQ being reverse complemented
	  # 0x20    SEQ of the next fragment in the template being reversed
	  # 0x40    the first fragment in the template
	  # 0x80    the last fragment in the template
	  pairhits     => [{ "3,35" => 1 }],
	  norc         => 1,
	  samflags_map => [{ 3 => (1 | 2 | 32 | 64), 35 => (1 | 2 | 16 | 128) }],
	  tlen_map     => [{ 3 => 53, 35 => -53 }] },

	{ name     => "SAM paired-end where both mates align 2",
	  ref      => [ "GCACTATCTACGCTTCGGCGTCGGCGAAAAAACGCACGACCGGGTGTGTGACAATCATATATAGCGCGC" ],
	  #              012345678901234567890123456789012345678901234567890123456789012345678
	  #              0         1         2         3         4         5         6
	  mate1s   => [    "TCACCGACGCCGAAGCGTAGATAG" ],
	  mate2s   =>                                    [ "ACGACCGGGTGTGTGACAATC" ],
	  #                 -----------------------------------------------------
	  #                 01234567890123456789012345678901234567890123456789012
	  #                 0         1         2         3         4         5
	  #  0x1    template having multiple fragments in sequencing
	  #  0x2    each fragment properly aligned according to the aligner
	  #  0x4    fragment unmapped
	  #  0x8    next fragment in the template unmapped
	  # 0x10    SEQ being reverse complemented
	  # 0x20    SEQ of the next fragment in the template being reversed
	  # 0x40    the first fragment in the template
	  # 0x80    the last fragment in the template
	  mate1fw      => 0,
	  mate2fw      => 1,
	  pairhits     => [{ "3,35" => 1 }],
	  norc         => 1,
	  samflags_map => [{ 3 => (1 | 2 | 16 | 64), 35 => (1 | 2 | 32 | 128) }],
	  tlen_map     => [{ 3 => 53, 35 => -53 }] },

	{ name     => "SAM paired-end where both mates align 3",
	  ref      => [ "GCACTATCTACGCTTCGGCGTCGGCGAAAAAACGCACGACCGGGTGTGTGACAATCATATATAGCGCGC" ],
	  #              012345678901234567890123456789012345678901234567890123456789012345678
	  #              0         1         2         3         4         5         6
	  mate1s   => [    "CTATCTACGCTTCGGCGTCGGTGA" ],
	  mate2s   =>                                    [ "ACGACCGGGTGTGTGACAATC" ],
	  #                 -----------------------------------------------------
	  #                 01234567890123456789012345678901234567890123456789012
	  #                 0         1         2         3         4         5
	  #  0x1    template having multiple fragments in sequencing
	  #  0x2    each fragment properly aligned according to the aligner
	  #  0x4    fragment unmapped
	  #  0x8    next fragment in the template unmapped
	  # 0x10    SEQ being reverse complemented
	  # 0x20    SEQ of the next fragment in the template being reversed
	  # 0x40    the first fragment in the template
	  # 0x80    the last fragment in the template
	  mate1fw      => 1,
	  mate2fw      => 1,
	  pairhits     => [{ "3,35" => 1 }],
	  norc         => 1,
	  samflags_map => [{ 3 => (1 | 2 | 64), 35 => (1 | 2 | 128) }],
	  tlen_map     => [{ 3 => 53, 35 => -53 }] },

	{ name     => "SAM paired-end where mate #1 aligns but mate #2 doesn't",
	  ref      => [ "GCACTATCTACGCTTCGGCGTCGGCGAAAAAACGCACGACCGGGTGTGTGACAATCATATATAGCGCGC" ],
	  #              012345678901234567890123456789012345678901234567890123456789012345678
	  #              0         1         2         3         4         5         6
	  mate1s   => [    "CTATCTACGCTTCGGCGTCGGCGA" ],
	  mate2s   =>                                    [ "GATTGTCTTTTCCCGGAAAAATCGT" ],
	  #  0x1    template having multiple fragments in sequencing
	  #  0x2    each fragment properly aligned according to the aligner
	  #  0x4    fragment unmapped
	  #  0x8    next fragment in the template unmapped
	  # 0x10    SEQ being reverse complemented
	  # 0x20    SEQ of the next fragment in the template being reversed
	  # 0x40    the first fragment in the template
	  # 0x80    the last fragment in the template
	  pairhits     => [{ "*,3" => 1 }],
	  norc         => 1,
	  samflags_map => [{ 3 => (1 | 8 | 64), "*" => (1 | 4 | 128) }] },

	{ name     => "SAM paired-end where mate #1 aligns but mate #2 doesn't because its repetitive",
	  ref      => [ "GCACTATCTACGCTTCGGCGTCGGCGAAAAAACGCACGACCGGGTGTGTGACAATCATATATAGCGCGCGCACGACCGGGTGTGTGACAATCATATATAGCGCGC" ],
	  #              012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234
	  #              0         1         2         3         4         5         6         7         8         9        10
	  mate1s   => [    "CTATCTACGCTTCGGCGTCGGCGA" ],
	  mate2s   =>                                    [ "GATTGTCTTTTCCCGGAAAAATCGT" ],
	  report   => "-a -m 1",
	  #  0x1    template having multiple fragments in sequencing
	  #  0x2    each fragment properly aligned according to the aligner
	  #  0x4    fragment unmapped
	  #  0x8    next fragment in the template unmapped
	  # 0x10    SEQ being reverse complemented
	  # 0x20    SEQ of the next fragment in the template being reversed
	  # 0x40    the first fragment in the template
	  # 0x80    the last fragment in the template
	  pairhits     => [{ "*,3" => 1 }],
	  norc         => 1,
	  samflags_map => [{ 3 => (1 | 8 | 64), "*" => (1 | 4 | 128) }] },

	{ name     => "SAM paired-end where neither mate aligns",
	  ref      => [ "GCACTATCTACGCTTCGGCGTCGGCGAAAAAACGCACGACCGGGTGTGTGACAATCATATATAGCGCGC" ],
	  #              012345678901234567890123456789012345678901234567890123456789012345678
	  #              0         1         2         3         4         5         6
	  mate1s   => [    "CTATATACGAAAAAGCGTCGGCGA" ],
	  mate2s   =>                                    [ "GATTGTCTTTTCCCGGAAAAATCGT" ],
	  #  0x1    template having multiple fragments in sequencing
	  #  0x2    each fragment properly aligned according to the aligner
	  #  0x4    fragment unmapped
	  #  0x8    next fragment in the template unmapped
	  # 0x10    SEQ being reverse complemented
	  # 0x20    SEQ of the next fragment in the template being reversed
	  # 0x40    the first fragment in the template
	  # 0x80    the last fragment in the template
	  pairhits     => [{ "*,*" => 1 }],
	  norc         => 1,
	  samoptflags_flagmap => [{
		(1 | 4 | 8 |  64) => { "YT:Z:UP" => 1 },
		(1 | 4 | 8 | 128) => { "YT:Z:UP" => 1 }
	}] },

	{ name     => "SAM paired-end where both mates align, but discordantly",
	  ref      => [ "GCACTATCTACGCTTCGGCGTCGGCGAAAAAACGCACGACCGGGTGTGTGACAATCATATATAGCGCGC" ],
	  #              012345678901234567890123456789012345678901234567890123456789012345678
	  #              0         1         2         3         4         5         6
	  mate1s   => [    "CTATCTACGCTTCGGCGTCGGCGA" ],
	  mate2s   =>                                    [ "ACGACCGGGTGTGTGACAATC" ],
	  #                 -----------------------------------------------------
	  #                 01234567890123456789012345678901234567890123456789012
	  #                 0         1         2         3         4         5
	  #  0x1    template having multiple fragments in sequencing
	  #  0x2    each fragment properly aligned according to the aligner
	  #  0x4    fragment unmapped
	  #  0x8    next fragment in the template unmapped
	  # 0x10    SEQ being reverse complemented
	  # 0x20    SEQ of the next fragment in the template being reversed
	  # 0x40    the first fragment in the template
	  # 0x80    the last fragment in the template
	  pairhits     => [{ "3,35" => 1 }],
	  norc         => 1,
	  samflags_map => [{ 3 => (1 | 64), 35 => (1 | 128) }],
	  tlen_map     => [{ 3 => 0, 35 => 0 }] },

	{ name   => "matchesRef regression 4",
	  ref    => [ "CCGGGTCGTCACGCCCCGCTTGCGTCANGCCCCTCACCCTCCCTTTGTCGGCTCCCACCCCTCCCCATCCGTTGTCCCCGCCCCCGCCCGCCGGGTCGTCACGCCCCGCTTGCGTCANGC",
	              "GCTCGGAATTCGTGCTCCGNCCCGTACGGTT" ],
	  #
	  #       NNNNNGA------A-------------------G-NTTT
	  #            ||||||||||||||||||||||||||||||||||
	  #       CCAAT-ATTTTTAATTTCCTCTATTTTTCTCTCGTCTTG
	  args   => "-P \"SNP=18.5740807665497;NP=Q;RDG=46.3220993654702;RFG=41.3796024365659;MIN=5.57015383125426,-3.28597145122829;NCEIL=0.263054599454459,0.130843661549367;SEED=2,29;IVAL=L,0.0169183264663712,3.75762168662522\" --overhang --trim5 6",
	  reads  => [ "CTTTGCACCCCTCCCTTGTCGGCTCCCACCCATCCCCATCCGTTGTCCCCGCCCCCGCCCGCCGGTCGTCACTCCCCGTTTGCGTCATGCCCCTCACCCTCCCTTTGTCGGCTCGCACCCCTCCCCATCCGTTGTCCCCGCCCCCGCTCTCGGGGTCTTCACGCCCCGCTTGCTTCATGCCCCTCACTCGCACCCCG" ],
	},
	
	{ name   => "matchesRef regression 3",
	  ref    => [ "GAAGNTTTTCCAATATTTTTAATTTCCTCTATTTTTCTCTCGTCTTGNTCTAC" ],
	  #
	  #       NNNNNGA------A-------------------G-NTTT
	  #            ||||||||||||||||||||||||||||||||||
	  #       CCAAT-ATTTTTAATTTCCTCTATTTTTCTCTCGTCTTG
	  args   => "-P \"MMP=R;MIN=8.8,-8.1\" --overhang",
	  reads  => [ "CAAGACGAGAGAAAAATAGAGGAAATTAAAAATATTGG" ],
	},

	{ name   => "matchesRef regression 2",
	  ref    => ["GTTGTCGGCAGCTCTGGATATGTGNTCTCGGGTTTATNTCGTTGTCG",
	             "CCTTGTTNTTAATGCTGCCTGGTTTNG"],
	  args   =>  "-P \"RDG=2.02030755427021,2.81949533273331;MIN=-6.52134769703939,-3.39889659588514;IVAL=L,0.127835912101927\" --overhang --trim5 5",
	  mate1s => ["TCTGGCGGTTGCGAAGGCCCCTGGCGGTTGCTATGTCCTCTGGCGGTTGCGTTGTCGGCAGCTCG"],
	  mate2s => ["AGAACACATATCCAGAGCTGCCGACAACGAAATGAACCCGAGAGCACAAATCCAGAG"] },

	# Regression test for an issue observed once
	{ name   => "matchesRef regression 1",
	  #            0         1         2         3         4         5         6         7
	  #            01234567890123456789012345678901234567890123456789012345678901234567890
	  ref    => [ "AGGTCGACCGAAAGGCCTAGAGGTCGACCGACAATCTGACCATGGGGCGAGGAGCGAGTAC" ],
	  #                       ||||||||||||||||||||||||||||||||||||||||||||||||||
	  reads  => [            "AAGGCCTAGAGGTCGACCGACAATCTGACCATGGGGCGAGGAGCGAGTACTGGTCTGGGG" ],
	  #                       012345678901234567890123456789012345678901234567890123456789
	  #                       0         1         2         3         4         5                  
	  args   => "--overhang" },

	# 1 discordant alignment and one concordant alignment.  Discordant because
	# the fragment is too long.

	{ name => "Discordant with different chromosomes",
	  ref    => [ "TTTATAAAAATATTTCCCCCCCC",
										 "CCCCCCTGTCGCTACCGCCCCCCCCCCC" ],
	#                 ATAAAAATAT                 GTCGCTACCG
	#                 ATAAAAATAT                TGTCGCTACC
	#              01234567890123456789012
	#              0         1         2  
	#                                     0123456789012345678901234567
	#                                     0         1         2
	  mate1s    => [ "ATAAAAATAT", "ATAAAAATAT" ],
	  mate2s    => [ "GTCGCTACCG", "TGTCGCTACC" ],
	  mate1fw   => 1,
	  mate2fw   => 1,
	  args      =>   "-I 0 -X 35",
	  # Not really any way to flag an alignment as discordant
	  pairhits  => [ { "3,7" => 1 }, { "3,6" => 1 } ],
	  rnext_map => [ { 3 => 1, 7 => 0 }, { 3 => 1, 6 => 0 } ],
	  pnext_map => [ { 3 => 7, 7 => 3 }, { 3 => 6, 6 => 3 } ] },

	{ name   => "Fastq 1",
	  ref    => [   "AGCATCGATCAGTATCTGA" ],
	  fastq  => "\@r0\nCATCGATCAGTATCTG\n+\nIIIIIIIIIIIIIIII",
	  hits   => [{ 2 => 1 }] },

	{ name   => "Tabbed 1",
	  ref    => [   "AGCATCGATCAGTATCTGA" ],
	  tabbed => "r0\tCATCGATCAGTATCTG\tIIIIIIIIIIIIIIII",
	  hits   => [{ 2 => 1 }] },

	{ name   => "Fasta 1",
	  ref    => [  "AGCATCGATCAGTATCTGA" ],
	  fasta  => ">r0\nCATCGATCAGTATCTG",
	  hits   => [{ 2 => 1 }] },

	{ name   => "Qseq 1",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  qseq   => join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "X", "Y",
						   "Index",
						   "0", # Mate
						   "CATCGATCAGTATCTG",
						   "IIIIIIIIIIIIIIII",
						   "1"),
	  hits   => [{ 2 => 1 }] },

	{ name   => "Raw 1",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  raw    =>     "CATCGATCAGTATCTG",
	  hits   => [{ 2 => 1 }] },

	# Like Fastq 1 but with extra newline
	{ name   => "Fastq 2",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fastq  => "\@r0\nCATCGATCAGTATCTG\n+\nIIIIIIIIIIIIIIII\n",
	  hits   => [{ 2 => 1 }] },

	{ name   => "Tabbed 1",
	  ref    => [   "AGCATCGATCAGTATCTGA" ],
	  tabbed => "r0\tCATCGATCAGTATCTG\tIIIIIIIIIIIIIIII\n",
	  hits   => [{ 2 => 1 }] },
	
	{ name   => "Fasta 2",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fasta  => ">r0\nCATCGATCAGTATCTG\n",
	  hits   => [{ 2 => 1 }] },

	{ name   => "Qseq 2",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  qseq   => join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "X", "Y",
						   "Index",
						   "0", # Mate
						   "CATCGATCAGTATCTG",
						   "IIIIIIIIIIIIIIII",
						   "1")."\n",
	  hits   => [{ 2 => 1 }] },

	{ name   => "Raw 2",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  raw    => "CATCGATCAGTATCTG\n",
	  hits   => [{ 2 => 1 }] },

	# Like Fastq 1 but with many extra newlines
	{ name   => "Fastq 3",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fastq  => "\n\n\r\n\@r0\nCATCGATCAGTATCTG\r\n+\n\nIIIIIIIIIIIIIIII\n\n",
	  hits   => [{ 2 => 1 }] },

	{ name   => "Tabbed 3",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  tabbed => "\n\n\r\nr0\tCATCGATCAGTATCTG\tIIIIIIIIIIIIIIII\n\n",
	  hits   => [{ 2 => 1 }] },

	{ name   => "Fasta 3",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fasta  => "\n\n\r\n>r0\nCATCGATCAGTATCTG\r\n\n",
	  hits   => [{ 2 => 1 }] },

	{ name   => "Qseq 3",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  qseq   => "\n\n\n".join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "X", "Y",
						   "Index",
						   "0", # Mate
						   "CATCGATCAGTATCTG",
						   "IIIIIIIIIIIIIIII",
						   "1")."\n\n",
	  hits   => [{ 2 => 1 }] },

	{ name   => "Raw 3",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  raw    => "\n\n\nCATCGATCAGTATCTG\n\n",
	  hits   => [{ 2 => 1 }] },

	# Quality string length doesn't match (too short by 1)
	{ name   => "Fastq 4",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fastq  => "\n\n\r\n\@r0\nCATCGATCAGTATCTG\r\n+\n\nIIIIIIIIIIIIIII\n\n",
	  should_abort => 1},

	{ name   => "Tabbed 4",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  tabbed => "\n\n\r\nr0\tCATCGATCAGTATCTG\tIIIIIIIIIIIIIII\n\n",
	  should_abort => 1},

	{ name   => "Qseq 4",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  qseq   => "\n\n\n".join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "X", "Y",
						   "Index",
						   "0", # Mate
						   "CATCGATCAGTATCTG",
						   "IIIIIIIIIIIIIII",
						   "1")."\n\n",
	  should_abort => 1},

	# Name line doesn't start with @
	{ name   => "Fastq 5",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fastq  => "\n\n\r\nr0\nCATCGATCAGTATCTG\r\n+\n\nIIIIIIIIIIIIIII\n\n",
	  should_abort => 1,
	  hits   => [{ }] },

	# Name line doesn't start with >
	{ name   => "Fasta 5",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fasta  => "\n\n\r\nr0\nCATCGATCAGTATCTG\r",
	  should_abort => 1,
	  hits   => [{ }] },

	# Name line doesn't start with @ (2)
	{ name   => "Fastq 6",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fastq  => "r0\nCATCGATCAGTATCTG\r\n+\n\nIIIIIIIIIIIIIII\n\n",
	  should_abort => 1,
	  hits   => [{ }] },

	# Name line doesn't start with > (2)
	{ name   => "Fasta 6",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fasta  => "r0\nCATCGATCAGTATCTG\r",
	  should_abort => 1,
	  hits   => [{ }] },

	# Part of sequence is trimmed
	{ name   => "Fastq 7",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fastq  => "\n\n\r\n\@r0\nCATCGATCAGTATCTG\r\n+\n\nIIIIIIIIIIIIIIII\n\n",
	  args   => "--trim3 4",
	  norc   => 1,
	  hits   => [{ 2 => 1 }] },

	{ name   => "Tabbed 7",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  tabbed => "\n\n\r\nr0\tCATCGATCAGTATCTG\tIIIIIIIIIIIIIIII\n\n",
	  args   => "--trim3 4",
	  norc   => 1,
	  hits   => [{ 2 => 1 }] },

	{ name   => "Fasta 7",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fasta  => "\n\n\r\n\>r0\nCATCGATCAGTATCTG\r\n",
	  args   => "--trim3 4",
	  norc   => 1,
	  hits   => [{ 2 => 1 }] },

	{ name   => "Qseq 7",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  qseq   => "\n\n\n".join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "X", "Y",
						   "Index",
						   "0", # Mate
						   "CATCGATCAGTATCTG",
						   "IIIIIIIIIIIIIIII",
						   "1")."\n\n",
	  args   => "--trim3 4",
	  norc   => 1,
	  hits   => [{ 2 => 1 }] },

	{ name   => "Raw 7",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  raw    => "\n\n\r\nCATCGATCAGTATCTG\r\n",
	  args   => "--trim3 4",
	  norc   => 1,
	  hits   => [{ 2 => 1 }] },

	# Whole sequence is trimmed
	{ name   => "Fastq 8",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fastq  => "\n\n\r\n\@r0\nCATCGATCAGTATCTG\r\n+\n\nIIIIIIIIIIIIIIII\n\n",
	  args   => "--trim5 16",
	  hits   => [{ "*" => 1 }] },

	{ name   => "Tabbed 8",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  tabbed => "\n\n\r\nr0\tCATCGATCAGTATCTG\tIIIIIIIIIIIIIIII\n\n",
	  args   => "--trim5 16",
	  hits   => [{ "*" => 1 }] },

	{ name   => "Fasta 8",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fasta  => "\n\n\r\n\>r0\nCATCGATCAGTATCTG\r\n",
	  args   => "--trim3 16",
	  hits   => [{ "*" => 1 }] },

	{ name   => "Qseq 8",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  qseq   => "\n\n\n".join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "X", "Y",
						   "Index",
						   "0", # Mate
						   "CATCGATCAGTATCTG",
						   "IIIIIIIIIIIIIIII",
						   "1")."\n\n",
	  args   => "--trim3 16",
	  hits   => [{ "*" => 1 }] },

	{ name   => "Raw 8",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  raw    => "\n\n\r\nCATCGATCAGTATCTG\r\n",
	  args   => "--trim3 16",
	  hits   => [{ "*" => 1 }] },

	# Sequence is skipped
	{ name   => "Fastq 9",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fastq  => "\n\n\r\n\@r0\nCATCGATCAGTATCTG\r\n+\n\nIIIIIIIIIIIIIIII\n\n",
	  args   => "-s 1",
	  hits   => [{ }] },

	{ name   => "Tabbed 9",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  tabbed => "\n\n\r\nr0\tCATCGATCAGTATCTG\tIIIIIIIIIIIIIIII\n\n",
	  args   => "-s 1",
	  hits   => [{ }] },

	{ name   => "Fasta 9",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fasta  => "\n\n\r\n>r0\nCATCGATCAGTATCTG\r\n",
	  args   => "-s 1",
	  hits   => [{ }] },

	{ name   => "Qseq 9",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  qseq   => "\n\n\n".join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "X", "Y",
						   "Index",
						   "0", # Mate
						   "CATCGATCAGTATCTG",
						   "IIIIIIIIIIIIIIII",
						   "1")."\n\n",
	  args   => "-s 1",
	  hits   => [{ }] },

	{ name   => "Raw 9",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  raw    => "CATCGATCAGTATCTG\n",
	  args   => "-s 1",
	  hits   => [{ }] },

	# Like Fastq 1 but with many extra newlines
	{ name   => "Fastq multiread 1",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fastq  => "\n\n\r\n\@r0\nCATCGATCAGTATCTG\r\n+\n\nIIIIIIIIIIIIIIII\n\n".
	            "\n\n\r\n\@r1\nATCGATCAGTATCTG\r\n+\n\nIIIIIIIIIIIIIII\n\n",
	  hits   => [{ 2 => 1 }, { 3 => 1 }] },

	{ name   => "Tabbed multiread 1",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  tabbed => "\n\n\r\nr0\tCATCGATCAGTATCTG\tIIIIIIIIIIIIIIII\n\n".
	            "\n\n\r\nr1\tATCGATCAGTATCTG\tIIIIIIIIIIIIIII\n\n",
	  hits   => [{ 2 => 1 }, { 3 => 1 }] },

	{ name   => "Fasta multiread 1",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  fasta  => "\n\n\r\n>r0\nCATCGATCAGTATCTG\n\n".
	            "\n\n\r\n>r1\nATCGATCAGTATCTG\n\n",
	  hits   => [{ 2 => 1 }, { 3 => 1 }] },

	{ name   => "Qseq multiread 1",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  qseq   => "\n\n\n".join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "10", "10",
						   "Index",
						   "0", # Mate
						   "CATCGATCAGTATCTG",
						   "IIIIIIIIIIIIIIII",
						   "1")."\n\n".
						 join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "12", "15",
						   "Index",
						   "0", # Mate
						   "ATCGATCAGTATCTG",
						   "IIIIIIIIIIIIIII",
						   "1")."\n\n",
	  idx_map => { "MachName_RunNum_Lane_Tile_10_10_Index" => 0,
	               "MachName_RunNum_Lane_Tile_12_15_Index" => 1 },
	  hits   => [{ 2 => 1 }, { 3 => 1 }] },
	
	{ name   => "Raw multiread 1",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  raw    => "\n\n\r\nCATCGATCAGTATCTG\n\n".
	            "\n\n\r\nATCGATCAGTATCTG\n\n",
	  hits   => [{ 2 => 1 }, { 3 => 1 }] },

	# Like Fastq multiread 1 but with -u 1
	{ name   => "Fastq multiread 2",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  args   =>   "-u 1",
	  fastq  => "\n\n\r\n\@r0\nCATCGATCAGTATCTG\r\n+\n\nIIIIIIIIIIIIIIII\n\n".
	            "\n\n\r\n\@r1\nATCGATCAGTATCTG\r\n+\n\nIIIIIIIIIIIIIII\n\n",
	  hits   => [{ 2 => 1 }] },

	{ name   => "Tabbed multiread 2",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  args   =>   "-u 1",
	  tabbed => "\n\n\r\nr0\tCATCGATCAGTATCTG\tIIIIIIIIIIIIIIII\n\n".
	            "\n\n\r\nr1\tATCGATCAGTATCTG\tIIIIIIIIIIIIIII\n\n",
	  hits   => [{ 2 => 1 }] },

	{ name   => "Fasta multiread 2",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  args   =>   "-u 1",
	  fasta  => "\n\n\r\n>r0\nCATCGATCAGTATCTG\r\n".
	            "\n\n\r\n>r1\nATCGATCAGTATCTG\r\n",
	  hits   => [{ 2 => 1 }] },

	{ name   => "Qseq multiread 2",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  args   =>   "-u 1",
	  qseq   => "\n\n\n".join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "10", "10",
						   "Index",
						   "0", # Mate
						   "CATCGATCAGTATCTG",
						   "IIIIIIIIIIIIIIII",
						   "1")."\n\n".
						 join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "12", "15",
						   "Index",
						   "0", # Mate
						   "ATCGATCAGTATCTG",
						   "IIIIIIIIIIIIIII",
						   "1")."\n\n",
	  idx_map => { "MachName_RunNum_Lane_Tile_10_10_Index" => 0,
	               "MachName_RunNum_Lane_Tile_12_15_Index" => 1 },
	  hits   => [{ 2 => 1 }] },

	{ name   => "Raw multiread 2",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  args   =>   "-u 1",
	  raw    => "\n\n\r\nCATCGATCAGTATCTG\r\n".
	            "\n\n\r\nATCGATCAGTATCTG\r\n",
	  hits   => [{ 2 => 1 }] },

	# Like Fastq multiread 1 but with -u 2
	{ name   => "Fastq multiread 3",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  args   =>   "-u 2",
	  fastq  => "\n\n\r\n\@r0\nCATCGATCAGTATCTG\r\n+\n\nIIIIIIIIIIIIIIII\n\n".
	            "\n\n\r\n\@r1\nATCGATCAGTATCTG\r\n+\n\nIIIIIIIIIIIIIII\n\n",
	  hits   => [{ 2 => 1 }, { 3 => 1 }] },

	{ name   => "Tabbed multiread 3",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  args   =>   "-u 2",
	  tabbed => "\n\n\r\nr0\tCATCGATCAGTATCTG\tIIIIIIIIIIIIIIII\n\n".
	            "\n\n\r\nr1\tATCGATCAGTATCTG\tIIIIIIIIIIIIIII\n\n",
	  hits   => [{ 2 => 1 }, { 3 => 1 }] },

	{ name   => "Fasta multiread 3",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  args   =>   "-u 2",
	  fasta  => "\n\n\r\n>r0\nCATCGATCAGTATCTG\r\n".
	            "\n\n\r\n>r1\nATCGATCAGTATCTG\r\n",
	  hits   => [{ 2 => 1 }, { 3 => 1 }] },

	{ name   => "Qseq multiread 3",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  args   =>   "-u 2",
	  qseq   => "\n\n\n".join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "10", "10",
						   "Index",
						   "0", # Mate
						   "CATCGATCAGTATCTG",
						   "IIIIIIIIIIIIIIII",
						   "1")."\n\n".
						 join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "12", "15",
						   "Index",
						   "0", # Mate
						   "ATCGATCAGTATCTG",
						   "IIIIIIIIIIIIIII",
						   "1")."\n\n",
	  idx_map => { "MachName_RunNum_Lane_Tile_10_10_Index" => 0,
	               "MachName_RunNum_Lane_Tile_12_15_Index" => 1 },
	  hits   => [{ 2 => 1 }, { 3 => 1 }] },
	
	{ name   => "Raw multiread 3",
	  ref    => [ "AGCATCGATCAGTATCTGA" ],
	  args   =>   "-u 2",
	  raw    => "\n\n\r\nCATCGATCAGTATCTG\r\n".
	            "\n\n\r\nATCGATCAGTATCTG\r\n",
	  hits   => [{ 2 => 1 }, { 3 => 1 }] },
	
	# Paired-end reads that should align
	{ name     => "Fastq paired 1",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  fastq1  => "\n\n\r\n\@r0\nAGCATCGATC\r\n+\n\nIIIIIIIIII\n\n".
	             "\n\n\@r1\nTCAGTTTTTGA\r\n+\n\nIIIIIIIIIII\n\n",
	  fastq2  => "\n\n\r\n\@r0\nTCAGTTTTTGA\n+\n\nIIIIIIIIIII\n\n".
	             "\n\n\r\n\@r1\nAGCATCGATC\r\n+\n\nIIIIIIIIII",
	  pairhits => [ { "0,8" => 1 }, { "0,8" => 1 } ] },

	{ name     => "Tabbed paired 1",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  tabbed  => "\n\n\r\nr0\tAGCATCGATC\tIIIIIIIIII\tTCAGTTTTTGA\tIIIIIIIIIII\n\n".
	             "\n\nr1\tTCAGTTTTTGA\tIIIIIIIIIII\tAGCATCGATC\tIIIIIIIIII\n\n",
	  paired => 1,
	  pairhits => [ { "0,8" => 1 }, { "0,8" => 1 } ] },

	{ name     => "Fasta paired 1",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  fasta1  => "\n\n\r\n>r0\nAGCATCGATC\r\n".
	             "\n\n>r1\nTCAGTTTTTGA\r\n",
	  fasta2  => "\n\n\r\n>r0\nTCAGTTTTTGA\n".
	             "\n\n\r\n>r1\nAGCATCGATC",
	  pairhits => [ { "0,8" => 1 }, { "0,8" => 1 } ] },

	{ name   => "Qseq paired 1",
	  ref    => [       "AGCATCGATCAAAAACTGA" ],
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  qseq1  => "\n\n\n".join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "10", "10",
						   "Index",
						   "1", # Mate
						   "AGCATCGATC",
						   "ABCBGACBCB",
						   "1")."\n\n".
						 join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "12", "15",
						   "Index",
						   "1", # Mate
						   "TCAGTTTTTGA",
						   "95849456875",
						   "1")."\n\n",
	  qseq2  => "\n\n\n".join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "10", "10",
						   "Index",
						   "2", # Mate
						   "TCAGTTTTTGA",
						   "ABCBGACBCBA",
						   "1")."\n\n".
						 join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "12", "15",
						   "Index",
						   "2", # Mate
						   "AGCATCGATC",
						   "AGGCBBGCBG",
						   "1")."\n\n",
	  idx_map => { "MachName_RunNum_Lane_Tile_10_10_Index" => 0,
	               "MachName_RunNum_Lane_Tile_12_15_Index" => 1 },
	  pairhits => [ { "0,8" => 1 }, { "0,8" => 1 } ] },
	
	{ name     => "Raw paired 1",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  raw1    => "\n\n\r\nAGCATCGATC\r\n".
	             "\n\nTCAGTTTTTGA\r\n",
	  raw2    => "\n\n\r\nTCAGTTTTTGA\n".
	             "\n\n\r\nAGCATCGATC",
	  pairhits => [ { "0,8" => 1 }, { "0,8" => 1 } ] },

	# Paired-end reads that should align
	{ name     => "Fastq paired 2",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  args     => "-s 1",
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  fastq1  => "\@r0\nAGCATCGATC\r\n+\n\nIIIIIIIIII\n\n".
	             "\n\n\@r1\nTCAGTTTTTGA\n+\n\nIIIIIIIIIII\n\n",
	  fastq2  => "\n\n\r\n\@r0\nTCAGTTTTTGA\n+\n\nIIIIIIIIIII\n\n".
	             "\n\n\r\n\@r1\nAGCATCGATC\r\n+\n\nIIIIIIIIII",
	  pairhits => [ { }, { "0,8" => 1 } ] },

	{ name     => "Tabbed paired 2",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  args     => "-s 1",
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  tabbed   => "r0\tAGCATCGATC\tIIIIIIIIII\tTCAGTTTTTGA\tIIIIIIIIIII\n\n".
	             "\nr1\tTCAGTTTTTGA\tIIIIIIIIIII\tAGCATCGATC\tIIIIIIIIII",
	  paired   => 1,
	  pairhits => [ { }, { "0,8" => 1 } ] },

	{ name     => "Fasta paired 2",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  args     => "-s 1",
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  fasta1  => ">r0\nAGCATCGATC\r\n".
	             "\n\n>r1\nTCAGTTTTTGA\n",
	  fasta2  => "\n\n\r\n>r0\nTCAGTTTTTGA\n".
	             "\n\n\r\n>r1\nAGCATCGATC",
	  pairhits => [ { }, { "0,8" => 1 } ] },

	{ name   => "Qseq paired 1",
	  ref    => [       "AGCATCGATCAAAAACTGA" ],
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  args     => "-s 1",
	  qseq1  => "\n\n\n".join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "10", "10",
						   "Index",
						   "1", # Mate
						   "AGCATCGATC",
						   "ABCBGACBCB",
						   "1")."\n\n".
						 join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "12", "15",
						   "Index",
						   "1", # Mate
						   "TCAGTTTTTGA",
						   "95849456875",
						   "1")."\n\n",
	  qseq2  => "\n\n\n".join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "10", "10",
						   "Index",
						   "2", # Mate
						   "TCAGTTTTTGA",
						   "ABCBGACBCBA",
						   "1")."\n\n".
						 join("\t", "MachName",
	                       "RunNum",
						   "Lane",
						   "Tile",
						   "12", "15",
						   "Index",
						   "2", # Mate
						   "AGCATCGATC",
						   "AGGCBBGCBG",
						   "1")."\n\n",
	  idx_map => { "MachName_RunNum_Lane_Tile_10_10_Index" => 0,
	               "MachName_RunNum_Lane_Tile_12_15_Index" => 1 },
	  pairhits => [ { }, { "0,8" => 1 } ] },

	{ name     => "Raw paired 2",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  args     => "-s 1",
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  raw1    => "AGCATCGATC\r\n".
	             "\n\nTCAGTTTTTGA\n",
	  raw2    => "\n\n\r\nTCAGTTTTTGA\n".
	             "\n\n\r\nAGCATCGATC",
	  pairhits => [ { }, { "0,8" => 1 } ] },

	# Paired-end reads that should align
	{ name     => "Fastq paired 3",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  args     => "-u 1",
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  fastq1  => "\n\n\r\n\@r0\nAGCATCGATC\r\n+\n\nIIIIIIIIII\n\n".
	             "\n\n\@r1\nTCAGTTTTTGA\r\n+\n\nIIIIIIIIIII\n\n",
	  fastq2  => "\n\n\r\n\@r0\nTCAGTTTTTGA\n+\n\nIIIIIIIIIII\n\n".
	             "\n\n\r\n\@r1\nAGCATCGATC\r\n+\nIIIIIIIIII",
	  pairhits => [ { "0,8" => 1 }, { } ] },

	{ name     => "Tabbed paired 3",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  args     => "-u 1",
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  tabbed   => "\n\n\r\nr0\tAGCATCGATC\tIIIIIIIIII\tTCAGTTTTTGA\tIIIIIIIIIII\n\n".
	              "\n\nr1\tTCAGTTTTTGA\tIIIIIIIIIII\tAGCATCGATC\tIIIIIIIIII",
	  paired   => 1,
	  pairhits => [ { "0,8" => 1 }, { } ] },

	{ name     => "Fasta paired 3",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  args     => "-u 1",
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  fasta1  => "\n\n\r\n>r0\nAGCATCGATC\r\n".
	             "\n\n>r1\nTCAGTTTTTGA\r\n",
	  fasta2  => "\n\n\r\n>r0\nTCAGTTTTTGA\n".
	             "\n\n\r\n>r1\nAGCATCGATC",
	  pairhits => [ { "0,8" => 1 }, { } ] },

	{ name     => "Raw paired 3",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  args     => "-u 1",
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  raw1    => "\n\n\r\nAGCATCGATC\r\n".
	             "\n\nTCAGTTTTTGA\r\n",
	  raw2    => "\n\n\r\nTCAGTTTTTGA\n".
	             "\n\n\r\nAGCATCGATC",
	  pairhits => [ { "0,8" => 1 }, { } ] },

	# Paired-end reads that should align
	{ name     => "Fastq paired 4",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  args     => "-s 1 -P \"SEED=2\"",
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  fastq1  => "\n\n\r\n\@r0\nAGCATCGATC\r\n+\n\nIIIIIIIIII\n\n".
	             "\n\n\@r1\nTC\r\n+\n\nII\n\n".
	             "\n\n\@r2\nTCAGTTTTTGA\r\n+\n\nIIIIIIIIIII\n\n",
	  fastq2  => "\n\n\r\n\@r0\nTCAGTTTTTGA\n+\n\nIIIIIIIIIII\n\n".
	             "\n\n\r\n\@r1\nAG\r\n+\nII".
	             "\n\@r2\nAGCATCGATC\r\n+\nIIIIIIIIII",
	  paired   => 1,
	  pairhits =>    [ { }, { "*,*" => 1 }, { "0,8" => 1 } ],
	  samoptflags_map => [
	  { },
	  { "*" => { "YT:Z:UP" => 1, "YF:i:1"  => 1 } },
	  { 0   => { "MD:Z:10" => 1, "YT:Z:CP" => 1 },
		8   => { "MD:Z:11" => 1, "YT:Z:CP" => 1 } }]
	},

	{ name     => "Tabbed paired 4",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  args     => "-s 1 -P \"SEED=2\"",
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  tabbed   => "\n\n\r\nr0\tAGCATCGATC\tIIIIIIIIII\tTCAGTTTTTGA\tIIIIIIIIIII\n\n".
	              "\n\nr1\tTC\tII\tAG\tII".
	              "\n\nr2\tTCAGTTTTTGA\tIIIIIIIIIII\tAGCATCGATC\tIIIIIIIIII\n\n",
	  paired   => 1,
	  pairhits =>    [ { }, { "*,*" => 1 }, { "0,8" => 1 } ],
	  samoptflags_map => [
	  { },
	  { "*" => { "YT:Z:UP" => 1, "YF:i:1"  => 1 } },
	  { 0   => { "MD:Z:10" => 1, "YT:Z:CP" => 1 },
		8   => { "MD:Z:11" => 1, "YT:Z:CP" => 1 } }]
	},

	{ name     => "Fasta paired 4",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  args     => "-s 1 -P \"SEED=2\"",
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  fasta1  => "\n\n\r\n>r0\nAGCATCGATC\r\n".
	             "\n\n>r1\nTC\r\n".
	             "\n\n>r2\nTCAGTTTTTGA\r\n",
	  fasta2  => "\n\n\r\n>r0\nTCAGTTTTTGA\n\n".
	             "\n\n\r\n>r1\nAG".
	             "\n>r2\nAGCATCGATC",
	  pairhits =>    [ { }, { "*,*" => 1 }, { "0,8" => 1 } ],
	  samoptflags_map => [
	  { },
	  { "*" => { "YT:Z:UP" => 1, "YF:i:1"  => 1 } },
	  { 0   => { "MD:Z:10" => 1, "YT:Z:CP" => 1 },
		8   => { "MD:Z:11" => 1, "YT:Z:CP" => 1 } }]
	},

	{ name     => "Raw paired 4",
	  ref      => [     "AGCATCGATCAAAAACTGA" ],
	  args     => "-s 1 -P \"SEED=2\"",
	  #                  AGCATCGATC
	  #                          TCAAAAACTGA
	  #                  0123456789012345678
	  raw1    => "\n\n\r\nAGCATCGATC\r\n".
	             "\n\nTC\r\n".
	             "\n\nTCAGTTTTTGA\r\n",
	  raw2    => "\n\n\r\nTCAGTTTTTGA\n\n".
	             "\n\n\r\nAG".
	             "\nAGCATCGATC",
	  pairhits =>    [ { }, { "*,*" => 1 }, { "0,8" => 1 } ],
	  samoptflags_map => [
	  { },
	  { "*" => { "YT:Z:UP" => 1, "YF:i:1"  => 1 } },
	  { 0   => { "MD:Z:10" => 1, "YT:Z:CP" => 1 },
		8   => { "MD:Z:11" => 1, "YT:Z:CP" => 1 } }]
	},

	#
	# Check that skipping of empty reads is handled correctly.  A read that is
	# empty or becomes empty after --trim3/--trim5 are applied should still
	# count as a first-class read that gets propagated up into the alignment
	# loop.  And it should be counted in the -s/-u totals.
	#
	
	{ ref      => [     "AGCATCGATCAGTATCTGA" ],
	  reads    => [ "",    "ATCGATCAGTA" ],
	  args     => "-s 1",
	  hits     => [ {}, { 3 => 1 }] },

	{ ref      => [     "AGCATCGATCAGTATCTGA" ],
	  mate1s   => [ "", "AGCATCGATC" ],
	  mate2s   => [ "",          "TCAGATACTG" ],
	  args     => "-s 1",
	  pairhits => [ {}, { "0,9" => 1 }] },

	{ ref      => [     "AGCATCGATCAGTATCTGA" ],
	  reads    => [ "",    "ATCGATCAGTA" ],
	  args     => "-s 2",
	  hits     => [ {}, {} ] },

	{ ref      => [     "AGCATCGATCAGTATCTGA" ],
	  mate1s   => [ "", "AGCATCGATC" ],
	  mate2s   => [ "",          "TCAGATACTG" ],
	  args     => "-s 2",
	  pairhits => [ {}, {} ] },

	{ ref    => [     "AGCATCGATCAGTATCTGA" ],
	  reads  => [ "",    "ATCGATCAGTA", "AGTATCTGA" ],
	  args   => "-s 1 -u 1",
	  hits   => [ {}, { 3 => 1 }] },

	{ ref    => [     "AGCATCGATCAGTATCTGA" ],
	  reads  => [ "AC",  "ATCGATCAGTA" ],
	  args   => "-s 1 --trim3 2",
	  norc   => 1,
	  hits   => [ {}, { 3 => 1 }] },

	{ ref    => [     "AGCATCGATCAGTATCTGA" ],
	  reads  => [ "AC",  "ATCGATCAGTA" ],
	  args   => "-s 1 --trim3 2",
	  nofw   => 1,
	  hits   => [ {}, { 5 => 1 }] },

	{ ref    => [     "AGCATCGATCAGTATCTGA" ],
	  reads  => [ "AC",  "ATCGATCAGTA" ],
	  args   => "-s 1 --trim5 2",
	  nofw   => 1,
	  hits   => [ {}, { 3 => 1 }] },

	{ ref    => [     "AGCATCGATCAGTATCTGA" ],
	  reads  => [ "AC",  "ATCGATCAGTA" ],
	  args   => "-s 1 --trim5 2",
	  norc   => 1,
	  hits   => [ {}, { 5 => 1 }] },

	#
	# Alignment with overhang
	#

	{ ref    => [ "TGC" ],
	  reads  => [ "ATGC" ],
	  args   => "--overhang -P \"SEED=0,3,1;NCEIL=1,0\"",
	  hits   => [ { 0 => 1 } ],
	  mapq   => [ 40 ],
	  cigar  => [ "1S3M" ],
	  samoptflags => [
		{ "AS:i:-1" => 1, "YT:Z:UU" => 1, "MD:Z:3" => 1, "XN:i:1" => 1 } ]
	},

	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TTGTTCG" ],
	  args   => "-P \"SEED=0,2,1;NCEIL=2,0\"",
	  hits   => [ { 0 => 1 } ],
	  mapq   => [ 40 ],
	  cigar  => [ "7M" ],
	  samoptflags => [ { "AS:i:0" => 1, "YT:Z:UU" => 1, "MD:Z:7" => 1 } ]
	},

	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TTGTTCG" ],
	  args   => "",
	  hits   => [ { 0 => 1 } ],
	  mapq   => [ 40 ],
	  flags => [ "XM:0,XP:0,XT:UU,XC:7=" ],
	  cigar  => [ "7M" ],
	  samoptflags => [ { "AS:i:0" => 1, "YT:Z:UU" => 1, "MD:Z:7" => 1 } ]
	},

	{ ref    => [ "TTGTTCGT" ],
	  reads  => [  "TGTTCGT", "TTGTTCG" ],
	  args   => "--overhang",
	  hits   => [ { 1 => 1 }, { 0 => 1 } ],
	  flags => [ "XM:0,XP:0,XT:UU,XC:7=", "XM:0,XP:0,XT:UU,XC:7=" ],
	  cigar  => [ "7M", "7M" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:7" => 1 },
		{ "YT:Z:UU" => 1, "MD:Z:7" => 1 }
	  ]
	},

	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TGTTCGT", "TTGTTCG" ],
	  args   => "",
	  hits   => [ { 1 => 1 }, { 0 => 1 } ],
	  flags => [ "XM:0,XP:0,XT:UU,XC:7=", "XM:0,XP:0,XT:UU,XC:7=" ],
	  cigar  => [ "7M", "7M" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:7" => 1 },
		{ "YT:Z:UU" => 1, "MD:Z:7" => 1 }
	  ]
	},

	# Reads 1 and 2 don't have overhang, reads 3 and 4 overhang opposite ends
	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TGTTCGT", "GTTCGTA", "ATTGTTC" ],
	  args   => "--overhang -P \"SEED=0,2,1;NCEIL=2,0\"",
	  hits   => [ { 1 => 1 }, { 2 => 1 }, { 0 => 1 } ],
	  cigar  => [ "7M", "6M1S", "1S6M" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:7" => 1 },
		{ "AS:i:-1" => 1, "XN:i:1" => 1, "YT:Z:UU" => 1, "MD:Z:6" => 1 },
		{ "AS:i:-1" => 1, "XN:i:1" => 1, "YT:Z:UU" => 1, "MD:Z:6" => 1 }
	  ]},

	# Same as previous case but --overhang not specified
	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TGTTCGT", "TTGTTCG", "GTTCGTA", "ATTGTTC" ],
	  args   => "-P \"SEED=0,2,1;NCEIL=2,0\"",
	  hits   => [ { 1 => 1 }, { 0 => 1 } ], # only the internal hits
	  cigar  => [ "7M", "7M", "*", "*" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:7" => 1 },
		{ "YT:Z:UU" => 1, "MD:Z:7" => 1 },
		{ "YT:Z:UU" => 1 },
		{ "YT:Z:UU" => 1 }
	  ]
	},

	# A simple case that should align with or without overhang, with or without
	# a special NCEIL setting.
	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TTGTTCG" ],
	  args   => "--overhang -P \"SEED=0,2,1;NCEIL=2,0\"",
	  hits   => [ { 0 => 1 } ]},

	{ ref    => [ "TTGTTCGT" ],
	  reads  => [ "TTGTTCG" ],
	  args   => "--overhang",
	  hits   => [ { 0 => 1 } ]},

	#
	# Testing the various -M/-m/-k/-a alignment modes in both unpaired and
	# paired-end modes.  Ensuring that SAM optional flags such as YM:i, YP:i
	# are set properly in all cases.
	#
	
	#
	# Paired-end
	#

	{ name   => "P.M.58.G.b Unpaired -M 5 w/ 8 hits global, but mate #1 has just 1",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2                                             0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               0123456789012345678901234567                                      0123456789012345678901234567                                                                                                                                               0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGT" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6                                                                                                   7                                                                                                   8                                                                                                   9                                                                                                   
	  #            0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  args   =>   "-M 5 -X 1000",
	  mapq_map   => [{ 12 => 40, 78 => 0, 249 => 0, 315 => 0, 486 => 0,
	                   552 => 0, 723 => 0, 789 => 0, 960 => 0 }],
	  pairhits   => [{ "12,78"  => 1, "12,249" => 1, "12,315" => 1,
	                   "12,486" => 1, "12,552" => 1, "12,723" => 1,
					   "12,789" => 1, "12,960" => 1 }],
	  hits_are_superset => [ 1 ],
	  cigar_map => [{
		12   => "33M",   78 => "28M",
		249  => "28M",  315 => "28M",
		486  => "28M",  552 => "28M",
		723  => "28M",  789 => "28M",
		960  => "28M"
	  }],
	  samoptflags_map => [ {
	    12 => {   "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:0" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    78 => {   "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    249 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    315 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    486 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    552 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    723 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    789 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    960 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	} ] },
	
	{ name   => "P.M.58.L.b Unpaired -M 5 w/ 8 hits local, but mate #1 has just 1",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2                                             0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               0123456789012345678901234567                                      0123456789012345678901234567                                                                                                                                               0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGT" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6                                                                                                   7                                                                                                   8                                                                                                   9                                                                                                   
	  #            0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  args   =>   "-M 5 -X 1000 --local",
	  mapq_map   => [{ 12 => 40, 78 => 0, 249 => 0, 315 => 0, 486 => 0,
	                   552 => 0, 723 => 0, 789 => 0, 960 => 0 }],
	  pairhits   => [{ "12,78"  => 1, "12,249" => 1, "12,315" => 1,
	                   "12,486" => 1, "12,552" => 1, "12,723" => 1,
					   "12,789" => 1, "12,960" => 1 }],
	  hits_are_superset => [ 1 ],
	  cigar_map => [{
		12   => "33M",   78 => "28M",
		249  => "28M",  315 => "28M",
		486  => "28M",  552 => "28M",
		723  => "28M",  789 => "28M",
		960  => "28M"
	  }],
	  samoptflags_map => [ {
	    12 => {   "AS:i:66" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:0" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:56"  => 1 },
	    78 => {   "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    249 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    315 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    486 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    552 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    723 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    789 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    960 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	} ] },

	{ name   => "P.m.58.G.b Unpaired -m 5 w/ 8 hits global, but mate #1 has just 1",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2                                             0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               0123456789012345678901234567                                      0123456789012345678901234567                                                                                                                                               0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGT" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6                                                                                                   7                                                                                                   8                                                                                                   9                                                                                                   
	  #            0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  args   =>   "-m 5 -X 1000",
	  mapq_map   => [{ 12 => 40, "*" => 0 }],
	  pairhits   => [{ "*,12" => 1 }],
	  cigar_map  => [{ 12 => "33M", "*" => "*" }],
	  samoptflags_map => [ {
	    12 => {   "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:0" => 1, "YP:i:1" => 1, "YT:Z:UP"  => 1, "YS:i:0"  => 1 },
	    "*" => {  "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:UP"  => 1 },
	} ] },
	
	{ name   => "P.m.58.L.b Unpaired -m 5 w/ 8 hits local, but mate #1 has just 1",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2                                             0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               0123456789012345678901234567                                      0123456789012345678901234567                                                                                                                                               0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGT" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6                                                                                                   7                                                                                                   8                                                                                                   9                                                                                                   
	  #            0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  args   =>   "-m 5 -X 1000 --local",
	  mapq_map   => [{ 12 => 40, "*" => 0 }],
	  pairhits   => [{ "*,12" => 1 }],
	  cigar_map  => [{ 12 => "33M", "*" => "*" }],
	  samoptflags_map => [ {
	    12 => {   "AS:i:66" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:0" => 1, "YP:i:1" => 1, "YT:Z:UP"  => 1, "YS:i:56"  => 1 },
	    "*" => {  "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:UP"  => 1 },
	} ] },

	{ name   => "P.k.58.G.b Unpaired -k 5 w/ 8 hits global, but mate #1 has just 1",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2                                             0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               0123456789012345678901234567                                      0123456789012345678901234567                                                                                                                                               0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGT" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6                                                                                                   7                                                                                                   8                                                                                                   9                                                                                                   
	  #            0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  args   =>   "-X 1000",
	  report =>   "-k 5",
	  mapq_map   => [{ 12 => 255, 78 => 0, 249 => 0, 315 => 0, 486 => 0,
	                   552 => 0, 723 => 0, 789 => 0, 960 => 0 }],
	  pairhits   => [{ "12,78"  => 1, "12,249" => 1, "12,315" => 1,
	                   "12,486" => 1, "12,552" => 1, "12,723" => 1,
					   "12,789" => 1, "12,960" => 1 }],
	  hits_are_superset => [ 1 ],
	  cigar_map => [{
		12   => "33M",   78 => "28M",
		249  => "28M",  315 => "28M",
		486  => "28M",  552 => "28M",
		723  => "28M",  789 => "28M",
		960  => "28M"
	  }],
	  samoptflags_map => [ {
	    12 => {   "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:0" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    78 => {   "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    249 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    315 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    486 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    552 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    723 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    789 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    960 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	} ] },
	
	{ name   => "P.k.58.L.b Unpaired -k 5 w/ 8 hits local, but mate #1 has just 1",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2                                             0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               0123456789012345678901234567                                      0123456789012345678901234567                                                                                                                                               0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG                                      ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACACACACCCCTATAGCTCGGAGCTGACTGGATCGACGACGT" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6                                                                                                   7                                                                                                   8                                                                                                   9                                                                                                   
	  #            0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  args   =>   "-X 1000 --local",
	  report =>   "-k 5",
	  mapq_map   => [{ 12 => 255, 78 => 0, 249 => 0, 315 => 0, 486 => 0,
	                   552 => 0, 723 => 0, 789 => 0, 960 => 0 }],
	  pairhits   => [{ "12,78"  => 1, "12,249" => 1, "12,315" => 1,
	                   "12,486" => 1, "12,552" => 1, "12,723" => 1,
					   "12,789" => 1, "12,960" => 1 }],
	  hits_are_superset => [ 1 ],
	  cigar_map => [{
		12   => "33M",   78 => "28M",
		249  => "28M",  315 => "28M",
		486  => "28M",  552 => "28M",
		723  => "28M",  789 => "28M",
		960  => "28M"
	  }],
	  samoptflags_map => [ {
	    12 => {   "AS:i:66" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:0" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:56"  => 1 },
	    78 => {   "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    249 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    315 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    486 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    552 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    723 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    789 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    960 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	} ] },
	
	{ name   => "P.M.22.G. Paired -M 2 w/ 2 paired hit, 2 unpaired hits each, global",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args   =>   "-M 2",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  mapq => [ 0 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		78  => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		249 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		315 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
	  }]
	},
	
	{ name   => "P.M.22.L. Paired -M 2 w/ 2 paired hit, 2 unpaired hits each, local",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args   =>   "-M 2 --local",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  mapq => [ 0 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:66" => 1, "XS:i:66" => 1, "MD:Z:33" => 1,
		         "YM:i:0"  => 1, "YP:i:0"  => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		78  => { "AS:i:56" => 1, "XS:i:56" => 1, "MD:Z:28" => 1,
		         "YM:i:0"  => 1, "YP:i:0"  => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
		249 => { "AS:i:66" => 1, "XS:i:66" => 1, "MD:Z:33" => 1,
		         "YM:i:0"  => 1, "YP:i:0"  => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		315 => { "AS:i:56" => 1, "XS:i:56" => 1, "MD:Z:28" => 1,
		         "YM:i:0"  => 1, "YP:i:0"  => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
	  }]
	},

	{ name   => "P.k.2.G. Paired -k 1 w/ 2 paired hit, 2 unpaired hits each, global",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args    => "",
	  report  => "-k 1",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  hits_are_superset => [ 1 ],
	  mapq => [ 255 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:0"  => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YT:Z:CP" => 1, "YS:i:0" => 1 },
		78  => { "AS:i:0"  => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YT:Z:CP" => 1, "YS:i:0" => 1 },
		249 => { "AS:i:0"  => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YT:Z:CP" => 1, "YS:i:0" => 1 },
		315 => { "AS:i:0"  => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YT:Z:CP" => 1, "YS:i:0" => 1 },
	  }]
	},

	{ name   => "P.k.2.L. Paired -k 1 w/ 2 paired hit, 2 unpaired hits each, local",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args    => "--local",
	  report  => "-k 1",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  hits_are_superset => [ 1 ],
	  mapq => [ 255 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:66" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YT:Z:CP" => 1, "YS:i:56" => 1 },
		78  => { "AS:i:56" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YT:Z:CP" => 1, "YS:i:66" => 1 },
		249 => { "AS:i:66" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YT:Z:CP" => 1, "YS:i:56" => 1 },
		315 => { "AS:i:56" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YT:Z:CP" => 1, "YS:i:66" => 1 },
	  }]
	},
	
	{ name   => "P.m.2.G. Paired -m 1 w/ 2 paired hit, 2 unpaired hits each, global",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s   => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s   => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw  => 1, mate2fw => 0,
	  args     =>   "-m 1",
	  pairhits => [ { "*,*" => 1 } ],
	  hits_are_superset => [ 1 ],
	  mapq  => [ 0 ],
	  cigar => [ "*" ],
	  samoptflags => [{ "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:UP" => 1 }],
	},

	{ name   => "P.m.2.L. Paired -m 1 w/ 2 paired hit, 2 unpaired hits each, local",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s   => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s   => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw  => 1, mate2fw => 0,
	  args     =>   "-m 1 --local",
	  pairhits => [ { "*,*" => 1 } ],
	  hits_are_superset => [ 1 ],
	  mapq  => [ 0 ],
	  cigar => [ "*" ],
	  samoptflags => [{ "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:UP" => 1 }],
	},
	
	{ name   => "P.M.2.G. Paired -M 1 w/ 2 paired hit, 2 unpaired hits each, global",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args   =>   "-M 1",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  hits_are_superset => [ 1 ],
	  mapq => [ 0 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		78  => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		249 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		315 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
	  }]
	},
	
	{ name   => "P.M.2.L. Paired -M 1 w/ 2 paired hit, 2 unpaired hits each, local",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args   =>   "-M 1 --local",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  hits_are_superset => [ 1 ],
	  mapq => [ 0 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:66" => 1, "XS:i:66" => 1, "MD:Z:33" => 1,
		         "YM:i:1"  => 1, "YP:i:1"  => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		78  => { "AS:i:56" => 1, "XS:i:56" => 1, "MD:Z:28" => 1,
		         "YM:i:1"  => 1, "YP:i:1"  => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
		249 => { "AS:i:66" => 1, "XS:i:66" => 1, "MD:Z:33" => 1,
		         "YM:i:1"  => 1, "YP:i:1"  => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		315 => { "AS:i:56" => 1, "XS:i:56" => 1, "MD:Z:28" => 1,
		         "YM:i:1"  => 1, "YP:i:1"  => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
	  }]
	},
	
	{ name   => "P.k.1.G. Paired -k w/ 1 paired hit, 1 unpaired hit each, global",
	  #                        0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1,  mate2fw => 0,
	  report  =>   "-k 1",
	  pairhits  => [ { "12,78" => 1 } ],
	  mapq => [ 255 ],
	  cigar_map => [{
		12 => "33M",
		78 => "28M"
	  }],
	  samoptflags_map => [{
		12 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		        "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		78 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		        "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
	  }]
	},
	
	{ name   => "P.k.1.L. Paired -k 1 w/ 1 paired hit, 1 unpaired hit each, local",
	  #                        0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1,  mate2fw => 0,
	  args    => "--local",
	  report  => "-k 1",
	  pairhits  => [ { "12,78" => 1 } ],
	  mapq => [ 255 ],
	  cigar_map => [{
		12 => "33M",
		78 => "28M"
	  }],
	  samoptflags_map => [{
		12 => { "AS:i:66" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		        "YM:i:0"  => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		78 => { "AS:i:56" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		        "YM:i:0"  => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
	  }]
	},
	
	{ name   => "P.m.1.G. Paired -m w/ 1 paired hit, 1 unpaired hit each, global",
	  #                        0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1,  mate2fw => 0,
	  args   =>   "-m 1",
	  pairhits  => [ { "12,78" => 1 } ],
	  mapq => [ 40 ],
	  cigar_map => [{
		12 => "33M",
		78 => "28M"
	  }],
	  samoptflags_map => [{
		12 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		        "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		78 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		        "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
	  }]
	},
	
	{ name   => "P.m.1.L. Paired -m w/ 1 paired hit, 1 unpaired hit each, local",
	  #                        0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1,  mate2fw => 0,
	  args   =>   "-m 1 --local",
	  pairhits  => [ { "12,78" => 1 } ],
	  mapq => [ 40 ],
	  cigar_map => [{
		12 => "33M",
		78 => "28M"
	  }],
	  samoptflags_map => [{
		12 => { "AS:i:66" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		        "YM:i:0"  => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		78 => { "AS:i:56" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		        "YM:i:0"  => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
	  }]
	},
	
	{ name   => "P.M.1.G. Paired -M w/ 1 paired hit, 1 unpaired hit each, global",
	  #                        0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1,  mate2fw => 0,
	  args   =>   "-M 1",
	  pairhits  => [ { "12,78" => 1 } ],
	  mapq => [ 40 ],
	  cigar_map => [{
		12 => "33M",
		78 => "28M"
	  }],
	  samoptflags_map => [{
		12 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		        "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		78 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		        "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
	  }]
	},
	
	{ name   => "P.M.1.L. Paired -M w/ 1 paired hit, 1 unpaired hit each, local",
	  #                          0         1         2         3                                   0         1         2       
	  #                          012345678901234567890123456789012                                 0123456789012345678901234567
	  #                          CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref      => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #              012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
	  #              0         1         2         3         4         5         6         7         8         9         0         1         2
	  mate1s   => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s   => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw  => 1,  mate2fw => 0,
	  args     =>   "-M 1 --local",
	  pairhits => [ { "12,78" => 1 } ],
	  mapq => [ 40 ],
	  cigar_map => [{
		12 => "33M",
		78 => "28M"
	  }],
	  samoptflags_map => [{
		12 => { "AS:i:66" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		        "YM:i:0"  => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		78 => { "AS:i:56" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		        "YM:i:0"  => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
	  }]
	},





	{ name   => "P.M.58.G. Unpaired -M 5 w/ 8 hits global",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                   
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0      
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6                                                                                                   7                                                                                                   8                                                                                                   9                                                                                                   0                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6                                                                                                   7                                                                                                   8                                                                                                   9                                                                                                   0                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6
	  #            0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   2                                                                                                   2                                                                                                   2                                                                                                   2                                                                                                   2                                                                                                   2                                                                                                   2
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  args   =>   "-M 5",
	  mapq   => [ 0 ],
	  pairhits  => [ { "12,78"     => 1, "249,315"   => 1, "486,552"   => 1,
	                   "723,789"   => 1, "960,1026"  => 1, "1197,1263" => 1,
					   "1434,1500" => 1, "1671,1737" => 1 } ],
	  hits_are_superset => [ 1 ],
	  cigar_map => [{
		12   => "33M",   78 => "28M",
		249  => "33M",  315 => "28M",
		486  => "33M",  552 => "28M",
		723  => "33M",  789 => "28M",
		960  => "33M", 1026 => "28M",
		1197 => "33M", 1263 => "28M",
		1434 => "33M", 1500 => "28M",
		1671 => "33M", 1737 => "28M"
	  }],
	  samoptflags_map => [ {
	    12 => {   "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    78 => {   "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    249 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    315 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    486 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    552 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    723 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    789 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    960 => {  "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    1026 => { "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    1197 => { "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    1263 => { "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    1434 => { "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    1500 => { "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    1671 => { "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	    1737 => { "AS:i:0" => 1, "XS:i:0" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:0"  => 1 },
	} ] },
	
	{ name   => "P.M.58.L. Unpaired -M 5 w/ 8 hits local",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                   
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0      
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6                                                                                                   7                                                                                                   8                                                                                                   9                                                                                                   0                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6                                                                                                   7                                                                                                   8                                                                                                   9                                                                                                   0                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6
	  #            0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   0                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   1                                                                                                   2                                                                                                   2                                                                                                   2                                                                                                   2                                                                                                   2                                                                                                   2                                                                                                   2
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  args   =>   "-M 5 --local",
	  mapq   => [ 0 ],
	  pairhits  => [ { "12,78"     => 1, "249,315"   => 1, "486,552"   => 1,
	                   "723,789"   => 1, "960,1026"  => 1, "1197,1263" => 1,
					   "1434,1500" => 1, "1671,1737" => 1 } ],
	  hits_are_superset => [ 1 ],
	  cigar_map => [{
		12   => "33M",   78 => "28M",
		249  => "33M",  315 => "28M",
		486  => "33M",  552 => "28M",
		723  => "33M",  789 => "28M",
		960  => "33M", 1026 => "28M",
		1197 => "33M", 1263 => "28M",
		1434 => "33M", 1500 => "28M",
		1671 => "33M", 1737 => "28M"
	  }],
	  samoptflags_map => [ {
	    12 => {   "AS:i:66" => 1, "XS:i:66" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:56"  => 1 },
	    78 => {   "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    249 => {  "AS:i:66" => 1, "XS:i:66" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:56"  => 1 },
	    315 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    486 => {  "AS:i:66" => 1, "XS:i:66" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:56"  => 1 },
	    552 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    723 => {  "AS:i:66" => 1, "XS:i:66" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:56"  => 1 },
	    789 => {  "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    960 => {  "AS:i:66" => 1, "XS:i:66" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:56"  => 1 },
	    1026 => { "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    1197 => { "AS:i:66" => 1, "XS:i:66" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:56"  => 1 },
	    1263 => { "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    1434 => { "AS:i:66" => 1, "XS:i:66" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:56"  => 1 },
	    1500 => { "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	    1671 => { "AS:i:66" => 1, "XS:i:66" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:33" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:56"  => 1 },
	    1737 => { "AS:i:56" => 1, "XS:i:56" => 1, "XN:i:0"   => 1, "XM:i:0"  => 1,
		          "XO:i:0" => 1, "XG:i:0" => 1, "NM:i:0"   => 1, "MD:Z:28" => 1,
		          "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP"  => 1, "YS:i:66"  => 1 },
	} ] },

	{ name   => "P.m.22.G. Paired -m 2 w/ 2 paired hit, 2 unpaired hits each, global",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args   =>   "-m 2",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  mapq => [ 0 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		78  => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		249 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		315 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
	  }]
	},
	
	{ name   => "P.m.22.L. Paired -m 2 w/ 2 paired hit, 2 unpaired hits each, local",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args   =>   "-m 2 --local",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  mapq => [ 0 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:66" => 1, "XS:i:66" => 1, "MD:Z:33" => 1,
		         "YM:i:0"  => 1, "YP:i:0"  => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		78  => { "AS:i:56" => 1, "XS:i:56" => 1, "MD:Z:28" => 1,
		         "YM:i:0"  => 1, "YP:i:0"  => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
		249 => { "AS:i:66" => 1, "XS:i:66" => 1, "MD:Z:33" => 1,
		         "YM:i:0"  => 1, "YP:i:0"  => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		315 => { "AS:i:56" => 1, "XS:i:56" => 1, "MD:Z:28" => 1,
		         "YM:i:0"  => 1, "YP:i:0"  => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
	  }]
	},

	{ name   => "P.M.22.G. Paired -M 2 w/ 2 paired hit, 2 unpaired hits each, global",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args   =>   "-M 2",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  mapq => [ 0 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		78  => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		249 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		315 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
	  }]
	},
	
	{ name   => "P.M.22.L. Paired -M 2 w/ 2 paired hit, 2 unpaired hits each, local",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args   =>   "-M 2 --local",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  mapq => [ 0 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:66" => 1, "XS:i:66" => 1, "MD:Z:33" => 1,
		         "YM:i:0"  => 1, "YP:i:0"  => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		78  => { "AS:i:56" => 1, "XS:i:56" => 1, "MD:Z:28" => 1,
		         "YM:i:0"  => 1, "YP:i:0"  => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
		249 => { "AS:i:66" => 1, "XS:i:66" => 1, "MD:Z:33" => 1,
		         "YM:i:0"  => 1, "YP:i:0"  => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		315 => { "AS:i:56" => 1, "XS:i:56" => 1, "MD:Z:28" => 1,
		         "YM:i:0"  => 1, "YP:i:0"  => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
	  }]
	},

	{ name   => "P.k.2.G. Paired -k 1 w/ 2 paired hit, 2 unpaired hits each, global",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args    => "",
	  report  => "-k 1",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  hits_are_superset => [ 1 ],
	  mapq => [ 255 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:0"  => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YT:Z:CP" => 1, "YS:i:0" => 1 },
		78  => { "AS:i:0"  => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YT:Z:CP" => 1, "YS:i:0" => 1 },
		249 => { "AS:i:0"  => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YT:Z:CP" => 1, "YS:i:0" => 1 },
		315 => { "AS:i:0"  => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YT:Z:CP" => 1, "YS:i:0" => 1 },
	  }]
	},

	{ name   => "P.k.2.L. Paired -k 1 w/ 2 paired hit, 2 unpaired hits each, local",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args    => "--local",
	  report  => "-k 1",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  hits_are_superset => [ 1 ],
	  mapq => [ 255 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:66" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YT:Z:CP" => 1, "YS:i:56" => 1 },
		78  => { "AS:i:56" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YT:Z:CP" => 1, "YS:i:66" => 1 },
		249 => { "AS:i:66" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YT:Z:CP" => 1, "YS:i:56" => 1 },
		315 => { "AS:i:56" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YT:Z:CP" => 1, "YS:i:66" => 1 },
	  }]
	},
	
	{ name   => "P.m.2.G. Paired -m 1 w/ 2 paired hit, 2 unpaired hits each, global",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s   => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s   => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw  => 1, mate2fw => 0,
	  args     =>   "-m 1",
	  pairhits => [ { "*,*" => 1 } ],
	  hits_are_superset => [ 1 ],
	  mapq  => [ 0 ],
	  cigar => [ "*" ],
	  samoptflags => [{ "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:UP" => 1 }],
	},

	{ name   => "P.m.2.L. Paired -m 1 w/ 2 paired hit, 2 unpaired hits each, local",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s   => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s   => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw  => 1, mate2fw => 0,
	  args     =>   "-m 1 --local",
	  pairhits => [ { "*,*" => 1 } ],
	  hits_are_superset => [ 1 ],
	  mapq  => [ 0 ],
	  cigar => [ "*" ],
	  samoptflags => [{ "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:UP" => 1 }],
	},
	
	{ name   => "P.M.2.G. Paired -M 1 w/ 2 paired hit, 2 unpaired hits each, global",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args   =>   "-M 1",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  hits_are_superset => [ 1 ],
	  mapq => [ 0 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		78  => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		249 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		         "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		315 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		         "YM:i:1" => 1, "YP:i:1" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
	  }]
	},
	
	{ name   => "P.M.2.L. Paired -M 1 w/ 2 paired hit, 2 unpaired hits each, local",
	  #                        0         1         2         3                                   0         1         2                                                                                                                                                      0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567                                                                                                                                               012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG                                                                                                                                               CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGATGCGTATCGACTATCGACAATATGACGCGTCGGTCACCCCATAATATGCAAAAATTATAGCTCACGACGCGTACTAATAGAAAACGCGCTATCAGCCTCCGACGCGGCGGTATCGAAGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1, mate2fw => 0,
	  args   =>   "-M 1 --local",
	  pairhits  => [ { "12,78" => 1, "249,315" => 1 } ],
	  hits_are_superset => [ 1 ],
	  mapq => [ 0 ],
	  cigar_map => [{
		12 => "33M", 249 => "33M",
		78 => "28M", 315 => "28M"
	  }],
	  samoptflags_map => [{
		12  => { "AS:i:66" => 1, "XS:i:66" => 1, "MD:Z:33" => 1,
		         "YM:i:1"  => 1, "YP:i:1"  => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		78  => { "AS:i:56" => 1, "XS:i:56" => 1, "MD:Z:28" => 1,
		         "YM:i:1"  => 1, "YP:i:1"  => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
		249 => { "AS:i:66" => 1, "XS:i:66" => 1, "MD:Z:33" => 1,
		         "YM:i:1"  => 1, "YP:i:1"  => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		315 => { "AS:i:56" => 1, "XS:i:56" => 1, "MD:Z:28" => 1,
		         "YM:i:1"  => 1, "YP:i:1"  => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
	  }]
	},
	
	{ name   => "P.k.1.G. Paired -k w/ 1 paired hit, 1 unpaired hit each, global",
	  #                        0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1,  mate2fw => 0,
	  report  =>   "-k 1",
	  pairhits  => [ { "12,78" => 1 } ],
	  mapq => [ 255 ],
	  cigar_map => [{
		12 => "33M",
		78 => "28M"
	  }],
	  samoptflags_map => [{
		12 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		        "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		78 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		        "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
	  }]
	},
	
	{ name   => "P.k.1.L. Paired -k 1 w/ 1 paired hit, 1 unpaired hit each, local",
	  #                        0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1,  mate2fw => 0,
	  args    => "--local",
	  report  => "-k 1",
	  pairhits  => [ { "12,78" => 1 } ],
	  mapq => [ 255 ],
	  cigar_map => [{
		12 => "33M",
		78 => "28M"
	  }],
	  samoptflags_map => [{
		12 => { "AS:i:66" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		        "YM:i:0"  => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		78 => { "AS:i:56" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		        "YM:i:0"  => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
	  }]
	},
	
	{ name   => "P.m.1.G. Paired -m w/ 1 paired hit, 1 unpaired hit each, global",
	  #                        0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1,  mate2fw => 0,
	  args   =>   "-m 1",
	  pairhits  => [ { "12,78" => 1 } ],
	  mapq => [ 40 ],
	  cigar_map => [{
		12 => "33M",
		78 => "28M"
	  }],
	  samoptflags_map => [{
		12 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		        "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		78 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		        "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
	  }]
	},
	
	{ name   => "P.m.1.L. Paired -m w/ 1 paired hit, 1 unpaired hit each, local",
	  #                        0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1,  mate2fw => 0,
	  args   =>   "-m 1 --local",
	  pairhits  => [ { "12,78" => 1 } ],
	  mapq => [ 40 ],
	  cigar_map => [{
		12 => "33M",
		78 => "28M"
	  }],
	  samoptflags_map => [{
		12 => { "AS:i:66" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		        "YM:i:0"  => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		78 => { "AS:i:56" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		        "YM:i:0"  => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
	  }]
	},
	
	{ name   => "P.M.1.G. Paired -M w/ 1 paired hit, 1 unpaired hit each, global",
	  #                        0         1         2         3                                   0         1         2       
	  #                        012345678901234567890123456789012                                 0123456789012345678901234567
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #            012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2
	  mate1s => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw => 1,  mate2fw => 0,
	  args   =>   "-M 1",
	  pairhits  => [ { "12,78" => 1 } ],
	  mapq => [ 40 ],
	  cigar_map => [{
		12 => "33M",
		78 => "28M"
	  }],
	  samoptflags_map => [{
		12 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		        "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
		78 => { "AS:i:0" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		        "YM:i:0" => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:0" => 1 },
	  }]
	},
	
	{ name   => "P.M.1.L. Paired -M w/ 1 paired hit, 1 unpaired hit each, local",
	  #                          0         1         2         3                                   0         1         2       
	  #                          012345678901234567890123456789012                                 0123456789012345678901234567
	  #                          CAGCGTACGGTATCTAGCTATGGGCATCGATCG                                 ACACACCCCTATAGCTCGGAGCTGACTG
	  ref      => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGAGCGGTATCTACAGCCACTCATCACACACCCCTATAGCTCGGAGCTGACTGGGTTACTGGGGGGGTATCGA" ],
	  #              012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
	  #              0         1         2         3         4         5         6         7         8         9         0         1         2
	  mate1s   => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  mate2s   => [ "CAGTCAGCTCCGAGCTATAGGGGTGTGT" ], # rev comped
	  mate1fw  => 1,  mate2fw => 0,
	  args     =>   "-M 1 --local",
	  pairhits => [ { "12,78" => 1 } ],
	  mapq => [ 40 ],
	  cigar_map => [{
		12 => "33M",
		78 => "28M"
	  }],
	  samoptflags_map => [{
		12 => { "AS:i:66" => 1, "XS:i:0" => 1, "MD:Z:33" => 1,
		        "YM:i:0"  => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:56" => 1 },
		78 => { "AS:i:56" => 1, "XS:i:0" => 1, "MD:Z:28" => 1,
		        "YM:i:0"  => 1, "YP:i:0" => 1, "YT:Z:CP" => 1, "YS:i:66" => 1 },
	  }]
	},

	#
	# Unpaired
	#

	{ name   => "U.M.1.G. Unpaired -M w/ 1 hit global",
	  #                        0         1         2         3
	  #                        012345678901234567890123456789012
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGA" ],
	  #            01234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5
	  reads  => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  args   =>   "-M 1",
	  report =>   "-a",
	  hits   => [ { 12 => 1 } ],
	  mapq   => [ 40 ],
	  cigar  => [ "33M" ],
	  samoptflags => [{
		"AS:i:0"   => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:33"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	}]
	},
	
	{ name   => "U.M.1.L. Unpaired -M w/ 1 hit local",
	  #                        0         1         2         3
	  #                        012345678901234567890123456789012
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGA" ],
	  #            01234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5
	  reads  => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  args   =>   "-M 1 --local",
	  report =>   "-a",
	  hits   => [ { 12 => 1 } ],
	  mapq   => [ 40 ],
	  cigar  => [ "33M" ],
	  samoptflags => [{
		"AS:i:66"  => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:33"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	}]
	},

	{ name   => "U.m.1.G. Unpaired -m w/ 1 hit global",
	  #                        0         1         2         3
	  #                        012345678901234567890123456789012
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGA" ],
	  #            01234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5
	  reads  => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  args   =>   "-m 1",
	  report =>   "-a",
	  mapq   => [ 40 ],
	  hits   => [ { 12 => 1 } ],
	  cigar  => [ "33M" ],
	  samoptflags => [ {
		"AS:i:0"   => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:33"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },
	
	{ name   => "U.m.1.L. Unpaired -m w/ 1 hit local",
	  #                        0         1         2         3
	  #                        012345678901234567890123456789012
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGA" ],
	  #            01234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5
	  reads  => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  args   =>   "-m 1 --local",
	  report =>   "-a",
	  mapq   => [ 40 ],
	  hits   => [ { 12 => 1 } ],
	  cigar  => [ "33M" ],
	  samoptflags => [ {
		"AS:i:66"  => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:33"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },

	{ name   => "U.k.1.G. Unpaired -k 1 w/ 1 hit global",
	  #                        0         1         2         3
	  #                        012345678901234567890123456789012
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGA" ],
	  #            01234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5
	  reads  => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  report =>   "-k 1",
	  mapq   => [ 255 ],
	  hits   => [ { 12 => 1 } ],
	  cigar  => [ "33M" ],
	  samoptflags => [ {
		"AS:i:0"   => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:33"  => 1, # mismatching positions/bases
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },
	
	{ name   => "U.M.1.L. Unpaired -m w/ 1 hit local",
	  #                        0         1         2         3
	  #                        012345678901234567890123456789012
	  #                        CAGCGTACGGTATCTAGCTATGGGCATCGATCG
	  ref    => [ "AGACGCAGTCACCAGCGTACGGTATCTAGCTATGGGCATCGATCGACGACGTACGA" ],
	  #            01234567890123456789012345678901234567890123456789012345
	  #            0         1         2         3         4         5
	  reads  => [ "CAGCGTACGGTATCTAGCTATGGGCATCGATCG" ],
	  args   => "--local",
	  report => "-k 1",
	  mapq   => [ 255 ],
	  hits   => [ { 12 => 1 } ],
	  cigar  => [ "33M" ],
	  samoptflags => [ {
		"AS:i:66"  => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:33"  => 1, # mismatching positions/bases
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },
	
	{ name   => "U.M.2.G. Unpaired -M 1 w/ 2 hit global",
	  #                  0         1         2                     0         1         2         
	  #                  012345678901234567890123456789            012345678901234567890123456789
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234
	  #            0         1         2         3         4         5         6         7         8
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   =>   "-M 1",
	  report =>   "-a",
	  mapq   => [ 0 ],
	  hits   => [ { 6 => 1, 48 => 1 } ],
	  hits_are_superset => [ 1 ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:0"   => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:1"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },
	
	{ name   => "U.M.2.L. Unpaired -M 1 w/ 2 hit local",
	  #                  0         1         2                     0         1         2         
	  #                  012345678901234567890123456789            012345678901234567890123456789
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234
	  #            0         1         2         3         4         5         6         7         8
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   =>   "-M 1 --local",
	  report =>   "-a",
	  mapq   => [ 0 ],
	  hits   => [ { 6 => 1, 48 => 1 } ],
	  hits_are_superset => [ 1 ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:60"  => 1, # alignment score
		"XS:i:60"  => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:1"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },

	{ name   => "U.m.2.G. Unpaired -m 1 w/ 2 hit global",
	  #                  0         1         2                     0         1         2         
	  #                  012345678901234567890123456789            012345678901234567890123456789
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234
	  #            0         1         2         3         4         5         6         7         8
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   =>   "-m 1",
	  report =>   "-a",
	  mapq   => [ 0 ],
	  hits   => [ { "*" => 1 } ],
	  cigar  => [ "*" ],
	  samoptflags => [ {
		"YM:i:1"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },
	
	{ name   => "U.m.2.L. Unpaired -m 1 w/ 2 hit local",
	  #                  0         1         2                     0         1         2         
	  #                  012345678901234567890123456789            012345678901234567890123456789
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234
	  #            0         1         2         3         4         5         6         7         8
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   =>   "-m 1 --local",
	  report =>   "-a",
	  mapq   => [ 0 ],
	  hits   => [ { "*" => 1 } ],
	  cigar  => [ "*" ],
	  samoptflags => [ {
		"YM:i:1"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },

	{ name   => "U.k.2.G. Unpaired -k 1 w/ 2 hit global",
	  #                  0         1         2                     0         1         2         
	  #                  012345678901234567890123456789            012345678901234567890123456789
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234
	  #            0         1         2         3         4         5         6         7         8
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  report =>   "-k 1",
	  mapq   => [ 255 ],
	  hits   => [ { 6 => 1, 48 => 1 } ],
	  hits_are_superset => [ 1 ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:0"   => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },

	{ name   => "U.k.2.L. Unpaired -k 1 w/ 2 hit local",
	  #                  0         1         2                     0         1         2         
	  #                  012345678901234567890123456789            012345678901234567890123456789
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234
	  #            0         1         2         3         4         5         6         7         8
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  report =>   "-k 1",
	  mapq   => [ 255 ],
	  hits   => [ { 6 => 1, 48 => 1 } ],
	  hits_are_superset => [ 1 ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:0"   => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },

	{ name   => "U.M.22.G. Unpaired -M 2 w/ 2 hit global",
	  #                  0         1         2                     0         1         2         
	  #                  012345678901234567890123456789            012345678901234567890123456789
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234
	  #            0         1         2         3         4         5         6         7         8
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   =>   "-M 2",
	  mapq   => [ 0 ],
	  hits   => [ { 6 => 1, 48 => 1 } ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:0"   => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },
	
	{ name   => "U.M.22.L. Unpaired -M 2 w/ 2 hit local",
	  #                  0         1         2                     0         1         2         
	  #                  012345678901234567890123456789            012345678901234567890123456789
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234
	  #            0         1         2         3         4         5         6         7         8
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   =>   "-M 2 --local",
	  mapq   => [ 0 ],
	  hits   => [ { 6 => 1, 48 => 1 } ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:60"  => 1, # alignment score
		"XS:i:60"  => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },

	{ name   => "U.m.22.G. Unpaired -m 2 w/ 2 hit global",
	  #                  0         1         2                     0         1         2         
	  #                  012345678901234567890123456789            012345678901234567890123456789
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234
	  #            0         1         2         3         4         5         6         7         8
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   =>   "-m 2",
	  mapq   => [ 0 ],
	  hits   => [ { 6 => 1, 48 => 1 } ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:0"   => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },
	
	{ name   => "U.m.22.L. Unpaired -m 2 w/ 2 hit local",
	  #                  0         1         2                     0         1         2         
	  #                  012345678901234567890123456789            012345678901234567890123456789
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234
	  #            0         1         2         3         4         5         6         7         8
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   =>   "-m 2 --local",
	  mapq   => [ 0 ],
	  hits   => [ { 6 => 1, 48 => 1 } ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:60"  => 1, # alignment score
		"XS:i:60"  => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },

	{ name   => "U.k.22.G. Unpaired -k 2 w/ 2 hit global",
	  #                  0         1         2                     0         1         2         
	  #                  012345678901234567890123456789            012345678901234567890123456789
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234
	  #            0         1         2         3         4         5         6         7         8
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  report =>   "-k 2",
	  mapq   => [ 0 ],
	  hits   => [ { 6 => 1, 48 => 1 } ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:0"   => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },
	
	{ name   => "U.k.22.L. Unpaired -k 2 w/ 2 hit local",
	  #                  0         1         2                     0         1         2         
	  #                  012345678901234567890123456789            012345678901234567890123456789
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234
	  #            0         1         2         3         4         5         6         7         8
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   => "--local",
	  report => "-k 2",
	  mapq   => [ 0 ],
	  hits   => [ { 6 => 1, 48 => 1 } ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:60"  => 1, # alignment score
		"XS:i:60"  => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },

	{ name   => "U.M.58.G. Unpaired -M 5 w/ 8 hits global",
	  #                  0         1         2                     0         1         2                      0         1         2                     0         1         2                      0         1         2                     0         1         2                      0         1         2                     0         1         2                
	  #                  012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789       
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA       
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   =>   "-M 5",
	  mapq   => [ 0 ],
	  hits   => [ { 6 => 1, 48 => 1, 91 => 1, 133 => 1, 176 => 1, 218 => 1, 261 => 1, 303 => 1 } ],
	  hits_are_superset => [ 1 ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:0"   => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:1"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },
	
	{ name   => "U.M.58.L. Unpaired -M 5 w/ 8 hits global",
	  #                  0         1         2                     0         1         2                      0         1         2                     0         1         2                      0         1         2                     0         1         2                      0         1         2                     0         1         2                
	  #                  012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789       
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA       
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   =>   "-M 5 --local",
	  mapq   => [ 0 ],
	  hits   => [ { 6 => 1, 48 => 1, 91 => 1, 133 => 1, 176 => 1, 218 => 1, 261 => 1, 303 => 1 } ],
	  hits_are_superset => [ 1 ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:60"  => 1, # alignment score
		"XS:i:60"  => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:1"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },

	{ name   => "U.m.58.G. Unpaired -m 5 w/ 8 hits global",
	  #                  0         1         2                     0         1         2                      0         1         2                     0         1         2                      0         1         2                     0         1         2                      0         1         2                     0         1         2                
	  #                  012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789       
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA       
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   =>   "-m 5",
	  mapq   => [ 0 ],
	  hits   => [ { "*" => 1 } ],
	  cigar  => [ "*" ],
	  samoptflags => [ {
		"YM:i:1"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },
	
	{ name   => "U.m.58.L. Unpaired -m 5 w/ 8 hits global",
	  #                  0         1         2                     0         1         2                      0         1         2                     0         1         2                      0         1         2                     0         1         2                      0         1         2                     0         1         2                
	  #                  012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789       
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA       
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   =>   "-m 5 --local",
	  mapq   => [ 0 ],
	  hits   => [ { "*" => 1 } ],
	  cigar  => [ "*" ],
	  samoptflags => [ {
		"YM:i:1"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },

	{ name   => "U.k.58.G. Unpaired -k 5 w/ 8 hits global",
	  #                  0         1         2                     0         1         2                      0         1         2                     0         1         2                      0         1         2                     0         1         2                      0         1         2                     0         1         2                
	  #                  012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789       
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA       
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  report =>   "-k 5",
	  mapq   => [ 0 ],
	  hits   => [ { 6 => 1, 48 => 1, 91 => 1, 133 => 1, 176 => 1, 218 => 1, 261 => 1, 303 => 1 } ],
	  hits_are_superset => [ 1 ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:0"   => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },
	
	{ name   => "U.k.58.L. Unpaired -k 5 w/ 8 hits local",
	  #                  0         1         2                     0         1         2                      0         1         2                     0         1         2                      0         1         2                     0         1         2                      0         1         2                     0         1         2                
	  #                  012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789             012345678901234567890123456789            012345678901234567890123456789       
	  #                  AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA             AGATTACGGATCTACGATTCGAGTCGGTCA            AGATTACGGATCTACGATTCGAGTCGGTCA       
	  ref    => [ "AGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGAAGACGCAGATTACGGATCTACGATTCGAGTCGGTCAGTCACCAGCGTAAGATTACGGATCTACGATTCGAGTCGGTCAAGTGCGA" ],
	  #            0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	  #            0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         
	  #            0                                                                                                   1                                                                                                   2                                                                                                   3
	  reads  => [ "AGATTACGGATCTACGATTCGAGTCGGTCA" ],
	  args   => "--local",
	  report => "-k 5",
	  mapq   => [ 0 ],
	  hits   => [ { 6 => 1, 48 => 1, 91 => 1, 133 => 1, 176 => 1, 218 => 1, 261 => 1, 303 => 1 } ],
	  hits_are_superset => [ 1 ],
	  cigar  => [ "30M" ],
	  samoptflags => [ {
		"AS:i:60"  => 1, # alignment score
		"XS:i:60"  => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:30"  => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YT:Z:UU"  => 1, # unpaired read aligned in unpaired fashion
	} ] },

	#
	# Local alignment
	#

	# Local alignment for a short perfect hit where hit spans the read
	{ name   => "Local alignment 1",
	  ref    => [ "TTGT" ],
	  reads  => [ "TTGT" ],
	  args   =>   "--local -P \"MIN=0.0,0.75\"",
	  hits   => [ { 0 => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU,XC:4=" ],
	  cigar  => [ "4M" ],
	  samoptflags => [ {
		"AS:i:8"   => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:4"   => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YP:i:0"   => 1, # read aligned repetitively in paired fashion
		"YT:Z:UU"  => 1, # type of alignment (concordant/discordant/etc)
	} ] },
	
	#   T T G A     T T G A
	# T x         T   x
	# T   x       T      
	# G     x     G        
	# T           T
	
	# Local alignment for a short hit where hit is trimmed at one end
	{ name   => "Local alignment 2",
	  ref    => [ "TTGA" ],
	  reads  => [ "TTGT" ],
	  args   =>   "--local -P \"MIN=0.0,0.75;SEED=0,3,1\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ],
	  flags => [ "XM:0,XP:0,XT:UU,XC:3=1S" ],
	  cigar  => [ "3M1S" ],
	  samoptflags => [ {
		"AS:i:6"   => 1, # alignment score
		"XS:i:0"   => 1, # suboptimal alignment score
		"XN:i:0"   => 1, # num ambiguous ref bases
		"XM:i:0"   => 1, # num mismatches
		"XO:i:0"   => 1, # num gap opens
		"XG:i:0"   => 1, # num gap extensions
		"NM:i:0"   => 1, # num edits
		"MD:Z:3"   => 1, # mismatching positions/bases
		"YM:i:0"   => 1, # read aligned repetitively in unpaired fashion
		"YP:i:0"   => 1, # read aligned repetitively in paired fashion
		"YT:Z:UU"  => 1, # type of alignment (concordant/discordant/etc)
	} ] },

	#     0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
	#     T T G T T C G T T T G T T C G T
	# 0 T x
	# 1 T   x
	# 2 G     x
	# 3 T       x
	# 4 T         x
	# 5 C           x
	# 6 G             x
	# 7 T               x
	# 8 T                 x
	# 9 T                   x
	# 0 G                     x
	# 1 T                       x
	# 2 T                         x
	#
	# Score=130

	#     0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
	#     T T G T T C G T T T G T T C G T
	# 0 T                 x
	# 1 T                   x
	# 2 G                     x
	# 3 T                       x
	# 4 T                         x
	# 5 C                           x
	# 6 G                             x
	# 7 T                               x
	# 8 T
	# 9 T
	# 0 G
	# 1 T
	# 2 T
	#
	# Score=80

	# Local alignment for a perfect hit
	{ name   => "Local alignment 3",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGTTTGTT" ],
	  args   =>   "--local -P \"POSF=0,1.0\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1, 8 => 1 } ],
	  flags_map => [{
		0 => "XM:0,XP:0,XT:UU,XC:13=",
		8 => "XM:0,XP:0,XT:UU,XC:8="
	  }],
	  cigar_map => [{
		0 => "13M",
		8 => "8M5S"
	  }],
	  samoptflags_map => [{
		0 => { "AS:i:26" => 1, "XS:i:16" => 1, "YT:Z:UU" => 1, "MD:Z:13" => 1 },
		8 => { "AS:i:16" => 1, "XS:i:16" => 1, "YT:Z:UU" => 1, "MD:Z:8"  => 1 }
	  }]
	},

	#                          1
	#      0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
	#      T T G T T C G T T T G T T C G T
	#  0 T                 x
	#  1 T                   x
	#  2 G                     x
	#  3 T                       x
	#  4 T                         x
	#  5 C                           x
	#  6 G                             x
	#  7 T                               x
	#  8 T
	#  9 T
	# 10 G
	#  1 T

	# Local alignment for a hit that should be trimmed from the right end
	{ name   => "Local alignment 4",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGTTTGT" ],
	  args   =>   "--local -P \"SEED=0,3,1;POSF=0,1.0\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1, 8 => 1 } ],
	  flags_map => [{
		0 => "XM:0,XP:0,XT:UU,XC:12=",
		8 => "XM:0,XP:0,XT:UU,XC:8=4S"
	  }],
	  cigar_map => [{
		0 => "12M",
		8 => "8M4S"
	  }],
	  samoptflags_map => [{
		0 => { "AS:i:24" => 1, "XS:i:16" => 1, "YT:Z:UU" => 1, "MD:Z:12" => 1 },
		8 => { "AS:i:16" => 1, "XS:i:16" => 1, "YT:Z:UU" => 1, "MD:Z:8"  => 1 }
	  }]
	},

	# Following cases depend on this being the case:
	#
	# static const float DEFAULT_CEIL_CONST = 3.0f;
	# static const float DEFAULT_CEIL_LINEAR = 3.0f;

	# Just enough budget for hits, so it should align
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCAT" ], # budget = 3 + 8 * 3 = 27
	  args   => "-P \"SEED=2;MMP=C27;MIN=-3,-3;RDG=25,15;RFG=25,15\"", # penalty = 27
	  report => "-a",
	  hits   => [ { 0 => 1, 8 => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU,XC:6=1X1=" ],
	  cigar  => [ "8M" ],
	  samoptflags => [ { "AS:i:-27" => 1, "XS:i:-27" => 1, "NM:i:1" => 1,
	                     "XM:i:1" => 1, "YT:Z:UU" => 1, "MD:Z:6G1" => 1 } ] },

	# Not quite enough budget for hits, so it should NOT align
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCAT" ], # budget = 3 + 8 * 3 = 27
	  args   =>   "-P \"SEED=2;MMP=C28;MIN=-3,-3;RDG=25,15;RFG=25,15\"", # penalty = 28
	  report =>   "-a",
	  hits   => [ { "*" => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU" ],
	  cigar  => [ "*" ],
	  samoptflags => [ { "YT:Z:UU" => 1 } ] },

	# Check that using a seed of length 1 with 1-mismatch doesn't crash.
	# Perhaps we should disallow it though?

	{ ref    => [ "AAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCC" ],
	  reads  => [ "AA", "AA", "AA", "AA", "CC", "CC", "CC", "CC", "AA", "AA", "AA", "AA", "CC", "CC", "CC", "CC" ],
	  names  => [ "r1", "r1", "r1", "r1", "r2", "r2", "r2", "r2", "r3", "r3", "r3", "r3", "r4", "r4", "r4", "r4" ],
	  args   => "-P \"SEED=1,1\"",
	  check_random => 1,
	  report => "-k 1" },

	#
	# Gap penalties
	#

	# Alignment with 1 read gap
	{ name   => "Gap penalties 1",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCTTTGTT" ], # budget = 3 + 12 * 3 = 39
	  args   =>   "-P \"MMP=C30;SEED=0,3,1;RDG=29,10;RFG=25,15;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU,XC:6=1D6=" ],
	  cigar  => [ "6M1D6M" ],
	  samoptflags => [{
		"AS:i:-39" => 1, "NM:i:1" => 1, "XO:i:1" => 1, "XG:i:1" => 1,
		"YT:Z:UU" => 1, "MD:Z:6^G6" => 1 }]
	},

	# Alignment with 1 read gap - colorspace
	{ name   => "Gap penalties 1 (colorspace)",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "01102200110" ], # budget = 3 + 11 * 3 = 36
	  #           "TTGTTCTTTGTT"
	  args   =>   "-P \"MMP=C30;SEED=0,3,1;RDG=26,10;RFG=25,15;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { 1 => 1 } ],
	  color  => 1,
	  flags  => [ "XM:0,XP:0,XT:UU,XC:5=1D5=" ],
	  cigar  => [ "5M1D5M" ],
	  samoptflags => [{
		"AS:i:-36" => 1, "NM:i:1" => 1, "XO:i:1" => 1, "XG:i:1" => 1,
		"YT:Z:UU" => 1, "MD:Z:5^G5" => 1 }]
	},

	# Alignment with 1 read gap, but not enough budget
	{ name   => "Gap penalties 2",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCTTTGTT" ], # budget = 3 + 12 * 3 = 39
	  args   =>   "-P \"MMP=C30;SEED=0,3,1;RDG=30,10;RFG=25,15;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { "*" => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU" ],
	  cigar  => [ "*" ],
	  samoptflags => [{ "YT:Z:UU" => 1 }]
	},

	# Alignment with 1 read gap, but not enough budget - colorspace
	{ name   => "Gap penalties 2 (colorspace)",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "01102200110" ], # budget = 3 + 11 * 3 = 36
	  #           "TTGTTCTTTGTT"
	  args   =>   "-P \"MMP=C30;SEED=0,3,1;RDG=27,10;RFG=25,15;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { "*" => 1 } ],
	  color  => 1,
	  flags  => [ "XM:0,XP:0,XT:UU" ],
	  cigar  => [ "*" ],
	  samoptflags => [{ "YT:Z:UU" => 1 }]
	},

	# Alignment with 1 reference gap
	{ name   => "Gap penalties 3",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGATTTGTT" ], # budget = 3 + 14 * 3 = 45
	  args   =>   "-P \"MMP=C30;SEED=0,3,1;RDG=25,15;RFG=30,15;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU,XC:7=1I6=" ],
	  cigar  => [ "7M1I6M" ],
	  samoptflags => [{ "AS:i:-45" => 1, "NM:i:1" => 1, "XO:i:1" => 1,
	                    "XG:i:1" => 1, "YT:Z:UU" => 1, "MD:Z:13" => 1 }]
	},

	#      0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
	#      T T G T T C G T T T G T T C G T
	#       0 1 1 0 2 3 1 0 0 1 1 0 2 3 1
	# 0 T  x
	# 0  0  x
	# 1 T    x
	# 1  1    x
	# 2 G      x
	# 2  1      x
	# 3 T        x
	# 3  0        x
	# 4 T          x
	# 4  2          x
	# 5 C            x
	# 5  3            x
	# 6 G
	# 6  2            x
	# 7 A
	# 7  3              x
	# 8 T                x
	# 8  0                x
	# 9 T                  x
	# 9  0                  x
	# 0 T                    x
	# 0  1                    x
	# 1 G                      x
	# 1  1                      x
	# 2 T                        x
	# 2  0                        x
	# 3 T
	#

	# Alignment with 1 reference gap - colorspace
	{ name   => "Gap penalties 3 (colorspace)",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "0110232300110" ], # budget = 3 + 13 * 3 = 42
	  #           "TTGTTCGATTTGTT"
	  args   =>   "-P \"SEED=0,3,1;RFG=27,15\"",
	  report =>   "-a",
	  hits   => [ { 1 => 1 } ],
	  color  => 1,
	  flags  => [ "XM:0,XP:0,XT:UU,XC:6=1I5=" ],
	  cigar  => [ "6M1I5M" ],
	  samoptflags => [ { "AS:i:-42" => 1, "NM:i:1" => 1, "XO:i:1" => 1,
	                     "XG:i:1" => 1, "YT:Z:UU" => 1, "MD:Z:11" => 1} ] },

	# Alignment with 1 reference gap, but not enough budget
	{ name   => "Gap penalties 4",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGATTTGTT" ], # budget = 3 + 14 * 3 = 45
	  args   =>   "-P \"MMP=C30;SEED=0,3,1;RDG=25,15;RFG=30,16;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { "*" => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU" ],
	  cigar  => [ "*" ],
	  samoptflags => [ { "YT:Z:UU" => 1 } ] },

	# Alignment with 1 reference gap, but not enough budget - colorspace
	{ name   => "Gap penalties 4 (colorspace)",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "0110232300110" ], # budget = 3 + 13 * 3 = 42
	  #           "TTGTTCGATTTGTT"
	  args   =>   "-P \"MMP=C30;SEED=0,3,1;RDG=25,15;RFG=27,16;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { "*" => 1 } ],
	  color  => 1,
	  flags  => [ "XM:0,XP:0,XT:UU" ],
	  cigar  => [ "*" ],
	  samoptflags => [ { "YT:Z:UU" => 1 } ] },

	# Alignment with 1 reference gap, but not enough budget
	{ name   => "Gap penalties 5",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGATTTGTT" ], # budget = 3 + 14 * 3 = 45
	  args   =>   "-P \"MMP=C30;SEED=0,3,1;RDG=25,15;RFG=31,15;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { "*" => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU" ],
	  cigar  => [ "*" ],
	  samoptflags => [ { "YT:Z:UU" => 1 } ] },

	# Alignment with 1 reference gap, but not enough budget - colorspace
	{ name   => "Gap penalties 5 (colorspace)",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "0110232300110" ], # budget = 3 + 13 * 3 = 42
	  #           "TTGTTCGATTTGTT"
	  args   =>   "-P \"MMP=C30;SEED=0,3,1;RDG=25,15;RFG=28,15;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { "*" => 1 } ],
	  color  => 1,
	  flags  => [ "XM:0,XP:0,XT:UU" ],
	  cigar  => [ "*" ],
	  samoptflags => [ { "YT:Z:UU" => 1 } ] },

	# Alignment with 1 reference gap and 1 read gap
	{ name   => "Gap penalties 6",
	  ref    => [ "ATTGTTCGTTTGTTCGTA" ],
	  reads  => [ "ATTGTTGTTTGATTCGTA" ], # budget = 3 + 18 * 3 = 57
	  args   =>   "-P \"MMP=C30;SEED=0,3,1;RDG=19,10;RFG=18,10;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU,XC:6=1D5=1I6=" ],
	  cigar  => [ "6M1D5M1I6M" ] },

	# Alignment with 1 reference gap and 1 read gap - colorspace
	{ name   => "Gap penalties 6 (colorspace)",
	  ref    => [ "ATTGTTCGTTTGTTCGTA" ],
	  reads  => [ "30110110012302313" ], # budget = 3 + 17 * 3 = 54
	  #           "ATTGTTGTTTGATTCGTA"
	  args   =>   "-P \"MMP=C30;SEED=0,3,1;RDG=16,10;RFG=18,10;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { 1 => 1 } ],
	  color  => 1,
	  flags  => [ "XM:0,XP:0,XT:UU,XC:5=1D5=1I5=" ],
	  cigar  => [ "5M1D5M1I5M" ] },

	# Alignment with 1 reference gap and 1 read gap, but not enough budget
	{ name   => "Gap penalties 7",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTGTTTGATTCGT" ], # budget = 3 + 16 * 3 = 51
	  args   =>   "-P \"MMP=C30;SEED=0,3,1;RDG=16,10;RFG=16,10;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { "*" => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU" ],
	  cigar  => [ "*" ] },

	# Alignment with 1 reference gap and 1 read gap, but not enough budget -
	# colorspace
	{ name   => "Gap penalties 7 (colorspace)",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "011011001230231" ], # budget = 3 + 15 * 3 = 48
	  #           "TTGTTGTTTGATTCGT"
	  args   =>   "-P \"MMP=C30;SEED=0,3,1;RFG=16,10;RDG=13,10;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { "*" => 1 } ],
	  color  => 1,
	  flags  => [ "XM:0,XP:0,XT:UU" ],
	  cigar  => [ "*" ] },

	# Experiment with N filtering
	
	{ name => "N filtering 1",
	  ref      => [ "GAGACTTTATACGCATCGAACTATCGCTCTA" ],
	  reads    => [         "ATACGCATCGAAC" ],
	  #              0123456789012345678901234567890
	  #                        1         2         3
	  args     =>   "-P \"NCEIL=0,0\"",
	  report   =>   "-a",
	  hits     => [ { 8 => 1 } ],
	  flags => [ "XM:0,XP:0,XT:UU,XC:13=" ] },

	{ name => "N filtering 2",
	  ref      => [ "GAGACTTTATNCGCATCGAACTATCGCTCTA" ],
	  reads    => [         "ATACGCATCGAAC" ],
	  #              0123456789012345678901234567890
	  #                        1         2         3
	  args     =>   "-P \"NCEIL=0,0\"",
	  report   =>   "-a",
	  hits     => [ { "*" => 1 } ] },

	{ name => "N filtering 3",
	  ref      => [ "GAGACTTTATACGCATCGAANTATCGCTCTA" ],
	  reads    => [         "ATACGCATCGAAC" ],
	  #              0123456789012345678901234567890
	  #                        1         2         3
	  args     =>   "-P \"NCEIL=0,0\"",
	  report   =>   "-a",
	  hits     => [ { "*" => 1 } ] },

	{ name => "N filtering 4",
	  ref      => [ "GAGACTTTNTACGCATCGAACTATCGCTCTA" ],
	  reads    => [         "ATACGCATCGAAC" ],
	  #              0123456789012345678901234567890
	  #                        1         2         3
	  args     =>   "-P \"NCEIL=0,0\"",
	  report   =>   "-a",
	  hits     => [ { "*" => 1 } ] },

	{ name => "N filtering 5",
	  ref      => [ "GAGACTTTATNCGCATCGAACTATCGCTCTA" ],
	  reads    => [         "ATACGCATCGAAC" ],
	  #              0123456789012345678901234567890
	  #                        1         2         3
	  args     =>   "-P \"NCEIL=0,0.1;SEED=0,10,1\"",
	  report   =>   "-a",
	  hits     => [ { 8 => 1 } ],
	  flags => [ "XM:0,XP:0,XT:UU,XC:2=1X10=" ] },

	{ name => "N filtering 6",
	  ref      => [ "GAGACTTTNTACGCATCGAANTATCGCTCTA" ],
	  reads    => [         "ATACGCATCGAAC" ],
	  #              0123456789012345678901234567890
	  #                        1         2         3
	  args     =>   "-P \"NCEIL=0,0.1;SEED=0,10,1\"",
	  report   =>   "-a",
	  hits     => [ { "*" => 1 } ] },

	# No discordant alignment because one mate is repetitive.

	{ name => "Simple paired-end 15",
	  ref    => [ "TTTATAAAAATATTTTTTATAAAAATATTTTCCCCCCCCCCCGCCGGCGCGCCCAACGCCGGCGCGCCCCC" ],
	#                 ATAAAAATAT     ATAAAAATAT             CGCCGGCGCG     CGCCGGCGCG
	#                                                     CGCGCCGGCG     CGCGCCGGCG
	#              01234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7
	#                 -------------------------------------------
	#                 0123456789012345678901234567890123456789012
	  mate1s   => [ "ATAAAAATAT" ],
	  mate2s   => [ "CGCCGGCGCG" ],
	  mate1fw  => 1, mate2fw => 1,
	  args     => "-I 0 -X 70 -m 2 -P \"MMP=C30\"",
	  mapq     => [ 0 ],
	  pairhits => [ { 3 => 1, 18 => 1, 41 => 1, 56 => 1 } ],
	  flags    => [ "XM:0,XP:1,XT:UP,XC:10=" ], # Pair aligns repetitively
	  samoptflags_map => [{
		  3  => { "AS:i:0" => 1, "YT:Z:UP" => 1, "MD:Z:10" => 1, "YP:i:1" => 1, "YM:i:0" => 1 },
		  18 => { "AS:i:0" => 1, "YT:Z:UP" => 1, "MD:Z:10" => 1, "YP:i:1" => 1, "YM:i:0" => 1 },
		  41 => { "AS:i:0" => 1, "YT:Z:UP" => 1, "MD:Z:10" => 1, "YP:i:1" => 1, "YM:i:0" => 1 },
		  56 => { "AS:i:0" => 1, "YT:Z:UP" => 1, "MD:Z:10" => 1, "YP:i:1" => 1, "YM:i:0" => 1 } }]
	},

	# No discordant alignment because one mate is repetitive.

	{ name => "Simple paired-end 14",
	  ref    => [ "TTTATAAAAATATTTCCCCCCCGCCCGCCCGCCCCCCGCCCGCCCGCCCCC" ],
	#                 ATAAAAATAT        CGCCCGCCCG     CGCCCGCCCG
	#              012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5
	#                 -------------------------------------------
	#                 0123456789012345678901234567890123456789012
	  mate1s   => [ "ATAAAAATAT" ],
	  mate2s   => [ "CGCCCGCCCG" ],
	  mate1fw  => 1,  mate2fw => 1,
	  args     => "-I 0 -X 50 -m 1",
	  mapq_map => [{ "*" => 0, 3 => 40 }],
	  pairhits => [ { "*,3" => 1 } ],
	  flags_map  => [{
		3   => "XM:0,XP:1,XT:UP,XC:10=",
		"*" => "XM:1,XP:1,XT:UP" } ],
	  samoptflags_map => [{
		3   => { "AS:i:0" => 1, "YT:Z:UP" => 1, "MD:Z:10" => 1, "YP:i:1" => 1, "YM:i:0" => 1 },
		"*" => { "YT:Z:UP" => 1, "YP:i:1" => 1, "YM:i:1" => 1 } }]
	},

	# Alignment with 1 reference gap
	{ ref    => [ "TTTTGTTCGTTTG" ],
	  reads  => [ "TTTTGTTCGATTTG" ], # budget = 3 + 14 * 3 = 45
	  args   =>   "-P \"SEED=0,8,1;MMP=C30;RDG=25,15;RFG=25,20;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ],
	  flags => [ "XM:0,XP:0,XT:UU,XC:9=1I4=" ],
	  cigar  => [ "9M1I4M" ],
	  samoptflags => [
		{ "AS:i:-45" => 1, "NM:i:1"  => 1, "XO:i:1" => 1, "XG:i:1" => 1,
		  "YT:Z:UU" => 1, "MD:Z:13" => 1 },
	  ]
	},

	#  TTGTTCGTTTGTT
	# Tx
	# T x
	# G  x
	# T   x
	# T    x
	# C     x
	# G      x
	# A      x
	# T       x
	# T        x
	# T         x
	# G          x
	# T           x
	# T            x
	
	# Alignment with 1 reference gap
	{ ref    => [ "TTGTTCGTTTGTT" ],
	  reads  => [ "TTGTTCGATTTGTT" ], # budget = 3 + 14 * 3 = 45
	  args   =>   "-P \"SEED=0,3,1;MMP=C30;RDG=25,15;RFG=25,20;MIN=-3,-3\"",
	  report =>   "-a",
	  hits   => [ { 0 => 1 } ],
	  flags => [ "XM:0,XP:0,XT:UU,XC:7=1I6=" ],
	  cigar  => [ "7M1I6M" ],
	  samoptflags => [
		{ "AS:i:-45" => 1, "NM:i:1"  => 1, "XO:i:1" => 1, "XG:i:1" => 1,
		  "YT:Z:UU" => 1, "MD:Z:13" => 1 },
	  ]
	},

	{ ref    => [ "ACNCA" ],
	  reads  => [ "CA" ],
	  args   => "",
	  report => "-a -P \"SEED=0,2,1;NCEIL=0,0\"",
	  hits   => [ { 3 => 1 } ],
	  edits  => [ ],
	  flags => [ "XM:0,XP:0,XT:UU,XC:2=" ],
	  cigar  => [ "2M" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:2" => 1 },
	  ]
	},

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
	  flags => [ "XM:0,XP:0,XT:UU,XC:2=", "XM:0,XP:0,XT:UU,XC:2=" ],
	  cigar  => [ "2M", "2M" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:2" => 1 },
		{ "YT:Z:UU" => 1, "MD:Z:2" => 1 },
	  ]
	},

	{ ref    => [ "ACNCANNNNNNAACGNNNNNNNACGAANNNNCGAAAN" ],
	#              0123456789012345678901234567890123456
	#              0         1         2         3
	  reads  => [ "CG" ],
	  args   => "",
	  report => "-a -P \"SEED=0,2,1;NCEIL=0,0;ROWM=3\"",
	  hits   => [ { 13 => 2, 23 => 2, 31 => 2 } ],
	  edits  => [ ],
	  flags => [ "XM:0,XP:0,XT:UU,XC:2=",
	             "XM:0,XP:0,XT:UU,XC:2=",
				 "XM:0,XP:0,XT:UU,XC:2=" ],
	  cigar  => [ "2M", "2M", "2M" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:2" => 1 },
		{ "YT:Z:UU" => 1, "MD:Z:2" => 1 },
		{ "YT:Z:UU" => 1, "MD:Z:2" => 1 },
	  ]
	},

	{ ref    => [ "ACNCANNNNNNAACGNNNNNNNACGAANNNNCGAAAN" ],
	#              0123456789012345678901234567890123456
	#              0         1         2         3
	  reads  => [ "CG" ],
	  args   => "",
	  report => "-a -P \"SEED=0,1,1;NCEIL=0,0;ROWM=3\"",
	  hits   => [ { 13 => 2, 23 => 2, 31 => 2 } ],
	  edits  => [ ],
	  flags => [ "XM:0,XP:0,XT:UU,XC:2=",
	             "XM:0,XP:0,XT:UU,XC:2=",
	             "XM:0,XP:0,XT:UU,XC:2=" ],
	  cigar  => [ "2M", "2M", "2M" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:2" => 1 },
		{ "YT:Z:UU" => 1, "MD:Z:2" => 1 },
		{ "YT:Z:UU" => 1, "MD:Z:2" => 1 },
	  ]
	},

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
	  flags => [ "XM:0,XP:0,XT:UU,XC:5=1X2=", "XM:0,XP:0,XT:UU,XC:5=1X2=" ],
	  cigar  => [ "8M", "8M" ],
	  samoptflags => [
		{ "AS:i:-1" => 1, "NM:i:1" => 1, "XM:i:1" => 1, "XN:i:1" => 1,
		  "YT:Z:UU" => 1, "MD:Z:5N2" => 1 },
		{ "AS:i:-1" => 1, "NM:i:1" => 1, "XM:i:1" => 1, "XN:i:1" => 1,
		  "YT:Z:UU" => 1, "MD:Z:5N2" => 1 },
	  ]
	},

	#
	# Alignment with multi-character read gap
	#

	# Relatively small example with a read gap extend
	{ ref    => [ "ATAACCTTCG" ],
	  reads  => [ "ATAATTCG" ], # 3 * 19 + 3 = 60
	  #                ^
	  #                4:CC>- 
	  args   => "",
	  report => "-a --overhang --gbar 3 -P \"MMP=C30;RDG=5,5;SEED=0,4,1;RFG=25,20;MIN=-3,-3\"",
	  hits   => [ { 0 => 1 } ],
	  edits  => [ "4:CC>-" ],
	  flags  => [ "XM:0,XP:0,XT:UU,XC:4=2D4=" ],
	  cigar  => [ "4M2D4M" ],
	  samoptflags => [
		{ "AS:i:-15" => 1, "NM:i:2" => 1,
		  "XO:i:1" => 1, "XG:i:2" => 3, "YT:Z:UU" => 1, "MD:Z:4^CC4" => 1 }
	  ]
	},

	# Reads 1 and 2 don't have overhang, reads 3 and 4 overhang opposite ends
	{ ref    => [ "ATATGCCCCATGCCCCCCTCCG" ],
	  reads  => [ "ATATGCCCCCCCCCCTCCG" ], # 3 * 19 + 3 = 60
	  #                     ^
	  #                     9:ATG>- 
	  args   =>   "-P \"SEED=0,8,1;MMP=C30;RDG=5,5;RFG=25,15;MIN=-3,-3\"",
	  hits   => [ { 0 => 1 } ],
	  edits  => [ "9:ATG>-" ],
	  norc   => 1,
	  flags => [ "XM:0,XP:0,XT:UU,XC:9=3D10=" ],
	  cigar  => [ "9M3D10M" ],
	  samoptflags => [
		{ "AS:i:-20" => 1, "NM:i:3" => 1,
		  "XO:i:1" => 1, "XG:i:3" => 3, "YT:Z:UU" => 1, "MD:Z:9^ATG10" => 1 }
	  ]
	},

	# Reads 1 and 2 don't have overhang, reads 3 and 4 overhang opposite ends
	{ ref    => [ "ATATGCCCCATGCCCCCCTCCG" ],
	  reads  => [ "CGGAGGGGGGGGGGCATAT" ],
	  #            ATATGCCCCCCCCCCTCCG
	  #                     ^         
	  #                     10:GTA>- 
	  args   => "",
	  report => "-a --overhang -P \"SEED=0,8,1;MMP=C30;RDG=5,5;RFG=25,20;MIN=-3,-3\"",
	  hits   => [ { 0 => 1 } ],
	  edits  => [ "10:GTA>-" ],
	  norc   => 1,
	  flags => [ "XM:0,XP:0,XT:UU,XC:9=3D10=" ],
	  cigar  => [ "9M3D10M" ],
	  samoptflags => [
		{ "AS:i:-20" => 1, "NM:i:3" => 1,
		  "XO:i:1" => 1, "XG:i:3" => 3, "YT:Z:UU" => 1, "MD:Z:9^ATG10" => 1 }
	  ]
	},

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
	  mate1fw => 1,  mate2fw => 1,
	  args     =>   "-I 0 -X 35",
	  # Not really any way to flag an alignment as discordant
	  pairhits => [ { "3,30" => 1 }, { "3,29" => 1 }, { "3,28" => 1 }, { "3,27" => 1 },
	                { "4,30" => 1 }, { "4,29" => 1 }, { "4,28" => 1 }, { "4,27" => 1 } ],
	  flags    => [ "XM:0,XP:0,XT:DP,XC:10=", "XM:0,XP:0,XT:DP,XC:10=",
	                "XM:0,XP:0,XT:CP,XC:10=", "XM:0,XP:0,XT:CP,XC:10=",
	                "XM:0,XP:0,XT:DP,XC:10=", "XM:0,XP:0,XT:CP,XC:10=",
					"XM:0,XP:0,XT:CP,XC:10=", "XM:0,XP:0,XT:CP,XC:10=" ] },

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
	  mate1fw => 1,  mate2fw => 1,
	  args   =>   "-I 0 -X 36",
	  # Not really any way to flag an alignment as discordant
	  pairhits => [ { "3,30" => 1 }, { "3,27" => 1 } ],
	  flags => [ "XM:0,XP:0,XT:DP,XC:10=", "XM:0,XP:0,XT:CP,XC:10=" ] },

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
	  mate1fw => 1,  mate2fw => 1,
	  args   =>   "-I 0 -X 36",
	  # Not really any way to flag an alignment as discordant
	  pairhits => [ { "3,30" => 1 } ],
	  flags => [ "XM:0,XP:0,XT:DP,XC:10=" ] },

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
	  mate1fw => 1,  mate2fw => 1,
	  args   =>   "-I 28 -X 80",
	  # Not really any way to flag an alignment as discordant
	  pairhits => [ { "3,20" => 1 } ],
	  flags => [ "XM:0,XP:0,XT:DP,XC:10=" ] },

	# Like 6, but with -M limit

	{ name => "Simple paired-end 9",
	  ref    => [ "CCCATATATATATCCTCCCATATATATATCCCTCCCCATATATATATCCCTTTTCCTTTCGCGCGCGCGTTTCCCCCCCCC" ],
	#                 ATATATATAT      ATATATATAT        ATATATATAT            CGCGCGCGCG
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s  => [ "ATATATATAT" ],
	  mate2s  => [ "CGCGCGCGCG" ],
	  mate1fw => 1,  mate2fw => 0,
	  args    => "-I 0 -X 80 -M 2 -P \"ROWM=3.0\"",
	  lines   => 2,
	  pairhits => [ { "3,59" => 1, "19,59" => 1, "37,59" => 1 } ],
	  hits_are_superset => [ 1 ],
	  flags  => [ "XM:1,XP:1,XT:CP,XC:10=", "XM:1,XP:1,XT:CP,XC:10=" ] },

	# Like 6, but without -m limit

	{ name => "Simple paired-end 8",
	  ref    => [ "CCCATATATATATCCTCCCATATATATATCCCTTCCCATATATATATCCCTTTTTTTTTCGCGCGCGCGTTTCCCCCCCCC" ],
	#                 ATATATATAT      ATATATATAT        ATATATATAT            CGCGCGCGCG
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CGCGCGCGCG" ],
	  mate1fw => 1,  mate2fw => 0,
	  args   =>   "-I 0 -X 80 -P \"ROWM=3.0\"",
	  pairhits => [ { "3,59" => 1, "19,59" => 1, "37,59" => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:CP,XC:10=" ] },

	# Like 6, but with lower -m limit

	{ name => "Simple paired-end 7",
	  ref    => [ "TTTATATATATATTTTTTTATATATATATTTTTTTTTATATATATATTTTTTTTCCCCCCGCGCGCGCGCCCCCCCCCCCC" ],
	#                 ATATATATAT      ATATATATAT        ATATATATAT            CGCGCGCGCG
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CGCGCGCGCG" ],
	  mate1fw => 1,  mate2fw => 0,
	  args   =>   "-I 0 -X 80 -m 1",
	  pairhits => [ { "*,*" => 1 } ],
	  flags  => [ "XM:1,XP:1,XT:CP", "XM:1,XP:1,XT:CP" ] },

	# Same but with mates reversed, first mate aligns 3 times

	{ name => "Simple paired-end 6",
	  ref    => [ 
	   "CCCATATATATATCCCCCCATATATATATCCCTTCCCATATATATATCCCTTTTCCTTTCGCGCGCGCGTTTCCCCCCCCC" ],
	#      ATATATATAT      ATATATATAT        ATATATATAT            CGCGCGCGCG
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CGCGCGCGCG" ],
	  mate1fw => 1,  mate2fw => 0,
	  args   =>   "-I 0 -X 80 -m 1 --norc",
	  pairhits => [ { "*,59" => 1 } ],
	  flags_map => [ { "*"  => "XM:1,XP:1,XT:UP",
	                   "59" => "XM:0,XP:1,XT:UP,XC:10=" } ] },

	# Same but with mates reversed, second mate doesn't align

	{ name => "Simple paired-end 5",
	  ref    => [ "CCCATATATATATCCCTCCATATATATATCCCTTTCCATATATATATCCCTTTTCCCTTCGCGCGCGCGTTTCCCCCCCCC" ],
	#                 ATATATATAT      ATATATATAT        ATATATATAT            "CCCCCGGGGG"
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CCCCCGGGGG" ],
	  mate1fw => 1,  mate2fw => 0,
	  args   =>   "-I 0 -X 80 -m 2",
	  pairhits => [ { "*,*" => 1 } ] },

	# Paired-end read, but only the first mate aligns within the -m 2 ceiling.
	# Second mate aligns 3 places.

	{ name => "Simple paired-end 4",
	  ref    => [
	   "CCCATATATATATCCCTTTTTTTCCCCCCCCCTTTCGCGCGCGCGTTTCTTTCGCGCGCGCGTTTCCCTTTCGCGCGCGCG" ],
	#      ATATATATAT                      CGCGCGCGCG       CGCGCGCGCG         CGCGCGCGCG
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CGCGCGCGCG" ],
	  mate1fw => 1,  mate2fw => 1,
	  args   =>   "-I 0 -X 80 -m 1 --norc",
	  pairhits => [ { "*,3" => 1 } ],
	  flags_map  => [ {  3  => "XM:i:0,XP:1,XT:UP,XC:10=",
	                    "*" => "XM:i:1,XP:1,XT:UP" } ] },

	# Paired-end read, but only the first mate aligns within the -m 2 ceiling.
	# Second mate doesn't align at all.

	{ name => "Simple paired-end 3",
	  ref    => [ 
	   "CCCATATATATATCCCTTTTTTTCCCCCCCCCTTTCGCGCGCGCGTTTCCTTCGCGCGCGCGTTTTCCTTTCGCGCGCGCG" ],
	#      ATATATATAT                      CGCGCGCGCG       CGCGCGCGCG         CGCGCGCGCG
	#              012345678901234567890123456789012345678901234567890123456789012345678901234567890
	#              0         1         2         3         4         5         6         7         8
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CCCCCGGGGG" ],
	  mate1fw => 1,  mate2fw => 1,
	  args   =>   "-I 0 -X 80 -m 1 --norc",
	  norc   => 1,
	  pairhits  => [ { "*,3"  => 1 } ],
	  flags_map => [ {  3     => "XM:0,XP:0,XT:UP,XC:10=",
	                   "*"    => "XM:1,XP:1,XT:UP" } ],
	  cigar_map => [{
		3 => "10M",
		"*" => "*"
	  }],
	  samoptflags_map => [{
		3 => {
			"MD:Z:10"  => 1, # mismatching positions/bases
			"YT:Z:UP"  => 1, # type of alignment (concordant/discordant/etc)
		},
		"*" => {
			"YT:Z:UP"  => 1, # type of alignment (concordant/discordant/etc)
		}
	  }]
	},

	# Paired-end read, but only one mate aligns

	{ name => "Simple paired-end 2; no --no-mixed",
	  ref    => [ "CCCATATATATATCCCTTTTTTTCCCCCCCCCCTTCGCGCGCGCGTTTCCCCC" ],
	#                 ATATATATAT                      CGCGCGCGCG
	#              01234567890123456789012345678901234567890123456789012
	#              0         1         2         3         4         5
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CCCCCGGGGG" ],
	  mate1fw => 1,  mate2fw => 1,
	  args   =>   "-I 0 -X 50 --nofw",
	  nofw   => 1,
	  pairhits  => [ { "*,3"  => 1 } ],
	  flags_map => [ {  3     => "XM:0,XP:0,XT:UP,XC:10=",
	                   "*"    => "XM:0,XP:0,XT:UP" } ],
	  cigar_map => [{
		3 => "10M",
		"*" => "*"
	  }],
	  samoptflags_map => [{
		3 => {
			"MD:Z:10"  => 1, # mismatching positions/bases
			"YT:Z:UP"  => 1, # type of alignment (concordant/discordant/etc)
		},
		"*" => {
			"YT:Z:UP"  => 1, # type of alignment (concordant/discordant/etc)
		}
	  }]
	},

	{ name => "Simple paired-end 2; --no-mixed",
	  ref    => [ "CCCATATATATATCCCTTTTTTTCCCCCCCCTTTTCGCGCGCGCGTTTCCCCC" ],
	#                 ATATATATAT                      CGCGCGCGCG
	#              01234567890123456789012345678901234567890123456789012
	#              0         1         2         3         4         5
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CCCCCGGGGG" ],
	  mate1fw => 1,  mate2fw => 1,
	  args   =>   "-I 0 -X 50 --no-mixed",
	  pairhits => [ { "*,*" => 1 } ] },

	# Simple paired-end alignment
	
	{ name => "Simple paired-end 1",
	  ref    => [ "CCCATATATATATCCCTTTTTTTCCCCCCCCTTTTCGCGCGCGCGTTTTCCCC" ],
	#                 ATATATATAT                      CGCGCGCGCG
	#              01234567890123456789012345678901234567890123456789012
	#              0         1         2         3         4         5
	  mate1s => [ "ATATATATAT" ],
	  mate2s => [ "CGCGCGCGCG" ],
	  mate1fw => 1,  mate2fw => 1,
	  args   =>   "-I 0 -X 50",
	  pairhits => [ { "3,35" => 1 } ],
	  flags => [ "XM:0,XP:0,XT:CP,XC:10=" ],
	  cigar_map => [{
		3  => "10M",
		35 => "10M"
	  }],
	  samoptflags_map => [{
		3 => {
			"MD:Z:10"  => 1, # mismatching positions/bases
			"YT:Z:CP"  => 1, # type of alignment (concordant/discordant/etc)
		},
		35 => {
			"MD:Z:10"  => 1, # mismatching positions/bases
			"YT:Z:CP"  => 1, # type of alignment (concordant/discordant/etc)
		}
	  }]
	},

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
	  hits     => [ { 1 => 1 } ],
	  color    => 1,
	  flags => [ "XM:0,XP:0,XT:UU,XC:3=1X1=1X" ]
	},

	{ ref      => [ "TTGTTC" ],
	  reads    => [ "01102"  ],
	  dec_seq  => [ "TTGTTC" ],
	  dec_qual => [ "IqqqqI" ],
	  args     =>   "--col-keepends",
	  hits     => [ { 0 => 1 } ],
	  color    => 1,
	  flags => [ "XM:0,XP:0,XT:UU,XC:6=" ],
	  cigar  => [ "6M" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:6" => 1 }
	  ]
	},

	{ ref      => [ "TTGTTC" ],
	  reads    => [ "01102"  ],
	  dec_seq  => [ "TTGTTC" ],
	  dec_qual => [ "IqqqqI" ],
	  args     =>   "--col-keepends -P \"SEED=0,3,1\"",
	  hits     => [ { 0 => 1 } ],
	  color    => 1,
	  flags => [ "XM:0,XP:0,XT:UU,XC:6=" ],
	  cigar  => [ "6M" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:6" => 1 }
	  ]
	},

	{ ref      => [ "TTGTTC" ],
	  reads    => [ "01102" ],
	  dec_seq  => [ "TGTT" ],
	  dec_qual => [ "qqqq" ],
	  args     =>   "-P \"SEED=0,3,1\"",
	  hits     => [ { 1 => 1 } ],
	  color    => 1,
	  flags => [ "XM:0,XP:0,XT:UU,XC:4=" ],
	  cigar  => [ "4M" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:4" => 1 }
	  ]
	},

	{ ref      => [ "TTGTTC" ],
	  reads    => [ "01112" ],
	  #                 ^mm
	  args     =>   "--col-keepends -P \"SEED=0,4,1;MMP=C5\"",
	  hits     => [ { "*" => 1 } ],
	  color    => 1,
	  cigar  => [ "*" ],
	  samoptflags => [{ "YT:Z:UU" => 1 }]
	},

	{ ref      => [ "TTGTTC" ],
	  reads    => [ "01112" ],
	  #              21110
	  #                 ^mm
	  dec_seq  => [ "TTGTTC" ],
	  dec_qual => [ "Iqq!!I" ],
	  args     =>   "--col-keepends -P \"SEED=0,3,1;MMP=C7\"",
	  hits     => [ { 0 => 1 } ],
	  color    => 1,
	  flags => [ "XM:0,XP:0,XT:UU,XC:6=" ],
	  cigar  => [ "6M" ],
	  samoptflags => [
		{ "AS:i:-7" => 1, "YT:Z:UU" => 1, "MD:Z:6" => 1 }
	  ]
	},

	# Two color mismatches carry a smaller penalty than 1 SNP
	{ ref      => [ "TTGTTC" ],
	  #              TTGTAC
	  #                  ^snp
	  reads    => [ "01131" ],
	  dec_seq  => [ "TTGTTC" ],
	  dec_qual => [ "Iqq!!!" ],
	  args     =>   "--col-keepends -P \"SEED=0,3,1;MMP=C5;SNP=11\"",
	  hits     => [ { 0 => 1 } ],
	  color    => 1,
	  flags => [ "XM:0,XP:0,XT:UU,XC:6=" ],
	  cigar  => [ "6M" ],
	  samoptflags => [
		{ "AS:i:-10" => 1, "YT:Z:UU" => 1, "MD:Z:6" => 1 }
	  ]
	},

	# 1 SNP carries a smaller penalty than 2 color mismatches
	{ ref      => [ "TTGTTC" ],
	  #              TTGTAC
	  #                  ^snp
	  reads    => [ "01131" ],
	  dec_seq  => [ "TTGTAC" ],
	  dec_qual => [ "IqqqqI" ],
	  args     =>   "--col-keepends -P \"SEED=0,3,1;MMP=C5;SNP=9\"",
	  hits     => [ { 0 => 1 } ],
	  color    => 1,
	  flags => [ "XM:0,XP:0,XT:UU,XC:4=1X1=" ],
	  cigar  => [ "6M" ],
	  samoptflags => [
		{ "NM:i:1" => 1, "XM:i:1" => 1, "AS:i:-9" => 1,
		  "YT:Z:UU" => 1, "MD:Z:4T1" => 1 }
	  ]
	},

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
	  args   => "-m 1",
	  hits   => [ { 0 => 1 }, { 0 => 1 }, { 0 => 1 }, { 0 => 1 }, { 0 => 1 } ],
	  edits  => [ "-",        "-",        "-",        "-",        "-",       ],
	  flags  => [ "XM:0,XP:0,XT:UU,XC:8=", "XM:0,XP:0,XT:UU,XC:8=",
	              "XM:0,XP:0,XT:UU,XC:8=", "XM:0,XP:0,XT:UU,XC:8=",
	              "XM:0,XP:0,XT:UU,XC:8=" ],
	  cigar  => [ "8M", "8M", "8M", "8M", "8M" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:8" => 1, "YM:i:0" => 1 },
		{ "YT:Z:UU" => 1, "MD:Z:8" => 1, "YM:i:0" => 1 },
		{ "YT:Z:UU" => 1, "MD:Z:8" => 1, "YM:i:0" => 1 },
		{ "YT:Z:UU" => 1, "MD:Z:8" => 1, "YM:i:0" => 1 },
		{ "YT:Z:UU" => 1, "MD:Z:8" => 1, "YM:i:0" => 1 },
	  ]
	},

	{ ref    => [ "TTGTTCGT" ],
	#              TTGTTCGT
	#              TAATTCGT
	#              TTGTTCGT
	  reads  => [ "TTGTTCGT", "TAAAACGT", "TTGTTCGT", "TAATTCGT",    "TTGTTCGT" ],
	  args   => "-m 1 -P \"SEED=0,3,1;MIN=-3,-3;MMP=C9;RFG=25,15;RDG=25,15\"",
	  hits   => [ { 0 => 1 }, { '*' => 1 }, { 0 => 1 }, { 0 => 1 },    { 0 => 1 } ],
	  edits  => [ "-",        undef,        "-",        "1:T>A,2:G>A", "-",        ],
	  flags  => [ "XM:0,XP:0,XT:UU,XC:8=", "XM:0,XP:0,XT:UU",
	              "XM:0,XP:0,XT:UU,XC:8=", "XM:0,XP:0,XT:UU,XC:1=2X5=",
				  "XM:0,XP:0,XT:UU,XC:8=" ],
	  cigar  => [ "8M", "*", "8M", "8M", "8M" ],
	  samoptflags => [
		{ "YT:Z:UU"  => 1, "MD:Z:8" => 1, "YM:i:0" => 1 },
		{ "YT:Z:UU"  => 1, "YM:i:0" => 1 },
		{ "YT:Z:UU"  => 1, "MD:Z:8" => 1, "YM:i:0" => 1 },
		{ "AS:i:-18" => 1, "NM:i:2" => 1, "XM:i:2" => 1, "YT:Z:UU" => 1, "MD:Z:1TG5" => 1, "YM:i:0" => 1 },
		{ "YT:Z:UU"  => 1, "MD:Z:8" => 1, "YM:i:0" => 1 },
	  ],
	  norc   => 1 },

	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => "-m 3",
	  mapq   => [ 0 ],
	  hits   => [ { 0 => 1, 8 => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU,XC:8=", "XM:0,XP:0,XT:UU,XC:8=" ],
	  cigar  => [ "8M" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:8" => 1, "YM:i:0" => 1 }
	  ],
	},

	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => "-m 2",
	  mapq   => [ 0 ],
	  hits   => [ { 0 => 1, 8 => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU,XC:8=", "XM:0,XP:0,XT:UU,XC:8=" ],
	  cigar  => [ "8M" ],
	  samoptflags => [
		{ "YT:Z:UU" => 1, "MD:Z:8" => 1, "YM:i:0" => 1 }
	  ],
	},

	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => "-m 1",
	  mapq   => [ 0 ],
	  hits   => [ { "*" => 1 } ],
	  flags  => [ "XM:1,XP:0,XT:UU" ],
	  cigar  => [ "*" ],
	  samoptflags => [
		{ "XM:i:1" => 1, "YT:Z:UU" => 1, "YM:i:1" => 1 }
	  ],
	},

	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => "-M 1",
	  mapq   => [ 0 ],
	  hits   => [ { 0 => 1, 8 => 1 } ],
	  flags  => [ "XM:1,XP:0,XT:UU,XC:8=" ],
	  hits_are_superset => [ 1 ],
	  cigar  => [ "8M" ],
	  samoptflags => [
		{ "YM:i:1" => 1, "YT:Z:UU" => 1, "MD:Z:8" => 1, "YM:i:1" => 1 }
	  ],
	},

	# Read 3 overhangs right end
	{ ref    => [ "TTGTTCGT"  ],
	  reads  => [   "GTTCGTA" ],
	  args   => "--overhang -P \"SEED=0,3,1;NCEIL=2,0\"",
	  hits   => [ { 2 => 1 } ],
	  flags => [ "XM:0,XP:0,XT:UU,XC:6=1X" ] },

	# Mess with arguments
	
	# Default should be 1-mismatch, so this shouldn't align
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTATTAGT" ],
	  args   => "",
	  hits   => [ { "*" => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU" ],
	  cigar  => [ "*" ],
	  samoptflags => [{ "YT:Z:UU" => 1 }],
	},

	# Shouldn't align with 0 mismatches either
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTATTAGT" ],
	  args   => "-P SEED=0",
	  hits   => [ { "*" => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU" ],
	  cigar  => [ "*" ],
	  samoptflags => [{ "YT:Z:UU" => 1 }],
	},

	# Should align with 2 mismatches, provided the mismatches are in the right spots
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TAGTTCAT" ],
	  #            TAGTTCAT
	  #                    TAGTTCAT
	  args   => "-P \"SEED=2;MMP=C1\"",
	  hits   => [ { 0 => 1, 8 => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU,XC:1=1X4=1X1=" ],
	  cigar  => [ "8M" ],
	  samoptflags => [{ "AS:i:-2" => 1, "XS:i:-2" => 1, "NM:i:2" => 1,
	                    "XM:i:2" => 1, "YT:Z:UU" => 1, "MD:Z:1T4G1" => 1 }],
	},

	# Should align with 2 mismatches, provided the mismatches are in the right spots
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTATTAGT" ],
	  args   =>   "-P \"SEED=2;MMP=C1\"",
	  hits   => [ { "*" => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU" ],
	  cigar  => [ "*" ],
	  samoptflags => [{ "YT:Z:UU" => 1 }],
	},

	# Should align with 0 mismatches if we can wedge a seed into the 2
	# matching characters between the two mismatches.  Here we fail to
	# wedge a length-3 seed in (there's no room)
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTATTAGT" ],
	  args   => "-P \"SEED=0,3,1;MMP=C1\"",
	  hits   => [ { "*" => 1 } ],
	  flags  => [ "XM:0,XP:0,XT:UU" ],
	  cigar  => [ "*" ],
	  samoptflags => [{ "YT:Z:UU" => 1 }],
	},

	# Should align with 0 mismatches if we can wedge a seed into the 2
	# matching characters between the two mismatches.  Here we wedge a
	# length-2 seed in
	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [      "TTATTAGT" ],
	  args   => "-P \"SEED=0,2,1;MMP=C1;ROWM=10\"",
	  #
	  # TTGTTCGTTTGTTCGT TTGTTCGTTTGTTCGT TTGTTCGTTTGTTCGT
	  # || || ||            ||  |             |  || ||
	  # TTATTAGT            TTATTAGT          TTATTAGT
	  #
	  # TTGTTCGTTTGTTCGT TTGTTCGTTTGTTCGT TTGTTCGTTTGTTCGT
	  #         ||  |           ||  |             || || ||
	  #      TTATTAGT           TTATTAGT          TTATTAGT
	  #
	  hits   =>   [ { 0 => 1, 3 => 1, 4 => 1,
	                  5 => 1, 7 => 1, 8 => 1} ],
	  flag_map => [ { 0 => "XM:0,XP:0,XT:UU,XC:2=1X2=1X2=",
	                  3 => "XM:0,XP:0,XT:UU,XC:2=2X1=3X",
					  4 => "XM:0,XP:0,XT:UU,XC:1=2X2=1X2=",
					  5 => "XM:0,XP:0,XT:UU,XC:3X2=2X1=",
					  7 => "XM:0,XP:0,XT:UU,XC:2=2X1=3X",
					  8 => "XM:0,XP:0,XT:UU,XC:2=1X2=1X2="} ],
	  cigar_map => [ { 0 => "8M", 3 => "8M", 4 => "8M",
	                   5 => "8M", 7 => "8M", 8 => "8M" } ],
	  samoptflags_map => [{
		0 => { "AS:i:-2" => 1, "XS:i:-2" => 1, "NM:i:2" => 1, "XM:i:2" => 1,
		       "YT:Z:UU" => 1, "MD:Z:2G2C2"    => 1 },
		3 => { "AS:i:-5" => 1, "XS:i:-2" => 1, "NM:i:5" => 1, "XM:i:5" => 1,
		       "YT:Z:UU" => 1, "MD:Z:2CG1TTG0" => 1 },
		4 => { "AS:i:-3" => 1, "XS:i:-2" => 1, "NM:i:3" => 1, "XM:i:3" => 1,
		       "YT:Z:UU" => 1, "MD:Z:1CG2T2"   => 1 },
		5 => { "AS:i:-5" => 1, "XS:i:-2" => 1, "NM:i:5" => 1, "XM:i:5" => 1,
		       "YT:Z:UU" => 1, "MD:Z:0CGT2GT1" => 1 },
		7 => { "AS:i:-5" => 1, "XS:i:-2" => 1, "NM:i:5" => 1, "XM:i:5" => 1,
		       "YT:Z:UU" => 1, "MD:Z:2TG1TCG0" => 1 },
		8 => { "AS:i:-2" => 1, "XS:i:-2" => 1, "NM:i:2" => 1, "XM:i:2" => 1,
		       "YT:Z:UU" => 1, "MD:Z:2G2C2"    => 1 },
	  }],
	},
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
sub runbowtie2($$$$$$$$$$$$$$$$$$) {
	my (
		$do_build,
		$args,
		$color,
		$fa,
		$reportargs,       #5
		$read_file_format,
		$read_file,
		$mate1_file,
		$mate2_file,
		$reads,
		$mate1s,
		$mate2s,
		$names,
		$ls,
		$rawls,
		$header_ls,
		$raw_header_ls,
		$should_abort) = @_;
	$args .= " --quiet";
	$args .= " --print-flags";
	$reportargs = "-a" unless defined($reportargs);
	$args .= " -C" if $color;
	$args .= " $reportargs";
	# Write the reference to a fasta file
	print "References:\n";
	open(FA, $fa) || die;
	while(<FA>) { print $_; }
	close(FA);
	if($do_build) {
		my $build_args = ($color ? "-C" : "");
		my $cmd = "$bowtie2_build --quiet --sanity $build_args $fa .simple_tests.tmp";
		print "$cmd\n";
		system($cmd);
		($? == 0) || die "Bad exitlevel from bowtie2-build: $?";
	}
	my $pe = (defined($mate1s) && $mate1s ne "");
	$pe = $pe || (defined($mate1_file));
	my $mate1arg;
	my $mate2arg;
	my $readarg;
	my $formatarg = "-c";
	my ($readstr, $m1str, $m2str) = (undef, undef, undef);
	$readstr = join(",", @$reads)  if defined($reads);
	$m1str   = join(",", @$mate1s) if defined($mate1s);
	$m2str   = join(",", @$mate2s) if defined($mate2s);
	if(defined($read_file) || defined($mate1_file)) {
		defined($read_file_format) || die;
		my $ext = "";
		if($read_file_format eq "fastq") {
			$formatarg = "-q";
			$ext = ".fq";
		} elsif($read_file_format eq "tabbed") {
			$formatarg = "--12";
			$ext = ".tab";
		} elsif($read_file_format eq "fasta") {
			$formatarg = "-f";
			$ext = ".fa";
		} elsif($read_file_format eq "qseq") {
			$formatarg = "--qseq";
			$ext = "_qseq.txt";
		} elsif($read_file_format eq "raw") {
			$formatarg = "-r";
			$ext = ".raw";
		} else {
			die "Bad format: $read_file_format";
		}
		if(defined($read_file)) {
			# Unpaired
			open(RD, ">.simple_tests$ext") || die;
			print RD $read_file;
			close(RD);
			$readarg = ".simple_tests$ext";
		} else {
			defined($mate1_file) || die;
			defined($mate2_file) || die;
			# Paired
			open(M1, ">.simple_tests.1$ext") || die;
			print M1 $mate1_file;
			close(M1);
			open(M2, ">.simple_tests.2$ext") || die;
			print M2 $mate2_file;
			close(M2);
			$mate1arg = ".simple_tests.1$ext";
			$mate2arg = ".simple_tests.2$ext";
		}
	} else {
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
	}
	my $cmd;
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
		if(substr($_, 0, 1) eq "@") {
			push @$header_ls, [ split(/\t/, $_, -1) ];
			push @$raw_header_ls, $_;
		} else {
			push @$ls, [ split(/\t/, $_, -1) ];
			push @$rawls, $_;
		}
	}
	close(BT);
	($? == 0 ||  $should_abort) || die "bowtie2 aborted with exitlevel $?\n";
	($? != 0 || !$should_abort) || die "bowtie2 failed to abort!\n";
}

##
# Compare a hash ref of expected SAM flags with a hash ref of observed SAM
# flags.
#
sub matchSamOptionalFlags($$) {
	my ($flags, $ex_flags) = @_;
	my %ex = ();
	for(keys %$ex_flags) {
		my ($nm, $ty, $vl) = split(/:/, $_);
		defined($vl) || die "Could not parse optional flag field \"$_\"";
		($ex{$nm}{ty}, $ex{$nm}{vl}) = ($ty, $vl);
	}
	for(keys %$flags) {
		my ($ex_ty, $ex_vl);
		if(defined($ex{$_})) {
			($ex_ty, $ex_vl) = ($ex{$_}{ty}, $ex{$_}{vl});
		} else {
			($ex_ty, $ex_vl) = ("i", "0");
		}
		defined($ex_ty) || die;
		defined($ex_vl) || die;
		my ($ty, $vl) = ($flags->{$_}{ty}, $flags->{$_}{vl});
		defined($ty) || die;
		defined($vl) || die;
		$ex_ty eq $ty ||
			die "Expected SAM optional flag $_ to have type $ex_ty, had $ty";
		$ex_vl eq $vl ||
			die "Expected SAM optional flag $_ to have value $ex_vl, had $vl";
	}
	return 1;
}

my $tmpfafn = ".simple_tests.pl.fa";
my $last_ref = undef;
for (my $ci = 0; $ci < scalar(@cases); $ci++) {
	my $c = $cases[$ci];
	# If there's any skipping of cases to be done, do it here prior to the
	# eq_deeply check
	my $color = 0;
	$color = $c->{color} if defined($c->{color});
	next if ($color && $skipColor);
	my $do_build = 0;
	unless(defined($last_ref) && eq_deeply($c->{ref}, $last_ref)) {
		writeFasta($c->{ref}, $tmpfafn);
		$do_build = 1;
	}
	$last_ref = $c->{ref};
	# For each set of arguments...
	my $case_args = $c->{args};
	$case_args = "" unless defined($case_args);
	my $first = 1; # did we build the index yet?
	# Forward, then reverse-complemented
	my $fwlo = ($c->{nofw} ? 1 : 0);
	my $fwhi = ($c->{norc} ? 0 : 1);
	for(my $fwi = $fwlo; $fwi <= $fwhi; $fwi++) {
		my $fw = ($fwi == 0);
		my $sam = 1;
		
		my $reads     = $c->{reads};
		my $m1s       = $c->{mate1s};
		my $m2s       = $c->{mate2s};
		
		my $read_file = undef;
		my $mate1_file = undef;
		my $mate2_file = undef;
		
		$read_file  = $c->{fastq}   if defined($c->{fastq});
		$read_file  = $c->{tabbed}  if defined($c->{tabbed});
		$read_file  = $c->{fasta}   if defined($c->{fasta});
		$read_file  = $c->{qseq}    if defined($c->{qseq});
		$read_file  = $c->{raw}     if defined($c->{raw});

		$mate1_file = $c->{fastq1}  if defined($c->{fastq1});
		$mate1_file = $c->{tabbed1} if defined($c->{tabbed1});
		$mate1_file = $c->{fasta1}  if defined($c->{fasta1});
		$mate1_file = $c->{qseq1}   if defined($c->{qseq1});
		$mate1_file = $c->{raw1}    if defined($c->{raw1});

		$mate2_file = $c->{fastq2}  if defined($c->{fastq2});
		$mate2_file = $c->{tabbed2} if defined($c->{tabbed2});
		$mate2_file = $c->{fasta2}  if defined($c->{fasta2});
		$mate2_file = $c->{qseq2}   if defined($c->{qseq2});
		$mate2_file = $c->{raw2}    if defined($c->{raw2});
		
		my $read_file_format = undef;
		if(!defined($reads) && !defined($m1s) && !defined($m2s)) {
			defined($read_file) || defined($mate1_file) || die;
			$read_file_format = "fastq"  if defined($c->{fastq})  || defined($c->{fastq1});
			$read_file_format = "tabbed" if defined($c->{tabbed}) || defined($c->{tabbed});
			$read_file_format = "fasta"  if defined($c->{fasta})  || defined($c->{fasta1});
			$read_file_format = "qseq"   if defined($c->{qseq})   || defined($c->{qseq1});
			$read_file_format = "raw"    if defined($c->{raw})    || defined($c->{raw1});
			next unless $fw;
		}
		# Run bowtie2
		my @lines = ();
		my @rawlines = ();
		my @header_lines = ();
		my @header_rawlines = ();
		print $c->{name}." " if defined($c->{name});
		print "(fw:".($fw ? 1 : 0).", sam:$sam)\n";
		my $mate1fw = 1;
		my $mate2fw = 0;
		$mate1fw = $c->{mate1fw} if defined($c->{mate1fw});
		$mate2fw = $c->{mate2fw} if defined($c->{mate2fw});
		if(!$fw) {
			# Reverse-complement the reads
			my @s = (); @s = @$reads if defined($reads);
			# Reverse-complement mates and switch mate1 with mate2
			my @m1 = (); @m1 = @$m1s if defined($m1s);
			my @m2 = (); @m2 = @$m2s if defined($m2s);
			for(0..scalar(@s )-1) { $s [$_] = DNA::revcomp($s [$_], $color); }
			if($mate1fw == $mate2fw) {
				for(0..$#m1) { $m1[$_] = DNA::revcomp($m1[$_], $color); }
				for(0..$#m2) { $m2[$_] = DNA::revcomp($m2[$_], $color); }
			}
			$reads = \@s if defined($reads);
			$m1s = \@m2 if defined($m1s);
			$m2s = \@m1 if defined($m2s);
		}
		my $a = $case_args;
		if(defined($m2s)) {
			$a .= " --";
			$a .= ($mate1fw ? "f" : "r");
			$a .= ($mate2fw ? "f" : "r");
		}
		runbowtie2(
			$do_build && $first,
			"$a",
			$color,
			$tmpfafn,
			$c->{report},
			$read_file_format, # formate of read/mate files
			$read_file,        # read file
			$mate1_file,       # mate #1 file
			$mate2_file,       # mate #2 file
			$reads,            # read list
			$m1s,              # mate #1 list
			$m2s,              # mate #2 list
			$c->{names},
			\@lines,
			\@rawlines,
			\@header_lines,
			\@header_rawlines,
			$c->{should_abort});
		$first = 0;
		my $pe = defined($c->{mate1s}) && $c->{mate1s} ne "";
		$pe = $pe || defined($mate1_file);
		$pe = $pe || $c->{paired};
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
			my ($readname, $orient, $chr, $off, $seq, $qual, $mapq, $oms,
				$editstr, $flagstr, $samflags, $cigar, $rnext, $pnext, $tlen);
			my %samoptflags = ();
			if($sam) {
				scalar(@$l) >= 11 ||
					die "Bad number of fields; expected at least 11 got ".
						scalar(@$l).":\n$rawlines[$li]\n";
				($readname, $samflags, $chr, $off) = @$l[0..3];
				($seq, $qual) = @$l[9..10];
				$orient = ((($samflags >> 4) & 1) == 0) ? "+" : "-";
				$mapq  = $l->[4]; # mapping quality
				$cigar = $l->[5]; # CIGAR string
				$rnext = $l->[6]; # ref seq of next frag in template
				$pnext = $l->[7]; # position of next frag in template
				$tlen  = $l->[8]; # template length
				if($pnext == 0) { $pnext = "*"; } else { $pnext--; }
				for(my $m = 11; $m < scalar(@$l); $m++) {
					next if $l->[$m] eq "";
					my ($nm, $ty, $vl) = split(/:/, $l->[$m]);
					defined($vl) ||
						die "Could not parse optional flag field $m: ".
							"\"$l->[$m]\"";
					$samoptflags{$nm}{ty} = $ty;
					$samoptflags{$nm}{vl} = $vl;
				}
				if($off > 0) { $off--; }
				else { $off = "*"; }
			} else {
				scalar(@$l) == 9 ||
					die "Bad number of fields; expected 9 got ".
						scalar(@$l).":\n$rawlines[$li]\n";
				($readname, $orient, $chr, $off, $seq,
				 $qual, $oms, $editstr, $flagstr) = @$l;
			}
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
			$rdi = $c->{idx_map}{$rdi} if defined($c->{idx_map}{$rdi});
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
			%hits = %{$c->{hits}->[$rdi]} if
				defined($c->{hits}->[$rdi]);
			# 'flags'
			my $flags = undef;
			$flags = $c->{flags}->[$rdi] if
				defined($c->{flags}->[$rdi]);
			# 'samflags'
			my $ex_samflags = undef;
			$ex_samflags = $c->{ex_samflags}->[$rdi] if
				defined($c->{ex_samflags}->[$rdi]);
			# 'samflags_map'
			my $ex_samflags_map = undef;
			$ex_samflags_map = $c->{samflags_map}->[$rdi] if
				defined($c->{samflags_map}->[$rdi]);
			# 'samoptflags'
			my $ex_samoptflags = undef;
			$ex_samoptflags = $c->{samoptflags}->[$rdi] if
				defined($c->{samoptflags}->[$rdi]);
			# 'cigar'
			my $ex_cigar = undef;
			$ex_cigar = $c->{cigar}->[$rdi] if
				defined($c->{cigar}->[$rdi]);
			# 'cigar_map'
			my $ex_cigar_map = undef;
			$ex_cigar_map = $c->{cigar_map}->[$rdi] if
				defined($c->{cigar_map}->[$rdi]);
			# 'mapq'
			my $ex_mapq = undef;
			$ex_mapq = $c->{mapq}->[$rdi] if
				defined($c->{mapq}->[$rdi]);
			# 'mapq_map'
			my $ex_mapq_map = undef;
			$ex_mapq_map = $c->{mapq_map}->[$rdi] if
				defined($c->{mapq_map}->[$rdi]);
			# 'rnext_map'
			my $ex_rnext_map = undef;
			$ex_rnext_map = $c->{rnext_map}->[$rdi] if
				defined($c->{rnext_map}) && defined($c->{rnext_map}->[$rdi]);
			# 'pnext_map'
			my $ex_pnext_map = undef;
			$ex_pnext_map = $c->{pnext_map}->[$rdi] if
				defined($c->{pnext_map}) && defined($c->{pnext_map}->[$rdi]);
			# 'tlen_map'
			my $ex_tlen_map = undef;
			$ex_tlen_map = $c->{tlen_map}->[$rdi] if
				defined($c->{tlen_map}) && defined($c->{tlen_map}->[$rdi]);
			# 'flags_fw'
			my $flags_fw = undef;
			$flags_fw = $c->{flags_fw}->[$rdi] if
				defined($c->{flags_fw}->[$rdi]);
			# 'flags_rc'
			my $flags_rc = undef;
			$flags_rc = $c->{flags_rc}->[$rdi] if
				defined($c->{flags_rc}->[$rdi]);
			# 'pairhits'
			my %pairhits = ();
			%pairhits = %{$c->{pairhits}->[$rdi]} if
				defined($c->{pairhits}->[$rdi]);
			# 'pairflags'
			my %pairflags = ();
			%pairflags = %{$c->{pairflags}->[$rdi]} if
				defined($c->{pairflags}->[$rdi]);
			# 'hits_are_superset'
			my $hits_are_superset = 0;
			$hits_are_superset = $c->{hits_are_superset}->[$rdi] if
				defined($ci);
			# edits
			my $ex_edits = undef;
			$ex_edits = $c->{edits}->[$rdi] if
				defined($c->{edits}->[$rdi]);
			if(!$sam) {
				# Bowtie flags
				if(defined($flags)) {
					$flagstr eq $flags ||
						die "Expected flags=\"$flags\", got \"$flagstr\"";
				}
				if(defined($flags_fw) && $fw) {
					$flagstr eq $flags_fw ||
						die "Expected flags=\"$flags_fw\", got \"$flagstr\"";
				}
				if(defined($flags_rc) && !$fw) {
					$flagstr eq $flags_rc ||
						die "Expected flags=\"$flags_rc\", got \"$flagstr\"";
				}
				if(defined($c->{flag_map})) {
					if(defined($c->{flag_map}->[$rdi]->{$off})) {
						$flagstr eq $c->{flag_map}->[$rdi]->{$off} ||
							die "Expected flags=\"$c->{flag_map}->[$rdi]->{$off}\"".
								" at offset $off, got \"$flagstr\"";
					}
				}
			}
			if($sam) {
				# SAM flags
				if(defined($ex_samflags)) {
					$samflags eq $ex_samflags ||
						die "Expected flags $ex_samflags, got $samflags";
				}
				if(defined($ex_samflags_map)) {
					if(defined($c->{samflags_map}->[$rdi]->{$off})) {
						my $ex = $c->{samflags_map}->[$rdi]->{$off};
						$samflags eq $ex || die
							"Expected FLAGS value $ex at offset $off, got $samflags"
					} else {
						die "Expected to see alignment with offset $off parsing samflags_map";
					}
				}
				# CIGAR string
				if(defined($ex_cigar)) {
					$cigar eq $ex_cigar ||
						die "Expected CIGAR string $ex_cigar, got $cigar";
				}
				if(defined($ex_cigar_map)) {
					if(defined($c->{cigar_map}->[$rdi]->{$off})) {
						my $ex = $c->{cigar_map}->[$rdi]->{$off};
						$cigar eq $ex || die
							"Expected CIGAR string $ex at offset $off, got $cigar"
					} else {
						die "Expected to see alignment with offset $off parsing cigar_map";
					}
				}
				# MAPQ
				if(defined($ex_mapq)) {
					$mapq eq $ex_mapq ||
						die "Expected MAPQ $ex_mapq, got $mapq";
				}
				if(defined($ex_mapq_map)) {
					if(defined($c->{mapq_map}->[$rdi]->{$off})) {
						my $ex = $c->{mapq_map}->[$rdi]->{$off};
						$mapq eq $ex || die
							"Expected MAPQ string $ex at offset $off, got $mapq"
					} else {
						die "Expected to see alignment with offset $off parsing mapq_map";
					}
				}
				# SAM optional flags
				if(defined($ex_samoptflags)) {
					matchSamOptionalFlags(\%samoptflags, $ex_samoptflags);
				}
				if(defined($c->{samoptflags_map})) {
					if(defined($c->{samoptflags_map}->[$rdi]->{$off})) {
						matchSamOptionalFlags(
							\%samoptflags,
							$c->{samoptflags_map}->[$rdi]->{$off});
					} else {
						die "Expected to see alignment with offset $off parsing samoptflags_map";
					}
				}
				if(defined($c->{samoptflags_flagmap})) {
					if(defined($c->{samoptflags_flagmap}->[$rdi]->{$samflags})) {
						matchSamOptionalFlags(
							\%samoptflags,
							$c->{samoptflags_flagmap}->[$rdi]->{$samflags});
					} else {
						die "Expected to see alignment with flag $samflags parsing samoptflags_flagmap";
					}
				}
				# RNEXT map
				if(defined($c->{rnext_map})) {
					if(defined($c->{rnext_map}->[$rdi]->{$off})) {
						my $ex = $c->{rnext_map}->[$rdi]->{$off};
						$rnext eq $ex || die
							"Expected RNEXT '$ex' at offset $off, got '$rnext'"
					} else {
						die "Expected to see alignment with offset $off parsing rnext_map".Dumper($c);
					}
				}
				# PNEXT map
				if(defined($c->{pnext_map})) {
					if(defined($c->{pnext_map}->[$rdi]->{$off})) {
						my $ex = $c->{pnext_map}->[$rdi]->{$off};
						$pnext eq $ex || die
							"Expected PNEXT '$ex' at offset $off, got '$pnext'"
					} else {
						die "Expected to see alignment with offset $off parsing pnext_map";
					}
				}
				# TLEN map
				if(defined($c->{tlen_map})) {
					if(defined($c->{tlen_map}->[$rdi]->{$off})) {
						my $ex = $c->{tlen_map}->[$rdi]->{$off};
						$tlen eq $ex || die
							"Expected TLEN '$ex' at offset $off, got '$tlen'"
					} else {
						die "Expected to see alignment with offset $off parsing tlen_map";
					}
				}
			}
			if($pe && $lastchr ne "") {
				my $offkey = $lastoff.",".$off;
				$offkey = $off.",".$lastoff if $off eq "*";
				if($lastoff ne "*" && $off ne "*") {
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
			if(!$sam && defined($ex_edits)) {
				my $eds = $l->[7];
				$eds eq $ex_edits ||
					die "For edit string, expected \"$ex_edits\" got \"$eds\"\n";
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
				die "Had $hitsLeft hit(s) left over at position $k";
			}
			my $pairhitsLeft = scalar(keys %pairhits);
			if($pairhitsLeft != 0 && !$hits_are_superset) {
				print Dumper(\%pairhits);
				die "Had $pairhitsLeft hit(s) left over at position $k";
			}
		}
		
		$c->{hits} = $hitstmp;
		$c->{pairhits} = $pairhitstmp;
	}
	$last_ref = undef if $first;
}
print "PASSED\n";
