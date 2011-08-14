#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 2/23/2011

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin); 
use lib "$Bin";
use lib "$Bin/contrib";
use Sim;
use ForkManager;

# Simulator configuration
my %conf = (
	bowtie2_build       => "bowtie2-build",
	bowtie2             => "bowtie2",
	bowtie2_build_debug => "bowtie2-build-debug",
	bowtie2_debug       => "bowtie2-debug",
	tempdir             => "/tmp"
);

# Number of parallel processes to use
my $cpus = 1;

my $usage = qq!
run.pl [options*]

Options:

  --bowtie2 <path>              Path to bowtie2 release binary
  --bowtie2-debug <path>        Path to bowtie2 debug binary
  --bowtie2-build <path>        Path to bowtie2-build release binary
  --bowtie2-build-debug <path>  Path to bowtie2-build debug binary
  --tempdir <path>              Put temporary files here
  --cases <int>                 Each thread runs around <int> cases (def: 5)
  --cpus <int> / -p <int>       Run test cases in <int> threads at once
  --maxreads <int>              Handle at most <int> reads per case
  --numrefs <int>               Generate <int> refs per case
  --die-with-child              Kill parent as soon as 1 child dies
  --no-die-with-child           Don\'t kill parent as soon as 1 child dies
  --small                       Make small test cases
  --help                        Print this usage message

!;

my $help = 0;
my $ncases = 5;
my $dieWithChild = 1;

GetOptions(
	"bowtie2=s"             => \$conf{bowtie2},
	"bowtie2-debug=s"       => \$conf{bowtie2_debug},
	"bowtie2-build=s"       => \$conf{bowtie2_build},
	"bowtie2-build-debug=s" => \$conf{bowtie2_build_debug},
	"tempdir|tmpdir=s"      => \$conf{tempdir},
	"cases-per-thread=i"    => \$ncases,
	"small"                 => \$conf{small},
	"no-paired"             => \$conf{no_paired},
	"no-color"              => \$conf{no_color},
	"help"                  => \$help,
	"die-with-child"        => \$dieWithChild,
	"no-die-with-child"     => sub { $dieWithChild = 0 },
	"p|cpus=i"              => \$cpus,
	"u|qupto|maxreads=i"    => \$conf{maxreads},
	"numrefs|num-refs=i"    => \$conf{numrefs},
) || die "Bad options;";

if($help) {
	print $usage;
	exit 0;
}

my $sim = Sim->new();
my $pm = new Parallel::ForkManager($cpus); 

# Callback for when a child finishes so we can get its exit code
my @childFailed = ();
my @childFailedPid = ();

$pm->run_on_finish(sub {
	my ($pid, $exit_code, $ident) = @_;
	if($exit_code != 0) {
		push @childFailed, $exit_code;
		push @childFailedPid, $pid;
		!$dieWithChild || die "Dying with child with PID $pid";
	}
});

my $totcases = $ncases * $cpus;
for(1..$totcases) {
	my $childPid = $pm->start;
	if($childPid != 0) {
		next; # spawn the next child
	}
	$sim->nextCase(\%conf);
	$pm->finish;
}
$pm->wait_all_children;
for(@childFailedPid) {
	print STDERR "Error message from child with pid $_:\n";
	my $fn = ".run.pl.child.$_";
	if(open(ER, $fn)) {
		print STDERR "---------\n";
		while(<ER>) {
			print STDERR $_;
		}
		print STDERR "---------\n";
		close(ER);
	} else {
		print STDERR "(could not open $fn)\n";
	}
}
print STDERR "PASSED\n";
