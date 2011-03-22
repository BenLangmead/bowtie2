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
  --cpus <int> / -p <int>       Run test cases in <int> threads at once
  --small                       Make small test cases
  --help                        Print this usage message

!;

my $help = 0;

GetOptions(
	"bowtie2=s"             => \$conf{bowtie2},
	"bowtie2-debug=s"       => \$conf{bowtie2_debug},
	"bowtie2-build=s"       => \$conf{bowtie2_build},
	"bowtie2-build-debug=s" => \$conf{bowtie2_build_debug},
	"tempdir|tmpdir=s"      => \$conf{tempdir},
	"small"                 => \$conf{small},
	"no-paired"             => \$conf{no_paired},
	"help"                  => \$help,
	"p|cpus=i"              => \$cpus
) || die "Bad options;";

# TODO: fix
$conf{no_paired} = 1;

if($help) {
	print $usage;
	exit 0;
}

my $sim = Sim->new();
my $pm = new Parallel::ForkManager($cpus); 

# Callback for when a child finishes so we can get its exit code
my $childFailed = 0;
my $childFailedPid = 0;
$pm->run_on_finish(sub {
	my ($pid, $exit_code, $ident) = @_;
	if($exit_code != 0) {
		$childFailed = $exit_code;
		$childFailedPid = $pid;
	}
});

my $case = $sim->nextCase(\%conf);
for(1..$cpus) {
	my $childPid = $pm->start;
	if($childPid != 0) {
		next; # spawn the next child
	}
	$sim->nextCase(\%conf);
	$pm->finish;
}
$pm->wait_all_children;
if($childFailed) {
	print STDERR "One or more children failed!\n";
	exit 1
}
print STDERR "PASSED\n";
