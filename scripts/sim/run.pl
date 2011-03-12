#!/usr/bin/perl -w

##
# Author: Ben Langmead
#   Date: 2/23/2011

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin); 
use lib $Bin;
use lib $Bin/contrib;
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

GetOptions(
	"bowtie2=s"             => \$conf{bowtie2},
	"bowtie2-debug=s"       => \$conf{bowtie2_debug},
	"bowtie2-build=s"       => \$conf{bowtie2_build},
	"bowtie2-build-debug=s" => \$conf{bowtie2_build_debug},
	"tempdir=s"             => \$conf{tempdir},
	"p|cpus=i"              => \$cpus
) || die "Bad options";

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

my $case = $sim->nextCase( { no_paired => 1 } );
for(1..$cpus) {
	my $childPid = $pm->start;
	if($childPid != 0) {
		next; # spawn the next child
	}
	$sim->nextCase( { no_paired => 1 } );
	$pm->finish;
}
$pm->wait_all_children;
if($childFailed) {
	print STDERR "One or more children failed!\n";
	exit 1
}
print STDERR "PASSED\n";
