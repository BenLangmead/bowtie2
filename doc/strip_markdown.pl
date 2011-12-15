#!/usr/bin/perl -w

##
# strip_markdown.pl
#
# Used to convert MANUAL.markdown to MANUAL.  Leaves all manual content, but
# strips away some of the clutter that makes it hard to read the markdown.
#

use strict;
use warnings;

my $lastBlank = 0;

while(<>) {
	# Skip comments
	next if /^\s*<!--/;
	next if /^\s*!/;
	next if /^\s*-->/;
	# Skip internal links
	next if /\[.*\]: #/;
	# Skip HTML
	next if /^\s?\s?\s?<.*>\s*$/;
	# Skip HTML
	next if /^\s*<table/;
	next if /^\s*<\/td/;
	next if /^\s*<.*>\s*$/;
	# Strip [`...`]
	s/\[`/`/g;
	s/`\]/`/g;
	# Strip [#...]
	#s/\[#[^\]]*\]//g;
	# Strip (#...)
	s/\(#[^\)]*\)//g;
	# Turn hashes into spaces
	#s/^####/   /;
	#s/^###/ /;
	if(/^\s*$/) {
		next if $lastBlank;
		$lastBlank = 1;
	} else {
		$lastBlank = 0;
	}
	print $_;
}
