#!/bin/sh

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

#  simple_tests.sh

export ASAN_OPTIONS=halt_on_error=1
export UBSAN_OPTIONS=halt_on_error=1

MAKE=`which gmake || which make`
$MAKE $* bowtie2-align-s \
		bowtie2-align-l \
		bowtie2-align-s-debug \
		bowtie2-align-l-debug \
		bowtie2-align-s-sanitized \
		bowtie2-align-l-sanitized \
		bowtie2-build-s \
		bowtie2-build-l \
		bowtie2-build-s-debug \
		bowtie2-build-l-debug \
		bowtie2-build-s-sanitized \
		bowtie2-build-l-sanitized && \
perl scripts/test/simple_tests.pl \
	--bowtie2=./bowtie2 \
	--bowtie2-build=./bowtie2-build \
	--compressed
