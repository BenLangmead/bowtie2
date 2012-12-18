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

CPUS=$1
shift
make -j $CPUS bowtie2-align bowtie2-align-debug bowtie2-build bowtie2-build-debug && \
perl scripts/sim/run.pl \
	--bowtie2=./bowtie2-align \
	--bowtie2-debug=./bowtie2-align-debug \
	--bowtie2-build=./bowtie2-build \
	--bowtie2-build-debug=./bowtie2-build-debug \
	--cpus=$CPUS \
	$*
