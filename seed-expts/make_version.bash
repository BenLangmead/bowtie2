#!/bin/bash

# Warning: destroys any existing bowtie2-align-s in parent directory

set -ex

cd ..
patch="seed-expts/versions/${1}/patch.diff"
test -f "${patch}"
test -s "${patch}" && git apply "${patch}"
make bowtie2-align-s
cp bowtie2-align-s "seed-expts/versions/${1}/"
test -s "${patch}" && git apply -R "${patch}"
cd seed-expts
