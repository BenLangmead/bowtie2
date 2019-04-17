
<!-- badges: start -->
[![Generic badge](https://img.shields.io/badge/version-2.3.5.1-green.svg)](https://shields.io/)
[![Build Status](https://travis-ci.org/BenLangmead/bowtie2.svg?branch=master)](https://travis-ci.org/BenLangmead/bowtie2)
[![License: GPL v3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!--badges: end -->

## Overview

*Bowtie 2* is an ultrafast and memory-efficient tool for aligning sequencing reads
to long reference sequences. It is particularly good at aligning reads of about 50
up to 100s or 1,000s of characters, and particularly good at aligning to relatively
long (e.g. mammalian) genomes. Bowtie 2 indexes the genome with an FM Index to keep
its memory footprint small: for the human genome, its memory footprint is typically
around 3.2 GB. Bowtie 2 supports gapped, local, and paired-end alignment modes.

## Obtaining Bowtie2

*Bowtie 2* is available from various package managers, notably
[Bioconda](https://anaconda.org/bioconda/bowtie2). With Bioconda installed, you
should be able to install Bowtie 2 with `conda install bowtie2`.

Containerized versions of Bowtie 2 are also available via the
[Biocontainers](https://biocontainers.pro/) project (e.g. via
[Docker Hub](https://hub.docker.com/r/biocontainers/bowtie2/)).

You can also download Bowtie 2 sources and binaries from the "releases" tab on this
page. Binaries are available for the x86_64 architecture running Linux, Mac OS X,
and Windows. We are planning on adding experimental support for ARM-64 in an
upcoming release. If you plan to compile Bowtie 2 yourself, make sure you have the
[TBB](https://www.threadingbuildingblocks.org/) and [zlib](https://www.zlib.net)
libraries installed. See the
[Building from source](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#building-from-source)
section of the manual for details.

## Getting started
Looking to try out Bowtie 2? Check out the [Bowtie 2 UI](http://bit.ly/bt2ui-beta) (currently in beta).

### Alignment
`bowtie2` takes a Bowtie 2 index and a set of sequencing read files and outputs a
set of alignments in SAM format.

"Alignment" is the process by which we discover how and where the read sequences are
similar to the reference sequence. An “alignment” is a result from this process,
specifically: an alignment is a way of “lining up” some or all of the characters in
the read with some characters from the reference in a way that reveals how they’re
similar. For example:

```
  Read:      GACTGGGCGATCTCGACTTCG
             |||||  |||||||||| |||
  Reference: GACTG--CGATCTCGACATCG
```
Where dash symbols represent gaps and vertical bars show where aligned characters match.

We use alignment to make an educated guess as to where a read originated with
respect to the reference genome. It’s not always possible to determine this with
certainty. For instance, if the reference genome contains several long stretches of
As (`AAAAAAAAA` etc.) and the read sequence is a short stretch of As (`AAAAAAA`), we
cannot know for certain exactly where in the sea of As the read originated.

**Examples**
```
# Aligning unpaired reads
bowtie2 -x example/index/lambda_virus -U example/reads/longreads.fq

# Aligning paired reads
bowtie2 -x example/index/lambda_virus -1 example/reads/reads_1.fq -2 example/reads/reads_2.fq
```

### Building an index

`bowtie2-build` builds a Bowtie index from a set of DNA sequences. `bowtie2-build`
outputs a set of 6 files with suffixes `.1.bt2`, `.2.bt2`, `.3.bt2`, `.4.bt2`,
`.rev.1.bt2`, and `.rev.2.bt2`. In the case of a large index these suffixes will
have a `bt2l` termination. These files together constitute the index: they are all
that is needed to align reads to that reference. The original sequence FASTA files
are no longer used by Bowtie 2 once the index is built.

Bowtie 2’s `.bt2` index format is different from Bowtie 1’s `.ebwt` format, and they
are not compatible with each other.

**Examples**
```
# Building a small index
bowtie2-build example/reference/lambda_virus.fa example/index/lambda_virus

# Building a large index
bowtie2-build --large-index example/reference/lambda_virus.fa example/index/lambda_virus
```

### Index inpection

`bowtie2-inspect` extracts information from a Bowtie 2 index about what kind of
index it is and what reference sequences were used to build it. When run without any
options, the tool will output a FASTA file containing the sequences of the original
references (with all non-A/C/G/T characters converted to Ns). It can also be used to
extract just the reference sequence names using the `-n/--names` option or a more
verbose summary using the `-s/--summary` option.

**Examples**
```
# Inspecting a lambda_virus index (small index) and outputting the summary
bowtie2-inspect --summary example/index/lambda_virus

# Inspecting the entire lambda virus index (large index)
bowtie2-inspect --large-index example/index/lambda_virus
```

## Publications

### Bowtie 2 Papers

- Langmead B, Wilks C., Antonescu V., Charles R. __[Scaling read aligners to hundreds of threads on general-purpose processors](https://doi.org/10.1093/bioinformatics/bty648)__. [Bioinformatics](https://academic.oup.com/bioinformatics). bty648.

- Langmead B, Salzberg S. __[Fast gapped-read alignment with Bowtie 2](http://www.nature.com/nmeth/journal/v9/n4/full/nmeth.1923.html)__. [Nature Methods](http://www.nature.com/nmeth). 2012, 9:357-359.

- Langmead B, Trapnell C, Pop M, Salzberg SL. __[Ultrafast and memory-efficient alignment of short DNA sequences to the human genome](http://genomebiology.com/2009/10/3/R25)__. [Genome Biology](http://genomebiology.com/) 10:R25.

### Related Publications

- P. Ferragina, G. Manzini __[Opportunistic data structures with applications](https://ieeexplore.ieee.org/document/892127)__. [IEEE Xplore](http://genomebiology.com/) 10.1109/SFCS.2000.892127

## Related Work

Check out the [Bowtie 2 UI](http://bit.ly/bt2ui-beta), a [shiny](https://shiny.rstudio.com), frontend to the Bowtie 2 command line.
