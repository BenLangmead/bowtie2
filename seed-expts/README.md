Rough order of operations:
* Download the indexes by changing to `indexes` subdir and running `get.sh`
* Download the SAM sources for the reads by changing to `reads` subdir and running `get.sh`
* Convert the SAMs to FASTQs by running `./reads.sh` in this directory
* Truncate the reads by running `truncate.bash` in `reads` subdir
* Run `aligner.sh` to run the suite of experiments
    * In the outermost loop, it will build and evaluate two different versions of `bowtie2-aligns`; the versions are described by patch files in the subdirectories of the `versions` directory
    * In the next nested loop, it iterates over `local` and `semiglobal`, since our Vargas-annotated reads are from both kinds of `bowtie2` runs
    * In the next nested loop, it iterates over all the reads files, which are distinguished by which index should be aligned to.  That mapping is recorded in `genome_index.txt`. 
* This produces a suite of `.csv` files that can be analyzed with `adseed.Rmd`
    * The plots need work

Inventory of scripts, not all of which have been needed/used lately:

* `vassess_aligner.py`/`aligner.sh`: described above
* `make_version.sh`: helper script for `vassess_aligner.py` that cds up a directory, applies appropriate patch, builds `bowtie2-align-s`, copies it to subdir, undoes the patch.  Assumes this working dir is at HEAD of adseed without staged changes.
* `vassess_reads.py`/`reads.sh`: these scripts work to convert Vargas SAM records into FASTQ files where the correct alignment scores are recorded in the read name
* `combine.py`: similar to previous, but combines a Vargas SAM with an aligner SAM
* `compare.py`/`compare_all.sh`: work together to search for "obvious winners" among all the configurations tried 
* `adseed.Rmd`: for data analysis of the csvs that come from running `aligner.sh`

Inventory of other files:

* `genome_index.txt`: list of names of genomes used for alignment, along with path to index files for them.  Used by `align.sh`
* `cmds_local.txt`: list of all the parameter-setting configurations to use for testing aligner's local mode
* `cmds_semiglobal.txt`: list of all the parameter-setting configurations to use for testing aligner's semi-global (end-to-end) mode
