# aliases for converting sample read files
# `source conversion_utilities.sh` to add them to your environment

alias fastq_to_fasta="sed 'N;x;N;N;x;s/@/>/"
alias paired_to_tab5="paste <(sed 'N;x;N;g;N;s/\n/	/g' reads_1.fq) <(sed  -n 'n;h;n;g;N;s/\n/	/g;p' reads_2.fq) > reads_12.tab5"
alias paired_to_tab6="paste <(sed 'N;x;N;g;N;s/\n/	/g' reads_1.fq) <(sed 'N;x;N;g;N;s/\n/	/g' reads_2.fq) > reads_12.tab6"
alias paired_to_interleaved="paste -d'\n' <(sed 'N;N;N;s/\n/	/g' reads_1.fq) <(sed 'N;N;N;s/\n/	/g' reads_2.fq) | tr '\t' '\n' > reads_12.fq"
