
./shuffle.py ERR055337_g_shuffle_1.fq ERR055337_g_shuffle_2.fq

# default params
../bowtie2 -x ~/data/hg19 -1 ERR055337_g_shuffle_1.fq.out -2 ERR055337_g_shuffle_2.fq.out --no-unal -S out_$1.sam 2>out_$1.err 

# Use only system time for random generator
#../bowtie2 -x ~/data/hg19 -1 ERR055337_g_shuffle_1.fq.out -2 ERR055337_g_shuffle_2.fq.out --no-unal --non-deterministic -S out_nd_$1.sam 2>out_nd_$1.err 

# Reporting only to check all 3 regions
#../bowtie2 -x ~/data/hg19 -1 ERR055337_g_shuffle_1.fq.out -2 ERR055337_g_shuffle_2.fq.out -k 3 --no-unal -S out_all3_$1.sam 2>out_all_3$1.err 
 
#
