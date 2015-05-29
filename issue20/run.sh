
./shuffle.py ERR055337_g_shuffle_1.fq ERR055337_g_shuffle_2.fq

# default params
#../bowtie2 -x ~/data/hg19 -1 ERR055337_g_shuffle_1.fq.out -2 ERR055337_g_shuffle_2.fq.out --no-unal -S out_$1.sam 2>out_$1.err 

#for (( i=1; i<18; i++)); do
#  ../bowtie2 -x ~/data/hg19 -1 ERR055337_g_shuffle_1.fq.out -2 ERR055337_g_shuffle_2.fq.out --no-unal -S out_$1_p$i.sam 2>out_$1_p$i.err 
#done

# Use only system time for random generator
#../bowtie2 -x ~/data/hg19 -1 ERR055337_g_shuffle_1.fq.out -2 ERR055337_g_shuffle_2.fq.out --no-unal --non-deterministic -S out_nd_$1.sam 2>out_nd_$1.err 

#for (( i=1; i<18; i++)); do
#  ../bowtie2 -x ~/data/hg19 -1 ERR055337_g_shuffle_1.fq.out -2 ERR055337_g_shuffle_2.fq.out --no-unal --non-deterministic -S out_nd_$1_p$i.sam 2>out_nd_$1_p$i.err 
#done

# Reporting only to check all 3 regions
#../bowtie2 -x ~/data/hg19 -1 ERR055337_g_shuffle_1.fq.out -2 ERR055337_g_shuffle_2.fq.out -k 3 --no-unal 2>out_all_3$1.err | samtools view -bS - -o out_all_3$1.bam 
 

# For igm machines
#mkdir -p default nondeterministic
#for (( i=1; i<18; i++)); do
#  ../bowtie2 -x ~/work/hg19/hg19 -1 ERR055337_g_shuffle_1.fq.out -2 ERR055337_g_shuffle_2.fq.out --no-unal -S default/out_$1_p$i.sam 2>default/out_$1_p$i.err 
#done

# Use only system time for random generator
#../bowtie2 -x ~/data/hg19 -1 ERR055337_g_shuffle_1.fq.out -2 ERR055337_g_shuffle_2.fq.out --no-unal --non-deterministic -S out_nd_$1.sam 2>out_nd_$1.err 

#for (( i=1; i<18; i++)); do
#  ../bowtie2 -x ~/work/hg19/hg19 -1 ERR055337_g_shuffle_1.fq.out -2 ERR055337_g_shuffle_2.fq.out --no-unal --non-deterministic -S nondeterministic/out_nd_$1_p$i.sam 2>nondeterministic/out_nd_$1_p$i.err
#done

