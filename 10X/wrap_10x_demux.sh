#! /bin/bash

#$ -cwd
#$ -q long
#$ -P regevlab

#$ -l h_vmem=32g
#$ -e demux.err
#$ -o demux.log

source /broad/software/scripts/useuse 
reuse -q .bcl2fastq2-2.17.1.14
export PATH=/seq/regev_genome_portal/SOFTWARE/10X/cellranger-VERSION:$PATH
cellranger demux --run=INDIR --maxjobs=10 --mempercore=16
