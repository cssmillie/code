#! /bin/bash

#$ -cwd
#$ -q long
#$ -P regevlab
#$ -l h_vmem=4g
#$ -e demux.err
#$ -o demux.log

source /broad/software/scripts/useuse 
reuse -q .bcl2fastq2-2.17.1.14
reuse UGER
export PATH=/seq/regev_genome_portal/SOFTWARE/10X/cellranger-VERSION:$PATH
cellranger mkfastq --run=INDIR --csv CSV --jobmode=/home/unix/csmillie/code/10X/sge.template
