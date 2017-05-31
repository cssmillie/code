#!/bin/sh

#$ -cwd
#$ -q long
#$ -P regevlab

#$ -m be
#$ -N run10X
#$ -t 1-4
#$ -l h_vmem=8g
#$ -e run.error
#$ -o run.log

source /broad/software/scripts/useuse
reuse UGER
export PATH=/seq/regev_genome_portal/SOFTWARE/10X/cellranger-VERSION:$PATH

SEED=$(awk "NR==$SGE_TASK_ID" BARCODES)

channel_id=$(echo "$SEED" | cut -d$'\t' -f1)
channel_indices=$(echo "$SEED" | cut -d$'\t' -f2)

echo $channel_id
echo $channel_indices

cellranger count --id=${channel_id} \
               --transcriptome=REF \
               --fastqs=FQPATH \
               --jobmode=sge \
               --maxjobs=10 \
               --mempercore=16 \
               --jobinterval=2500 \
               --indices=${channel_indices}
