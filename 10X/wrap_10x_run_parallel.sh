#!/bin/sh
#$ -cwd
#$ -P regevlab
#$ -N run10X
#$ -l h_vmem=4g
#$ -e run.error
#$ -o run.log
#$ -l h_rt=30:00:00
#$ -t 1-NUM
#$ -l os=RedHat7

source /broad/software/scripts/useuse
reuse UGER
use .cellranger-3.0.2
#export PATH=/seq/regev_genome_portal/SOFTWARE/10X/cellranger-VERSION:$PATH

SEED=$(awk "NR==($SGE_TASK_ID+1)" CSV)
channel_lane=$(echo "$SEED" | cut -d$',' -f1)
channel_id=$(echo "$SEED" | cut -d$',' -f2)
channel_indices=$(echo "$SEED" | cut -d$',' -f3)

command="
cellranger count --id=${channel_id} \
	   --transcriptome=REF \
	   --fastqs=FASTQ_PATH \
	   --jobmode=/home/unix/csmillie/code/10X/sge.template \
	   --mempercore=8 \
	   --nosecondary \
	   --force-cells=6000
"

# Add lanes
if [[ $channel_lane  != '*' ]]; then
    # Get comma-separated lanes
    lane1=$(echo "$channel_lane" | cut -d$'-' -f 1)
    lane2=$(echo "$channel_lane" | cut -d$'-' -f 2)
    lanes=$(seq -s, $lane1 $lane2)
    #command="${command} --lanes=${channel_lane}"
    command="${command} --lanes=${lanes}"
fi

# If demultiplexed with --demux, add indices
if [[ -e ./SERIAL/BCL_PROCESSOR_CS ]]; then
    command="${command} --indices=${channel_indices}"
# If demultiplexed with --mkfastq, add sample
elif [[ -e ./SERIAL/MAKE_FASTQS_CS ]]; then
    command="${command} --sample=${channel_id}"
else
    exit 1
fi

${command}
