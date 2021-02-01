
FR1=$1 # FASTQ file (R1)
FR2=$2 # FASTQ file (R2)
BDB=$3 # bowtie database
SAM=$4 # SAM file
GTF=$5 # GTF file
OUT=$6 # output file

# run bowtie
bowtie2 -x $BDB --very-sensitive -a --no-unal -1 $FR1 -2 $FR2 -S $SAM

# sort SAM by position
samtools sort $SAM -O sam > ${SAM%.*}.sort.sam

# count CDS in sorted SAM
htseq-count -r pos -s no -a 0 -t CDS ${SAM%.*}.sort.sam $GTF > $OUT
