#!/bin/bash

refFastaPath=/broad/mccarroll/software/metadata/merged_reference/hg19_mm10_transgenes/hg19_mm10_transgenes.fasta
numCells=(3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000)

rm -rf ../Analysis ../qsub_logs ../bam* ../*DGE ../temp* ../QC*
mkdir -p ../qsub_logs
mkdir -p ../Analysis
mkdir -p ../QC_files
mkdir -p ../QC_reports
mkdir -p ../tempQC

l=0
for fq1 in ../Data/*R1*;
do
#STEP 1 : ALIGNMENT TO REFERENCE (CHECK CODE TO ENSURE THAT LEFT AND RIGHT FASTQ FILES ARE SUPPLIED CORRECTLY)
fq2=${fq1/R1/R2}
echo ${fq1}
echo ${fq2}

bfq1=`basename ${fq1}`
bfq2=`basename ${fq2}`

b=`echo ${bfq1} | grep -P '^[a-zA-Z0-9\_]*_S' -o`
b=${b/_S/}
echo $b

sed "s|fName|${b}|g;s|fastq1|${bfq1}|g;s|fastq2|${bfq2}|g;s|bamFileName|${b}_bq10_star|g;s|numCellsNum|${numCells[l]}|g;s|refFasta|${refFastaPath}|g" < run_Alignment.sh > run_Alignment_${b}.sh
chmod +x run_Alignment_${b}.sh
#STEP 2 : BAMTagHistogram 
mkdir -p ../bam_reads
#STEP 3 : Cumulative Plot
Rfile=DropSeqCumuPlot_${b}.R
bamReadsFile=../bam_reads/${b}_bq10_star.reads.txt
temp1=`echo ${bamReadsFile} | sed 's:\/:\\\/:g'`
temp2=..\\\/bam_reads\\\/${b}
sed "s/fileName/${temp1}/g;s/figName/${temp2}/g" < DropSeqCumuPlot.R > DropSeqCumuPlot_${b}.R

#STEP 4 : Collect Cell Barcodes

#STEP 5 : Split Bam files
mkdir -p ../bams_HUMAN_MOUSE

#STEP 6 : DGE
mkdir -p ../UMI_DGE
mkdir -p ../reads_DGE

ssub -q long -m 150 -o ../qsub_logs/${b}.log  "./run_Alignment_${b}.sh"

echo ${numCells[l]}
l=`expr $l + 1`
done
