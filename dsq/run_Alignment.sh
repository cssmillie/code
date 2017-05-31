#!/bin/bash

source ~/.bashrc

b=fName
bamName=bamFileName
fq1=../Data/fastq1
fq2=../Data/fastq2
numCells=numCellsNum
reference_fasta=refFasta

#STEP 1 : Alignment
./align_star_human_mouse.sh ${fq1} ${fq2} ${b}

mv ../bams/${bamName}.bam ../bams/${bamName}_old.bam
mv ../bams/${bamName}.bam.bai ../bams/${bamName}_old.bam.bai

#STEP 1.1 : Detect Bead Synthesis errors
/broad/mccarroll/software/dropseq/prod/DetectBeadSynthesisErrors I=../bams/${bamName}_old.bam O=../bams/${bamName}_unmerged.bam OUTPUT_STATS=../bam_reads/${bamName}.synthesis_stats.txt SUMMARY=../bam_reads/${bamName}.synthesis_stats.summary.txt NUM_BARCODES=$((numCells*2)) PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC MAX_NUM_ERRORS=2 CELL_BARCODE_TAG=XC

samtools index ../bams/${bamName}_unmerged.bam

/broad/mccarroll/software/dropseq/prod/CollapseBarcodesInPlace I=../bams/${bamName}_unmerged.bam O=../bams/${bamName}.bam PRIMARY_BARCODE=XC OUT_BARCODE=ZC MIN_NUM_READS_CORE=3000 MIN_NUM_READS_NONCORE=30000000 EDIT_DISTANCE=2

rm ../bams/${bamName}.bam.bai
samtools index ../bams/${bamName}.bam

rm ../bams/${bamName}_old.bam*
rm ../bams/${bamName}_unmerged.bam*

#STEP 2: Bam Tag Histogram
/broad/mccarroll/software/dropseq/prod/BAMTagHistogram I=../bams/${bamName}.bam O=../bam_reads/${bamName}.reads.txt.gz TAG=ZC
#/broad/mccarroll/software/dropseq/prod/BAMTagHistogram I=../bams/${bamName}_old.bam O=../bam_reads/${bamName}_old.reads.txt.gz TAG=XC
#/broad/mccarroll/software/dropseq/prod/BAMTagHistogram I=../bams/${bamName}_unmerged.bam O=../bam_reads/${bamName}_unmerged.reads.txt.gz TAG=XC

gunzip ../bam_reads/${bamName}.reads.txt.gz
#gunzip ../bam_reads/${bamName}_old.reads.txt.gz
#gunzip ../bam_reads/${bamName}_unmerged.reads.txt.gz

#STEP 3: Cumulative plot to compute cell barcodes
R CMD BATCH DropSeqCumuPlot_${b}.R ../qsub_logs/DropSeqCumuPlot_${b}.out
readsTable=${bamName}.reads.txt
numCells=`cat ../bam_reads/${b}_numCells.txt` 
sed "s/\<filename_input\>/${readsTable}/g;s/\<Ncells_input\>/${numCells}/g" < collect_cell_barcodes.R > collect_cell_barcodes.${b}.R
R CMD BATCH collect_cell_barcodes.${b}.R

#STEP 4 : Segregate Human and Mouse BAM files
/broad/mccarroll/software/dropseq/prod/FilterBAM I=../bams/${bamName}.bam O=../bams_HUMAN_MOUSE/${bamName}.HUMAN.bam REF_SOFT_MATCHED_RETAINED=HUMAN
/broad/mccarroll/software/dropseq/prod/FilterBAM I=../bams/${bamName}.bam O=../bams_HUMAN_MOUSE/${bamName}.MOUSE.bam REF_SOFT_MATCHED_RETAINED=MOUSE

#STEP 5: DGE UMIs
/broad/mccarroll/software/dropseq/prod/DigitalExpression I=../bams_HUMAN_MOUSE/${bamName}.HUMAN.bam O=../UMI_DGE/${bamName}.HUMAN.umi.dge.txt.gz SUMMARY=../UMI_DGE/${bamName}.HUMAN.umi.dge.summary.txt CELL_BARCODE_TAG=ZC CELL_BC_FILE=../bam_reads/${bamName}.barcodes_use.txt
/broad/mccarroll/software/dropseq/prod/DigitalExpression I=../bams_HUMAN_MOUSE/${bamName}.MOUSE.bam O=../UMI_DGE/${bamName}.MOUSE.umi.dge.txt.gz SUMMARY=../UMI_DGE/${bamName}.MOUSE.umi.dge.summary.txt CELL_BARCODE_TAG=ZC CELL_BC_FILE=../bam_reads/${bamName}.barcodes_use.txt

#STEP 6: DGE READS
/broad/mccarroll/software/dropseq/prod/DigitalExpression I=../bams_HUMAN_MOUSE/${bamName}.HUMAN.bam O=../reads_DGE/${bamName}.HUMAN.reads.dge.txt.gz SUMMARY=../reads_DGE/${bamName}.HUMAN.reads.dge.summary.txt CELL_BARCODE_TAG=ZC CELL_BC_FILE=../bam_reads/${bamName}.barcodes_use.txt OUTPUT_READS_INSTEAD=true
/broad/mccarroll/software/dropseq/prod/DigitalExpression I=../bams_HUMAN_MOUSE/${bamName}.MOUSE.bam O=../reads_DGE/${bamName}.MOUSE.reads.dge.txt.gz SUMMARY=../reads_DGE/${bamName}.MOUSE.reads.dge.summary.txt CELL_BARCODE_TAG=ZC CELL_BC_FILE=../bam_reads/${bamName}.barcodes_use.txt OUTPUT_READS_INSTEAD=true

#STEP 7: QC
BASEDIR=/broad/mccarroll/software/dropseq/prod
source $BASEDIR/configDropSeqRNAEnvironment.bash

$BASEDIR/DropSeqStandardAnalysis \
--BAMFile ../bams/${bamName}.bam \
--reference $reference_fasta \
--numCells ${numCells} \
--estimatedNumCells ${numCells} \
--estimatedNumBeads $((numCells*20)) \
--report_dir ../QC_files \
--pointSize=0.75 \
--batchSystem local \
--beadSynthesisErrorDetail ../bam_reads/${bamName}.synthesis_stats.txt \
--cellTag XC \
--cellTagCollapsed ZC \
--outPDF ../QC_reports/${bamName}_QC.pdf \
--tempDir ../tempQC \
--use_threads \
--verbose 1 

