#!/bin/bash

. /broad/tools/scripts/useuse

use Samtools
use Java-1.8
use BWA

fastQ1=$1
fastQ2=$2
outName=$3

set -x

rootDir=$(dirname $(readlink -e `pwd`/))
metaDataDir=/broad/mccarroll/software/metadata/merged_reference/hg19_mm10_transgenes
metaDataName=hg19_mm10_transgenes
refSequence=${metaDataDir}/${metaDataName}.fasta
baseQuality=10
picardLoc=/seq/software/picard/current/bin
gapToolsLoc=/broad/mccarroll/software/dropseq/prod

# IF THE NUMBER OF CORE BARCODES IS SET to a non-zero number, then use that for barcode collapse.  OTHERWISE, use the NUM_READS_CORE to set the min number of reads for a barcode to be considered core.
NUM_CORE_BARCODES=0

#derived variables
genomeDir=${metaDataDir}/STAR
geneIntervals=${metaDataDir}/${metaDataName}.genes.intervals
exonIntervals=${metaDataDir}/${metaDataName}.exons.intervals
refFlat=${metaDataDir}/${metaDataName}.refFlat
tempDir=${rootDir}/temp/star_fast_${baseQuality}_${outName}
outDir=${rootDir}/bams

#probably never changing variables
FivePAdapter=AAGCAGTGGTATCAACGCAGAGTGAATGGG
NUM_READS_CORE=100
MAX_RECORDS_IN_RAM=10000000
java_mem=4g

set -x

mkdir -p ${outDir}
mkdir -p ${tempDir}

java -Xmx${java_mem} -jar ${picardLoc}/picard.jar FastqToSam \
FASTQ=$fastQ1 \
FASTQ2=$fastQ2 \
TMP_DIR=${tempDir} \
QUALITY_FORMAT=Standard \
OUTPUT=/dev/stdout \
COMPRESSION_LEVEL=0 \
MAX_RECORDS_IN_RAM=${MAX_RECORDS_IN_RAM} \
SAMPLE_NAME=${outName} SORT_ORDER=queryname | \
${gapToolsLoc}/TagBamWithReadSequenceExtended \
INPUT=/dev/stdin \
OUTPUT=/dev/stdout \
COMPRESSION_LEVEL=0 \
SUMMARY=${tempDir}/unaligned_tagged1.bam_summary.txt \
BASE_RANGE=13-20 \
BASE_QUALITY=${baseQuality} \
BARCODED_READ=1 \
DISCARD_READ=False \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1 | \
${gapToolsLoc}/TagBamWithReadSequenceExtended \
INPUT=/dev/stdin \
OUTPUT=/dev/stdout \
COMPRESSION_LEVEL=0 \
SUMMARY=${tempDir}/unaligned_tagged3.bam_summary.txt \
BASE_RANGE=1-12 \
BASE_QUALITY=${baseQuality} \
BARCODED_READ=1 \
DISCARD_READ=True \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1 | \
${gapToolsLoc}/FilterBAM \
TAG_REJECT=XQ \
INPUT=/dev/stdin \
OUTPUT=/dev/stdout \
COMPRESSION_LEVEL=0 | \
${gapToolsLoc}/TrimStartingSequence \
INPUT=/dev/stdin \
OUTPUT=/dev/stdout \
COMPRESSION_LEVEL=0 \
OUTPUT_SUMMARY=${tempDir}/adapter_trimming_report.txt \
SEQUENCE=${FivePAdapter} \
MISMATCHES=0 \
NUM_BASES=5 | \
${gapToolsLoc}/PolyATrimmer \
INPUT=/dev/stdin \
OUTPUT=${tempDir}/unaligned_mc_tagged_polyA_filtered.bam \
OUTPUT_SUMMARY=${tempDir}/polyA_trimming_report.txt \
MISMATCHES=0 \
NUM_BASES=6

java -jar -Xmx2g ${picardLoc}/picard.jar SamToFastq \
INPUT=${tempDir}/unaligned_mc_tagged_polyA_filtered.bam \
FASTQ=${tempDir}/unaligned_tagged.fastq.gz \
tmp_dir=${tempDir} CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2

#align with star
/fg/software/gap/gap_analysis/3rdParty/STAR_2.4.0a/STAR \
--genomeDir ${genomeDir} \
--readFilesIn <(gunzip -c ${tempDir}/unaligned_tagged.fastq.gz) \
--outFileNamePrefix ${tempDir}/star

java -Djava.io.tmpdir=${tempDir} -Xmx${java_mem} -jar ${picardLoc}/picard.jar SortSam \
I=${tempDir}/starAligned.out.sam \
O=${tempDir}/starAligned.out.bam \
SO=queryname

rm ${tempDir}/starAligned.out.sam

java -Djava.io.tmpdir=${tempDir} -Xmx${java_mem} -jar ${picardLoc}/picard.jar MergeBamAlignment \
REFERENCE_SEQUENCE=${refSequence} \
UNMAPPED_BAM=${tempDir}/unaligned_mc_tagged_polyA_filtered.bam \
ALIGNED_BAM=${tempDir}/starAligned.out.bam \
MAX_RECORDS_IN_RAM=${MAX_RECORDS_IN_RAM} \
OUTPUT=/dev/stdout \
COMPRESSION_LEVEL=0 \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false VALIDATION_STRINGENCY=SILENT | \
${gapToolsLoc}/TagReadWithInterval \
I=/dev/stdin \
O=/dev/stdout \
COMPRESSION_LEVEL=0 \
LOCI=${geneIntervals} \
TAG=XG | \
${gapToolsLoc}/TagReadWithInterval \
I=/dev/stdin \
O=/dev/stdout \
COMPRESSION_LEVEL=0 \
LOCI=${exonIntervals} \
TAG=XE | \
${gapToolsLoc}/TagReadWithGeneExon \
I=/dev/stdin \
O=${tempDir}/${outName}_star_gene_exon_tagged2.bam \
ANNOTATIONS_FILE=${refFlat} \
TAG=GE



#if [[ $NUM_CORE_BARCODES -eq 0 ]]; then
#${gapToolsLoc}/CollapseBarcodesInPlace \
#INPUT=${tempDir}/${outName}_star_gene_exon_tagged2.bam \
#OUTPUT=${tempDir}/${outName}_star_gene_exon_tagged_cell.bam \
#PRIMARY_BARCODE=XC \
#OUT_BARCODE=ZC \
#READ_QUALITY=10 \
#FILTER_PCR_DUPLICATES=false \
#FIND_INDELS=true \
#MAX_RECORDS_IN_RAM=${MAX_RECORDS_IN_RAM} \
#MIN_NUM_READS_CORE=${NUM_READS_CORE} \
#NUM_THREADS=${num_threads} \
#MIN_NUM_READS_NONCORE=10
#else 
#${gapToolsLoc}/CollapseBarcodesInPlace \
#INPUT=${tempDir}/${outName}_star_gene_exon_tagged2.bam \
#OUTPUT=${tempDir}/${outName}_star_gene_exon_tagged_cell.bam \
#PRIMARY_BARCODE=XC \
#OUT_BARCODE=ZC \
#READ_QUALITY=10 \
#FILTER_PCR_DUPLICATES=false \
#FIND_INDELS=true \
#MAX_RECORDS_IN_RAM=${MAX_RECORDS_IN_RAM} \
#NUM_CORE_BARCODES=${NUM_CORE_BARCODES} \
#NUM_THREADS=${num_threads} \
#MIN_NUM_READS_NONCORE=10  
#fi

#cp ${tempDir}/${outName}_star_gene_exon_tagged_cell.bam ${outDir}/${outName}_bq${baseQuality}_star.bam
cp ${tempDir}/${outName}_star_gene_exon_tagged2.bam ${outDir}/${outName}_bq${baseQuality}_star.bam
samtools index ${outDir}/${outName}_bq${baseQuality}_star.bam
