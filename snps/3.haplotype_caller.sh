use Java-1.7

bam=$1
ref=$2
sid=$3
out=$4
ivl=$5

java -jar ~/bin/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref -I $bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --filter_reads_with_N_cigar --emitRefConfidence GVCF --sample_name $sid -o $out -L $ivl
