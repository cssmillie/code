v = ['-V interval%d.vcf' %(i) for i in range(100)]
print 'java -cp ~/bin/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R /broad/mccarroll/software/metadata/merged_reference/hg19_mm10_transgenes/hg19_mm10_transgenes.fasta -out all.vcf ' + ' '.join(v)
