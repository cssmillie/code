import glob

for i in range(100):
    fns = glob.glob('interval%d.*.g.vcf' %(i))
    v = ' '.join(['--variant %s' %(fn) for fn in fns])
    cmd = 'java -jar ~/bin/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /broad/mccarroll/software/metadata/merged_reference/hg19_mm10_transgenes/hg19_mm10_transgenes.fasta -o interval%d.vcf %s' %(i, v)
    print cmd
