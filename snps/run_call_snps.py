

collections = [line.rstrip() for line in open('collections.lst')]

for collection in collections:
    
    bam = './bams/MISTRG_%s_bq10_star.HUMAN.bam' %(collection)
    dge = './dges/MISTRG_%s_bq10_star.HUMAN.umi.dge.txt.gz' %(collection)
    out = collection
    
    print 'python ~/aviv/box/snps/0.call_snps_from_bam.py --bam %s --dge %s --prefix %s --no_replace' %(bam, dge, out)
