# Use HaplotypeCaller to call SNPs in a BAM file

import argparse, os, os.path, ssub, sys


# Command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--bam', help='bam file', required=True)
parser.add_argument('--fst', help='fasta file', default='/broad/mccarroll/software/metadata/merged_reference/hg19_mm10_transgenes/hg19_mm10_transgenes.fasta')
parser.add_argument('--dge', help='dge file', required=True)
parser.add_argument('--btag', help='barcode tag in bam file', default='ZC:Z')
parser.add_argument('--utag', help='umi tag in bam file', default='XM:Z')
parser.add_argument('--gtag', help='gene tag in bam file', default='XG:Z')
parser.add_argument('--prefix', help='output prefix', required=True)
parser.add_argument('--no_replace', help='skip finished vcf files', default=False, action='store_true')
args = parser.parse_args()


# Job submitter
submitter = ssub.Submitter()

# Set filenames
base_dir = os.path.dirname(os.path.realpath(__file__))
barcodes_fn = '%s.barcodes.out' %(args.prefix)
sam_fn = '%s.sam' %(args.prefix)
bam_fn = '%s.bam' %(args.prefix)


# Get good barcodes from DGE (> 300 genes)
cmd = 'Rscript %s/1.get_barcodes.r %s %s' %(base_dir, args.dge, barcodes_fn)
submitter.submit_and_wait(cmd)
#print cmd


# Convert BAM to SAM
# 1. Add @RG header
# 2. Remove reads with bad barcodes
cmd = 'samtools view -h %s | %s/2.process_bam.py --barcodes %s --btag %s --utag %s --gtag %s --prefix %s > %s' %(args.bam, base_dir, barcodes_fn, args.btag, args.utag, args.gtag, args.prefix, sam_fn)
submitter.submit_and_wait(cmd)
#print cmd


# Convert back to BAM
cmd = 'samtools view -b -S %s > %s; rm %s' %(sam_fn, bam_fn, sam_fn)
submitter.submit_and_wait(cmd)
#print cmd


# Make BAM index
cmd = 'samtools index %s' %(bam_fn)
submitter.submit_and_wait(cmd)
#print cmd

# Create folders for analysis
if not os.path.exists('./vcf'):
    os.mkdir('./vcf')
if not os.path.exists('./vcf/%s' %(args.prefix)):
    os.mkdir('./vcf/%s' %(args.prefix))

# Call variants separately for each sample using the -ERC GVCF option
barcodes = [line.rstrip() for line in open(barcodes_fn)]
cmds = []

for barcode in barcodes:
    sample = '%s_%s' %(args.prefix, barcode)
    if not os.path.exists('./vcf/%s/%s' %(args.prefix, barcode)):
        os.mkdir('./vcf/%s/%s' %(args.prefix, barcode))
    # Split into 100 separate intervals to speed things up
    for i in range(100):
        vcf_fn = './vcf/%s/%s/%s.%d.g.vcf' %(args.prefix, barcode, sample, i)
        if args.no_replace == True:
            idx_fn = vcf_fn + '.idx'
            if os.path.exists(idx_fn):
                continue
        int_fn = '%s/intervals/hg19.exons.%d.intervals' %(base_dir, i)
        cmd = '%s/3.haplotype_caller.sh %s %s %s %s %s' %(base_dir, bam_fn, args.fst, sample, vcf_fn, int_fn)
        cmds.append(cmd)
submitter.submit(cmds)
#print cmds

'''
# Get list of GVCFs to combine (list of lists)
fns_list = []
for i in range(100):
    fns = []
    for barcode in barcodes:
        sample = '%s_%s' %(args.prefix, barcode)
        vcf_fn = './vcf/%s/%s/%s.%d.g.vcf' %(args.prefix, barcode, sample, i)
        fns.append(vcf_fn)
    fns_list.append(fns)

# Hierarchically combine GVCFs
# TEST - NEED TO BREAK INTO SEPARATE SCRIPT
iteration = 0
while True:
    iteration += 1
    if sum([len(fns) for fns in fns_list]) == 0:
        break
    cmds = []
    for i, fns in enumerate(fns_list):
        count = 0
        if len(fns) <= 1:
            fns_list[i] = []
        else:
            new_fns = []
            while len(fns) > 0:
                count += 1
                j = min(200, len(fns))
                v = ' '.join(['-V %s' %(fn) for fn in fns[:j]])
                fns = fns[j:]
                new_fn = './vcf/region%d.iter%d.%d.combine.g.vcf' %(i, iteration, count)
                new_fns.append(new_fn)
                cmd = 'java -jar ~/bin/GenomeAnalysisTK.jar -T CombineGVCFs -R %s %s -o %s' %(args.fst, v, new_fn)
                cmds.append(cmd)
                fns_list[i] = new_fns
    if len(cmds) > 0:
        submitter.submit_and_wait(cmds)
'''
