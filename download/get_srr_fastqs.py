import os, sys

# steps
# 1. download srr ids ("Select All" in Run Selector -> copy URL)

# input arguments
fn = sys.argv[1]

# get FASTQ urls
for line in open(fn):
    for srr in line.rstrip().split(','):
        if not os.path.exists('%s.fastq.gz' %(srr)):
            url = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/00%s/%s/*' %(srr[:6], srr[-1], srr)
            cmd = 'wget "%s"' %(url)
            print(cmd)
