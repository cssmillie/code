#!/usr/bin/env python

import argparse, re, sys

# get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--barcodes', help='list of barcodes', required=True)
parser.add_argument('--btag', help='barcode tag', default='ZC')
parser.add_argument('--gtag', help='gene tag', default='XG')
parser.add_argument('--utag', help='umi tag', default='XM')
parser.add_argument('--prefix', help='sample prefix', required=True)
args = parser.parse_args()

# compile regular expressions
barcode_regex = re.compile('%s:(\S*)' %(args.btag))
gene_regex = re.compile('%s:(\S*)' %(args.gtag))
umi_regex = re.compile('%s:(\S*)' %(args.utag))

# initialize variables
bcode2rgz = {} # map barcode to rgz number
data = {} # store barcode-gene-umi data
match = 0 # count regex matches
total = 0 # count total lines

# write RG:Z header
barcodes = [line.rstrip() for line in open(args.barcodes)]

# header flag
header = 1

# modify each line of stdin
for line in sys.stdin:
    
    total += 1
    
    # write header
    if line.startswith('@RG'):
        if header == 1:
            for i, barcode in enumerate(barcodes):
                print '@RG\tID:%d\tSM:%s_%s' %(i, args.prefix, barcode)
            header = 0
        continue
    if line.startswith('@'):
        print line.rstrip()
        continue
    
    # parse line in samfile
    try:
        barcode = barcode_regex.search(line).group(1)
        umi = umi_regex.search(line).group(1)
        gene = gene_regex.search(line).group(1)
    except:
        if 'INTERGENIC' not in line:
            print line.rstrip()
        continue
    
    match += 1
    
    # get barcode number
    if barcode not in barcodes:
        continue
    if barcode not in bcode2rgz:
        bcode2rgz[barcode] = len(bcode2rgz)
    rgz = bcode2rgz[barcode]
    
    # skip previously seen umis
    if barcode not in data:
        data[barcode] = {}
    if gene not in data[barcode]:
        data[barcode][gene] = {}
    if umi in data[barcode][gene]:
        continue
    else:
        data[barcode][gene][umi] = 1

    # replace rgz
    line = re.sub('RG:Z:A', 'RG:Z:%d' %(rgz), line)
    
    print line.rstrip()
