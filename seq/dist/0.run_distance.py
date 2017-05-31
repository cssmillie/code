import pandas as pd
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fst', help='Input fasta file', default='')
parser.add_argument('--aln', help='Input alignment', default='')
parser.add_argument('--subset', help='Subset (list)', default='')
parser.add_argument('--cpus', help='Number of CPUs', required=True, type=int)
parser.add_argument('--prefix', help='Output prefix', default='D')
args = parser.parse_args()

if args.subset != '':
    len_otus = len(open(args.subset).readlines())
else:
    len_otus = 0
    fiche = args.fst
    if not fiche:
        fiche = args.aln
    for line in open(fiche):
        if line.startswith('>'):
            len_otus += 1

size = int((2*args.cpus)**.5)
basename = os.getcwd()

for i in range(size):
    for j in range(i+1):
        
        beg1 = (1.*len_otus/size)*i
        end1 = (1.*len_otus/size)*(i+1)

        beg2 = (1.*len_otus/size)*j
        end2 = (1.*len_otus/size)*(j+1)

        ooo_fn = '%s/%s.%d_%d_%d_%d.dist' %(basename, args.prefix, beg1, end1, beg2, end2)
        if os.path.exists(ooo_fn):
            continue

        cmd = 'python /home/csmillie/sbin/dist/get_distance.py --beg1 %d --end1 %d --beg2 %d --end2 %d' %(beg1, end1, beg2, end2)
        if args.fst != '':
            cmd += ' --fst %s' %(args.fst)
        if args.aln != '':
            cmd += ' --aln %s' %(args.aln)
        if args.subset != '':
            cmd += ' --subset %s' %(args.subset)
        cmd += ' > %s' %(ooo_fn)
        print cmd
