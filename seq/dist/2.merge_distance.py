import glob, sys, os.path, re
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fst', required=True, help='Fasta file')
parser.add_argument('--subset', default='', help='Subset file (list)')
parser.add_argument('--prefix', required=True, help='Output prefix')
parser.add_argument('--out', required=True, help='Outfile (distance matrix)')
args = parser.parse_args()

fst_fn = args.fst
prefix = args.prefix
out_fn = args.out

if args.subset != '':
    otus = [line.rstrip() for line in open(args.subset).readlines()]
else:
    otus = []
    for line in open(fst_fn):
        if line.startswith('>'):
            otus.append(line.rstrip().split()[0][1:])
otus = sorted(otus)
otus = [otu.split(';')[0] for otu in otus]


basename = os.getcwd()
corr_fns = glob.glob('%s/%s.*.dist' %(basename, prefix))

R = pd.DataFrame(index=otus, columns=otus)

for corr_fn in corr_fns:

    correlations = open(corr_fn).readline().rstrip().split()

    coords = re.search('.*\/%s.(.*?).dist' %(prefix),corr_fn).group(1)
    [beg1, end1, beg2, end2] = [int(xi) for xi in coords.split('_')]

    for i in range(beg1, end1):
        if i<0 or i>len(otus):
            continue
        k = 0
        for j in range(beg2, end2):
            if j<0 or j>=i:
                continue
            k += 1
        v = correlations[:k][:]
        correlations[:k] = []
        R.iloc[i,beg2:beg2+k] = R.iloc[beg2:beg2+k,i] = v
    
for i in range(len(R.index)):
    R.iloc[i,i] = 0

R.to_csv(out_fn, sep='\t')
