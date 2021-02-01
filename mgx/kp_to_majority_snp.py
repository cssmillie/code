import numpy as np
import argparse, pickle, glob, re, sys

# Read input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--kpileup', help='list of kpileup files (newline-separated)')
parser.add_argument('--gene_file', help='kpileup gene file (.txt)')
args = parser.parse_args()

# Read kpileup gene file
fst = {}
for line in open(args.gene_file):
    line = line.rstrip().split()
    contig = line[0]
    gene = line[1]
    beg = line[2]
    end = line[3]
    seq = line[5]
    fst[contig] = seq

# Read kpileup output file
for line in open(args.kpileup):
    line = line.rstrip().split()
    if len(line) == 10 and line[0] != 'Sample':
        sample = line[0]
        




    # Add kpileup results to numpy arrays
for sample in sample2index:
    for line in open('%s.kp.txt' %(sample)):
        line = line.rstrip().split()
        if len(line) == 10 and line[0] != 'Sample':
            sample = line[0]
            i = sample2index[sample]
            contig = line[1]
            j = int(line[2])
            nt = line[7]
            k = nts.index(nt)
            count = int(line[8])
            x[contig][i,j-1,k] = count

# Write numpy arrays to file
pickle.dump(x, open(args.out, 'w'))
