import numpy as np
from scipy import sparse
import argparse, pickle, glob, re, sys

# Merge kpileups from multiple samples. Write dictionary of (M, N, 4) numpy arrays where:
# M = samples
# N = alignment sites
# 4 = nucleotides (ACGT)
# Entry (i,j,k) of this array corresponds to the count of nucleotide k at position j of sample i

# Read input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--samples', help='sample list (newline-delimited)')
parser.add_argument('--gene_file', help='kpileup gene file')
parser.add_argument('--out', help='output file (.pickle)')
args = parser.parse_args()

# Map samples to indices
sample2index = {}
i = 0
for line in open(args.samples):
    sample = line.rstrip()
    sample2index[sample] = i
    i += 1
M = len(sample2index)

# Initialize numpy arrays for each genome
gene2coords = {}
x = {'A':{}, 'C':{}, 'G':{}, 'T':{}}
for line in open(args.gene_file):
    line = line.rstrip().split()
    contig = line[0]
    gene = line[1]
    beg = int(line[2])
    end = int(line[3])
    gene2coords[gene] = [beg, end]
    for nt in ['A', 'C', 'G', 'T']:
        x[nt][gene] = sparse.csr_matrix((M, (end - beg + 1)))
        print(np.shape(x[nt][gene]))
        
# Add kpileup results to numpy arrays
for sample in sample2index:
    for line in open('%s.kp.txt' %(sample)):
        line = line.rstrip().split()
        if len(line) == 10 and line[0] != 'Sample':
            sample = line[0]
            i = sample2index[sample]
            contig = line[1]
            pos = int(line[2])
            gene = line[3]
            nt = line[7]
            count = int(line[8])
            j = pos - gene2coords[gene][0]
            x[nt][gene][i,j] = count

# Write numpy arrays to file
pickle.dump(x, open(args.out, 'w'))
