import numpy as np
import argparse, pickle, glob, os, re, sys, time

# Merge kpileups from multiple samples. Write dictionary of (M, N, 4) numpy arrays where:
# M = samples
# N = alignment sites
# 4 = nucleotides (ACGT)
# Entry (i,j,k) of this array corresponds to the count of nucleotide k at position j of sample i

# Read input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--kp_dir', help='Directory for kpileup files', default='.')
parser.add_argument('--samples', help='Sample list (newline-delimited)')
parser.add_argument('--gene', help='Gene name')
parser.add_argument('--gene_file', help='kpileup gene file')
parser.add_argument('--out', help='Output file (.pickle)')
args = parser.parse_args()

if os.path.exists(args.out):
    quit()

# Initialize data
nts = 'ACGT'

# Map samples to indices
sample2index = {}
i = 0
for line in open(args.samples):
    sample = line.rstrip()
    sample2index[sample] = i
    i += 1
M = len(sample2index)

# Initialize numpy arrays for each genome
gene2beg = {}
x = {}
for line in open(args.gene_file):
    line = line.rstrip().split('\t')
    contig = line[0]
    gene = line[1]
    beg = int(line[2])
    end = int(line[3])
    if gene == args.gene:
        gene2beg[gene] = beg
        x = np.zeros([M, (end - beg + 1), 4])

# Add kpileup results to numpy arrays
for sample in sample2index:
    
    # read kpileup index
    [beg, end] = [-np.inf, np.inf]
    ifn = '%s/%s.kp.index' %(args.kp_dir, sample)
    if os.path.exists(ifn):
        for line in open(ifn):
            line = line.rstrip().split('\t')
            if line[0] == args.gene:
                beg = int(line[1])
                end = int(line[2])
                print(sample, beg, end)
                break
        else:
            print('skipping %s' %(sample))
            continue
    
    with open('%s/%s.kp.out' %(args.kp_dir, sample)) as fh:
        for i, line in enumerate(fh):
            
            if i <= 1:
                continue
            
            if beg <= i and i <= end:
                
                line = line.rstrip().split()
                sample = line[0]
                contig = line[1]
                pos = int(line[2])
                gene = line[3]
                nt = line[7]
                count = int(line[8])

                if args.gene != gene:
                    continue
                
                i = sample2index[sample]
                j = pos - gene2beg[gene]
                k = nts.index(nt)
                
                x[i,j,k] = count
                
            elif i >= end:
                break
            else:
                continue

# Filter alignment
I = np.arange(x.shape[0])
J = np.arange(gene2beg[gene], gene2beg[gene] + x.shape[1])

print(x.shape)

# Remove samples with less than 50% coverage
if x.shape[0] > 0:
    i = (x.sum(axis=2) > 0).mean(axis=1) >= 0.50
    x = x[i,:,:]
    I = I[i]

print(x.shape)

# Remove positions with less than 5% coverage
if x.shape[1] > 0:
    j = (x.sum(axis=2) > 0).mean(axis=0) >= 0.05
    x = x[:,j,:]
    J = J[j]

print(x.shape)

# Remove monomorphic positions
if x.shape[1] > 0:
    j = ((x >= 10).sum(axis=0) > 0).sum(axis=1) > 1
    x = x[:,j,:]
    J = J[j]
    
print(x.shape)

print(J)

# Write numpy arrays to file
res = {'x':x, 'i':I, 'j':J}
pickle.dump(res, open(args.out, 'wb'))
