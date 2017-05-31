import pandas as pd
import argparse
import numpy as np

# Convert a derep mapping file to an OTU table

parser = argparse.ArgumentParser()
parser.add_argument('--map', help='mapping file')
parser.add_argument('--min_count', help = 'minimum read count', type = int)
parser.add_argument('--min_samples', help='minimum number of samples', type = int)
parser.add_argument('--out', help='output table')

args = parser.parse_args()

samples = {}
otus = {}

for line in open(args.map):
    otu, table = line.rstrip().split('\t')
    entries = table.split(' ')
    count = sum([int(entry.split(':')[1]) >= args.min_count for entry in entries])
    if count < args.min_samples:
        continue
    if otu not in otus:
        otus[otu] = len(otus)
    for entry in entries:
        sample, count = entry.split(':')
        if sample not in samples:
            samples[sample] = len(samples)

x = np.zeros([len(samples), len(otus)])

for line in open(args.map):
    otu, table = line.rstrip().split('\t')
    if otu in otus:
        for entry in table.split(' '):
            sample, count = entry.split(':')
            i = samples[sample]
            j = otus[otu]
            x[i,j] += int(count)

x = pd.DataFrame(x)

sort_otus = sorted(otus.keys(), key=lambda a: otus[a])
sort_samp = sorted(samples.keys(), key=lambda a: samples[a])

keep = ((x > 0).sum(axis=0) > 10)

x = x.ix[:, keep]
sort_otus = np.array(sort_otus)[keep]

out = open(args.out, 'w')
out.write( 'sample\t' + '\t'.join(sort_otus) + '\n')
for sample in sort_samp:
    i = samples[sample]
    out.write( '%s\t%s' %(sample, '\t'.join(['%d' %(xi) for xi in x.ix[i,:]])) + '\n')
