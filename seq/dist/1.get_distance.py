import sys
import pandas as pd
import scipy.stats
import nwalign as nw
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fst', default = '', help='Input fasta file')
parser.add_argument('--aln', default = '', help='Input aln file')
parser.add_argument('--subset', default = '', help='Subset file (list)')
parser.add_argument('--beg1', required=True, type=int, help='beg1')
parser.add_argument('--end1', required=True, type=int, help='end1')
parser.add_argument('--beg2', required=True, type=int, help='beg2')
parser.add_argument('--end2', required=True, type=int, help='end2')
args = parser.parse_args()

fst = {}

fst_fn = args.fst
aln_fn = args.aln
beg1 = args.beg1
end1 = args.end1
beg2 = args.beg2
end2 = args.end2

if aln_fn != '':
    fst_fn = aln_fn

subset = []
if args.subset != '':
    subset = [line.rstrip() for line in open(args.subset).readlines()]

for line in open(fst_fn):
    if line.startswith('>'):
        seq_id = line.rstrip().split()[0][1:]
        if len(subset) == 0:
            fst[seq_id] = ''
        elif seq_id in subset:
            fst[seq_id] = ''
    else:
        try:
            fst[seq_id] += line.rstrip()
        except:
            continue

data = []

def calc_distance(x,y):
    if aln_fn == '':
        x, y = nw.global_align(x, y, gap_open=-5, gap_extend=-1, match=1)
    dxy = 0
    lxy = 0
    for i in range(min(len(x),len(y))):
        if x[i] != '-' and y[i] != '-':
            lxy += 1
            if x[i] != y[i]:
                dxy += 1
    return 1.*dxy/lxy

otus = sorted(fst.keys())

for i in range(beg1, end1):
    if (i < 0) or (i > len(otus)):
        continue
    oi = otus[i]
    xi = fst.get(oi, '')
    for j in range(beg2, end2):
        if (j < 0) or (j >= i):
            continue
        oj = otus[j]
        xj = fst.get(oj, '')
        if xi != '' and xj != '':
            dij = calc_distance(xi,xj)
            data.append('%.4f' %(dij))
        else:
            data.append('NA')

print '\t'.join(data)
