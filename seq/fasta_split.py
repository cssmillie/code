'''
split fasta or fastq into k files

usage:
  python fasta_split.py in.fst 10

creates 3 files:
  in.1.fst
  in.2.fst
  in.3.fst
  
'''

import argparse, itertools, os.path, sys
from util import *

parser = argparse.ArgumentParser()
parser.add_argument('-f', help='FASTA file')
parser.add_argument('-q', help='FASTQ file')
parser.add_argument('-k', help='Number of files to split into', type=int)
parser.add_argument('-o', help='Output prefix')
args = parser.parse_args()

if args.f:
    iter_seq = iter_fst
    fn = args.f
if args.q:
    iter_seq = iter_fsq
    fn = args.q

# Open files
suffix = fn.split('.')[-1]
fns = ['%s.%d.%s' %(args.o, i, suffix) for i in range(args.k)]
for fn in fns:
    if os.path.exists(fn):
        exit('file %s exists' %(fn))
fhs = cycle([open(fn, 'w') for fn in fns])


# Write files
for record in iter_seq(fn):
    fh = fhs.next()
    fh.write('\n'.join(record) + '\n')

