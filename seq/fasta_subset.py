import argparse
import util

# parse args
parser = argparse.ArgumentParser()
parser.add_argument('-f', help='FASTA file')
parser.add_argument('-q', help='FASTQ file')
parser.add_argument('-s', help='Subset ids')
args = parser.parse_args()

# load subset
subset = [line.rstrip() for line in open(args.s)]

# get iterator
iter_seq = ''
if args.f:
    iter_seq = util.iter_fst(args.f)
if args.q:
    iter_seq = util.iter_fsq(args.q)

# subset file
for record in iter_seq:
    sid = record[0][1:].split(';')[0]
    if sid in subset:
        print '\n'.join(record)
    
