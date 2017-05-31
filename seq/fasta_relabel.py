import argparse
import re
import util

parser = argparse.ArgumentParser()
parser.add_argument('-f', help='fasta file')
parser.add_argument('-q', help='fastq file')
parser.add_argument('-p', help='prefix')
args = parser.parse_args()

if args.f:
    fn = args.f
    iter_seq = util.iter_fst
elif args.q:
    fn = args.q
    iter_seq = util.iter_fsq

for record in iter_seq(fn):
    old_name = record[0][1:].split(' ')[0]
    new_name = args.p + old_name
    record[0] = re.sub(old_name, new_name, record[0])
    print '\n'.join(record)
