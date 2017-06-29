import argparse, sys
from Bio import SeqIO

# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--abi', help='ABI file', required=True)
parser.add_argument('--fmt', help='Format (fasta, fastq, fastq-illumina)', default='fastq-illumina')
parser.add_argument('--out', help='Output file', required=True)
args = parser.parse_args()

# Convert abi to output format
out = open(args.out, 'w')
for record in SeqIO.parse(open(args.abi, 'rb'), 'abi'):
    out.write('%s\n' %(record.format(args.fmt)))
out.close()
