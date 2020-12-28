import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='input sequence (fasta)', required=True)
parser.add_argument('-o', help='rRNA output (fasta)', required=True)
args = parser.parse_args()

if '.gz' in args.i:
    cmd = 'zcat %s | barrnap -o %s' %(args.i, args.o)
else:
    cmd = 'cat %s | barrnap -o %s' %(args.i, args.o)

print(cmd)
