import argparse

# get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--fst', help='input fasta file')
parser.add_argument('--out', help='output file prefix')
args = parser.parse_args()

# run kraken
cmd = '/home/unix/csmillie/bin/kraken2/kraken2 --db /home/unix/csmillie/aviv/db/kraken/fungi %s --report %s.kraken.report.txt > %s.kraken.out' %(args.fst, args.out, args.out)
print(cmd)
