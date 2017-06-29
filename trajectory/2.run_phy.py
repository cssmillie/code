import argparse, random

parser = argparse.ArgumentParser()
parser.add_argument('--aln', help='phylip or fasta alignment')
parser.add_argument('--boot', help='bootstrap', default=False, action='store_true')
parser.add_argument('--fast', help='fasttree', default=False, action='store_true')
args = parser.parse_args()

r1 = random.randint(0, 32767)
r2 = random.randint(0, 32767)

prefix = '.'.join(args.aln.split('.')[:-1])

if args.fast == False:
    
    if args.aln.split('.')[-1] != 'phy':
        quit('Check alignment format')

    if args.boot:
        cmd = '~/bin/standard-RAxML/raxmlHPC-SSE3 -s %s -n %s -m MULTIGAMMA -p %d -f a -N 100 -x %d' %(args.aln, prefix, r1, r2)
    else:
        cmd = '~/bin/standard-RAxML/raxmlHPC-SSE3 -s %s -n %s -m MULTIGAMMA -p %d -f D' %(args.aln, prefix, r1)
else:
    
    if args.aln.split('.')[-1] == 'phy':
        quit('Check alignment format')

    cmd = '~/bin/FastTree -nt -gtr -gamma < %s > %s.tree' %(args.aln, prefix)

print cmd
