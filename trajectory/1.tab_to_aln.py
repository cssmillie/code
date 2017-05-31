import argparse, sys

# Get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--tab', help='alignment matrix (generated in R)', required=True)
parser.add_argument('--sep', help='sample separator', default=';')
parser.add_argument('--phy', help='phylip output file', default='')
parser.add_argument('--fst', help='fasta output file (converted to dna)', default='')
parser.add_argument('--rmv', help='samples to remove from alignment', default='')
args = parser.parse_args()

# Initialize parameters
remove = args.rmv.split(',')

# Write PHYLIP file
if args.phy:    
    m = 0
    n = 0
    l = 9
    for line in open(args.tab):
        line = line.rstrip().split(args.sep)
        if line[0] in remove:
            continue
        m += 1
        if n == 0:
            n = len(line[1])
        else:
            if n != len(line[1]):
                quit('alignment lengths not equal')
        l = max(l, len(line[0]))

    out = open(args.phy, 'w')
    out.write('%s %s\n' %(m, n))    
    for line in open(args.tab):
        line = line.rstrip().split(args.sep)
        if line[0] in remove:
            continue
        out.write(line[0] + ' '*(l + 1 - len(line[0])) + line[1] + '\n')
    out.close()

# Write FASTA file
if args.fst:
    
    from string import maketrans
    nums = '0123'
    dnas = 'ACGT'
    ttab = maketrans(nums, dnas)
    
    out = open(args.fst, 'w')
    for line in open(args.tab):
        line = line.rstrip().split(args.sep)
        if line[0] in remove:
            continue
        out.write('>%s\n%s\n' %(line[0], line[1].translate(ttab)))
    out.close()
