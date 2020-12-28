import argparse, math, util

# read input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--fst', help='FASTA file', required=True)
parser.add_argument('--map', help='Sequence ID map', default='')
parser.add_argument('--bin', help='Divide genome into bins', default=1, type=int)
args = parser.parse_args()

# get sequence alphabet
alphabet = sorted('G A L M F W K Q E S P V I C Y H R N D T'.split())

# divide sequence into bins
if args.bin == 1:
    bins = []
else:
    bins = range(1, args.bin + 1)


# construct matrix header
if len(bins) == 0:
    header = ['gene'] + ['count_%s' %(aa) for aa in alphabet] + ['freq_%s' %(aa) for aa in alphabet]
else:
    header = ['gene'] + ['count_full_%s' %(aa) for aa in alphabet] + ['freq_full_%s' %(aa) for aa in alphabet]

for bin in bins:
    header = header + ['count_q%d_%s' %(bin, aa) for aa in alphabet] + ['freq_q%d_%s' %(bin, aa) for aa in alphabet]

print('\t'.join(header))
        
# read sequence id map
smap = {}
if args.map:
    for line in open(args.map):
        line = line.rstrip().split()
        smap[line[0]] = line[1]

# iterate over sequences
for record in util.iter_fst(args.fst):
    sid, seq = record

    # map sequence ids
    sid = sid.split()[0][1:]
    if args.map:
        if sid in smap:
            sid = smap[sid]
        else:
            continue

    # count amino acids
    count = {}
    for i in range(len(seq)):
        aa = seq[i]
        count[aa] = count.get(aa, 0) + 1
        if args.bin > 1:
            aq = '%s_Q%d' %(aa, i / math.ceil(len(seq) / args.bin) + 1)
            count[aq] = count.get(aq, 0) + 1
    
    # calculate frequencies
    freq = {}

    total = sum([count.get(aa, 0) for aa in alphabet])
    for aa in alphabet:
        freq[aa] = count.get(aa, 0)/total
        
    for bin in bins:
        total = sum([count.get('%s_Q%d' %(aa, bin), 0) for aa in alphabet])
        for aa in alphabet:
            aq = '%s_Q%d' %(aa, bin)
            freq[aq] = count.get(aq, 0)/total
    
    # construct output
    out = [sid] + [count.get(aa, 0) for aa in alphabet] + [freq.get(aa, 0) for aa in alphabet]
    for bin in bins:
        out = out + [count.get('%s_Q%d' %(aa, bin), 0) for aa in alphabet] + [freq.get('%s_Q%d' %(aa, bin), 0) for aa in alphabet]
    print('\t'.join(map(str, out)))
