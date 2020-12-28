import argparse, translate, util

# read input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--fst', help='FASTA file', required=True)
parser.add_argument('--map', help='Sequence ID map', default='')
args = parser.parse_args()

# get sequence alphabet
nts = sorted('A C G T'.split())
cds = sorted(translate.codon_table.keys())
aas = sorted('G A L M F W K Q E S P V I C Y H R N D T'.split())
alphabet = nts + cds + aas
header = ['gene'] + ['count_nt_%s' %(li) for li in nts] + ['count_cd_%s' %(li) for li in cds] + ['count_aa_%s' %(li) for li in aas] + \
    ['freq_nt_%s' %(li) for li in nts] + ['freq_cd_%s' %(li) for li in cds] + ['freq_aa_%s' %(li) for li in aas]
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
    
    if len(seq) % 3 != 0:
        continue
    
    # map sequence ids
    sid = sid.split()[0][1:]
    if args.map:
        if sid in smap:
            sid = smap[sid]
        else:
            continue

    # count fna features
    # ------------------
    count = {}
    freq = {}
    
    for nt in seq:
        count[nt] = count.get(nt, 0) + 1
    
    for cd in translate.get_codons(seq):
        count[cd] = count.get(cd, 0) + 1

    for aa in translate.translate(seq):
        count[aa] = count.get(aa, 0) + 1

    # calculate frequencies
    # ---------------------
    
    tot_nt = sum([count.get(nt, 0) for nt in nts])
    for nt in nts:
        freq[nt] = count.get(nt, 0)/tot_nt
    
    tot_cd = sum([count.get(cd, 0) for cd in cds])
    for cd in cds:
        freq[cd] = count.get(cd, 0)/tot_cd
    
    tot_aa = sum([count.get(aa, 0) for aa in aas])
    for aa in aas:
        freq[aa] = count.get(aa, 0)/tot_aa

    out = map(str, [sid] + [count.get(li, 0) for li in alphabet] + [freq.get(li, 0) for li in alphabet])
    print('\t'.join(out))
    quit()
        

'''
nucleotides = sort('A C G T'.split())
amino_acids = sorted('G A L M F W K Q E S P V I C Y H R N D T'.split())



print('gene\t' + '\t'.join(alphabet))

# read sequence id map
smap = {}
if args.map:
    for line in open(args.map):
        line = line.rstrip().split()
        smap[line[0]] = line[1]

# iterate over sequences
for record in util.iter_fst(args.fst):
    sid, seq = record
    sid = sid.split()[0][1:]
    if args.map:
        if sid in smap:
            sid = smap[sid]
        else:
            continue
    count = {}
    for i in range(len(seq)):
        letter = seq[i]
        if letter not in count:
            count[letter] = 0
        count[letter] += 1
    out = map(str, [sid] + [count[li] if li in count else 0 for li in alphabet])
    print('\t'.join(out))
'''
