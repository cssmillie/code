import argparse, sys, util

# Remove primers from sequences

def parse_args():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', default='', help='Input FASTA file')
    parser.add_argument('-q', default='', help='Input FASTQ file')
    parser.add_argument('-p', default='', help='Primer sequence')
    parser.add_argument('-l', default='', help='Primer list')
    parser.add_argument('-d', default=1, type=int, help='Max primer differences')
    parser.add_argument('-w', default=10, type=int, help='Search positions 1-w for primer')
    parser.add_argument('-b', default=False, help='Remove primers from both ends', action='store_true')
    args = parser.parse_args()
    
    # Check for consistency
    if args.f and args.q:
        quit('Error: cannot specify both FASTA and FASTQ file')
    if args.p and args.l:
        quit('Error: cannot specify both primer and primer list')
    
    return args


def mismatches(seq, subseq, w):
    # Calculate the number of mismatches between a sequence and a given subsequence
    # Searches in a sliding window that starts at position 1 and ends at position w
    best_i = 0 # index (start position)
    best_d = len(seq) # edit distance
    # for every start position
    for i in range(w):
        # calculate edit distance to the given subsequence
        d = primer.MatchPrefix(seq[i:], subseq)
        # keep track of the best index and edit distance
        if d < best_d:
            best_i = i
            best_d = d
    return [best_i, best_d]


def find_best_match(seq, primers, w, max_diff):
    # Find the sample with the best matching barcode
    best_i = ''
    best_p = '' # barcode
    best_d = len(seq) # edit distance
    # Calculate edit distance to every barcode
    for p in primers:
        [i,d] = mismatches(seq, p, w) # index, edit distance
        if d < best_d:
            best_i = i
            best_p = p
            best_d = d
    # Return [index, edit distance, barcode, sample id] of best match
    if best_d <= max_diff:
        return [best_i, best_d, best_p]
    else:
        return ['', '', '']


def run():
    # Remove primers from FASTA/FASTQ file
    
    # Get command line arguments
    args = parse_args()
    
    # Get primer sequences
    if args.p:
        primers = [args.p]
    elif args.l:
        primers = [line.rstrip() for line in open(args.l)]
    else:
        quit('Error: must specify primer or primer list')

    # Get FASTA/FASTQ iterators
    if args.f:
        fn = args.f
        iter_fst = util.iter_fst
    elif args.q:
        fn = args.q
        iter_fst = util.iter_fsq
    else:
        quit('Error: must specify FASTA or FASTQ file')
    
    # Iterate through FASTA/FASTQ file
    n = 0
    t = 0
    for record in iter_fst(fn):
        t += 1
        seq = record[1]
        [i,d,p] = find_best_match(seq, primers, args.w, args.d)
        if i != '':
            record[1] = record[1][i + len(p):]
            record[3] = record[3][i + len(p):]
            if args.b == False:
                print '\n'.join(record)
                n += 1
            if args.b == True:
                SEQ = util.reverse_complement(record[1])
                [I,D,P] = find_best_match(SEQ, primers, args.w, args.d)
                if I != '':
                    record[1] = record[1][:-(I + len(P))]
                    record[3] = record[3][:-(I + len(P))]
                    print '\n'.join(record)
                    n += 1
    
    sys.stderr.write('%d total sequences\n' %(t))
    sys.stderr.write('%d primers removed (%.2f%%)\n' %(n, 100.*n/t))
                    
run()
