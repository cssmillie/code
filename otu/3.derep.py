# Dereplicate sequences in fasta file

import util, argparse, re, os, seqdb

def parse_args():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='Input FASTA file', default='')
    parser.add_argument('-q', help='Input FASTQ file', default='')
    parser.add_argument('-s', help='Sample ID separator', required=True)
    parser.add_argument('-o', help='Output mapping file', required=True)
    parser.add_argument('-d', help='OTU database', default='', required=True)
    parser.add_argument('-M', help='Min count', default=10, type=int)
    parser.add_argument('-S', help='Min samples', default=3, type=int)
    parser.add_argument('-l', help='Trim length', type=int, default=0)
    args = parser.parse_args()
    return args


def dereplicate(fst='', fsq='', sep='_', trim_len=''):
    # Dereplicate sequences
    x = {}
    if fst:
        fn = fst
        iter_fst = util.iter_fst
    if fsq:
        fn = fsq
        iter_fst = util.iter_fsq
    for record in iter_fst(fn):
        [sid, seq] = record[:2]
        sa = re.search('(.*?)%s' %(sep), sid).group(1)
        if trim_len:
            if len(seq) >= trim_len:
                seq = seq[:trim_len]
            else:
                continue
        if seq not in x:
            x[seq] = {}
        if sa not in x[seq]:
            x[seq][sa] = 0
        x[seq][sa] += 1
    return x


def write_output(x, map_fn, db_fn, min_size=1, min_samples=1):
    # Write output (database + mapping file)
    # Load SeqDB
    db = seqdb.SeqDB(fn=db_fn)
    out = open(map_fn, 'w')
    for seq in x:
        if sum(x[seq].values()) < min_size:
            continue
        if len(x[seq]) < min_samples:
            continue
        sid = db.get_otu(seq)
        out.write('%s\t%s\n' %(sid, ' '.join(['%s:%d' %(sa, x[seq][sa]) for sa in x[seq]])))
    out.close()
    db.write(db_fn)


args = parse_args()
if args.l == 0:
    args.l = ''
x = dereplicate(fst=args.f, fsq=args.q, sep=args.s, trim_len=args.l)
write_output(x, map_fn=args.o, db_fn=args.d, min_size=args.M, min_samples=args.S)
