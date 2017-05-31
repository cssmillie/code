# Dereplicate sequences in fasta file
# Option to dereplicate with a database to avoid OTU name collisions

import util, argparse, re, os

def parse_args():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='Input FASTA file', default='')
    parser.add_argument('-q', help='Input FASTQ file', default='')
    parser.add_argument('-d', help='Input OTU database', default='')
    parser.add_argument('-u', help='Output OTU database', required=True)
    parser.add_argument('-m', help='Output mapping file', required=True)
    parser.add_argument('-s', help='Sample ID separator', required=True)
    parser.add_argument('-M', help='Min size', default=10, type=int)
    parser.add_argument('-S', help='Min samples', default=3, type=int)
    parser.add_argument('-l', help='Trim length', type=int, default=0)
    args = parser.parse_args()
    return args


def load_db(fn, trim_len):
    # Load OTU database (otu id -> sequence)
    if not fn:
        return {}
    db = {}
    for [sid, seq] in util.iter_fst(fn):
        sid = int(sid)
        if trim_len:
            if len(seq) >= trim_len:
                seq = seq[:trim_len]
            else:
                continue
        db[sid] = seq
    return db


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


def write_output(x, db, ufn, mfn, min_size=1, min_samples=1):
    # Write output (dereplicated fasta file + mapping file)
    fst_out = open(ufn, 'a')
    map_out = open(mfn, 'w')
    if len(db) == 0:
        otu_id = 0
    else:
        otu_id = max(db.keys())
    for seq in x:
        if sum(x[seq].values()) < min_size:
            continue
        if len(x[seq]) < min_samples:
            continue
        if seq in db:
            seqid = db[seq]
        else:
            otu_id += 1
            seqid = otu_id
            fst_out.write('>%s\n%s\n' %(seqid, seq))
        map_out.write('%s\t%s\n' %(seqid, ' '.join(['%s:%d' %(sa, x[seq][sa]) for sa in x[seq]])))
    fst_out.close()
    map_out.close()


args = parse_args()
if args.l == 0:
    args.l = ''
db = load_db(args.d, args.l)
if args.d and args.u:
    os.system('cp %s %s' %(args.d, args.u))
x = dereplicate(fst = args.f, fsq = args.q, sep = args.s, trim_len = args.l)
write_output(x, db, args.u, args.m, min_size=args.M, min_samples=args.S)
