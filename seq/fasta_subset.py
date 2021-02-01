#!/usr/bin/env python

import argparse, re, sys, util

# parse args
parser = argparse.ArgumentParser()
parser.add_argument('--fst', help='FASTA file', default='')
parser.add_argument('--fsq', help='FASTQ file', default='')
parser.add_argument('--FST', help='FASTA flag (for sys.stdin)', action='store_true', default=False)
parser.add_argument('--FSQ', help='FASTQ flag (for sys.stdin)', action='store_true', default=False)
parser.add_argument('--regex', help='Sequence ID regex', default='')
parser.add_argument('--keep', help='IDs to keep (.txt)')
parser.add_argument('--remove', help='IDs to remove (.txt)')
parser.add_argument('--minlen', help='Minimum length', type=int, default=0)
parser.add_argument('--num', help='Number to keep', type=int, default=0)
parser.add_argument('--prefix', help='Prefix to add', type=str, default='')
parser.add_argument('--prefix_sep', help='Prefix separator', type=str, default='.')
parser.add_argument('--debug', help='Debug mode', action='store_true', default=False)
args = parser.parse_args()

# get iterator
if args.fst:
    iter_seq = util.iter_fst(args.fst)
elif args.fsq:
    iter_seq = util.iter_fsq(args.fsq)
elif args.FST:
    iter_seq = util.iter_fst(sys.stdin)
elif args.FSQ:
    iter_seq = util.iter_fsq(sys.stdin)
else:
    quit('error: must specify fst, fsq, FST, or FSQ')

# initialize variables
keep = {}
remove = {}

# load IDs/coordinates to keep
if args.keep:
    for line in open(args.keep):
        line = line.rstrip().split('\t')
        sid = line[0]
        if len(line) > 1:
            beg = min(map(int, line[1:]))
            end = max(map(int, line[1:]))
        else:
            beg = None
            end = None
        if sid not in keep:
            keep[sid] = []
        keep[sid].append([beg, end])

# load IDs/coordinates to remove
if args.remove:
    for line in open(args.remove):
        line = line.rstrip().split('\t')
        sid = line[0]
        if len(line) > 1:
            beg = min(map(int, line[1:]))
            end = max(map(int, line[1:]))
        else:
            beg = None
            end = None
        if sid not in remove:
            remove[sid] = []
        remove[sid].append([beg, end])

# subset FASTA/FASTQ sequences

count = {'total': 0, 'keep': 0, 'remove': 0, 'trim': 0, 'pass': 0, 'fail': 0}

for record in iter_seq:
    
    if args.num > 0:
        if count['total'] >= args.num:
            break
    
    # extract record info
    sid = record[0][1:]
    seq = record[1]
    if args.fsq or args.FSQ:
        qua = record[3]
    
    # filter by regex
    if args.regex:
        m = re.search(args.regex, sid)
        if not m:
            continue
    
    # select regions to keep
    if args.keep:
        ind = [0]*len(seq)
        if sid.split()[0] in keep:
            count['keep'] += 1
            for element in keep[sid.split()[0]]:
                [beg, end] = element
                if beg is None:
                    beg = 1
                if end is None:
                    end = len(seq)
                ind[(beg-1):end] = [1]*(end-beg+1)
        else:
            continue
    
    # select regions to remove
    if args.remove:
        ind = [1]*len(seq)
        if sid in remove:
            count['remove'] += 1
            for element in remove[sid]:
                [beg, end] = element
                if beg is None:
                    beg = 1
                if end is None:
                    end = len(seq)
                ind[(beg-1):end] = [0]*(end-beg+1)

    if not args.keep and not args.remove:
        ind = [1]*len(seq)
    
    # trim sequence 
    seq = ''.join([seq[i] if ind[i] == 1 else 'N' for i in range(len(seq))])
    if args.fsq or args.FSQ:
        qua = ''.join([qua[i] if ind[i] == 1 else '#' for i in range(len(qua))])
    for i in range(len(seq)):
        if seq[i] != 'N':
            break
    for j in range(len(seq))[::-1]:
        if seq[j] != 'N':
            break
    seq = seq[i:(j+1)]
    if args.fsq or args.FSQ:
        qua = qua[i:(j+1)]
    
    if i > 0 or j < (len(seq) - 1):
        count['trim'] += 1
    
    # print sequence
    if len(seq) - seq.count('N') >= args.minlen:
        count['pass'] += 1
        record[1] = seq
        if args.fsq or args.FSQ:
            record[3] = qua
        record[0] = record[0][0] + args.prefix + args.prefix_sep + record[0][1:]
        print('\n'.join(record))
        count['total'] += 1
    else:
        count['fail'] += 1

if args.debug:
    for k in count:
        print(k, count[k])
