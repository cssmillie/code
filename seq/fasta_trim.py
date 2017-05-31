import util, sys

fn = sys.argv[1]
k = int(sys.argv[2])


for [sid, seq] in util.iter_fst(fn):

    if len(seq) >= k:
        print '>%s\n%s' %(sid, seq[:k])
