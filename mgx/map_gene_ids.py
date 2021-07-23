import glob, re

imap = {}

# map gene ids
for line in open('GCF_000482265.1.gtf'):
    line = line.rstrip().split('\t')
    try:
        old = re.search('gene_id "(.*?)"', line[8]).group(1)
        new = re.search('gene "(.*?)"', line[8]).group(1)
        imap[old] = new
    except:
        pass

# map counts
fns = glob.glob('*gene_count.txt')
for old_fn in fns:
    new_fn = re.sub('.txt', '.map.txt', old_fn)
    out = open(new_fn, 'w')
    for line in open(old_fn):
        line = line.rstrip().split('\t')
        if line[0] in imap:
            line[0] = imap[line[0]]
        out.write('\t'.join(line) + '\n')
    out.close()
