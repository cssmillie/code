import argparse, glob, re

# input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--regex', help='regex')
parser.add_argument('--level', help='taxonomic level (G,S)')
args = parser.parse_args()

# initialize data
data = {}
rows = []
cols = []

for fn in glob.glob(args.regex):
    sample = re.sub('.bracken', '', fn)
    data[sample] = {}
    rows.append(sample)
    for line in open(fn):
        if line.startswith('name'):
            continue
        line = line.rstrip().split('\t')
        name = line[0]
        count = int(line[5])
        data[sample][name] = count
        if name not in cols:
            cols.append(name)

rows = sorted(rows)
cols = sorted(cols)
print('Sample\t' + '\t'.join(cols))
for ri in rows:
    out = [ri]
    for ci in cols:
        out.append(str(data[ri].get(ci, 0)))
    print('\t'.join(out))

