import argparse
import glob
import re

# input arguments
# ---------------
parser = argparse.ArgumentParser()
parser.add_argument('--infile', help='KEGG pathways (file)')
parser.add_argument('--inlist', help='KEGG pathways (list)')
args = parser.parse_args()

# get filenames
# -------------
if args.infile:
    fns = [args.infile]
else:
    fns = [line.rstrip() for line in open(args.inlist)]

# parse kegg files
# ----------------
c2num = {}
c2lab = {}
cflag = 0
for fn in fns:
    for line in open(fn):
        if line.startswith('COMPOUND'):
            line = re.sub('COMPOUND', '', line)
            cflag = 1
        if line.startswith('REFERENCE'):
            cflag = 0
        m = re.search(' *(C[0-9]*)(.*)', line)
        if cflag and m:
            cid = m.group(1)
            lab = m.group(2)
            if cid != 'C':
                c2num[cid] = c2num.get(cid, 0) + 1
                c2lab[cid] = lab

for cid in c2num:
    print('%s\t%s\t%s' %(cid, c2num[cid], c2lab[cid]))
