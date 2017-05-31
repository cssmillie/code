import argparse

# -----------------------------------------------------------------------
# inputs:
# 1. tree file (newick)
# 2. metadata file (map of labels -> continuous trait)
#
# this script assumes that the metadata ranges from negative to positive
# it labels every internal node and calculates the mean at each node
#
# outputs:
# 1. tree file (newick) - with node labels
# 2. itol file - metadata for all nodes. +/- labeled separately
# -----------------------------------------------------------------------



# read command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--in_tree', help='Input tree (newick)', required=True, default='')
parser.add_argument('--in_meta', help='Metadata (continuous)', required=True, default='')
parser.add_argument('--out_tree', help='Output tree (newick)', required=True, default='')
parser.add_argument('--out_meta', help='Output iTol file', required=True, default='')
parser.add_argument('--out_popup', help='Output popup file', required=True, default='')
args = parser.parse_args()

import ete2, re
import numpy as np

# load tree
t = ete2.Tree(args.in_tree, format=1)

# get map
b2r = {}
for line in open(args.in_meta):
    try:
        line = line.rstrip().split()
        b2r[line[0]] = float(line[1])
    except:
        # skip header
        continue


# write itol file
out = open(args.out_meta, 'w')
popup = open(args.out_popup, 'w')

# for every node
i = 0 # node label
for node in t.traverse():
    i += 1
    if node.is_leaf():
        try:
            res = b2r[node.name]
        except:
            continue
        out.write('%s\t%f\n' %(node.name, res))
    else:
        if node.name == '':
            node.name = 'n%d' %(i)
        if node.name:
            popup.write('%s\tlabel\t%s\n' %(node.name, re.sub('[^0-9a-zA-Z]+', '_', node.name)))
        leaves = node.get_leaves()
        res = [b2r[leaf.name] for leaf in leaves if leaf.name in b2r]
        if len(res) == 0:
            continue
        else:
            res = np.mean(res)
        out.write('%s\t%f\n' %(node.name, res))
out.close()

# write relabeled tree
t.write(outfile=args.out_tree, format=1)
