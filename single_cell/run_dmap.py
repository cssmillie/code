import argparse
import numpy as np
import pandas as pd
import scanpy.api as sc
sc.settings.verbosity = 3


# -------------------------------------------------------------------------------
# Runs scanpy: diffusion map, diffusion pseudotime, approximate graph abstraction
# -------------------------------------------------------------------------------


# ---------------
# input arguments
# ---------------

parser = argparse.ArgumentParser()

# general
group1 = parser.add_argument_group('general')
group1.add_argument('--data', help='Input matrix [cells x feats]')
group1.add_argument('--out', help='Output prefix', default='test')

# dmap
group2 = parser.add_argument_group('dmap')
group2.add_argument('--n_pcs', help='Number of PCs (0 = no PCA)', type=int, default=0)

# dpt
group3 = parser.add_argument_group('dpt')
group3.add_argument('--dpt', help='Run DPT?', default=False, action='store_true')
group3.add_argument('--iroot', help='Index of root cell', type=int, default=0)
group3.add_argument('--n_dcs', help='Number of DCs', type=int, default=10)
group3.add_argument('--n_branches', help='Number of branchings', type=int, default=0)
group3.add_argument('--min_size', help='Min group size (fraction)', type=float, default=.001)

# aga
group4 = parser.add_argument_group('aga')
group4.add_argument('--aga', help='Run AGA?', default=False, action='store_true')
group4.add_argument('--clusters', help='Clusters [cells x 1]', default='')

args = parser.parse_args()


# ---------
# load data
# ---------

# matrix
x = pd.read_csv(args.data, sep='\t', header=None, index_col=False)
x = np.array(x)

# clusters
if args.clusters:
    clusters = pd.read_csv(args.clusters, sep='\t', header=None, index_col=False)
    clusters = np.array(clusters)
    g = 'groups'
else:
    clusters = ''
    g = 'louvain_groups'


# -------------
# diffusion map
# -------------

adata = sc.AnnData(X=x)
adata.uns['iroot'] = args.iroot

if type(clusters) != type(''):
    adata.obs['groups'] = clusters

try:
    sc.tl.diffmap(adata, n_pcs=args.n_pcs, n_comps=50)
except:
    pass


# --------------------
# diffusion pseudotime
# --------------------

if args.dpt:
    try:
        sc.tl.dpt(adata, n_pcs=args.n_pcs, n_dcs=args.n_dcs, n_branchings=args.n_branches, min_group_size=args.min_size)
    except:
        pass


# -----------------------------    
# approximate graph abstraction
# -----------------------------

if args.aga:
    try:
        sc.tl.aga(adata, groups=g, n_pcs=args.n_pcs, n_dcs=args.n_dcs)
    except:
        pass


# ------------------
# write output files
# ------------------

out_map = {'X_diffmap': adata.obsm, 'dpt_pseudotime': adata.obs, 'dpt_groups': adata.obs,
           'aga_adjacency_tree_confidence': adata.uns, 'aga_adjacency_full_attachedness': adata.uns, 'aga_adjacency_full_confidence': adata.uns, 'aga_pseudotime': adata.obs}

for obj in out_map:
    try:
        out = out_map[obj][obj]
        if obj == 'aga_adjacency_tree_confidence':
            out = out.toarray()
        np.savetxt('%s.%s.txt' %(args.out, obj), out, delimiter='\t', fmt='%s')
    except:
        pass

try:
    np.savetxt('%s.categories.txt' %(args.out), adata.obs['groups'].cat.categories, delimiter='\t', fmt='%s')
except:
    pass
