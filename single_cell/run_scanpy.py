import argparse
import copy
import numpy as np
import pandas as pd
import re
import scanpy
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
group1.add_argument('--counts', help='Input matrix [cells x feats] (sparse = .mtx)', default='pcd-transpose.matrix.mtx')
group1.add_argument('--data', help='Input matrix [cells x feats] (sparse = .mtx)', default='')
group1.add_argument('--transpose', help='Transpose matrix', default=False, action='store_true')
group1.add_argument('--n_pcs', help='Number of PCs (0=none) (leiden, dmap)', type=int, default=25)
group1.add_argument('--n_neighbors', help='Number of neighbors (leiden, dmap)', type=int, default=10)
group1.add_argument('--out', help='Output prefix', default='test')

# cluster
group2 = parser.add_argument_group('cluster')
group2.add_argument('--leiden', help='Run Leiden? (n_neighbors, n_pcs, resolution)', default=False, action='store_true')
group2.add_argument('--resolution', help='Clustering resolution', default=1.0, type=float)

# dmap
group3 = parser.add_argument_group('dmap')
group3.add_argument('--dmap', help='Run DMAP? (n_neighbors, n_pcs)', default=False, action='store_true')

# dpt
group4 = parser.add_argument_group('dpt')
group4.add_argument('--dpt', help='Run DPT?', default=False, action='store_true')
group4.add_argument('--iroot', help='Index of root cell', type=int, default=0)

# paga
group5 = parser.add_argument_group('aga')
group5.add_argument('--paga', help='Run PAGA?', default=False, action='store_true')
group5.add_argument('--clusters', help='Clusters [cells x 1]', default='')

# dca
group6 = parser.add_argument_group('dca')
group6.add_argument('--dca', help='Run DCA?', default=False, action='store_true')

# output
group7 = parser.add_argument_group('out')
group7.add_argument('--h5ad', help='Write h5ad filename', default='')
group7.add_argument('--npz', help='Write npz filename', default='')

args = parser.parse_args()

print(args)

# ---------
# load data
# ---------

adata = None

# counts
if args.counts:
    if '.mtx' in args.counts:
        print('Reading sparse matrix %s' %(args.counts), flush=True)
        x = sc.read_mtx(args.counts).X
    else:
        print('Reading dense matrix %s' %(args.data), flush=True)
        x = pd.read_csv(args.counts, sep='\t', header=None, index_col=False)
        x = np.array(x)
    if args.transpose:
        x = x.T
    adata = sc.AnnData(X=x, layers={'counts': copy.deepcopy(x)})
    print('Normalization (total = 10,000 UMIs per cell)')
    sc.pp.normalize_total(adata, target_sum=10000, inplace=True)
    sc.pp.log1p(adata, base=2)

# data
if args.data:
    if '.mtx' in args.data:
        print('Reading sparse matrix %s' %(args.data), flush=True)
        x = sc.read_mtx(args.data).X
    else:
        print('Reading dense matrix %s' %(args.data), flush=True)
        x = pd.read_csv(args.data, sep='\t', header=None, index_col=False)
        x = np.array(x)
    if args.transpose:
        x = x.T
    adata = sc.AnnData(X=x)

# rows/columns
barcodes = features = []
fn = ''
if '.mtx' in args.counts:
    fn = args.counts
if '.mtx' in args.data:
    fn = args.data
if fn:
    print('Adding features and barcodes')
    temp = [line.rstrip() for line in open(re.sub('matrix.mtx', 'barcodes.tsv', fn))]
    if len(temp) == x.shape[0]:
        barcodes = temp
        print('Read %d barcodes' %(len(temp)), flush=True)
    else:
        features = temp
        print('Read %d features' %(len(temp)), flush=True)
    temp = [line.rstrip() for line in open(re.sub('matrix.mtx', 'genes.tsv', fn))]
    if len(temp) == x.shape[0]:
        barcodes = temp
        print('Read %d barcodes' %(len(temp)), flush=True)
    else:
        features = temp
        print('Read %d features' %(len(temp)), flush=True)

print('Data = [%d cells x %d genes]' %(np.shape(x)), flush=True)

# clusters
if args.clusters:
    clusters = pd.read_csv(args.clusters, sep='\t', header=None, index_col=False)
    clusters = np.array(clusters)
    g = 'groups'
else:
    clusters = ''
    g = 'louvain_groups'

# outputs
out_map = {}


# -------------
# setup anndata
# -------------

adata.uns['iroot'] = args.iroot

if type(clusters) != type(''):
    adata.obs['groups'] = clusters

if len(barcodes) > 0:
    print('Adding barcodes: %s, ...' %(', '.join(barcodes[:3])), flush=True)
    adata.obs_names = barcodes

if len(features) > 0:
    print('Adding features: %s, ...' %(', '.join(features[:3])), flush=True)
    adata.var_names = features


# ----------
# clustering
# ----------

if args.leiden:
    try:
        print('Nearest neighbor graph', flush=True)
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)
        print('Clustering with leiden', flush=True)
        sc.tl.leiden(adata, resolution=args.resolution)
        out_map['leiden'] = adata.obs
    except:
        pass


# -------------
# diffusion map
# -------------

if args.dmap:
    try:
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)
        sc.tl.diffmap(adata)
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep='X_diffmap')
        out_map['X_diffmap'] = adata.obsm
        #sc.tl.diffmap(adata, n_pcs=args.n_pcs, n_comps=50, n_neighbors=args.n_neighbors)
    except:
        pass


# --------------------
# diffusion pseudotime
# --------------------

if args.dpt:
    try:
        sc.tl.dpt(adata)
        out_map['dpt_pseudotime'] = adata.obs
    except:
        pass


# -----------------------------    
# approximate graph abstraction
# -----------------------------

if args.paga:
    try:
        sc.tl.paga(adata, groups=g)
        out_map['connectivities'] = adata.uns['paga']
        out_map['connectivities_tree'] = adata.uns['paga']
    except:
        pass


# ----------------------
# deep count autoencoder
# ----------------------

if args.dca:
    try:
        sc.pp.dca(adata)
        np.savetxt('%s.%s.txt' %(args.out, 'dca'), adata.X, delimiter='\t', fmt='%s')
    except:
        pass


# ------------------
# write output files
# ------------------

print('Writing out', flush=True)
for obj in out_map:
    try:
        out = out_map[obj][obj]
        if obj in ['connectivities', 'connectivities_tree']:
            out = out.toarray()
        np.savetxt('%s.%s.txt' %(args.out, obj), out, delimiter='\t', fmt='%s')
    except:
        pass

try:
    np.savetxt('%s.categories.txt' %(args.out), adata.obs['groups'].cat.categories, delimiter='\t', fmt='%s')
except:
    pass

if args.h5ad:
    print('Writing h5ad', flush=True)
    adata.write_h5ad(args.h5ad)
