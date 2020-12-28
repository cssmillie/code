import argparse, sys
import scanpy as sc
import scipy.io
import numpy as np

# read input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', help='input file', default='NoCiteSeq_demux.annotated.h5ad')
parser.add_argument('-o', help='output prefix', default='NoCiteSeq')
args = parser.parse_args()

# read h5ad file
adata = sc.read(args.i)

# find counts key
ks = []
for ki in adata.layers.keys():
    if ki.contains('count'):
        ks.append(ki)
print(ks)
if len(ks) == 0:
    try:
        scipy.io.mmwrite('%s.matrix.mtx' %(args.o), adata.raw.X)
    except:
        scipy.io.mmwrite('%s.matrix.mtx' %(args.o), adata.X)
elif len(ks) == 1:
    ki = ks[0]
    scipy.io.mmwrite('%s.matrix.mtx' %(args.o), adata.layers[ki])
else:
    sys.exit('multiple keys contain "count"')

# write data
adata.obs.to_csv('%s.meta.csv' %(args.o))
adata.var_names.to_frame().to_csv('%s.genes.csv' %(args.o))
adata.obs_names.to_frame().to_csv('%s.barcodes.csv' %(args.o))

if 'X_umap' in adata.obsm:
    np.savetxt('%s.umap.csv' %(args.o), adata.obsm['X_umap'], delimiter=',')

if 'X_pca' in adata.obsm:
    np.savetxt('%s.pca.csv' %(args.o), adata.obsm['X_pca'], delimiter=',')
    
