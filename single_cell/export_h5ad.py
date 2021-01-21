import argparse, sys
import scanpy as sc
import scipy.io
import numpy as np
from PIL import Image

# read input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', help='input file', default='NoCiteSeq_demux.annotated.h5ad')
parser.add_argument('-o', help='output prefix', default='NoCiteSeq')
parser.add_argument('--dense', help='dense format', default=False, action='store_true')
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
        if args.dense:
            np.savetxt('%s.matrix.csv' %(args.o), adata.raw.X, delimiter=',')
        else:
            scipy.io.mmwrite('%s.matrix.mtx' %(args.o), adata.raw.X)
    except:
        if args.dense:
            np.savetxt('%s.matrix.csv' %(args.o), adata.X, delimiter=',')
        else:
            scipy.io.mmwrite('%s.matrix.mtx' %(args.o), adata.X)
elif len(ks) == 1:
    ki = ks[0]
    if args.dense:
        np.savetxt('%s.matrix.csv' %(args.o), adata.layers[ki], delimiter=',')
    else:
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
    
if 'image_hires' in adata.uns:
    im = Image.fromarray(adata.uns['image_hires'])
    im.save('%s.image_hires.png' %(args.o))
