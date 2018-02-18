import argparse
import magic
import sys

# Run magic imputation
# input = TPM or log2TPM (authors use TPM)
# output = imputed data

# Arguments
# ---------

# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--txt', help='TXT file (rows = cells, columns = features)', default='')
parser.add_argument('--mtx', help='MTX file (rows = cells, columns = features)', default='')
parser.add_argument('--genes', help='Genes file', default='')
parser.add_argument('-p', help='Number of PCs', default=20, type=int)
parser.add_argument('-k', help='Number of nearest neighbors', default=30, type=int)
parser.add_argument('--ka', help='Knn autotune parameter', default=10, type=int)
parser.add_argument('-e', help='Epsilon', default=1, type=float)
parser.add_argument('-r', help='Rescale percent', default=99, type=float)
parser.add_argument('--out', help='Outfile (clusters)')
args = parser.parse_args()

# Check arguments
if args.txt == '' and args.mtx == '':
    sys.exit('must specify txt or mtx')
if args.mtx != '' and args.genes == '':
    sys.exit('must specify mtx and genes')


# Run magic
# ---------

def run_magic(args):
    
    if args.txt != '':
        print('Reading TXT')
        scdata = magic.mg.SCData.from_csv(args.txt, data_type='sc-seq', normalize=False)
    
    if args.mtx != '':
        print('Reading MTX')
        scdata = magic.mg.SCData.from_mtx(args.mtx, args.genes, normalize=False)
    
    print('Running magic')
    scdata.run_magic(n_pca_components=args.p, random_pca=True, t=None, k=args.k, ka=args.ka, epsilon=args.e, rescale_percent=args.r)
    
    print('Writing data')
    scdata.magic.to_csv(args.out)


run_magic(args)
