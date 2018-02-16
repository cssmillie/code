import argparse
import magic

# This script runs Phenograph 

# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--txt', help='TXT file (rows = cells, columns = features)', default='')
parser.add_argument('--mtx', help='MTX file (rows = cells, columns = features)', default='')
parser.add_argument('--genes', help='Genes file', default='')
parser.add_argument('-n', help='Number of PCs', default=20, type=int)
parser.add_argument('-k', help='Number of nearest neighbors', default=30, type=int)
parser.add_argument('--out', help='Outfile (clusters)')
args = parser.parse_args()

# Run phenograph
if args.txt != '':
    print('Reading TXT')
    scdata = magic.mg.SCData.from_csv(args.txt, data_type='sc-seq', normalize=False)

if args.mtx != '':
    print('Reading MTX')
    scdata = magic.mg.SCData.from_mtx(args.mtx, args.genes, normalize=False)

print('Running magic')
scdata.run_magic(n_pca_components=args.n, random_pca=True, t=6, k=args.k, ka=10, epsilon=1, rescale_percent=99)

print('Writing data')
scdata.magic.to_csv(args.out)
