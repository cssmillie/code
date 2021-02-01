import argparse

# read input arguments
parser = argparser.ArgumentParser()
parser.add_argument('--sam', help='SAM file')
parser.add_argument('--gtf', help='GTF file')
parser.add_argument('--pct', help='min %id', default=90, type=float)
parser.add_argument('--len', help='min length', default=25, type=float)
parser.add_argument('--out', help='outfile')
args = parser.parse_args()


