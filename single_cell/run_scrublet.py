import argparse
import scipy.io
import scrublet

# This script runs Scrublet

# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--data', help='Sparse counts matrix (rows = cells, columns = features)')
parser.add_argument('--out', help='Scrublet outfile')
args = parser.parse_args()

# Read data
data = scipy.io.mmread(args.data)

# Run scrublet
scrub = scrublet.Scrublet(data)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

# Write results
out = open(args.out, 'w')
for i in range(len(doublet_scores)):
    x = ['NA', 'NA']
    try:
        x[0] = doublet_scores[i]
    except:
        pass
    try:
        x[1] = predicted_doublets[i]
    except:
        pass
    out.write('%s\t%s\n' %(x[0], x[1]))
out.close()
