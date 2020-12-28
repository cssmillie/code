import sys, util

# input arguments
fna = sys.argv[1]
gtf = sys.argv[2]

# iterate over gtf
for line in open(gtf):
    line = line.rstrip().split('\t')
    
