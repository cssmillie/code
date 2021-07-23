import argparse, re, util

# Make kpileup "gene file" from FASTA reference
# Output columns:
# 1 = Contig ID
# 2 = Gene ID
# 3 = Beg position (1-indexed)
# 4 = End position (1-indexed)
# 5 = Strand (+ or -)
# 6 = Sequence

# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--fst', help='Input FST file')
parser.add_argument('--gtf', help='Input GTF file')
parser.add_argument('--tag', help='GTF tag', default='gene')
parser.add_argument('--names', help='use gene names (default = gene ids)', default=False, action='store_true')
args = parser.parse_args()

# Read FST sequence
fst = {}
for record in util.iter_fst(args.fst):
    sid = record[0].split()[0][1:]
    seq = record[1]
    fst[sid] = seq

# Parse GTF file
for line in open(args.gtf):

    # Skip header
    if line.startswith('#'):
        continue
    
    # Extract fields    
    line = line.rstrip().split('\t')
    contig = line[0]
    tag = line[2]
    beg = int(line[3])
    end = int(line[4])
    strand = line[6]
    anno = line[8]
    gene = ''

    # Check tag
    if tag != args.tag:
        continue
    
    # Get gene name
    gene = re.search('gene_id "(.*?)"', anno).group(1)
    if args.names:
        try:
            gene = re.search('gene "(.*?)"', anno).group(1)
        except:
            pass
    
    # Get sequence
    seq = fst[contig][(beg-1):end]
    
    # Print line
    print('\t'.join(map(str, [contig, gene, beg, end, strand, seq])))
