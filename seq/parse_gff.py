import argparse, re, util

def nice_name(x):
    return util.get_valid_filename(x)

# input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--gff', help='gff file')
parser.add_argument('--type', default='', help='type (CDS, gene, etc)')
args = parser.parse_args()

print('\t'.join(['id', 'contig', 'beg', 'end', 'gene', 'product']))

for line in util.open_gz(args.gff):

    # parse gff record
    # ----------------
    line = line.rstrip().split('\t')

    if args.type:
        if len(line) < 3 or line[2] != args.type:
            continue

    # coords
    contig = line[0]
    beg = min(int(line[3]), int(line[4]))
    end = max(int(line[3]), int(line[4]))
    
    # gene id
    try:
        gid = re.search('ID=(.*?);', line[8]).group(1)
    except:
        gid = ''
    
    # gene name
    try:
        gene = re.search('gene=(.*?);', line[8]).group(1)
    except:
        gene = ''
    
    # gene product
    try:
        product = re.search('product=(.*?);', line[8]).group(1)
    except:
        product = ''
        
    # print record
    print('\t'.join(map(str, [gid, contig, beg, end, gene, product])))
