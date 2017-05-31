import argparse, re

# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--gene', help='gene id')
parser.add_argument('--genes', help='list of gene ids')
parser.add_argument('--col', help='column', type=int, default=0)
parser.add_argument('--replace', help='replace?', default=False, action='store_true')
parser.add_argument('--mouse', help='use mouse target & db', default=False, action='store_true')
parser.add_argument('--multi', help='multiple mappings', default=False, action='store_true')
parser.add_argument('--full_na', help='full NAs', default=False, action='store_true')
parser.add_argument('--skip_na', help='skip NAs', default=False, action='store_true')
parser.add_argument('--db', help='gene info database', default='/home/unix/csmillie/aviv/db/human.gene_info.txt')
parser.add_argument('--target', help='target list', default='/home/unix/csmillie/aviv/db/human_genes.txt')
args = parser.parse_args()

# Get mouse files
if args.mouse:
    args.db = '/home/unix/csmillie/aviv/db/mouse.gene_info.txt'
    args.target = '/home/unix/csmillie/aviv/db/mouse_genes.txt'

# Fix case
tcase = {} # map gene.upper() -> gene
qcase = {} # map gene.upper() -> gene
def fix_gene(gene):
    gene = gene.upper()
    return gene


# Get target genes
target = {}
for line in open(args.target):
    gene = line.rstrip()
    GENE = fix_gene(gene)
    target[GENE] = 1
    tcase[GENE] = gene

# Get query genes
query = {}
if args.gene:
    gene = args.gene
    GENE = fix_gene(gene)
    query[GENE] = []
    qcase[GENE] = gene

if args.genes:
    for line in open(args.genes):
        gene = line.rstrip().split()[args.col]
        GENE = fix_gene(gene)
        query[GENE] = []
        qcase[GENE] = gene

# Populate query
for GENE in query:
    if GENE in target:
        query[GENE] = [GENE]

# Read gene info
for line in open(args.db):

    if line.startswith('#'):
        continue

    line = line.rstrip().split('\t')
    genes = [line[2]] + line[4].split('|')
    GENES = [fix_gene(gene) for gene in genes]

    for QGENE in GENES:
        if QGENE in query:
           
            if args.multi == False and len(query[QGENE]) > 0:
                continue
            
            for TGENE in GENES:
                if TGENE in target:
                    query[QGENE].append(TGENE)
                    
                    if args.multi == False:
                        break

# Print map
if args.gene:
    qgene = args.gene
    QGENE = fix_gene(args.gene)
    if len(query[QGENE]) > 0:
        TGENES = query[QGENE]
        tgenes = ','.join([tcase[TGENE] for TGENE in TGENES])
        print '%s\t%s' %(qgene, tgenes)
    else:
        if args.skip_na == False:
            if args.full_na == True:
                tgene = '%s:NA' %(qgene)
            else:
                tgene = 'NA'
            print '%s\t%s' %(qgene, tgene)

if args.genes:
    for line in open(args.genes):
        line = line.rstrip().split()
        qgene = line[args.col]
        QGENE = fix_gene(qgene)
        if len(query[QGENE]) > 0:
            TGENES = query[QGENE]
            tgenes = ','.join(set([tcase[TGENE] for TGENE in TGENES]))
            if args.replace == True:
                line[args.col] = tgenes
            else:
                line.insert(args.col+1, tgenes)
        else:
            if args.skip_na == True:
                continue
            if args.full_na == True:
                if args.replace == True:
                    line[args.col] = '%s:NA' %(qgene)
                else:
                    line.insert(args.col+1, '%s:NA' %(qgene))
            else:
                if args.replace == True:
                    line[args.col] = 'NA'
                else:
                    line.insert(args.col+1, 'NA')

        print '\t'.join(line)
