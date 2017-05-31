import argparse, copy, re

# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--gene', help='Gene ID')
parser.add_argument('--genes', help='List of gene IDs (txt)')
parser.add_argument('--col', help='Column number (origin = 1)', type=int, default=1)
parser.add_argument('--replace', help='Replace original ID (default = insert)?', default=False, action='store_true')
parser.add_argument('--db1', help='Species 1 gene info database', default='/home/unix/csmillie/aviv/db/human.gene_info.txt')
parser.add_argument('--db2', help='Species 2 gene info database', default='/home/unix/csmillie/aviv/db/mouse.gene_info.txt')
parser.add_argument('--target', help='List of target genes', default='/home/unix/csmillie/aviv/db/mouse_genes.txt')
parser.add_argument('--ortho', help='Orthologs map (tsv)', default='/home/unix/csmillie/aviv/db/orthologs.h2m.txt')
parser.add_argument('--multi', help='Allow multiple mappings (comma-separated)', default=False, action='store_true')
parser.add_argument('--skip_na', help='Skip NAs', default=False, action='store_true')
parser.add_argument('--m2h', help='Map mouse to human', default=False, action='store_true')
args = parser.parse_args()

# Get mapping files
if args.m2h:
    args.db1 = '/home/unix/csmillie/aviv/db/mouse.gene_info.txt'
    args.db2 = '/home/unix/csmillie/aviv/db/human.gene_info.txt'
    args.target = '/home/unix/csmillie/aviv/db/human_genes.txt'
    args.ortho = '/home/unix/csmillie/aviv/db/orthologs.m2h.txt'

# Fix case
qcase = {} # map gene.upper() -> gene
tcase = {} # map gene.upper() -> gene
def fix_gene(gene):
    gene = gene.upper()
    return gene

# Get query genes (default=human)
query = {}
if args.gene:
    GENE = fix_gene(gene)
    query[GENE] = []
    qcase[GENE] = gene
if args.genes:
    for line in open(args.genes):
        gene = line.rstrip().split()[(args.col-1)]
        GENE = fix_gene(gene)
        query[GENE] = []
        qcase[GENE] = gene

# Get target genes (default=mouse)
target = {}
for line in open(args.target):
    gene = line.rstrip()
    GENE = fix_gene(gene)
    target[GENE] = 1
    tcase[GENE] = gene

# Get ortholog map
ortho = {}
for line in open(args.ortho):
    [gene1,gene2] = line.rstrip().split()
    GENE1 = fix_gene(gene1)
    GENE2 = fix_gene(gene2)
    ortho[GENE1] = GENE2

def map_genes(query, target, db):
    # Map query to target using synonyms in db
    
    m = copy.deepcopy(query)
    
    # Check for matching names
    for GENE in query:
        if GENE in target:
            m[GENE] = [GENE]
    
    # Use synonyms from db
    for line in open(db):
        line = line.rstrip().split('\t')
        
        # Skip header
        if line[0].startswith('#'):
            continue
        
        # Get gene synonyms
        if len(line) == 2:
            genes = line
        else:
            genes = [line[2]] + line[4].split('|')
        GENES = [fix_gene(gene) for gene in genes]
        
        # Map query to target
        for QGENE in GENES:
            if QGENE in query:
                for TGENE in GENES:
                    if TGENE in target:
                        if args.multi == False and len(m[QGENE]) > 0:
                            continue
                        m[QGENE].append(TGENE)
    
    # Return gene map
    return m


# 1: map query genes to human ortholog keys
map1 = map_genes(query, ortho, args.db1)

# 2: map human orthologs to mouse orthologs
map2 = {}
for u in map1:
    for v in map1[u]:
        map2[v] = ortho[v]

# 3: map mouse orthologs to target genes
map3 = dict([(v, []) for v in map2.values()])
map3 = map_genes(map3, target, args.db2)

# Merge maps: query -> target
m = {}
for qgene in map1:
    m[qgene] = []
    hgenes = map1[qgene]
    for hgene in hgenes:
        mgene = map2[hgene]
        tgenes = map3[mgene]
        for tgene in tgenes:
            m[qgene].append(tgene)

# Print map

if args.gene:
    qgene = args.gene
    QGENE = fix_gene(args.gene)
    if len(m[QGENE]) > 0:
        TGENES = m[QGENE]
        tgenes = ','.join([tcase[TGENE] for TGENE in TGENES])
        print '%s\t%s' %(qgene, tgenes)
    else:
        if args.skip_na == False:
            tgene = 'NA'
            print '%s\t%s' %(qgene, tgene)

if args.genes:
    for line in open(args.genes):
        line = line.rstrip().split()
        qgene = line[(args.col-1)]
        QGENE = fix_gene(qgene)
        if len(m[QGENE]) > 0:
            TGENES = m[QGENE]
            tgenes = ','.join(set([tcase[TGENE] for TGENE in TGENES]))
            if args.replace == True:
                line[(args.col-1)] = tgenes
            else:
                line.insert(args.col, tgenes)
        else:
            if args.skip_na == False:
                if args.replace == True:
                    line[(args.col-1)] = 'NA'
                else:
                    line.insert(args.col, 'NA')

        print '\t'.join(line)
