
import argparse, gzip, os, re, sys, util


# input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--fna', help='fna file', default='')
parser.add_argument('--gtf', help='gtf file', default='')
parser.add_argument('--pre', help='output prefix', default='')
parser.add_argument('--map', help='map file (1=prefix, 2=fna, 3=gtf)', default='')
parser.add_argument('--tag', help='feature tag regex (gene|CDS|tRNA)', default='CDS')
parser.add_argument('--regex', help='gene regex', default='')
parser.add_argument('--fmt', help='output format', choices=['org', 'gene'], default='org')
parser.add_argument('--gmap', help='gene map (for fmt="gene")', default='')
parser.add_argument('--ming', help='min # organisms per gene', default=0, type=int)
parser.add_argument('--topg', help='select top genes (ranked by # organisms)', default=0, type=int)
parser.add_argument('--exon', help='exons only?', default=False, action='store_true')
args = parser.parse_args()


# check arguments
if not args.fna or not args.gtf:
    if not args.map:
        quit('error: must provide fna/gtf or map')
if args.regex or args.fmt == 'gene':
    args.exon = True


# read input files
imap = {}
if args.map:
    for line in open(args.map):
        pre, fna, gtf = line.rstrip().split()
        imap[pre] = [fna, gtf]
else:
    imap[args.pre] = [args.fna, args.gtf]

# read gene map
gmap = {}
if args.gmap:
    for line in open(args.gmap):
        line = line.rstrip().split('\t')
        gene = line[0]
        for gid in line[1].split(','):
            gmap[gid] = gene

# nice filenames
def nice_name(x):
    return util.get_valid_filename(x)
    
    
# gtf fields
# ----------

def extract_gene_id(line):
    m = re.search('gene_id "(.*?)"', line[8])
    if not m:
        m = re.search('ID=(.*?);', line[8])
    if m:
        return nice_name(m.group(1))
    else:
        return ''


def extract_gene_name(line):
    m = re.search('gene "(.*?)"', line[8])
    if not m:
        m = re.search('gene=(.*?);', line[8])
    if m:
        return nice_name(m.group(1))
    else:
        return ''


# filter genes
def filter_genes():
    count = {}
    if args.topg == 0 and args.ming == 0:
        return([])
    for pre in imap:
        gtf = imap[pre][1]
        for line in util.open_gz(gtf):
            line = line.rstrip().split('\t')
            try:
                #gid = nice_name(re.search('gene_id "(.*?)"', line[8]).group(1))
                #gname = nice_name(re.search('gene "(.*?)"', line[8]).group(1))
                gid = extact_gene_id(line)
                gname = extract_gene_name(line)
                if gname not in count:
                    count[gname] = []
                count[gname].append(pre)
            except:
                continue
    for gi in count:
        count[gi] = len(set(count[gi]))
    genes = sorted(count.keys(), key=lambda x: count[x], reverse=True)
    if args.topg > 0:
        genes = genes[:args.topg]
    if args.ming > 0:
        genes = [x for x in genes if count[x] >= args.ming]
    print([count[gi] for gi in genes])
    print('Min: %s organisms' %(count[genes[-1]]))
    print('Max: %s organisms' %(count[genes[0]]))
    return(genes)

def filter_gmap():
    genes = []
    count = {}
    if args.topg == 0 and args.ming == 0:
        return([])
    for line in open(args.gmap):
        line = line.rstrip().split('\t')
        gene = line[0]
        gids = line[1].split(',')
        if len(gids) >= args.ming:
            genes.append(gene)
            count[gene] = len(gids)
    genes = sorted(count.keys(), key=lambda x: count[x], reverse=True)
    print(len(genes))
    if args.topg > 0:
        genes = genes[:args.topg]
    print(len(genes))
    if args.ming > 0:
        genes = [x for x in genes if count[x] >= args.ming]
    print(len(genes))
    print([count[gi] for gi in genes])
    print('Min: %s organisms' %(count[genes[-1]]))
    print('Max: %s organisms' %(count[genes[0]]))
    return(genes)

if args.gmap:
    genes = filter_gmap()
else:
    genes = filter_genes()

if len(genes) > 0:
    print('Selecting %s genes:\n%s' %(len(genes), ', '.join(genes)), flush=True)

# write output
fhs = {'': sys.stdout}
def writeout(seq, org, gene):
    ki = ''
    if args.fmt == 'org':
        ki = org
    elif args.fmt == 'gene':
        ki = gene
    else:
        quit('error: invalid --fmt "%s" argument' %(args.fmt))
    if ki not in fhs:
        fn = '%s.fna' %(ki)
        if os.path.exists(fn):
            quit('error: file %s exists' %(fn))
        fhs[ki] = open(fn, 'w')
    fh = fhs[ki]
    fh.write('%s\n' %(seq))


# iterate over files
for pre in imap:
    fna, gtf = imap[pre]
    
    # initialize vars
    ibeg = 0
    ipre = ''
    ifin = ''
    
    # read fasta file
    # keys = contig names
    fna = util.read_fst(fna, split=True)
    
    # iterate over gtf entries
    for line in util.open_gz(gtf):
        
        # skip header rows
        if line.startswith('#'):
            continue
        
        # extract gtf record
        line = line.rstrip().split('\t')
        #if line[2] in ['gene', 'CDS', 'tRNA']:
        #if line[2] in ['gene', 'CDS', 'tRNA']:
        m = re.search(args.tag, line[2])
        if m:
            
            # gene information
            contig = line[0]
            beg = int(line[3])
            end = int(line[4])
            strand = line[6]
            #gid = nice_name(re.search('gene_id "(.*?)"', line[8]).group(1))
            #try:
            #    gname = nice_name(re.search('gene "(.*?)"', line[8]).group(1))
            #except:
            #    gname = ''
            gid = extract_gene_id(line)
            if args.gmap:
                if gid in gmap:
                    gname = gmap[gid]
                else:
                    continue
            else:
                gname = extract_gene_name(line)
                
            # filter genes
            if len(genes) > 0:
                if gname not in genes:
                    continue
            
            # get intron
            intron = fna[contig][ibeg:(beg-1)]
            ibeg = end
            if ipre and not args.exon and len(intron) > 0:
                intron_seq = re.sub('^>\.', '>', '>%s.intron_%s %s\n%s' %(pre, ipre, gname, intron))
                writeout(intron_seq, org=pre, gene=gname)
            ipre = gid
            if ifin == '':
                ifin = beg
            
            # get gene
            gene = fna[contig][(beg-1):end]
            if strand == '-':
                gene = util.reverse_complement(gene)
            if args.regex and not re.search(args.regex, gname):
                continue
            if len(gene) > 0:
                exon_seq = re.sub('^>\.', '>', '>%s.%s %s\n%s' %(pre, gid, gname, gene))
                writeout(exon_seq, org=pre, gene=gname)

# final intron
intron = fna[contig][ibeg:] + fna[contig][:(ifin-1)]
if not args.exon and len(intron) > 0:
    intron_seq = re.sub('^>\.', '>', '>%s.intron_%s %s\n%s' %(pre, ipre, gname, intron))
    writeout(intron_seq, org=pre, gene=gname)

for fh in fhs:
    try:
        fh.close()
    except:
        pass
