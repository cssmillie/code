import argparse, re, util

def nice_name(x):
    return util.get_valid_filename(x)

# input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--gtf', help='gtf file')
parser.add_argument('--pre', help='prefix')
parser.add_argument('--tab', help='gtf list (1=prefix, 2=gtf)')
parser.add_argument('--order', help='print gene order', default=False, action='store_true')
args = parser.parse_args()

# gtf filenames
gtfs = {}
if args.gtf:
    gtfs[args.pre] = args.gtf
elif args.tab:
    for line in open(args.tab):
        pre, gtf = line.rstrip().split()
        gtfs[pre] = gtf
else:
    quit('error: must provide --gtf or --tab')

# iterate over gtfs
if args.order:
    res = {}
    for pre in gtfs:
        gtf = gtfs[pre]
        res[pre] = []
        for line in util.open_gz(gtf):
            line = line.rstrip().split('\t')
            try:
                gene = nice_name(re.search('gene "(.*?)"', line[8]).group(1))
                if gene not in res[pre]:
                    res[pre].append(gene)
            except:
                continue
    for pre in res:
        if len(res[pre]) > 100:
            print('%s\t%s' %(pre, ','.join(res[pre])), flush=True)
        
