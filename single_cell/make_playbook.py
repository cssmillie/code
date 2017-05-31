import argparse, os

parser = argparse.ArgumentParser()
parser.add_argument('--name', help='Project name', default='test')
parser.add_argument('--dge', help='Input DGE matrix', default='', required=True)
parser.add_argument('--regex', help='Select cells by regex', default='.*')
parser.add_argument('--ming', help='Min genes per cell', default=500, type=int)
parser.add_argument('--minc', help='Min cells per gene', default=0, type=int)
parser.add_argument('--maxc', help='Max cells per group', default=0)
parser.add_argument('--varmet', help='Variable gene selection', default='karthik')
parser.add_argument('--sep', help='Identity field separator', default=r'\\\\.')
parser.add_argument('--field', help='Identity field', default=1, type=int)
parser.add_argument('--batch', help='Batch correction', default=0, type=int)
parser.add_argument('--cores', help='Number of cores', default=1, type=int)
parser.add_argument('--template', help='Playbook template', default='~/code/single_cell/playbook_template.r')
args = parser.parse_args()

import os, re, ssub
import numpy as np

if args.dge != '':
    dge = os.path.abspath(args.dge)
dge = re.sub('/', '\\/', dge)

cmd = '''
sed -e 's/NAME/%s/' \
    -e 's DGE_FILE %s ' \
    -e 's/REGEX/%s/' \
    -e 's/MING/%s/' \
    -e 's/MINC/%s/' \
    -e 's/MAXC/%s/' \
    -e 's/VARMET/%s/' \
    -e 's/SEP/%s/' \
    -e 's/FIELD/%s/' \
    -e 's/BATCH/%s/' \
    -e 's/CORES/%s/' \
    %s | grep -v '^#'
''' %(args.name, args.dge, args.regex, args.ming, args.minc, args.maxc, args.varmet, args.sep, args.field, args.batch, args.cores, args.template)

os.system(cmd)
