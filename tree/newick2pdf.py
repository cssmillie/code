import argparse, os, re
from itolapi import Itol
from itolapi import ItolExport

# Read input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--tree', help='newick tree file', required=True)
parser.add_argument('--boot', help='collapse nodes with bootstrap < boot', default=0)
parser.add_argument('--root', help='root outgroup (or midpoint)', default='')
parser.add_argument('--mode', help='display mode', choices=['normal', 'circular'], default='normal')
parser.add_argument('--out', help='pdf outfile', required=True)
args = parser.parse_args()

# Get root
if args.root == 'midpoint':
    root = '1'
else:
    root = args.root

# Modify tree
cmd = 'java -Xmx16m -jar ~/bin/TreeCollapseCL4.jar -f %s -b %s -r %s' %(args.tree, args.boot, root)
os.system(cmd)
prefix, suffix = re.search('(.*?)\.(.*?)$', args.tree).groups()
tree = '%s_%scoll.%s' %(prefix, args.boot, suffix)
args.tree = '%s.tmp' %(args.tree)
cmd = 'mv %s %s' %(tree, args.tree)
os.system(cmd)

# Upload tree
uploader = Itol.Itol()
uploader.add_variable('treeFile', args.tree)
uploader.add_variable('treeFormat', 'newick')
uploader.add_variable('treeName', args.out)
if args.root == 'midpoint':
    uploader.add_variable('midpointRoot', '1')
uploader.upload()
tree_id = uploader.comm.tree_id

# Export tree
exporter = ItolExport.ItolExport()
exporter.set_export_param_value('tree', tree_id)
exporter.set_export_param_value('format', 'pdf')
exporter.set_export_param_value('displayMode', args.mode)
exporter.export(args.out)
