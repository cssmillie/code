import argparse, ete3

# parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-t', help='newick tree')
parser.add_argument('-l', help='labels')
parser.add_argument('-r', help='remove internal nodes', action='store_true', default=False)
parser.add_argument('-o', help='output newick')
args = parser.parse_args()

# read input tree
tree = ete3.Tree(args.t, format=1)

# get labels to keep
labels = [line.rstrip().split()[0] for line in open(args.l)]
labels = dict(zip(labels, [1]*len(labels)))

if args.r == False:
    
    # subset tree
    remove = []
    for node in tree.traverse('levelorder'):
        descendants = [node.name] + [desc.name for desc in node.get_descendants()]
        n_keep = sum(map(int, [desc in labels for desc in descendants]))
        if n_keep == 0:
            remove.append(node)
    
    # remove nodes
    for node in remove:
        node.detach()

if args.r == True:
    
    # intersect labels
    descendants = [desc.name for desc in tree.get_descendants()]
    labels = [desc for desc in descendants if desc in labels]
    
    # prune tree
    tree.prune(labels)


# write tree
tree.write(outfile=args.o, format=1)
