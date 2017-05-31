import re

# Convert a list of RDP IDs to a newick tree

def rdp2tree(x):
    # Convert a list of RDP IDs to a python tree structure
	# Input: list of RDP IDs, e.g.
	# - Bacteria|Proteobacteria|Gammaproteobacteria|Pasteurellales|Pasteurellaceae|Actinobacillus
	# - Bacteria|Proteobacteria|Gammaproteobacteria|Pseudomonadales|Moraxellaceae|Acinetobacter
	# Output: tree (python dictionary)
	# - tree[node] = {child1:{...}, child2:{...}, ...}
	# - tree[leaf] = {}
    tree = {}
    for otu in x:
        nodes = otu.split('|')
        ctree = tree
        for node in nodes:
            node = re.sub('[ :;\(\)\[\]]', '', node)
            if node not in ctree:
                ctree[node] = {}
            ctree = ctree[node]
    return tree


def add_leaf_ids(tree, leaf_id=[]):
    # Add leaf ids to the tips of a tree
    # Input: tree (python dictionary)
    # - tree[node] = {child1:{...}, child2:{...}, ...}
    # - tree[node] = {}
    # Output: relabeled tree (python dictionary)
    # - tree[node] = {child1:{...}, child2:{...}, ...}
    # - tree[leaf] = RDP ID
    for child in tree:
        child_leaf_id = leaf_id + [child]
        if len(tree[child]) == 0:
            tree[child] = '|'.join(child_leaf_id)
        else:
            add_leaf_ids(tree[child], child_leaf_id)
    return tree


def tree2newick(tree, leaf_id=[], format='normal'):
    # Convert tree (python dictionary) to Newick string
    # Input: tree (python dictionary)
    # Output: newick string
    if type(tree) == str:
        return '%s' %(tree)
    newick = []
    for child in tree:
        child_leaf_id = leaf_id + [child]
        newick.append(tree2newick(tree[child], leaf_id=child_leaf_id, format=format))
    if len(leaf_id) == 0:
        node_id = 'Root'
    else:
        if format == 'normal':
            node_id = leaf_id[-1]
        if format == 'full':
            node_id = '|'.join(leaf_id)
    newick = '(%s)%s' %(','.join(newick), node_id)
    return newick


def run():
    # Convert RDP list to newick file
    
    # Get command line arguments
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='RDP list (input)')
    parser.add_argument('-o', help='Newick tree (output)')
    parser.add_argument('-f', help='Format of internal node labels', choices=['normal', 'full'])
    args = parser.parse_args()
    
    # Get RDP list
    rdp = [line.rstrip() for line in open(args.i)]
    
    # Convert to tree
    tree = rdp2tree(rdp)
    tree = add_leaf_ids(tree)
    
    # Write to newick
    newick = tree2newick(tree, format=args.f)
    out = open(args.o, 'w')
    out.write(newick+';\n')
    out.close()


if __name__ == '__main__':
    run()
