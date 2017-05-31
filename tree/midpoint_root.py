import dendropy, sys

t = dendropy.Tree.get_from_path(sys.argv[1], 'newick')
t.reroot_at_midpoint()
t.write_to_path(sys.argv[2], 'newick')
