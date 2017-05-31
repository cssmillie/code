library(ape)
library(optparse)

read_input = function(){
  option_list= list(
    make_option(c('--tree'), default='', help='Infile (newick tree)'),
    make_option(c('--otus'), default='', help='OTUs to keep (list'),
    make_option(c('--out'), default='', help='Outfile (newick tree)')    
  )
  args = parse_args(OptionParser(option_list = option_list), args=commandArgs(trailingOnly=TRUE))
  return(args)
}

args = read_input()
tree = read.tree(args$tree)
otus = as.character(read.csv(args$otus, header=F)[,1])
otus.drop = setdiff(tree$tip.label, otus)
tree = drop.tip(tree, otus.drop)
write.tree(tree, file=args$out)
