
# command line arguments
infile = commandArgs(trailingOnly=T)[[1]]
outfile = commandArgs(trailingOnly=T)[[2]]

min_genes = 1000

# read dge
x = read.table(infile, sep='\t', header=T, row.names=1)
g = colSums(x > 0)
x = x[, g >= min_genes]

# write umis
cat(colnames(x), file=outfile, sep='\n')
