library(argparse)
source('~/code/util/mtx.r')

if(interactive()){

    # set default arguments
    # ---------------------
    args = list(kraken='CD7017.Epi_A.kraken.out')

} else {

    # read input arguments
    # --------------------
    parser = ArgumentParser()
    parser$add_argument('--kraken', help='kraken outfile')
    parser$add_argument('--report', help='kraken report')
    parser$add_argument('--prefix', help='sample prefix')
    parser$add_argument('--out', help='dge output file')
    args = parser$parse_args()

}

# map taxa to names
y = ffread(args$report, sep='\t', header=F, as.dt=T)
y = setNames(y[[6]], y[[5]])

# load kraken data
x = ffread(args$kraken, sep='\t', header=F, as.dt=T)
x = x[x[[1]] == 'C']
u = sapply(strsplit(x[[2]], ';'), '[[', 2) # barcodes
v = as.factor(x[[3]]) # taxa
x = do.call(rbind, tapply(v, u, table))

# fix duplicates
print(any(duplicated(rownames(x))))
print(any(duplicated(colnames(x))))

# fix names
rownames(x) = paste(args$prefix, rownames(x), sep='.')
colnames(x) = paste(colnames(x), gsub(' ', '_', y[colnames(x)]), sep='.')

# write output file
write.table(x, file=args$out, sep='\t', quote=F)
