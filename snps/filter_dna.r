library(optparse)

# get input arguments
options = list(
    make_option('--in', help='input dna matrix'),
    make_option('--M', help='number of total alleles per sample', default=0, type='numeric'),
    make_option('--m', help='number of minor alleles per sample', default=0, type='numeric'),
    make_option('--N', help='number of total alleles per snp', default=0, type='numeric'),
    make_option('--n', help='number of minor alleles per snp', default=0, type='numeric'),
    make_option('--out', help='output dna matrix')
)
args = parse_args(OptionParser(option_list=options))

# read data (rows = snps, columns = samples)
x = read.table(args$i, sep='\t', header=T)
x = t(x)

# filter matrix (rows = samples, columns = snps)
f = data.frame(x)
f = sapply(f, function(a){
    a[a == '-'] = NA
    a = droplevels(a)
    order(table(a), decreasing=T)[as.numeric(a)]
})
rownames(f) = colnames(x)

# filter samples
x = x[rowSums(f > 0, na.rm=T) >= args$M & rowSums(f > 1, na.rm=T) >= args$m,]

# filter snps
x = x[,colSums(f > 0, na.rm=T) >= args$N & colSums(f > 1, na.rm=T) >= args$n]

# fix rownames
rownames(x) = paste(rownames(x), ';', sep='')

# write data
write.table(x, file=args$out, sep='', quote=F, col.names=F)
