library(optparse)

option_list = list(make_option('--folder', help='directory name'),
                   make_option('--ming', help='genes per cell cutoff', type='integer'),
                   make_option('--out', help='output prefix'),
                   make_option('--num', help='number of downsamples', default=10, type='integer')
                   )
args = parse_args(OptionParser(option_list=option_list))

library(vegan)

# Reads per collection
x.r = c()

# Cells per collection
y.h = c()
y.m = c()
y.t = c()

# Get folder name
folder = args$folder

print('Loading human DGE')
h.dge = read.table(paste('~/aviv/data/', folder, '/dge/human.dge.txt.gz', sep=''), header=T, row.names=1, sep='\t')

print('Loading human reads')
h.reads = read.table(paste('~/aviv/data/', folder, '/dge/human.reads.dge.txt.gz', sep=''), header=T, row.names=1, sep='\t')

print('Loading mouse DGE')
m.dge = read.table(paste('~/aviv/data/', folder, '/dge/mouse.dge.txt.gz', sep=''), header=T, row.names=1, sep='\t')

print('Loading mouse reads')
m.reads = read.table(paste('~/aviv/data/', folder, '/dge/mouse.reads.dge.txt.gz', sep=''), header=T, row.names=1, sep='\t')

# Count total reads
total_reads = 0
if(ncol(h.dge) > 1){total_reads = total_reads + sum(h.reads)}
if(ncol(m.dge) > 1){total_reads = total_reads + sum(m.reads)}

# Get read counts
min_reads = 10**6
max_reads = total_reads
x.i = seq(log10(min_reads), log10(max_reads), length.out=args$num)

# Count collections
ncoll = 0

# Subsample human reads
if(ncol(h.dge) > 1){
    h.coll = length(unique(sapply(strsplit(colnames(h.dge), '\\.'), '[', 1)))
    ncoll = h.coll
    for(xi in x.i){
        new_reads = rrarefy(h.reads, round((10**xi)*rowSums(h.reads)/sum(h.reads)))
        new_dge = matrix(0, nrow=nrow(h.dge), ncol=ncol(h.dge))
        subsample_dge = function(i,j){
            x1 = h.dge[i,j]
            y1 = h.reads[i,j]
            y2 = new_reads[i,j]
            r1 = rmultinom(1, y1-x1, rep(1,x1)) + 1
            r2 = rrarefy(r1, y2)
            x2 = sum(r2 > 0)
            return(x2)
        }
        subsample_dge = Vectorize(subsample_dge)
        ijs = which(new_reads > 0, arr.ind=T)
        new_dge[ijs] = apply(ijs, 1, function(a){subsample_dge(a[[1]], a[[2]])})
        h.cells = sum(colSums(new_dge > 0) >= args$ming)
        y.h = c(y.h, h.cells/h.coll)
    }
} else {
    y.h = c(y.h, rep(0, length(x.i)))
}

# subsample mouse reads
if(ncol(m.dge) > 1){
    m.coll = length(unique(sapply(strsplit(colnames(m.dge), '\\.'), '[', 1)))
    ncoll = m.coll
    for(xi in x.i){
        new_reads = rrarefy(m.reads, round((10**xi)*rowSums(m.reads)/sum(m.reads)))
        new_dge = matrix(0, nrow=nrow(m.dge), ncol=ncol(m.dge))
        subsample_dge = function(i,j){
            x1 = m.dge[i,j]
            y1 = m.reads[i,j]
            y2 = new_reads[i,j]
            r1 = rmultinom(1, y1-x1, rep(1,x1)) + 1
            r2 = rrarefy(r1, y2)
            x2 = sum(r2 > 0)
            return(x2)
        }
        subsample_dge = Vectorize(subsample_dge)
        ijs= which(new_reads > 0, arr.ind=T)
        new_dge[ijs] = apply(ijs, 1, function(a){subsample_dge(a[[1]], a[[2]])})
        m.cells = sum(colSums(new_dge > 0) >= args$ming)
        y.m = c(y.m, m.cells/m.coll)
    }
} else{
    y.m = c(y.m, rep(0, length(x.i)))
}

# reads per collection
x.r = c(x.r, x.i/ncoll)

# total genes per collection
y.t = y.h + y.m

# save output to file
z = list(x.r=x.r, y.h=y.h, y.m=y.m, y.t=y.t)
saveRDS(z, paste(args$out, 'overseq.rds', sep='.'))
