library(Seurat)
library(methods)
library(optparse)
library(Hmisc)

# Get command line arguments
if(interactive()){
    args = list()
    args$seur = '~/aviv/human/tsne/human.1000.10.seur.rds'
    args$log = '~/aviv/human/tsne/human.1000.10.log.rds'
    args$ming = 1000
    args$data = TRUE
    args$pca = FALSE
    args$go = FALSE
    args$goa = FALSE
    args$scale = FALSE
    args$acut = -1
    args$gcut = -1
    args$icut = 4
    args$cells = 100
    args$clust = TRUE
    args$tab = 'test.tab'
    args$ident = '~/aviv/human/cluster/human.1000.10.infomap.n25.clust.rds'

} else {
options = list(make_option('--seur', help='seurat file'),
               make_option('--log', help='log file'),
               make_option('--ming', help='minimum number of genes', default=500, type='integer'),
               make_option('--data', help='use data', action='store_true', default=FALSE),
               make_option('--pca', help='use pca', action='store_true', default=FALSE),
               make_option('--go', help='use go', action='store_true', default=FALSE),
               make_option('--goa', help='use go + genes', action='store_true', default=FALSE),
               make_option('--scale', help='scale data', action='store_true', default=FALSE),
               make_option('--acut', help='binary cutoff', default=-1, type='numeric'),
               make_option('--gcut', help='number of groups (discretize each sample)', default=-1, type='integer'),
               make_option('--icut', help='number of groups (discretize matrix)', default=-1, type='integer'),
               make_option('--cells', help='number of cells a gene must be present in', default=0, type='integer'),
               make_option('--clust', help='aggregate by clusters', action='store_true', default=FALSE),
               make_option('--tab', help='name of alignment file'),
               make_option('--ident', help='clusters list (1=cell, 2=group)', default='')
               )
args = parse_args(OptionParser(option_list=options))
}

# Read Seurat object
seur = readRDS(args$seur)
log = readRDS(args$log)

# Select cells
cells.use = colnames(seur@raw.data)[colSums(seur@raw.data > 0) >= args$ming]
cells.use = cells.use[cells.use %in% colnames(seur@data)]

# Select data
if(args$pca == TRUE){
    npcs = min(log$pc_sig, ncol(seur@pca.rot))
    data = seur@pca.rot[cells.use, 1:npcs]
} else if(args$go == TRUE){
    go = read.table('/home/unix/csmillie/aviv/db/gopca/go_annotations_human.tsv', sep='\t', row.names=1, stringsAsFactors=F)
    go = go[go[,2] == 'BP',]
    data = sapply(strsplit(go[,4],','), function(a){colSums(seur@data[a,], na.rm=T)})
} else if(args$goa == TRUE){
    go = read.table('/home/unix/csmillie/aviv/db/gopca/go_annotations_human.tsv', sep='\t', row.names=1, stringsAsFactors=F)
    go = go[go[,2] == 'BP',]
    data = sapply(strsplit(go[,4],','), function(a){colSums(seur@data[a,], na.rm=T)})
    genes.go = unique(unlist(strsplit(go[,4], ',')))
    genes.ad = rownames(seur@data)[! rownames(seur@data) %in% genes.go]
    data = cbind(data, t(seur@data[genes.ad,]))
} else {
    data = t(seur@data[,cells.use])
}

# Filter zeros
data = data[,colSums(data > 0) >= args$cells]

# Get alignment
if(args$ident != ''){
    ident.use = readRDS(args$ident)$membership
    seur = set.ident(seur, ident.use=ident.use)
}
if(args$clust == FALSE){
    aln = data
} else{
    aln = data.frame(aggregate(data, list(seur@ident[cells.use]), mean), row.names=1)
}

# Scale
if(args$scale == TRUE){
    aln = t(scale(t(aln)))
}

# Discretize
if(args$acut != -1){
    aln = aln > args$acut
    aln = data.frame(t(apply(aln, 1, as.integer)))
}
if(args$gcut != -1){
    aln = data.frame(t(apply(aln, 1, function(a){as.integer(cut2(a, g=args$gcut))-1})))
}
if(args$icut != -1){
    if(min(aln) < 0){stop()}
    r = rownames(aln)
    x = as.numeric(as.matrix(aln))
    x = x[x > 0]
    p = seq(0, 1, 1/(args$icut-1))
    q = c(-1, quantile(x, prob=p))
    aln = data.frame(apply(aln, 2, function(a){as.integer(cut(a, breaks=q))-1}))
    rownames(aln) = r
}

# Format
rownames(aln) = paste(rownames(aln), ";", sep='')

# Write alignment
write.table(aln, file=args$tab, sep='', quote=F, col.names=F)
