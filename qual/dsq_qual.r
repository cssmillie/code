library(optparse)

# Dropseq quality metrics

# Input arguments
option_list = list(make_option('--folder', help='project folder'),
                   #make_option('--cutoff', help='genes per cell', type='numeric', default=500),
		   make_option('--out', help='output prefix')
              )
args = parse_args(OptionParser(option_list=option_list))

# Get folders
qc = file.path('~/aviv/data', args$folder, 'raw', 'QC_files')
print(qc)
dge = file.path('~/aviv/data', args$folder, 'raw', 'UMI_DGE')
print(dge)

# Get collections
collections = gsub('_bq.*', '', list.files(dge, pattern='*.HUMAN.*summary.txt'))
print(collections)

# Make filenames
get_filename = function(collection, type=''){
    if(type == 'hcells'){
        return(file.path(dge, paste(collection, '_bq10_star.HUMAN.umi.dge.summary.txt', sep='')))
    }
    if(type == 'hreads'){
        return(file.path(qc, paste(collection, '_bq10_star_numReads_perCell_ZC_mq_10.txt.gz', sep='')))
    }
    if(type == 'mcells'){
        return(file.path(dge, paste(collection, '_bq10_star.MOUSE.umi.dge.summary.txt', sep='')))
    }
    if(type == 'mreads'){
        return(file.path(qc, paste(collection, '_bq10_star_numReads_perCell_ZC_mq_10.txt.gz', sep='')))
    }
    if(type == 'reads'){
        return(file.path(qc, paste(collection, '_bq10_star_ReadQualityMetrics.txt', sep='')))
    }
    if(type == 'genic'){
        return(file.path(qc, paste(collection, '_bq10_star_fracIntronicExonic.txt', sep='')))
    }
    if(type == 'doublet'){
        return(list.files(qc, pattern=paste(collection, '_.*_categorized_cellTypes.txt', sep=''), full.names=T))
    }
}

#print('human cells')
#h.cells = sapply(collections, function(a){
#    q = get_filename(a, 'hcells')
#    q = read.table(q, sep='\t', header=T, row.names=1)
#    q = q[q$NUM_GENES >= args$cutoff,]
#    num_cells = nrow(q)
#    num_genes = q[,1]
#    num_umis = q[,2]
#    return(list(num_cells=num_cells, num_genes=num_genes, num_umis=num_umis))
#})
#h.cells.names = c('Human cells', 'Human genes per cell (median)', 'Human UMIs per cell (median)')
#rownames(h.cells) = h.cells.names
#print(h.cells)
#
#print('mouse cells')
#m.cells = sapply(collections, function(a){
#    q = get_filename(a, 'mcells')
#    q = read.table(q, sep='\t', header=T, row.names=1)
#    q = q[q$NUM_GENES >= args$cutoff,]
#    num_cells = nrow(q)
#    num_genes = q[,1]
#    num_umis = q[,2]
#    return(list(num_cells=num_cells, num_genes=num_genes, num_umis=num_umis))
#})
#m.cells.names = c('Mouse cells', 'Mouse genes per cell (median)', 'Mouse UMIs per cell (median)')
#rownames(m.cells) = m.cells.names
#print(m.cells)

print('reads')
reads = sapply(collections, function(a){
    q = get_filename(a, 'reads')
    q = readLines(q)
    q = unlist(strsplit(q[grep('^all', q)], '\t'))
    total_reads = as.numeric(q[[2]])
    mapped_reads = as.numeric(q[[3]])
    hq_reads = as.numeric(q[[4]])
    return(c(total_reads, mapped_reads, hq_reads))
})
reads.names = c('Total reads', 'Mapped reads', 'Mapped HQ reads')
rownames(reads) = reads.names

print('genic')
genic = sapply(collections, function(a){
    q = get_filename(a, 'genic')
    q = readLines(q)
    q = unlist(strsplit(q[grep('^PF_', q)+1], '\t'))
    total = as.numeric(q[[1]])
    aligned = as.numeric(q[[2]])
    ribosomal = as.numeric(q[[3]])
    coding = as.numeric(q[[4]])
    utr = as.numeric(q[[5]])
    intronic = as.numeric(q[[6]])
    intergenic = as.numeric(q[[7]])
    return(c(total, aligned, ribosomal, coding, utr, intronic, intergenic))
})
genic.names = c('Total', 'Mapped', 'Ribosomal', 'Coding', 'UTR', 'Intronic', 'Intergenic')
rownames(genic) = genic.names

print('doublets')
doublets = sapply(collections, function(a){
    q = get_filename(a, 'doublet')
    q = read.table(q, sep='\t', header=T, row.names=1)
    return(as.character(q$organism))
})

# Per-collection

doublet_rates = sapply(doublets, function(a){
    q = table(a)
    if(length(q) == 3){
        h = q[['HUMAN']]
        m = q[['MOUSE']]
        b = q[['Mixed']]
        t = b/(2*(h/(h+m))*(m/(h+m)))
	return(t/(h+m+b))
    }else{return(NA)}
})
doublet_rates = t(data.frame(Doublets=doublet_rates))

# Cell statistics
#h.cells.tmp = h.cells
#m.cells.tmp = m.cells
#t.cells.tmp = matrix(NA, nrow=nrow(h.cells), ncol=ncol(h.cells))
#t.cells.tmp[1,] = as.numeric(h.cells.tmp[1,]) + as.numeric(m.cells.tmp[1,])
#for(i in 2:3){
#    h.cells.tmp[i,] = sapply(h.cells.tmp[i,], median)
#    m.cells.tmp[i,] = sapply(m.cells.tmp[i,], median)
#    for(j in 1:ncol(h.cells)){
#    	  u = unlist(h.cells[i,j])
#	  v = unlist(m.cells[i,j])
#	  t.cells.tmp[i,j] = median(c(u,v))
#    }
#}
#t.cells.names = c('Total cells', 'Total genes per cell (median)', 'Total UMIs per cell (median)')
#rownames(t.cells.tmp) = t.cells.names
#colnames(t.cells.tmp) = colnames(h.cells.tmp)
#print(t.cells.tmp)

# Calculate pcts
pct_reads = reads
for(i in 2:3){pct_reads[i,] = pct_reads[i,]/pct_reads[1,]}
pct_reads = pct_reads[2:nrow(pct_reads),]*100
pct_reads.names = gsub('$', ' (%)', rownames(pct_reads))
rownames(pct_reads) = pct_reads.names

pct_genic = genic
for(i in 3:7){pct_genic[i,] = pct_genic[i,]/pct_genic[2,]}
pct_genic = pct_genic[3:nrow(pct_genic),]*100
pct_genic.names = gsub('$', ' (%)', rownames(pct_genic))
rownames(pct_genic) = pct_genic.names

# Aggregate data
ncoll = ncol(reads)
#data = t(rbind(h.cells.tmp, m.cells.tmp, t.cells.tmp, reads, pct_reads, pct_genic, doublet_rates))
data = t(rbind(reads, pct_reads, pct_genic, doublet_rates))
data = cbind(rep(as.character(args$folder), nrow(data)), data)
colnames(data)[[1]] = 'Project'
write.table(data, file=paste(args$out, '.collections.dsq_qual.txt', sep=''), sep='\t', quote=F)

# Cell statistics
#h1 = sum(as.numeric(h.cells[1,]))/ncol(h.cells)
#h2 = median(unlist(c(h.cells[2,])))
#h3 = median(unlist(c(h.cells[3,])))
#h.cells.tmp = c(h1, h2, h3)
#
#m1 = sum(as.numeric(m.cells[1,]))/ncol(m.cells)
#m2 = median(unlist(c(m.cells[2,])))
#m3 = median(unlist(c(m.cells[3,])))
#m.cells.tmp = c(m1, m2, m3)
#
#t1 = h1 + m1
#t2 = median(c(unlist(c(h.cells[2,])), unlist(c(m.cells[2,]))))
#t3 = median(c(unlist(c(h.cells[3,])), unlist(c(m.cells[3,]))))
#t.cells.tmp = c(t1, t2, t3)

reads = rowSums(reads)
genic = rowSums(genic)
doublets = unlist(doublets)

doublet_rate = NA
q = table(doublets)
if(length(q) == 3){
    h = q[['HUMAN']]
    m = q[['MOUSE']]
    b = q[['Mixed']]
    t = b/(2*(h/(h+m))*(m/(h+m)))
    doublet_rate = t/(h+m+b)
}
doublet_rate = list(Doublet = doublet_rate)

pct_reads = reads
for(i in 2:3){pct_reads[[i]] = pct_reads[[i]]/pct_reads[[1]]}
pct_reads = pct_reads[2:length(pct_reads)]*100

pct_genic = genic
for(i in 3:7){pct_genic[[i]] = pct_genic[[i]]/pct_genic[[2]]}
pct_genic = pct_genic[3:length(pct_genic)]*100

#data = c(args$folder, ncoll, h.cells.tmp, m.cells.tmp, t.cells.tmp, reads, pct_reads, pct_genic, doublet_rate)
data = c(args$folder, ncoll, reads, pct_reads, pct_genic, doublet_rate)
data = data.frame(x=data)

#h.cells.names[[1]] = 'Human cells per collection'
#m.cells.names[[1]] = 'Mouse cells per collection'
#t.cells.names[[1]] = 'Total cells per collection'

#names = c('Project', 'Number of collections', h.cells.names, m.cells.names, t.cells.names, reads.names, pct_reads.names, pct_genic.names, 'Doublet rate')
names = c('Project', 'Number of collections', reads.names, pct_reads.names, pct_genic.names, 'Doublet rate')
write.table(data, file=paste(args$out, '.experiments.dsq_qual.txt', sep=''), sep='\t', col.names=names, quote=F)
