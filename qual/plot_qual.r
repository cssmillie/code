library(ggplot2)
library(optparse)
library(plyr)
library(proxy)

multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    require(grid)
    plots = c(list(...), plotlist)
    numPlots = length(plots)
    if (is.null(layout)) {
        layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots==1) {
        print(plots[[1]])
    } else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        for (i in 1:numPlots) {
            matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

if(interactive()){
    args = list()
    args$map = 'test.in'
    args$out = 'test'
} else{
# Read input arguments
option_list = list(make_option('--map', help='map (1=id, 2=folder, 3=pattern'),
                   make_option('--out', help='output prefix')
                   )
args = parse_args(OptionParser(option_list=option_list))
}

# Load input data
data = read.table(args$map, stringsAsFactors=F)
data = cbind(data, paste('/home/unix/csmillie/aviv/data/', data[,2], '/raw', sep=''))
colnames(data) = c('id', 'project_id', 'pattern', 'path')

# Get files
groups = data[,1]
q_fns = list() # read quality metrics
h_fns = list() # human dge summaries
m_fns = list() # mouse dge summaries
c_fns = list() # cell types
for(i in 1:nrow(data)){
    group_id = data[i,1]
    
    path = paste(data[i,4], '/QC_files', sep='')
    pattern = paste(data[i,3], '.*ReadQualityMetrics.txt', sep='')
    files = list.files(path=path, pattern=pattern, full.names=TRUE)
    q_fns[[group_id]] = files
    file_ids = gsub('.*/', '', gsub('_star.*', '', files))
    
    path = paste(data[i,4], '/UMI_DGE', sep='')
    pattern = paste(data[i,3], '.*HUMAN.*summary.txt', sep='')
    files = list.files(path=path, pattern=pattern, full.names=TRUE)
    new_ids = gsub('.*/', '', gsub('_star.*', '', files))
    h_fns[[group_id]] = files
    if(sum(file_ids != new_ids) > 0){stop('files not aligned')}
    
    path = paste(data[i,4], '/UMI_DGE', sep='')
    pattern = paste(data[i,3], '.*MOUSE.*summary.txt', sep='')
    files = list.files(path=path, pattern=pattern, full.names=TRUE)
    new_ids = gsub('.*/', '', gsub('_star.*', '', files))
    m_fns[[group_id]] = files
    if(sum(file_ids != new_ids) > 0){stop('files not aligned')}
    
    path = paste(data[i,4], '/QC_files', sep='')
    pattern = paste(data[i,3], '.*categorized_cellTypes.txt', sep='')
    files = list.files(path=path, pattern=pattern, full.names=TRUE)
    new_ids = gsub('.*/', '', gsub('_star.*', '', files))
    c_fns[[group_id]] = files
    if(sum(file_ids != new_ids) > 0){stop('files not aligned')}
}

# Get read counts (rows = group ids, columns: 1=total, 2=mapped, 3=hq)
reads = sapply(q_fns, function(a){
    sapply(a, function(b){
        q = readLines(b)
        q = unlist(strsplit(q[grep('^all', q)], '\t'))
        q = as.integer(q[c(2,3,5)])
        return(q)
    })
})
# Calculate number of collections (names=group ids)
collections = sapply(reads, ncol)
reads = t(sapply(reads, rowSums))

# Get gene and transcript counts
# columns = group ids
# row 1 = gene counts
# row 2 = transcript counts
h.counts = lapply(h_fns, function(a){
    lapply(a, function(b){
        q = read.table(b, sep='\t', header=T, row.names=1)
        q = q[q[,1] >= 500,]
        return(q)
    })
})
h.counts = sapply(h.counts, function(a){do.call('rbind', a)})

m.counts = lapply(m_fns, function(a){
    lapply(a, function(b){
        q = read.table(b, sep='\t', header=T, row.names=1)
        q = q[q[,1] >= 500,]
        return(q)
    })
})
m.counts = sapply(m.counts, function(a){do.call('rbind', a)})

# Get categorized cell types
cell_types = sapply(c_fns, function(a){
    lapply(a, function(b){
        q = read.table(b, sep='\t', header=T, row.names=1)
        return(q$organism)
    })
})

# Plot variables
# num_genes = number of genes per cell
# num_cells = number of cells per collection
# num_cells.norm = number of cells per collection per read
# num_reads.norm = number of reads per collection per read

# human plot variables
h.group_ids = c()
h.num_genes = c()
h.num_cells = c()
h.num_cells.norm = c()
h.num_reads.norm = c()

# mouse plot variables
m.group_ids = c()
m.num_genes = c()
m.num_cells = c()
m.num_cells.norm = c()
m.num_reads.norm = c()

# total plot variables
t.group_ids = c()
t.num_genes = c()
t.num_cells = c()
t.num_cells.norm = c()
t.num_reads.norm = c()

for(i in 1:nrow(data)){

    group = data[i,1]
    ncoll = collections[[group]]
    
    h.genes = h.counts[1,group][[1]]
    h.cells = rep(1, length(h.genes))
    h.reads = h.counts[2,group][[1]]

    m.genes = m.counts[1,group][[1]]
    m.cells = rep(1, length(m.genes))
    m.reads = m.counts[2,group][[1]]

    t.genes = c(h.genes, m.genes)
    t.reads = c(h.reads, m.reads)
    t.cells = rep(1, length(t.genes))

    h.o = order(h.genes, decreasing=T)
    h.genes = h.genes[h.o]
    h.reads = h.reads[h.o]

    m.o = order(m.genes, decreasing=T)
    m.genes = m.genes[m.o]
    m.reads = m.reads[m.o]

    t.o = order(t.genes, decreasing=T)
    t.genes = t.genes[t.o]
    t.reads = t.reads[t.o]

    total_reads = reads[group,1]
    mapped_reads = reads[group,2]
    hq_reads = reads[group,3]
    
    h.group_ids = c(h.group_ids, rep(group, length(h.genes)))
    m.group_ids = c(m.group_ids, rep(group, length(m.genes)))
    t.group_ids = c(t.group_ids, rep(group, length(t.genes)))
    
    h.num_genes = c(h.num_genes, h.genes)
    m.num_genes = c(m.num_genes, m.genes)
    t.num_genes = c(t.num_genes, t.genes)

    h.num_cells = c(h.num_cells, cumsum(h.cells)/ncoll)
    m.num_cells = c(m.num_cells, cumsum(m.cells)/ncoll)
    t.num_cells = c(t.num_cells, cumsum(t.cells)/ncoll)
    
    h.num_cells.norm = c(h.num_cells.norm, cumsum(h.cells)/reads[group,1]/ncoll)
    m.num_cells.norm = c(m.num_cells.norm, cumsum(m.cells)/reads[group,1]/ncoll)
    t.num_cells.norm = c(t.num_cells.norm, cumsum(t.cells)/reads[group,1]/ncoll)

    h.num_reads.norm = c(h.num_reads.norm, cumsum(h.reads)/reads[group,1]/ncoll)
    m.num_reads.norm = c(m.num_reads.norm, cumsum(m.reads)/reads[group,1]/ncoll)
    t.num_reads.norm = c(t.num_reads.norm, cumsum(t.reads)/reads[group,1]/ncoll)

}

# Convert to dataframe
h.data = data.frame(group_ids=h.group_ids, num_genes=h.num_genes, num_cells=h.num_cells, num_cells.norm=h.num_cells.norm, num_reads.norm=h.num_reads.norm)
m.data = data.frame(group_ids=m.group_ids, num_genes=m.num_genes, num_cells=m.num_cells, num_cells.norm=m.num_cells.norm, num_reads.norm=m.num_reads.norm)
t.data = data.frame(group_ids=t.group_ids, num_genes=t.num_genes, num_cells=t.num_cells, num_cells.norm=t.num_cells.norm, num_reads.norm=t.num_reads.norm)

# Plot cells per collection
p1 = ggplot(h.data) + geom_line(aes(x=num_genes, y=num_cells, colour=group_ids)) + xlim(c(500,2000)) + xlab('Genes per cell (cutoff)') + ylab('Cells per collection') + theme_minimal() + ggtitle('Human')
p2 = ggplot(m.data) + geom_line(aes(x=num_genes, y=num_cells, colour=group_ids)) + xlim(c(500,2000)) + xlab('Genes per cell (cutoff)') + ylab('Cells per collection') + theme_minimal() + ggtitle('Mouse')
p3 = ggplot(t.data) + geom_line(aes(x=num_genes, y=num_cells, colour=group_ids)) + xlim(c(500,2000)) + xlab('Genes per cell (cutoff)') + ylab('Cells per collection') + theme_minimal() + ggtitle('Total')
pdf(paste(args$out, '.cells.pdf', sep=''), width=14, height=9)
multiplot(p1, p2, p3, cols=2)
dev.off()

# Plot cells per collection per read
p1 = ggplot(h.data) + geom_line(aes(x=num_genes, y=num_cells.norm, colour=group_ids)) + xlim(c(500,2000)) + xlab('Genes per cell (cutoff)') + ylab('Cells per collection per read') + theme_minimal() + ggtitle('Human')
p2 = ggplot(m.data) + geom_line(aes(x=num_genes, y=num_cells.norm, colour=group_ids)) + xlim(c(500,2000)) + xlab('Genes per cell (cutoff)') + ylab('Cells per collection per read') + theme_minimal() + ggtitle('Mouse')
p3 = ggplot(t.data) + geom_line(aes(x=num_genes, y=num_cells.norm, colour=group_ids)) + xlim(c(500,2000)) + xlab('Genes per cell (cutoff)') + ylab('Cells per collection per read') + theme_minimal() + ggtitle('Total')
pdf(paste(args$out, '.cells.norm.pdf', sep=''), width=14, height=9)
multiplot(p1, p2, p3, cols=2)
dev.off()

# Plot transcripts per collection per read
p1 = ggplot(h.data) + geom_line(aes(x=num_genes, y=num_reads.norm, colour=group_ids)) + xlim(c(500,2000)) + xlab('Genes per cell (cutoff)') + ylab('Transcripts per collection per read') + theme_minimal() + ggtitle('Human')
p2 = ggplot(m.data) + geom_line(aes(x=num_genes, y=num_reads.norm, colour=group_ids)) + xlim(c(500,2000)) + xlab('Genes per cell (cutoff)') + ylab('Transcripts per collection per read') + theme_minimal() + ggtitle('Mouse')
p3 = ggplot(t.data) + geom_line(aes(x=num_genes, y=num_reads.norm, colour=group_ids)) + xlim(c(500,2000)) + xlab('Genes per cell (cutoff)') + ylab('Transcripts per collection per read') + theme_minimal() + ggtitle('Total')
pdf(paste(args$out, '.reads.norm.pdf', sep=''), width=14, height=9)
multiplot(p1, p2, p3, cols=2)
dev.off()

# Plot total, mapped, and hq reads
pdf(paste(args$out, '.read_stats.pdf', sep=''))
par(mfrow=c(2,1), mar=c(8,5,1,1))
barplot(t(reads/collections), las=2, beside=T, col=rainbow(3), ylab='Reads per collection', cex.names=.6, cex.axis=.6, legend=TRUE)
legend('topright', legend=c('Total', 'Mapped', 'HQ'), fill=rainbow(3))
barplot(t(reads[,2:3]/reads[,1]), las=2, beside=T, col=rainbow(2), ylim=c(0,1), ylab='Fraction of total reads', cex.names=.6, cex.axis=.6, legend=TRUE)
legend('topright', legend=c('Mapped', 'HQ'), fill=rainbow(2))
dev.off()

# Estimate doublet rates
doublet_rates = list()
for(group in names(cell_types)){
    counts = table(unlist(cell_types[[group]]))
    if(length(counts) == 3){
        h = counts[['HUMAN']] # human cells
        m = counts[['MOUSE']] # mouse cells
        b = counts[['Mixed']] # doublets
        # observed doublets = 2*Pr(h)*Pr(m)*total_doublets
        t = b/(2*(h/(h+m))*(m/(h+m)))
        doublet_rates[[group]] = t/(h+m+b)
    }
}
pdf(paste(args$out, '.doublet_rates.pdf', sep=''))
barplot(t(as.matrix(doublet_rates)), col='orange', ylab='Estimated doublet rate')
dev.off()
