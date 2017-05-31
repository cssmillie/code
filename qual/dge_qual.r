library(optparse)

option_list = list(make_option('--folder', help='folder name'),
	           make_option('--field', help='collection field', type='integer'),
		   make_option('--out', help='output prefix'))

args = parse_args(OptionParser(option_list=option_list))

# Get filenames
h.dge = paste(args$folder, '/dge/human.dge.txt.gz', sep='')
m.dge = paste(args$folder, '/dge/mouse.dge.txt.gz', sep='')

# Get dges
print('Loading Human DGE')
if(file.exists(h.dge)){
    h.dge = read.table(h.dge, sep='\t', header=T, row.names=1)
} else {
    h.dge=data.frame()
}

print('Loading Mouse DGE')
if(file.exists(m.dge)){
    m.dge = read.table(m.dge, sep='\t', header=T, row.names=1)
} else {
    m.dge=data.frame()
}

if(ncol(h.dge) == 0 & ncol(m.dge) == 0){stop('Check folder name')}

# Get collections
collections = c()

if(ncol(h.dge) > 1){
    h.collections = sapply(strsplit(colnames(h.dge), '\\.'), '[', args$field)
    print(paste('Found', length(unique(h.collections)), 'human collections'))
    collections = c(collections, h.collections)
} else {h.collections = c()}

if(ncol(m.dge) > 1){
    m.collections = sapply(strsplit(colnames(m.dge), '\\.'), '[', args$field)
    print(paste('Found', length(unique(m.collections)), 'mouse collections'))
    collections = c(collections, m.collections)
} else {m.collections = c()}

collections = unique(collections)
print(collections)

# Calculate quality statistics
quals1 = c() # per-collection
quals2 = c() # per-experiment
for(cutoff in c(500, 750, 1000)){
    
    qual1 = matrix(0, nrow=length(collections), ncol=9)
    qual2 = rep(0, 12)
    rownames(qual1) = collections

    if(ncol(h.dge) > 1){
        h.num_genes = colSums(h.dge > 0)
	h.num_reads = colSums(h.dge)
	h.cells = tapply(h.num_genes >= cutoff, h.collections, sum)
	h.genes = tapply(h.num_genes, h.collections, median)
	h.reads = tapply(h.num_reads, h.collections, median)
	qual1[collections, 1:3] = cbind(h.cells[collections], h.genes[collections], h.reads[collections])
	qual2[1:4] = c(sum(h.num_genes >= cutoff), sum(h.num_genes >= cutoff)/length(unique(h.collections)), median(h.num_genes), median(h.num_reads))
    } else {h.num_genes = c(); h.num_reads = c()}

    if(ncol(m.dge) > 1){
        m.num_genes = colSums(m.dge > 0)
	m.num_reads = colSums(m.dge)
	m.cells = tapply(m.num_genes >= cutoff, m.collections, sum)
        m.genes = tapply(m.num_genes, m.collections, sum)
        m.reads = tapply(m.num_reads, m.collections, median)
	qual1[collections,4:6] = cbind(m.cells[collections], m.genes[collections], m.reads[collections])
	qual2[5:8] = c(sum(m.num_genes >= cutoff), sum(m.num_genes >= cutoff)/length(unique(m.collections)), median(m.num_genes), median(m.num_reads))
    } else {m.num_genes = c(); m.num_reads = c()}

    t.num_genes = c(h.num_genes, m.num_genes)
    t.num_reads = c(h.num_reads, m.num_reads)
    t.collections = c(h.collections, m.collections)
    
    t.cells = tapply(t.num_genes >= cutoff, t.collections, sum)
    t.genes = tapply(t.num_genes, t.collections, median)
    t.reads = tapply(t.num_reads, t.collections, median)
    qual1[collections, 7:9] = cbind(t.cells[collections], t.genes[collections], t.reads[collections])
    qual2[9:12] = c(sum(t.num_genes >= cutoff), sum(t.num_genes >= cutoff)/length(collections), median(t.num_genes), median(t.num_reads))

    quals1 = cbind(quals1, qual1)
    quals2 = c(quals2, qual2)
}

q1 = expand.grid(c('cells', 'median genes per cell', 'median reads per cell'), c('Human', 'Mouse', 'Total'), c('(cutoff=500)', '(cutoff=750)', '(cutoff=1000)'))[,c(2,1,3)]
q1 = c('Project', paste(as.character(q1[,1]), as.character(q1[,2]), as.character(q1[,3])))

q2 = expand.grid(c('cells', 'cells per collection', 'median genes per cell', 'median reads per cell'), c('Human', 'Mouse', 'Total'), c('(cutoff=500)', '(cutoff=750)', '(cutoff=1000)'))[,c(2,1,3)]
q2 = c('Project', paste(as.character(q2[,1]), as.character(q2[,2]), as.character(q2[,3])))

quals1 = cbind(rep(as.character(args$folder), nrow(quals1)), quals1)
colnames(quals1) = q1

quals2 = c(args$folder, quals2)
quals2 = t(data.frame(cbind(q2, quals2), row.names=1))

write.table(quals1, file=paste(args$out, '.collections.dge_qual.txt', sep=''), sep='\t', quote=F)
write.table(quals2, file=paste(args$out, '.experiments.dge_qual.txt', sep=''), sep='\t', quote=F, row.names=F)
