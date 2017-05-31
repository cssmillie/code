library(Matrix)

args = commandArgs(trailingOnly=TRUE)

folder = args[[1]]
out = args[[2]]

genes = file.path(folder, 'genes.tsv')
barcodes = file.path(folder, 'barcodes.tsv')
counts = file.path(folder, 'matrix.mtx')

counts = as.matrix(readMM(counts))
genes = read.table(genes, sep='\t')[,2]
barcodes = readLines(barcodes)

rownames(counts) = genes
colnames(counts) = barcodes

write.table(counts, file=out, sep='\t', quote=F, col.names=NA)
