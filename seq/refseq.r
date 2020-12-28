seurat()

# load metadata
info = as.data.table(read.table('~/aviv/db/refseq/bacteria.assembly_summary.txt', quote='', sep='\t', comment.char='', fill=T, header=T, stringsAsFactors=F))
meta = as.data.table(read.table('~/aviv/db/refseq/prokaryotes.meta.csv', sep=',', comment.char='', fill=T, header=T, stringsAsFactors=F))

# merge metadata
meta$Assembly = gsub('GCA', 'GCF', meta$Assembly)
meta = merge(info, meta, by.x='assembly_accession', by.y='Assembly', all=T)




# map strain to host
s2h = setNames(meta$Host, meta$Strain)
s2s = setNames(gsub(' .*', '', meta$X.Organism.Name), meta$Strain)
