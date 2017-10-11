library(optparse)
library(tibble)
source('~/code/single_cell/seurat.r')

if(interactive()){
args = list(name='./outs/Naive_DC', dge='./outs/Naive_DC.dge.txt.gz', ming=500, sigdb='')
} else {

option_list = list(make_option('--name', help='filename prefix'),
	           make_option('--dge', help='dge file'),
		   make_option('--ming', help='minimum number of genes', default=500, type='integer'),
		   make_option('--out', help='write metadata and tpm', default=FALSE, action='store_true'),
		   make_option('--sigdb', help='gene signature database'))
args = parse_args(OptionParser(option_list=option_list))
}

seur = run_seurat(name=args$name, dge=args$dge, ming=args$ming, verbose=T, write_out=F, ncores=4)

format_portal = function(x){

    # Fix factor columns
    i = sapply(x, is.factor)
    x[,i] = sapply(x[,i], as.character)
    
    # Add data types
    h = c('TYPE', ifelse(sapply(x, is.numeric), 'numeric', 'group'))
    
    # Combine rows
    x = data.frame(x, stringsAsFactors=F) %>% rownames_to_column('NAME')
    h = data.frame(t(h), stringsAsFactors=F)
    colnames(h) = colnames(x)
    x = rbind(h, x)
    
    return(x)
}

# Write cluster information
clusters = seur@tsne.rot
colnames(clusters) = c('X', 'Y')
clusters = format_portal(clusters)
write.table(clusters, file=paste0(args$name, '.clusters.txt'), quote=F, sep='\t', row.names=F)

if(args$out){

    # Write metadata
    metadata = subset(seur@data.info, select=-c(num_pcs))
    metadata = format_portal(metadata)
    write.table(metadata, file=paste0(args$name, '.metadata.txt'), quote=F, sep='\t', row.names=F)
    
    # Write log2(TPM+1)
    tpm_fn = paste0(args$name, '.log2tpm.txt')
    tpm = seur@data %>% rownames_to_column('GENE')
    write.table(tpm, file=tpm_fn, sep='\t', quote=F, row.names=F)
}
