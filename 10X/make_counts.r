library(optparse)
library(methods)
library(tidyverse)
library(tibble)

# Get input arguments
option_list = list(
	    make_option('--folder', help='folder name'),
	    make_option('--prefix', help='prepend prefix to barcode', default=''),
	    make_option('--out', help='output file prefix')
	    )
args = parse_args(OptionParser(option_list=option_list))

# Get filenames
genes_fn = file.path(args$folder, 'genes.tsv')
barcodes_fn = file.path(args$folder, 'barcodes.tsv')
counts_fn = file.path(args$folder, 'matrix.mtx')

# Make sure files exist
if(!file.exists(genes_fn)){stop(paste(genes_fn, 'not found'))}
if(!file.exists(barcodes_fn)){stop(paste(barcodes_fn, 'not found'))}
if(!file.exists(counts_fn)){stop(paste(counts_fn, 'not found'))}

# Load data
counts = as.matrix(readMM(counts_fn))
genes = read.table(genes_fn, sep='\t')[,2]
barcodes = readLines(barcodes_fn)

# Add prefix to cell barcodes
if(args$prefix != ''){
    colnames(counts) = paste0(prefix, '.', barcodes)
} else {
    colnames(counts) = barcodes
}

# Fix duplicate gene names
counts = data.frame(aggregate(counts, list(genes), sum), row.names=1)

# Fix cell barcodes
colnames(counts) = gsub('\\.1$', '', colnames(counts))

# Write output
counts = counts %>% rownames_to_column('GENE')
write.table(counts, file=args$out, sep='\t', quote=F, row.names=F)
