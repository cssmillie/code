require(optparse)


# ----------------------------------------------------------------------------------
# Combine columns of sparse or dense expression matrices with fast rbind/cbind joins
# ----------------------------------------------------------------------------------
#
# If "path" is given, join genes.tsv, barcodes.tsv, and matrix.mtx within that folder
# If "path" is null, then use the sample mapping file in "map"
#
# If "path" points to a directory, it assumes CellRanger directory format
# If "path" points to a file and ".mtx" suffix is found, then it reads with readMM
# Otherwise, it reads the file with ffread(path, row.names=1)
#
# Optionally:
# - "ming" (min genes per cell) and "minc" (min cells per gene) filters
# - "rename" removes all prefixes from cell barcodes
# - "sparse" saves output as a sparse matrix with writeMM
#
# Usage:
# Rscript merge_dges.r --path /path/to/cellranger/hg19 --out out.txt
# Rscript merge_dges.r --map /path/to/sample/map --out out.mtx --sparse


# ---------------
# Input arguments
# ---------------

option_list = list(make_option('--name', help='prefix added to cell names', default=NULL),
                   make_option('--path', help='sparse matrix path', default=NULL),
		   make_option('--pattern', help='regex pattern for cells to keep', default=NULL),
	           make_option('--map', help='input mapping file (1=prefix, 2=path, 3=pattern)', default=NULL),
		   make_option('--minc', help='cells per gene cutoff', type='integer', default=1),
                   make_option('--ming', help='genes per cell cutoff', type='integer', default=1),
                   make_option('--out', help='output file'),
		   make_option('--rename', help='use base cell names', action='store_true', default=FALSE),
		   make_option('--sparse', help='save as sparse matrix', action='store_true', default=FALSE)		   
                   )
args = parse_args(OptionParser(option_list=option_list))

require(data.table)
require(Matrix)
require(plyr)
require(tidyverse)


# ------------
# Map datasets
# ------------

if(!is.null(args$path)){
    if(is.null(args$pattern)){args$pattern='.*'}
    map = matrix(c(args$name, args$path, args$pattern), nrow=1, ncol=3)
} else {
    map = read.table(args$map, stringsAsFactors=F)
}
colnames(map) = c('name', 'path', 'pattern')
map$name[is.na(map$name)] = ''


# -------------
# Read matrices
# -------------


# CellRanger matrix filenames
mtx_filenames = function(path){
    
    if(dir.exists(path)){
        path = paste0(path, '/')
    } else {
        path = paste0(path, '.')
    }
    
    path = gsub('//', '/', path)
    path = gsub('\\.\\.', '\\.', path)
    path = gsub('/\\.', '/', path)
    
    list(genes=paste0(path, 'genes.tsv'),
         barcodes=paste0(path, 'barcodes.tsv'),
	 matrix=paste0(path, 'matrix.mtx'))
}


# Sparse matrix
read_mtx = function(path){
    
    # Get filenames
    fns = mtx_filenames(path)
    genes_fn = fns$genes
    barcodes_fn = fns$barcodes
    counts_fn = fns$matrix
    
    # Read matrix
    counts = readMM(counts_fn)
    genes = read.table(genes_fn, sep='\t')
    genes = genes[,ncol(genes)]
    barcodes = readLines(barcodes_fn)
    
    # Fix duplicate gene names
    counts = rowsum(as.matrix(counts), genes)
    colnames(counts) = barcodes
    
    as(counts, 'sparseMatrix')
}


# Text matrix
read_txt = function(path){
    counts = ffread(path, row.names=TRUE)
    as(as.matrix(counts), 'sparseMatrix')
}


# Read path (sparse or text)
read_path = function(path, prefix='', pattern='.*', ming=1, minc=1, rename=FALSE){
    
    # Read sparse or text matrix
    if(file.exists(path) & !dir.exists(path)){
        counts = read_txt(path)
    } else {
        counts = read_mtx(path)
    }
    
    # Fix formatting
    rownames(counts) = gsub('^mm10_', '', rownames(counts))
    rownames(counts) = gsub('^hg19_', '', rownames(counts))
    colnames(counts) = gsub('\\.1$', '', colnames(counts))
    if(rename == TRUE){colnames(counts) = gsub('.*\\.', '', colnames(counts))}
    colnames(counts) = paste(prefix, colnames(counts), sep='.')
    colnames(counts) = gsub('^\\.', '', colnames(counts))
    
    # Filter matrix
    genes.use = rowSums(counts > 0) >= minc
    cells.use = grepl(pattern, colnames(counts)) & (colSums(counts > 0) >= ming)
    counts = counts[genes.use, cells.use]
    
    # Print and return
    cat(paste0('\nRead ', path, ' [', nrow(counts), ' x ', ncol(counts), ']\n'))
    return(counts)
}

dges = list()
for(i in 1:nrow(map)){
    name = map[i,1]
    path = map[i,2]
    pattern = map[i,3]
    dges[[path]] = read_path(path, prefix=name, pattern=pattern, ming=args$ming, minc=args$minc, rename=args$rename)
}


# --------------
# Merge matrices
# --------------

cat('\nInitializing DGE\n')

# Make merged matrix
genes = unique(sort(unlist(sapply(dges, rownames))))
x = Matrix(0, nrow=length(genes), sparse=T)
rownames(x) = genes

# Fast merge using rbind/cbind
for(name in names(dges)){
    cat(paste('\nMerging', name, '\n'))
    
    # Select DGE
    dge = dges[[name]]

    # Fast align    
    new = Matrix(0, nrow=length(genes)-nrow(dge), ncol=ncol(dge), sparse=T)
    rownames(new) = setdiff(genes, rownames(dge))
    dge = rbind(dge, new)
    dge = dge[genes,]

    # Fast combine
    x = cbind(x, dge)
}

# Remove first column
x = x[,2:ncol(x)]


# ------------
# Write output
# ------------

# Free up memory
rm(dges)

# Filter data
genes.use = rowSums(x > 0) >= args$minc
cells.use = colSums(x > 0) >= args$ming
x = x[genes.use, cells.use]
cat(paste0('\nWriting DGE [', nrow(x), ' x ', ncol(x), ']\n'))

# Write to file
if(args$sparse == TRUE){
    
    # Get filenames
    fns = mtx_filenames(args$out)
    genes_fn = fns$genes
    barcodes_fn = fns$barcodes
    counts_fn = fns$matrix

    # Write sparse matrix
    writeLines(rownames(x), genes_fn)
    writeLines(colnames(x), barcodes_fn)
    writeMM(x, file=counts_fn)
    
} else {
    cat('\nWriting dense matrix\n')
    require(tibble)
    x = as.data.frame(as.matrix(x)) %>% rownames_to_column('GENE')
    write.table(x, file=args$out, quote=F, sep='\t', row.names=F)
}
