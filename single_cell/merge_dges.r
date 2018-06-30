require(optparse)


# ----------------------------------------------------------------------------------
# Combine columns of sparse or dense expression matrices with fast rbind/cbind joins
# ----------------------------------------------------------------------------------
#
# Must specify sample mapping file ("map") or ("name", "path", and "pattern")
# If "path" points to a file, reads as dense matrix with ffread
# Otherwise, reads as sparse matrix with read_mtx(prefix = path)
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
		   make_option('--pattern', help='regex pattern for cells to keep', default='.*'),
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
source('~/code/util/mtx.r')


# ------------
# Map datasets
# ------------

if(!is.null(args$path)){
    map = matrix(c(args$name, args$path, args$pattern), nrow=1, ncol=3)
} else {
    map = read.table(args$map, stringsAsFactors=F)
}
colnames(map) = c('name', 'path', 'pattern')
map$name[is.na(map$name)] = ''


# -------------
# Read matrices
# -------------


# Read DGE (sparse or text)
read_dge = function(path, prefix='', pattern='.*', ming=1, minc=1, rename=FALSE){
    
    # Read sparse or text matrix
    if(file.exists(path) & !dir.exists(path)){
        counts = ffread(path, row.names=TRUE)
	counts = as(as.matrix(counts), 'sparseMatrix')
    } else {
        counts = read_mtx(prefix=path, data='matrix.mtx', rows='genes.tsv', cols='barcodes.tsv', fix_duplicates=TRUE)
    }
    
    # Fix formatting
    rownames(counts) = gsub('^mm10_', '', rownames(counts))
    rownames(counts) = gsub('^hg19_', '', rownames(counts))
    colnames(counts) = gsub('\\.1$', '', colnames(counts))
    colnames(counts) = gsub('-1$', '', colnames(counts))
    if(rename == TRUE){colnames(counts) = gsub('.*\\.', '', colnames(counts))}
    colnames(counts) = paste(prefix, colnames(counts), sep='.')
    colnames(counts) = gsub('^\\.', '', colnames(counts))
    
    # Filter matrix
    print(dim(counts))
    genes.use = rowSums(counts > 0) >= minc
    cells.use = grepl(pattern, colnames(counts)) & (colSums(counts > 0) >= ming)
    counts = counts[genes.use, cells.use]
    print(dim(counts))
    
    # Print and return
    cat(paste0('\nRead ', path, ' [', nrow(counts), ' x ', ncol(counts), ']\n'))
    return(counts)
}

dges = list()
for(i in 1:nrow(map)){
    name = map[i,1]
    path = map[i,2]
    pattern = map[i,3]
    dge = read_dge(path, prefix=name, pattern=pattern, ming=args$ming, minc=args$minc, rename=args$rename)
    if(!is.null(dim(dge))){dges = c(dges, dge)}
}


# --------------
# Merge matrices
# --------------

cat('\nInitializing DGE\n')
x = sparse_cbind(dges)

# Check for replicate cell barcodes
if(max(table(colnames(x))) > 1){
    cat(paste('\nMerging', ncol(x) - length(unique(colnames(x))), 'replicate cell barcodes'))
    print(mean(colSums(x[,grep('ColFr0_Mye_3', colnames(x))])))
    x = as(t(rowsum(t(as.matrix(x)), colnames(x))), 'sparseMatrix')
    print(mean(colSums(x[,grep('ColFr0_Mye_3', colnames(x))])))
}


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
    write_mtx(x, prefix=args$out, data='matrix.mtx', rows='genes.tsv', cols='barcodes.tsv')
} else {
    cat('\nWriting dense matrix\n')
    require(tibble)
    x = as.data.frame(as.matrix(x)) %>% rownames_to_column('GENE')
    write.table(x, file=args$out, quote=F, sep='\t', row.names=F)
}
