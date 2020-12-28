# notes:
# 1. run this script with "conda deactivate" in the beginning
#    example:
#    for x in {1..10}; do echo "conda deactivate; Rscript ~/code/tree/treedist.parallel.r --trees test.trees.txt -i $x -n 10 --out test.trees.qdist --ntip 100 --norm"; done | ssub3
# 2. after running in parallel, merge results by re-running command with --merge flag
#    example:
#    Rscript ~/code/tree/treedist.parallel.r --trees test.trees.txt -i $x -n 10 --out test.trees.qdist --ntip 100 --norm --merge test.trees.qdist.rds


library(argparse)

# input arguments
parser = ArgumentParser()
parser$add_argument('--trees', help='list of tree files')
parser$add_argument('-i', help='index', type='integer')
parser$add_argument('-n', help='total', type='integer')
parser$add_argument('--out', help='output prefix (.rds)')
parser$add_argument('--regex', help='leaf regex', default='\\..*')
parser$add_argument('--type', help='type of distance', default='quartet')
parser$add_argument('--norm', help='normalize?', default=FALSE, action='store_true')
parser$add_argument('--ntip', help='# tips to subsample', default=0, type='integer')
parser$add_argument('--merge', help='merge results (.rds)', default=FALSE, action='store_true')
parser$add_argument('--boot', help='collapse nodes with low bootstrap support', default=0, type='integer')
args = parser$parse_args()

if(args$merge == FALSE){

library(ape)
library(phangorn)
library(phytools)
source('~/code/tree/treedist.r')

# load trees
print('loading trees')
trees = sapply(readLines(args$trees), read.tree, simplify=F)
if(args$boot > 0){
    library(ggtree)
    trees = sapply(trees, function(a) as.polytomy(a, feature='node.label', fun=function(x) as.numeric(x) < args$boot), simplify=F)
}

# distance matrix
print('initializing data')
D = matrix(NA, nrow=length(trees), ncol=length(trees))
diag(D) = 0
rownames(D) = colnames(D) = names(trees)
N = D

# fix ntip
if(args$type == 'quartet'){
    if(args$ntip == 0 | args$ntip > 450){
        args$ntip = 450
	print('warning: quartet method requires ntip < 477. setting to 450')
    }
}

# get i,j indices
ind = which(upper.tri(D), arr.ind=T)
ind = ind[split(1:nrow(ind), ceiling(1:nrow(ind)/(nrow(ind)/args$n)))[[args$i]],]

# calculate distances
print('calculating distances')
for(k in 1:nrow(ind)){
    i = ind[k,1]
    j = ind[k,2]
    res = tdist_pair(trees[[i]], trees[[j]], label_fxn=function(a) gsub(args$regex, '', a), type=args$type, normalize=args$norm, ntip=args$ntip, ret_ntip=TRUE)
    D[i,j] = res$dist
    N[i,j] = res$ntip
}

# save results
print('saving results')
saveRDS(D, file=paste0(args$out, '.', args$i, '.dist.rds'))
saveRDS(N, file=paste0(args$out, '.', args$i, '.ntip.rds'))

} else {

dfns = paste0(args$out, '.', 1:args$n, '.dist.rds')
nfns = paste0(args$out, '.', 1:args$n, '.ntip.rds')

D = N = NULL

print('Loading results')
for(fn in dfns){
    xi = readRDS(fn)
    xi = replace(xi, is.na(xi), 0)
    if(is.null(D)){D = xi} else {D = D + xi}
}
for(fn in nfns){
    xi = readRDS(fn)
    xi = replace(xi, is.na(xi), 0)
    if(is.null(N)){N = xi} else {N = N + xi}    
}

#dres = sapply(dfns, readRDS, simplify=F)
#nres = sapply(nfns, readRDS, simplify=F)
print('Merging results')
#sum.na = function(x,y) {replace(x, is.na(x), 0) + replace(y, is.na(y), 0)}
#D = Reduce(sum.na, dres)
#N = Reduce(sum.na, nres)
print('Saving results')
saveRDS(D, file=paste0(args$out, '.dist.rds'))
saveRDS(N, file=paste0(args$out, '.ntip.rds'))

}
