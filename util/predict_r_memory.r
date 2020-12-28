# input arguments
args = commandArgs(trailingOnly=T)
rows = as.integer(args[[1]])
cols = as.integer(args[[2]])

# calculate memory
byte = rows*cols*8
mega = round(1.0*byte/(2^{20}),2)
giga = round(1.0*mega/1024,2)

print(paste('bytes:', byte))
print(paste('MB:', mega))
print(paste('GB:', giga))
