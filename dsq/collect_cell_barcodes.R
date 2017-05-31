

collect_cell_barcodes <- function(filename,Ncells){ 

data = read.table(file=filename)
barcodes = data$V2[1:Ncells]

filename = gsub("reads.txt","barcodes_use.txt",filename)
write.table(as.data.frame(barcodes), file=filename, col.names=FALSE, row.names=FALSE,quote=FALSE)
}

collect_cell_barcodes("../bam_reads/filename_input", Ncells_input)
