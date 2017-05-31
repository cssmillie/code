library(ineq)

getNumReadsPerCellBarcode<-function (bamFile, cellTagCollapsed="XC", readMapQuality=10, organism=NULL, tempDir, reportDir, gapAnalysisDir="/broad/mccarroll/software/dropseq/prod" ) {
	
	o=paste(reportDir, "/", basename(bamFile), "_numReads_perCell_", cellTagCollapsed, "_mq_", readMapQuality, ".txt.gz", sep="")	
	command=paste(gapAnalysisDir, "/BAMTagHistogram I=", bamFile, " OUTPUT=", o, " TAG=", cellTagCollapsed, " READ_QUALITY=", readMapQuality, " FILTER_PCR_DUPLICATES=false", " TMP_DIR=", tempDir, sep="")
	system(command)
	return (o)
}

find_elbow <- function(a){

x = c(1:length(a))
a <- c(0,a)
x <- c(0,x)

x = x[1:5000]
a = a[1:5000]


dist <- c()

z = x[length(x)]
w = a[length(a)]

for (i in 1:length(x)){
    
    p = x[i]; q = a[i] 
    r = abs(z*q - w*p) / sqrt(z^2 + w^2)
    print(r)
    dist = c(dist,r)
    
}

print(dist[1:100])
return(x[which.max(dist)])


}

plotNumCellBarcodes<-function (cellBCCountsFile, cellBCCollapsedCountsFile, xlimit=NULL, selectedCellsFile=NULL) {
	
	a=read.table(cellBCCountsFile, header=F, stringsAsFactors=F)[,1:2]
	cell_barcodes=a[order(a$V1,decreasing=T),]
	
	b=read.table(cellBCCollapsedCountsFile, header=F, stringsAsFactors=F)[,1:2]
	cell_barcodes_collapsed=b[order(b$V1,decreasing=T),]
	
	cum=cumsum(cell_barcodes$V1)
	cum=cum/cum[length(cum)]

	cumB=cumsum(cell_barcodes_collapsed $V1)
	cumB=cumB/cumB[length(cumB)]
	
	
	numCells = find_elbow(cum)
	print(numCells)
        numCellsFound = NULL	
	
	
	pdf("figName_cumplot.pdf",width=10,height=10)
	if (is.null(xlimit)==F) {
		plot(1:length(cum), cum, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads", xlim=range(xlimit))
		points(1:length(cumB), cumB, type='l', col="green")
		points(c(numCells,numCells),c(0,1),type='l',col="black")
		text(x=numCells+10,y=0.9, labels=paste0("Ncells = ",numCells))
                text(x=numCells+10,y=0.75, labels=paste0("Mapping Rate = ", cum[numCells]/cum[length(cum)]*100, "%"))
                text(x=numCells+10, y=0.6, labels=paste0("Gini = ", round(Gini(cell_barcodes$V1[1:numCells]),2))) 
	} else {
		plot(1:length(cum), cum, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads")
		points(1:length(cumB), cumB, type='l', col="green")
		points(c(numCells,numCells),c(0,1),type='l',col="black")
		text(x=numCells+10,y=0.9, labels=paste0("Ncells = ",numCells))
                text(x=numCells+10,y=0.75, labels=paste0("Mapping Rate = ", round(cum[numCells]/cum[length(cum)]*100,2), "%"))
                text(x=numCells+10, y=0.6, labels=paste0("Gini = ", round(Gini(cell_barcodes$V1[1:numCells]),2))) 

	}
	
	if (!is.null(selectedCellsFile)) {
		z=read.table(selectedCellsFile, header=F, stringsAsFactors=F)
		abline(v=dim(z)[1], lwd=2)
		numCellsFound=dim(z)[1]
	}
	
	if (!is.null(numCellsFound)) {
		title(paste("Cumulative fraction of reads per cell barcode, cells found:", numCellsFound))
	} else {
		title("Cumulative fraction of reads per cell barcode")
	}
	
	legend('bottomright', legend=c("cell barcodes", "collapsed cell barcodes"), fill=c("blue", "green"))
	dev.off()
	
	#READS
	
	pdf("figName_NreadsHiToLo.pdf",width=10,height=10)
	barplot(cell_barcodes$V1,col="blue", xlab="Cells [descending]", ylab="# of reads per cell barcode", xlim=c(0,100), ylim=c(0, 1.05*max(cell_barcodes$V1)))
	title(paste0("Num of reads per cell (top 100), Total reads = ", sum(cell_barcodes$V1) ))
	dev.off()

	#Write output
	write(numCells, file="figName_numCells.txt")
	
	
}

cellBCCountsFile_Name = "fileName"
cellBCCollapsedCountsFile_Name = "fileName"
plotNumCellBarcodes(cellBCCountsFile_Name, cellBCCollapsedCountsFile_Name, xlimit=c(0,7000))

