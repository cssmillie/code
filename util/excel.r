options(java.parameters = "-Xmx8g")
library(xlsx)

write_excel = function(data, out, max_rows=500, row.names=FALSE){

    # Write list of data frames to Excel spreadsheet
    # data = input data (named list)
    # out = output file (.xls)
    
    for(name in names(data)){print(paste('Writing', name))
        di = data[[name]]
	di = di[1:min(max_rows, nrow(di))]
        write.xlsx(di, file=out, sheetName=name, append=TRUE, row.names=row.names)
    }
}
