options(java.parameters = "-Xmx8g")
library(xlsx)

write_excel = function(data, out, row.names=FALSE){

    # Write list of data frames to Excel spreadsheet
    # data = input data (named list)
    # out = output file (.xls)
    
    for(name in names(data)){
        write.xlsx(data[[name]], file=out, sheetName=name, append=TRUE, row.names=FALSE)
    }
}
