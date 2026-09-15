readExcel = function(filename, tibble=FALSE) {
    # https://stackoverflow.com/questions/12945687/read-all-worksheets-in-an-excel-workbook-into-an-r-list-with-data-frames
    require(readxl)
    sheets = readxl::excel_sheets(filename)
    x = lapply(sheets, function(.) {
        readxl::read_excel(filename, sheet=., na=c("", "NA", "NaN"))
    })
    if (!tibble) {x = lapply(x, as.data.frame)}
    names(x) = sheets
    return(x)
}
