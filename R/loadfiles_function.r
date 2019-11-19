#' generate folder function
#'
#' This loads count and target file info (internal function only) 
#' @param count_file count matrix raw or normalized
#' @param target_file file with experimental information
#' @keywords loadsfiles 
#' @export

loadfiles <- function(count_file, target_file) {
    #Load in count data  and target information from csv'
    counts    <- read.table(count_file, header = TRUE, sep = "\t", row.names=1, as.is = TRUE)
    targets   <- read.table(target_file, header = TRUE, sep = ",", row.names = 1, as.is = TRUE)
    results   <- list("counts" = counts, "targets" = targets)
    return(results)
}
