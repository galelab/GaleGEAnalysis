#' generate folder function
#'
#' This loads count and target file info (internal function only) 
#' @param count_file count matrix raw or normalized
#' @param target_file file with experimental information
#' @keywords loadsfiles
#' @export

loadfiles <- function(count_file, target_file) {
    #Load in count data  and target information from csv'
    if (typeof(count_file) == "character") {
        counts    <- read.table(count_file, header = TRUE,
                            sep = "\t", row.names = 1,
                            as.is = TRUE, check.names = FALSE)
        counts    <- counts[, order(names(counts))]
    } else if ((isFALSE(count_file)) || (is.null(count_file))) {
        counts <- FALSE
    }

    if (typeof(target_file) == "character") {
        targets   <- read.table(target_file, header = TRUE,
                                sep = ",", row.names = 1,
                                as.is = TRUE, check.names = FALSE)
        targets   <- targets[order(rownames(targets)), ]
    } else if (((isFALSE(target_file)) || is.null(count_file))) {
        targets <- FALSE
    }

    results   <- list("counts" = counts, "targets" = targets)
    return(results)
}
