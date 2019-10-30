#' s3 differential gene analysis 
#'
#' This function runs differential gene analysis (linear modeling on normalized counts)
#' 
#' 
#' 
#' @keywords gene differential analysis
#' @export
#' @import edgeR
#' @import limma
#' @examples
#' 

s3_DE_analysis <- function() { 
    design     <- readRDS('./s1_norm_raw_counts_results/1.designobject.rds')
    eset_voom  <- readRDS('./s1_norm_raw_counts_results/1.voomobject.rds')


    
}