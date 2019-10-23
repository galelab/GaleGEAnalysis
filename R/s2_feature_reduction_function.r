#' s2 feature reduction 
#'
#' This function does feature reduction and visualization on normalized voom counts
#' @param countfile raw counts table (generally output by htseq).
#' @param targetfile target file.
#' @param  figres resolution at which to output figures (default is 300).
#' @keywords gene expression feature reduction
#' @export
#' @import limma
#' @import edgeR
#' @examples
#' s2_feature_reduction(count_file.txt, target_file.csv, vizualize_data=TRUE, FilterGenesWithCounts=100)

s2_feature_reduction <- function(countfile, targetfile, figres=100) { 
    V.CPM <- readRDS('./s1_norm_raw_counts_results/1.voomobject.rds')
    files <- loadfiles(count_file=countfile, target_file=targetfile)
    results_path <- generate_folder('s2_feature_reduction_results')
    vizualize_feature_reduction_data(V.CPM, files$target$treatment,results_path, figres)

}

vizualize_feature_reduction_data <- function(data, labels, results_path, figres) { 
    ###MDS (multidimensional scaling) uses log fold changes between genes as distances
    pdf(NULL)
    MDS <- plotMDS(data, labels = labels,  cex= .8)
    minx<-min(MDS$x)
    maxx<-max(MDS$x)
    miny<-min(MDS$y)
    maxy<-max(MDS$y)
    png(file.path(results_path,  '2.mds_vnorm_matrix.png'), res = figres)
    plot(MDS$x, MDS$y, cex=1, xlim=c(minx-1, maxx+1), ylim=c(miny-1, maxy+1), xlab=paste0(MDS$axislabel,' 1'), ylab=paste0(MDS$axislabel,' 2'), frame = FALSE)
    text(MDS$x, MDS$y, labels, cex=0.6, pos=4)
    dev.off()

}

