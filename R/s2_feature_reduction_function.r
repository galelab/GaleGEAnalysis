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
#' @import Biobase
#' @import ExpressionNormalizationWorkflow
#' @import stringr
#' @import stats
#' @import factoextra
#' @examples
#' s2_feature_reduction(count_file.txt, target_file.csv, vizualize_data=TRUE, FilterGenesWithCounts=100)

s2_feature_reduction <- function(countfile, targetfile, figres=100) { 
    # V.CPM <- readRDS('./s1_norm_raw_counts_results/1.voomobject.rds')
    pdf(NULL)
    files <- loadfiles(count_file=countfile, target_file=targetfile)
    results_path <- generate_folder('s2_feature_reduction_results')
    vizualize_feature_reduction_data(files$counts, files$target$treatment,results_path, figres)
    pvca(files$counts, files$target, results_path)
    pca(files$counts, files$target, results_path)

}

vizualize_feature_reduction_data <- function(data, labels, results_path, figres=100) { 
    ###MDS (multidimensional scaling) uses log fold changes between genes as distances
    MDS <- plotMDS(data, gene.selection='pairwise', cex= .8)
    minx<-min(MDS$x)
    maxx<-max(MDS$x)
    miny<-min(MDS$y)
    maxy<-max(MDS$y)
    png(file.path(results_path,  '2.mds_vnorm_matrix.png'), res = figres)
    plot(MDS$x, MDS$y, cex=1, xlim=c(minx-1, maxx+1), ylim=c(miny-1, maxy+1), xlab=paste0(MDS$axislabel,' 1'), ylab=paste0(MDS$axislabel,' 2'), frame = FALSE)
    text(MDS$x, MDS$y, labels, cex=0.6, pos=4)
    dev.off()
}

pvca <- function(exprs, covrts, results_path, figres=100) {
    #Principal Component analysis of variation (PVCA)

    inpData <- expSetobj(exprs, covrts)
    cvrts_eff_var <- colnames(covrts) ## Set the covariates whose effect size on the data needs to be calculated
    cvrts_eff_var_f <- c()
    for (val in cvrts_eff_var) {
        val1 <- lapply(val, tolower)
        test <- grep('id', val1)
        if (length(test) == 1) { 
        } else { 
            cvrts_eff_var_f <- append(cvrts_eff_var_f, val, after = length(cvrts_eff_var_f))
        }
    }
    # cvrts_eff_var <- cvrts_eff_var[-1]
    pct_thrsh <- 0.75
    png(file.path(results_path,  '2.pvca_vnorm_matrix.png'), res = figres)
    pvcAnaly(inpData, pct_thrsh, cvrts_eff_var_f)
    dev.off()
}


pca<-function(exprs, labels, results_path, figres=100) {
    normcounts <-exprs
    normcounts300<-normcounts[sample(1:nrow(normcounts),300),]
    pca <- prcomp(t(normcounts))
    cx <- sweep(t(normcounts), 2, colMeans(t(normcounts)), "-")
    sv <- svd(cx)
    vizualize_pca1(file.path(results_path, '2.svd_vnorm_matrix.png'), sv$u, labels$Sex, labels$Animal_Outcome, figres)
    vizualize_pca1(file.path(results_path, '2.pca_vnorm_matrix.png'), pca$x, labels$Sex, labels$Animal_Outcome, figres)    
    vizualize_scree_plot(file.path(results_path, '2.scree_vnorm_matrix.png'), pca, figres)
    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[with(loadingscores, order(-PC1)),]
    save_loading_scores(file.path(results_path, '2.loadingscores_pc1.txt'), head(loadingscores['PC1'],20), figres)
    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[with(loadingscores, order(-PC2)),] 
    save_loading_scores(file.path(results_path, '2.loadingscores_pc2.txt'), head(loadingscores['PC2'],20), figres)

}

# vizualize_pca<-function(plot_file, PCA, class1, class2, figres) { 
#     minx<-min(PCA[[1]]$PC1)
#     maxx<-max(PCA[[1]]$PC1)
#     miny<-min(PCA[[1]]$PC2)
#     maxy<-max(PCA[[1]]$PC2)
#     png(plot_file, res = figres)
#     plot(PCA[[1]]$PC1, PCA[[1]]$PC2, frame=FALSE, ylim=c(miny-5, maxy+5),
#          xlim=c(minx-5, maxx+5), pch=as.numeric(as.factor(class1)), col=as.numeric(as.factor(class2)) )

#     legend("topright", bty = "n", pch=as.numeric(as.factor(as.numeric(as.factor(class1)))),
#            legend= levels(as.factor(class1)))
#     legend("bottomright", bty = "n", pch='-', col=levels(as.factor(as.numeric(as.factor(class2)))),
#             legend= c(levels(as.factor(class2))))
#     dev.off()
# }
vizualize_pca1<-function(plot_file, PCA, class1, class2, figres) { 
    minx<-min(PCA[,1])
    maxx<-max(PCA[,1])
    miny<-min(PCA[,2])
    maxy<-max(PCA[,2])
    png(plot_file, res = figres)
    plot(PCA[,1], PCA[,2], frame=FALSE, ylim=c(miny, maxy), xlim=c(minx, maxx),
         pch=as.numeric(as.factor(class1)), col=as.numeric(as.factor(class2)) )
    legend("topright", inset=c(-1,-0.2), bty = "n", pch=as.numeric(as.factor(as.numeric(as.factor(class1)))),
           legend= levels(as.factor(class1)))
    legend("bottomright", bty = "n", pch='-', col=levels(as.factor(as.numeric(as.factor(class2)))),
           legend= c(levels(as.factor(class2))))
    dev.off()
}

vizualize_scree_plot<-function(plot_file, PCA, figres) {
    scree.plot<-fviz_eig(PCA, addlabels=TRUE, hjust = -0.3)
    png(plot_file, res = figres)
    print(scree.plot)
    dev.off()
}

save_loading_scores<-function(write_file, df, figres) { 
    write.table(df, file=write_file)
}
