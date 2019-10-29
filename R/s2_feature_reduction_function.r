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

s2_feature_reduction <- function(countfile, targetfile, figres=100, pcva=TRUE, pca=TRUE, umap=TRUE) { 
    # V.CPM <- readRDS('./s1_norm_raw_counts_results/1.voomobject.rds')
    pdf(NULL)
    files <- loadfiles(count_file=countfile, target_file=targetfile)
    results_path <- generate_folder('s2_feature_reduction_results')
    vizualize_feature_reduction_data(files$counts, files$target$treatment,results_path, figres)
    if (pcva == TRUE) {
        pvca_fun(files$counts, files$target, results_path)
    }
    if (pca == TRUE) {
        pca_fun(files$counts, files$target, results_path)
    }
    umap_fun(files$counts, files$target, results_path)
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

pvca_fun <- function(exprs, covrts, results_path, figres=100) {
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


pca_fun<-function(exprs, labels, results_path, figres=100) {
    normcounts <-exprs
    normcounts300<-normcounts[sample(1:nrow(normcounts),300),]

    pca <- prcomp(t(normcounts))
    E <- get_eig(pca)
    cx <- sweep(t(normcounts), 2, colMeans(t(normcounts)), "-")
    sv <- svd(cx)

    vizualize_pca1(file.path(results_path, '2.svd_vnorm_matrix.png'), sv$u, labels$Sex, labels$Animal_Outcome, figres, E)
    vizualize_pca1(file.path(results_path, '2.pca_vnorm_matrix.png'), pca$x, labels$Sex, labels$Animal_Outcome, figres, E)    
    vizualize_scree_plot(file.path(results_path, '2.scree_vnorm_matrix.png'), pca, figres)

    vizualize_pca1(file.path(results_path, '2.svd_vnorm_matrix1.png'), sv$u, labels$Sex, labels$animal.ID, figres, E)
    vizualize_pca1(file.path(results_path, '2.pca_vnorm_matrix1.png'), pca$x, labels$Sex, labels$animal.ID, figres, E)

    vizualize_pca1(file.path(results_path, '2.svd_vnorm_matrix2.png'), sv$u, labels$Sex, labels$Time_Point, figres, E)
    vizualize_pca1(file.path(results_path, '2.pca_vnorm_matrix2.png'), pca$x, labels$Sex, labels$Time_Point, figres, E)

    vizualize_pca1(file.path(results_path, '2.svd_vnorm_matrix3.png'), sv$u, labels$Study_Group, labels$treatment, figres, E)
    vizualize_pca1(file.path(results_path, '2.pca_vnorm_matrix3.png'), pca$x, labels$Study_Group, labels$treatment, figres, E)

    vizualize_pca1(file.path(results_path, '2.svd_vnorm_matrix4.png'), sv$u, labels$Study_Group, labels$Animal_Outcome, figres, E)
    vizualize_pca1(file.path(results_path, '2.pca_vnorm_matrix4.png'), pca$x, labels$Study_Group, labels$Animal_Outcome, figres, E)

    vizualize_pca1(file.path(results_path, '2.svd_vnorm_matrix5.png'), sv$u, labels$animal.ID, labels$Animal_Outcome, figres, E)
    vizualize_pca1(file.path(results_path, '2.pca_vnorm_matrix5.png'), pca$x, labels$animal.ID, labels$Animal_Outcome, figres, E)

    vizualize_pca1(file.path(results_path, '2.svd_vnorm_matrix6.png'), sv$u, labels$Number_Challenges, labels$Animal_Outcome, figres, E)
    vizualize_pca1(file.path(results_path, '2.pca_vnorm_matrix6.png'), pca$x, labels$Number_Challenges, labels$Animal_Outcome, figres, E)

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[with(loadingscores, order(-PC1)),]
    save_loading_scores(file.path(results_path, '2.loadingscores_pc1.txt'), head(loadingscores['PC1'],20), figres)
    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[with(loadingscores, order(-PC2)),] 
    save_loading_scores(file.path(results_path, '2.loadingscores_pc2.txt'), head(loadingscores['PC2'],20), figres)

}

umap_fun<-function(exprs, labels, results_path, figres=100) { 
    U <- umap(t(exprs))
    vizualize_umap(file.path(results_path, '2.umap_reduction.png'), U$layout, labels$Sex, labels$Animal_Outcome, figres)
    vizualize_umap(file.path(results_path, '2.umap_reduction1.png'), U$layout, labels$Sex, labels$animal.ID, figres)
    vizualize_umap(file.path(results_path, '2.umap_reduction2.png'), U$layout, labels$Sex, labels$Time_Point, figres)
    vizualize_umap(file.path(results_path, '2.umap_reduction3.png'), U$layout, labels$Study_Group, labels$treatment, figres)
    vizualize_umap(file.path(results_path, '2.umap_reduction4.png'), U$layout, labels$Study_Group, labels$Animal_Outcome, figres)
    vizualize_umap(file.path(results_path, '2.umap_reduction5.png'), U$layout, labels$animal.ID, labels$Animal_Outcome, figres)
    vizualize_umap(file.path(results_path, '2.umap_reduction6.png'), U$layout, labels$Number_Challenges, labels$Animal_Outcome, figres)

}

vizualize_umap<-function(plot_file, U, class1, class2, figres) { 
    minx<-min(U[,1])
    maxx<-max(U[,1])
    miny<-min(U[,2])
    maxy<-max(U[,2])
    png(plot_file, res = figres)
    par(mar = c(5, 4, 2, 4), xpd=TRUE)
    plot(U[,1], U[,2], frame=FALSE, ylim=c(miny, maxy), xlim=c(minx, maxx),  pch=as.numeric(as.factor(class1)), col=as.numeric(as.factor(class2)), 
    xlab='Dim 1', ylab='Dim 2')
    legend("topright", inset=c(-0.25,-0.1), bty = "n", pch=as.numeric(levels(as.factor(as.numeric(as.factor(class1))))),
           legend= levels(as.factor(class1)))
    legend("bottomright", inset=c(-0.3, 0), bty = "n", pch='-', col=levels(as.factor(as.numeric(as.factor(class2)))),
           legend= c(levels(as.factor(class2))))
    dev.off()
}


vizualize_pca1<-function(plot_file, PCA, class1, class2, figres, E) { 
    minx<-min(PCA[,1])
    maxx<-max(PCA[,1])
    miny<-min(PCA[,2])
    maxy<-max(PCA[,2])
    png(plot_file, res = figres)
    par(mar = c(5, 4, 2, 5.5), xpd=TRUE)
    plot(PCA[,1], PCA[,2], frame=FALSE, ylim=c(miny, maxy), xlim=c(minx, maxx),
         pch=as.numeric(as.factor(class1)), col=as.numeric(as.factor(class2)), xlab=paste0('PC1 ',round(E$variance.percent[1],  digits = 2), '%'),
         ylab=paste0('PC2 ',round(E$variance.percent[2],  digits = 2), '%'  ))
    legend("topright", inset=c(-0.25,-0.1), bty = "n", pch=as.numeric(levels(as.factor(as.numeric(as.factor(class1))))),
           legend= levels(as.factor(class1)))
    legend("bottomright", inset=c(-0.4, 0), bty = "n", pch='-', col=levels(as.factor(as.numeric(as.factor(class2)))),
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
