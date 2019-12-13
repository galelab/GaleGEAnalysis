#' s2 feature reduction 
#'
#' This function does feature reduction and visualization on normalized voom counts
#' @param countfile normalized counts table (generally should have been generated by first step s1_normalize_raw_counts, but also this function can be run on raw cout file).
#' @param targetfile target file.
#' @param target_columns columns from the target file to label samples on the pca/umap plots (has to be 2)
#' @param figres resolution at which to output figures (default is 300).
#' @param pca run principal component analysis (and singular vector decomposition) (default set to TRUE)
#' @param UMAP run Uniform Manifold Aproximation Projection (default set to TRUE)
#' @param base_file_name file name for all figures
#' @keywords gene expression feature reduction
#' @export
#' @import Biobase
#' @import ExpressionNormalizationWorkflow
#' @import stringr
#' @import stats
#' @import factoextra
#' @import umap
#' @import Rtsne
#' @import png
#' @examples
#' s2_feature_reduction(countfile='./s1_norm_raw_counts_results/1.norm_matrix.txt',  targetfile='./p1_modified_count_matrix_results/target_file.csv', target_columns=c(2,5), base_file_name='vnorm.png')

s2_feature_reduction <- function(countfile='./s1_norm_raw_counts_results/1.norm_matrix.txt', targetfile='./p1_modified_count_matrix_results/targets_mod.csv', target_columns=c(2,5), figres=100, pca=TRUE, UMAP=TRUE, tsne=TRUE, base_file_name='vnorm.png') { 
    pdf(NULL)
    files <- loadfiles(count_file=countfile, target_file=targetfile)
    results_path <- generate_folder('s2_feature_reduction_results')
    
    print ('STATUS: Running MDS feature reduction')
    vizualize_feature_reduction_data(files$counts, files$targets[,target_columns[1]], results_path, base_file_name, figres)
    if (all.equal(rownames(files$targets), colnames(files$counts)) != TRUE) { 
        print ('WARNING: order of samples in target and count file is not the same, this needs fixing before code can proceed')
    } else {
        if (pca == TRUE) {
            print ('STATUS: Running PCA feature reduction')
            pca_fun(files$counts, files$target, results_path, base_file_name, target_columns, figres)
        }
        
        if (UMAP == TRUE) {
            print ('STATUS: Running UMAP feature reduction')
            umap_fun(files$counts, files$target, results_path, base_file_name, target_columns, figres)
        }
        if (tsne == TRUE) {
           print ('STATUS: Running tSNE feature reduction')
           tsne_fun(files$counts, files$target, results_path, base_file_name, target_columns, figres)
        }
    }
}

vizualize_feature_reduction_data <- function(data, labels, results_path, base_file_name, figres=100) { 
    ###MDS (multidimensional scaling) uses log fold changes between genes as distances
    MDS <- plotMDS(data, gene.selection='pairwise', cex= .8)
    minx<-min(MDS$x)
    maxx<-max(MDS$x)
    miny<-min(MDS$y)
    maxy<-max(MDS$y)
    png(file.path(results_path,  paste0('2.mds_',base_file_name)), res = figres)
    plot(MDS$x, MDS$y, cex=1, xlim=c(minx-1, maxx+1), ylim=c(miny-1, maxy+1), xlab=paste0(MDS$axislabel,' 1'), ylab=paste0(MDS$axislabel,' 2'), frame = FALSE)
    text(MDS$x, MDS$y, labels, cex=0.6, pos=4)
    dev.off()
}

pvca_fun <- function(exprs, covrts, results_path, base_file_name, target_columns, figres=100) {
    #Principal Component analysis of variation (PVCA)

    inpData <- expSetobj(exprs, covrts)
    cvrts_eff_var <- colnames(covrts) ## Set the covariates whose effect size on the data needs to be calculated
 
    # for (val in cvrts_eff_var) {
    #     val1 <- lapply(val, tolower)
    #     test <- grep('id', val1)
    #     if (length(test) == 1) { 
    #     } else { 
    #         cvrts_eff_var_f <- append(cvrts_eff_var_f, val, after = length(cvrts_eff_var_f))
    #     }
    # }
    cvrts_eff_var <- cvrts_eff_var[-1]
    pct_thrsh <- 0.75
    png(file.path(results_path, paste0('2.pvca_', base_file_name)), res = figres)
    pvcAnaly(inpData, pct_thrsh, cvrts_eff_var)
    dev.off()
}


pca_fun<-function(exprs, labels, results_path, base_file_name, target_columns, figres=100) {
    #Run PCA/SVD reduction 

    pca <- prcomp(t(exprs))
    E <- get_eig(pca)
    cx <- sweep(t(exprs), 2, colMeans(t(exprs)), "-")
    sv <- svd(cx)

    vizualize_pca(file.path(results_path, paste0('2.svd_',base_file_name)), sv$u, labels[,target_columns[1]], labels[,target_columns[2]], figres, E)
    vizualize_pca(file.path(results_path, paste0('2.pca_',base_file_name)), pca$x, labels[,target_columns[1]], labels[,target_columns[2]], figres, E)    
    vizualize_scree_plot(file.path(results_path, paste0('2.scree_',base_file_name)), pca, figres)

    loadingscores <- as.data.frame(pca$rotation)
    is_pc1_0 <- loadingscores$PC1 > 0
    is_pc2_0 <- loadingscores$PC2 > 0

    loadingscores <- loadingscores[is_pc1_0,]
    loadingscores <- loadingscores[with(loadingscores, order(-PC1)),]
    save_loading_scores(file.path(results_path, '2.loadingscores_pc1.txt'), loadingscores['PC1'], figres)

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[is_pc2_0,]
    loadingscores <- loadingscores[with(loadingscores, order(-PC2)),] 
    save_loading_scores(file.path(results_path, '2.loadingscores_pc2.txt'), loadingscores['PC2'], figres)

}

umap_fun<-function(exprs, labels, results_path, base_file_name, target_columns, figres=100) {
    #Runs default paramaters of umap 
    U <- umap(t(exprs))
    vizualize_umap(file.path(results_path, paste0('2.umap_',base_file_name)), U$layout, labels[,target_columns[1]], labels[,target_columns[2]], figres)
}

vizualize_umap<-function(plot_file, U, class1, class2, figres) { 
    #Vizualize umap reduction
    minx<-min(U[,1])
    maxx<-max(U[,1])
    miny<-min(U[,2])
    maxy<-max(U[,2])
    png(plot_file, res = figres)
    par(mar = c(5, 4, 2, 4), xpd=TRUE)
    plot(U[,1], U[,2], frame=FALSE, ylim=c(miny-1, maxy+1), xlim=c(minx-1, maxx+1),  pch=as.numeric(as.factor(class1)), col=as.numeric(as.factor(class2)), 
    xlab='Dim 1', ylab='Dim 2')
    legend("topright", inset=c(-0.25,-0.1), bty = "n", pch=as.numeric(levels(as.factor(as.numeric(as.factor(class1))))),
           legend= levels(as.factor(class1)))
    legend("bottomright", inset=c(-0.25, 0), bty = "n", pch='-', col=levels(as.factor(as.numeric(as.factor(class2)))),
           legend= c(levels(as.factor(class2))))
    dev.off()
}


vizualize_pca<-function(plot_file, PCA, class1, class2, figres, E) { 
    #Vizualize PCA  results 
    minx<-min(PCA[,1])
    maxx<-max(PCA[,1])
    miny<-min(PCA[,2])
    maxy<-max(PCA[,2])
    png(plot_file, res = figres)
    par(mar = c(5, 4, 2, 5.5), xpd=TRUE)
    plot(PCA[,1], PCA[,2], frame=FALSE, ylim=c(miny, maxy), xlim=c(minx, maxx),
         pch=as.numeric(as.factor(class1)), col=as.numeric(as.factor(class2)), xlab=paste0('PC1 ',round(E$variance.percent[1],  digits = 2), '%'),
         ylab=paste0('PC2 ',round(E$variance.percent[2],  digits = 2), '%'  ))
    legend("topright", inset=c(-0.35,-0.1), bty = "n", pch=as.numeric(levels(as.factor(as.numeric(as.factor(class1))))),
           legend= levels(as.factor(class1)))
    legend("bottomright", inset=c(-0.37, 0), bty = "n", pch='-', col=levels(as.factor(as.numeric(as.factor(class2)))),
           legend= c(levels(as.factor(class2))))
    dev.off()
}

vizualize_scree_plot<-function(plot_file, PCA, figres) {
    #Vizualize principle component variation results 
    scree.plot<-fviz_eig(PCA, addlabels=TRUE, hjust = -0.3)
    png(plot_file, res = figres)
    print(scree.plot)
    dev.off()
}

save_loading_scores<-function(write_file, df, figres) { 
    #Save list of genes that have a positive effect on variation of principle component 1 and 2 sorted from most influential 
    write.table(df, file=write_file)
}

tsne_fun<-function(exprs, labels, results_path, base_file_name, target_columns, figres=100) {
    #Runs default paramaters of umap 
    T <- Rtsne(t(exprs), perplexity=1)
    vizualize_tSNE(file.path(results_path, paste0('2.tsne_',base_file_name)), T$Y, labels[,target_columns[1]], labels[,target_columns[2]], figres)
}

vizualize_tSNE<-function(plot_file, U, class1, class2, figres) { 
    #Vizualize umap reduction
    minx<-min(U[,1])
    maxx<-max(U[,1])
    miny<-min(U[,2])
    maxy<-max(U[,2])
    png(plot_file, res = figres)
    par(mar = c(5, 4, 2, 4), xpd=TRUE)
    plot(U[,1], U[,2], frame=FALSE, ylim=c(miny-1, maxy+1), xlim=c(minx-1, maxx+1),  pch=as.numeric(as.factor(class1)), col=as.numeric(as.factor(class2)), 
    xlab='Dim 1', ylab='Dim 2')
    legend("topright", inset=c(-0.25,-0.1), bty = "n", pch=as.numeric(levels(as.factor(as.numeric(as.factor(class1))))),
           legend= levels(as.factor(class1)))
    legend("bottomright", inset=c(-0.25, 0), bty = "n", pch='-', col=levels(as.factor(as.numeric(as.factor(class2)))),
           legend= c(levels(as.factor(class2))))
    dev.off()
}