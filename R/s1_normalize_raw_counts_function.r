#' s1 normalize raw counts function following Gale lab protocol
#'
#' This function allows the normalization of raw counts following a general protocol developed by previous members of the gale lab
#' @param countfile raw counts table (generally output by htseq).
#' @param targetfile target file.
#' @param vizualize_data whether or not to generate figures (default set to true).
#' @param  FilterGenesWithCounts filter out genes with counts below a certain value, (default set to 0).
#' @param  figres resolution at which to output figures (default is 300).
#' @keywords gene expression normalization
#' @export
#' @import limma
#' @import edgeR
#' @examples
#' normalize_raw_counts(count_file.txt, target_file.csv, vizualize_data=TRUE, FilterGenesWithCounts=100)

s1_normalize_raw_counts <- function(countfile, targetfile, visualize_data = TRUE, FilterGenesWithCounts=0, figres=100) { 
    ###READ IN FILES
    print("STATUS: loading files")
    files <- loadfiles(count_file=countfile, target_file=targetfile)
    DE_DF <- DGEList(counts = files$counts)

    #FILTER OUT GENES WITH LOW COUNTS
    print("STATUS: filtering out genes with low counts")
    DE_DF_fl = DE_DF[apply(FUN=max, X=DE_DF, MARGIN=1)>FilterGenesWithCounts,]

    ###get biological coefficients of variation
    print("STATUS: getting biological coefficient of variation (takes time...)")
    count_matrix_flv <- biological_coefficents_variation(DE_DF_fl)
    
    ###NORMALIZE VIA TMM
    print("STATUS: getting normalizing factors (method TMM)")
    DE_DF_fl_norm <- calcNormFactors(DE_DF_fl)

    ###SET UP MODEL DESIGN
    print("STATUS: setting up model design")

    if (length(files$targets$treatment) != length(colnames(DE_DF))) {
        print ('WARNING: different number of treatments and column names in count file (needs fixing before we can proceed)')
        print (paste0("Length of treatments:", length(files$targets$treatment)))
        print (paste0("Length of column names in count/normalized matrix:", length(colnames(DE_DF))))
    }
    else {
        Treatment <- factor(files$targets$treatment, levels=unique(files$targets$treatment))
        design <- model.matrix(~0 + Treatment)
        colnames(design) <- levels(Treatment)
        results_path <- generate_folder('s1_norm_raw_counts_results')
        
        ###RUN VOOM
        print("STATUS: running voom")
        png(file.path(results_path,'1.voomplot.png'), res=figres)
        par(mar=c(1,1,1,1))
        V.CPM = voom(DE_DF_fl_norm, design=design, plot=T, span=0.1)
        dev.off()

        ###save normalized counts and design variable used for linear modeling later
        write.csv(V.CPM$E, file=file.path(results_path,"1.matrix_norm.csv"))
        saveRDS(design, file=file.path(results_path, "1.designfile.rds"))

        if (visualize_data == TRUE) { 
            print("STATUS: generating figures")
            vizualize_counts(files$counts, V.CPM$E, files$targets$treatment, count_matrix_flv, figres=figres, results_path=results_path)
        }


        results_norm <- list("norm_exprs_voom" = V.CPM, "design" = design)

        return (results_norm)
    }
}

###OTHER FUNCTIONS USED BY normalize_raw_counts
vizualize_counts <- function(countsmatrix, norm_exprs, labels, count_matrix_flv, figres=100, results_path) {
    #Generate figures for counts
    print('STATUS: generating log2 boxplot of counts')
    png(file.path(results_path, "1.boxplot_raw_count_matrix.png"), res=figres)
    # par(mar=c(1,1,1,1))
    minvalue <-min(log2(countsmatrix+1))
    maxvalue <- max(log2(countsmatrix+1))
    boxplot(log2(countsmatrix+1), labels = labels, ylim=c(minvalue-5, maxvalue+5), ylab = "log2 Expression", main = "Raw count matrix", cex.axis=.6, las=2, frame=FALSE)
    dev.off()

    print('STATUS: generating boxplot of normalized voom counts')
    png(file.path(results_path, "1.boxplot_vnorm_matrix.png"), res=figres)
    # par(mar=c(1,1,1,1))
    minvalue <-min(norm_exprs)
    maxvalue <- max(norm_exprs)
    boxplot(norm_exprs, labels = labels, ylim=c(minvalue-1, maxvalue+1), ylab = "log2 Expression", main = "Raw count matrix", cex.axis=.6, las=2, frame=FALSE)
    dev.off()


    print('STATUS: generating density plot of all sample counts')
    png(file.path(results_path, "1.densities_raw_count_matrix.png"), res=figres)
    # par(mar=c(1,1,1,1))
    if (length(labels) > 10) {
        plotDensities(log2(countsmatrix+1), legend = FALSE)    
    } else {
        plotDensities(log2(countsmatrix+1), legend = "topright", levels(labels))
    }
    dev.off()

    print('STATUS: generating density plot of normalized voom counts')
    png(file.path(results_path, "1.densities_vnorm_matrix.png"), res=figres)
    # par(mar=c(1,1,1,1))
    if (length(labels) > 10) {
        plotDensities(norm_exprs, legend = FALSE)    
    } else {
        plotDensities(norm_exprs, legend = "topright", levels(labels))
    }
    dev.off()


    print('STATUS: generating biological varation vs abundance')
    png(file.path(results_path, "1.biologicalcoefficentvariation_raw.png"), res=figres)
    # par(mar=c(1,1,1,1))
    plotBCV(count_matrix_flv, cex=0.4, main="Biological coefficient of variation (BCV) vs abundance")
    dev.off()
}

# Biological coefficients of variation

biological_coefficents_variation <- function(count_matrix_fl) {
    count_matrix_flv <- estimateCommonDisp(count_matrix_fl,verbose=T) #print the BCV value
    count_matrix_flv <- estimateTrendedDisp(count_matrix_flv)
    count_matrix_flv <- estimateTagwiseDisp(count_matrix_flv)
    return(count_matrix_flv)
}




