#' s1 normalize raw counts function following Gale lab protocol
#'
#' This function allows the normalization of raw counts following a general protocol developed by previous members of the gale lab
#' @param countfile raw counts table (generally output by htseq).
#' @param targetfile target file.
#' @param vizualize_data whether or not to generate figures (default set to true).
#' @param  FilterGenesWithCounts filter out genes with counts below a certain value, (default set to 100).
#' @param  figres resolution at which to output figures (default is 300).
#' @keywords gene expression normalization
#' @export
#' @import limma
#' @import edgeR
#' @examples
#' normalize_raw_counts(count_file.txt, target_file.csv, vizualize_data=TRUE, FilterGenesWithCounts=100)

s1_normalize_raw_counts <- function(countfile, targetfile, visualize_data = TRUE, FilterGenesWithCounts=100, figres=100) { 
    ###READ IN FILES
    print("STATUS: loading files")
    files <- loadfiles(count_file=countfile, targets_file=targetfile)
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

        ###RUN VOOM
        print("STATUS: running voom")
        png('1.voomplot.png', res=figres)
        par(mar=c(1,1,1,1))
        V.CPM = voom(DE_DF_fl_norm, design=design, plot=T, span=0.1)
        dev.off()

        ###GENERATEOUTPUT DIRECTORY IN CURRENT WORKING DIRECTORY 
        workDir <-getwd()
        subDir='s1_norm_raw_counts_results'
        results_path = file.path(workDir, subDir)
        if (file.exists(subDir)){
        } else { 
            dir.create(results_path)
        }
        if (visualize_data == TRUE) { 
            print("STATUS: generating figures")
            vizualize_counts(files$counts, files$targets$treatment, count_matrix_flv, figres=figres, results_path=results_path)
        }

        ###save normalized counts and design variable used for linear modeling later
        write.csv(V.CPM$E, file=file.path(results_path,"1.matrix_norm.csv"))
        saveRDS(design, file=file.path(results_path, "1.designfile.rds"))

        results_norm <- list("norm_exprs_voom" = V.CPM, "design" = design)

        return (results_norm)
    }
}

###OTHER FUNCTIONS USED BY normalize_raw_counts
loadfiles <- function(count_file, targets_file) {
    #Load in count data  and target information from csv'
    counts <- read.table(count_file, header = TRUE, sep = "\t", row.names=1, as.is = TRUE)
    targets <- read.table(targets_file, header = TRUE, sep = ",", row.names = 1, as.is = TRUE)
    results <- list("counts" = counts, "targets" = targets)
    return(results)
}

vizualize_counts <- function(countsmatrix, labels, count_matrix_flv, figres=100, results_path) {
    #Generate figures for counts
    print('STATUS: generating log2 boxplot of counts')
    png(file.path(results_path, "1.boxplot_count_matrix.png"), res=figres)
    # par(mar=c(1,1,1,1))
    boxplot(log2(countsmatrix+1), labels = labels, ylab = "log2 Expression", main = "Raw count_matrix", cex.axis=.6, las=2)
    dev.off()

    print('STATUS: generating density plot of all sample counts')
    png(file.path(results_path, "1.densities_count_matrix.png"), res=figres)
    # par(mar=c(1,1,1,1))
    if (length(labels) > 10) {
        plotDensities(log2(countsmatrix+1), legend = FALSE)    
    } else {
        plotDensities(log2(countsmatrix+1), legend = "topright", levels(labels))
    }
    dev.off()

    print('STATUS: generating biological varation vs abundance')
    png(file.path(results_path, "1.biologicalcoefficentvariation.png"), res=figres)
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




