#' s1 normalize raw counts function following Gale lab protocol
#'
#' This function allows the normalization of raw counts following a general protocol developed by previous members of the gale lab
#' @param countfile raw counts table (generally output by htseq).
#' @param targetfile target file.
#' @param target_class columns from the target file to build design matrix for future DE analysis 
#' @param vizualize_data whether or not to generate figures (default set to true).
#' @param filter_genes_below_counts filter out genes with counts below a certain value, (default set to 0).
#' @param figres resolution at which to output figures (default is 300).
#' @keywords gene expression normalization
#' @export
#' @import limma
#' @import edgeR
#' @import edgeR
#' @examples
#' normalize_raw_counts(count_file.txt, target_file.csv, vizualize_data=TRUE, FilterGenesWithCounts=100)

s1_normalize_raw_counts <- function(countfile, targetfile, gene_conversion_file=FALSE, target_class=10, blocking_column=FALSE, visualize_data=TRUE, filter_genes_below_counts=0, figres=100) { 
    ###READ IN FILES
    print("STATUS: loading files")
    files            <- loadfiles(count_file=countfile, target_file=targetfile)
    DE_DF            <- DGEList(counts = files$counts)

    #FILTER OUT GENES WITH LOW COUNTS
    print("STATUS: filtering out genes with low counts")
    A                <- rowSums(DE_DF$counts)
    isexpr           <- A > filter_genes_below_counts
    DE_DF            <- DE_DF[isexpr,]
    DE_DF_fl         <- calcNormFactors(DE_DF)

    ###get biological coefficients of variation
    print("STATUS: getting biological coefficient of variation (takes time...)")
    count_matrix_flv <- biological_coefficents_variation(DE_DF_fl)
    
    ###NORMALIZE VIA TMM
    print("STATUS: getting normalizing factors (method TMM)")
    DE_DF_fl_norm    <- calcNormFactors(DE_DF_fl)

    ###SET UP MODEL DESIGN
    print("STATUS: setting up model design")

    if (length(files$targets$treatment) != length(colnames(DE_DF)) | all.equal(rownames(files$targets), colnames(files$counts)) != TRUE) {
        print ('WARNING: different number of treatments and column names in count file or order of samples in target file does not match order in count file (needs fixing before we can proceed)')
        print (paste0("Length of treatments:", length(files$targets$treatment)))
        print (paste0("Length of column names in count/normalized matrix:", length(colnames(DE_DF))))
    }
    else {
        results_path <- generate_folder('s1_norm_raw_counts_results')
        unlink('./s1_norm_raw_counts_results/*')
        factors<-list()
        # for (i in target_class) {
        treatment    <- factor(files$targets[,target_class], levels=unique(files$targets[,target_class]))
            # factors <- list.append(factors, i = F)
        # }
        design       <- model.matrix(~0 + treatment)

        if (is.fullrank(design) == TRUE & is.null(nonEstimable(design))) { 
            # colnames(design) <- levels(CLASS1)
            if (blocking_column != FALSE) {
                BLOCKID  <- factor(files$targets[,blocking_column], levels=unique(files$targets[,blocking_column]))
                corfit   <- duplicateCorrelation(DE_DF_fl_norm$counts,design,block=BLOCKID)
            }
            ###RUN VOOM
            print("STATUS: running voom")
            png(file.path(results_path,'1.voomplot.png'), res=figres)
            # par(mar=c(1,1,1,1))
            V.CPM        <- voomWithQualityWeights(DE_DF_fl_norm, design=design, plot=T, span=0.1)
            dev.off()

            ###save normalized counts and design variable used for linear modeling later
            write.table(data.frame(V.CPM$E), sep='\t', col.names=NA, file=file.path(results_path,"1.norm_matrix.txt"))
            norm_matrix  <- V.CPM$E
            
            if (typeof(gene_conversion_file) == 'character') {
                rhesus2human <- read.csv(file=gene_conversion_file, header=TRUE, stringsAsFactors = FALSE)
                nm_hgnc      <- merge(rhesus2human, norm_matrix, by.x='Gene.stable.ID', by.y='row.names')
                nm_hgnc      <- avereps(nm_hgnc, ID = nm_hgnc$Gene.stable.ID)
                write.table(nm_hgnc, sep='\t', file=file.path(results_path, "1.norm_matrix_HGNC.txt"))
            }

            saveRDS(V.CPM, file.path(results_path, "1.voomobject.rds"))

            if (blocking_column != FALSE) {
                saveRDS(corfit, file.path(results_path, "1.corfit.rds"))
            }
            saveRDS(design, file=file.path(results_path, "1.designobject.rds"))

            if (visualize_data == TRUE) { 
                print("STATUS: generating figures")
                vizualize_counts(files$counts, V.CPM$E, files$targets, count_matrix_flv, figres=figres, results_path=results_path)
            }

            # results_norm <- list("norm_exprs_voom" = V.CPM, "design" = design)

            # return (results_norm)
        } else { 
            print ('WARNING: error with design matrix... rethink how it is being set up')

        }

    }
}

###OTHER FUNCTIONS USED BY normalize_raw_counts
vizualize_counts <- function(countsmatrix, norm_exprs, labels, count_matrix_flv, figres=100, results_path) {
    #Generate figures for counts
    print('STATUS: generating log2 boxplot of counts')
    generate_boxplots(log2(countsmatrix+1), labels[,1], file.path(results_path, "1.boxplot_raw_count_matrix.png"), figres, maintitle='Raw count matrix', ylabtitle='log2 Expression')

    print('STATUS: generating boxplot of normalized voom counts')
    generate_boxplots(norm_exprs, labels[,1], file.path(results_path, "1.boxplot_vnorm_matrix.png"), figres, maintitle='Normalized count matrix', ylabtitle='voom normalized expression')

    print('STATUS: generating density plot of all log counts')
    generate_density_plot(log2(countsmatrix+1), labels[,1], file.path(results_path, "1.densities_raw_log_counts.png"), figres)

    print('STATUS: generating density plot of raw counts')
    generate_density_plot(countsmatrix, labels[,1], file.path(results_path, "1.densities_raw_counts.png"), figres)

    print('STATUS: generating density plot of normalized voom counts')
    generate_density_plot(norm_exprs, labels[,1], file.path(results_path, "1.densities_vnorm_matrix.png"), figres)


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

generate_density_plot <- function(data, labels, filename, figres){ 
    png(filename, res=figres)
    par( xpd=TRUE)
    if (length(labels) > 10) {
        plotDensities(data, legend = FALSE)    
    } else {
        plotDensities(data, legend = "topright", inset=c(-0.2,0), levels(labels))
    }
    dev.off()
}



generate_boxplots <- function(data, labels, filename, figres, maintitle, ylabtitle){ 
    png(filename, res=figres)
    # par(mar=c(1,1,1,1))
    minvalue <-min(data)
    maxvalue <- max(data)
    boxplot(data, labels = labels, ylim=c(minvalue-1, maxvalue+1), ylab = ylabtitle, main = maintitle, cex.axis=.6, las=2, frame=FALSE)
    dev.off()

}    