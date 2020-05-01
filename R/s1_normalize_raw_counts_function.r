#' s1 normalize raw counts function following Gale lab protocol
#'
#' This function allows the normalization of raw counts following a general protocol developed by previous members of the gale lab
#' @param countfile raw counts table (generally output by htseq).
#' @param targetfile target file.
#' @param target_column columns from the target file to build design matrix for future DE analysis
#' @param blocking_column column to account sampling from the same animal multiple times
#' @param vizualize_data whether or not to generate figures (default set to true).
#' @param filter_genes_below_counts filter out genes with counts below a certain value, (default set to 0).
#' @param norm_method normalization method to use when normalizing LogCPM
#' @param results_folder folder to store results in defualt is s1_norm_raw_counts_results
#' @param figres resolution at which to output figures (default is 300).
#' @keywords gene expression normalization
#' @export
#' @import limma
#' @import edgeR
#' @examples
#' normalize_raw_counts(countfile="./p1_modified_count_matrix_results/count_file.txt", targetfile="./p1_modified_count_matrix_results/target_file.csv", gene_conversion_file="rhesus2human.csv", target_column=10, blocking_column=2, vizualize_data=TRUE, filter_genes_below_counts=0, figres=100)

s1_normalize_raw_counts <- function(countfile, targetfile,
                                    gene_conversion_file=FALSE,
                                    target_column=10, batch_column=FALSE,
                                    blocking_column=FALSE,
                                    visualize_data=TRUE,
                                    filter_genes_below_counts=0,
                                    norm_method="none",
                                    results_folder = "s1_norm_raw_counts_results",
                                    figres=150) {
    ###READ IN FILES
    print("STATUS: loading files")
    files <- loadfiles(count_file = countfile,
                      target_file = targetfile)
    DE_DF <- DGEList(counts = files$counts)

    #FILTER OUT GENES WITH LOW COUNTS
    print("STATUS: filtering out genes with low counts")
    A <- rowSums(DE_DF$counts)
    isexpr <- A >= filter_genes_below_counts
    DE_DF <- DE_DF[isexpr, ]

    ###NORMALIZE VIA TMM
    print("STATUS: getting normalizing factors (method TMM)")
    DE_DF_fl <- calcNormFactors(DE_DF)

    ###get biological coefficients of variation
    # print("STATUS: getting biological coefficient of variation (takes time...)")
    # count_matrix_flv <- biological_coefficents_variation(DE_DF_fl)

    ###SET UP MODEL DESIGN
    print("STATUS: setting up model design")

    if (length(files$targets[, target_column]) != length(colnames(DE_DF))) {
        print("WARNING: different number of treatments and column names in count file ")
        print(paste0("Length of treatments:",
                     length(files$targets[, target_column])))
        print (paste0("Length of column names in count/normalized matrix:", length(colnames(DE_DF))))
    }
    else if ( all.equal(rownames(files$targets), colnames(files$counts)) != TRUE) {
       print ("WARNING: Order of samples in target file does not match order in count file (needs fixing before we can proceed)")
    #    files$counts <- files$counts[, rownames(files$targets)]
    #    if (all.equal(rownames(files$targets), colnames(files$counts)) == TRUE) {
    #        print("WARNING: Order of samples has been corrected!")
    #    }
    }
    else {
        results_path <- generate_folder(results_folder)
        unlink(paste0(results_folder, "/*"))
        factors <- list()

        if (typeof(batch_column) == "logical") {
            treatment   <- factor(files$targets[, target_column],
                            levels = unique(files$targets[, target_column]))
            design      <- model.matrix(~0 + treatment)
            rownames(design) <- colnames(DE_DF$counts)
            colnames(design) <- make.names(colnames(design))
            design <- design[, colnames(design)[order(tolower(colnames(design[, ])))]]

        } else {
            treatment    <- factor(files$targets[, target_column],
                            levels = unique(files$targets[, target_column]))
            batch        <- factor(files$targets[, batch_column],
                            levels = unique(files$targets[, batch_column]))
            design       <- model.matrix(~0 + treatment + batch)
            rownames(design) <- colnames(DE_DF$counts)
            colnames(design) <- make.names(colnames(design))
            design <- design[, colnames(design)[order(tolower(colnames(design[, ])))]]

        }

        if (is.fullrank(design) == TRUE & is.null(nonEstimable(design))) {
            # colnames(design) <- levels(CLASS1)
            if (blocking_column != FALSE) {
                BLOCKID <- factor(files$targets[, blocking_column],
                                  levels = unique(files$targets[,blocking_column]))
                corfit <- duplicateCorrelation(DE_DF_fl$counts, design,
                                               block = BLOCKID)
            }
            ###RUN VOOM
            print("STATUS: running voom")
            png(file.path(results_path, "1.voomplot.png"), res = figres)

            #Transform count data to log2-counts per million
            V.CPM <- voom(DE_DF_fl, normalize.method = norm_method, design = design, plot = T, span = 0.1)
            dev.off()

            ###save normalized counts and design variable used for linear modeling later
            orig.cols    <- colnames(files$counts)
            # orig.cols    <- append(orig.cols, "Name", after=0)
            orig.rows    <- rownames(V.CPM$E)
            write.table(data.frame(V.CPM$E), sep = "\t",
                       row.names = orig.rows, col.names = orig.cols,
                       file = file.path(results_path, "1.norm_matrix.txt"))
            norm_matrix  <- V.CPM$E

            if (typeof(gene_conversion_file) == "character") {
                rhesus2human <- read.csv(file = gene_conversion_file,
                                         header = TRUE,
                                         stringsAsFactors = FALSE)
                nm_hgnc      <- merge(rhesus2human, norm_matrix,
                                      by.x = "Gene.stable.ID", by.y = "row.names")
                nm_hgnc      <- avereps(nm_hgnc, ID = nm_hgnc$Gene.stable.ID)
                write.table(nm_hgnc, sep = "\t",
                            file = file.path(results_path,
                                             "1.norm_matrix_HGNC.txt"))
            }

            saveRDS(V.CPM, file.path(results_path, "1.voomobject.rds"))

            if (blocking_column != FALSE) {
                saveRDS(corfit, file.path(results_path, "1.corfit.rds"))
            }
            saveRDS(design, file = file.path(results_path,
                                            "1.designobject.rds"))

            if (visualize_data == TRUE) {
                print("STATUS: generating figures")
                vizualize_counts(files$counts, V.CPM$E, files$targets,
                                 figres = figres,
                                 results_path = results_path)
            }

            # results_norm <- list("norm_exprs_voom" = V.CPM, "design" = design)

            # return (results_norm)
        } else {
            print("WARNING: error with design matrix... rethink how it is being set up")
            print(paste0("WARNING: full rank check is ", is.fullrank(design)))
            print(paste0("WARNING: nonEstimatable check is ",
                         is.null(nonEstimable(design))))

        }

    }
}

###OTHER FUNCTIONS USED BY normalize_raw_counts
vizualize_counts <- function(countsmatrix, norm_exprs, labels,
                             figres=100, results_path) {
    #Generate figures for counts
    print("STATUS: generating log2 boxplot of counts")
    generate_boxplots(log2(countsmatrix + 1), labels[, 1],
                      file.path(results_path, "1.boxplot_raw_count_matrix.png"),
                      figres, maintitle = "Raw count matrix",
                      ylabtitle = "log2 Expression")

    print("STATUS: generating boxplot of normalized voom counts")
    generate_boxplots(norm_exprs, labels[, 1],
                      file.path(results_path, "1.boxplot_vnorm_matrix.png"),
                      figres, maintitle = "Normalized count matrix",
                      ylabtitle = "voom normalized expression")

    print("STATUS: generating density plot of all log counts")
    generate_density_plot(log2(countsmatrix + 1), labels[, 1],
                         file.path(results_path,
                                   "1.densities_raw_log_counts.png"),
                         figres)

    print("STATUS: generating density plot of raw counts")
    generate_density_plot(countsmatrix, labels[, 1],
                          file.path(results_path, "1.densities_raw_counts.png"),
                          figres)

    print("STATUS: generating density plot of normalized voom counts")
    generate_density_plot(norm_exprs, labels[, 1],
                          file.path(results_path,
                                    "1.densities_vnorm_matrix.png"),
                          figres)


    # print("STATUS: generating biological varation vs abundance")
    # png(file.path(results_path, "1.biologicalcoefficentvariation_raw.png"),
    #      res = figres)
    # # par(mar=c(1,1,1,1))
    # plotBCV(count_matrix_flv, cex = 0.4,
    #         main = "Biological coefficient of variation (BCV) vs abundance")
    # dev.off()
}

# Biological coefficients of variation

biological_coefficents_variation <- function(count_matrix_fl) {
    count_matrix_flv <- estimateCommonDisp(count_matrix_fl,
                                           verbose = T) #print the BCV value
    count_matrix_flv <- estimateTrendedDisp(count_matrix_flv)
    count_matrix_flv <- estimateTagwiseDisp(count_matrix_flv)
    return(count_matrix_flv)
}

generate_density_plot <- function(data, labels, filename, figres) {
    png(filename, res = figres)
    par( xpd = TRUE)
    if (length(labels) > 10) {
        plotDensities(data, legend = FALSE)
    } else {
        plotDensities(data, legend = "topright",
                      inset = c(-0.2,0), levels(labels))
    }
    dev.off()
}



generate_boxplots <- function(data, labels, filename, figres, maintitle, ylabtitle) { 
    png(filename, res=figres)
    # par(mar=c(1,1,1,1))
    minvalue <- min(data)
    maxvalue <- max(data)
    boxplot(data, labels = labels, ylim = c(minvalue - 1, maxvalue + 1),
            ylab = ylabtitle, main = maintitle, cex.axis = .6, las = 2,
            frame = FALSE)
    dev.off()

}