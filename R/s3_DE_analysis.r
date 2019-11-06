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
#' @import gplots 
#' @examples
#' 

s3_DE_analysis <- function(countfile, targetfile,  gene_conversion_file=FALSE, target_class=10, blocking_column=FALSE, matrixfile=FALSE, p_value=0.05, log_fold_change=1.5) { 

    results_path <- generate_folder('s3_DE_results')
    unlink('./s3_DE_results/*')

    files <- loadfiles(count_file=countfile, target_file=targetfile)

    if (file.exists('./s1_norm_raw_counts_results/1.designobject.rds')) { 
        design     <- readRDS('./s1_norm_raw_counts_results/1.designobject.rds')
    } else {
        print ('WARNING: Could not find design file...make sure you are in correct working directory or repeat step (s1_normalize_raw_counts)')
    }
    if (file.exists('./s1_norm_raw_counts_results/1.voomobject.rds')) {
        eset_voom  <- readRDS('./s1_norm_raw_counts_results/1.voomobject.rds')
    } else { 
        print ('WARNING: Could not find voom file...make sure you are in correct working directory or repeat step (s1_normalize_raw_counts)')
    }
    if (file.exists('./s1_norm_raw_counts_results/1.corfit.rds')) {
        corfit     <- readRDS('./s1_norm_raw_counts_results/1.corfit.rds')
    } else { 
        print ('WARNING: Could not find corfit file...not necessary for analysis but just warning the user')
    }


    if (exists("design") == FALSE | exists("eset_voom") == FALSE) { 
        print ('WARNING: could not find design and voom model... make sure user is in the correct working directory')
    } else { 
        if (typeof(blocking_column) != 'logical') {
            BLOCKID  <- factor(files$targets[,blocking_column], levels=unique(files$targets[,blocking_column]))

            if (exists('corfit') == TRUE) {
                fit <- lmFit(eset_voom,design,block=BLOCKID, correlation=corfit$consensus)
            } else { 
                fit <- lmFit(eset_voom,design,block=BLOCKID)
            }
        } else { 
            if (exists('corfit') == TRUE) {
                fit <- lmFit(eset_voom,design, correlation=corfit$consensus)
            } else { 
                fit <- lmFit(eset_voom,design)
            }
        }
        if (typeof(matrixfile) == 'character') {
            matrix_contrast <- scan(matrixfile,  character(), quote = "")

            cont.matrix <- makeContrasts(contrasts=matrix_contrast, levels=design)

            fit <- contrasts.fit(fit, cont.matrix)
            fit <-eBayes(fit)

            results <- decideTests(fit, lfc=log2(log_fold_change), method="separate", adjust.method="BH", p.value=p_value)
            a <- vennCounts(results)

            write.fit(fit, file=file.path(results_path, "3.DE_orig_fit.txt"), digits=3, method="separate", adjust="BH")
            DE_HGNC <- read.csv(file.path(results_path, "3.DE_orig_fit.txt"), header = T,row.names = 1, check.names=FALSE,sep = "\t")
            if (typeof(gene_conversion_file) == 'character') {
                rhesus2human <- read.csv(file=gene_conversion_file, header=TRUE, stringsAsFactors = FALSE)
                DE_HGNC <- merge(rhesus2human, DE_HGNC, by.x='Gene.stable.ID', by.y='row.names')
                DE_HGNC <- avereps(DE_HGNC, ID = DE_HGNC$HGNC.symbol)
                write.csv(DE_HGNC, file=file.path(results_path, "3.DE_orig_fit_HGNC.csv"))
            }

            #Pull out signifcantly expressed genes 
            dataMatrix <- fit$coefficients # Extract results of differential expression #LogFold change is the coefficients
            sigMask <- dataMatrix * (results**2) # 1 if significant, 0 otherwise
            ExpressMatrix <- subset(dataMatrix, rowSums(sigMask) != 0) # filter for significant genes
            sigMask <- subset(sigMask, rowSums(sigMask) != 0)
            write.csv(ExpressMatrix, file=file.path(results_path,"3.ExpressMatrix_separate_LFC.csv"))
            global_modules <- vizualize_DE_genes_HM(ExpressMatrix, file.path(results_path, "3.heatmap_djn.png"))
            write.csv(global_modules, file=file.path(results_path, "3.modules.csv"))

            if (typeof(gene_conversion_file) == 'character') {
                DE_HGNC_LFC <- read.csv(file.path(results_path, "3.ExpressMatrix_separate_LFC.csv"), header = T,row.names = 1, check.names=FALSE, sep = ",")
                DE_HGNC_LFC <- merge(rhesus2human, DE_HGNC_LFC, by.x='Gene.stable.ID', by.y='row.names')
                DE_HGNC_LFC <- avereps(DE_HGNC_LFC, ID = DE_HGNC_LFC$HGNC.symbol)
                write.csv(DE_HGNC_LFC, file=file.path(results_path, "3.ExpressMatrix_separate_LFC_HGNC_AV.csv"))
                global_modulesM <- as.matrix(global_modules)
                GM_HGNC <- merge(rhesus2human, global_modulesM, by.x='Gene.stable.ID', by.y='row.names',all.X=T,all.Y=T)
                write.csv(GM_HGNC, file=file.path(results_path, "3.modules_HGNC.csv"))
            }
        } else { 
            print ('WARNING: need to specify matrix file')
        }
    }
}

vizualize_DE_genes_HM<-function(data, plot_file) {
    png(plot_file, width = 8, height = 10, units = 'in', res = 300)
    global_modules <- heatmap.F.4(data, cutoff = 1, distmethod = "euclidean", clustermethod = "ward.D", clusterdim='row')
    dev.off()
    return (global_modules)
}