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

s3_DE_analysis <- function(countfile, targetfile, target_class=10, blocking_column=FALSE, matrixfile=FALSE, p_value=0.05, log_fold_change=1.5) { 
    # design = FALSE
    # eset_voom = FALSE
    # corfit = FALSE
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
            print (summary(results))
            a <- vennCounts(results)
            # print(a)
            write.fit(fit, file="SG_Watkins_DE.txt", digits=3, method="separate", adjust="BH")
            DE_HGNC <- read.csv("SG_Watkins_DE.txt", header = T,row.names = 1, check.names=FALSE,sep = "\t")
            rhesus2human <- read.csv(file="rhesus2human.csv", header=TRUE, stringsAsFactors = FALSE)
            print (' 1')
            DE_HGNC <- merge(rhesus2human, DE_HGNC, by.x='Gene.stable.ID', by.y='row.names')
            print (' 2')
            DE_HGNC <- avereps(DE_HGNC, ID = DE_HGNC$HGNC.symbol)
            # print (' 3')
            # print (DE_HGNC)
            # row.names(DE_HGNC) <- unique(DE_HGNC[,1])
            # print (' 4')
            # write.csv(DE_HGNC, file="DEgenes_HGNC.csv")
            dataMatrix <- fit$coefficients # Extract results of differential expression

            sigMask <- dataMatrix * (results**2) # 1 if significant, 0 otherwise

            ExpressMatrix <- subset(dataMatrix, rowSums(sigMask) != 0) # filter for significant genes

            # Filter sigMask to use for selecting DE genes from ExpressMatrix
            sigMask <- subset(sigMask, rowSums(sigMask) != 0)

            dim(sigMask)

            length(sigMask)

            write.csv(ExpressMatrix, file="ExpressMatrix_separate_Watkins.csv")
            # source('heatmap2.F.R')
            png("Watkins_heatmap_djn.png",width = 8, height = 10, units = 'in', res = 300)
            global_modules <- heatmap.F.4(ExpressMatrix, cutoff = 1, distmethod = "euclidean", clustermethod = "ward.D", clusterdim='row')
            dev.off()
 
        } else { 
            print ('WARNING: need to specify matrix file')
        }
    }
}

