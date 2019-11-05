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
        print (' 1')
        if (typeof(blocking_column) != 'logical') {
            BLOCKID  <- factor(files$targets[,blocking_column], levels=unique(files$targets[,blocking_column]))
            print (' 2')
            if (exists('corfit') == TRUE) {
                fit <- lmFit(eset_voom,design,block=BLOCKID, correlation=corfit$consensus)
            } else { 
                fit <- lmFit(eset_voom,design,block=BLOCKID)
            }
        } else { 
            print (' 3')
            if (exists('corfit') == TRUE) {
                fit <- lmFit(eset_voom,design, correlation=corfit$consensus)
            } else { 
                fit <- lmFit(eset_voom,design)
            }
        }
        if (typeof(matrixfile) == 'character') {
            print (' 4')
            matrix_contrast <- scan(matrixfile,  character(), quote = "")
            print (matrix_contrast)
            print (' 5')
            # cont.matrix <- makeContrasts(contrasts=matrix_contrast, levels=design)
            cont.matrix <- makeContrasts(
            (treatmentInfected_BL - treatmentUninfected_BL),
            (treatmentInfected_w1 - treatmentUninfected_w1),
            (treatmentInfected_w18 - treatmentUninfected_w18),
            (treatmentInfected_w19 - treatmentUninfected_w19),
            (treatmentInfected_w25 - treatmentUninfected_w25),
            (treatmentInfected_w26 - treatmentUninfected_w26),
            (treatmentInfected_w27 - treatmentUninfected_w27),
            (treatmentInfected_w41 - treatmentUninfected_w41),
            levels=design)
            print (' 6')
            fit <- contrasts.fit(fit, cont.matrix)
            fit <-eBayes(fit)
            results <- decideTests(fit, lfc=log2(log_fold_change), method="separate", adjust.method="BH", p.value=p_value)
            summary(results)
            a <- vennCounts(results)
            print(a)
        } else { 
            print ('WARNING: need to specify matrix file')
        }
    }
}

