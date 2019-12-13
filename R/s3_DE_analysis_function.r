#' s3 differential gene analysis 
#'
#' This function runs differential gene analysis (linear modeling on normalized counts)
#' @param countfile normalized counts table (generally should have been generated by first step s1_normalize_raw_counts, but also this function can be run on raw cout file).
#' @param targetfile target file.
#' @param gene_conversion_file file with alternative gene names (we usually convert Rhesus Ensembl genes to HGNC) 
#' @param blocking_column column to account sampling from the same animal multiple times (needs to be the same as whatever was specified in step s1_normalize_raw_counts)
#' @param matrixfile text file outlining different contrasts to do in DE analysis 
#' @param pvalue parameter for determining significantly DE genes (default 0.05)
#' @param logfoldchange parameter for determining significantly DE genes (default 1.5)
#' @keywords differential expression linear modeling 
#' @export
#' @import edgeR
#' @import limma
#' @import gplots 
#' @import ggplot2
#' @import data.table
#' @examples
#' s3_DE_analysis(countfile='./s1_norm_raw_counts_results/1.norm_matrix.txt', targetfile='./p1_modified_count_matrix_results/target_file.csv', gene_conversion_file='rhesus2human.csv', blocking_column=2, matrixfile='./MATRIX.txt')
 

s3_DE_analysis <- function(countfile='./s1_norm_raw_counts_results/1.norm_matrix.txt', targetfile='./p1_modified_count_matrix_results/target_file.csv',  gene_conversion_file=FALSE, blocking_column=FALSE, matrixfile=FALSE, pvalue=0.05, logfoldchange=1.5) { 

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
            print ('STATUS: getting DE genes...')
 
            matrix_contrast <- scan(matrixfile,  character(), quote = "")

            cont.matrix <- makeContrasts(contrasts=matrix_contrast, levels=design)

            fit <- contrasts.fit(fit, cont.matrix)
            fit <-eBayes(fit)

            results <- decideTests(fit, lfc=log2(logfoldchange), method="separate", adjust.method="BH", p.value=pvalue)
            a <- vennCounts(results)

            write.fit(fit, file=file.path(results_path, "3.All_LFC.txt"), digits=3, method="separate", adjust="BH")
            DE_HGNC <- read.csv(file.path(results_path, "3.All_LFC.txt"), header = T,row.names = 1, check.names=FALSE,sep = "\t")
            if (typeof(gene_conversion_file) == 'character') {
                rhesus2human <- read.csv(file=gene_conversion_file, header=TRUE, stringsAsFactors = FALSE)
                DE_HGNC <- merge(rhesus2human, DE_HGNC, by.x='Gene.stable.ID', by.y='row.names')
                DE_HGNC <- avereps(DE_HGNC, ID = DE_HGNC$HGNC.symbol)
                write.csv(DE_HGNC, file=file.path(results_path, "3.All_LFC_HGNC.csv"))
            }

            write.table(fit$p.value, file=file.path(results_path, "3.All_Pvalues.txt"), sep='\t')#, digits=3, method="separate", adjust="BH")
            DE_HGNC <- read.csv(file.path(results_path, "3.All_Pvalues.txt"), header = T,row.names = 1, check.names=FALSE,sep = "\t")
            if (typeof(gene_conversion_file) == 'character') {
                rhesus2human <- read.csv(file=gene_conversion_file, header=TRUE, stringsAsFactors = FALSE)
                DE_HGNC <- merge(rhesus2human, DE_HGNC, by.x='Gene.stable.ID', by.y='row.names')
                DE_HGNC <- avereps(DE_HGNC, ID = DE_HGNC$HGNC.symbol)
                write.csv(DE_HGNC, file=file.path(results_path, "3.All_Pvalues_HGNC.csv"))
            }            
            write.table(fit$t, file=file.path(results_path, "3.All_tvalues.txt"), sep='\t')#, digits=3, method="separate", adjust="BH")
            DE_HGNC <- read.csv(file.path(results_path, "3.All_tvalues.txt"), header = T,row.names = 1, check.names=FALSE,sep = "\t")
            if (typeof(gene_conversion_file) == 'character') {
                rhesus2human <- read.csv(file=gene_conversion_file, header=TRUE, stringsAsFactors = FALSE)
                DE_HGNC <- merge(rhesus2human, DE_HGNC, by.x='Gene.stable.ID', by.y='row.names')
                DE_HGNC <- avereps(DE_HGNC, ID = DE_HGNC$HGNC.symbol)
                write.csv(DE_HGNC, file=file.path(results_path, "3.All_tvalues_HGNC.csv"))
            }            


            #Pull out signifcantly expressed genes

            ##SIGNIFICANT LOGVALUES
            dataMatrix <- fit$coefficients # Extract results of differential expression #LogFold change is the coefficients
            sigMask <- dataMatrix * (results**2) # 1 if significant, 0 otherwise
            ExpressMatrixLFC <- subset(dataMatrix, rowSums(sigMask) != 0) # filter for significant genes
            sigMask <- subset(sigMask, rowSums(sigMask) != 0)
            write.csv(ExpressMatrixLFC, file=file.path(results_path,"3.Significant_separate_LFC.csv"))
            convert2HGNC(gene_conversion_file, "3.Significant_separate_LFC.csv", "3.Significant_separate_LFC_HGNC_AV.csv", results_path)

            ##SIGNIFICANT T values
            dataMatrix <- fit$t # Extract results of differential expression #LogFold change is the coefficients
            sigMask <- dataMatrix * (results**2) # 1 if significant, 0 otherwise
            ExpressMatrixtvalue <- subset(dataMatrix, rowSums(sigMask) != 0) # filter for significant genes
            sigMask <- subset(sigMask, rowSums(sigMask) != 0)
            write.csv(ExpressMatrixtvalue, file=file.path(results_path,"3.Significant_separate_tvalues.csv"))            
            convert2HGNC(gene_conversion_file, "3.Significant_separate_tvalues.csv", "3.Significant_separate_tvalues_HGNC_AV.csv", results_path)
    
             ##SIGNIFICANT P values
            dataMatrix <- fit$p.value # Extract results of differential expression #LogFold change is the coefficients
            sigMask <- dataMatrix * (results**2) # 1 if significant, 0 otherwise
            ExpressMatrixPvalue <- subset(dataMatrix, rowSums(sigMask) != 0) # filter for significant genes
            sigMask <- subset(sigMask, rowSums(sigMask) != 0)
            write.csv(ExpressMatrixPvalue, file=file.path(results_path,"3.Significant_separate_Pvalues.csv"))         
            convert2HGNC(gene_conversion_file, "3.Significant_separate_Pvalues.csv", "3.Significant_separate_Pvalues_HGNC_AV.csv", results_path)            
 
            results_path2 <- generate_folder('s3_DE_results/enrichfiles')
            unlink('./s3_DE_results/enrichfiles/*')

            counter <- 0
            for (i in colnames(fit$coefficients)) {
                counter <- counter + 1
                results_single <-  results[, c(i)]
                results_single <- results_single[(results_single[,1] != 0),]
                significantgenes <- fit$coefficients[rownames(results_single), counter]
                allgenes         <- fit$coefficients[, counter]
                if (typeof(gene_conversion_file) == 'character') {
                    sig_HGNC <- merge(rhesus2human, significantgenes, by.x='Gene.stable.ID', by.y='row.names',all.X=T,all.Y=T)
                    sig_HGNC <- sig_HGNC[ , !(names(sig_HGNC) %in% c('Gene.stable.ID'))]
                    sig_HGNC <- avereps(sig_HGNC, ID = sig_HGNC$HGNC.symbol)

                    all_HGNC <- merge(rhesus2human, allgenes, by.x='Gene.stable.ID', by.y='row.names',all.X=T,all.Y=T)
                    all_HGNC <- all_HGNC[ , !(names(all_HGNC) %in% c('Gene.stable.ID'))]
                    all_HGNC <- avereps(all_HGNC, ID = all_HGNC$HGNC.symbol)
                    write.table(sig_HGNC, file=file.path(results_path2, paste0(i,'_sig.rnk')), row.names=FALSE, col.names = FALSE, sep='\t', quote=FALSE)
                    write.table(all_HGNC, file=file.path(results_path2, paste0(i,'_all.rnk')), row.names=FALSE, col.names = FALSE, sep='\t', quote=FALSE)

                    all_HGNCnodup           <- all_HGNC[!duplicated(all_HGNC[,2]), ]
                    sig_HGNCnodup           <- sig_HGNC[!duplicated(sig_HGNC[,2]), ]
                    write.table(sig_HGNCnodup, file=file.path(results_path2, paste0(i,'_sig4GSEA.rnk')), row.names=FALSE, col.names = FALSE, sep='\t', quote=FALSE)
                    write.table(all_HGNCnodup, file=file.path(results_path2, paste0(i,'_all4GSEA.rnk')), row.names=FALSE, col.names = FALSE, sep='\t', quote=FALSE)
                } else { 
                    write.table(significantgenes, file=file.path(results_path2, paste0(i,'_sig.rnk')), col.names = FALSE, sep='\t', quote=FALSE)
                    write.table(allgenes, file=file.path(results_path2, paste0(i,'_all.rnk')), col.names = FALSE, sep='\t', quote=FALSE)
                }   

                # datamatrix    <- fit$p.value[, c(i)]
                # sigMask       <- dataMatrix * (results_single**2) # 1 if significant, 0 otherwise
                # ExpressMatrix <- subset(dataMatrix, rowSums(sigMask) != 0) # filter for significant genes
                # sigMask       <- subset(sigMask, rowSums(sigMask) != 0
                # P.values      <-p.adjust(ExpressMatrix, method="BH");

            }


            hm_results      <- vizualize_DE_genes_HM(ExpressMatrixLFC, file.path(results_path, "3.heatmap_djn.png"))
            global_modules  <- hm_results$modules
            write.csv(global_modules, file=file.path(results_path, "3.modules.csv"))
            global_modulesM <- as.matrix(global_modules)
            GM_HGNC         <- merge(rhesus2human, global_modulesM, by.x='Gene.stable.ID', by.y='row.names',all.X=T,all.Y=T)
            write.csv(GM_HGNC, file=file.path(results_path,  "3.modules_HGNC.csv"))
            clustermatrix  <- hm_results$clustermatrix
            clustermatrix   <- clustermatrix[order(nrow(clustermatrix):1),] #invert row order
            write.csv(clustermatrix, file=file.path(results_path, "3.Clustered_LFC.csv"), quote=FALSE)
            # clustermatrixM <- as.matrix(clustermatrix)
            # clustermatrixM_hgnc    <- merge(clustermatrixM, rhesus2human, by.y='Gene.stable.ID', by.x='row.names',all.X=T,all.Y=T)#, all.X=T,all.Y=T)
            # clustermatrixM  <- as.data.frame(clustermatrix)
            # clustermatrixM$Gene.Stable.ID <- rownames(clustermatrixM)
            # rhesus2human    <- as.data.frame(rhesus2human)
            # print (head(rhesus2human$Gene.stable.ID))
            # print (head(clustermatrixM_hgnc))
            # clustermatrixM$HGNC.symbol <- rhesus2human$HGNC.symbol[match(rhesus2human$Gene.stable.ID, rownames(clustermatrixM))]
            # clustermatrixM_hgnc    <- inner_join(rhesus2human, clustermatrixM, by='Gene.Stable.ID')
            # print (head(clustermatrixM_hgnc))
            # write.csv(clustermatrixM_hgnc, file=file.path(results_path,  "3.Clustered_LFC_HGNC.csv"))     
            vizualize_DE_genes_bp(results, file.path(results_path,'3.barplot_NumDEgenes.png'))

        } else { 
            print ('WARNING: need to specify matrix file')
        }
    }
}

convert2HGNC<-function(gene_conversion_file, input_file, output_file, results_path) { 
    if (typeof(gene_conversion_file) == 'character') {
        rhesus2human <- read.csv(file=gene_conversion_file, header=TRUE, stringsAsFactors = FALSE)
        DE_HGNC_LFC <- read.csv(file.path(results_path, input_file), header = T,row.names = 1, check.names=FALSE, sep = ",")
        DE_HGNC_LFC <- merge(rhesus2human, DE_HGNC_LFC, by.x='Gene.stable.ID', by.y='row.names')
        DE_HGNC_LFC <- avereps(DE_HGNC_LFC, ID = DE_HGNC_LFC$HGNC.symbol)
        write.csv(DE_HGNC_LFC, file=file.path(results_path, output_file))

    } else { 
        print ('WARNING: need to specify conversion file to convert Ensembls to HGNCs')
    }
}

vizualize_DE_genes_HM<-function(data, plot_file) {
    print ('STATUS: Generating heatmap of DE genes...')
    png(plot_file, width = 8, height = 10, units = 'in', res = 300)
    global_modules <- heatmap.F.4(data, cutoff = 1, distmethod = "euclidean", clustermethod = "ward.D", clusterdim='row')
    dev.off()
    return (global_modules)
}

vizualize_DE_genes_bp<-function(results, plot_file) {
    print ('STATUS: Generating bar plot of number of DE genes...')
    results_t <- t(summary(results))
    results_t <- results_t[,-2]

    for (i in 1:(length(row.names(results_t)))) {
        results_t[i, 1] <- results_t[i, 1] * -1
    }

    DE <- as.data.frame(results_t)
    DE <- setnames(DE, old=c("Var1","Var2", "Freq"), new=c("Time_Point", "group", "DE_genes"))

    #write.csv(DE , file="barplot_TP186.csv",row.names=FALSE)

    #Create plot
    ggplot(DE, aes(x=Time_Point, y=DE_genes, fill=group, label = DE$DE_genes))+  geom_bar(stat="identity", position="identity")+
    # geom_text(size = 5, position = position_stack(vjust = 0) )+
    scale_fill_manual(values = c("#9d9dff", "#ff4d4d")) +
    ylab("Number of Differentially Expressed Genes") +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    ggsave(plot_file, dpi = 300)
}