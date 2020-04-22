#' s4 gene enrichment analysis
#'
#' This function runs gene enrichment analysis
#'
#' @keywords gene enrichment analysis
#' @param go_enrich_type type of GO enrichment to do (BP, CC, or MF) (default BP)
#' @param universe Will specify whether you want to use the whole genome as background (FALSE) or set of total expressed genes(TRUE)
#' @param universe_file specifies path to universe file(defualt ./s1_raw_counts_results/1.norm_matrix_HGNC.txt)
#' @param DEgenes file containg log fold change for DE genes (default FALSE, needs to be given with log_values_column if used)
#' @param log_fold_column column to pull from DEgenes
#' @param rnkfile file containing list of users own genes (not generated in step s3).  Must be in rank file format (gene name & logfold value seperated by tab)
#' @param results_folder user specified output folder (default is s4_gene_enrichment_results)
#' @param comparison specify name of time point comparison to perform over enrichment analysis 
#' @param modules perform over enrichment analysis on modules generated in step s3
#' @param module_results_folder folder where DE results are (default s3_DE_results/)
#' @param NumTopGoTerms top GO terms to show (default 10)
#' @param figres resolution of output figures (default 300)
#' @param ensembl_retrieve whether or not to retrieve gene name descriptions
#' @param base_file_name name to save files under (default ge.png)
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import AnnotationHub
#' @import biomaRt
#' @import data.table
#' @import ggplot2
#' @import stringr
#' @export
#' @examples
#' s4_gene_enrichment_analysis(DEgenes=FALSE, go_enrich_type="BP", log_values_column=FALSE, modules=TRUE, pvalue=0.05, qvalue=0.05, NumTopGoTerms=30)

s4_gene_enrichment_analysis <- function(go_enrich_type="BP", universe=TRUE,
                                        universe_file =" ./s1_raw_counts_results/1.norm_matrix_HGNC.txt",
                                        DEgenes=FALSE, log_values_column=FALSE,
                                        rnkfile=FALSE, results_folder='s4_gene_enrichment_results',
                                        comparison=FALSE, modules=TRUE,
                                        module_results_folder="s3_DE_results",
                                        NumTopGoTerms=30, figres=300,
                                        ensembl_retrieve=TRUE,
                                        gene_name_type="SYMBOL",
                                        base_file_name="ge.png") {
    results_path <- generate_folder(results_folder)
    unlink(paste0(results_path, "/*"))
    if (typeof(comparison) == "character") {
        results_path  <- generate_folder(paste0(results_folder,
                                                "/", comparison))
        unlink(paste0(results_path, "/*"))
    } else if (typeof(log_values_column) != "logical") {
        results_path  <- generate_folder(paste0(results_folder,
                                                "/column",
                                                log_values_column))
        unlink(paste0(results_path, "/*"))
    }
    if (modules == TRUE) {
        results_path_mod  <- generate_folder(paste0(results_folder,
                                                    "/modules"))
        unlink(paste0(results_path_mod, "/*"))
    }

    if (isTRUE(universe)) {
        all_genes_table <- read.table(universe_file,
                                      row.names = 1, sep = "\t")
        all_universe_genes <- all_genes_table$HGNC.symbol
        all_universe_genes <- as.character(all_universe_genes)

    } else {
        print("STATUS: Will be using whole genome as background in over enrichment analysis")
    }
    ###Connect to biomart
    if (ensembl_retrieve == TRUE) {
        ensembl <- useEnsembl(biomart = "ensembl", mirror="uswest",
                              dataset = "hsapiens_gene_ensembl")
    } else {
        ensembl <- FALSE
    }
    ###Perform over enrichment on modules
    ###generating from cluster heatmap in step s3
    if (modules == TRUE) {
        module <- read.csv(paste0(module_results_folder, "/3.modules_HGNC.csv"), row.names = 1)
        unique_modules <- unique(module$V1)

        for (m in unique_modules) {
            mod <-  module[(module[, 3] == m),]
            mod <- mod[!(is.na(mod$HGNC.symbol) | mod$HGNC.symbol == ""), ]
            geneList <- cbind(mod$HGNC.symbol, c(rep(1, length(mod$HGNC.symbol))))
            head(geneList)
            rownames(geneList) <- geneList[, 1]
            head(geneList)
            geneList <- geneList[, -1]
            head(geneList)
            class(geneList) <- "numeric"
            # print(geneList)
            if (isTRUE(universe)) {
                ego <- run_over_enrichment(mod$HGNC.symbol,
                                           go_enrich_type = go_enrich_type,
                                           gene_name_type = gene_name_type,
                                           universe = all_universe_genes)
            } else {
                ego <- run_over_enrichment(mod$HGNC.symbol,
                                          go_enrich_type = go_enrich_type,
                                          gene_name_type = gene_name_type)
            }

            write.csv(as.data.frame(ego), file = file.path(results_path_mod,
                                                 paste0("total_enrichment_",
                                                         m, ".csv")))

            if (is.null(ego)) {
                print ("WARNNING: No enrichments found so no output figures or files will be generated")
            } else {
                barplot(ego, showCategory = NumTopGoTerms)
                ggsave(file.path(results_path_mod,
                       paste0("over_enrich_barplot_", m, "_", base_file_name)),
                        width = 10, height = 8, dpi = figres)
                cnetplot(ego, foldChange = NULL)
                ggsave(file.path(
                    results_path_mod,
                    paste0("over_enrich_cnetplot_", m, "_", base_file_name)
                    ), width = 6, height = 6, dpi = 150
                )
                if (ensembl_retrieve == TRUE) {
                    extract_genesego(ego, rnk = FALSE,
                                     results_path_mod, ensembl,
                                     enrich_type = "ora", direction = m)
                }
            }
        }
    }

    if ((typeof(comparison) != "logical") | (typeof(log_values_column) != "logical") | (typeof(rnkfile) != "logical")) {
        if (typeof(comparison) == "character") {
            genes <- read.table(file.path(paste0(module_results_folder, "enrichfiles",
                                          comparison, "_sig.rnk")),
                                sep = "\t")
            genesall <- read.table(file.path(paste0(module_results_folder, "/enrichfiles",
                                             comparison, "_all.rnk")),
                                   sep = "\t")

        } else if (typeof(log_values_column) != "logical") {
            genefile <- read.csv(DEgenes)
            cnames <- colnames(genefile)
            genes <- read.table(file.path(paste0(module_results_folder, "/enrichfiles",
                                cnames[log_values_column], "_sig.rnk")),
                                sep = "\t")
        } else if (typeof(rnkfile) != "logical") {
            genes <- read.table(rnkfile, sep = "\t")
        }

        gene_up <-  genes[(genes[, 2] > 0), ]
        gene_down <-  genes[(genes[, 2] < 0), ]

        geneList_up <- gene_up[, 2]
        names(geneList_up) <- as.character(gene_up[, 1])
        geneList_up <- sort(geneList_up, decreasing = TRUE)
        genenames_up <- names(geneList_up)

        geneList_down <- gene_down[, 2]
        names(geneList_down) <- as.character(gene_down[, 1])
        geneList_down <- sort(geneList_down, decreasing = TRUE)
        genenames_down <- names(geneList_down)
        geneList <- genes[, 2]

        geneList <- genes[, 2]
        names(geneList) <- as.character(genes[, 1])
        genenames <- names(geneList)
        geneList <- sort(geneList, decreasing = TRUE)



        if (isTRUE(universe)) {
            print ("STATUS: running over enrichment analysis with expressed genes as background")
            egoup <- run_over_enrichment(genenames_up,
                                         go_enrich_type = go_enrich_type,
                                         universe = all_universe_genes,
                                         gene_name_type = gene_name_type)
            egodown <- run_over_enrichment(genenames_down,
                                           go_enrich_type = go_enrich_type,
                                           universe = all_universe_genes,
                                           gene_name_type = gene_name_type)
            ego <- run_over_enrichment(genenames,
                                       go_enrich_type = go_enrich_type,
                                       universe = all_universe_genes,
                                      gene_name_type = gene_name_type)

        } else {
            print("STATUS: running over enrichment analysis with whole genome as background")
            egoup  <- run_over_enrichment(genenames_up,
                                          go_enrich_type = go_enrich_type,
                                         gene_name_type = gene_name_type)
            egodown <- run_over_enrichment(genenames_down,
                                           go_enrich_type = go_enrich_type,
                                           gene_name_type = gene_name_type)
            ego <- run_over_enrichment(genenames,
                                       go_enrich_type = go_enrich_type,
                                       gene_name_type = gene_name_type)

        }
        ego_dim <- dim(as.data.frame(ego))
        if (ego_dim[1] != 0) {
            cnetplot(ego, foldChange = NULL)
            ggsave(file.path(results_path,
                             paste0("over_enrich_cnetplotall_",
                                    base_file_name)), dpi = figres)
            barplot(ego, showCategory = NumTopGoTerms)
            ggsave(file.path(results_path, paste0("over_enrich_barplotall_",
                                                 base_file_name)),
                             width = 10, height = 8, dpi = figres)
            if (ensembl_retrieve == TRUE) {
                extract_genesego(ego, genes, results_path, ensembl,
                                 enrich_type = "ora", direction = "all")
            }
            write.csv(as.data.frame(ego), file = file.path(results_path,
                                          paste0("total_enrichment_all.csv")))
        } else {
            print ("WARNING: not enough data to generate cnetplot for all differentially expressed genes over enrichment")
        }

        egoup_dim <- dim(as.data.frame(egoup))
        if (egoup_dim[1] != 0) {
            cnetplot(egoup, foldChange = NULL)
            ggsave(file.path(results_path, paste0("over_enrich_cnetplotup_",
                                                  base_file_name)),
                   dpi = figres)
            barplot(egoup, showCategory = NumTopGoTerms)
            ggsave(file.path(results_path, paste0("over_enrich_barplotup_",
                                                  base_file_name)),
                   width = 10, height = 8, dpi = figres)
            if (ensembl_retrieve == TRUE) {
                extract_genesego(egoup, genes, results_path, ensembl,
                                 enrich_type = "ora", direction = "up")
            }
            write.table(as.data.frame(egoup), file = file.path(results_path,
                                            paste0("total_enrichment_up.csv")))
        } else {
            print ("WARNING: not enough data to generate cnetplot for upregulated genes over enrichment")
        }

        egodown_dim <- dim(as.data.frame(egodown))
        if (egodown_dim[1] != 0) {
            cnetplot(egodown, foldChange = NULL)
            ggsave(file.path(results_path, paste0("over_enrich_cnetplotdown_",
                                                  base_file_name)),
                    dpi = figres)
            barplot(egodown, showCategory = NumTopGoTerms)
            ggsave(file.path(results_path, paste0("over_enrich_barplotdown_",
                                                  base_file_name)),
                    width = 10, height = 8, dpi = figres)
            if (ensembl_retrieve == TRUE) {
                extract_genesego(egodown, genes, results_path, ensembl,
                                 enrich_type = "ora", direction = "down")
            }
            write.csv(as.data.frame(egodown), file = file.path(results_path,
                                       paste0("total_enrichment_down.csv")))
        } else {
            print("WARNING: not enough data to generate cnetplot for downregulated genes over enrichment")
        }
    }
}

run_over_enrichment <- function(genes,  go_enrich_type, gene_name_type,
                                universe=FALSE) {
    # print (go_enrich_type)
    if (typeof(universe) != "logical") {
        ego <- enrichGO(genes, universe = universe,
                OrgDb         = org.Hs.eg.db,
                keyType       = gene_name_type,
                ont           = go_enrich_type,
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1)

    } else {
        ego <- enrichGO(genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = gene_name_type,
                ont           = go_enrich_type,
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1)
    }
    return(ego)

}


extract_genesego <- function(enrichment, rnk, results_path, ensembl,
                             enrich_type="ora", direction="up", sig_cutoff=0.1) {


    enrichment <- enrichment[enrichment$p.adjust <= sig_cutoff, ]

    for (i in 1:dim(enrichment)[1]) {
        description <- str_replace(enrichment$Description[i], "/", "_")
        print (paste0("STATUS: getting genes in ", description, " term"))
        genes <- unlist(strsplit(enrichment$geneID[i], "/"))

        if (length(genes)>0) {
            genedesc <- getBM(attributes = c("external_gene_name", "description"),
                            filters = "external_gene_name",
                            values = genes, mart = ensembl)

            ###Remove unecessary extra information from description
            count <- 1
            for (i in genedesc$description) {
                i <- gsub("\\s+\\[Source\\S+\\s+\\S+$", "", i, perl = TRUE)
                genedesc$description[count] <- i
                count <- count + 1
            }

            if (typeof(rnk) != "logical") {
                tabl <- data.table(genedesc)
                tabl1 <- data.table(rnk)

                genedescfinal <- merge(tabl, tabl1,
                                    by.x = "external_gene_name",
                                    by.y = "V1")
                write.table(genedescfinal, file = file.path(results_path,
                                                    paste0(description,
                                                        "_", direction,
                                                        "_", enrich_type,
                                                        "_genes.csv")),
                            sep = ",", row.names = FALSE)
            } else {
                tabl <- data.table(genedesc)
                write.table(tabl, file = file.path(results_path,
                                        paste0(description, "_",
                                                direction, "_",
                                                enrich_type,
                                                "_genes.csv")),
                            sep = ",", row.names = FALSE)
            }
        }
    }
}
