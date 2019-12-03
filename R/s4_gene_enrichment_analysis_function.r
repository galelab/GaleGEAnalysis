#' s4 gene enrichment analysis
#'
#' This function runs gene enrichment analysis
#' 
#' @keywords gene enrichment analysis 
#' @import clusterProfiler 
#' @import org.Hs.eg.db
#' @import AnnotationHub
#' @import biomaRt
#' @import data.table
#' @import ggplot2
#' @export
#' @examples
#' 

s4_gene_enrichment_analysis <-function(DEgenes='./s3_DE_results/3.ExpressMatrix_separate_LFC_HGNC_AV.csv', go_enrich_type='BP', result_folder=FALSE, gene_name_column=3, log_values_column=7, pvalue=0.1, NumTopGoTerms=10, figres=100, base_file_name='ge.png') {
    if (typeof(result_folder) == 'logical') {
        results_path  <- generate_folder('s4_gene_enrichment_results')
        results_path  <- generate_folder(paste0('s4_gene_enrichment_results/column',log_values_column))
        unlink(paste0(results_path,'/*'))
    } else { 
        results_path  <- generate_folder(result_folder)
        results_path  <- generate_folder(paste0(result_folder,'/column',log_values_column))
        unlink(paste0(results_path,'/*'))
    }


    genefile      <- read.csv(DEgenes)
    genernk       <- genefile[,c(gene_name_column, log_values_column)]
    write.table(genernk, sep='\t',quote = FALSE, col.names=FALSE, row.names=FALSE, file=file.path(results_path,'DEgenes.rnk'))

    genes         <- read.table(file.path(results_path,'DEgenes.rnk'), sep='\t')
    geneList      <- genes[,2]
    names(geneList) = as.character(genes[,1])
    x             <-names(geneList)
    geneList      <- sort(geneList,decreasing=TRUE)


    ego <- enrichGO(gene          = x,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'SYMBOL',
                    ont           = go_enrich_type,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = pvalue,
                    qvalueCutoff  = 0.05)

    gse   <- gseGO(geneList=geneList, ont=go_enrich_type, pvalueCutoff  = pvalue, keyType='SYMBOL', OrgDb=org.Hs.eg.db, verbose=F)   
    
    cnetplot(ego, foldChange=geneList)
    ggsave(file.path(results_path, paste0('over_enrich_cnetplot_',base_file_name)), dpi=300)

    barplot(ego, showCategory=NumTopGoTerms)
    ggsave(file.path(results_path, paste0('over_enrich_barplot_',base_file_name)), dpi=300)

    dotplot(gse,showCategory=NumTopGoTerms)
    ggsave(file.path(results_path, paste0('gse_enrich_doplot_',base_file_name)), dpi=300)

    extract_genesego(ego, genes, results_path, enrich_type='ora', NumGOterms=NumTopGoTerms)
    extract_genes(gse, genes, results_path, enrich_type='gsea', NumGOterms=NumTopGoTerms)
}

extract_genes <- function(enrichment, rnk, results_path, enrich_type='ora', NumGOterms=10) {

    ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")    
    
    for(i in 1:NumGOterms) {
    
        print (paste0('STATUS: getting genes in ', enrichment[i]$Description, ' term'))
        genes <- unlist(strsplit(enrichment[i]$core_enrichment, '/'))

        genedesc <- getBM(attributes=c('external_gene_name', 'description'), filters = 'external_gene_name', values = genes, mart =ensembl)
        tabl <- data.table(genedesc)
        tabl1 <- data.table(rnk)

        genedescfinal <- merge(tabl, tabl1, by.x='external_gene_name', by.y='V1')
        write.table(genedescfinal, file=file.path(results_path, paste0(enrichment[i]$Description,'_', enrich_type, '_genes.csv')), sep=',', row.names=FALSE)
     }
}

extract_genesego <- function(enrichment, rnk, results_path, enrich_type='ora', NumGOterms=10) {

    ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")    
    
    for(i in 1:NumGOterms) {
    
        print (paste0('STATUS: getting genes in ', enrichment[i]$Description, ' term'))
        genes <- unlist(strsplit(enrichment[i]$geneID, '/'))

        genedesc <- getBM(attributes=c('external_gene_name', 'description'), filters = 'external_gene_name', values = genes, mart =ensembl)
        tabl <- data.table(genedesc)
        tabl1 <- data.table(rnk)

        genedescfinal <- merge(tabl, tabl1, by.x='external_gene_name', by.y='V1')
        write.table(genedescfinal, file=file.path(results_path, paste0(enrichment[i]$Description,'_', enrich_type, '_genes.csv')), sep=',', row.names=FALSE)
     }
}
