#' s4 gene enrichment analysis
#'
#' This function runs gene enrichment analysis
#' 
#' @keywords gene enrichment analysis 
#' @import topGO 
#' @import org.Hs.eg.db
#' @import AnnotationHub
#' @export
#' @examples
#' 

s4_gene_enrichment_analysis<-function(DEgenes='./s3_DE_results/3.ExpressMatrix_separate_LFC_HGNC_AV.csv', result_folder=FALSE, gene_name_column=1, log_values_column=7, pvalue=0.1, NumTopGoTerms=10, figres=100) {
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

    topGOdata     <- new("topGOdata", description="topGO BP", ontology="BP",nodeSize=5, allGenes=geneList, geneSel=function(p) p>pvalue, annot=annFUN.org, mapping="org.Hs.eg.db", ID="SYMBOL")
    # resultFisher  <- runTest(topGOdata, algorithm="classic", statistic="fisher")
    resultKS      <- runTest(topGOdata, algorithm ="classic", statistic='ks')
    allRes        <- GenTable(topGOdata, classicKS=resultKS, orderBy ="classicFisher", topNodes=NumTopGoTerms)
    write.table(allRes, file=file.path(results_path, 'GOTermTable.csv'),sep=',', row.names=FALSE)
    vizualize_enrichment(file.path(results_path,'GoEnrichmentFig.png'), allRes, figres=figres)
    goIDs         <- allRes$GO.ID
    
     for (gt in goIDs) {
         print (paste0('STATUS: getting genes in ', gt, ' term'))
         term.genes <- genesInTerm(topGOdata, gt)
         affID      <- term.genes[[gt]]
         pval       <- sort(geneScore(topGOdata, affID, use.names = TRUE))
         # pval <- pval[pval <= pvalue]
         affID      <- names(pval)
         tablex     <-select(org.Hs.eg.db, affID, c("GENENAME"), "ALIAS")
         for (i in 1:length(pval)) {
             tablex[['logFC']][i]<-pval[i]
         }
         write.table(tablex, file=file.path(results_path, paste0(gt,'_genes.csv')), sep=',', row.names=FALSE)
     }
}

vizualize_enrichment<-function(plot_file, goEnrichment, figres) {
    ###NOTE ggplot FUNCTION WAS PULLED FROM THE WEB https://www.biostars.org/p/350710/
    pdf(NULL)
    print (goEnrichment)
    png(plot_file)
    par(mar = c(2, 18, 2, 2))
    barplot(goEnrichment$Annotated, names.arg=goEnrichment$Term, las=2, col='blue', horiz=TRUE)
    dev.off()
}