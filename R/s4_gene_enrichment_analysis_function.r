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

s4_gene_enrichment_analysis<-function(DEgenes='./s3_DE_results/3.ExpressMatrix_separate_LFC_HGNC_AV.csv', log_values_column=7, pvalue=0.1, NumTopGoTerms=10) {
    results_path <- generate_folder('s4_gene_enrichment_results')
    results_path <- generate_folder(paste0('s4_gene_enrichment_results/column',log_values_column))
    unlink(paste0(results_path,'/*'))

    genefile <- read.csv(DEgenes)
    
    genernk <- genefile[,c(3,log_values_column)]

    write.table(genernk, sep='\t',quote = FALSE, col.names=FALSE, row.names=FALSE, file=file.path(results_path,'DEgenes.rnk'))

    genes <- read.table(file.path(results_path,'DEgenes.rnk'), sep='\t')
    geneList <- genes[,2]
    names(geneList) = as.character(genes[,1])
    x<-names(geneList)

    topGOdata     <- new("topGOdata", description="topGO BP", ontology="BP",nodeSize=5, allGenes=geneList, geneSel=function(p) p>pvalue, annot=annFUN.org, mapping="org.Hs.eg.db", ID="SYMBOL")
    resultFisher  <- runTest(topGOdata, algorithm="classic", statistic="fisher")
    resultKS <- runTest(topGOdata, algorithm ="classic", statistic='ks')
    allRes        <- GenTable(topGOdata, classicFisher=resultFisher, classicKS=resultKS, orderBy ="classicFisher", topNodes=NumTopGoTerms)
    write.table(allRes, file=file.path(results_path, 'GOTermTable.csv'),sep=',', row.names=FALSE)
    goIDs          <- allRes$GO.ID
    
    for (gt in goIDs) {
        print (paste0('STATUS: getting genes in ', gt, ' term'))
        term.genes <- genesInTerm(topGOdata, gt)
        affID <- term.genes[[gt]]
        pval <- sort(geneScore(topGOdata, affID, use.names = TRUE))
        # pval <- pval[pval <= pvalue]
        affID <- names(pval)
        tablex<-select(org.Hs.eg.db, affID, c("GENENAME"), "ALIAS")
        for (i in 1:length(pval)) {
            tablex[['logFC']][i]<-pval[i]
        }
        write.table(tablex, file=file.path(results_path, paste0(gt,'_genes.csv')), sep=',', row.names=FALSE)
    }

    # WebGestaltR(enrichMethod="GSEA",
    #         organism="hsapiens",
    #         interestGeneFile="./s4_gene_enrichment_results/DEgenes.rnk",
    #         minNum=5,
    #         interestGeneType="genesymbol",
    #         enrichDatabase = c("geneontology_Biological_Process","geneontology_Cellular_Component","geneontology_Molecular_Function"),
    #         sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8, outputDirectory='./s4_gene_enrichment_results'
    #         )

    # x<-gseGO(geneList='s4_gene_enrichment_results/DEgenes.rnk', pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb=org.Hs.eg.db, keyType = "SYMBOL")
}