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

s4_gene_enrichment_analysis<-function(DEgenes='./s3_DE_results/3.ExpressMatrix_separate_LFC_HGNC_AV.csv', log_values_column=7, pvalue=0.1, NumTopGoTerms=10, figres=100) {
    results_path  <- generate_folder('s4_gene_enrichment_results')
    results_path  <- generate_folder(paste0('s4_gene_enrichment_results/column',log_values_column))
    unlink(paste0(results_path,'/*'))

    genefile      <- read.csv(DEgenes)
    
    genernk       <- genefile[,c(3,log_values_column)]

    write.table(genernk, sep='\t',quote = FALSE, col.names=FALSE, row.names=FALSE, file=file.path(results_path,'DEgenes.rnk'))

    genes         <- read.table(file.path(results_path,'DEgenes.rnk'), sep='\t')
    geneList      <- genes[,2]
    names(geneList) = as.character(genes[,1])
    x             <-names(geneList)

    topGOdata     <- new("topGOdata", description="topGO BP", ontology="BP",nodeSize=5, allGenes=geneList, geneSel=function(p) p>pvalue, annot=annFUN.org, mapping="org.Hs.eg.db", ID="SYMBOL")
    resultFisher  <- runTest(topGOdata, algorithm="classic", statistic="fisher")
    resultKS      <- runTest(topGOdata, algorithm ="classic", statistic='ks')
    allRes        <- GenTable(topGOdata, classicFisher=resultFisher, classicKS=resultKS, orderBy ="classicFisher", topNodes=NumTopGoTerms)
    write.table(allRes, file=file.path(results_path, 'GOTermTable.csv'),sep=',', row.names=FALSE)
    # vizualize_enrichment(file.path(results_path,'GoEnrichmentFig.pdf'), allRes, figres=figres)
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
    # ggsave(plot_file,height=6.75,width=9)
    # gp<-ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
    #     stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    #     xlab("Biological process") +
    #     ylab("Enrichment") +
    #     ggtitle("Title") +
    #     scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
    #     theme_bw(base_size=24) +
    #     theme(
    #         legend.position='none',
    #         legend.background=element_rect(),
    #         plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
    #         axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
    #         axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
    #         axis.title=element_text(size=24, face="bold"),
    #         legend.key=element_blank(),     #removes the border
    #         legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    #         legend.text=element_text(size=18),  #Text size
    #         title=element_text(size=18)) +
    #     guides(colour=guide_legend(override.aes=list(size=2.5))) +
    #     coord_flip()
    pdf(plot_file)
    barplot(goEnrichment$classicKS, names.arg=goEnrichment$Term)
    dev.off()
}