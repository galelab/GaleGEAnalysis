
# Gale Lab differential gene analysis pipeline

Goal of this software is to provide a pipeline for lab members to analyze gene expression data


## Install 
First load library('devtools'). If you do not have devtools run install.packages('devtools). To download in R please use install_github('galelab/GaleGEAnalysis', dependencies=TRUE).  

## How to execute pipeline
GaleGEAnalysis currently has 5 primary functions
1. p1_modify_count_matrix - function will remove samples specified by the user from the count matrix (generally produced by HTSeq) and target file.  The target file is a csv with information about each sample (treatment vs no treatment etc).  This step is not necessary and can be skipped. Results in folder p1_modified_count_matrix_results/ in working directory. Samples to remove must be specified in a csv or txt file.
2. s1_normalize_raw_counts - function will normalize counts using limma and the trimmed mean of M-values (TMM)  normalization method.  The main inputs needed are the count file from htseq, the target file, the column in the target file specifying whether the sample under went treatment or non treatment (target_column) Can put column name or number for this parameter.  Other options are provided here as to whether to filter out genes that have low counts across all samples and to account for a batch effect (batch_column).  All results will be outputed to in a folder named s1_norm_raw_count_results/ in your working directory.
3. s2_feature_reduction - function will perform a variety of feature reduction alorithms on the normalized count data. As a default as long as you are in the same folder that you ran the s1_normalize_raw_counts step the default count file is ./s1_norm_raw_count_results/1.norm_matrix.txt.  If this file is located elsewhere you will have to specify its location in the countfile parameter.  Additionally, if you didn't run p1_modify_count_matrix you will have to specify the lcoation of the target file with the targetfile parameter in this funtion.  Finally, the user must specify the columns in the target file to label samples by in the pca and umap reduction plots with the target_columns parameter (must have 2 columns in vector format (i.e c(2,5))). Output will go to folder s2_feature_reduction in your working directory. (Note: this is not a ncessary step and can be skipped)
4. s3_DE_analysis - function will perform the DE analysis.  Like in step s2_feature_reduction specify count and target files in the same manner.  Additionaly user will need to generate matrix file which tells the pipeline what samples to do the DE analysis on.  See example MATRIXEXAMPLE.txt tells which treatments vs non treatments to compare to each other. Output files from this step will be stored in folder s3_DE_analysis_results/. IMPORTANT OUTPUT FILES:

    3.All_LFC.txt - contains log fold changes for all genes for each comparison in the matrix file 

    3.All_LFC_HGNC.txt - contains log fold changes for all genes for each comparison in the matrix file with HGNC gene labels 

    3.All_Pvalues.txt - contains P values for all genes for each comparison in the matrix file 

    3.All_Pvalues_HGNC.txt - contains P values for all genes for each comparison in the matrix file with HGNC gene labels 

    3.All_tvalues.txt - contains t values for all genes for each comparison in the matrix file 

    3.All_tvalues_HGNC.txt - contains t values for all genes for each comparison in the matrix file with HGNC gene labels 

    3.Significant_separate_LFC.csv - contains log fold for significant differentially expressed genes at atleast in one comparison 

    3.Significant_separate_LFC_HGNC_AV.csv - contains log fold for significant differentially expressed genes at atleast in one comparison with HGNC gene labels (AV means that genes with the same HGNC labels are averaged)

    3.Significant_separate_Pvalues.csv - contains p values for significant differentially expressed genes at atleast in one comparison 

    3.Significant_separate_Pvalues_HGNC_AV.csv - contains p values for significant differentially expressed genes at atleast in one comparison with HGNC gene labels (AV means that genes with the same HGNC labels are averaged)

    3.Significant_separate_tvalues.csv - contains t values for significant differentially expressed genes at atleast in one comparison 

    3.Significant_separate_tvalues_HGNC_AV.csv - contains t values for significant differentially expressed genes at atleast in one comparison with HGNC gene labels (AV means that genes with the same HGNC labels are averaged)

    3.heatmap_djn.png - heatmap of log fold change values of genese that are significantly differentially expressed at atleast one comparison.
    enrichfiles/ - this folder contains list of rank files that are used in step s4_gene_enrichment_analysis

        *_sig.rnk - contain list of significant differnetially expressed genes and log fold change values at given comparison to be used for GSEA analysis in other programs such as WebGestalt

        *_all4GSEA.rnk - contain list of all genes and log fold change values at given comparison to be used for GSEA analysis in other programs such as WebGestalt (Note: No duplicate genes in this list)    

5. s4_gene_enrichment_analysis - function will perform over representation analysis (ORA) on genes significantly expressed identified in step s3_DE_anlysis.  Step 3 clusters these genes according to expression patterns (see heatmap generated in step s3_DE_anlysis) into modules.  ORA is done on each module.  If user wants to look at ORA of a specific comparison specify the comparison in the comparison option of the s4_gene_enrichment_analysis function.  Also the user can input their own rank file and perform ORA on that instead of gene lists in step s3_DE_anlysis (see examples below). For each comparison an output folder will be generated with a number of output files such as:

    over_enrich_barplotall_ge.png - bar plot of significantly enriched GO terms using all differentially expressed genes

    over_enrich_barplotup_ge.png - bar plot of significantly enriched GO terms using up regulated differentially expressed genes

    over_enrich_barplotdown_ge.png - bar plot of significantly enriched GO terms using down regulated differentially expressed genes

    over_enrich_cneplotall_ge.png - cnet plot (network) of significantly enriched GO terms using all differentially expressed genes

    over_enrich_cneplotup_ge.png - cnet plot (network) of significantly enriched GO terms using up regulated differentially expressed genes

    over_enrich_cneplotdown_ge.png - cnet plot (network) of significantly enriched GO terms using down regulated differentially expressed genes

    *all_ora_genes.csv tables of genes and corresponding log fold changes in GO terms when all differentially expresed genes were used for analysis

    *up_ora_genes.csv tables of genes and corresponding log fold changes in GO terms when up regulated differentially expresed genes were used for analysis

    *down_ora_genes.csv tables of genes and corresponding log fold changes in GO terms when down regulated differentially expresed genes were used for analysis

## Example of pipeline run
Commands for running pipeline

1. p1_modify_count_matrix(countfile='count_matrix.txt', targetfile='target.csv', samples_to_remove_count_matrix='samplestoremove.csv')

2. s1_normalize_raw_counts(countfile='./p1_modified_count_matrix_results/count_matrix_mod.txt', targetfile='./p1_modified_count_matrix_results/targets_mod.csv', gene_conversion_file='rhesus2human.csv', target_column=3, batch_column=2,filter_genes_below_counts=50)

3. s2_feature_reduction(countfile='./s1_norm_raw_counts_results/1.norm_matrix.txt', targetfile='./p1_modified_count_matrix_results/targets_mod.csv', target_columns=c(2,5))

4. s3_DE_analysis(countfile='./s1_norm_raw_counts_results/1.norm_matrix.txt', targetfile='./p1_modified_count_matrix_results/target_file.csv', gene_conversion_file='rhesus2human.csv',  'matrixfile=MATRIXEXAMPLE.txt', pvalue=0.05, logfoldchange=1.5)

* This shows couple examples of how s4_gene_enrichment_analysis could be run
5. s4_gene_enrichment_analysis(go_enrich_type='BP', universe=TRUE, modules=TRUE, NumTopGoTerms=30) - This will perform ORA on genes classified in specific modules from step s3_DE_analysis ONLY as well as use all expressed genes as background (collected from ./s1_norm_raw_counts_results/1.norm_matrix_HGNC.txt) instead of whole genome (because universe is specified as true) and create tables and figures including the top 30 enriched go terms.

5. s4_gene_enrichment_analysis(go_enrich_type='BP', comparison='treatmentInfected_w19-treatmentUninfected_w19' universe=TRUE, modules=TRUE, NumTopGoTerms=30) - This will perform ORA on genes classified in specific modules from step s3_DE_analysis and genes significantly differentially expressed between the comparison 'treatmentInfected_w19-treatmentUninfected_w19' as well as use all expressed genes as background (collected from ./s1_norm_raw_counts_results/1.norm_matrix_HGNC.txt) instead of whole genome (because universe is specified as true) and create tables and figures including the top 30 enriched go terms.

5. s4_gene_enrichment_analysis(go_enrich_type='BP', comparison='treatmentInfected_w19-treatmentUninfected_w19' universe=FALSE, modules=FALSE, NumTopGoTerms=30) - This will perform ORA on genes classified in specific modules from step s3_DE_analysis and genes significantly differentially expressed between the comparison 'treatmentInfected_w19-treatmentUninfected_w19' as well as use whole genome (because universe is specified as false) as background genes and create tables and figures including the top 30 enriched go terms.

## Other important notes 
1. matrix file must be layed out in the exact same manner in as the example file (MATRIXEXAMPLE.txt).  The word 'treatment' must be in front of each comparison as in the example file.
2. different conversion file can be used but it must be in csv format layed out in the same manner as rhesus2human.csv
3. target file must be layed out in a similar manner as the example file (targetexample.csv) but amount of information in the target file will vary with each experiment