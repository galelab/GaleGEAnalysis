
# Gale Lab Differential gene analysis

Goal of this software is to provide a pipeline for lab members to analyze gene expression data


## Install 
To download in R please use install_github('lwhitmore/GaleGEAnalysis').  First load library('devtools')

## How to execute pipeline
GaleGEAnalysis currently has 4 primary functions
1. p1_modify_count_matrix - function will remove samples specified by the user from the count matrix (generally produced by HTSeq) and target file.  The target file is a csv with information about each sample (treatment vs no treatment etc).  This step is not necessary and can be skipped. Results in folder p1_modified_count_matrix_results/ in working directory.
2. s1_normalize_raw_counts - function will normalize counts using limma and the trimmed mean of M-values (TMM)  normalization method.  The main inputs needed are the count file from htseq, the target file, the column in the target file specifying whether the sample under went treatment or non treatment (target_column).  Other options are provided here as to whether to filter out genes that have low counts across all samples and to account for a batch effect (batch_column).  All results will be outputed to in a folder named s1_norm_raw_count_results/ in your working directory.
3. s2_feature_reduction - function will perform a variety of feature reduction alorithms on the normalized count data. As a default as long as you are in the same folder that you ran the s1_normalize_raw_counts step the default count file is ./s1_norm_raw_count_results/1.norm_matrix.txt.  If this file is located elsewhere you will have to specify its location in the countfile parameter.  Additionally, if you didn't run p1_modify_count_matrix you will have to specify the lcoation of the target file with the targetfile parameter in this funtion.  Finally, the user must specify the columns in the target file to label samples by in the pca and umap reduction plots with the target_columns parameter (must have 2 columns in vector format (i.e c(2,5))). Output will go to folder s2_feature_reduction in your working directory.
4. s3_DE_analysis - function will perform the DE analysis.  Like in step s2_feature_reduction specify count and target files in the same manner.  Additionaly user will need to generate matrix file which tells the pipeline what samples to do the DE analysis on.  See example MATRIXEXAMPLE.txt tells which treatments vs non treatments to compare to each other. Output files from this step will be stored in folder s3_DE_analysis_results/. Values in files 3.ExpressMatrix_separate_LFC_HGNC_AV.csv and 3.ExpressMatrix_separate_LFC.csv are log fold change.  
5. s4_gene_enrichment_analysis - function will perform over representation analysis (ORA) and gene set enrichment analysis (GSEA) on list of genes provide by the user or from a list of differentially expressed genes from step s3_DE_analysis. The column with the log fold change the user wants to use to run GSEA or ORA on must be specified by parameter log_values_column as well as the column where the gene symbols are by parameter gene_name_column.  Note that gene names must be HGNC symbol and must be human.  This means previous steps must have used the gene_conversion_parameter unless gene IDs were orginally HGNCs.  Defaults for this may not work with users analysis so be careful.  Results are put in folder s4_gene_enrichment_results/ unless otherwise specified by user.

## Example of pipeline run
p1_modify_count_matrix(countfile='count_matrix.txt', targetfile='target.csv', samples_to_remove_count_matrix=c(sample1,sample2))

s1_normalize_raw_counts(countfile=./p1_modified_count_matrix_results/count_matrix_mod.txt, targetfile='./p1_modified_count_matrix_results/targets_mod.csv', gene_conversion_file=rhesus2human.csv, target_column=3, batch_column=2,filter_genes_below_counts=50)

s2_feature_reduction(countfile='./s1_norm_raw_counts_results/1.norm_matrix.txt', targetfile='./p1_modified_count_matrix_results/targets_mod.csv', target_columns=c(2,5))

s3_DE_analysis(countfile='./s1_norm_raw_counts_results/1.norm_matrix.txt', targetfile='./p1_modified_count_matrix_results/target_file.csv', gene_conversion_file=rhesus2human.csv,  'matrixfile=MATRIXEXAMPLE.txt', pvalue=0.05, logfoldchange=1.5)

s4_gene_enrichment_analysis(DEgenes='./s3_DE_results/3.ExpressMatrix_separate_LFC_HGNC_AV.csv',  go_enrich_type='BP')