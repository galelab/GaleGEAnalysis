#' p1 modifying count matrix modifies count matrix 
#'
#' This function allows remove samples from count matrix if they are not to be included in downstream gene expression analysis
#' @param countfile raw counts table (generally output by htseq).
#' @param targetfile target file.
#' @param samples_to_remove_count_matrix text file (no header) or r list of samples to remove (defualt is FALSE)
#' @keywords remove samples
#' @export
#' @examples
#' p1_modify_count_matrix(count_file.txt, target_file.csv, samples_to_remove_count_matrix=c(sample1,sample2))

p1_modify_count_matrix <- function(countfile, targetfile, samples_to_remove_count_matrix=FALSE) {
    #remove samples to not include analyses
    if (samples_to_remove_count_matrix != FALSE) {
        if ((grep('.txt', samples_to_remove_count_matrix) == 1) || (grep('.csv', samples_to_remove_count_matrix) == 1)) {
            samples_to_remove <- read.table(samples_to_remove_count_matrix, header=FALSE, row.names=1)
            samples_to_remove_count_matrix <- rownames(samples_to_remove)
        }
        if (length(samples_to_remove_count_matrix) > 0) {

            files <- loadfiles(countfile, targetfile)
            # print (dim(files$counts))
            # print (dim(files$targets))

            modcounts <- files$counts[ , !(names(files$counts) %in% samples_to_remove_count_matrix)]
            modtargets <- files$targets[!(rownames(files$targets) %in% samples_to_remove_count_matrix), ]
            # print(dim(modcounts))
            # print(dim(modtargets))

            print(paste0('STATUS: The number of samples removed were ', length(samples_to_remove_count_matrix)))
            results_path <- generate_folder('p1_modified_count_matrix_results')

            write.table(data.frame("Name"=rownames(modcounts), modcounts), sep = "\t", row.names=FALSE, file=file.path(results_path,"count_matrix_mod.txt"))
            write.csv(data.frame("Name"=rownames(modtargets), modtargets), row.names=FALSE, file=file.path(results_path,"targets_mod.csv"))
            
            results <- list("counts" = modcounts, "targets" = modtargets)
            
            return (results)
        }
    }
    print ('STATUS: NO SAMPLES REMOVED')
}


loadfiles <- function(count_file, targets_file) {
    #Load in count data  and target information from csv'
    counts <- read.table(count_file, header = TRUE, sep = "\t", row.names=1, as.is = TRUE)
    targets <- read.table(targets_file, header = TRUE, sep = ",", row.names = 1, as.is = TRUE)
    results <- list("counts" = counts, "targets" = targets)
    return(results)
}

generate_folder <- function(foldername) {
    workDir <-getwd()
    subDir = foldername
    results_path = file.path(workDir, subDir)
    if (file.exists(subDir)){
    } else { 
        dir.create(results_path)
    }
    return (results_path)
}