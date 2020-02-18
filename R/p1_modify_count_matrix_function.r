#' p1 modifying count matrix modifies count matrix
#'
#' This function allows remove samples from count matrix if they
#'  are not to be included in downstream gene expression analysis
#' @param countfile raw counts table (generally output by htseq).
#' @param targetfile target file.
#' @param samples_to_remove_count_matrix text file (no header) or r list of samples to remove (defualt is FALSE)
#' @keywords remove samples
#' @export
#' @examples
#' p1_modify_count_matrix(count_file.txt, target_file.csv, samples_to_remove_count_matrix=c(sample1,sample2))

p1_modify_count_matrix <- function(countfile, targetfile,
                                   samples_to_remove_count_matrix=FALSE) {
    #remove samples to not include analyses
    if (typeof(samples_to_remove_count_matrix) == "character") {
        if ((grep(".txt", samples_to_remove_count_matrix) == 1) ||
            (grep(".csv", samples_to_remove_count_matrix) == 1)) {
            samples_to_remove <- read.table(samples_to_remove_count_matrix,
                                            header = FALSE, row.names = 1,
                                            check.names = FALSE)
            samples_to_remove_count_matrix <- rownames(samples_to_remove)
        }
        if (length(samples_to_remove_count_matrix) > 0) {

            files <- loadfiles(countfile, targetfile)

            #Remove unecessaary samples from count matrix
            modcounts <- files$counts[ , !(names(files$counts) %in% samples_to_remove_count_matrix)]
            #Sort column names in count matrix
            modcounts <- modcounts[, order(names(modcounts))]

            #Remove unecessaary samples from target file
            modtargets <- files$targets[!(rownames(files$targets) %in% samples_to_remove_count_matrix), ]
            ##Sorts targets by row names so hopefully it matches ordered column rows above in count matrix
            modtargets <- modtargets[order(rownames(modtargets)), ]


            print(paste0("STATUS: The number of samples removed were ",
                         length(samples_to_remove_count_matrix)))
            results_path <- generate_folder("p1_modified_count_matrix_results")

            write.table(data.frame("Name" = rownames(modcounts),
                                   modcounts, check.names = FALSE),
                        sep = "\t", col.names = TRUE,  row.names = FALSE,
                        file = file.path(results_path, "count_matrix_mod.txt"))

            write.csv(data.frame("Name" = rownames(modtargets),
                                 modtargets, check.names = FALSE),
                     row.names = FALSE,
                    file = file.path(results_path, "targets_mod.csv"))

            results <- list("counts" = modcounts, "targets" = modtargets)
        }
    } else {
       print("STATUS: NO SAMPLES REMOVED")
    }
}
