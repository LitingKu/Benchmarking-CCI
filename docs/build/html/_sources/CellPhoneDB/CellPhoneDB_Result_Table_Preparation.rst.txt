CellPhoneDB Result Table Preparation
=========================================== 

The following R script reads CellPhoneDB statistical result files (means and p-values) and reshapes the data for downstream analysis.

1. Dependencies
-------------------------

.. code-block:: r


   library(tidyr)
   library(dplyr)

2. Create a function to find and read the text file into a data frame
---------------------------------------------------------------------------

.. code-block:: r

   
   read_statistical_file <- function(base_dir, pattern) {
     file_list <- list.files(path = base_dir, pattern = pattern, full.names = TRUE)
     if (length(file_list) == 1) {
       file_path <- file_list[1]
       tryCatch({
         df <- read.table(file_path, header = TRUE, sep = "\t")
         cat("Successfully read", file_path, "\n")
         return(df)
       }, error = function(e) {
         cat("Error reading", file_path, ":", e$message, "\n")
         return(NULL)
       })
     } else if (length(file_list) == 0) {
       cat("No files matching the pattern were found in", base_dir, "\n")
     } else {
       cat("Multiple files matching the pattern were found in", base_dir, "\n")
     }
   }


3. Generating table for downstream analysis
--------------------------------------------------

.. code-block:: r

   base_dir <- paste0("/Volumes/lku/ILIBD/data/slice1")

   means <- read_statistical_file(base_dir, "statistical_analysis_means")
   means <- means[which(means$directionality %in% c("Ligand-Receptor")),]
   ind <- which(colnames(means)=="classification")
   coln <- c("interacting_pair","partner_a","partner_b","gene_a","gene_b", colnames(means)[(ind+1):length(colnames(means))])
   dat <- means[, c(coln)]
   dat_means <- dat %>% gather(key = "detect", value = "value", -c(interacting_pair, partner_a, partner_b, gene_a, gene_b))
   dat_means <- dat_means[which(dat_means$value > 0.0), ]
   colnames(dat_means) <- c("interacting_pair", "partner_a", "partner_b", "gene_a", "gene_b", "detect", "means")

   pval <- read_statistical_file(base_dir, "statistical_analysis_pvalues")
   pval <- pval[which(pval$directionality %in% c("Ligand-Receptor")),]
   ind <- which(colnames(pval)=="classification")
   coln <- c("interacting_pair", colnames(pval)[(ind+1):length(colnames(pval))])
   dat <- pval[, c(coln)]
   dat_pval <- dat %>% gather(key = "detect", value = "value", -interacting_pair)

   df <- left_join(dat_means, dat_pval, by=c("interacting_pair", "detect"))
   df <- df[!is.na(df$value),]

   df$ligand <- sapply(strsplit(df$detect, "\\."), `[`, 1)
   df$receptor <- sapply(strsplit(df$detect, "\\."), `[`, 2)

   df <- df[,c("ligand","receptor","means","value","interacting_pair","partner_a", "partner_b", "gene_a", "gene_b")]
   colnames(df) <- c("ligand","receptor","means","pvalue","interaction_name","partner_a", "partner_b", "gene_a", "gene_b")

   write.csv(df, paste0("/Volumes/lku/ILIBD/data/slice1/cellphonedb_result.csv"))








