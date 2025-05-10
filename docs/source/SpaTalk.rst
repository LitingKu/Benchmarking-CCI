SpaTalk Analysis of CCI
=========================================== 

The following R script demonstrates how we use SpaTalk conduct the CCI analysis.

- The SpaTalk documentation and package could be found here: https://github.com/ZJUFanLab/SpaTalk

- The paper could be found here: https://doi.org/10.1038/s41467-022-32111-8


1. Dependencies
-------------------------

.. code-block:: r


   library(SpaTalk)
   library(NNLM)
   library(dplyr)



2. Run the analysis
-------------------------

- We follow the SpaTalk guidance and tutorial for the analysis.

.. code-block:: r

   folder_names <- c("slice1")
   for (num in folder_names) {
     print(paste0("doing dataset: ", num))
     path <- paste0("/rsrch5/home/biostatistics/lku/ILIBD/data/", num, "/exprsn_df.csv")
     data <- read.csv(path)
     data <- data[!is.na(data$cell_type), ]
     index <- which(colnames(data) == "cell_type")
     gene_spot_df <- data[,c((index+1) : ncol(data))]
     rownames(gene_spot_df) <- data$cell
     meta_data <- data[,c("cell","x","y","cell_type")]
     meta_data$cell <- paste0("cell",c(1:nrow(meta_data)))
     rownames(gene_spot_df) <- meta_data$cell
     


     ###### Option 1: Doing the deconvolution #######################
     # create spatalk object
     obj <- createSpaTalk(st_data = as.matrix(t(gene_spot_df)),
                          st_meta = meta_data[, -4],
                          species = "Human",
                          if_st_is_sc = F,
                          spot_max_cell = 10)
     
     # scRNA reference data: include the sc_count and sc_meta
     load("/rsrch5/home/biostatistics/lku/ILIBD/data/sc_ref_data.rda")

     obj <- dec_celltype(object = obj,
                         sc_data = as.matrix(sc_count),
                         sc_celltype = sc_meta$cell_type)
  
     ###### Option 2: If single cell level ST, direct assign cell type ##
     obj <- createSpaTalk(st_data = as.matrix(t(gene_spot_df)),
                          st_meta = meta_data[, -4],
                          species = "Human",
                          if_st_is_sc = T,
                          spot_max_cell = 1,
                          celltype = meta_data$cell_type)
     
     # Find possible L-R database
     obj <- find_lr_path(object = obj, lrpairs = lrpairs, pathways = pathways)
     obj <- dec_cci_all(object = obj, pvalue =1.1)
     df <- obj@lrpair
     df$interaction_name <- paste0(df$ligand,"_",df$receptor)
     final_table <- df[,c("celltype_sender","celltype_receiver","score","lr_co_ratio_pvalue","interaction_name")]
     write.csv(final_table,paste0("/rsrch5/home/biostatistics/lku/ILIBD/data/",num,"/spatalk_result.csv"))

   }




