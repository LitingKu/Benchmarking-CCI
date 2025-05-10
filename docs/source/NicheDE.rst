NicheDE Analysis of CCI
=========================================== 

The following R script demonstrates how we use NicheDE conduct the CCI analysis.

- The NicheDE documentation and package could be found here: https://github.com/kaishumason/NicheDE

- The paper could be found here: https://doi.org/10.1186/s13059-023-03159-6


1. Dependencies
-------------------------

.. code-block:: r


   Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
   options(timeout=9999999)
   #devtools::install_github("kaishumason/NicheDE") # install
   library(nicheDE)
   library(dplyr)

2. Load the default L-R database from the package
--------------------------------------------------

.. code-block:: r

   
   data("niche_net_ligand_target_matrix")
   data("ramilowski_ligand_receptor_list")


3. Prepare the function and scRNA-seq reference for `CreateLibraryMatrix` 
---------------------------------------------------------------------------

We extracted and modified this function because the default implementation provided in the original package resulted in an error. The corrected version is shown below:

.. code-block:: r

   # scRNA reference data: include the sc_count and sc_meta
   load("/rsrch5/home/biostatistics/lku/ILIBD/data/sc_ref_data.rda")

   CreateLibraryMatrix = function(data, cell_type) {
     if (mean(rownames(data) == cell_type[, 1]) != 1) {
       stop('data rownames and Cell type matrix names do not match')
     }

     CT = unique(as.vector(cell_type[, 2]))
     n_CT = length(CT)
     L = matrix(NA, n_CT, ncol(data))
     rownames(L) = CT
     colnames(L) = colnames(data)

     for (j in c(1:n_CT)) {
       cells = which(cell_type[, 2] == CT[j])
       if (length(cells) > 1000) {
         print(paste0("Too many cell of type ", CT[j], " downsampling to 1000."))
         cells = sample(cells, 1000, replace = FALSE)
       } else if (length(cells) == 1) {
         cells = data[cells, ]
         L[j, ] = mean(cells)
       }
       cells = data[cells, ]
       L[j, ] = apply(cells, 2, mean)
     }

     print('Average expression matrix computed')
     return(L)
   }






4. Run the analysis and generate output table for downstream analysis
---------------------------------------------------------------------------

.. code-block:: r

   folder_names <- c("slice1")
   for (num in folder_names) {
     print(paste0("doing dataset: ", num))
     path <- paste0("/rsrch5/home/biostatistics/lku/ILIBD/data/", num, "/exprsn_df.csv")
     data <- read.csv(path)
     data <- data[!is.na(data$cell_type), ]
     index <- which(colnames(data) == "cell_type")

     # Load counts matrix
     gene_spot_df <- data[, c((index + 1):ncol(data))]
     rownames(gene_spot_df) <- data$cell
     gene_spot_df <- as.matrix(gene_spot_df)

     # Load coordinate matrix
     coord_data <- data[, c("y", "x")]
     colnames(coord_data) <- c("imagerow", "imagecol")
     rownames(coord_data) <- rownames(gene_spot_df)

     # Extract cell proportions
     cell_index <- which(colnames(data) == "cell")
     celltype_index <- which(colnames(data) == "cell_type")
     cell_prop_df <- as.matrix(data[, c((cell_index + 1):(celltype_index - 1))])
     rownames(cell_prop_df) <- rownames(gene_spot_df)


     # Load expression profile matrix using CreateLibraryMatrix function
     our <- CreateLibraryMatrix(t(sc_count), sc_meta[, c("TAG", "cell_type")])
     expression_profile <- our[colnames(cell_prop_df), ]

     # Construct NicheDE object and run analysis
     NDE_obj <- CreateNicheDEObject(gene_spot_df, coord_data, expression_profile, cell_prop_df, sigma = c(1, 100, 250))
     NDE_obj <- CalculateEffectiveNiche(NDE_obj)
     NDE_obj <- niche_DE(NDE_obj, num_cores = 4, outfile = "", C = 50, M = 10, gamma = 0.8, print = TRUE, Int = TRUE, batch = TRUE, self_EN = FALSE, G = 1)

     # Run niche_LR analysis for all cell type pairs
     print(paste0("- - - - - - - - -doing niche_LR"))
     LR_results <- list()
     celltypes <- colnames(cell_prop_df)

     for (type in celltypes) {
       for (type2 in celltypes) {
         tryCatch({
           LR <- niche_LR_spot(
             NDE_obj,
             ligand_cell = type,
             receptor_cell = type2,
             ligand_target_matrix = niche_net_ligand_target_matrix,
             lr_mat = ramilowski_ligand_receptor_list,
             K = 25,
             M = 50,
             alpha = 0.05,
             truncation_value = 3
           )

           result_name <- paste(type, type2, sep = "_")
           LR <- as.data.frame(LR)
           LR$ligand_celltype <- type
           LR$receptor_celltype <- type2
           LR$interaction_name <- paste0(LR$ligand, "_", LR$receptor)
           LR_results[[result_name]] <- LR
         }, error = function(e) {
           cat("Error in ligand_cell:", type, "receptor_cell:", type2, "\n")
         })
       }
     }

     # Combine and save final result
     print(paste0("- - - - - - - - - - - - - - writing dataframe"))
     final_df <- do.call(rbind, LR_results)
     write.csv(final_df, paste0("/rsrch5/home/biostatistics/lku/ILIBD/data/", num, "/nichede_result_new.csv"))
   }




