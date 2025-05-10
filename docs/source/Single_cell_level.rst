Single-cell-level Simulated Data
=========================================== 

The following R script demonstrates how we generate simulated single-cell-level ST data.

- We use the scRNA-seq to generate the ST data according to `GSM2230759 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi>`__, where the file could be downloaded from `Here <https://drive.google.com/file/d/1ME3fhDi4muYvDoSQWUVyJoGoZryV0M2k/view?usp=sharing>`__.


- We design that the single-cell-level ST data is well-clustered.



1. Dependencies
-------------------------

.. code-block:: r

   library(MASS)
   library(FNN)
   library(dplyr)
   library(ggplot2)
   library(deldir)
   library(grDevices)

2. Design the spatial coordinates and cell counts for each cell type based on scRNA reference
---------------------------------------------------------------------------------------------- 

.. code-block:: r

   # Read the scRNA-seq reference
   data3 <- read.csv("/rsrch5/home/biostatistics/lku/spot_level/test/GSM2230759_human3_umifm_counts.csv")

   cell_counts <- c(
     acinar = 843, # give it 7 cells 843+7
     activated_stellate = 100, 
     alpha = 1130, 
     beta = 787, 
     delta = 161, 
     ductal = 376, 
     endothelial = 92, 
     gamma = 36, 
     macrophage = 14, 
     quiescent_stellate = 54
   )
   # Total number of cells
   num_cells <- sum(cell_counts)

   # Generate spatial grid
   grid_size <- ceiling(sqrt(num_cells))


3. Voronoi tessellation to create the well-clustered ST data
--------------------------------------------------------------

- Here we generate different simulated spatial coordinates and save them into ``spatial_cluster_coords.csv``.

- Then we also visualize the well-clustered spatial organization by saving plots to ``spatial_cluster_coords.pdf``.

- We set random seed for generating random simulated data.


.. code-block:: r

   for (num in 1:50) {
     set.seed(num + 41)

     # Create grid points
     grid_points <- expand.grid(
       X = 1:grid_size,
       Y = 1:grid_size
     )

     # Seed points for Voronoi tessellation
     seed_points <- data.frame(
       X = sample(1:grid_size, length(cell_counts), replace = FALSE),
       Y = sample(1:grid_size, length(cell_counts), replace = FALSE),
       cell_type = names(cell_counts)
     )

     # Assign each grid point to the nearest seed (Voronoi region)
     grid_points <- grid_points %>%
       rowwise() %>%
       mutate(region = which.min((X - seed_points$X)^2 + (Y - seed_points$Y)^2)) %>%
       ungroup()

     # Add cell type labels
     grid_points <- grid_points %>%
       mutate(cell_type = seed_points$cell_type[region])

     # Sample desired number of cells per type
     final_points <- data.frame()
     remaining_grid <- grid_points

     for (cell_type in names(cell_counts)) {
       subset <- remaining_grid %>% filter(cell_type == cell_type)

       if (nrow(subset) < cell_counts[cell_type]) {
         stop(paste("Not enough points for cell type:", cell_type))
       }

       sampled_points <- subset %>% slice_sample(n = cell_counts[cell_type])
       final_points <- bind_rows(final_points, sampled_points)
       remaining_grid <- remaining_grid %>%
         filter(!paste(X, Y) %in% paste(sampled_points$X, sampled_points$Y))
     }

     if (nrow(final_points) != num_cells) {
       stop("Total number of points doesn't match 3593.")
     }

     if (any(duplicated(final_points[, c("X", "Y")]))) {
       stop("Duplicate coordinates found!")
     }

     specific_celltypes <- c("acinar", "activated_stellate", "alpha", "beta", "delta",
                             "ductal", "endothelial", "gamma", "macrophage", "quiescent_stellate")
     colors <- scPalette(length(specific_celltypes))
     names(colors) <- specific_celltypes

     p1 <- ggplot(final_points, aes(x = X, y = Y, color = cell_type)) +
       geom_point(size = 2) +
       scale_color_manual(values = colors) +
       labs(title = "", x = " ", y = " ", color = "Cell Type") +
       theme_minimal()

     path1 <- paste0("/rsrch5/home/biostatistics/lku/simdata/new/", num, "/spatial_cluster_coords.pdf")
     ggsave(path1, plot = p1, width = 6, height = 5)

     final_points <- final_points[, c("X", "Y", "cell_type")]
     path2 <- paste0("/rsrch5/home/biostatistics/lku/simdata/new/", num, "/spatial_cluster_coords.csv")
     write.csv(final_points, path2)
   }


4. Pre-assign the expression counts into each single cell
-------------------------------------------------------------

- We read the previous ``spatial_cluster_coords.csv``.

- Then we derive the gene expression distribution from a scRNA-seq reference and remodel it based on the designed cell type composition.

- Finally, we assign the expression counts into each single cell.

.. code-block:: r


   scRNA <- data3[which(data3$assigned_cluster %in% c(
     "acinar", "activated_stellate", "alpha", "beta", "delta",
     "ductal", "endothelial", "gamma", "macrophage", "quiescent_stellate"
   )), ]

   index <- which(colnames(scRNA) == "A1BG")
   allgene <- names(which(colSums(scRNA[, index:ncol(scRNA)]) > 0))
   scRNA <- scRNA[, c(colnames(scRNA[, 1:(index - 1)]), allgene)]

   for (num in 1:50) {
     print(paste0("doing: num ", num))
     path <- paste0("/rsrch5/home/biostatistics/lku/simdata/new/", num, "/spatial_cluster_coords.csv")
     coords <- read.csv(path)[, -1]
     celltype <- unique(coords$cell_type)

     column_names <- c("cell_type", "cell", allgene)
     simulated_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
     colnames(simulated_df) <- column_names

     for (type in celltype) {
       print(paste0("doing: ", type))
       RNA <- scRNA[which(scRNA$assigned_cluster == type), ]
       n <- nrow(coords[which(coords$cell_type == type), ])

       simulated_gene <- list()
       for (gene in allgene) {
         exprsn <- RNA[, gene]
         simulated_exprsn <- tryCatch({
           fit <- fitdistr(exprsn, densfun = "Negative Binomial")
           rnbinom(n, mu = fit$estimate["mu"], size = fit$estimate["size"])
         }, error = function(e) {
           rep(0, n)
         })
         simulated_gene[[gene]] <- simulated_exprsn
       }

       simulated_gene <- as.data.frame(simulated_gene)
       cell_type <- data.frame("cell_type" = rep(type, n))
       ID <- data.frame("cell" = paste0(type, ".", 1:n))
       dat <- cbind(cell_type, ID, simulated_gene)
       simulated_df <- rbind(simulated_df, dat)
     }

     ################ Assigning cell IDs to spatial coordinates ################

     simulated_df$cell_type <- as.character(simulated_df$cell_type)
     coords$cell_type <- as.character(coords$cell_type)
     celltypes <- unique(simulated_df$cell_type)
     colnames(coords) <- c("x", "y", "cell_type")

     dat_list <- list()
     set.seed(100)
     for (i in seq_along(celltypes)) {
       types <- celltypes[i]
       print(paste0("doing: ", types))
       select_coords <- coords[which(coords$cell_type == types), ]
       scrna <- simulated_df[which(simulated_df$cell_type == types), ]
       id <- sample(scrna$cell)
       select_coords$cell <- id
       dat_list[[i]] <- select_coords
     }

     assign <- do.call(rbind, dat_list)

     scdata_with_coordinates <- left_join(assign, simulated_df, by = c("cell_type", "cell"))
     rownames(scdata_with_coordinates) <- scdata_with_coordinates$cell

     write.csv(scdata_with_coordinates,
       paste0("/rsrch5/home/biostatistics/lku/simdata/new/", num, "/assigned_exprsn_df.csv"))
   }




5. Generate final expression counts for single-cell ST data
----------------------------------------------------------------

- Based on our design, we use different radius and expression fold change to make cell type specific L-R interaction more intense for certain detection. For different radius, expression values for specific cell type pairs are scaled based on fold-change settings (e.g., 5×, 10×, 15×).

- We select overlap L-R pairs according to different database from all the tools we benchmarked, where the L-R pairs file could be download from `Here <https://drive.google.com/file/d/1oG0s9-N9ufpSTjZlp9BFFdwWKLmRXLl0/view?usp=sharing>`__.

- We save the ground truth cell type specific interactions into ``_truth_LR.csv``

- Finally, we save different radius and fold change expression single-cell-level ST simulated data into ``paste0("k_",k,"_fc_",fc,"_assigned_exprsn_df.csv")``.

.. code-block:: r

   final <- read.csv("/rsrch5/home/biostatistics/lku/simdata/possible_LRpair.csv")
   final <- final[,-1]

   sim_name <- c(1:50)

   for (num in sim_name) {
     file <- read.csv(paste0("/rsrch5/home/biostatistics/lku/simdata/new/", num, "/assigned_exprsn_df.csv"))
     index <- which(colnames(file) == "cell")
     gene_spot_df <- file[, (index + 1):ncol(file)]
     rownames(gene_spot_df) <- file$cell
     meta_data <- file[, c("x", "y", "cell", "cell_type")]

     final$ligand_Exprsn <- 0
     final$receptor_Exprsn <- 0

     for (i in 1:nrow(final)) {
       if (final$ligand[i] %in% colnames(gene_spot_df)) {
         final$ligand_Exprsn[i] <- mean(gene_spot_df[, final$ligand[i]])
       }
       if (final$receptor[i] %in% colnames(gene_spot_df)) {
         final$receptor_Exprsn[i] <- mean(gene_spot_df[, final$receptor[i]])
       }
     }

     possible <- final[final$ligand_Exprsn > 0 & final$receptor_Exprsn > 0, ]

     specific_list <- list()
     for (celltype in unique(file$cell_type)) {
       sub_data <- file[file$cell_type == celltype, ]
       index <- which(colnames(sub_data) == "cell")
       sub_data_df <- sub_data[, (index + 1):ncol(sub_data)]
       rownames(sub_data_df) <- sub_data$cell
       sub_data_df <- as.data.frame(t(sub_data_df))
       sub_data_df[is.na(sub_data_df)] <- 0
       specific_list[[celltype]] <- sub_data_df
     }

     #### Set k range: this is for interacting radius ####
     k_values <- c(18, 36, 60)
     contact_results <- list()
 
     for (k in k_values) {
       knn_result <- get.knn(meta_data[, c("x", "y")], k = k)
       knn_contacts <- data.frame(cell = meta_data$cell, cell_type = meta_data$cell_type)
       knn_contacts$neighbor_cell_types <- apply(knn_result$nn.index, 1, function(neighbors) {
         paste(meta_data$cell_type[neighbors], collapse = ",")
       })
       long_contacts <- knn_contacts %>%
         tidyr::separate_rows(neighbor_cell_types, sep = ",") %>%
         count(cell_type, neighbor_cell_types, name = "contact_count")
       long_contacts$ligand_receptor_pairs <- paste0(long_contacts$cell_type, "-", long_contacts$neighbor_cell_types)
       contact_results[[paste0("k_", k)]] <- long_contacts
     }
     
     #### Then in each radius, we find the cell type involved in this radius ####
     #### This could help us define the cell type specific interactions ####
     LR <- possible[, c("ligand", "receptor", "interaction_name")]

     for (k in k_values) {
       print(paste0("doing Knn K = ", k))
       all_vec <- c()
       dff <- contact_results[[paste0("k_", k)]]

       for (n in 1:nrow(LR)) {
         ligand_gene <- LR[n, ]$ligand
         receptor_gene <- LR[n, ]$receptor

         ligand_avg_vector <- sapply(names(specific_list), function(celltype) {
           rowMeans(specific_list[[celltype]][ligand_gene, , drop = FALSE])
         })
         receptor_avg_vector <- sapply(names(specific_list), function(celltype) {
           rowMeans(specific_list[[celltype]][receptor_gene, , drop = FALSE])
         })

         ligand_avg_vector <- sort(ligand_avg_vector, decreasing = TRUE)
         receptor_avg_vector <- sort(receptor_avg_vector, decreasing = TRUE)

         ligand_positive <- names(ligand_avg_vector)[ligand_avg_vector > 0][1]
         receptor_positive <- names(receptor_avg_vector)[receptor_avg_vector > 0][1]

         if (!is.na(ligand_positive) && !is.na(receptor_positive)) {
           pair <- paste0(ligand_positive, "-", receptor_positive)
           if (pair %in% dff$ligand_receptor_pairs) {
             vec <- paste0(LR[n, ]$interaction_name, ".", pair)
             all_vec <- c(all_vec, vec)
           }
         }
       }

       all_vecdf <- data.frame(celltype = all_vec)
       all_vecdf$LR <- sapply(strsplit(as.character(all_vecdf$celltype), "\\."), `[`, 1)
       all_vecdf$cell_pair <- sapply(strsplit(as.character(all_vecdf$celltype), "\\."), `[`, 2)
       all_vecdf$ligand <- sapply(strsplit(all_vecdf$LR, "_"), `[`, 1)
       all_vecdf$receptor <- sapply(strsplit(all_vecdf$LR, "_"), `[`, 2)
       all_vecdf$celltype_A <- sapply(strsplit(all_vecdf$cell_pair, "-"), `[`, 1)
       all_vecdf$celltype_B <- sapply(strsplit(all_vecdf$cell_pair, "-"), `[`, 2)

       path1 <- paste0("/rsrch5/home/biostatistics/lku/simdata/new/", num, "/k_", k, "_truth_LR.csv")
       write.csv(all_vecdf, path1)

       ##### This is for different expression fold change ####
       ##### This means within the radius, the expression will be, for example 5 fold higher then others ####
       for (fc in c(5, 10, 15)) {
         print(paste0("doing fc = ", fc))
         knn_result <- get.knn(as.matrix(meta_data[, c("x", "y")]), k = k)[[1]]
         gene_count_new <- gene_spot_df

         for (celltype_A in unique(all_vecdf$celltype_A)) {
           cellid_A <- which(file$cell_type == celltype_A)
           for (cellid in cellid_A) {
             nearest_cells <- knn_result[cellid, ]
             for (nearid in nearest_cells) {
               nearest_celltype_B <- file[nearid, "cell_type"]
               tmp_LR <- all_vecdf[all_vecdf$celltype_A == celltype_A & all_vecdf$celltype_B == nearest_celltype_B, ]
               genes_to_scale <- unique(c(tmp_LR$ligand, tmp_LR$receptor))
               gene_ind <- na.omit(match(genes_to_scale, colnames(gene_count_new)))

               if (length(gene_ind) > 0) {
                 gene_count_new[c(cellid, nearid), gene_ind] <- fc * gene_spot_df[c(cellid, nearid), gene_ind]
               }
             }
           }
         }

         info <- file[, 1:(which(colnames(file) == "cell"))]
         count_new <- cbind(info, gene_count_new)[, -1]
         path2 <- paste0("/rsrch5/home/biostatistics/lku/simdata/new/", num, "/k_", k, "_fc_", fc, "_assigned_exprsn_df.csv")
         write.csv(count_new, path2)
       }
     }
   }






