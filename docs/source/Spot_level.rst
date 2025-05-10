Spot-level Simulated Data
=========================================== 

The following R script demonstrates how we generate simulated spot-level ST data.

- We use the scRNA-seq to generate the ST data where we put 12 single cell inside each spot.

- The scRNA-seq cell type and cell type counts were according to `GSM2230759 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi>`__, where the file could be downloaded from `Here <https://drive.google.com/file/d/1ME3fhDi4muYvDoSQWUVyJoGoZryV0M2k/view?usp=sharing>`__

- We design that the center region will have intense cell-cell interactions (high expression of ligand-receptor) for a more certain global detection, the center region spot IDs could be downloaded from `Here <https://drive.google.com/file/d/1W13I22aToE8BSYn6a4yi2xdpMEoFuwfI/view?usp=sharing>`__

- The spatial coordinates for constructing the ST data could be downloaded from `Here <https://drive.google.com/file/d/1_YYniZa2BZXZLcPStnjtccHSWAWYnk-W/view?usp=sharing>`__

1. Dependencies
-------------------------

.. code-block:: r


   library(MASS)
   library(dplyr)


2. Design the cell counts for each cell type based on scRNA reference
---------------------------------------------------------------------------

- Cell types contain activated_stellate, beta, delta, ductal, macrophage  quiescent_stellate 

.. code-block:: r

   # Read the scRNA-seq reference
   data3 <- read.csv("/rsrch5/home/biostatistics/lku/spot_level/test/GSM2230759_human3_umifm_counts.csv")

   # High expression cell type ids
   ma_ids <- paste0("macrophage.",c(1:10))
   dc_ids <- paste0("ductal.",c(1:400))
   de_ids <- paste0("delta.",c(1:188))
   be_ids <- paste0("beta.",c(1:841))
   as_ids <- paste0("activated_stellate.",c(1:102))
   qs_ids <- paste0("quiescent_stellate.",c(1:55))
   high_ids <- c(ma_ids,dc_ids,de_ids,be_ids,as_ids, qs_ids  )
   
   # Low expression cell type ids
   ma_ids <- paste0("macrophage.",c(11: (78+11-1) ))
   dc_ids <- paste0("ductal.",c(401:(1965+401-1) ))
   de_ids <- paste0("delta.",c(189:(825+189-1) ))
   be_ids <- paste0("beta.",c(842:(4108+842-1) ))
   as_ids <- paste0("activated_stellate.",c(103:(527+103-1) ))
   qs_ids <- paste0("quiescent_stellate.",c(56:(285+56-1) ))
   low_ids <- c(ma_ids,dc_ids,de_ids,be_ids,as_ids, qs_ids  )


3. Read the spatial coordinates, center region spot ids
---------------------------------------------------------

- Here we separate the spot ids into center region ids (where we out intense CCI , high expressed L-R), and the remain spot ids

.. code-block:: r

   interest_region_Spot_IDs <- read.csv("/rsrch5/home/biostatistics/lku/spot_level/test/interest_region_Spot_IDs.csv")
   exp_id <- interest_region_Spot_IDs$interest_region_Spot_IDs
   coord <- read.csv("/rsrch5/home/biostatistics/lku/spot_level/test/coord.csv")
   rownames(coord) <- coord$X
   coord <- coord[,-1]
   selected_region <- coord[exp_id,]
   center_row <- mean(selected_region$y)
   center_col <- mean(selected_region$x)
   selected_region$distance <- sqrt((selected_region$y - center_row)^2 + (selected_region$x - center_col)^2)
   sorted_region <- selected_region[order(selected_region$distance), ]
   sorted_ids <- rownames(sorted_region)
   exp_id <- sorted_ids
   unique_cell_types <- c("beta","delta" ,"ductal","macrophage","activated_stellate","quiescent_stellate")
   remain_id <- setdiff(rownames(coord),exp_id)


4. Generate the simulated spot-level ST data
--------------------------------------------------

- We first derive the gene expression distribution from a scRNA-seq reference and remodel it based on the designed cell type composition..

- Then we place the higher expressed counts into the center region, and the low expressed counts into the remaining region. 

-  We use different random seed to generate random simulated data.


.. code-block:: r

   for (si in 1:100) {
     set.seed(si)
     random_cell_ids <- sample(high_ids)
     high_index <- matrix(random_cell_ids, nrow = 12)

     set.seed(si)
     random_cell_ids <- sample(low_ids)
     low_index <- matrix(random_cell_ids, nrow = 12)

     all_index <- cbind(high_index, low_index)
     colnames(all_index) <- c(1:782)
     spot_prop <- data.frame(matrix(ncol = length(unique_cell_types), nrow = 782))
     colnames(spot_prop) <- unique_cell_types

     for (i in 1:ncol(all_index)) {
       print(paste("doing column", i))
       simplified_index <- data.frame("celltype" = gsub("(\\..*)", "", all_index[, i]))
       tt51 <- simplified_index %>% group_by(celltype) %>% summarise(n = n()) %>%
               mutate(freq = n / sum(n)) %>% select(celltype, freq)

       for (cluster in tt51$celltype) {
         if (cluster %in% colnames(spot_prop)) {
           spot_prop[i, cluster] <- tt51[tt51$celltype == cluster, "freq"]
         }
       }
     }

     spot_prop[is.na(spot_prop)] <- 0
     rownames(spot_prop) <- c(exp_id, remain_id)
     spot_prop <- spot_prop[rownames(coord), ]
     spot_prop$cell_type <- colnames(spot_prop)[max.col(spot_prop, ties.method = "first")]
     spot_prop$cell <- rownames(spot_prop)

     unique_cell_types <- c("beta", "delta", "ductal", "macrophage", "activated_stellate", "quiescent_stellate")
     dd <- data3[data3$assigned_cluster %in% unique_cell_types, ]
     genes <- colnames(dd[, 4:ncol(dd)])

     sample_length <- data.frame(
       "celltype" = unique_cell_types,
       "high" = c(841, 188, 400, 10, 102, 55),
       "low" = c(4108, 825, 1965, 78, 527, 285)
     )
     rownames(sample_length) <- unique_cell_types

     gene_list <- list()

     for (gene in genes) {
       print(paste("doing:", gene))
       df_big <- data.frame(matrix(ncol = 2, nrow = 0))
       colnames(df_big) <- c("expression", "gene_symbol")

       for (celltype in unique_cell_types) {
         dd_data <- data3[data3$assigned_cluster == celltype, ]
         high_len <- sample_length[celltype, "high"]
         low_len <- sample_length[celltype, "low"]
         exp <- sort(dd_data[, gene], decreasing = TRUE)

         if ((sum(exp) <= 10) | (sum(exp > 0) < 7)) {
           new_samples <- rep(0, (high_len + low_len))
           names(new_samples) <- paste0(celltype, ".", seq_along(new_samples))
         } else {
           fit_result <- tryCatch({
             fit <- fitdistr(exp, "negative binomial")
             size <- fit$estimate["size"]
             mu <- fit$estimate["mu"]
             new_samples <- rnbinom(high_len + low_len, size = size, mu = mu)
             new_samples <- sort(new_samples, decreasing = TRUE)
             names(new_samples) <- paste0(celltype, ".", seq_along(new_samples))
           }, error = function(e) {
             print(paste("Error fitting NB for gene:", gene, "cell type:", celltype))
             return(NULL)
           })

           if (is.null(fit_result)) next
         }

         high_samples <- new_samples[1:high_len]
         low_samples <- new_samples[(high_len + 1):(high_len + low_len)]
         df <- data.frame(expression = c(high_samples, low_samples), gene_symbol = gene)
         df_big <- rbind(df_big, df)
       }

       gene_list[[gene]] <- df_big
     }

     gs_exp <- data.frame(matrix(ncol = 782, nrow = length(genes)))
     rownames(gs_exp) <- names(gene_list)

     for (i in 1:ncol(all_index)) {
       print(paste("doing row", i))
       for (gene in genes) {
         tt5 <- gene_list[[gene]]
         tt51 <- tt5[all_index[, i], ]
         gs_exp[gene, i] <- sum(tt51$expression)
       }
     }

     colnames(gs_exp) <- c(exp_id, remain_id)
     gs_exp <- gs_exp[, rownames(coord)]
     gs_exp[is.na(gs_exp)] <- 0
     gene_spot <- as.data.frame(t(gs_exp))
     gene_spot$cell <- rownames(gene_spot)

     sel_sums <- function(df, selected_rows, selected_cols) {
       selected_df <- df[selected_rows, selected_cols, drop = FALSE]
       return(colSums(selected_df))
     }

     specific_list <- list()
     for (celltype in unique_cell_types) {
       print(paste("doing celltype", celltype))
       cell_exp <- data.frame(matrix(ncol = 782, nrow = length(genes)))
       rownames(cell_exp) <- names(gene_list)

       for (i in 1:ncol(all_index)) {
         cc <- sub("\\..*", "", all_index[, i])
         ind <- all_index[, i][cc == celltype]
         mm <- lapply(gene_list, sel_sums, selected_rows = ind, selected_cols = "expression")
         cell_exp[, i] <- as.numeric(unlist(mm))
       }

       colnames(cell_exp) <- c(exp_id, remain_id)
       cell_exp <- cell_exp[, rownames(coord)]
       cell_exp[is.na(cell_exp)] <- 0
       specific_list[[celltype]] <- cell_exp
     }

     df <- left_join(coord, spot_prop, by = "cell")
     df <- left_join(df, gene_spot, by = "cell")

     save(specific_list, spot_prop, gs_exp, file = "/rsrch5/home/biostatistics/lku/spot_level/test/spotdata.rda")

     path_name <- paste0("/Users/lku/Desktop/CCI/big_sim/sim", si, "/sim_df.csv")
     write.csv(df, path_name)
   }
