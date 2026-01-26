CellNEST Analysis of CCI
=========================================== 

The following script demonstrates how we use CellNEST to conduct the CCI analysis.

- The CellNEST documentation, tutorial and package could be found here: https://github.com/schwartzlab-methods/CellNEST

- The paper could be found here: https://doi.org/10.1038/s41592-025-02721-3

CellNEST could only be run on GPU node. Here we demonstrate how to run on High Computing Cluster (HPC) using singularity by submitting LSF jobs.

1. Preprocessing Step
-------------------------

.. code-block:: r

   #BSUB -J CellNEST_IPMN_preprocess
   #BSUB -W 12:00
   #BSUB -q gpu
   #BSUB -gpu "num=1:gmem=16"
   #BSUB -n 4
   #BSUB -R "rusage[mem=16GB]"
   #BSUB -M 16
   #BSUB -o logs/cellnest_preprocess.%J.out
   #BSUB -e logs/cellnest_preprocess.%J.err
   #BSUB -cwd /home/lku/CellNEST
   #BSUB -B 
   #BSUB -N 
   #BSUB -u lku@mdanderson.org

   module load singularity

   singularity exec --nv /home/lku/cellnest_image.sif \
     bash cellnest preprocess \
     --data_name='IPMN_LG2' \ 
     --data_from='data/IPMN_LG2/'


2. Run the analysis 
-------------------------

- We follow the SpaCCI guidance and tutorial for the analysis

.. code-block:: r

   folder_names <- c("slice1")
   for (num in folder_names) {
     print(paste0("doing dataset: ", num))
     path <- paste0("/rsrch5/home/biostatistics/lku/ILIBD/data/", num, "/exprsn_df.csv")
     data <- read.csv(path)
     data <- data[!is.na(data$cell_type), ]

     # Prepare the gene expression dataframe
     index <- which(colnames(data)=="cell_type")
     gene_spot_df <- data[,c((index+1) : ncol(data))]
     rownames(gene_spot_df) <- data$cell
     gene_spot_df <- as.data.frame(t(gene_spot_df))
     library.size <- Matrix::colSums(gene_spot_df)
     gene_spot_df <- as.data.frame(Matrix::t(Matrix::t(gene_spot_df)/library.size) *10000) # normalized data matrix
     gene_spot_df[is.na(gene_spot_df)] <- 0
     
     # Prepare the cell type proportion data frame
     cell_index <- which(colnames(data)=="cell")
     celltype_index <- which(colnames(data)=="cell_type")
     cell_prop_df <- data[,c((cell_index+1):(celltype_index-1))]
     rownames(cell_prop_df) <- data$cell
     spatial_coords_df <- data[,c("x","y"),drop=FALSE]
     colnames(spatial_coords_df) <- c("imagecol", "imagerow")  
     rownames(spatial_coords_df) <- file$cell
     spatial_coords_df$Spot_ID <- rownames(spatial_coords_df)


     # Extract the L-R database, here we use the CellChat human database (you could explore other options)
     LRDB <- LR_database("Human", "CellChat", gene_spot_df)

     # Run the global analysis
     result_global <- run_SpaCCI(gene_spot_expression_dataframe = gene_spot_df,
                                 spot_cell_proportion_dataframe = cell_prop_df,
                                 spatial_coordinates_dataframe = spatial_coords_df,
                                 LR_database_list = LRDB,
                                 analysis_scale = "global")
      
     df <- result_global1$pvalue_df
     df <- df[,c("Cell_type_Ligand" ,"Cell_type_Receptor", "strength" ,"adjusted.PValue" ,"interaction_name","ligand.symbol","receptor.symbol")]
     colnames(df) <- c("ligand" ,  "receptor" ,"strength"  ,  "pvalue"  , "interaction_name" ,"ligand_gene","receptor_gene")
  
     write.csv(df,paste0("/rsrch5/home/biostatistics/lku/ILIBD/data/",num,"/spacci_result.csv"))
 

   }




3. Post-Processing Step   
-------------------------





