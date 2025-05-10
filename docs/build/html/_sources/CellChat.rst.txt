CellChat_v2 Analysis of CCI
=========================================== 

The following R script demonstrates how we use CellChat_v2 conduct the CCI analysis.

- The CellChat_v2 documentation and package could be found here: https://github.com/jinworks/CellChat

- The paper could be found here: https://doi.org/10.1038/s41596-024-01045-4


1. Dependencies
-------------------------

.. code-block:: r


   library(CellChat)
   library(patchwork)
   options(stringsAsFactors = FALSE)


2. Run the analysis
-------------------------

- We follow the CellChat_v2 guidance and tutorial for the analysis.


.. code-block:: r

   num <- "slice1"
   path <- paste0("/rsrch5/home/biostatistics/lku/ILIBD/data/",num,"/exprsn_df.csv")
   file <- read.csv(path)
   file <- file[!is.na(file$cell_type), ]
   index <- which(colnames(file)=="cell_type")
   data.input <- file[,c((index+1) : ncol(file))]
   rownames(data.input) <- file$cell
   data.input <- as.data.frame(t(data.input))
   library.size <- Matrix::colSums(data.input)
   data.input <- log1p(Matrix::t(Matrix::t(data.input)/library.size) *10000) # normalized data matrix
   data.input[is.na(data.input)] <- 0
   meta <- file[,c("x","y","cell","cell_type")]
   rownames(meta) <- meta$cell
   spatial.locs <- as.matrix(meta[,c("x","y")])
   spot.size = 65 # the theoretical spot size (um) in 10X Visium
   conversion.factor = spot.size/24.56651
   spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

   cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type", datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
   CellChatDB <- CellChatDB.human
   CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
   cellchat@DB <- CellChatDB.use
   cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
   future::plan("multisession", workers = 4) 
   cellchat <- identifyOverExpressedGenes(cellchat)
   cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
   cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                                 distance.use = TRUE, interaction.range = 250, scale.distance = 0.1,
                                 contact.dependent = TRUE, contact.range = 100)
   cellchat <- filterCommunication(cellchat, min.cells = 10)
   cellchat <- computeCommunProbPathway(cellchat)
   cellchat <- aggregateNet(cellchat)
   df.net <- subsetCommunication(cellchat1,thresh = 1.1)
   df.net$ligand <- df.net$source
   df.net$receptor <- df.net$target
   df.net$pvalue <- df.net$pval
   df.net$value <- df.net$prob
   df.net <- df.net[,c("interaction_name","pathway_name", "ligand","receptor","value","pvalue")]
   df.net <- df.net[,c("ligand" ,  "receptor" ,"value"  ,  "pvalue"  , "interaction_name" )]
   colnames(df.net) <- c("ligand" ,  "receptor" ,"score"  ,  "pvalue"  , "interaction_name" )
   write.csv(df.net ,paste0("/rsrch5/home/biostatistics/lku/ILIBD/data/",num,"/cellchat_result.csv"))






