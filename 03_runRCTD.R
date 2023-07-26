library(Matrix)
library(spacexr)
library(Seurat)
sampcor <- read.csv('/home/whou10/scratch16/whou10/boldrinin/doc/visium_snmultiome_sample.csv',as.is = T)

jid <- as.numeric(commandArgs(trailingOnly = T))
  print(jid)
  mat <- readRDS(paste0('/home/whou10/scratch16/whou10/boldrinin/visium/deconvolute/data/',tolower(sampcor[jid,2]),'/count.rds'))
  ct <- readRDS(paste0('/home/whou10/scratch16/whou10/boldrinin/visium/deconvolute/data/',tolower(sampcor[jid,2]),'/celltype.rds'))
  tab <- table(ct)
  id <- which(ct %in% names(tab)[tab >= 25])
  ct <- ct[id]
  mat <- mat[,id]
  nct <- names(ct)
  ct <- factor(as.character(ct))
  names(ct) <- nct
  
  reference <- Reference(mat, ct, colSums(mat))
  
  spa <- Load10X_Spatial(paste0('/home/whou10/scratch16/whou10/boldrinin/visium/data/',sampcor[jid,1]))
  mat <- spa@assays$Spatial@counts
  coord <- spa@images$slice1@coordinates[,c('imagerow','imagecol')]
  puck <- SpatialRNA(coord, mat, colSums(mat))
  
  myRCTD <- create.RCTD(puck, reference, max_cores = 1)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
  
  results <- myRCTD@results
  norm_weights = normalize_weights(results$weights) 
  saveRDS(norm_weights,file=paste0('/home/whou10/scratch16/whou10/boldrinin/visium/deconvolute/res/',sampcor[jid,1],'.rds'))


