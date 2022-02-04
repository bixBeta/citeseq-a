source("packages.R")

sobj.list <- lapply(X = sobjs.list.rna.filtered, FUN = function(x){
  x = FindVariableFeatures(x)
  x = SCTransform(object = x, method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = FALSE)
  
  x = RunPCA(x, verbose = FALSE, dims = 1:50)
  x = RunUMAP(x, dims = 1:50)
  
  x = FindNeighbors(x, dims = 1:50)
  x = FindClusters(x, resolution = 0.4)
  
})

#-------------------------------------------------------------------------------------------------------------
m= 7.6652E-6
b= -3.57897E-07
sobj.list = lapply(sobj.list, function(x){
  AddMetaData(object = x, metadata = round(m * nrow(x@meta.data) + b, 2), col.name = "expectedRate"  )
})

#-------------------------------------------------------------------------------------------------------------
saveRDS(sobj.list, "data/citeSeq__INPUT_DOUBLET_FINDER_split__by__orig.ident.Rds")

