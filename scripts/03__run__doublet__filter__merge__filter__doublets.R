source("packages.R")
# ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
# sweep.res.list <- paramSweep_v3(R2339, PCs = 1:10, sct = TRUE)
# sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
# bcmvn <- find.pK(sweep.stats)
# 
# ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
# annotations = R2339@meta.data$SCT_snn_res.0.6
# homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- R2339@meta.data$ClusteringResults
# nExp_poi <- round(R2339$expectedRate[[1]]*nrow(R2339@meta.data))  ## Assuming x% doublet formation rate - tailor for your dataset
# nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# 
# ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
# pK_value <- as.numeric( as.vector( bcmvn$pK )[ which.max( bcmvn$BCmetric ) ] )
# R2339 <- doubletFinder_v3(R2339, PCs = 1:10, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
# R2339.2 <- doubletFinder_v3(R2339, PCs = 1:10, pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.18_106", sct = T)
# 
# 
# DimPlot( R2339 , reduction = "umap" , group.by = colnames( R2339@meta.data )[ length(R2339@meta.data) ] , cols=c("khaki","blue") , order=c('Doublet','Singlet') )
# 

#____________________________________________________________________________

#---------------------------------------------------------------------------------------------------
for ( i in 1:length(sobj.list) ){
  
  seuratObject = sobj.list[[i]]
  library_id <- names(sobj.list)[[i]]
  print( paste( i , ' ==> ' , 'Library ID: ' , library_id ) )
  
  ######## pK Identification (no ground-truth) #############################################################
  sweep.res.list <- paramSweep_v3( seuratObject , PCs = 1:10 , sct = T )
  sweep.stats <- summarizeSweep( sweep.res.list , GT = FALSE )
  bcmvn <- find.pK( sweep.stats )
  pK_value <- as.numeric( as.vector( bcmvn$pK )[ which.max( bcmvn$BCmetric ) ] )
  
  ######## Homotypic Doublet Proportion Estimate ###########################################################
  
  nExp_poi     <- round(seuratObject$expectedRate[[1]]*nrow(seuratObject@meta.data))  ## Assuming doublet formation rate - tailor for your dataset
  ######## Run DoubletFinder with varying classification stringencies ######################################
  
  seuratObject <- doubletFinder_v3( seuratObject , PCs = 1:10 , pN = 0.25 , pK = pK_value , nExp = nExp_poi     , reuse.pANN = FALSE , sct = T )
  ######## Write down ######################################################################################
  
  OutputTab <- paste( "results/doubletFinder-tsvs/" , library_id , '.DoubletFinder.tsv' , sep='' )
  write.table( seuratObject@meta.data  ,  file = OutputTab  ,  sep='\t' )
  OutputPdf <- paste( "results/" , library_id , '.DoubletFinder.umap.pdf' , sep='' )
  png( file = OutputPdf , width = 8.5 , height = 8 )
  DimPlot( seuratObject , reduction = "umap" , group.by = colnames( seuratObject@meta.data )[ length(seuratObject@meta.data) ] , cols=c("khaki","blue") , order=c('Doublet','Singlet') )
  dev.off()
  Sys.sleep(3)
  ######## Clear objects ###################################################################################
  
}

df.list= list()
df.paths = list.files("results/doubletFinder-tsvs/", pattern = ".tsv$", full.names = T)

for (i in 1:length(df.paths)) {
  df.list[[i]] = read.delim(df.paths[i], header = T, sep = "\t")
  names(df.list)[[i]] = basename(df.paths)[[i]]
}


doubletBarcodes = sapply(df.list, function(x){
  
  rownames(x[which(x[14] == "Doublet"),])
  
})

bcs = unname(unlist(doubletBarcodes))

sobj.sct <- merge(x = sobj.list[[1]], y = sobj.list[2:length(sobj.list)], merge.data=TRUE)
sobj.sct.meta = sobj.sct@meta.data


cells.use = rownames(sobj.sct.meta)[which(!rownames(sobj.sct.meta) %in% bcs)]

adt.obj = subset(sobj.sct, cells = cells.use)
DimPlot(adt.obj)

save.image("data/doublet-filtered-adt-state.Rdata")
