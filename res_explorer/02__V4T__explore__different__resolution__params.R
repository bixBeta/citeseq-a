v4t <- readRDS("/workdir/DOCKER_DUMP_LM31/data/TCELLS_sct__harmony__ndim50__res0.4_0.5_0.6_0.8.RDS")
sanity.meta =v4t@meta.data

all.res = seq(0.4,2,.1)
DefaultAssay(v4t) <- "SCT"
v4t.obj.added.clusts = FindClusters(v4t, resolution = all.res)
v4t.obj.added.clusts.meta = v4t.obj.added.clusts@meta.data

table(v4t.obj.added.clusts$SCT_snn_res.0.6 == sanity.meta$SCT_snn_res.0.6)

dotplots.list = list()

lvfeats = c("CD3E","CD247","CD4","CD3","CD8","CD8A","TRBC1","TRGC1","TRDC","CCR7","CD28","SELL","CD69","CCL5","IFNG","EOMES","TBX21","CXCR3","IL2","GATA3","IL4R","CCR4","PTGDR2","IRF4","SPI1","KLRB1","CCR6","RORC","AHR","FOXO4","CCR10","CTLA4","FOXP3","IL2RA","LAG3","IL10","NCAM1","FCGR3A","KLRK1","NCR3","LILRB1","TRPM3")

getUmapDotPlots = function(sobj, res, feats){
  
  Idents(object = sobj) <- sobj@meta.data[,paste0("SCT_snn_res.", res)]
  uplot = DimPlot(sobj , group.by = paste0("SCT_snn_res.", res), label = T, label.size = 6, repel = T) + ggtitle(paste0("sv4 Tcells UMAP -- res ", res))
  
  dplot = DotPlot(object = sobj, assay="RNA", features=feats, col.min=-1.5,
                  group.by = paste0("SCT_snn_res.", res), cols="Spectral", dot.scale=3, cluster.idents = T, scale.by = 'size') + 
    RotatedAxis() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggtitle(paste0("sv4 Tcells DotPlot -- res ", res))
  
  cplot = uplot + dplot
  
  
  png(filename = paste0("res_explorer/umap_dotplots_sv4t/sv4-Tcells__umap__dotPlot__res__",res, ".png"), width = 1800, height = 700, res = 100)
  print(cplot)
  #Sys.sleep(3)
  dev.off()
  
  # return(cplot)
  
}

for (i in 1:length(all.res)) {
  getUmapDotPlots(sobj = v4t.obj.added.clusts, res = all.res[i], feats = lvfeats)
}

saveRDS(object = v4t.obj.added.clusts, file = "res_explorer/v4t.object.allRES.RDS")


DefaultAssay(v4t.obj.added.clusts) <- "RNA"

png("res_explorer/featurePlots-SV4-Tcells-1.png", width = 1200, height = 1080, res = 100)
FeaturePlot(v4t.obj.added.clusts, features = lvfeats[1:20], raster = F)
dev.off()


png("res_explorer/featurePlots-SV4-Tcells-2.png", width = 1200, height = 1080, res = 100)
FeaturePlot(v4t.obj.added.clusts, features = lvfeats[21:42], raster = F)
dev.off()





