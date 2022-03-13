# read in cite-seq object with both GEX and ADT clustering
adt.tcells <- readRDS("/home/MAIN2/February_03_CITEseq_fixed_AB_annot/data/ADT__TCELLS__ADT__Clustering.Rds")
sanity.meta = adt.tcells@meta.data

all.res = seq(0.4,2,.1)
DefaultAssay(adt.tcells) <- "SCT"
adt.obj.added.clusts = FindClusters(adt.tcells, resolution = all.res)
adt.obj.added.clusts.meta = adt.obj.added.clusts@meta.data

table(adt.obj.added.clusts$SCT_snn_res.0.6 == sanity.meta$SCT_snn_res.0.6)

dotplots.list = list()

lvfeats = c("CD3E","CD247","CD4","CD8A","TRBC1","TRGC1","TRDC","CCR7","CD28","SELL","CD69","CCL5","IFNG","EOMES","TBX21","CXCR3","IL2","GATA3","IL4R","CCR4","PTGDR2","IRF4","SPI1","KLRB1","CCR6","RORC","AHR","FOXO4","CCR10","CTLA4","FOXP3","IL2RA","LAG3","IL10","NCAM1","FCGR3A","KLRK1","NCR3","LILRB1","TRPM3")

getUmapDotPlots = function(sobj, res, feats){
  
  Idents(object = sobj) <- sobj@meta.data[,paste0("SCT_snn_res.", res)]
  uplot = DimPlot(sobj , group.by = paste0("SCT_snn_res.", res), label = T, label.size = 6, repel = T) + ggtitle(paste0("CITE-seq GEX UMAP -- res ", res))
  
  dplot = DotPlot(object = sobj, assay="RNA", features=feats, col.min=-1.5,
          group.by = paste0("SCT_snn_res.", res), cols="Spectral", dot.scale=3, cluster.idents = T, scale.by = 'size') + 
          RotatedAxis() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
          ggtitle(paste0("CITE-seq GEX DotPlot -- res ", res))
  
  cplot = uplot + dplot
  
  
  png(filename = paste0("res_explorer/umap_dotplots/CITE-seq__umap__dotPlot__res__",res, ".png"), width = 1800, height = 600, res = 100)
  print(cplot)
  #Sys.sleep(3)
  dev.off()
  
  # return(cplot)
  
}

for (i in 1:length(all.res)) {
  getUmapDotPlots(sobj = adt.obj.added.clusts, res = all.res[i], feats = lvfeats)
}

saveRDS(object = adt.obj.added.clusts, file = "res_explorer/adt.tcells.object.allRES.RDS")

