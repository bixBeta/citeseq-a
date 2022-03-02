DefaultAssay(sobj.sct) <- "RNA"
DimPlot(sobj.sct)
sobj.sct$seurat_clusters <- sobj.sct$SCT_snn_res.0.4
Idents(sobj.sct) <- sobj.sct$seurat_clusters

sobj.sct.all.markers <- FindAllMarkers(sobj.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(sobj.sct.all.markers, file = "results/ADT__TCELLS__findAllMarkers_minpct.25_lfc.25.RDS") 


DefaultAssay(adt.tcells) <- "RNA"
Idents(object = adt.tcells) <- adt.tcells$SCT_snn_res.0.4
res0.4 = DimPlot(adt.tcells , group.by = "SCT_snn_res.0.4", label = T, label.size = 5) + ggtitle("res 0.4")

Idents(object = adt.tcells) <- adt.tcells$SCT_snn_res.0.5
res0.5 = DimPlot(adt.tcells , group.by = "SCT_snn_res.0.5", label = T, label.size = 5) + ggtitle("res 0.5")

Idents(object = adt.tcells) <- adt.tcells$SCT_snn_res.0.6
res0.6 = DimPlot(adt.tcells , group.by = "SCT_snn_res.0.6", label = T, label.size = 5) + ggtitle("res 0.6")

Idents(object = adt.tcells) <- adt.tcells$SCT_snn_res.0.8
res0.8 = DimPlot(adt.tcells , group.by = "SCT_snn_res.0.8", label = T, label.size = 5) + ggtitle("res 0.8")

png(filename = "figures/TCELLS__sct__harmony__umaps.png", width = 1920, height = 1080)
res0.4 + res0.5 + res0.6 + res0.8
dev.off()


adt.tcells.meta = adt.tcells@meta.data
adt.tcells$seurat_clusters = adt.tcells$SCT_snn_res.0.8

DimPlot(adt.tcells, label = T)
lv3 = c("CD3E","CD247","CD4","CD8A","TRBC1","TRGC1","TRDC","CCR7","CD28","SELL","CD69","CCL5","IFNG","GATA3","IRF4","SPI1","KLRB1","AHR","CTLA4","FOXP3","NCAM1","FCGR3A","KLRK1","NCR3","LILRB1","CD14","LGALS2","S100A8","PF4","MS4A1","CD1C","HLA-DRA","ITGAX","CD34")

png(filename = "figures_ADT_TCELLS/res0.8-dotPlot-lv2-MEGA.png", width = 1200, height = 800, res = 150)
DotPlot(adt.tcells, assay="RNA", features=lv3, col.min=-1.5,
        group.by = "SCT_snn_res.0.8",cols="Spectral", dot.scale=3, cluster.idents = T, scale.by = 'size') + RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("res0.8 dotPlot-gex-clusters") 
dev.off()

DefaultAssay(adt.tcells) <- "ADT"
png(filename = "figures_ADT_TCELLS/res0.8-ADTdotPlot-lv2-MEGA.png", width = 1200, height = 800, res = 150)
DotPlot(adt.tcells, assay="ADT", features=rownames(adt.tcells), col.min=-1.5,
        group.by = "SCT_snn_res.0.8",cols="Spectral", dot.scale=3, cluster.idents = T, scale.by = 'size') + RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("res0.8 ADT dotPlot gex clusters") 
dev.off()

