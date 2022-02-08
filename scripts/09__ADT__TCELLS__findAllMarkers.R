DefaultAssay(sobj.sct) <- "RNA"
DimPlot(sobj.sct)
sobj.sct$seurat_clusters <- sobj.sct$SCT_snn_res.0.4
Idents(sobj.sct) <- sobj.sct$seurat_clusters

sobj.sct.all.markers <- FindAllMarkers(sobj.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(sobj.sct.all.markers, file = "results/ADT__TCELLS__findAllMarkers_minpct.25_lfc.25.RDS") 


