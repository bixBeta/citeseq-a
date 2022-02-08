# ---------------------------------------
# ADT  features from gene annotation LV excel sheet
lv.feats = c("CD3E","CD4","CD8A","TRAC","TRBV2","TRBC1","TRBC2","TRGC2","TRGC1","TRDC","FCGR3B","NCAM1","
              PTPRC","FAS","CCR7","CXCR3","CCR4","KLRB1","CCR6","CCR10","IL2RA", "IL7R", "LEF1", "SELL", "FOXP1", "GZMK", "GZMM")

DefaultAssay(sobj.sct) <- "ADT"
sobj.sct = NormalizeData(sobj.sct, assay = "ADT", normalization.method = "CLR", margin = 2)
sobj.sct = FindVariableFeatures(sobj.sct)
sobj.sct = ScaleData(sobj.sct)
rownames(sobj.sct)

adt.data <- GetAssayData(sobj.sct, slot = "data",assay = "ADT")
adt.dist <- dist(t(adt.data))
sobj.sct = RunPCA(sobj.sct, reduction.name = "pca_adt", reduction.key = "pca_adt_")

FeaturePlot(sobj.sct, features = rownames(sobj.sct), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)

sobj.sct[["rnaClusterID"]] <- Idents(sobj.sct)

# Now, we rerun tSNE using our distance matrix defined only on ADT (protein) levels.
sobj.sct[["umap_adt"]] <- RunUMAP(adt.dist, assay = "ADT", reduction.key = "adtUMAP_")
sobj.sct[["adt_snn"]] <- FindNeighbors(adt.dist)$snn
sobj.sct <- FindClusters(sobj.sct, resolution = 0.2, graph.name = "adt_snn")

clustering.table <- table(Idents(sobj.sct), sobj.sct$rnaClusterID)
clustering.table

png(filename = "figures_ADT_TCELLS/dotPlot-adtAssay-groupBy-gex-clusters.png", width = 1000, height = 1500)
DotPlot(sobj.sct, assay="ADT", features=rownames(sobj.sct), col.min=-1.5,
        group.by = "rnaClusterID", split.by="orig.ident",cols="Spectral", dot.scale=3) + RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("dotPlot-adtAssay-groupBy-gex-clusters") 
dev.off()

DefaultAssay(sobj.sct)  <- "RNA"

png(filename = "figures_ADT_TCELLS/dotPlot-gex-clusters.png", width = 1000, height = 1500)
DotPlot(sobj.sct, assay="RNA", features=lv.feats, col.min=-1.5,
        group.by = "rnaClusterID", split.by="orig.ident",cols="Spectral", dot.scale=3) + RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("dotPlot-gex-clusters") 
dev.off()


###--seurat vignette
tsne_rnaClusters <- DimPlot(sobj.sct, reduction = "umap_adt", group.by = "rnaClusterID") + NoLegend()
tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Clustering based on scRNA-seq") + theme(plot.title = element_text(hjust = 0.5))
tsne_rnaClusters <- LabelClusters(plot = tsne_rnaClusters, id = "rnaClusterID", size = 7)
tsne_rnaClusters

tsne_adtClusters <- DimPlot(sobj.sct, reduction = "umap_adt", pt.size = 0.5) + NoLegend()
tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on ADT signal") + theme(plot.title = element_text(hjust = 0.5))
tsne_adtClusters <- LabelClusters(plot = tsne_adtClusters, id = "ident", size = 7)
wrap_plots(list(tsne_rnaClusters, tsne_adtClusters), ncol = 2)

fplot.RNA = FeaturePlot(sobj.sct, features = lv.feats, min.cutoff = "q05", max.cutoff = "q95", ncol = 5)

DefaultAssay(sobj.sct) <- "ADT"
fplot.adt = FeaturePlot(sobj.sct, features = rownames(sobj.sct), min.cutoff = "q05", max.cutoff = "q95", ncol = 4, reduction = "umap")


saveRDS(sobj.sct, "data/ADT__TCELLS__ADT__Clustering.Rds")

Idents(sobj.sct) <- sobj.sct$adt_snn_res.0.2
FeaturePlot(sobj.sct, features = rownames(sobj.sct), min.cutoff = "q05", max.cutoff = "q95",
            ncol = 5, reduction = "umap_adt", label = T, label.size = 5) 

DefaultAssay(sobj.sct) <- "RNA"
FeaturePlot(sobj.sct, features = lv.feats, min.cutoff = "q05", max.cutoff = "q95",
            ncol = 5, reduction = "umap_adt", label = T, label.size = 5) 


Idents(sobj.sct) <- sobj.sct$SCT_snn_res.0.4
FeaturePlot(sobj.sct, features = lv.feats, min.cutoff = "q05", max.cutoff = "q95",
            ncol = 5, reduction = "umap", label = T, label.size = 5) 

DefaultAssay(sobj.sct) <- "ADT"
FeaturePlot(sobj.sct, features = rownames(sobj.sct), min.cutoff = "q05", max.cutoff = "q95",
            ncol = 5, reduction = "umap", label = T, label.size = 5) 



lv.feats2 = c("CD3E","CD247","CD4","CD8A","TRBC1","TRGC1","TRDC",
              "CCR7","CD28","SELL","CD69","CCL5","IFNG","GATA3","IRF4",
              "SPI1","KLRB1","AHR","CTLA4","FOXP3","NCAM1","FCGR3A","KLRK1",
              "NCR3","LILRB1","CD14","LGALS2","S100A8","PF4","MS4A1","CD1C",
              "HLA-DRA","ITGAX","CD34")

s = DimPlot(sobj.sct, label = T, label.size = 5, ncol = 4, split.by = "SCT_snn_res.0.4")
s1 = DimPlot(sobj.sct, label = T, label.size = 8, ncol = 1)
png(filename = "ADT_TCELLS_GEX_UMAP_0.4.png", width = 3200, height = 1080, res = 150)
s1 + s 
dev.off()

# new marker list 
DefaultAssay(sobj.sct) <- "RNA"
Idents(sobj.sct) <- sobj.sct$SCT_snn_res.0.4
FeaturePlot(sobj.sct, features = lv.feats2, min.cutoff = "q05", max.cutoff = "q95",
            ncol = 5, reduction = "umap", label = T, label.size = 5) 



DefaultAssay(sobj.sct)  <- "RNA"
Idents(sobj.sct)
png(filename = "figures_ADT_TCELLS/lv2__dotPlot-gex-clusters.png", width = 1000, height = 1500)
DotPlot(sobj.sct, assay="RNA", features=lv.feats2, col.min=-1.5,
        group.by = "rnaClusterID", split.by="orig.ident",cols="Spectral", dot.scale=3) + RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("dotPlot-gex-clusters") 
dev.off()

DefaultAssay(sobj.sct)  <- "ADT"
Idents(sobj.sct)
png(filename = "figures_ADT_TCELLS/lv2__dotPlot-adtAssay-groupBy-gex-clusters.png", width = 1000, height = 1500)
DotPlot(sobj.sct, assay="ADT", features=rownames(sobj.sct), col.min=-1.5,
        group.by = "rnaClusterID", split.by="orig.ident",cols="Spectral", dot.scale=3) + RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("dotPlot-adtAssay-groupBy-gex-clusters") 
dev.off()
