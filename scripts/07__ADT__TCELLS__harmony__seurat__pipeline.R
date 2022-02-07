sobjs.list.rna.filtered = SplitObject(object = adt.Tcells, split.by = "orig.ident")

sobjs.list.rna.filtered <- lapply(X = sobjs.list.rna.filtered, 
                                  FUN = function(x){SCTransform(object = x, method = "glmGamPoi", vars.to.regress = "percent.mt", 
                                                                return.only.var.genes = FALSE)})

var.features <- SelectIntegrationFeatures(object.list = sobjs.list.rna.filtered, nfeatures = 3000)

sobj.sct <- merge(x = sobjs.list.rna.filtered[[1]], y = sobjs.list.rna.filtered[2:length(sobjs.list.rna.filtered)], merge.data=TRUE)

VariableFeatures(sobj.sct) <- var.features

sobj.sct

sobj.sct <- RunPCA(sobj.sct, verbose = FALSE, dims = 1:50)

sobj.sct <- RunHarmony(sobj.sct, assay.use="SCT", group.by.vars = "orig.ident", plot_convergence = T)

sobj.sct <- RunUMAP(sobj.sct, reduction = "harmony", dims = 1:50)

sobj.sct <- FindNeighbors(sobj.sct, reduction = "harmony", dims = 1:50) %>% FindClusters(resolution = c(0.4, 0.5, 0.6, 0.8))

Idents(object = sobj.sct) <- sobj.sct$SCT_snn_res.0.4
res0.4 = DimPlot(sobj.sct , group.by = "SCT_snn_res.0.4", label = T, label.size = 5) + ggtitle("res 0.4")

Idents(object = sobj.sct) <- sobj.sct$SCT_snn_res.0.5
res0.5 = DimPlot(sobj.sct , group.by = "SCT_snn_res.0.5", label = T, label.size = 5) + ggtitle("res 0.5")

Idents(object = sobj.sct) <- sobj.sct$SCT_snn_res.0.6
res0.6 = DimPlot(sobj.sct , group.by = "SCT_snn_res.0.6", label = T, label.size = 5) + ggtitle("res 0.6")

Idents(object = sobj.sct) <- sobj.sct$SCT_snn_res.0.8
res0.8 = DimPlot(sobj.sct , group.by = "SCT_snn_res.0.8", label = T, label.size = 5) + ggtitle("res 0.8")


res0.4 + res0.5 + res0.6 + res0.8

Idents(object = sobj.sct) <- sobj.sct$SCT_snn_res.0.4
sobj.sct$seurat_clusters = sobj.sct$SCT_snn_res.0.4

saveRDS(sobj.sct, file = "data/ADT__TCELLS__harmony__ndim50__res0.4_0.5_0.6_0.8.RDS")

sct.meta = sobj.sct@meta.data
write.csv(sct.meta, file = "results/ADT__TCELLS__sct.meta.postIntegration.csv", quote = F)

DimPlot(sobj.sct, label = T)
